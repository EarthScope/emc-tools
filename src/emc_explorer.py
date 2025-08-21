#!/usr/bin/env python
import sys, os, getopt, traceback
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from textwrap import fill
from pyproj import Geod

# ---- Imports from your project layout ----
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
PROP_DIR = os.path.join(ROOT_DIR, "prop")
LIB_DIR = (
    os.path.join(ROOT_DIR, "lib")
    if os.path.exists(os.path.join(ROOT_DIR, "lib"))
    else os.path.join(ROOT_DIR, "src")
)
plot_var = None

sys.path.extend([PROP_DIR, LIB_DIR])

import shared_prop as prop
import tool_prop as slicer_prop
import shared_lib as lib
import tool_lib as slicer_lib

from typing import Mapping

# ---- Constants / Globals trimmed ----
VALID_FILE_TYPES = {"csv": ".csv", "geocsv": ".csv", "netcdf": ".nc"}
logger = lib.get_logger()
DASH = "-" * 10

from contextlib import contextmanager
import time


# ---------- Messaging Helpers ----------
def header(text: str) -> None:
    logger.info("\n" + "=" * 72)
    logger.info(text)
    logger.info("=" * 72)


def subheader(text: str) -> None:
    logger.info("\n" + text)
    logger.info("-" * len(text))


def tip(text: str) -> None:
    logger.info(f"[TIP] {text}")


def note(text: str) -> None:
    logger.info(f"[INFO] {text}")


def ok(text: str) -> None:
    logger.info(f"[OK] {text}")


def warn(text: str) -> None:
    logger.warning(f"[WARN] {text}")


def err(text: str) -> None:
    logger.error(f"[ERR] {text}")


def kv(key: str, value) -> None:
    logger.info(f"  • {key}: {value}")


def fmt_range(lo, hi) -> str:
    return f"[{lo:.3g}, {hi:.3g}]"


def list_menu(title: str, items: List[str]) -> None:
    subheader(title)
    for it in items:
        logger.info(f"  - {it}")


@contextmanager
def timed(section: str, enabled: bool):
    t0 = time.perf_counter()
    try:
        if enabled:
            subheader(section)
            tip("Working…")
        yield
    finally:
        if enabled:
            dt = (time.perf_counter() - t0) * 1000
            ok(f"{section} done in {dt:.0f} ms")


def _format_units(u: Optional[str]) -> Optional[str]:
    """
    - None  -> truly missing (return None so caller can decide to warn/error)
    - ""    -> dimensionless (return a friendly label)
    - other -> pass through
    """
    if u is None:
        return None
    if u == "":
        return "(dimensionless)"
    return u


# ---------- Hint system ----------
from typing import Iterable


def show_hints(hints: Dict[str, str], title: str = "Hints") -> None:
    """Pretty-print a hint table."""
    logger.info("\n" + title)
    logger.info("-" * len(title))
    # keep stable order: commands first, then alpha
    for key in sorted(hints.keys()):
        logger.info(f"  {key:<10} {hints[key]}")


def choose_with_hints(
    prompt: str,
    valid: Iterable[str],
    hints: Optional[Dict[str, str]] = None,
    allow_blank: bool = False,
) -> str:
    """
    Input loop that supports 'h'/'?' to display hints, and validates choice.
    Returns the user's choice (may be '' if allow_blank=True).
    """
    valid_set = set(v.lower() for v in valid)
    while True:
        choice = input(prompt).strip()
        low = choice.lower()
        if low in ("h", "?"):
            if hints:
                show_hints(hints)
            else:
                logger.info("No hints available.")
            continue
        if allow_blank and choice == "":
            return choice
        if choice and low in valid_set:
            return low
        # If you want to echo choices:
        logger.warning(f"[WARN] Unrecognized option: '{choice}'. Type 'h' for help.")


MAIN_HINTS = {
    "meta": "View global attributes, coords, and variables",
    "range": "Show min/max for coords and variables",
    "subset": "Enter sub-menu: volume/slice/xsection/surface",
    "map": "Show coverage map (if lat/lon available)",
    "help": "Print usage text",
    "exit": "Quit the tool",
}


SUBSET_HINTS_2D = {
    "volume": "Extract a subvolume (limits for each dimension)",
    "surface": "Surface plot of a 2D variable",
    "back": "Return to main menu",
    "exit": "Quit the tool",
}

SUBSET_HINTS_3D = {
    "volume": "Extract a subvolume (limits for each dimension)",
    "slice": "Slice along one coordinate (e.g., depth)",
    "xsection": "Interpolated vertical cross-section along a path",
    "surface": "Surface plot of a 2D variable",
    "back": "Return to main menu",
    "exit": "Quit the tool",
}

SLICE_ACTION_HINTS = {
    "plot2d": "2D image plot of the slice",
    "plot3d": "Interactive 3D surface of the slice",
    "gmap": "Geographic 2D plot (when applicable)",
    "cmap": "Change colormap; supports 'name' or 'name,vmin,vmax'",
    "save": "Save sliced data (csv/geocsv/netcdf)",
    "back": "Back to slice setup",
    "exit": "Quit the tool",
}


# -----------------------
# Dataclasses / Config
# -----------------------
@dataclass
class AppConfig:
    verbose: bool
    input_file: str
    cmap: str = prop.cmap
    interpolation_methods: List[str] = tuple(slicer_prop.interpolation_method)
    xsection_steps: int = slicer_prop.steps
    vertical_exaggeration: float = 0.0


@dataclass
class DsContext:
    ds: xr.Dataset
    dataset_type: str
    base_title: str
    coords: List[str]
    data_vars: List[str]
    coord_values: Dict[str, List[float]]


class _Once:
    printed = set()


def _subset_with_xarray(
    ds: xr.Dataset, limits: Dict[str, Tuple[float, float]]
) -> xr.Dataset:
    """
    Robust subsetting:
      1) Try .sel with slice(...) for each 1D coord (handles asc/desc).
      2) If .sel fails (e.g., non-monotonic coords), fallback to .where(..., drop=True).
    """
    # First, try building a selectors dict for .sel
    selectors = {}
    for dim, (lo, hi) in limits.items():
        if dim not in ds.coords:
            continue
        coord = ds[dim]
        if coord.ndim != 1:
            # skip non-1D coordinates in selectors; we'll handle them in the where-fallback
            continue
        # If coordinate runs descending, flip slice bounds so xarray can walk along it
        first, last = float(coord.values[0]), float(coord.values[-1])
        if first <= last:
            selectors[dim] = slice(lo, hi)
        else:
            selectors[dim] = slice(hi, lo)

    try:
        sub = ds.sel(selectors)
        # Also filter any leftover dims (non-1D or dims we didn’t sel on) via where
        for dim, (lo, hi) in limits.items():
            if dim in selectors:
                continue
            if dim in sub.coords and np.issubdtype(sub[dim].dtype, np.number):
                lo2, hi2 = (min(lo, hi), max(lo, hi))
                sub = sub.where((sub[dim] >= lo2) & (sub[dim] <= hi2), drop=True)
        return sub
    except Exception:
        # Fallback: numeric where on all provided limits
        sub = ds
        for dim, (lo, hi) in limits.items():
            if dim in sub.coords and np.issubdtype(sub[dim].dtype, np.number):
                lo2, hi2 = (min(lo, hi), max(lo, hi))
                sub = sub.where((sub[dim] >= lo2) & (sub[dim] <= hi2), drop=True)
        return sub


def _normalize_limits(lims: Mapping) -> Dict[str, Tuple[float, float]]:
    norm: Dict[str, Tuple[float, float]] = {}
    for k, v in lims.items():
        key = str(k)
        if isinstance(v, (list, tuple, np.ndarray)):
            if len(v) != 2:
                raise ValueError(f"Bad limit for {key}: {v}")
            lo, hi = float(v[0]), float(v[1])
        else:
            # Catch weird cases (e.g., a scalar slipped in)
            raise ValueError(f"Bad limit type for {key}: {type(v).__name__}")
        # enforce (min, max) and tuples
        norm[key] = (min(lo, hi), max(lo, hi))
    return norm


def print_once(key: str, fn):
    if key not in _Once.printed:
        fn()
        _Once.printed.add(key)


def _parse_pair(raw: str, what: str) -> Optional[Tuple[float, float]]:
    try:
        a, b = (float(x.strip()) for x in raw.split(",", 1))
        return a, b
    except Exception:
        err(f"Invalid {what} '{raw}'. Expected 'a,b'.")
        return None


def _norm_lon(lon: float) -> float:
    """Normalize to [-180, 180]."""
    lon = float(lon)
    while lon > 180:
        lon -= 360
    while lon <= -180:
        lon += 360
    return lon


def get_spatial_bounds(ctx: DsContext) -> Optional[Dict[str, Tuple[float, float]]]:
    """
    Return {'latitude': (lat_min, lat_max), 'longitude': (lon_min, lon_max)}.
    Use CF global attrs when present; fall back to coords.
    """
    ds = ctx.ds
    if "latitude" not in ds.coords or "longitude" not in ds.coords:
        return None

    lat_min = ds.attrs.get("geospatial_lat_min")
    lat_max = ds.attrs.get("geospatial_lat_max")
    lon_min = ds.attrs.get("geospatial_lon_min")
    lon_max = ds.attrs.get("geospatial_lon_max")

    if (
        lat_min is not None
        and lat_max is not None
        and lon_min is not None
        and lon_max is not None
    ):
        lat_min, lat_max = float(lat_min), float(lat_max)
        lon_min, lon_max = _norm_lon(float(lon_min)), _norm_lon(float(lon_max))
    else:
        # fall back to coord arrays
        lat = np.asarray(ds["latitude"].values).astype(float)
        lon = np.asarray(ds["longitude"].values).astype(float)
        lat_min, lat_max = float(np.nanmin(lat)), float(np.nanmax(lat))
        # normalize and then bound
        lon = np.vectorize(_norm_lon)(lon)
        lon_min, lon_max = float(np.nanmin(lon)), float(np.nanmax(lon))

    # ensure (min <= max)
    if lat_min > lat_max:
        lat_min, lat_max = lat_max, lat_min
    if lon_min > lon_max:
        lon_min, lon_max = lon_max, lon_min

    return {"latitude": (lat_min, lat_max), "longitude": (lon_min, lon_max)}


def _set_axes_from_data(ax: plt.Axes, da: xr.DataArray) -> None:
    """
    Set x/y limits using the actual coordinate values on the plotted DataArray.
    Works for 2D plots where the last two dims are the plotted axes.
    """
    if da.ndim < 2:
        return
    ydim, xdim = da.dims[-2], da.dims[-1]

    # Coordinate values -> numeric arrays
    x = np.asarray(da[xdim].values, dtype=float)
    y = np.asarray(da[ydim].values, dtype=float)

    # handle degenerate or NaN-only cases defensively
    if x.size and np.isfinite(x).any():
        ax.set_xlim(np.nanmin(x), np.nanmax(x))
    if y.size and np.isfinite(y).any():
        ax.set_ylim(np.nanmin(y), np.nanmax(y))


def _ask_point(tag: str) -> Optional[Tuple[float, float]]:
    raw = input(f"[xsection] {tag} point lon,lat (back/exit, h)? ").strip()
    if raw.lower() in ("h", "?"):
        tip("Enter longitude,latitude. Example: -120.5, 35.2")
        return _ask_point(tag)
    if raw == "exit":
        sys.exit(0)
    if raw == "back" or not raw:
        return None
    pair = _parse_pair(raw, f"{tag} point")
    return pair


def _append_history(ds: xr.Dataset, note: str) -> None:
    """Append a UTC-stamped note to global history."""
    from datetime import datetime, timezone

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    prev = ds.attrs.get("history", "")
    ds.attrs["history"] = (prev + ("; " if prev else "") + f"{ts} {note}").strip()


def update_geospatial_attrs_inplace(ds: xr.Dataset) -> None:
    """
    Update CF-style geospatial global attributes based on the dataset's current coords.
    Sets: geospatial_lat/lon_(min|max|units), geospatial_vertical_(min|max|units|positive).
    No-op for missing coords.
    """

    def set_if(k, v):
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            ds.attrs[k] = v

    # Latitude
    if "latitude" in ds.coords:
        lat = np.asarray(ds["latitude"].values, dtype=float)
        set_if("geospatial_lat_min", float(np.nanmin(lat)))
        set_if("geospatial_lat_max", float(np.nanmax(lat)))
        set_if("geospatial_lat_units", ds["latitude"].attrs.get("units", "degrees"))

    # Longitude (keep native range; don’t force wrap unless you want to)
    if "longitude" in ds.coords:
        lon = np.asarray(ds["longitude"].values, dtype=float)
        set_if("geospatial_lon_min", float(np.nanmin(lon)))
        set_if("geospatial_lon_max", float(np.nanmax(lon)))
        set_if("geospatial_lon_units", ds["longitude"].attrs.get("units", "degrees"))

    # Vertical (prefer common names)
    zvar = None
    for cand in ("depth", "height", "z"):
        if cand in ds.coords:
            zvar = cand
            break
    if zvar:
        z = np.asarray(ds[zvar].values, dtype=float)
        set_if("geospatial_vertical_min", float(np.nanmin(z)))
        set_if("geospatial_vertical_max", float(np.nanmax(z)))
        set_if("geospatial_vertical_units", ds[zvar].attrs.get("units", ""))
        pos = ds[zvar].attrs.get("positive")
        if pos:
            set_if("geospatial_vertical_positive", pos)


def prompt_cmap(
    cfg: AppConfig,
    prompt_prefix: str = "[cmap]",
    allow_back: bool = True,
    allow_exit: bool = True,
) -> Optional[str]:
    """
    Interactively choose a colormap.
    - 'list' shows available colormaps and re-prompts.
    - '' (blank) keeps current cmap and returns None.
    - 'back' returns None (no change).
    - 'exit' exits the app.
    - 'name' or 'name,vmin,vmax' sets cfg.cmap (vmin/vmax message only).
    Returns the chosen cmap name, or None if unchanged.
    """
    while True:
        raw = input(
            f"{prompt_prefix} name or 'name,vmin,vmax' "
            f"(type 'list' to see available cmaps; blank=keep '{cfg.cmap}', back/exit)? "
        ).strip()

        low = raw.lower()
        if allow_exit and low == "exit":
            sys.exit(0)
        if allow_back and (low == "back" or raw == ""):
            # back or blank → no change
            return None

        if low == "list":
            cmaps = plt.colormaps()
            logger.info(f"Available colormaps ({len(cmaps)}):")
            cols = 4
            for i in range(0, len(cmaps), cols):
                logger.info("  ".join(f"{c:<20}" for c in cmaps[i : i + cols]))
            # stay in this loop and re-prompt
            continue

        parts = [p.strip() for p in raw.split(",")]
        if len(parts) == 1:
            name = parts[0]
            if name in plt.colormaps():
                cfg.cmap = name
                note(f"cmap set to {cfg.cmap}")
                return name
            else:
                err(f"invalid cmap: {name}")
                # re-prompt
                continue

        if len(parts) == 3:
            name, vmin, vmax = parts
            if name in plt.colormaps():
                cfg.cmap = name
                note(
                    f"cmap={name}, vmin={vmin}, vmax={vmax} "
                    "(note: per-plot vmin/vmax wiring can be added)"
                )
                return name
            else:
                err(f"invalid cmap: {name}")
                continue

        err("invalid cmap spec")  # re-prompt


# -----------------------
# CLI / Usage
# -----------------------
def usage() -> None:
    logger.info(
        f"""
A tool to extract data from an EMC netCDF file. Inspect metadata, slice, plot, and save.
Slicing supports existing axes and cross-sections for gridded 3D data.

Options:
  -h, --help        Show this help
  -v, --verbose     Verbose output
  -m, --meta        Print metadata and exit
  -i, --input FILE  Input netCDF file (required)
        """
    )


def parse_args(argv: List[str]) -> Tuple[AppConfig, bool]:
    try:
        opts, _ = getopt.getopt(argv, "hvmi:", ["help", "verbose", "meta", "input="])
    except getopt.GetoptError as err:
        logger.error(err)
        usage()
        sys.exit(2)

    verbose = False
    meta_only = False
    input_file = None

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o in ("-v", "--verbose"):
            verbose = True
        elif o in ("-m", "--meta"):
            meta_only = True
        elif o in ("-i", "--input"):
            input_file = a

    if not input_file or not os.path.isfile(input_file):
        logger.error("[ERR] missing or invalid -i/--input")
        usage()
        sys.exit(1)

    cfg = AppConfig(verbose=verbose, input_file=input_file)
    return cfg, meta_only


# -----------------------
# Dataset utilities
# -----------------------
def get_netcdf_engine(filename: str) -> Optional[str]:
    for ext, engine in prop.netcdf_engines.items():
        if filename.endswith(ext):
            return engine
    return None


def compute_dataset_type(ds: xr.Dataset) -> str:
    has2d = any(len(da.dims) == 2 for da in ds.data_vars.values())
    has3d = any(len(da.dims) == 3 for da in ds.data_vars.values())
    if has3d:
        return "3D"
    if has2d:
        return "2D"
    return "Unknown"


def load_dataset(cfg: AppConfig) -> DsContext:
    engine = get_netcdf_engine(cfg.input_file)

    with timed(f"Loading dataset '{os.path.basename(cfg.input_file)}'", cfg.verbose):
        with xr.open_dataset(cfg.input_file, engine=engine) as ds:
            ds = ds.load()  # load into memory for interactive use

    base_title = ds.attrs.get("model") or ds.attrs.get("id", "")
    dataset_type = compute_dataset_type(ds)

    coords = list(ds.coords)
    data_vars = list(ds.data_vars)
    coord_values: Dict[str, List[float]] = {}

    for var in coords:
        if var in ds and ds[var].size > 0:
            coord_values[var] = list(np.asarray(ds[var].data).ravel())
        else:
            warn(f"Missing or empty coord '{var}'")

    print("\n\n\n")
    note(f"Loaded {cfg.input_file}")
    kv("engine", engine or "<auto>")
    kv("dataset type", dataset_type)
    kv("coords", ", ".join(coords))
    kv("data vars", ", ".join(data_vars) if data_vars else "<none>")
    if base_title:
        kv("title", base_title)

    return DsContext(
        ds=ds,
        dataset_type=dataset_type,
        base_title=base_title,
        coords=coords,
        data_vars=data_vars,
        coord_values=coord_values,
    )


# -----------------------
# Presentation helpers
# -----------------------
def _to_python_type(val):
    if isinstance(val, np.generic):
        return val.item()
    if isinstance(val, np.ndarray):
        return val.tolist()
    return val


def show_metadata(ctx: DsContext) -> None:
    header(f"Metadata ({ctx.dataset_type})")
    subheader("Global attributes (geospatial*)")
    for k, v in ctx.ds.attrs.items():
        if "geospatial" in k:
            kv(k, _to_python_type(v))

    subheader("Coordinate Variables")
    indent = 14
    for name in ctx.ds.dims:
        logger.info(f"\n  {name}:")
        for ak, av in ctx.ds[name].attrs.items():
            logger.info(f"{' '*indent}{ak}: {av}")
        vals = np.asarray(ctx.ds[name].values).astype(float).ravel().tolist()
        val_str = ", ".join(f"{v:g}" for v in vals)
        logger.info(f"{' '*indent}Values:")
        logger.info(
            f"{' '*indent}{fill(val_str, width=88, subsequent_indent=' '*indent)}"
        )

    subheader("Data Variables")
    slicer_lib.display_var_meta(ctx.ds, ctx.data_vars, indent=14, values=False)


def prompt_for_variable(
    valid_vars: List[str],
    prompt_msg: str = "variable",
    allow_back: bool = True,
    allow_exit: bool = True,
    case_insensitive: bool = True,
) -> Optional[str]:
    """
    Keep prompting until a valid variable is entered.
    Returns the canonical variable name, or None if user types 'back' (when allowed).
    Calls sys.exit(0) on 'exit' (when allowed).
    """
    if not valid_vars:
        warn("No variables available.")
        return None

    # Map lowercase -> canonical (preserve original capitalization)
    canon = {v.lower(): v for v in valid_vars} if case_insensitive else None

    while True:
        choice = input(f"{prompt_msg} {valid_vars} (back/exit)? ").strip()
        low = choice.lower()

        if allow_exit and low == "exit":
            sys.exit(0)
        if allow_back and (low == "back" or choice == ""):
            return None

        # validate
        if case_insensitive:
            if low in canon:
                return canon[low]
        else:
            if choice in valid_vars:
                return choice

        warn(f"Invalid variable '{choice}'. Please choose one from the list.")


def show_ranges(ctx: DsContext) -> None:
    header(f"Ranges ({ctx.dataset_type})")

    subheader("Coordinate Variables")
    for var in ctx.ds.dims:
        raw_units = ctx.ds[var].attrs.get("units", None)
        units = _format_units(raw_units)
        if units is None:  # only truly missing (no 'units' attr) is an error
            warn(f"Missing 'units' attribute for '{var}' (treating as unknown)")
            units = "(unknown)"
        vmin, vmax = np.nanmin(ctx.ds[var].data), np.nanmax(ctx.ds[var].data)
        kv(var, f"{vmin:0.2f} to {vmax:0.2f} {units}")

    subheader("Data Variables")
    for var in ctx.data_vars:
        raw_units = ctx.ds[var].attrs.get("units", None)
        units = _format_units(raw_units)
        if units is None:  # only truly missing is an issue; warn, don't kill the app
            warn(f"Missing 'units' attribute for '{var}' (treating as unknown)")
            units = "(unknown)"
        vmin, vmax = np.nanmin(ctx.ds[var].data), np.nanmax(ctx.ds[var].data)
        kv(var, f"{vmin:0.2f} to {vmax:0.2f} {units}")


# -----------------------
# Output / saving
# -----------------------
def save_dataset(subset: xr.Dataset, filename: str, messages: List[str]) -> List[str]:
    precision = 3
    data_vars = list(subset.data_vars)
    main_coords = list(subset.dims)

    # ---- UPDATE GEO ATTRS based on the subset ----
    try:
        update_geospatial_attrs_inplace(subset)
        # Optional: record the subset ranges in history
        pretty = ", ".join(
            f"{d}=[{float(subset[d].min()):.3g},{float(subset[d].max()):.3g}]"
            for d in main_coords
            if d in subset.coords
        )
        _append_history(subset, f"subset saved with {pretty}")
    except Exception as ex:
        logger.warning(f"[WARN] could not update geospatial attrs: {ex}")

    # round aux coords
    for c in subset.coords:
        if c not in main_coords:
            subset[c].values = np.round(subset[c].values, precision)

    if filename.endswith(".gcsv"):
        meta = lib.get_geocsv_metadata_from_ds(subset)
        meta = f"# dataset: GeoCSV2.0\n# delimiter: ,\n{meta}"
        with open(filename, "w") as fp:
            fp.write(meta)
        subset.to_dataframe().to_csv(filename, mode="a")
    elif filename.endswith(".csv"):
        subset.to_dataframe().to_csv(filename, mode="w")
    elif filename.endswith(".nc"):
        enc = {v: {"zlib": True, "complevel": 4} for v in data_vars}
        for c in main_coords:
            enc.setdefault(c, {})["_FillValue"] = None
        subset.to_netcdf(filename, mode="w", format="NETCDF4_CLASSIC", encoding=enc)
    elif filename.endswith(".h5"):
        enc = {v: {"zlib": True, "complevel": 4} for v in data_vars}
        for c in main_coords:
            enc.setdefault(c, {})["_FillValue"] = None
        subset.to_netcdf(
            filename, mode="w", engine=prop.netcdf_engines[".h5"], encoding=enc
        )
    else:
        messages.append(f"[ERR] invalid file: {filename}")
        return messages

    messages.append(f"[INFO] Saved as {filename}")
    return messages


# -----------------------
# Interactive flows
# -----------------------
def prompt_save(subset: xr.Dataset) -> None:
    """
    Keep prompting the user to save the given Dataset until they either save,
    type back, or exit. Supports 'h' for help.
    """
    header("Save")
    tip("Format: <filename>,<type>  (types: csv|geocsv|netcdf)")
    tip("Examples:  out,netcdf   → writes out.nc (NETCDF4_CLASSIC)")
    tip("           out,geocsv   → writes out.csv with GeoCSV header")
    tip("           out,csv      → writes out.csv (no GeoCSV header)")

    messages: List[str] = []

    while True:
        # flush, then prompt
        try:
            sys.stdout.flush()
        except Exception:
            pass

        raw = input("[save] filename,type (e.g., out,netcdf | back/exit, h)? ").strip()

        low = raw.lower()
        if low == "exit":
            sys.exit(0)
        if low == "back" or raw == "":
            # blank means “don’t save now”; return to the action menu
            return
        if low == "h" or low == "?":
            logger.info("Help:")
            logger.info("  out,netcdf  → writes out.nc (NETCDF4_CLASSIC)")
            logger.info("  out,geocsv  → writes out.csv with GeoCSV header")
            logger.info("  out,csv     → writes out.csv (no GeoCSV header)")
            continue

        if "," not in raw:
            err(f"Bad input: '{raw}' (expected 'name,type')")
            continue

        name, kind = (p.strip() for p in raw.split(",", 1))
        if kind not in VALID_FILE_TYPES:
            err(f"Bad type: '{kind}' (valid: {list(VALID_FILE_TYPES)})")
            continue

        filename = f"{name}{VALID_FILE_TYPES[kind]}"
        note(f"Saving → {filename}")
        with timed("Writing file", True):
            messages = save_dataset(subset, filename, messages)

        for m in messages:
            # save_dataset already prefixes with [INFO]; keep output tidy
            note(m.replace("[INFO] ", ""))
        messages.clear()

        ok("Save complete")
        # after a successful save, either let them save another,
        # or type back/Enter to return to the action menu


def prompt_menu(verbose: bool, ctx: Optional[DsContext] = None) -> str:
    # compute which options to show
    has_geo = bool(ctx and "latitude" in ctx.ds.coords and "longitude" in ctx.ds.coords)
    options = (
        ["meta", "range", "subset"] + (["map"] if has_geo else []) + ["help", "exit"]
    )

    if verbose:
        header("Main Menu")
        tip("Type 'h' for hints at any prompt.")
        items = [
            "meta   — view file metadata",
            "range  — show variable ranges",
            "subset — volume / slice / xsection / surface",
        ]
        if has_geo:
            items.append("map    — show coverage map")
        items += ["help   — usage", "exit   — quit"]
        list_menu("Choose an option", items)

    choice = choose_with_hints(
        f"[data] select option [{', '.join(options)} | h]? ",
        valid=options,
        hints={k: v for k, v in MAIN_HINTS.items() if k in options},
        allow_blank=False,
    )
    return choice


def handle_subset(cfg: AppConfig, ctx: DsContext) -> None:
    while True:
        if ctx.dataset_type == "2D":
            valid = ["volume", "surface", "back", "exit"]
            if cfg.verbose:
                header("Subset Menu (2D)")
                tip("Type 'h' for hints at any prompt.")

            choice = choose_with_hints(
                "[subset] select [volume, surface, back, exit | h]? ",
                valid=valid,
                hints=SUBSET_HINTS_2D,
            )
        else:
            if (
                ctx.dataset_type != "2D"
                and "longitude" in ctx.ds.coords
                and "latitude" in ctx.ds.coords
            ):
                valid = ["volume", "slice", "xsection", "back", "exit"]
            else:
                valid = ["volume", "slice", "back", "exit"]
            if cfg.verbose:
                header("Subset Menu (2D/3D)")
                tip("Type 'h' for hints at any prompt.")
            choice = choose_with_hints(
                "[subset] select [volume, slice, xsection, back, exit | h]? ",
                valid=valid,
                hints=SUBSET_HINTS_3D,
            )

        if choice in ("back", ""):
            return
        if choice == "exit":
            sys.exit(0)

        if choice == "volume":
            do_volume_flow(cfg, ctx)
        elif choice == "slice":
            do_slice_flow(plot_var, cfg, ctx)
        elif choice == "xsection":
            do_xsection_flow(cfg, ctx)
        elif choice == "surface":
            do_surface_flow(cfg, ctx)


# ---- Volume flow (trimmed example) ----
def do_volume_flow(cfg: AppConfig, ctx: DsContext) -> None:
    ds = ctx.ds
    if cfg.verbose:
        header("Volume")
        list_menu(
            "Steps",
            [
                "1) Enter min,max for each dimension (blank = full range)",
                "2) Preview/Save",
            ],
        )

    # start with tuple defaults
    subset_limits: Dict[str, Tuple[float, float]] = {
        dim: (float(np.nanmin(ds[dim].data)), float(np.nanmax(ds[dim].data)))
        for dim in ds.dims
    }

    for dim, default in subset_limits.items():
        tip(f"{dim} valid range = {fmt_range(*default)}")
        raw = input(f"[subset-volume] {dim} range (min,max | back/exit, h)? ").strip()
        if raw.lower() in ("h", "?"):
            logger.info("Enter two numbers as 'min,max'. Blank = full range.")
            raw = ""  # fall back to full
        if raw == "exit":
            sys.exit(0)
        if raw == "back":
            return

        if raw:
            try:
                lo, hi = map(float, raw.split(",", 1))
                subset_limits[dim] = (min(lo, hi), max(lo, hi))  # tuple
                ok(f"{dim} = {fmt_range(*subset_limits[dim])}")
            except Exception:
                warn(f"Invalid '{raw}', using full {fmt_range(*default)}")
                subset_limits[dim] = default  # ensure tuple
        else:
            # user pressed Enter → explicitly reassign the default (tuple)
            subset_limits[dim] = default
            note(f"{dim} = full {fmt_range(*default)}")

    # normalize before calling subsetter (defensive)
    try:
        subset_limits = _normalize_limits(subset_limits)
    except Exception as ex:
        err(f"Invalid limits before subsetting: {ex}")
        return

    with timed("Subsetting volume", cfg.verbose):
        # Defensive normalize (if you kept _normalize_limits)
        # subset_limits = _normalize_limits(subset_limits)
        subset_volume = _subset_with_xarray(ds, subset_limits)

    # No 'warnings' object now; we produced a dataset (possibly empty)
    if any(size == 0 for size in subset_volume.dims.values()):
        warn("Subset produced an empty selection. Check your ranges.")
        return

    ok("Volume ready")
    for d, (lo, hi) in subset_limits.items():
        kv(d, fmt_range(lo, hi))

    prompt_save(subset_volume)


def show_coverage_map(cfg: AppConfig, ctx: DsContext) -> None:
    """
    Show a quick coverage map bounded by spatial lat/lon limits.
    Tries Cartopy; falls back to a simple Matplotlib rectangle.
    """
    bounds = get_spatial_bounds(ctx)
    if not bounds:
        err("Coverage map requires latitude/longitude coordinates.")
        return

    (lat_min, lat_max) = bounds["latitude"]
    (lon_min, lon_max) = bounds["longitude"]

    title = ctx.base_title or "Coverage"
    subtitle = f"Lat [{lat_min:.2f}, {lat_max:.2f}], Lon [{lon_min:.2f}, {lon_max:.2f}]"

    # Try Cartopy first
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        proj = ccrs.PlateCarree()
        with timed("Drawing coverage map (Cartopy)", cfg.verbose):
            fig = plt.figure(figsize=(7.5, 4.5))
            ax = plt.axes(projection=proj)
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=proj)
            ax.coastlines(resolution="110m", linewidth=0.6)
            ax.add_feature(cfeature.BORDERS.with_scale("110m"), linewidth=0.4)
            ax.add_feature(cfeature.LAND, facecolor="0.95")
            ax.add_feature(cfeature.OCEAN, facecolor="0.98")
            # draw the bounding box explicitly
            xs = [lon_min, lon_max, lon_max, lon_min, lon_min]
            ys = [lat_min, lat_min, lat_max, lat_max, lat_min]
            ax.plot(xs, ys, transform=proj)
            ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5, linestyle="--")
            plt.title(f"{title}\nCoverage area — {subtitle}")
            plt.tight_layout()
            plt.show()
        ok("Coverage map rendered")

    except Exception:
        # Fallback: plain Matplotlib
        with timed("Drawing coverage map (Matplotlib fallback)", cfg.verbose):
            fig, ax = plt.subplots(figsize=(7.0, 4.0))
            # bounding box
            xs = [lon_min, lon_max, lon_max, lon_min, lon_min]
            ys = [lat_min, lat_min, lat_max, lat_max, lat_min]
            ax.plot(xs, ys)
            ax.set_xlim(lon_min - 0.5, lon_max + 0.5)
            ax.set_ylim(lat_min - 0.5, lat_max + 0.5)
            ax.set_xlabel("longitude")
            ax.set_ylabel("latitude")
            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, linestyle="--", linewidth=0.3, alpha=0.6)
            plt.title(f"{title}\nCoverage area — {subtitle}")
            plt.tight_layout()
            plt.show()
        ok("Coverage map (fallback) rendered")


# ---- Slice/xsection/surface flows (stubs showing separation) ----
def do_slice_flow(plot_var, cfg: AppConfig, ctx: DsContext) -> None:
    ds = ctx.ds
    coords_main = [c for c in ds.coords if c in ds.dims]
    aux_coords = [c for c in ds.coords if c not in ds.dims]
    coord_values = ctx.coord_values
    data_vars = [v for v in ctx.data_vars if v not in aux_coords]

    while True:
        if cfg.verbose:
            print_once(
                "slice-header",
                lambda: (
                    header("Slice"),
                    list_menu(
                        "Steps",
                        [
                            "1) Choose direction (coord to slice along)",
                            "2) Enter target value (closest will be used)",
                            "3) Enter limits for remaining dimensions",
                            "4) Plot / save",
                        ],
                    ),
                    tip("Type 'h' for hints at any prompt."),
                ),
            )

        valid_dirs = coords_main + ["back", "exit"]
        direction = input(f"[slice] direction {coords_main} (back/exit)? ").strip()
        if direction == "exit":
            sys.exit(0)
        if direction == "back" or direction not in coords_main:
            return

        vmin = float(np.nanmin(ds[direction].data))
        vmax = float(np.nanmax(ds[direction].data))
        tip(f"{direction} valid range = {fmt_range(vmin, vmax)}")

        slice_value = None
        while slice_value is None:
            raw = input(f"[slice-{direction}] value (back/exit)? ").strip()
            if raw == "exit":
                sys.exit(0)
            if raw in ("", "back"):
                break
            try:
                slice_value = float(raw)
            except Exception:
                err(f"Invalid number: '{raw}'")
                slice_value = None

        if slice_value is None or raw == "back":
            continue

        closest = lib.closest(coord_values[direction], slice_value)
        ok(f"Using closest {direction} = {closest:g}")
        if abs(closest - slice_value) > 0:
            note(f"Requested {slice_value:g} differs by {abs(closest - slice_value):g}")

        slice_dims = [d for d in coords_main if d != direction]
        limits: Dict[str, Tuple[float, float]] = {}
        for dim in slice_dims:
            dmin = float(np.nanmin(ds[dim].data))
            dmax = float(np.nanmax(ds[dim].data))
            tip(
                f"{dim} valid range = {fmt_range(dmin, dmax)} (Enter 'min,max' or blank for full)"
            )
            raw = input(f"[slice-{direction}] {dim} limits (back/exit)? ").strip()
            if raw == "exit":
                sys.exit(0)
            if raw == "back":
                limits = None
                break
            if not raw:
                limits[dim] = (dmin, dmax)
                note(f"{dim} = full {fmt_range(dmin, dmax)}")
                continue
            try:
                lo, hi = map(float, raw.split(","))
                lo, hi = min(lo, hi), max(lo, hi)
                limits[dim] = (lo, hi)
                ok(f"{dim} = {fmt_range(lo, hi)}")
            except Exception:
                warn(f"Invalid '{raw}', using full range")
                limits[dim] = (dmin, dmax)

        if limits is None:
            continue

        if limits is None:
            continue  # back to direction selection

        # perform slice
        with timed("Computing slice", cfg.verbose):
            sliced = slicer_lib.slicer(ds, direction, closest, limits)

        ok("Slice ready")
        kv("direction", direction)
        for d, (lo, hi) in limits.items():
            kv(d, fmt_range(lo, hi))
        kv("variables", ", ".join(data_vars) if data_vars else "<none>")

        logger.info(f"\nSlice summary:\n{sliced}")
        # actions: plot2d/plot3d/gmap/cmap/save/back
        while True:
            gmap_option = ""
            if "latitude" in ds.coords and "longitude" in ds.coords:
                # only if direction not in lat/lon dims:
                if (direction not in ds["latitude"].dims) and (
                    direction not in ds["longitude"].dims
                ):
                    gmap_option = ", gmap"

            valid_actions = ["plot2d", "plot3d", "cmap", "save", "back", "exit"]
            if gmap_option:
                valid_actions.insert(2, "gmap")  # after plot3d

            action = choose_with_hints(
                f"[slice-{direction}] Action [{', '.join(valid_actions)} | h]? ",
                valid=valid_actions,
                hints=SLICE_ACTION_HINTS,
            )

            if action == "exit":
                sys.exit(0)
            if action in ("back", ""):
                break

            # inside the 'plot2d' branch
            if action == "plot2d":
                repeat = True
                while repeat:
                    if not data_vars:
                        err("No plottable variables.")
                        break
                    if len(data_vars) == 1:
                        repeat = False
                    plot_var = (
                        data_vars[0]
                        if len(data_vars) == 1
                        else prompt_for_variable(
                            data_vars, f"[slice-{direction}-plot2d] variable"
                        )
                    )
                    if plot_var is None:  # user typed back or blank
                        break  # return to action menu

                    plot_data = sliced[plot_var].copy()
                    if "depth" in plot_data.dims and (
                        "geospatial_vertical_positive" in sliced.attrs
                    ):
                        if sliced.attrs["geospatial_vertical_positive"] == "down":
                            plot_data["depth"] = -plot_data["depth"]
                    plot_data.plot(cmap=cfg.cmap)
                    title = f"{ctx.base_title}\n{direction} slice of {plot_var} at {closest} {ds[direction].attrs.get('units','')}"
                    plt.title(title)
                    _set_axes_from_data(plt.gca(), plot_data)

                    plt.tight_layout()
                    plt.show()

            # inside the 'plot3d' branch
            elif action == "plot3d":
                repeat = True
                while repeat:
                    if not data_vars:
                        err("No plottable variables.")
                        break
                    if len(data_vars) == 1:
                        repeat = False
                    plot_var = (
                        data_vars[0]
                        if len(data_vars) == 1
                        else prompt_for_variable(
                            data_vars, f"[slice-{direction}-plot3d] variable"
                        )
                    )
                    if plot_var is None:
                        break
                    sliced[plot_var].plot.surface(cmap=cfg.cmap)
                    title = f"{ctx.base_title}\n{direction} slice (surface) of {plot_var} at {closest} {ds[direction].attrs.get('units','')}"
                    plt.title(title)
                    plt.tight_layout()
                    plt.show()

            # inside the 'gmap' branch
            elif action == "gmap" and gmap_option:
                repeat = True
                while repeat:
                    if not data_vars:
                        err("No plottable variables.")
                        break

                    if len(data_vars) == 1:
                        repeat = False
                    plot_var = (
                        data_vars[0]
                        if len(data_vars) == 1
                        else prompt_for_variable(
                            data_vars, f"[slice-{direction}-gmap] variable"
                        )
                    )
                    if plot_var is None:
                        break

                    gmap_limits = {}
                    gmap_limits["longitude"] = limits.get(
                        "longitude",
                        (
                            float(np.nanmin(ds["longitude"])),
                            float(np.nanmax(ds["longitude"])),
                        ),
                    )
                    gmap_limits["latitude"] = limits.get(
                        "latitude",
                        (
                            float(np.nanmin(ds["latitude"])),
                            float(np.nanmax(ds["latitude"])),
                        ),
                    )
                    title = f"{ctx.base_title}\n{direction} slice of {plot_var} at {closest} {ds[direction].attrs.get('units','')}"
                    slicer_lib.gmap(
                        plot_var,
                        cfg.cmap,
                        gmap_limits,
                        sliced,
                        vmin=None,
                        vmax=None,
                        title=title,
                    )

            elif action == "cmap":
                prompt_cmap(cfg, prompt_prefix=f"[slice-{direction}-cmap]")
            elif action == "save":
                if plot_var is None:
                    plot_var = (
                        data_vars[0]
                        if len(data_vars) == 1
                        else prompt_for_variable(data_vars, f"[save] variable")
                    )
                prompt_save(sliced[[plot_var]])  # only that var
                continue


def do_xsection_flow(cfg: AppConfig, ctx: DsContext) -> None:
    ds = ctx.ds

    # Require geographic coords
    if "longitude" not in ds.coords or "latitude" not in ds.coords:
        err("xsection requires geographic coordinates: longitude, latitude")
        return

    # Pick Z coordinate
    zvar = "depth" if "depth" in ds.coords else None
    if not zvar:
        # choose first 3rd dimension-like coord that isn't lon/lat
        candidates = [c for c in ds.coords if c not in ("longitude", "latitude")]
        if not candidates:
            err("No vertical coordinate available for xsection.")
            return
        zvar = candidates[0]
        note(f"Using '{zvar}' as vertical coordinate")

    # Units & metadata
    z_units = ds[zvar].attrs.get("units", "")
    z_label = zvar
    positive = ds[zvar].attrs.get("positive", "down" if zvar == "depth" else "up")
    model_attrs = ds.attrs
    grid_ref = model_attrs.get("grid_ref", "latitude_longitude")
    utm_zone = model_attrs.get("utm_zone", None)
    ellipsoid = model_attrs.get("ellipsoid", "WGS84")

    if cfg.verbose:
        header("Cross-section")
        list_menu(
            "Steps",
            [
                "1) Enter start lon,lat",
                "2) Enter end lon,lat",
                f"3) Enter {zvar} min,max (blank = full)",
                "4) Interpolation method & number of points",
                "5) Plot / save",
            ],
        )
        tip("Type 'h' for hints at any prompt.")

    # ---- Start / End points
    while True:
        start = _ask_point("start")
        if start is None:  # back
            return
        end = _ask_point("end")
        if end is None:
            return
        if start == end:
            warn("Start and end are identical. Choose different points.")
            continue
        break

    # ---- Depth (or vertical) range
    zmin = float(np.nanmin(ds[zvar]))
    zmax = float(np.nanmax(ds[zvar]))
    tip(f"{zvar} valid range = {fmt_range(zmin, zmax)} {z_units or ''}")
    raw = input(f"[xsection] {zvar} min,max (blank=full | back/exit, h)? ").strip()
    if raw.lower() in ("h", "?"):
        tip(f"Example: {zmin:.2f},{zmax:.2f}")
        raw = input(f"[xsection] {zvar} min,max (blank=full | back/exit)? ").strip()
    if raw == "exit":
        sys.exit(0)
    if raw == "back":
        return
    if raw:
        zlohi = _parse_pair(raw, f"{zvar} range")
        if not zlohi:
            return
        zlo, zhi = min(zlohi), max(zlohi)
    else:
        zlo, zhi = zmin, zmax
    ok(f"{zvar} = {fmt_range(zlo, zhi)} {z_units or ''}")

    # ---- Interpolation method & steps
    default_method = cfg.interpolation_methods[0]
    methods_list = ", ".join(cfg.interpolation_methods)
    raw = input(
        f"[xsection] interpolation method [{methods_list}] (default {default_method} | back/exit, h)? "
    ).strip()
    if raw.lower() in ("h", "?"):
        tip("Common choices: 'linear', 'nearest'.")
        raw = input(
            f"[xsection] method [{methods_list}] (default {default_method} | back/exit)? "
        ).strip()
    if raw == "exit":
        sys.exit(0)
    if raw == "back" or not raw:
        interp_method = default_method
    else:
        interp_method = raw if raw in cfg.interpolation_methods else default_method
        if interp_method == default_method and raw != default_method:
            warn(f"Unknown method '{raw}', using default '{default_method}'")

    raw = input(
        f"[xsection] number of points (default {cfg.xsection_steps} | back/exit, h)? "
    ).strip()
    if raw.lower() in ("h", "?"):
        tip("Points along the path, including endpoints. Higher = smoother line.")
        raw = input(
            f"[xsection] number of points (default {cfg.xsection_steps} | back/exit)? "
        ).strip()
    if raw == "exit":
        sys.exit(0)
    if raw == "back" or not raw:
        steps = cfg.xsection_steps
    else:
        steps = int(raw) if raw.isnumeric() else cfg.xsection_steps
        if not raw.isnumeric():
            warn(f"Invalid steps '{raw}', using default {steps}")

    # ---- Subset vertically first to speed up interpolation
    plot_data = ds.where((ds[zvar] >= zlo) & (ds[zvar] <= zhi), drop=True)

    # ---- Interpolate cross-section
    with timed("Interpolating cross-section", cfg.verbose):
        try:
            xsec, lats, lons = slicer_lib.interpolate_path(
                plot_data,
                start,
                end,
                num_points=steps,
                method=interp_method,
                grid_ref=grid_ref,
                utm_zone=utm_zone,
                ellipsoid=ellipsoid,
            )
        except Exception as ex:
            err(f"interpolate_path failed: {ex}")
            return

    # ---- Distance coordinate
    geod = Geod(ellps=ellipsoid)
    _, _, dists = geod.inv(lons[:-1], lats[:-1], lons[1:], lats[1:])
    cumdist = np.concatenate(([0.0], np.cumsum(dists)))  # meters
    # Convert to km if sensible
    if np.nanmax(cumdist) > 1000:
        cumdist = cumdist / 1000.0
        dist_label = "distance"
        dist_unit = "km"
    else:
        dist_label = "distance"
        dist_unit = "m"

    xsec = xsec.assign_coords(distance=("points", cumdist))
    xsec = xsec.swap_dims({"points": "distance"})

    # set attributes for the distance coordinate
    xsec["distance"].attrs["long_name"] = dist_label
    xsec["distance"].attrs["units"] = dist_unit

    # ---- Plot one variable quickly, or let user choose
    plot_vars = [v for v in ctx.data_vars if v not in ("longitude", "latitude", zvar)]
    if not plot_vars:
        err("No plottable data variables found.")
        return

    while True:
        plot_var = (
            plot_vars[0]
            if len(plot_vars) == 1
            else prompt_for_variable(plot_vars, "[xsection] variable to plot")
        )
        if plot_var is None:
            return  # back to previous menu

        pdata = xsec.copy()
        pdata[zvar].attrs["long_name"] = z_label
        pdata[zvar].attrs["units"] = z_units

        pdata = xsec.copy()
        # flip first
        if zvar in pdata.dims and zvar == "depth":
            pdata[zvar] = (-1 if positive == "down" else 1) * pdata[zvar]
        # then set attrs
        pdata[zvar].attrs.update({"long_name": z_label, "units": z_units})

        with timed(f"Plotting xsection of {plot_var}", cfg.verbose):
            pdata[plot_var].plot.contourf(cmap=cfg.cmap)
            plt.gca().set_xlabel(dist_label)
            plt.title(f"{ctx.base_title}\nCross-section of {plot_var}")
            _set_axes_from_data(plt.gca(), pdata[plot_var])

            plt.tight_layout()
            plt.show()
        # loop continues so user can plot another var or type 'back'

        # ---- Save option
        raw = input("[xsection] save data (y/n)? ").strip().lower()
        if raw == "y":
            _data = pdata.copy()
            # Undo the depth sign change, if needed.
            _data = pdata.copy()
            if zvar in _data.dims and zvar == "depth":
                _data[zvar] = (-1 if positive == "down" else 1) * _data[zvar]
            # restore attrs again
            _data[zvar].attrs.update(
                {"long_name": z_label, "units": z_units, "positive": positive}
            )
            prompt_save(_data[[plot_var]])  # only that var


def do_surface_flow(cfg: AppConfig, ctx: DsContext) -> None:
    """
    Surface plotting for 2D variables.
    - choose a 2D data variable
    - choose limits for its two dims (blank = full)
    - plot2d / gmap / cmap / save
    """
    ds = ctx.ds

    # Find 2D variables
    vars_2d = [v for v, da in ds.data_vars.items() if len(da.dims) == 2]
    if not vars_2d:
        warn("No 2D variables found. Try 'slice' or 'xsection' instead.")
        return

    if cfg.verbose:
        header("Surface")
        list_menu(
            "Steps",
            [
                "1) Choose a 2D variable",
                "2) Enter limits for its two dimensions (blank = full)",
                "3) Plot / gmap / cmap / save",
            ],
        )
        tip("Type 'h' for hints at any prompt.")

    # Choose variable
    if not vars_2d:
        warn("No 2D variables found. Try 'slice' or 'xsection' instead.")
        return

    var = (
        vars_2d[0]
        if len(vars_2d) == 1
        else prompt_for_variable(
            vars_2d,
            prompt_msg="[surface] variable",
            allow_back=True,
            allow_exit=True,
        )
    )
    if var is None:
        return  # user typed back

    da = ds[var]
    dims = list(da.dims)  # exactly 2
    d0, d1 = dims

    # Collect limits per dim
    limits: Dict[str, Tuple[float, float]] = {}
    for dim in dims:
        dmin = float(np.nanmin(ds[dim]))
        dmax = float(np.nanmax(ds[dim]))
        tip(
            f"{dim} valid range = {fmt_range(dmin, dmax)} (Enter 'min,max' or blank for full)"
        )
        raw = input(f"[surface-{var}] {dim} limits (back/exit, h)? ").strip()
        if raw.lower() in ("h", "?"):
            tip("Enter two numbers like: 10,20   (blank = full range)")
            raw = input(f"[surface-{var}] {dim} limits (back/exit)? ").strip()
        if raw == "exit":
            sys.exit(0)
        if raw == "back":
            return
        if not raw:
            limits[dim] = (dmin, dmax)
            note(f"{dim} = full {fmt_range(dmin, dmax)}")
            continue
        try:
            lo, hi = map(float, raw.split(",", 1))
            lo, hi = min(lo, hi), max(lo, hi)
            limits[dim] = (lo, hi)
            ok(f"{dim} = {fmt_range(lo, hi)}")
        except Exception:
            warn(f"Invalid '{raw}', using full range")
            limits[dim] = (dmin, dmax)

    # Subselect
    with timed(f"Subselecting {var}", cfg.verbose):
        sub = ds.sel(
            {
                d0: slice(limits[d0][0], limits[d0][1]),
                d1: slice(limits[d1][0], limits[d1][1]),
            }
        )

    # Determine whether gmap is available
    gmap_available = (
        "latitude" in ds.coords
        and "longitude" in ds.coords
        and d0 in ds["latitude"].dims
        and d1 in ds["longitude"].dims
    ) or (d0 in ("latitude", "longitude") and d1 in ("latitude", "longitude"))

    # Build action menu + hints dynamically
    valid_actions = ["plot2d", "cmap", "save", "back", "exit"]
    if gmap_available:
        valid_actions.insert(1, "gmap")

    # Filter hints to only what we show
    local_hints = {k: v for k, v in SLICE_ACTION_HINTS.items() if k in valid_actions}

    while True:
        action = choose_with_hints(
            f"[surface-{var}] Action [{', '.join(valid_actions)} | h]? ",
            valid=valid_actions,
            hints=local_hints,
        )

        if action == "exit":
            sys.exit(0)
        if action in ("back", ""):
            return

        # Loop through coordinates in sub. Check for positive direction.
        for coord_name, coord in sub.coords.items():
            if "positive" in coord.attrs and coord.attrs["positive"].lower() == "down":
                # Multiply by -1 and keep attributes
                flipped = -1 * coord
                flipped.attrs = dict(coord.attrs)  # preserve existing attrs
                flipped.attrs["positive"] = "up"  # update to reflect the change
                sub = sub.assign_coords({coord_name: flipped})

        if action == "plot2d":
            with timed(f"Plotting {var}", cfg.verbose):
                sub[var].plot(cmap=cfg.cmap)
                title = f"{ctx.base_title}\nSurface plot of {var}"
                plt.title(title)
                ax = plt.gca()
                _set_axes_from_data(ax, sub[var])
                plt.tight_layout()
                plt.show()
            ok("Plot rendered")

        elif action == "gmap" and gmap_available:
            # gmap limits based on chosen limits if applicable
            gmap_limits = {}
            if "longitude" in limits:
                gmap_limits["longitude"] = limits["longitude"]
            else:
                gmap_limits["longitude"] = (
                    float(np.nanmin(ds["longitude"])),
                    float(np.nanmax(ds["longitude"])),
                )
            if "latitude" in limits:
                gmap_limits["latitude"] = limits["latitude"]
            else:
                gmap_limits["latitude"] = (
                    float(np.nanmin(ds["latitude"])),
                    float(np.nanmax(ds["latitude"])),
                )
            title = f"{ctx.base_title}\n{var} surface (geographic)"
            slicer_lib.gmap(
                var, cfg.cmap, gmap_limits, sub, vmin=None, vmax=None, title=title
            )

        elif action == "cmap":
            prompt_cmap(cfg, prompt_prefix=f"[surface-{var}-cmap]")

        elif action == "save":
            # save only the chosen 2D variable
            prompt_save(sub[[var]])  # Dataset with just that var


# -----------------------
# App loop
# -----------------------
def run(cfg: AppConfig, ctx: DsContext) -> None:
    if cfg.verbose:
        tip("Use 'meta' to verify geospatial attrs, then 'subset' to explore.")
    while True:
        choice = prompt_menu(cfg.verbose, ctx)
        if choice in ("exit", "back", ""):
            sys.exit(0)
        if choice == "help":
            usage()
            continue
        if choice == "meta":
            show_metadata(ctx)
            continue
        if choice == "range":
            show_ranges(ctx)
            continue
        if choice == "subset":
            handle_subset(cfg, ctx)
            continue
        if choice == "map":
            show_coverage_map(cfg, ctx)
            continue
        logger.info("Unrecognized option.")


# -----------------------
# Entrypoint
# -----------------------
def main():
    cfg, meta_only = parse_args(sys.argv[1:])
    if lib.check_file_type(cfg.input_file).get("engine") != "netcdf":
        logger.error(f"[ERR] {cfg.input_file} is not a valid netCDF file")
        sys.exit(1)

    ctx = load_dataset(cfg)
    if meta_only:
        show_metadata(ctx)
        return

    run(cfg, ctx)


if __name__ == "__main__":
    main()
