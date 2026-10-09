#!/usr/bin/env python3
# =============================================================================
# 🌾 WV Tree Counter – R19h Final (Basemap Edition)
# Author: Ashley R. Mitchell – U.S. Dept. of the Interior, OSMRE (2025)
#
# DESCRIPTION:
#   A fully automated vegetation and canopy analysis workflow for WV SMCRA permits.
#   This tool integrates NAIP imagery, DeepForest canopy detection, and NDVI/GRVI
#   computation to assess vegetation coverage and canopy density within a mining
#   permit boundary. It also generates visual NDVI maps, DeepForest shapefiles,
#   and summary CSVs across runs for long-term trend analysis.
#
# FEATURES:
#   ✅ Auto-installs required Python libraries (Rich, Rasterio, GeoPandas, etc.; DeepForest is optional)
#   ✅ Queries permit boundaries directly from WVDEP mining shapefile
#   ✅ Computes NDVI (or GRVI if NIR band unavailable)
#   ✅ Adaptive thresholding if no vegetation detected at default 0.25 NDVI
#   ✅ Generates per-tile and permit-wide statistics and maps
#   ✅ Adds basemap, scalebar, and north arrow for context
#   ✅ DeepForest canopy detection (optional)
#   ✅ Logs all operations and appends summary to summary.csv
#   ✅ Produces well-formatted reports in `data/<PERMIT>/results`
#
# OUTPUTS:
#   📂 data/
#       ├── permits/                  # WVDEP mining boundary shapefile
#       ├── naip_files/
#       │   ├── index/                # WV NAIP index shapefile
#       │   ├── downloads/            # Temporary zip/tif ingestion folder
#       │   └── *.tif                 # Processed NAIP tiles
#       ├── <PERMIT>/
#       │   ├── results/
#       │   │   ├── maps/             # NDVI rasters + PNGs
#       │   │   ├── deepforest/       # Canopy shapefiles
#       │   │   └── results.txt       # Summary + plain English explanation
#       └── summary.csv               # Cumulative record across all runs
#
# DEPENDENCIES (Auto-installed on first run):
#   rich, tqdm, geopandas, rasterio, shapely, numpy, matplotlib, pandas,
#   contextily, matplotlib-scalebar
#   Optional: deepforest>=1.4,<3 (install it yourself to enable canopy detection)
#
# HOW TO RUN:
#   1️⃣ Activate your Python or Conda environment:
#       conda activate smcra_wv
#   2️⃣ Run the script:
#       python WV_TreeCounter_R19h_Final_Basemap.py
#   3️⃣ Enter a WV permit ID (e.g. S500806)
#   4️⃣ Follow on-screen prompts to place NAIP tiles in `data/naip_files/downloads/`
#   5️⃣ Results will be automatically generated and summarized.
#
# NOTE:
#   • The tool can process both NDVI (red/NIR) and fallback GRVI (red/green)
#   • DeepForest requires a compatible GPU or Apple MPS for best performance.
#   • Default coordinate system: EPSG:26917 (UTM Zone 17N)
#
# Version History:
#   R19a – Permit mosaic prototype
#   R19c – NDVI w/ Legend + Permit Acreage
#   R19e – Safe NDVI + DeepForest Recovery
#   R19h – Final Basemap Edition (this release)
# =============================================================================

import os, sys, re, zipfile, datetime, subprocess, importlib, warnings, shutil, inspect
from pathlib import Path
import numpy as np, geopandas as gpd, rasterio, rasterio.mask, rasterio.features
from rasterio.merge import merge as rio_merge
from rasterio.vrt import WarpedVRT
from rasterio.enums import ColorInterp, Resampling
from shapely.geometry import box, mapping
from tqdm import tqdm
import textwrap

# ---------- Auto-install dependencies ----------
# DeepForest (PyTorch) is NOT installed here: it is heavy and only needed when the
# user opts in to canopy detection. It is imported lazily in run_deepforest().
def ensure_package(pkg, pip_name=None):
    try:
        importlib.import_module(pkg)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pip_name or pkg])
        importlib.invalidate_caches()

for pkg, pip_name in [("rich", None), ("tqdm", None), ("geopandas", None), ("rasterio", None),
                      ("shapely", None), ("numpy", None), ("matplotlib", None),
                      ("matplotlib_scalebar", "matplotlib-scalebar"), ("contextily", None),
                      ("pandas", None)]:
    ensure_package(pkg, pip_name)

# ---------- Imports after installation ----------
from rich.console import Console
from rich.panel import Panel
from rich.prompt import Prompt
from rich.theme import Theme
from rich.table import Table
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import contextily as ctx
import pandas as pd

try:
    from wv_tree_counter.states import STATES, get_state
    from wv_tree_counter import sources
except ImportError:  # running main.py directly as a script
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from states import STATES, get_state
    import sources

# ---------- UI ----------
theme = Theme({
    "info": "steel_blue3",
    "warn": "dark_orange3",
    "good": "yellow3",
    "hdr": "bold steel_blue1",
    "gold": "bold yellow3"
})
console = Console(theme=theme)

# ---------- Config ----------
ROOT = Path("data")
PERMIT_DIR = ROOT / "permits"
NAIP_DIR = ROOT / "naip_files"
DL_DIR = NAIP_DIR / "downloads"
INDEX_DIR = NAIP_DIR / "index"

PERMIT_URL = "https://tagis.dep.wv.gov/data/vector/SHP_WVDEP_GIS_data_mining_reclamation_permit_boundary.zip"
NAIP_INDEX_URL = "https://www.fpacbc.usda.gov/sites/default/files/2024-10/wv_naip22qq.zip"
SUMMARY_CSV = "summary.csv"
LOGFILE = "runlog.txt"

# Coordinate reference systems
WGS84 = "EPSG:4326"
AREA_CRS = "EPSG:5070"           # NAD83 / CONUS Albers – equal-area, valid across all of Appalachia
SQM_PER_ACRE = 4046.8564224

# Vegetation thresholds
VEG_THRESHOLD = 0.25             # fixed threshold used for Vegetation cover (%) so permits stay comparable
ADAPTIVE_THRESHOLDS = [0.20, 0.15, 0.10]
MIN_VEG_PIXELS = 500

# Files the tool itself writes; never treat these as raw imagery
DERIVED_SUFFIXES = ("_NDVI", "_GRVI", "_rgb_mosaic", "_mosaic")
IMAGERY_EXTS = (".tif", ".tiff")


def ensure_dirs():
    for d in [ROOT, PERMIT_DIR, NAIP_DIR, DL_DIR, INDEX_DIR]:
        d.mkdir(parents=True, exist_ok=True)

# ---------- Log ----------
def log(msg):
    ts = datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    line = f"{ts} {msg}"
    console.print(f"[info]{line}[/]")
    with open(LOGFILE, "a") as f: f.write(line + "\n")

# ---------- Download/Extract ----------
def unzip(zip_path, extract_to):
    with zipfile.ZipFile(zip_path, "r") as z:
        for n in tqdm(z.namelist(), desc=f"Extracting {os.path.basename(zip_path)}"):
            z.extract(n, extract_to)

def download(url, dest):
    import requests
    try:
        r = requests.get(url, stream=True, timeout=120)
    except Exception as e:
        log(f"⚠️ Download error for {url}: {e}")
        return False
    if r.status_code != 200:
        return False
    with open(dest, "wb") as f:
        for c in r.iter_content(1 << 20):
            if c:
                f.write(c)
    return True

def find_shp(folder):
    for f in folder.glob("*.shp"):
        return f
    return None

def ensure_permits():
    shp = find_shp(PERMIT_DIR)
    if shp:
        return shp
    z = PERMIT_DIR / "wv_permits.zip"
    log("📡 Downloading WVDEP permits…")
    if not download(PERMIT_URL, z):
        sys.exit("❌ Permit shapefile download failed.")
    unzip(z, PERMIT_DIR)
    z.unlink(missing_ok=True)
    shp = find_shp(PERMIT_DIR)
    if not shp:
        sys.exit("❌ Permit shapefile missing after extraction.")
    return shp

def ensure_index():
    shp = find_shp(INDEX_DIR)
    if shp:
        return shp
    z = INDEX_DIR / "wv_naip22qq.zip"
    log("📥 Downloading WV NAIP index…")
    if not download(NAIP_INDEX_URL, z):
        console.print(f"[warn]Manual download required:\n{NAIP_INDEX_URL}[/warn]")
        input(f"Place the zip (or the extracted shapefile) into {INDEX_DIR} and press Enter…")
    if z.exists():
        unzip(z, INDEX_DIR)
        z.unlink(missing_ok=True)
    shp = find_shp(INDEX_DIR)
    if not shp:
        sys.exit("❌ NAIP index missing after extraction.")
    return shp

# ---------- CRS helpers ----------
def to_crs(geom, src_crs, dst_crs):
    """Reproject a single Shapely geometry."""
    return gpd.GeoSeries([geom], crs=src_crs).to_crs(dst_crs).iloc[0]

def utm_crs_for(geom, crs):
    """NAD83 UTM zone (EPSG:269xx) containing the geometry's centroid.
    Appalachian coal country spans zones 16–18, so this can't be hard-coded."""
    c = to_crs(geom, crs, WGS84).centroid
    zone = int((c.x + 180) // 6) + 1
    return rasterio.crs.CRS.from_epsg(26900 + zone)

# ---------- Lookup info ----------
def make_lookup_info(perm_gdf, idx_path, pid, outdir):
    geom = perm_gdf.geometry.iloc[0]
    idx = gpd.read_file(idx_path)
    if idx.crs != perm_gdf.crs:
        idx = idx.to_crs(perm_gdf.crs)
    hits = idx[idx.intersects(geom)]
    tiles = []
    for _, r in hits.iterrows():
        for key in ["FileName","filename","QQNAME","FQUAD_NAME"]:
            if key in r.index and isinstance(r[key], str):
                tiles.append(Path(r[key]).stem.strip())
                break
    centroid_4326 = to_crs(geom.centroid, perm_gdf.crs, WGS84)
    info = outdir / "lookup_info.txt"
    with open(info, "w") as f:
        f.write(f"Permit: {pid}\nCentroid (4326): {centroid_4326.y:.6f}, {centroid_4326.x:.6f}\n\n")
        f.write("Suggested NAIP tiles:\n")
        for t in tiles:
            f.write(f" - {t}.tif\n")
        f.write(f"\nDrop .zip or .tif tiles into: {DL_DIR}\n")
    console.print(Panel.fit(
        f"[gold]Centroid:[/gold] {centroid_4326.y:.6f}, {centroid_4326.x:.6f}\n"
        f"[gold]Tiles found:[/gold] {len(tiles)} (see lookup_info.txt)",
        title=f"📍 {pid} Tile Suggestions", border_style="yellow3"))
    return tiles

# ---------- Auto-ingest ----------
def auto_ingest():
    count = 0
    zips = list(DL_DIR.glob("*.zip"))
    for z in zips:
        unzip(z, DL_DIR)
        z.unlink()
    for root, _, files in os.walk(DL_DIR):
        for f in files:
            if f.lower().endswith((".tif",".tiff",".tfw",".xml",".aux.xml",".ovr")):
                src = Path(root) / f
                dst = NAIP_DIR / f
                if not dst.exists():
                    shutil.move(src, dst)
                    count += 1
    log(f"🎯 Auto-ingested {count} files.")
    return count

# ---------- Imagery discovery ----------
def tile_name(t):
    """File name of a local Path or a (signed) COG URL."""
    return os.path.basename(str(t).split("?")[0])

def list_imagery(folder):
    """Raw imagery tiles in `folder`, excluding products this tool wrote."""
    out = []
    for p in sorted(folder.iterdir()):
        if p.suffix.lower() in IMAGERY_EXTS and not p.stem.endswith(DERIVED_SUFFIXES):
            out.append(p)
    return out

def tiles_for_permit(tifs, geom, geom_crs):
    """Tiles whose footprint intersects the permit. The permit geometry is reprojected
    into each tile's own CRS, so tiles in different UTM zones are handled correctly."""
    hits = []
    for t in tifs:
        try:
            with rasterio.open(t) as src:
                g = to_crs(geom, geom_crs, src.crs) if src.crs != geom_crs else geom
                if box(*src.bounds).intersects(g):
                    hits.append(t)
        except Exception as e:
            log(f"⚠️ Could not read {tile_name(t)}: {e}")
    return hits

def warn_mixed_years(tiles):
    """NAIP file names end in the acquisition date (…_YYYYMMDD.tif). Mixing years in one
    mosaic makes the result a blend of two different points in time."""
    years = set()
    for t in tiles:
        m = re.search(r"_((?:19|20)\d{2})\d{4}(?:_\d{8})?$", Path(tile_name(t)).stem)
        if m:
            years.add(m.group(1))
    if len(years) > 1:
        log(f"⚠️ Tiles span multiple acquisition years {sorted(years)} – keep one NAIP year per run.")
    return sorted(years)

def has_alpha(src):
    """Drone orthomosaics (Site Scan, Pix4D, WebODM…) are usually RGBA: band 4 is
    transparency, not near-infrared. Some tools also *label* a real NIR band as alpha,
    so the label is only trusted when band 4 actually looks like a mask (0/255 only)."""
    try:
        if src.count < 4 or ColorInterp.alpha not in (src.colorinterp or []):
            return False
        b = src.colorinterp.index(ColorInterp.alpha) + 1
        h = max(1, min(src.height, 512)); w = max(1, min(src.width, 512))
        sample = src.read(b, out_shape=(h, w))
        return bool(np.isin(np.unique(sample), [0, 255]).all())
    except Exception:
        return False

# ---------- Permit mosaic ----------
def build_mosaic(tiles, geom, geom_crs, work_crs, out_path):
    """Merge every intersecting tile into ONE seamless raster on a common grid
    (in `work_crs`), clipped to the permit. Overlapping tile buffers are taken
    once, so pixels are never double-counted.

    Returns (array[bands, rows, cols], transform, valid_mask, n_bands_used, alpha_dropped)."""
    g_work = to_crs(geom, geom_crs, work_crs)
    srcs, vrts = [], []
    alpha_dropped = False
    try:
        res = None
        for t in tiles:
            s = rasterio.open(t)
            srcs.append(s)
            if has_alpha(s):
                alpha_dropped = True
            ds = s
            if s.crs != work_crs:
                ds = WarpedVRT(s, crs=work_crs, resampling=Resampling.nearest)
                vrts.append(ds)
            r = min(abs(ds.res[0]), abs(ds.res[1]))
            res = r if res is None else min(res, r)
        # keep tile order, using the reprojected VRT where one was made
        ordered = []
        vi = 0
        for s in srcs:
            if s.crs != work_crs:
                ordered.append(vrts[vi]); vi += 1
            else:
                ordered.append(s)
        n_bands = min(d.count for d in ordered)
        if alpha_dropped:
            n_bands = min(n_bands, 3)
        arr, transform = rio_merge(ordered, bounds=g_work.bounds, res=res,
                                   nodata=0, indexes=list(range(1, n_bands + 1)))
    finally:
        for v in vrts: v.close()
        for s in srcs: s.close()

    inside = rasterio.features.geometry_mask([mapping(g_work)], out_shape=arr.shape[1:],
                                             transform=transform, invert=True)
    has_data = np.any(arr != 0, axis=0)
    valid = inside & has_data
    arr[:, ~valid] = 0

    profile = dict(driver="GTiff", height=arr.shape[1], width=arr.shape[2], count=arr.shape[0],
                   dtype=arr.dtype, crs=work_crs, transform=transform, nodata=0, compress="LZW")
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(arr)
    return arr, transform, valid, n_bands, alpha_dropped

# ---------- NDVI / GRVI ----------
def compute_index(arr, valid, transform, crs, outdir, pid, nir_available=True):
    """NDVI (R=band1, NIR=band4) when a real NIR band exists, otherwise GRVI (G-R)/(G+R).
    Pixels outside the permit or without data are NaN and excluded from every statistic."""
    if arr.shape[0] >= 4 and nir_available:
        index_type = "NDVI"
        red = arr[0].astype("float32"); nir = arr[3].astype("float32")
        index = (nir - red) / (nir + red + 1e-6)
    elif arr.shape[0] >= 3:
        index_type = "GRVI"
        red = arr[0].astype("float32"); green = arr[1].astype("float32")
        index = (green - red) / (green + red + 1e-6)
    else:
        raise ValueError(f"Unsupported band count: {arr.shape[0]}")
    index[~valid] = np.nan

    n_valid = int(valid.sum())
    # Adaptive threshold (reported for reference); pix is counted AT the threshold chosen.
    thr = VEG_THRESHOLD
    pix = int(np.count_nonzero(index > thr))
    if pix < MIN_VEG_PIXELS:
        for t in ADAPTIVE_THRESHOLDS:
            n = int(np.count_nonzero(index > t))
            if n > MIN_VEG_PIXELS:
                thr, pix = t, n
                log(f"⚠️  Adaptive threshold {t:.2f} applied to {pid}")
                break

    out_tif = outdir / f"{pid}_{index_type}.tif"
    with rasterio.open(out_tif, "w", driver="GTiff", height=index.shape[0], width=index.shape[1],
                       count=1, dtype="float32", crs=crs, transform=transform,
                       nodata=np.nan, compress="LZW") as dst:
        dst.write(index.astype("float32"), 1)

    mean = float(np.nanmean(index)) if n_valid else 0.0
    veg_pct = float(np.count_nonzero(index > VEG_THRESHOLD) / n_valid * 100) if n_valid else 0.0
    return {"tile": pid, "thr": thr, "pix": pix, "index": out_tif, "type": index_type,
            "mean": mean, "veg_pct": veg_pct, "valid_pixels": n_valid}

# ---------- Index PNG ----------
def render_index_png(index_tif, pid, mapsdir, index_type="NDVI"):
    with rasterio.open(index_tif) as src:
        m = src.read(1)
        b = src.bounds
    png = mapsdir / f"{pid}_index_mean.png"
    fig, ax = plt.subplots(figsize=(8, 6))
    show = ax.imshow(m, cmap="YlGn", vmin=0, vmax=1, extent=[b.left, b.right, b.bottom, b.top])
    plt.colorbar(show, ax=ax, label=index_type)
    plt.title(f"Vegetation Index ({index_type}) – {pid}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(png, dpi=150)
    plt.close()
    return png

# ---------- DeepForest ----------
def boxes_to_geo(df, transform, crs):
    """Convert DeepForest pixel boxes (xmin/ymin/xmax/ymax, row 0 at top) to map-coordinate polygons."""
    geoms = []
    for xmin, ymin, xmax, ymax in df[["xmin", "ymin", "xmax", "ymax"]].itertuples(index=False):
        x0 = transform.c + xmin * transform.a + ymin * transform.b
        y0 = transform.f + xmin * transform.d + ymin * transform.e
        x1 = transform.c + xmax * transform.a + ymax * transform.b
        y1 = transform.f + xmax * transform.d + ymax * transform.e
        geoms.append(box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)))
    keep = [c for c in ("label", "score") if c in df.columns]
    return gpd.GeoDataFrame(df[keep].reset_index(drop=True), geometry=geoms, crs=crs)

def run_deepforest(rgb_path, geom, geom_crs, outdir, pid):
    """Detect tree crowns on the permit RGB mosaic. Runs once on the merged mosaic
    (no double counting in tile overlaps) using tiled prediction, then georeferences
    the boxes and keeps crowns whose centre falls inside the permit.
    Returns (shapefile path or None, crown count)."""
    try:
        from deepforest import main as df_main
    except ImportError:
        log("❌ DeepForest is not installed. Run: pip install 'deepforest>=1.4,<3'")
        return None, 0

    with rasterio.open(rgb_path) as src:
        transform, crs, gsd = src.transform, src.crs, abs(src.res[0])
    if gsd > 0.3:
        log(f"⚠️ Imagery is {gsd:.2f} m/pixel. DeepForest's tree model was trained on ~0.1 m imagery; "
            "at NAIP resolution small/young crowns are missed, so treat counts as a lower bound.")

    log("🌳 Running DeepForest")
    model = df_main.deepforest()
    model.load_model("weecology/deepforest-tree")   # use_release() was removed in DeepForest 2.x

    # DeepForest 1.x names the argument raster_path; 2.x names it path.
    params = inspect.signature(model.predict_tile).parameters
    path_kw = "path" if "path" in params else "raster_path"
    try:
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None   # large permits exceed PIL's default size guard
    except ImportError:
        pass
    preds = model.predict_tile(**{path_kw: str(rgb_path)}, patch_size=400, patch_overlap=0.05)

    if preds is None or len(preds) == 0:
        log("DeepForest found no crowns.")
        return None, 0

    crowns = boxes_to_geo(preds, transform, crs)
    g = to_crs(geom, geom_crs, crs)
    crowns = crowns[crowns.geometry.centroid.within(g)].reset_index(drop=True)
    shp = outdir / f"{pid}_canopy.shp"
    if len(crowns):
        crowns.to_file(shp)
    log(f"🌳 DeepForest: {len(crowns):,} crowns inside the permit.")
    return (shp if len(crowns) else None), len(crowns)

def write_rgb(arr, transform, crs, out_path):
    rgb = arr[:3]
    with rasterio.open(out_path, "w", driver="GTiff", height=rgb.shape[1], width=rgb.shape[2],
                       count=3, dtype=rgb.dtype, crs=crs, transform=transform, nodata=0,
                       compress="LZW") as dst:
        dst.write(rgb)
    return out_path

# ---------- Interpretive summary ----------
def interpret_site(veg_pct, meanval, canopy_per_acre):
    """
    Converts numeric vegetation and canopy metrics into an Appalachian-style
    plain-language interpretation of site condition.
    """

    if veg_pct < 5 and meanval < 0.15:
        return ("The site shows very low vegetation cover and NDVI values, "
                "suggesting barren or recently disturbed mine lands with little "
                "active regrowth. Sparse canopy recovery indicates early-stage "
                "reclamation or continuing disturbance.")
    elif veg_pct < 25:
        return ("Vegetation cover remains patchy, with low to moderate NDVI. "
                "This typically represents grass-dominated reclamation or partial "
                "natural recovery. Tree canopy density is minimal but beginning "
                "to establish in sheltered areas.")
    elif veg_pct < 60:
        return ("The site exhibits moderate vegetation recovery, with NDVI "
                "values consistent with maturing herbaceous cover and early "
                "woody growth. Canopy presence suggests mixed-age succession "
                "or planted reclamation plots.")
    else:
        return ("The area shows strong vegetative recovery with high NDVI and "
                "substantial canopy density, indicating successful long-term "
                "reclamation or natural forest regeneration on former mine lands.")

# ---------- summary.csv ----------
SUMMARY_COLUMNS = ["date", "state", "permit", "permit_acres", "veg_pix", "thr", "mean_ndvi",
                   "veg_pct", "canopy", "canopy_per_acre", "lon", "lat",
                   "imagery_source", "imagery_year", "imagery_dates"]

def append_summary(row, path=None):
    """Append one row to summary.csv. Older WV-only files (no state/imagery columns)
    are upgraded in place first, with state = WV for existing rows."""
    path = Path(path or SUMMARY_CSV)
    if path.exists():
        with open(path) as f:
            header = f.readline().strip().split(",")
        if header != SUMMARY_COLUMNS:
            old = pd.read_csv(path, dtype=str)
            if "state" not in old.columns:
                old["state"] = "WV"
            for c in SUMMARY_COLUMNS:
                if c not in old.columns:
                    old[c] = ""
            extra = [c for c in old.columns if c not in SUMMARY_COLUMNS]
            old[SUMMARY_COLUMNS + extra].to_csv(path, index=False)
            log(f"📄 Upgraded {path} to the multi-state column layout.")
            header = SUMMARY_COLUMNS + extra
    else:
        header = SUMMARY_COLUMNS
        pd.DataFrame(columns=header).to_csv(path, index=False)
    pd.DataFrame([{c: row.get(c, "") for c in header}]).to_csv(path, mode="a", header=False, index=False)

# ---------- Summary + Appalachian Interpretation ----------
def summarize(pid, res, rec, canopy, geom, crs, permit_area, state="WV", imagery=None):
    """
    Summarizes NDVI, vegetation, canopy, and interpretive context for one permit.
    Writes results.txt, updates summary.csv, and prints formatted Rich table.
    Runs whether or not DeepForest was used (canopy = 0 when skipped).
    """
    veg = rec["pix"] if rec else 0
    thr = rec["thr"] if rec else 0
    meanval = rec["mean"] if rec else 0
    veg_pct = rec["veg_pct"] if rec else 0

    canopy_per_acre = canopy / permit_area if permit_area > 0 else 0

    # --- Console output table ---
    imagery = imagery or {}
    tab = Table(title=f"🌾 Tree Counter – {state} {pid}", title_style="gold")
    for n, v in [
        ("State", state),
        ("Imagery", f"{imagery.get('source', 'local')} {imagery.get('year', '')}".strip()),
        ("Permit acres", f"{permit_area:,.1f}"),
        ("Veg pixels", f"{veg:,}"),
        ("NDVI thresh (adaptive)", f"{thr:.2f}"),
        ("Mean NDVI", f"{meanval:.2f}"),
        (f"Veg. cover (% > {VEG_THRESHOLD:.2f})", f"{veg_pct:.1f}"),
        ("Canopy crowns", f"{canopy:,}")
    ]:
        tab.add_row(n, v)
    console.print(tab)

    # --- Generate human-readable interpretation ---
    interp_text = interpret_site(veg_pct, meanval, canopy_per_acre)

    # --- Write results.txt ---
    results_file = res / "results.txt"
    with open(results_file, "w") as f:
        f.write(f"State: {state}\n")
        f.write(f"Permit: {pid}\n")
        f.write(f"Imagery: {imagery.get('source', 'local')} {imagery.get('year', '')} "
                f"{', '.join(imagery.get('dates', []))}\n")
        f.write(f"Permit area (acres): {permit_area:.1f}\n")
        f.write(f"Veg pixels (at adaptive threshold): {veg:,}\n")
        f.write(f"NDVI threshold (adaptive): {thr:.2f}\n")
        f.write(f"Mean NDVI: {meanval:.2f}\n")
        f.write(f"Vegetation cover (% of permit pixels > {VEG_THRESHOLD:.2f}): {veg_pct:.1f}\n")
        f.write(f"Canopy crowns: {canopy:,}\n")
        f.write(f"Canopy per acre: {canopy_per_acre:.4f}\n\n")
        f.write("Plain-language interpretation:\n")
        f.write(textwrap.fill(interp_text, width=80))
        f.write("\n")

    # --- Append to running summary CSV ---
    cent = to_crs(geom, crs, WGS84).centroid
    append_summary({
        "date": str(datetime.date.today()), "state": state, "permit": pid,
        "permit_acres": f"{permit_area:.1f}", "veg_pix": veg, "thr": f"{thr:.2f}",
        "mean_ndvi": f"{meanval:.2f}", "veg_pct": f"{veg_pct:.1f}", "canopy": canopy,
        "canopy_per_acre": f"{canopy_per_acre:.4f}", "lon": f"{cent.x:.6f}", "lat": f"{cent.y:.6f}",
        "imagery_source": imagery.get("source", "local"), "imagery_year": imagery.get("year", ""),
        "imagery_dates": ";".join(imagery.get("dates", [])),
    })

    # --- Log and colorized completion panel ---
    console.print(Panel.fit(
        f"[gold]✅ Complete![/gold]\n"
        f"[steel_blue3]Permit:[/steel_blue3] {state} {pid}\n"
        f"[steel_blue3]Area:[/steel_blue3] {permit_area:,.1f} acres\n"
        f"[steel_blue3]NDVI mean:[/steel_blue3] {meanval:.2f}\n"
        f"[steel_blue3]Vegetation cover:[/steel_blue3] {veg_pct:.1f}%\n"
        f"[steel_blue3]Canopy crowns:[/steel_blue3] {canopy:,}\n\n"
        f"[gold]Summary → {res}/results.txt[/gold]\n\n"
        f"[italic]{interp_text}[/italic]",
        title="🌾 Tree Counter – Appalachian Summary",
        border_style="yellow3", width=80
    ))
    return {"veg_pix": veg, "thr": thr, "mean_ndvi": meanval, "veg_pct": veg_pct,
            "canopy": canopy, "canopy_per_acre": canopy_per_acre}

# ---------- Map Composer ----------
def make_map(index_tif, geom, geom_crs, mapsdir, pid, index_type="NDVI", basemap=True):
    """NDVI map drawn in map coordinates so the raster, permit outline, basemap and
    scalebar all line up."""
    if not index_tif or not Path(index_tif).exists():
        return None
    with rasterio.open(index_tif) as src:
        arr = src.read(1)
        b = src.bounds
        rcrs = src.crs
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.set_xlim(b.left, b.right); ax.set_ylim(b.bottom, b.top)
    ink = "black"
    if basemap:
        try:
            ctx.add_basemap(ax, crs=rcrs, source=ctx.providers.Esri.WorldImagery, zorder=0)
            ink = "white"
        except Exception as e:
            log(f"⚠️ Basemap unavailable ({e}); map drawn without it.")
    show = ax.imshow(arr, cmap="YlGn", vmin=0, vmax=1, alpha=0.85, zorder=1,
                     extent=[b.left, b.right, b.bottom, b.top])
    plt.colorbar(show, ax=ax, label=index_type)
    plt.title(f"{index_type} Map – {pid}")
    g = to_crs(geom, geom_crs, rcrs)
    gpd.GeoSeries([g], crs=rcrs).boundary.plot(ax=ax, color=ink, linewidth=1.5, zorder=2)
    ax.add_artist(ScaleBar(dx=1, units="m", dimension="si-length", location="lower left"))
    ax.text(0.95, 0.1, "N\n↑", transform=ax.transAxes, ha="center", va="center",
            fontsize=12, color=ink, zorder=3)
    ax.set_axis_off()
    plt.tight_layout()
    outpng = mapsdir / f"{pid}_NDVI_map.png"
    plt.savefig(outpng, dpi=150)
    plt.close()
    return outpng

# ---------- Permit Geometry and Area ----------
def get_permit_area(permits, pid):
    """
    Returns the Shapely geometry of the selected permit and its area in acres.
    Area is measured in an equal-area projection (EPSG:5070). Web Mercator
    (EPSG:3857) overstates area by ~1.4–1.7x at Appalachian latitudes.
    """
    sel = permits[permits["permit_id"].astype(str).str.upper() == pid]
    if sel.empty:
        raise ValueError(f"Permit {pid} not found in permit shapefile.")
    return permit_geometry(sel, pid)

def permit_geometry(sel, pid=""):
    """Dissolve every polygon for a permit (GeoMine often returns several) into one
    geometry and measure it in acres (equal-area EPSG:5070)."""
    geom = sources.union(sel.geometry)
    if not geom.is_valid:
        geom = geom.buffer(0)
    area_m2 = gpd.GeoSeries([geom], crs=sel.crs).to_crs(AREA_CRS).area.iloc[0]
    permit_area = area_m2 / SQM_PER_ACRE
    log(f"📐 Calculated permit area for {pid}: {permit_area:.1f} acres ({len(sel)} polygon(s))")
    return geom, permit_area

# ---------- Analysis (no prompts; reusable for batch runs) ----------
def analyze_permit(pid, geom, crs, permit_area, tifs, res, maps, deep, run_df=False, basemap=True,
                   state="WV", imagery=None):
    """Mosaic → index → maps → optional DeepForest → summary, for one permit."""
    tiles = tiles_for_permit(tifs, geom, crs)
    if not tiles:
        log(f"❌ None of the {len(tifs)} imagery tiles intersect permit {pid}.")
        return None
    log(f"🧩 {len(tiles)} tile(s) intersect {pid}: {', '.join(tile_name(t) for t in tiles)}")
    warn_mixed_years(tiles)

    work_crs = utm_crs_for(geom, crs)
    mosaic_path = res / f"{pid}_mosaic.tif"
    arr, transform, valid, n_bands, alpha_dropped = build_mosaic(tiles, geom, crs, work_crs, mosaic_path)
    if alpha_dropped:
        log("⚠️ Band 4 is an alpha (transparency) band, not NIR – using GRVI on RGB.")
    if not valid.any():
        log(f"❌ Imagery has no valid pixels inside permit {pid}.")
        return None

    rec = compute_index(arr, valid, transform, work_crs, res, pid, nir_available=not alpha_dropped)
    render_index_png(rec["index"], pid, maps, rec["type"])
    make_map(rec["index"], geom, crs, maps, pid, rec["type"], basemap=basemap)

    canopy = 0
    if run_df:
        rgb_path = write_rgb(arr, transform, work_crs, res / f"{pid}_rgb_mosaic.tif")
        _, canopy = run_deepforest(rgb_path, geom, crs, deep, pid)

    return summarize(pid, res, rec, canopy, geom, crs, permit_area, state=state, imagery=imagery)

# ---------- Permit lookup (any state) ----------
def resolve_permit(state, pid, permit_source=None, permit_file=None, id_field=None):
    """Return (geometry, crs, acres) for a permit in any supported state.

    Sources: WVDEP shapefile (WV default), OSMRE GeoMine (PA, OH, VA, KY, TN, AL, or WV
    with --permit-source geomine), or a local boundary file (--permit-file, needed for MD)."""
    state, cfg = get_state(state)
    src = "file" if permit_file else (permit_source or cfg["permit_source"])
    if src == "file":
        if not permit_file or not id_field:
            raise ValueError(f"{cfg['name']} needs --permit-file and --id-field "
                             f"(permit boundaries from {cfg['regulator']}).")
        sel = sources.load_permit_file(permit_file, id_field, pid)
        if sel.empty:
            raise ValueError(f"Permit {pid} not found in {permit_file} ({id_field}).")
    elif src == "wvdep":
        permits = gpd.read_file(ensure_permits())
        sel = permits[permits["permit_id"].astype(str).str.upper() == pid]
        if sel.empty:
            raise ValueError(f"Permit {pid} not found in the WVDEP permit shapefile.")
    elif src == "geomine":
        if cfg["geomine_contact"] is None:
            raise ValueError(f"{cfg['name']} is not in GeoMine; use --permit-file.")
        sel, suggestions = sources.fetch_permit_geomine(cfg["geomine_contact"], pid)
        if sel.empty:
            hint = f" Similar IDs: {', '.join(suggestions)}" if suggestions else ""
            raise ValueError(f"Permit {pid} not found in GeoMine for {cfg['name']}.{hint}")
        if "permittee" in sel.columns:
            names = sorted(set(sel["permittee"].dropna().astype(str)))
            log(f"ℹ️ GeoMine: {len(sel)} polygon(s); permittee: {', '.join(names) or 'n/a'}")
    else:
        raise ValueError(f"Unknown permit source '{src}'")
    geom, acres = permit_geometry(sel, pid)
    return geom, sel.crs, acres

# ---------- Imagery acquisition ----------
def get_imagery(mode, state, pid, geom, crs, pdir, year=None, interactive=False):
    """Return (tiles, imagery_meta). mode 'stac' streams NAIP COGs from the Planetary
    Computer catalog (works for every state, no downloads); 'local' uses .tif files
    in data/naip_files/."""
    if mode == "stac":
        g4326 = to_crs(geom, crs, WGS84)
        info = sources.naip_for_permit(g4326, year=year)
        log(f"🛰️ NAIP {info['year']} ({', '.join(info['states']).upper()}): {len(info['hrefs'])} tile(s), "
            f"{'/'.join(f'{g:g}' for g in info['gsd'])} m, dates {', '.join(info['dates'])}")
        log(f"   Other NAIP years here: {', '.join(y for y in info['available_years'] if y != info['year'])}")
        if info["off_season_dates"]:
            log(f"⚠️ Catalog date(s) {', '.join(info['off_season_dates'])} fall outside May–Sept. "
                "NAIP catalog dates are sometimes delivery dates – confirm the imagery is leaf-on.")
        with open(pdir / "imagery_info.txt", "w") as f:
            f.write(f"NAIP year: {info['year']}\nDates: {', '.join(info['dates'])}\n"
                    f"GSD (m): {info['gsd']}\nItems:\n" + "\n".join(f" - {i}" for i in info["ids"]) + "\n")
        return info["hrefs"], {"source": "naip-stac", "year": info["year"], "dates": info["dates"]}

    # local files
    if interactive and state == "WV":
        make_lookup_info(gpd.GeoDataFrame(geometry=[geom], crs=crs), ensure_index(), pid, pdir)
    if interactive:
        console.print(Panel.fit(f"Place NAIP .zip/.tif files for this permit into:\n{DL_DIR}\n\nPress Enter to continue.",
                                title="⏸️ Pause for Downloads", border_style="yellow3"))
        input()
    auto_ingest()
    tifs = list_imagery(NAIP_DIR)
    years = warn_mixed_years(tiles_for_permit(tifs, geom, crs)) if tifs else []
    return tifs, {"source": "local", "year": years[0] if len(years) == 1 else "", "dates": []}

# ---------- One permit, end to end ----------
def run_permit(state, pid, imagery="stac", year=None, run_df=False, basemap=True,
               permit_source=None, permit_file=None, id_field=None, interactive=False):
    state, cfg = get_state(state)
    pid = pid.strip().upper()
    console.print(f"[hdr]── {cfg['name']} permit {pid} ──[/hdr]")
    geom, crs, permit_area = resolve_permit(state, pid, permit_source, permit_file, id_field)

    pdir = ROOT / state / pid
    res = pdir / "results"
    maps = res / "maps"
    deep = res / "deepforest"
    for d in [pdir, res, maps, deep]:
        d.mkdir(parents=True, exist_ok=True)

    tiles, meta = get_imagery(imagery, state, pid, geom, crs, pdir, year=year, interactive=interactive)
    if not tiles:
        raise RuntimeError("No imagery found.")
    console.print("[info]🌿 Building permit mosaic and computing vegetation index…[/info]")
    out = analyze_permit(pid, geom, crs, permit_area, tiles, res, maps, deep, run_df=run_df,
                         basemap=basemap, state=state, imagery=meta)
    console.print(Panel.fit(f"[good]✅ Complete![/good]\nResults → {res}", title="🌾 Done", border_style="yellow3"))
    return out

def read_permit_list(path, default_state):
    """One permit per line: 'PERMIT' or 'STATE,PERMIT'. Blank lines and # comments skipped."""
    jobs = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#")[0].strip()
        if not line:
            continue
        parts = [p.strip() for p in line.split(",")]
        if len(parts) >= 2 and parts[0].upper() in STATES:
            jobs.append((parts[0].upper(), parts[1]))
        else:
            if not default_state:
                raise ValueError(f"'{line}' has no state; add 'ST,PERMIT' or pass --state.")
            jobs.append((default_state, parts[0]))
    return jobs

# ---------- CLI ----------
def build_parser():
    import argparse
    ap = argparse.ArgumentParser(
        prog="wv-tree-counter",
        description="Vegetation (NDVI) and canopy analysis on Appalachian coal mining permits. "
                    "Run with no arguments for the interactive prompts.")
    ap.add_argument("--state", help=f"Two-letter state: {', '.join(STATES)}")
    ap.add_argument("--permit", action="append", help="Permit ID (repeatable)")
    ap.add_argument("--permit-list", help="Text file: one 'PERMIT' or 'STATE,PERMIT' per line")
    ap.add_argument("--imagery", choices=["stac", "local"], default="stac",
                    help="stac = stream NAIP from the Planetary Computer catalog (default); "
                         "local = .tif files in data/naip_files/")
    ap.add_argument("--year", help="NAIP year (default: newest year that fully covers the permit)")
    ap.add_argument("--deepforest", action="store_true", help="Also run DeepForest crown detection")
    ap.add_argument("--no-basemap", action="store_true", help="Skip the Esri basemap on maps")
    ap.add_argument("--permit-source", choices=["wvdep", "geomine"],
                    help="Override the state's default permit source")
    ap.add_argument("--permit-file", help="Local permit boundary file (shapefile/GeoPackage/GeoJSON)")
    ap.add_argument("--id-field", help="Permit ID column in --permit-file")
    ap.add_argument("--list-states", action="store_true", help="Show supported states and exit")
    return ap

def show_states():
    tab = Table(title="Supported states", title_style="gold")
    for c in ["State", "Permit source", "Regulator", "Coal fields"]:
        tab.add_column(c)
    for code, cfg in STATES.items():
        tab.add_row(code, cfg["permit_source"], cfg["regulator"], cfg["coal_fields"])
    console.print(tab)

def run_cli(argv=None):
    """Parse arguments and run; returns the list of per-permit summaries."""
    ensure_dirs()
    args = build_parser().parse_args(argv)
    console.print(Panel.fit("[hdr]🌾 Tree Counter – Appalachian Coal Country (WV · PA · OH · MD · VA · KY · TN · AL)[/hdr]",
                            border_style="yellow3", width=80))
    if args.list_states:
        show_states()
        return []

    opts = dict(imagery=args.imagery, year=args.year, run_df=args.deepforest,
                basemap=not args.no_basemap, permit_source=args.permit_source,
                permit_file=args.permit_file, id_field=args.id_field)

    jobs = []
    if args.permit_list:
        jobs += read_permit_list(args.permit_list, args.state.upper() if args.state else None)
    if args.permit:
        if not args.state:
            sys.exit("❌ --permit needs --state.")
        jobs += [(args.state.upper(), p) for p in args.permit]

    if not jobs:  # interactive mode
        show_states()
        state = Prompt.ask("🗺️ State", choices=list(STATES), default="WV")
        cfg = STATES[state]
        if cfg["permit_source"] == "file" and not opts["permit_file"]:
            opts["permit_file"] = Prompt.ask(f"Path to {cfg['name']} permit boundary file")
            opts["id_field"] = Prompt.ask("Permit ID field name")
        opts["imagery"] = Prompt.ask("Imagery", choices=["stac", "local"], default="stac")
        opts["run_df"] = Prompt.ask("Run DeepForest? (y/N)", default="n").lower().startswith("y")
        for attempt in range(3):
            pid = Prompt.ask(f"🔢 {cfg['name']} permit number").strip().upper()
            try:
                return [run_permit(state, pid, interactive=True, **opts)]
            except ValueError as e:
                console.print(f"[warn]{e}[/warn]")
        sys.exit("❌ No matching permit.")

    results, failures = [], []
    for state, pid in jobs:
        try:
            results.append(run_permit(state, pid, **opts))
        except Exception as e:
            log(f"❌ {state} {pid}: {e}")
            failures.append((state, pid, str(e)))
    console.print(Panel.fit(f"[gold]Batch finished:[/gold] {len(results)} succeeded, {len(failures)} failed"
                            + "".join(f"\n  • {s} {p}: {m}" for s, p, m in failures),
                            title="🌾 Batch", border_style="yellow3"))
    if failures and not results:
        raise SystemExit(1)
    return results

def main(argv=None):
    """Console-script entry point (returns None so the exit code stays 0 on success)."""
    run_cli(argv)


# ---------- Entrypoint ----------
if __name__ == "__main__":
    main()
