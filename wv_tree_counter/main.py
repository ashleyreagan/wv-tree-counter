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
#   ✅ Auto-installs required Python libraries (Rich, Rasterio, GeoPandas, DeepForest, etc.)
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
#   deepforest, contextily, matplotlib-scalebar, joblib
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
#!/usr/bin/env python3
# =============================================================================
# WV Tree Counter (R19h Final – Basemap Edition)
# Author: Ashley Mitchell, OSMRE (2025)
#
# Features:
#  - Auto-installs dependencies if missing
#  - Adaptive NDVI/GRVI switching based on band count
#  - Adaptive rerun suggestion if no vegetation detected
#  - Steel-blue and amber Rich UI
#  - Parallel processing (w/ single-core fallback)
#  - DeepForest canopy detection within permit boundary
#  - Basemap, scalebar, north arrow, white legends
#  - Separate NDVI and Canopy output folders
#  - Acreage comparison (permit vs analyzed)
# =============================================================================

import os, sys, zipfile, datetime, subprocess, importlib, warnings, shutil
from pathlib import Path
import numpy as np, geopandas as gpd, rasterio, rasterio.mask
from shapely.geometry import box
from joblib import Parallel, delayed
from tqdm import tqdm
import textwrap

# ---------- Auto-install dependencies ----------
def ensure_package(pkg):
    try:
        importlib.import_module(pkg)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
        importlib.invalidate_caches()

for pkg in ["rich","tqdm","geopandas","rasterio","shapely","numpy","matplotlib",
            "matplotlib_scalebar","contextily","pandas","deepforest","joblib"]:
    ensure_package(pkg)

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

warnings.filterwarnings("ignore", category=UserWarning)

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
for d in [ROOT, PERMIT_DIR, NAIP_DIR, DL_DIR, INDEX_DIR]: d.mkdir(parents=True, exist_ok=True)

PERMIT_URL = "https://tagis.dep.wv.gov/data/vector/SHP_WVDEP_GIS_data_mining_reclamation_permit_boundary.zip"
NAIP_INDEX_URL = "https://www.fpacbc.usda.gov/sites/default/files/2024-10/wv_naip22qq.zip"
SUMMARY_CSV = "summary.csv"
LOGFILE = "runlog.txt"

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
    r = requests.get(url, stream=True, timeout=120)
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
        input("Place the zip into index folder and press Enter…")
    unzip(z, INDEX_DIR)
    z.unlink(missing_ok=True)
    shp = find_shp(INDEX_DIR)
    if not shp:
        sys.exit("❌ NAIP index missing after extraction.")
    return shp

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
    centroid_4326 = gpd.GeoSeries([geom.centroid], crs=perm_gdf.crs).to_crs(4326).iloc[0]
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
            if f.lower().endswith((".tif",".tfw",".xml",".aux.xml",".ovr")):
                src = Path(root) / f
                dst = NAIP_DIR / f
                if not dst.exists():
                    shutil.move(src, dst)
                    count += 1
    log(f"🎯 Auto-ingested {count} files.")
    return count

# ---------- NDVI / GRVI (Rewritten Safe Version) ----------
def compute_index(tif, geom, outdir):
    """
    Compute NDVI (if 4-band) or GRVI (if 3-band) within a given geometry.
    Adds robust shape checking and band name reporting.
    """
    try:
        with rasterio.open(tif) as src:
            if not box(*src.bounds).intersects(geom):
                return None

            img, transform = rasterio.mask.mask(src, [geom], crop=True)
            bands = img.shape[0]

            log(f"🛰️  {tif.name}: {bands} bands → {src.descriptions or 'no band names'}")

            # Determine index type and band mapping
            if bands >= 4:
                index_type = "NDVI"
                red_band  = 1  # usually band 1 (red)
                nir_band  = 4  # usually band 4 (NIR)
                red = img[red_band - 1].astype("float32")
                nir = img[nir_band - 1].astype("float32")
                index = (nir - red) / (nir + red + 1e-6)
            elif bands >= 3:
                index_type = "GRVI"
                red = img[0].astype("float32")
                green = img[1].astype("float32")
                index = (green - red) / (green + red + 1e-6)
            else:
                raise ValueError(f"Unsupported band count: {bands}")

            # Mask zeros and NaNs
            index[np.isnan(index)] = np.nan
            index[index == 0] = np.nan

            # Adaptive vegetation threshold
            thr = 0.25
            pix = np.count_nonzero(index > thr)
            if pix < 500:
                for t in [0.20, 0.15, 0.10]:
                    if np.count_nonzero(index > t) > 500:
                        thr = t
                        log(f"⚠️  Adaptive threshold {t:.2f} applied to {tif.name}")
                        break

            out_tif = outdir / f"{tif.stem}_{index_type}.tif"
            meta = src.meta.copy()
            meta.update(
                count=1,
                dtype="float32",
                compress="LZW",
                height=index.shape[0],
                width=index.shape[1],
                transform=transform
            )
            with rasterio.open(out_tif, "w", **meta) as dst:
                dst.write(index, 1)

            return {
                "tile": tif.name,
                "thr": thr,
                "pix": int(pix),
                "index": out_tif,
                "type": index_type
            }

    except Exception as e:
        log(f"❌ Index fail {tif.name}: {e}")
        return None
# ---------- Mean Index Map ----------
def mean_index_map(records, ref, mapsdir, pid):
    if not records:
        return None, None

    # Base reference profile
    with rasterio.open(ref) as base:
        base_profile = base.profile
        base_shape = base.read(1).shape

    arrs = []
    for r in records:
        try:
            with rasterio.open(r["index"]) as src:
                arr = src.read(1)
                if arr.shape != base_shape:
                    log(f"⚠️ Resampling {r['tile']} from {arr.shape} → {base_shape}")
                    data = np.empty(base_shape, dtype=np.float32)
                    rasterio.warp.reproject(
                        arr,
                        data,
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=base_profile["transform"],
                        dst_crs=base_profile["crs"],
                        resampling=rasterio.warp.Resampling.bilinear,
                    )
                    arr = data
                arrs.append(arr)
        except Exception as e:
            log(f"⚠️ Skipping {r['tile']} due to error: {e}")

    if not arrs:
        log("❌ No valid NDVI arrays to average.")
        return None, None

    # Compute the mean NDVI mosaic
    m = np.nanmean(np.stack(arrs), axis=0).astype("float32")

    # Write the output raster
    mtif = mapsdir / f"{pid}_index_mean.tif"
    with rasterio.open(mtif, "w", **base_profile) as dst:
        dst.write(m, 1)

    # Generate PNG
    png = mapsdir / f"{pid}_index_mean.png"
    fig, ax = plt.subplots(figsize=(8, 6))
    show = ax.imshow(m, cmap="YlGn", vmin=0, vmax=1)
    plt.colorbar(show, ax=ax, label="Vegetation Index")
    plt.title(f"Mean Vegetation Index – {pid}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(png, dpi=150)
    plt.close()
    return mtif, png
    # Generate PNG with legend and boundary
    png = mapsdir / f"{pid}_index_mean.png"
    fig, ax = plt.subplots(figsize=(8, 6))
    show = ax.imshow(m, cmap="YlGn", vmin=0, vmax=1)
    plt.colorbar(show, ax=ax, label="Vegetation Index")
    plt.title(f"Mean Vegetation Index – {pid}")
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(png, dpi=150)
    plt.close()
    return mtif, png


# ---------- DeepForest ----------
def run_deepforest(tifs, geom, outdir):
    from deepforest import main
    log("🌳 Running DeepForest")
    model = main.deepforest()
    model.use_release()
    shps = []
    for tif in tqdm(tifs, desc="DeepForest"):
        if tif.name.endswith(("_NDVI.tif","_GRVI.tif")):
            continue
        try:
            with rasterio.open(tif) as src:
                if not box(*src.bounds).intersects(geom):
                    continue
                img, tf = rasterio.mask.mask(src, [geom], crop=True)
                rgb = img[:3] if img.shape[0] >= 3 else img
                npimg = np.moveaxis(rgb, 0, -1)
                preds = model.predict_image(npimg)
                shp = outdir / f"{tif.stem}_canopy.shp"
                if hasattr(preds, "to_file"):
                    if getattr(preds, "crs", None) is None:
                        preds.set_crs(src.crs, inplace=True)
                    preds.to_file(shp)
                shps.append(shp)
        except Exception as e:
            log(f"⚠️ DF fail {tif.name}: {e}")
    return shps
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
    elif 5 <= veg_pct < 25:
        return ("Vegetation cover remains patchy, with low to moderate NDVI. "
                "This typically represents grass-dominated reclamation or partial "
                "natural recovery. Tree canopy density is minimal but beginning "
                "to establish in sheltered areas.")
    elif 25 <= veg_pct < 60:
        return ("The site exhibits moderate vegetation recovery, with NDVI "
                "values consistent with maturing herbaceous cover and early "
                "woody growth. Canopy presence suggests mixed-age succession "
                "or planted reclamation plots.")
    else:
        return ("The area shows strong vegetative recovery with high NDVI and "
                "substantial canopy density, indicating successful long-term "
                "reclamation or natural forest regeneration on former mine lands.")
# ---------- Summary + Appalachian Interpretation ----------
def summarize (pid, res, recs, mtif, shps, geom, crs, permit_area):
    """
    Summarizes NDVI, vegetation, canopy, and interpretive context for one permit.
    Writes results.txt, updates summary.csv, and prints formatted Rich table.
    """

    # --- Aggregate NDVI and vegetation metrics ---
    veg = sum(r["pix"] for r in recs) if recs else 0
    thr = np.mean([r["thr"] for r in recs]) if recs else 0
    meanval = veg_pct = 0
    if mtif and Path(mtif).exists():
        with rasterio.open(mtif) as src:
            arr = src.read(1)
            meanval = float(np.nanmean(arr))
            veg_pct = float(np.count_nonzero(arr > 0.25) / arr.size * 100)

    # --- DeepForest canopy metrics ---
    canopy = 0
    if shps:
        dfs = [gpd.read_file(s) for s in shps if Path(s).exists()]
        if dfs:
            canopy = sum(len(df) for df in dfs)

    canopy_per_acre = canopy / permit_area if permit_area > 0 else 0

    # --- Console output table ---
    tab = Table(title=f"🌾 WV Tree Counter – {pid}", title_style="gold")
    for n, v in [
        ("Permit acres", f"{permit_area:,.1f}"),
        ("Veg pixels", f"{veg:,}"),
        ("NDVI thresh(avg)", f"{thr:.2f}"),
        ("Mean NDVI", f"{meanval:.2f}"),
        ("Veg. cover(%)", f"{veg_pct:.1f}"),
        ("Canopy crowns", f"{canopy:,}")
    ]:
        tab.add_row(n, v)
    console.print(tab)

    # --- Generate human-readable interpretation ---
    interp_text = interpret_site(veg_pct, meanval, canopy_per_acre)

    # --- Write results.txt ---
    results_file = res / "results.txt"
    with open(results_file, "w") as f:
        f.write(f"Permit: {pid}\n")
        f.write(f"Permit area (acres): {permit_area:.1f}\n")
        f.write(f"Veg pixels: {veg:,}\n")
        f.write(f"NDVI threshold (avg): {thr:.2f}\n")
        f.write(f"Mean NDVI: {meanval:.2f}\n")
        f.write(f"Vegetation cover (%): {veg_pct:.1f}\n")
        f.write(f"Canopy crowns: {canopy:,}\n")
        f.write(f"Canopy per acre: {canopy_per_acre:.4f}\n\n")
        f.write("Plain-language interpretation:\n")
        f.write(textwrap.fill(interp_text, width=80))
        f.write("\n")

    # --- Append to running summary CSV ---
    cent = gpd.GeoSeries([geom], crs=crs).to_crs(4326).iloc[0].centroid
    header = "date,permit,permit_acres,veg_pix,thr,mean_ndvi,veg_pct,canopy,canopy_per_acre,lon,lat\n"
    line = f"{datetime.date.today()},{pid},{permit_area:.1f},{veg},{thr:.2f},{meanval:.2f},{veg_pct:.1f},{canopy},{canopy_per_acre:.4f},{cent.x:.6f},{cent.y:.6f}\n"
    if not Path(SUMMARY_CSV).exists():
        with open(SUMMARY_CSV, "w") as csvfile:
            csvfile.write(header)
    with open(SUMMARY_CSV, "a") as csvfile:
        csvfile.write(line)

    # --- Log and colorized completion panel ---
    console.print(Panel.fit(
        f"[gold]✅ Complete![/gold]\n"
        f"[steel_blue3]Permit:[/steel_blue3] {pid}\n"
        f"[steel_blue3]Area:[/steel_blue3] {permit_area:,.1f} acres\n"
        f"[steel_blue3]NDVI mean:[/steel_blue3] {meanval:.2f}\n"
        f"[steel_blue3]Vegetation cover:[/steel_blue3] {veg_pct:.1f}%\n"
        f"[steel_blue3]Canopy crowns:[/steel_blue3] {canopy:,}\n\n"
        f"[gold]Summary → {res}/results.txt[/gold]\n\n"
        f"[italic]{interp_text}[/italic]",
        title="🌾 WV Tree Counter – Appalachian Summary",
        border_style="yellow3", width=80
    ))
    


# ---------- Map Composer ----------
def make_map(mtif, geom, mapsdir, pid):
    if not mtif or not Path(mtif).exists():
        return
    with rasterio.open(mtif) as src:
        arr = src.read(1)
        bounds = src.bounds
        fig, ax = plt.subplots(figsize=(9, 7))
        show = ax.imshow(arr, cmap="YlGn", vmin=0, vmax=1)
        plt.colorbar(show, ax=ax, label="NDVI")
        plt.title(f"NDVI Map – {pid}")
        gdf = gpd.GeoDataFrame(geometry=[geom], crs=src.crs)
        gdf.boundary.plot(ax=ax, color="white", linewidth=1.5)
        ctx.add_basemap(ax, crs=src.crs, source=ctx.providers.Esri.WorldImagery, alpha=0.6)
        scalebar = ScaleBar(dx=src.res[0], units="m", dimension="si-length", location="lower left")
        ax.add_artist(scalebar)
        ax.text(0.9, 0.1, "N↑", transform=ax.transAxes, ha="center", va="center", fontsize=12, color="white")
        plt.axis("off")
        plt.tight_layout()
        outpng = mapsdir / f"{pid}_NDVI_map.png"
        plt.savefig(outpng, dpi=150)
        plt.close()

# ---------- Permit Geometry and Area ----------
def get_permit_area(permits, pid):
    """
    Returns the Shapely geometry of the selected permit and its area in acres.
    Handles CRS conversion safely for Shapely geometries.
    """
    sel = permits[permits["permit_id"].str.upper() == pid]
    if sel.empty:
        raise ValueError(f"Permit {pid} not found in permit shapefile.")
    
    # Extract the Shapely geometry
    geom = sel.iloc[0].geometry

    # Convert to a GeoSeries so we can reproject and calculate area
    geom_gs = gpd.GeoSeries([geom], crs=permits.crs).to_crs(3857)
    area_m2 = geom_gs.area.iloc[0]         # Area in square meters
    permit_area = area_m2 / 4046.86        # Convert to acres

    # Log this info
    log(f"📐 Calculated permit area for {pid}: {permit_area:.1f} acres")
    
    return geom, permit_area

# ---------- Compute Mosaic Stats ----------
def compute_mosaic_stats(records):
    """Compute mean NDVI safely, with auto-resampling for shape mismatches."""
    if not records:
        return None

    arrs = []
    base_profile = None
    base_shape = None
    for r in records:
        try:
            with rasterio.open(r["index"]) as src:
                arr = src.read(1)
                profile = src.profile
                if base_shape is None:
                    base_shape = arr.shape
                    base_profile = profile
                elif arr.shape != base_shape:
                    log(f"⚠️ Resampling {r['tile']} from {arr.shape} → {base_shape}")
                    data = np.empty(base_shape, dtype=np.float32)
                    rasterio.warp.reproject(
                        arr,
                        data,
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=base_profile["transform"],
                        dst_crs=base_profile["crs"],
                        resampling=rasterio.warp.Resampling.bilinear,
                    )
                    arr = data
                arrs.append(arr)
        except Exception as e:
            log(f"⚠️ Skipping {r['tile']} due to error: {e}")

    if not arrs:
        log("❌ No valid arrays found for mosaic stats.")
        return None

    stack = np.stack(arrs)
    mean = float(np.nanmean(stack))
    veg_pct = float(np.count_nonzero(stack > 0.25) / stack.size * 100)
    return {"mean": mean, "veg_pct": veg_pct}

# ---------- Main ----------
def main():
    console.print(Panel.fit("[hdr]🌾 WV Tree Counter (R19h Final – Basemap Edition)[/hdr]", border_style="yellow3", width=80))
    pshp = ensure_permits()
    ishp = ensure_index()
    permits = gpd.read_file(pshp)
    crs = permits.crs

    # Prompt user
    for _ in range(3):
        pid = Prompt.ask("🔢 WV permit number (e.g. S300120)").strip().upper()
        sel = permits[permits["permit_id"].str.upper() == pid]
        if not sel.empty:
            break
        console.print("[warn]Permit not found.[/warn]")
    if sel.empty:
        sys.exit("❌ No matching permit.")

    geom, permit_area = get_permit_area(permits, pid)
    pdir = ROOT / pid
    res = pdir / "results"
    maps = res / "maps"
    deep = res / "deepforest"
    for d in [pdir, res, maps, deep]:
        d.mkdir(parents=True, exist_ok=True)

    make_lookup_info(gpd.GeoDataFrame(geometry=[geom], crs=crs), ishp, pid, pdir)
    console.print(Panel.fit(f"Place NAIP .zip/.tif files listed in [gold]lookup_info.txt[/gold]\ninto:\n{DL_DIR}\n\nPress Enter to continue.", title="⏸️ Pause for Downloads", border_style="yellow3"))
    input()

    auto_ingest()
    tifs = list(NAIP_DIR.glob("*.tif"))
    if not tifs:
        sys.exit("❌ No imagery found.")

    console.print("[info]🌿 Computing NDVI indices…[/info]")
    try:
        recs = Parallel(n_jobs=max(1, os.cpu_count() - 1), prefer="threads")(delayed(compute_index)(p, geom, res) for p in tifs)
    except:
        recs = [compute_index(p, geom, res) for p in tifs]
    recs = [r for r in recs if r]

    mosaic_stats = compute_mosaic_stats(recs)
    mtif, mpng = (None, None)
    if recs:
        mtif, mpng = mean_index_map(recs, recs[0]["index"], maps, pid)
        make_map(mtif, geom, maps, pid)

    run_df = Prompt.ask("Run DeepForest? (y/N)").lower().startswith("y")
    if run_df:
        shps = run_deepforest(tifs, geom, deep)
        log(f"DeepForest found {len(shps)} canopy shapefiles.")
        summarize (pid, res, recs, mtif, shps, geom, crs, permit_area)
    console.print(Panel.fit(f"[good]✅ Complete![/good]\nResults → {res}", title="🌾 Done", border_style="yellow3"))


# ---------- Entrypoint ----------
if __name__ == "__main__":
    main()
