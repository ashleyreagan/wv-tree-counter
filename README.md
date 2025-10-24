# 🌾 WV Tree Counter – R19h (Basemap Edition)
**Author:** Ashley Mitchell, OSMRE • **Year:** 2025

A Python workflow to assess vegetation recovery on West Virginia mining permits using:
- **WVDEP permit boundaries**
- **NAIP 4-band imagery (RGB + NIR)**
- **NDVI/GRVI** indices (adaptive thresholding)
- **DeepForest** canopy crown detection (optional)
- Permit-clipped mosaics, maps (with **legend**, **scalebar**, **north arrow**)
- Human-readable summaries + a running **summary.csv** across permits

---

## ✨ What it does (end-to-end)
1. **Loads WVDEP permits** (auto-downloads if missing).
2. Prompts for a **permit number** (e.g., `S300120`) and computes **true GIS acres**.
3. Loads the **WV NAIP 2022+ tile index**, intersects with the permit, and writes:
   - **Centroid** (WGS84 + UTM 17N)
   - **Suggested NAIP tile names**
   - → `data/<PERMIT_ID>/lookup_info.txt`
4. **Pauses** so you can download NAIP tiles; place **`.tif` or `.zip`** into:
   - `data/naip_files/downloads/`
   - The tool **auto-ingests** (unzips, moves `.tif/.tfw/.xml/.ovr` into `data/naip_files/`).
5. Computes **NDVI** (or **GRVI** fallback if no NIR/band issues), **clipped to the permit**.
   - **Adaptive NDVI threshold**: starts at **0.25**, lowers to **0.20/0.15/0.10** only if there are no returns.
   - Skips any pre-computed `*_NDVI.tif` in the imagery folder.
6. (Optional) Runs **DeepForest** on **RGB** data **inside the permit**.
   - Auto-downloads the **default canopy model** on first run if missing.
7. Builds per-tile products **and** a permit-level **mosaic** and **maps** with:
   - Classic vegetation ramp (NDVI)
   - Canopy overlays (DeepForest)
   - **Permit boundary**, **north arrow**, **scalebar**
8. Writes a plain-English summary and appends **summary.csv** for cross-permit rollups.

## 📁 Project layout
```
WestVirginia/
├─ WV_TreeCounter_R19h_Final_Basemap.py      # main script
├─ runlog.txt                                # global log (also per-permit notes)
├─ summary.csv                               # running cross-permit summary
└─ data/
├─ permits/                                  # WVDEP shapefile (auto-downloaded)
├─ naip_files/                               # final NAIP .tif, .tfw, .ovr, etc (auto-ingested)
│  └─ downloads/                             # drop .zip or .tif here
├─ naip_index/ | index/                      # WV NAIP index (auto or manual)
├─ <PERMIT_ID>/
│  ├─ lookup_info.txt                        # centroid + suggested NAIP tiles
│  └─ results/
│     ├─ maps/                               # NDVI/Canopy PNGs with legend/scalebar
│     ├─ deepforest/                      	 # DF per-tile shapefiles + merged
│     ├─ *_NDVI.tif / *_NDVI.png           	 # per-tile clipped rasters
│     ├─ NDVI_mean.tif/.png       			 # permit mosaic + render
│     ├─ canopy_merged.shp        			 # merged canopy detections
│     ├─ results.txt                         # plain-English results
│     └─ runlog.txt               			 # permit-specific log (start/end stamps)
└─ results/ (legacy, if created by older runs)
```
---

## 🗂 Data sources
- **Permits:** WVDEP Mining & Reclamation Permit Boundary shapefile  
  `https://tagis.dep.wv.gov/data/vector/SHP_WVDEP_GIS_data_mining_reclamation_permit_boundary.zip`
- **NAIP tile index (WV, 2022+):**  
  `https://www.fpacbc.usda.gov/sites/default/files/2024-10/wv_naip22qq.zip`  
  (Place the shapefile in `data/naip_index/` or let the script download.)

---

## 🛠 Installation (Mac M2 friendly)
We strongly recommend a fresh **conda** env.

```bash
# From repo root (WestVirginia/)
conda create -n smcra_wv -c conda-forge python=3.10 gdal geos proj
conda activate smcra_wv

# Core deps
pip install rich tqdm deepforest geopandas rasterio shapely numpy matplotlib pandas joblib requests matplotlib-scalebar contextily

DeepForest on Apple Silicon: It uses PyTorch with Metal (MPS). If MPS behaves oddly, it falls back to CPU automatically in the script.

⸻

▶️ Running the tool

conda activate smcra_wv
python WV_TreeCounter_R19h_Final_Basemap.py

You’ll be prompted for a permit ID. The flow is:
	1.	Tool prints the centroid and suggested NAIP tiles, then pauses.
	2.	Download those tiles (raw 4-band NAIP .tif).
	3.	Drop them into: data/naip_files/downloads/
	4.	Press Enter — the tool auto-ingests (unzips/moves) and continues.

You can re-run the same permit any time; the script skips work it already did and only processes what’s missing.

⸻

🧮 How NDVI/GRVI are computed
	•	Band auto-detect: Uses GDAL/Rasterio metadata to infer bands.
	•	Preferred: Red = band 1, NIR = band 4 (typical NAIP).
	•	If NIR missing/unreadable, falls back to GRVI = (G–R)/(G+R).
	•	Adaptive threshold:
	•	Start at 0.25; if no vegetation pixels, it tries 0.20 → 0.15 → 0.10.
	•	The chosen threshold is reported per tile and in the summary.
	•	Clipping & mosaics: All indices are clipped to the permit; mosaics resample tiles to a common grid automatically.

⸻

🌳 DeepForest (optional)
	•	Runs inside the permit only, on RGB (bands 1–3).
	•	Auto-downloads the default canopy model on first run if missing.
	•	Outputs:
	•	per-tile *_canopy.shp
	•	permit-level <PERMIT>_canopy_merged.shp
	•	Summary shows total crowns and crowns per acre.

⸻

🗺 Map outputs
	•	NDVI (classic vegetation ramp) and Canopy maps:
	•	Permit outline
	•	North arrow & scalebar
	•	Legend
	•	Saved under data/<PERMIT>/results/maps/.

⸻

📄 Outputs & logs
	•	results.txt: Plain-English interpretation + statistics (NDVI mean, veg cover %, crowns, acres).
	•	summary.csv: Appended on every run → quick cross-permit rollups.
	•	lookup_info.txt: Centroid(s) + tile names (for manual imagery download).
	•	runlog.txt (root) and runlog_<PERMIT>.txt (per-permit): detailed progress with start/end stamps.

⸻

🧪 Quick sanity checks

Use rio info (rasterio) or gdalinfo to confirm your tiles are 4-band NAIP:

rio info data/naip_files/<tile>.tif
# Look for: "count": 4 and "colorinterp": ["red","green","blue","undefined"]

If you see 1 band it’s likely an NDVI you previously produced — move it out of naip_files/ so the tool doesn’t try to index it as raw data.

⸻

🧯 Troubleshooting

“tuple index out of range” / “Unsupported band count: 1”
	•	Cause: You placed *_NDVI.tif in data/naip_files/ or a non-4-band tile.
	•	Fix: Move NDVI products elsewhere. Keep only raw 4-band NAIP tiles in naip_files/.

“all input arrays must have the same shape” (mosaic)
	•	Cause: Tiles have different shapes/resolutions.
	•	Fix: The script now resamples on the fly. If you see it again, remove partial leftovers and re-run.

Very low crowns / low NDVI on mine lands
	•	Often normal for recently disturbed/reclaimed areas. Use the GRVI fallback and read the results.txt interpretation.

zsh: bad pattern when pasting log lines
	•	Your shell treats [...] as globs. Wrap the log snippet in quotes or just open runlog.txt.

Indentation/TabError
	•	If you edited the script, run a formatter (e.g., black) or make sure your editor uses spaces consistently.

⸻

🔍 Examples

Run for biggest permits (one-liner in Python):

python - <<'EOF'
import geopandas as gpd
g = gpd.read_file("data/permits/WVDEP_GIS_data_mining_reclamation_permit_boundary.shp")
g["acres_geom"] = g.geometry.to_crs(26917).area/4046.86
print(g.nlargest(2,"acres_geom")[["permit_id","operator","acres_geom"]].to_string(index=False))
EOF

Site you tested: S501397 centroid ≈ 37.876980, -81.799671 (WGS84) — the tool writes this to lookup_info.txt automatically.

⸻

✅ GitHub hygiene
	•	.gitignore excludes large geodata and outputs (already included).
	•	Include this README, your script, and optionally:
	•	requirements.txt (and/or a environment.yml)
	•	LICENSE (MIT)
	•	CHANGELOG.md for version bumps (R19h etc.)

⸻

📜 License

MIT — © 2025 Ashley Mitchell

---
