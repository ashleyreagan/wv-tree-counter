"""Regression tests for the permit mosaic, acreage, CRS, and DeepForest georeferencing fixes.

Uses small synthetic 4-band "NAIP-like" rasters, so no downloads are needed.
Run with:  pytest -q
"""
import os
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import pandas as pd
import geopandas as gpd
import pytest
import rasterio
from rasterio.transform import from_origin
from rasterio.enums import ColorInterp
from shapely.geometry import box, Polygon

from wv_tree_counter import main as wtc

UTM17 = "EPSG:26917"
UTM18 = "EPSG:26918"
VEG = (50, 80, 40, 200)      # NDVI ≈ 0.60
BARE = (120, 110, 100, 130)  # NDVI ≈ 0.04


def write_tile(path, crs, left, top, w, h, px=0.6, bands=VEG, alpha=False):
    arr = np.stack([np.full((h, w), b, np.uint8) for b in bands])
    with rasterio.open(path, "w", driver="GTiff", height=h, width=w, count=arr.shape[0],
                       dtype="uint8", crs=crs, transform=from_origin(left, top, px, px)) as d:
        d.write(arr)
        if alpha:
            d.colorinterp = [ColorInterp.red, ColorInterp.green, ColorInterp.blue, ColorInterp.alpha]
    return path


@pytest.fixture(autouse=True)
def _cwd(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # runlog.txt / summary.csv land in the temp dir


def test_permit_area_is_equal_area_not_web_mercator():
    sq = box(450000, 4190000, 451000, 4191000)  # exactly 1 km² ≈ 247.1 ac
    permits = gpd.GeoDataFrame({"permit_id": ["S000001"]}, geometry=[sq], crs=UTM17)
    _, acres = wtc.get_permit_area(permits, "S000001")
    assert acres == pytest.approx(247.1, abs=0.5)   # old code reported ~397


def test_multi_tile_permit_covers_both_tiles(tmp_path):
    # Two adjacent tiles: west is vegetated, east is bare. Permit spans both equally.
    a = write_tile(tmp_path / "a.tif", UTM17, 450000, 4190060, 100, 100, bands=VEG)
    b = write_tile(tmp_path / "b.tif", UTM17, 450060, 4190060, 100, 100, bands=BARE)
    permit = box(450000, 4190000, 450120, 4190060)
    tiles = wtc.tiles_for_permit([a, b], permit, UTM17)
    assert len(tiles) == 2
    arr, tf, valid, _, _ = wtc.build_mosaic(tiles, permit, UTM17, UTM17, tmp_path / "m.tif")
    assert arr.shape[2] == 200                       # full extent, not just tile A
    rec = wtc.compute_index(arr, valid, tf, UTM17, tmp_path, "P", nir_available=True)
    assert rec["veg_pct"] == pytest.approx(50.0, abs=1.0)
    assert rec["mean"] == pytest.approx((0.6 + 0.04) / 2, abs=0.02)   # old code: tile A diluted, B dropped


def test_overlapping_tiles_are_not_double_counted(tmp_path):
    # NAIP quarter-quads overlap their neighbours; overlap must be counted once.
    a = write_tile(tmp_path / "a.tif", UTM17, 450000, 4190060, 100, 100)
    b = write_tile(tmp_path / "b.tif", UTM17, 450030, 4190060, 100, 100)   # 50 % overlap
    permit = box(450000, 4190000, 450090, 4190060)
    arr, tf, valid, _, _ = wtc.build_mosaic([a, b], permit, UTM17, UTM17, tmp_path / "m.tif")
    assert valid.sum() == 150 * 100


def test_tiles_in_another_utm_zone_are_found_and_merged(tmp_path):
    # Straddle the zone 17/18 line at 78°W (WV Eastern Panhandle, PA, VA).
    line = gpd.GeoSeries.from_xy([-78.0], [39.5], crs="EPSG:4326")
    x17, y17 = line.to_crs(UTM17).iloc[0].coords[0]
    x18, y18 = line.to_crs(UTM18).iloc[0].coords[0]
    a = write_tile(tmp_path / "z17.tif", UTM17, x17 - 120, y17 + 60, 200, 200)   # west of line
    b = write_tile(tmp_path / "z18.tif", UTM18, x18 - 2, y18 + 60, 200, 200)     # east of line
    permit = box(x17 - 60, y17 - 30, x17 + 60, y17 + 30)  # permit in zone-17 coordinates
    tiles = wtc.tiles_for_permit([a, b], permit, UTM17)
    assert {t.name for t in tiles} == {"z17.tif", "z18.tif"}
    arr, tf, valid, _, _ = wtc.build_mosaic(tiles, permit, UTM17, UTM17, tmp_path / "m.tif")
    # 120 m x 60 m permit at 0.6 m ≈ 20,000 px; all of it should be covered
    assert valid.sum() == pytest.approx(200 * 100, rel=0.03)


def test_veg_cover_excludes_pixels_outside_permit_polygon(tmp_path):
    t = write_tile(tmp_path / "a.tif", UTM17, 450000, 4190060, 100, 100)
    tri = Polygon([(450000, 4190000), (450060, 4190000), (450000, 4190060)])  # half the bbox
    arr, tf, valid, _, _ = wtc.build_mosaic([t], tri, UTM17, UTM17, tmp_path / "m.tif")
    rec = wtc.compute_index(arr, valid, tf, UTM17, tmp_path, "P")
    assert rec["veg_pct"] == pytest.approx(100.0, abs=0.5)   # old code ≈ 50 %


def test_adaptive_threshold_pixel_count_matches_threshold(tmp_path):
    # NDVI ≈ 0.18: nothing above 0.25, everything above 0.15
    t = write_tile(tmp_path / "a.tif", UTM17, 450000, 4190060, 100, 100, bands=(100, 90, 80, 145))
    permit = box(450000, 4190000, 450060, 4190060)
    arr, tf, valid, _, _ = wtc.build_mosaic([t], permit, UTM17, UTM17, tmp_path / "m.tif")
    rec = wtc.compute_index(arr, valid, tf, UTM17, tmp_path, "P")
    assert rec["thr"] == 0.15
    assert rec["pix"] == 10000                      # old code reported 0 (count at 0.25)


def test_rgba_drone_ortho_uses_grvi_not_alpha_as_nir(tmp_path):
    t = write_tile(tmp_path / "ortho.tif", UTM17, 450000, 4190006, 100, 100, px=0.06,
                   bands=(50, 80, 40, 255), alpha=True)
    permit = box(450000, 4190000, 450006, 4190006)
    arr, tf, valid, nb, alpha = wtc.build_mosaic([t], permit, UTM17, UTM17, tmp_path / "m.tif")
    assert alpha and nb == 3
    rec = wtc.compute_index(arr, valid, tf, UTM17, tmp_path, "P", nir_available=not alpha)
    assert rec["type"] == "GRVI"


def test_deepforest_boxes_are_georeferenced():
    tf = from_origin(450000, 4190000, 0.6, 0.6)
    df = pd.DataFrame({"xmin": [0], "ymin": [0], "xmax": [10], "ymax": [10],
                       "label": ["Tree"], "score": [0.9]})
    g = wtc.boxes_to_geo(df, tf, UTM17)
    assert g.crs.to_epsg() == 26917
    assert tuple(round(v, 3) for v in g.geometry.iloc[0].bounds) == (450000, 4189994, 450006, 4190000)


def test_end_to_end_without_deepforest_writes_summary(tmp_path):
    a = write_tile(tmp_path / "a.tif", UTM17, 450000, 4190060, 100, 100)
    b = write_tile(tmp_path / "b.tif", UTM17, 450060, 4190060, 100, 100, bands=BARE)
    permit = box(450000, 4190000, 450120, 4190060)
    res = tmp_path / "res"; maps = res / "maps"; deep = res / "deepforest"
    for d in (res, maps, deep):
        d.mkdir(parents=True, exist_ok=True)
    out = wtc.analyze_permit("S000001", permit, UTM17, 1.78, [a, b], res, maps, deep,
                             run_df=False, basemap=False)
    assert out is not None and out["canopy"] == 0
    assert (res / "results.txt").exists()
    assert (tmp_path / "summary.csv").exists()
    assert (maps / "S000001_NDVI_map.png").exists()


def test_utm_zone_selection():
    pt = lambda lon, lat: gpd.GeoSeries.from_xy([lon], [lat], crs="EPSG:4326").iloc[0]
    assert wtc.utm_crs_for(pt(-81.8, 37.9), "EPSG:4326").to_epsg() == 26917   # southern WV
    assert wtc.utm_crs_for(pt(-77.9, 39.4), "EPSG:4326").to_epsg() == 26918   # Eastern Panhandle
    assert wtc.utm_crs_for(pt(-87.0, 33.5), "EPSG:4326").to_epsg() == 26916   # Alabama
