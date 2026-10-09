"""Multi-state tests: GeoMine permit lookup, NAIP STAC selection, batch CLI, and
summary.csv migration. Network calls are replaced by a fake HTTP session that
serves GeoJSON/STAC responses pointing at small local rasters."""
import os
os.environ.setdefault("MPLBACKEND", "Agg")

import json
import numpy as np
import pandas as pd
import geopandas as gpd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box, mapping

from wv_tree_counter import main as wtc
from wv_tree_counter import sources

UTM16 = "EPSG:26916"
LON, LAT = -84.2, 36.3                      # Cumberland Plateau, TN (UTM zone 16)


class Resp:
    def __init__(self, payload, status=200):
        self._p, self.status_code = payload, status

    def json(self):
        return self._p


class FakeHTTP:
    """Minimal stand-in for requests.Session used by sources.py."""

    def __init__(self, permits, items):
        self.permits = permits          # {(contact, PERMIT): GeoDataFrame in 4326}
        self.items = items              # STAC features
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append(("GET", url, params))
        if url == sources.SAS_TOKEN_URL:
            return Resp({"token": "sv=fake&sig=abc", "msft:expiry": "2099-01-01T00:00:00Z"})
        if url == sources.GEOMINE_URL:
            where = params["where"]
            feats = []
            for (contact, pid), gdf in self.permits.items():
                if f"contact = {contact}" in where and (f"'{pid}'" in where or ("LIKE" in where and pid[:3] in where)):
                    feats += json.loads(gdf.to_json())["features"]
            return Resp({"type": "FeatureCollection", "features": feats})
        raise AssertionError(f"unexpected GET {url}")

    def post(self, url, json=None, timeout=None):
        self.calls.append(("POST", url, json))
        assert url == sources.STAC_SEARCH_URL and json["collections"] == ["naip"]
        return Resp({"type": "FeatureCollection", "features": self.items, "links": []})


def utm_box(cx, cy, w, h):
    return box(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)


def write_tile(path, crs, bounds, px=0.6, bands=(50, 80, 40, 200)):
    left, bottom, right, top = bounds
    w, h = int(round((right - left) / px)), int(round((top - bottom) / px))
    arr = np.stack([np.full((h, w), b, np.uint8) for b in bands])
    with rasterio.open(path, "w", driver="GTiff", height=h, width=w, count=4, dtype="uint8",
                       crs=crs, transform=from_origin(left, top, px, px)) as d:
        d.write(arr)
    return str(path)


def stac_item(item_id, year, date, href, footprint_utm):
    geom = gpd.GeoSeries([footprint_utm], crs=UTM16).to_crs(4326).iloc[0]
    return {"type": "Feature", "id": item_id, "geometry": mapping(geom), "bbox": list(geom.bounds),
            "properties": {"naip:year": year, "datetime": f"{date}T00:00:00Z", "gsd": 0.6,
                           "naip:state": "tn", "proj:epsg": 26916},
            "assets": {"image": {"href": href}}}


@pytest.fixture
def tn_world(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(wtc, "ROOT", tmp_path / "data")
    for name in ["PERMIT_DIR", "NAIP_DIR", "DL_DIR", "INDEX_DIR"]:
        monkeypatch.setattr(wtc, name, tmp_path / "data" / name.lower())
    sources._token["value"] = None

    c = gpd.GeoSeries.from_xy([LON], [LAT], crs=4326).to_crs(UTM16).iloc[0]
    # Permit = two polygons (GeoMine often splits a permit), 120 m x 60 m overall
    parts = [utm_box(c.x - 30, c.y, 60, 60), utm_box(c.x + 30, c.y, 60, 60)]
    permit = gpd.GeoDataFrame({"permit_id": ["3270", "3270"], "permittee": ["Test Coal Co", "Test Coal Co"],
                               "contact": [4, 4]}, geometry=parts, crs=UTM16).to_crs(4326)

    # 2022: only the west half is imaged (partial)   2020: two tiles, full coverage
    west = (c.x - 100, c.y - 50, c.x, c.y + 50)
    east = (c.x, c.y - 50, c.x + 100, c.y + 50)
    t22 = write_tile(tmp_path / "m_2022_w.tif", UTM16, west)
    t20w = write_tile(tmp_path / "m_2020_w.tif", UTM16, west)
    t20e = write_tile(tmp_path / "m_2020_e.tif", UTM16, east, bands=(120, 110, 100, 130))
    items = [stac_item("tn_2022_w", "2022", "2022-06-10", t22, box(*west)),
             stac_item("tn_2020_w", "2020", "2020-07-01", t20w, box(*west)),
             stac_item("tn_2020_e", "2020", "2020-07-01", t20e, box(*east))]
    fake = FakeHTTP({(4, "3270"): permit}, items)
    monkeypatch.setattr(sources, "_http", fake)
    return {"fake": fake, "permit": permit, "tmp": tmp_path}


def test_geomine_lookup_dissolves_all_polygons(tn_world):
    geom, crs, acres = wtc.resolve_permit("TN", "3270")
    assert crs.to_epsg() == 4326
    assert acres == pytest.approx(120 * 60 / 4046.856, rel=0.01)   # both polygons counted
    where = tn_world["fake"].calls[0][2]["where"]
    assert "contact = 4" in where and "'3270'" in where


def test_geomine_missing_permit_reports_suggestions(tn_world):
    with pytest.raises(ValueError, match="not found in GeoMine for Tennessee"):
        wtc.resolve_permit("TN", "9999")


def test_maryland_requires_permit_file(tn_world, tmp_path):
    with pytest.raises(ValueError, match="--permit-file"):
        wtc.resolve_permit("MD", "SM-00-123")
    f = tmp_path / "md.gpkg"
    gpd.GeoDataFrame({"PERMIT_NO": ["SM-00-123"]}, geometry=[box(-79.1, 39.55, -79.099, 39.551)],
                     crs=4326).to_file(f)
    geom, crs, acres = wtc.resolve_permit("MD", "sm-00-123", permit_file=str(f), id_field="PERMIT_NO")
    assert acres > 0


def test_stac_picks_newest_year_with_full_coverage(tn_world):
    geom, crs, _ = wtc.resolve_permit("TN", "3270")
    info = sources.naip_for_permit(wtc.to_crs(geom, crs, "EPSG:4326"))
    assert info["year"] == "2020"                      # 2022 only covers half the permit
    assert len(info["hrefs"]) == 2
    assert info["available_years"] == ["2022", "2020"]


def test_stac_respects_requested_year(tn_world):
    geom, crs, _ = wtc.resolve_permit("TN", "3270")
    info = sources.naip_for_permit(wtc.to_crs(geom, crs, "EPSG:4326"), year="2022")
    assert info["year"] == "2022" and len(info["hrefs"]) == 1
    with pytest.raises(sources.SourceError, match="Available"):
        sources.naip_for_permit(wtc.to_crs(geom, crs, "EPSG:4326"), year="2016")


def test_sign_appends_token_only_to_urls(tn_world):
    assert sources.sign("/local/file.tif") == "/local/file.tif"
    assert sources.sign("https://naipeuwest.blob.core.windows.net/naip/x.tif") \
        == "https://naipeuwest.blob.core.windows.net/naip/x.tif?sv=fake&sig=abc"


def test_off_season_dates_are_flagged(tn_world):
    for it in tn_world["fake"].items:
        it["properties"]["datetime"] = "2020-12-11T00:00:00Z"
    geom, crs, _ = wtc.resolve_permit("TN", "3270")
    info = sources.naip_for_permit(wtc.to_crs(geom, crs, "EPSG:4326"))
    assert info["off_season_dates"] == ["2020-12-11"]


def test_cli_end_to_end_tennessee_zone16(tn_world):
    tmp = tn_world["tmp"]
    out = wtc.run_cli(["--state", "TN", "--permit", "3270", "--no-basemap"])
    assert len(out) == 1
    assert out[0]["veg_pct"] == pytest.approx(50.0, abs=2)   # west half vegetated, east bare
    res = tmp / "data" / "TN" / "3270" / "results"
    assert (res / "results.txt").exists() and "State: TN" in (res / "results.txt").read_text()
    with rasterio.open(res / "3270_mosaic.tif") as m:
        assert m.crs.to_epsg() == 26916                      # permit's own UTM zone
    df = pd.read_csv(tmp / "summary.csv", dtype=str)
    assert list(df.columns) == wtc.SUMMARY_COLUMNS
    assert df.iloc[0][["state", "imagery_source", "imagery_year"]].tolist() == ["TN", "naip-stac", "2020"]


def test_batch_continues_after_a_failure(tn_world, tmp_path):
    lst = tmp_path / "permits.txt"
    lst.write_text("# test batch\nTN,3270\nTN,0000\n")
    out = wtc.run_cli(["--permit-list", str(lst), "--no-basemap"])
    assert len(out) == 1   # 0000 fails, 3270 still runs


def test_permit_list_parsing(tmp_path):
    f = tmp_path / "p.txt"
    f.write_text("S300120\nPA, 56110101  # comment\n\nky,8360123\n")
    assert wtc.read_permit_list(f, "WV") == [("WV", "S300120"), ("PA", "56110101"), ("KY", "8360123")]
    with pytest.raises(ValueError):
        wtc.read_permit_list(f, None)


def test_old_wv_summary_csv_is_upgraded(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    p = tmp_path / "summary.csv"
    p.write_text("date,permit,permit_acres,veg_pix,thr,mean_ndvi,veg_pct,canopy,canopy_per_acre,lon,lat\n"
                 "2025-10-20,S300120,397.3,100,0.25,0.40,30.0,0,0.0000,-81.8,37.9\n")
    wtc.append_summary({"date": "2026-10-08", "state": "PA", "permit": "X1", "permit_acres": "10.0"}, p)
    df = pd.read_csv(p, dtype=str)
    assert list(df.columns) == wtc.SUMMARY_COLUMNS
    assert df["state"].tolist() == ["WV", "PA"]
    assert df.loc[0, "permit_acres"] == "397.3"
