"""Network data sources: permit boundaries (OSMRE GeoMine) and NAIP imagery
(Microsoft Planetary Computer STAC). Kept separate from main.py so they can be
mocked in tests and swapped out later (e.g. AWS, a DOI image server)."""
import os
import json
import time
from collections import defaultdict

import geopandas as gpd
from shapely.geometry import shape, mapping

# Faster, quieter remote COG reads (GDAL /vsicurl)
os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
os.environ.setdefault("GDAL_HTTP_MAX_RETRY", "3")
os.environ.setdefault("GDAL_HTTP_RETRY_DELAY", "2")

GEOMINE_URL = ("https://geoservices.osmre.gov/arcgis/rest/services/GeoMine/"
               "AllCoalmineOperations/MapServer/0/query")
STAC_SEARCH_URL = "https://planetarycomputer.microsoft.com/api/stac/v1/search"
SAS_TOKEN_URL = "https://planetarycomputer.microsoft.com/api/sas/v1/token/naip"

LEAF_ON_MONTHS = range(5, 10)   # May–September

_http = None


def http():
    """Shared requests session (replaced with a fake in tests)."""
    global _http
    if _http is None:
        import requests
        _http = requests.Session()
    return _http


class SourceError(RuntimeError):
    pass


def union(geoseries):
    """Dissolve a GeoSeries to one geometry (geopandas 0.x and 1.x)."""
    return geoseries.union_all() if hasattr(geoseries, "union_all") else geoseries.unary_union


# ---------------------------------------------------------------------------
# Permits – OSMRE GeoMine
# ---------------------------------------------------------------------------
def _sql_str(s):
    return "'" + str(s).replace("'", "''") + "'"


def geomine_query(where, out_fields="*", timeout=120):
    params = {"where": where, "outFields": out_fields, "outSR": 4326,
              "returnGeometry": "true", "f": "geojson"}
    try:
        r = http().get(GEOMINE_URL, params=params, timeout=timeout)
    except Exception as e:
        raise SourceError(
            f"Could not reach OSMRE GeoMine ({e}). The service may only be reachable "
            "on the DOI network/VPN; use --permit-file to supply boundaries instead.") from e
    if r.status_code != 200:
        raise SourceError(f"GeoMine returned HTTP {r.status_code}")
    data = r.json()
    if "error" in data:
        raise SourceError(f"GeoMine error: {data['error']}")
    feats = data.get("features", [])
    if not feats:
        return gpd.GeoDataFrame(geometry=[], crs="EPSG:4326")
    return gpd.GeoDataFrame.from_features(feats, crs="EPSG:4326")


def fetch_permit_geomine(contact_code, permit_id):
    """All GeoMine polygons for one permit from one regulatory authority.
    Returns (GeoDataFrame, suggestions). Suggestions are near-matches when the
    exact ID isn't found (permit ID formats differ by state)."""
    pid = permit_id.strip().upper()
    gdf = geomine_query(f"UPPER(permit_id) = {_sql_str(pid)} AND contact = {int(contact_code)}")
    if len(gdf):
        return gdf, []
    near = geomine_query(
        f"UPPER(permit_id) LIKE {_sql_str('%' + pid + '%')} AND contact = {int(contact_code)}",
        out_fields="permit_id,mine_name,permittee")
    sugg = sorted(set(near["permit_id"].dropna().astype(str))) if len(near) else []
    return gdf, sugg[:15]


# ---------------------------------------------------------------------------
# Permits – local file (any state, offline, or Maryland)
# ---------------------------------------------------------------------------
def load_permit_file(path, id_field, permit_id):
    gdf = gpd.read_file(path)
    if id_field not in gdf.columns:
        raise SourceError(f"Field '{id_field}' not in {path}. Fields: {', '.join(map(str, gdf.columns))}")
    sel = gdf[gdf[id_field].astype(str).str.strip().str.upper() == permit_id.strip().upper()]
    return sel


# ---------------------------------------------------------------------------
# Imagery – NAIP on Microsoft Planetary Computer
# ---------------------------------------------------------------------------
_token = {"value": None, "expires": 0}


def sign(href):
    """Append the Planetary Computer SAS token to a blob-storage URL."""
    if not str(href).startswith("http"):
        return href
    if _token["value"] is None or time.time() > _token["expires"] - 300:
        r = http().get(SAS_TOKEN_URL, timeout=60)
        if r.status_code != 200:
            raise SourceError(f"Could not get NAIP access token (HTTP {r.status_code})")
        j = r.json()
        _token["value"] = j["token"]
        _token["expires"] = time.time() + 45 * 60
    sep = "&" if "?" in href else "?"
    return f"{href}{sep}{_token['value']}"


def stac_search(geom_4326, limit=100):
    body = {"collections": ["naip"], "intersects": mapping(geom_4326), "limit": limit}
    feats, url, method = [], STAC_SEARCH_URL, "POST"
    for _ in range(50):  # pagination guard
        try:
            if method == "POST":
                r = http().post(url, json=body, timeout=120)
            else:
                r = http().get(url, timeout=120)
        except Exception as e:
            raise SourceError(f"Could not reach the NAIP catalog ({e}). Use --imagery local.") from e
        if r.status_code != 200:
            raise SourceError(f"NAIP catalog returned HTTP {r.status_code}")
        page = r.json()
        feats.extend(page.get("features", []))
        nxt = next((l for l in page.get("links", []) if l.get("rel") == "next"), None)
        if not nxt:
            break
        url = nxt["href"]
        method = nxt.get("method", "GET").upper()
        body = nxt.get("body", body)
    return feats


def choose_naip_year(items, geom_4326, year=None):
    """Pick ONE NAIP year whose tiles fully cover the permit (newest first, or the
    requested year). Returns (year, [items])."""
    by_year = defaultdict(list)
    for it in items:
        y = str(it["properties"].get("naip:year") or it["properties"]["datetime"][:4])
        by_year[y].append(it)
    if not by_year:
        raise SourceError("No NAIP imagery found for this permit.")
    years = [str(year)] if year else sorted(by_year, reverse=True)
    for y in years:
        if y not in by_year:
            raise SourceError(f"No NAIP for {year} here. Available: {', '.join(sorted(by_year, reverse=True))}")
        footprint = union(gpd.GeoSeries([shape(it["geometry"]) for it in by_year[y]], crs="EPSG:4326"))
        covered = footprint.buffer(1e-6).covers(geom_4326)
        if covered or year:
            return y, by_year[y]
    # nothing fully covers – fall back to the newest year
    y = sorted(by_year, reverse=True)[0]
    return y, by_year[y]


def naip_for_permit(geom_4326, year=None):
    """Signed COG URLs for one NAIP year covering the permit, plus metadata."""
    items = stac_search(geom_4326)
    y, chosen = choose_naip_year(items, geom_4326, year)
    dates = sorted({it["properties"]["datetime"][:10] for it in chosen})
    gsds = sorted({float(it["properties"].get("gsd", 0)) for it in chosen})
    hrefs = [sign(it["assets"]["image"]["href"]) for it in chosen]
    off_season = [d for d in dates if int(d[5:7]) not in LEAF_ON_MONTHS]
    return {
        "year": y,
        "hrefs": hrefs,
        "ids": [it["id"] for it in chosen],
        "dates": dates,
        "gsd": gsds,
        "states": sorted({it["properties"].get("naip:state", "") for it in chosen}),
        "available_years": sorted({str(it["properties"].get("naip:year")) for it in items}, reverse=True),
        "off_season_dates": off_season,
    }
