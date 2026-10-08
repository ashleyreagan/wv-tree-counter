"""State registry for Appalachian coal country.

Each entry says where permit boundaries come from:
  - "wvdep"   : WVDEP TAGIS permit shapefile (the original WV source)
  - "geomine" : OSMRE GeoMine "Surface Coalmine Boundary" service, filtered by the
                regulatory authority's GeoMine contact code
  - "file"    : a permit boundary file you supply (--permit-file / --id-field)

GeoMine contact codes come from the layer's "Contact" coded-value domain:
https://geoservices.osmre.gov/arcgis/rest/services/GeoMine/AllCoalmineOperations/MapServer/0
"""

STATES = {
    "WV": {
        "name": "West Virginia",
        "regulator": "WVDEP Division of Mining and Reclamation",
        "permit_source": "wvdep",
        "geomine_contact": 2,
        "coal_fields": "Northern and Central Appalachian",
        "example_permit": "S300120",
    },
    "PA": {
        "name": "Pennsylvania",
        "regulator": "PA DEP Bureau of Mining Programs",
        "permit_source": "geomine",
        "geomine_contact": 20,
        "coal_fields": "Northern Appalachian (bituminous) and Anthracite",
    },
    "OH": {
        "name": "Ohio",
        "regulator": "ODNR Division of Mineral Resources Management",
        "permit_source": "geomine",
        "geomine_contact": 21,
        "coal_fields": "Northern Appalachian",
    },
    "MD": {
        "name": "Maryland",
        "regulator": "MDE Mining Program",
        # Maryland does not appear in GeoMine's contact list; supply MDE's permit layer.
        "permit_source": "file",
        "geomine_contact": None,
        "coal_fields": "Northern Appalachian (Garrett and Allegany counties)",
    },
    "VA": {
        "name": "Virginia",
        "regulator": "Virginia Energy (formerly DMME), Big Stone Gap",
        "permit_source": "geomine",
        "geomine_contact": 1,
        "coal_fields": "Central Appalachian",
    },
    "KY": {
        "name": "Kentucky",
        "regulator": "KY Energy and Environment Cabinet, Division of Mine Permits",
        "permit_source": "geomine",
        "geomine_contact": 3,
        "coal_fields": "Central Appalachian (eastern KY); Illinois Basin (western KY)",
    },
    "TN": {
        "name": "Tennessee",
        "regulator": "OSMRE Knoxville Field Office (federal program state)",
        "permit_source": "geomine",
        "geomine_contact": 4,
        "coal_fields": "Central Appalachian (Cumberland Plateau)",
    },
    "AL": {
        "name": "Alabama",
        "regulator": "Alabama Surface Mining Commission",
        "permit_source": "geomine",
        "geomine_contact": 8,
        "coal_fields": "Southern Appalachian (Warrior, Cahaba, Coosa fields)",
    },
}


def get_state(code):
    code = (code or "").strip().upper()
    if code not in STATES:
        raise ValueError(f"Unknown state '{code}'. Supported: {', '.join(STATES)}")
    return code, STATES[code]
