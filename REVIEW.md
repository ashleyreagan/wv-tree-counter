# Codebase Review: WV Tree Counter (R19h – Basemap Edition)

**Reviewer:** Claude (automated review)
**Date:** 2026-02-09
**Version reviewed:** 0.19.1

## Project Summary

This is a Python-based geospatial analysis tool for assessing vegetation recovery on West Virginia SMCRA mining permits. It processes NAIP aerial imagery to compute NDVI/GRVI vegetation indices and optionally runs DeepForest-based canopy detection, producing maps, statistics, and plain-English summaries.

**Author:** Ashley Mitchell, U.S. Dept. of the Interior, OSMRE
**Version:** 0.19.1 | **License:** MIT | **Python:** 3.9+
**Structure:** Single-module package (`wv_tree_counter/main.py`, ~660 lines)

---

## Critical Bugs

### 1. `summarize()` is only called when DeepForest is enabled (`main.py:650-654`)

```python
run_df = Prompt.ask("Run DeepForest? (y/N)").lower().startswith("y")
if run_df:
    shps = run_deepforest(tifs, geom, deep)
    log(f"DeepForest found {len(shps)} canopy shapefiles.")
    summarize(pid, res, recs, mtif, shps, geom, crs, permit_area)
```

The `summarize()` call is nested inside the `if run_df:` block. If the user declines DeepForest (which is the default -- "y/N"), **no summary is generated** -- no `results.txt`, no `summary.csv` entry, no Rich console table. This is almost certainly unintended. The `summarize()` call should be outside the conditional, passing an empty list for `shps` when DeepForest is skipped.

### 2. Dead code after early `return` in `mean_index_map()` (`main.py:363-374`)

The function returns at line 363, making lines 364-374 (a duplicate PNG generation block) completely unreachable:

```python
    plt.savefig(png, dpi=150)
    plt.close()
    return mtif, png       # <-- returns here (line 363)
    # Generate PNG with legend and boundary   <-- dead code starts
    png = mapsdir / f"{pid}_index_mean.png"
    ...
    return mtif, png        # never reached
```

This appears to be a leftover from a refactor. The dead block should be removed.

### 3. `pix` count returned with wrong threshold (`main.py:277-301`)

```python
thr = 0.25
pix = np.count_nonzero(index > thr)  # counted at 0.25
if pix < 500:
    for t in [0.20, 0.15, 0.10]:
        if np.count_nonzero(index > t) > 500:
            thr = t                    # threshold lowered
            break
# ... but `pix` is still the count at 0.25, not at the new threshold
return {"tile": tif.name, "thr": thr, "pix": int(pix), ...}
```

When adaptive thresholding kicks in, `thr` is updated but `pix` still holds the count from the original 0.25 threshold. The downstream `summarize()` function sums `pix` from all records to report total vegetation pixels, which will be incorrect for tiles where the threshold was lowered.

---

## Moderate Issues

### 4. Bare `except` clause (`main.py:640`)

```python
try:
    recs = Parallel(...)(...)
except:
    recs = [compute_index(p, geom, res) for p in tifs]
```

Bare `except:` catches everything including `KeyboardInterrupt` and `SystemExit`. This should be `except Exception:` at minimum.

### 5. `compute_mosaic_stats()` result is unused (`main.py:644`)

```python
mosaic_stats = compute_mosaic_stats(recs)
```

The return value is computed (involving full raster I/O and reprojection) but never referenced again. This is wasted work.

### 6. Duplicate header block / shebang (`main.py:1 and 62`)

The file has two `#!/usr/bin/env python3` lines and two separate header comment blocks. Only the first shebang has any effect; the second is noise.

### 7. Missing `__init__.py` (`wv_tree_counter/` directory)

The package directory lacks an `__init__.py`. While `find:` in setup.cfg may still discover it, this can cause import issues in some environments and deviates from standard Python packaging practice.

### 8. `.gitignore` is minimal

The `.gitignore` only excludes `.DS_Store`. It should also exclude:
- `data/` (large GeoTIFFs, shapefiles, downloads)
- `__pycache__/`, `*.pyc`, `*.pyo`
- `*.egg-info/`, `dist/`, `build/`
- `runlog.txt`, `summary.csv`
- `.env`, `*.whl`

This risks accidentally committing large geospatial data files.

### 9. `ensure_index()` unzips unconditionally after download failure (`main.py:185-188`)

```python
if not download(NAIP_INDEX_URL, z):
    console.print(...)
    input("Place the zip into index folder and press Enter...")
unzip(z, INDEX_DIR)  # will crash if user didn't place the file
```

If the download fails and the user presses Enter without placing the file, `unzip()` will throw a `FileNotFoundError`. Contrast with `ensure_permits()` which properly exits on download failure.

### 10. Relative `ROOT` path (`main.py:123`)

```python
ROOT = Path("data")
```

This is relative to `os.getcwd()`, not the script location. If the user invokes `wv-tree-counter` from a different directory, all data will be written to an unexpected location.

---

## Code Quality Observations

- **Auto-install at import time (`main.py:88-97`)** -- Running `pip install` at import time is surprising behavior and a potential security concern. Should be handled during setup/installation.
- **No type hints** -- The entire codebase lacks type annotations.
- **No tests** -- No test files, no pytest configuration, and no CI/CD pipeline.
- **Global warning suppression (`main.py:110`)** -- `warnings.filterwarnings("ignore", category=UserWarning)` blanket-suppresses all UserWarnings.
- **Hardcoded CRS values** -- EPSG codes (26917, 3857, 4326) are scattered without named constants.
- **PEP 8 violation (`main.py:431`)** -- `def summarize (pid, ...)` has a space before the parenthesis.

---

## Architecture Assessment

### Strengths
- Clear linear workflow that's easy to follow
- Adaptive NDVI thresholding is a thoughtful domain-specific feature
- Plain-English interpretation function adds real user value
- Good use of Rich for terminal UI
- Auto-downloading of reference data reduces friction

### Areas for Improvement
- Single 660-line file would benefit from modular decomposition (e.g., `io.py`, `indices.py`, `visualization.py`, `summary.py`)
- No separation of concerns between I/O, computation, and presentation
- Functions like `mean_index_map()` both compute and visualize, violating single-responsibility
- No configuration file support -- all paths and URLs are hardcoded

---

## Priority Recommendations

| Priority | Issue | Impact |
|----------|-------|--------|
| **P0** | `summarize()` only called with DeepForest enabled | No output for typical runs |
| **P0** | Dead code in `mean_index_map()` | Confusion, possible lost feature |
| **P0** | `pix` count uses wrong threshold | Incorrect vegetation statistics |
| **P1** | Bare `except` catches `KeyboardInterrupt` | Can't Ctrl+C during parallel compute |
| **P1** | Missing `__init__.py` | Package import issues |
| **P1** | `.gitignore` incomplete | Risk of committing large data files |
| **P2** | Unused `mosaic_stats` result | Wasted computation |
| **P2** | Relative `ROOT` path | Wrong data directory if invoked elsewhere |
| **P2** | No tests or CI | No regression safety net |
