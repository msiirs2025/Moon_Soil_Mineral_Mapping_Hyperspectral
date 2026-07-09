# Runbook

How to run, rerun, and troubleshoot the pipeline.

## 0. Prerequisites

- **Python:** `C:\Users\HP\anaconda3\python.exe` — must have `numpy`, `scipy`, `rasterio`
  (bundled GDAL 3.10). There is **no** `osgeo` binding; do not `import osgeo`.
- **QGIS 3.40** at `C:\Program Files\QGIS 3.40.12\bin` — provides `gdalwarp.exe` and
  `gdal_translate.exe` used by the CH1 warp. Paths are hardcoded in `scripts/common.py`
  (`QGIS_BIN`) and `scripts/ch1_pipeline.py` (`GDAL_DATA` / `PROJ_LIB` / `PATH`).
  If QGIS moved, update `QGIS_BIN` in `common.py`.
- **Disk:** outputs total ~87 GB. Ensure free space on `D:` before a full rerun.
- **Raw inputs must be present** (already in place):
  - CH2: `ch2\...\_d_img_d18.qub` / `.xml` / `.spm` (see `IIRS_QUB` in `common.py`).
  - CH1: `ch1\cartOrder (1)\cartorder\m3g20090731t045352_v01_rfl.img` + `_v03_loc.img`.

All commands below assume this working directory:

```powershell
cd "D:\new project moon\pipeline\scripts"
$py = "C:\Users\HP\anaconda3\python.exe"
```

## 1. Run CH2 (Chandrayaan‑2 IIRS) — full chain

Steps are ordered and each reads the previous step's `.img`. Run in sequence:

```powershell
& $py ch2_pipeline.py 1          # step1 radiance (256 bands)      -> outputs/ch2/step1_radiance
& $py ch2_pipeline.py angle      # incidence-angle aux + angle.npy -> outputs/ch2/aux
& $py ch2_pipeline.py 2          # thermal removal -> reflectance  -> outputs/ch2/step2_thermal + aux/Temp
& $py ch2_pipeline.py 3          # photometric correction          -> outputs/ch2/step3_photometric
& $py ch2_seleno_destripe.py 4   # selenoreference (Moon eqc)      -> outputs/ch2/step4_seleno
& $py ch2_seleno_destripe.py 5   # destripe + selenoreference      -> outputs/ch2/step5_destripe (FINAL)
```

Shortcuts: `ch2_pipeline.py all` runs 1+angle+2+3; `ch2_seleno_destripe.py all` runs 4+5.
Note `step3` requires `aux/angle.npy` (created by the `angle` step) — run `angle` before `3`.

**If step 5 georeferencing was interrupted** (the sensor `.img` exists but the
selenoref one is partial):

```powershell
& $py finish_step5.py            # deletes the partial and re-warps the destriped cube
```

## 2. Run CH1 (Chandrayaan‑1 M3) — selenoreference

```powershell
& $py ch1_pipeline.py            # writes sensor/lon/lat rasters (reused if present)
& $py ch1_fix_warp.py            # GCP + TPS warp -> FINAL product in outputs/ch1/seleno
```

`ch1_fix_warp.py` is the definitive step and reuses the sensor raster written by
`ch1_pipeline.py`. The `-geoloc` warp inside `ch1_pipeline.py` was superseded (it
collapsed longitude) — `ch1_fix_warp.py` overwrites its output. If the sensor rasters
already exist, `ch1_pipeline.py` reuses them instead of rewriting.

## 3. Validate against ISRO reference products

```powershell
& $py validate_angle.py          # incidence angle vs ISRO _angle  (expect near-exact)
& $py validate_step2.py          # temperature vs _Temp, reflectance vs _inc_corrRef
```

Expected: incidence mean|Δ| ~0; temperature mean|Δ| ≈ 0.1 K; reflectance median rel < ~1%.

## 4. Inspect a raw cube

```powershell
& $py inspect_raw.py             # prints raw header/dimension info
```

## 5. Open the products

- **ENVI:** open the `.hdr` (or `.img`) directly. Wavelengths and, for georeferenced
  steps, `map info` are embedded. CH1 nodata = −999; CH2 georeferenced nodata = 0.
- **QGIS / ArcGIS:** open the `.tif`.

## Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `command failed` / `cmd failed` from a warp | QGIS path wrong — check `QGIS_BIN` in `common.py`; confirm `gdalwarp.exe` exists there. |
| `ModuleNotFoundError: rasterio` (or numpy/scipy) | Wrong interpreter — use `C:\Users\HP\anaconda3\python.exe`, not system Python. |
| CH1 longitude looks collapsed/wrong | You ran only `ch1_pipeline.py`. Run `ch1_fix_warp.py` (GCP+TPS) — that is the correct product. |
| Step 2/3/4/5 errors "file not found" | A prior step's `.img` is missing — steps are sequential; rerun the earlier step. |
| `step3` fails loading `angle.npy` | Run `ch2_pipeline.py angle` first. |
| Out of disk | `outputs/` is ~87 GB and fully regenerable; clear old steps or free space on `D:`. |
| Step 5 stopped mid‑georeference | Run `finish_step5.py`. |

## Rerun / clean notes

- Every step overwrites its own output folder, so reruns are safe and idempotent.
- `outputs/` is regenerable from `scripts/` + raw inputs; it can be deleted to reclaim space.
- To regenerate everything from scratch: run section 1, then section 2 (order matters
  within CH2; CH1 is independent of CH2).
