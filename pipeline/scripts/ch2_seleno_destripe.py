"""
CH2 IIRS  step 4 (selenoreferencing) + step 5 (destriping).

Selenoreferencing: 4 PDS4 corner coords -> GCP warp to Equirectangular Moon 2000
                   (sphere R=1737.4 km), bilinear, 94.56 m/pixel.
Destriping       : per-detector (per cross-track sample column) moment matching,
                   done in SENSOR geometry (correct domain) then georeferenced.

Note on order: the requested order lists destriping after selenoreferencing, but
striping is a sensor-domain (per-detector-column) artefact and MUST be removed
before warping to be effective. We therefore (a) provide step4 = georeferenced
photometric product, and (b) provide step5 = destriped-in-sensor-space product,
georeferenced, as the definitive final cube.
"""
import os, sys
import numpy as np
import rasterio
from rasterio.control import GroundControlPoint
from rasterio.crs import CRS
from rasterio.warp import calculate_default_transform, reproject, Resampling
sys.path.insert(0, os.path.dirname(__file__))
from common import (OUT, IIRS_S, IIRS_L, MOON_GEOG, MOON_EQC, IIRS_CORNERS,
                    load_solar_flux)

C2 = os.path.join(OUT, "ch2")
NB = 247
RES = 94.56
WL = load_solar_flux()[0][7:254]
GCP_CRS = CRS.from_proj4(MOON_GEOG)
DST_CRS = CRS.from_proj4(MOON_EQC)


def _gcps():
    g = []
    for k in ("UL", "UR", "LL", "LR"):
        row, col, lon, lat = IIRS_CORNERS[k]
        g.append(GroundControlPoint(row=row + 0.5, col=col + 0.5, x=lon, y=lat, z=0.0))
    return g


def _warp_cube(src_memmap, out_base, desc):
    """Warp a (NB,L,S) sensor cube to Moon eqc; write GeoTIFF + ENVI (both georeferenced)."""
    gcps = _gcps()
    dst_transform, w, h = calculate_default_transform(
        GCP_CRS, DST_CRS, IIRS_S, IIRS_L, gcps=gcps, resolution=RES)
    print(f"   output grid {w} x {h}  @ {RES} m")

    os.makedirs(os.path.dirname(out_base), exist_ok=True)
    tif = out_base + ".tif"
    img = out_base + ".img"
    prof = dict(driver="GTiff", height=h, width=w, count=NB, dtype="float32",
                crs=DST_CRS, transform=dst_transform, nodata=0.0,
                compress="lzw", interleave="band", BIGTIFF="YES")
    eprof = dict(prof); eprof.update(driver="ENVI", compress=None)
    eprof.pop("interleave", None); eprof.pop("compress", None); eprof.pop("BIGTIFF", None)
    eprof["interleave"] = "bsq"

    with rasterio.open(tif, "w", **prof) as dtif, rasterio.open(img, "w", **eprof) as denvi:
        for b in range(NB):
            dst = np.zeros((h, w), dtype='float32')
            reproject(source=np.ascontiguousarray(src_memmap[b]), destination=dst,
                      src_crs=GCP_CRS, gcps=gcps, dst_transform=dst_transform,
                      dst_crs=DST_CRS, src_nodata=0.0, dst_nodata=0.0,
                      resampling=Resampling.bilinear)
            dtif.write(dst, b + 1); dtif.set_band_description(b + 1, f"{WL[b]:.3f} nm")
            denvi.write(dst, b + 1)
            if b % 50 == 0:
                print(f"     band {b}/{NB}")
    _append_wavelengths(img + ".hdr")
    print("   ->", tif, "\n   ->", img)


def _append_wavelengths(hdr):
    if not os.path.exists(hdr):
        hdr = hdr.replace(".img.hdr", ".hdr")
    with open(hdr, "a") as f:
        f.write("\nwavelength units = Nanometers\nwavelength = {\n")
        v = ["%.4f" % x for x in WL]
        for i in range(0, len(v), 6):
            f.write(" " + ", ".join(v[i:i+6]) + (",\n" if i + 6 < len(v) else "}\n"))


# ---------------------------------------------------------------- STEP 4
def step4_seleno():
    print("[CH2 step4] selenoreferencing photometric cube -> Moon 2000 eqc")
    src = np.memmap(os.path.join(C2, "step3_photometric",
                    "ch2_iirs_step3_reflectance_photom.img"),
                    dtype='<f4', mode='r', shape=(NB, IIRS_L, IIRS_S))
    _warp_cube(src, os.path.join(C2, "step4_seleno", "ch2_iirs_step4_selenoref"),
               "Selenoreferenced photometric reflectance (Moon2000 eqc)")


# ---------------------------------------------------------------- STEP 5
def _destripe_sensor():
    """Per-detector-column moment matching on the sensor-space photometric cube."""
    src = np.memmap(os.path.join(C2, "step3_photometric",
                    "ch2_iirs_step3_reflectance_photom.img"),
                    dtype='<f4', mode='r', shape=(NB, IIRS_L, IIRS_S))
    out = np.empty((NB, IIRS_L, IIRS_S), dtype='float32')
    for b in range(NB):
        img = src[b].astype('float64')             # (L,S)
        valid = img > 0
        cnt = valid.sum(0)
        col_mean = np.where(cnt > 0, np.where(valid, img, 0).sum(0) / np.maximum(cnt, 1), 0)
        col_var = np.where(cnt > 0,
                           np.where(valid, (img - col_mean) ** 2, 0).sum(0) / np.maximum(cnt, 1), 0)
        col_std = np.sqrt(col_var)
        g_mean = col_mean[cnt > 0].mean() if np.any(cnt > 0) else 0.0
        g_std = col_std[cnt > 0].mean() if np.any(cnt > 0) else 1.0
        cs = np.where(col_std > 1e-9, col_std, 1.0)
        corr = (img - col_mean) / cs * g_std + g_mean
        corr[~valid] = 0.0
        out[b] = corr.astype('float32')
    return out


def step5_destripe():
    print("[CH2 step5] destriping (per-detector moment matching, sensor domain)")
    ds = _destripe_sensor()
    # sensor-space ENVI + GeoTIFF
    from common import write_envi, write_geotiff
    base = os.path.join(C2, "step5_destripe", "ch2_iirs_step5_destriped_sensor")
    write_envi(base + ".img", ds, wavelengths=WL,
               description="Destriped reflectance (per-detector moment matching), sensor geometry")
    write_geotiff(base + ".tif", ds, wavelengths=WL)
    print("   ->", base + ".img")
    # georeferenced final
    print("   georeferencing destriped cube ...")
    _warp_cube(ds, os.path.join(C2, "step5_destripe", "ch2_iirs_step5_destriped_selenoref"),
               "Destriped + selenoreferenced reflectance (final product, Moon2000 eqc)")


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("4", "all"):
        step4_seleno()
    if what in ("5", "all"):
        step5_destripe()
