"""CH1 M3 selenoreference - robust GCP + thin-plate-spline warp (replaces the
geoloc-array attempt, which collapsed longitude). Reuses the sensor raster already
written by ch1_pipeline.py."""
import os, sys, glob, subprocess
import numpy as np
import rasterio
from rasterio.control import GroundControlPoint
from rasterio.crs import CRS
sys.path.insert(0, os.path.dirname(__file__))
from common import OUT, MOON_GEOG, MOON_EQC, QGIS_BIN, GDALWARP, GDALTRANS
import ch1_pipeline as P

C1 = os.path.join(OUT, "ch1", "seleno")
sensor_tif = os.path.join(C1, "m3_sensor_rfl.tif")
out_tif = os.path.join(C1, "m3g20090731t045352_seleno_rfl.tif")
out_img = os.path.join(C1, "m3g20090731t045352_seleno_rfl.img")
ENV = P.ENV


def main():
    wl = P.parse_wavelengths()
    with rasterio.open(sensor_tif) as ds:
        cut, S = ds.height, ds.width
    print(f"[CH1-fix] sensor {S} x {cut}")
    loc = np.memmap(P.LOC, dtype='<f8', mode='r', shape=(P.L, 3, P.S))

    # dense GCP grid from LOC (every ~180 lines x 5 samples)
    rows = list(range(0, cut, 180)) + [cut - 1]
    cols = [8, 90, 152, 214, 295]
    gcps = []
    for r in rows:
        for c in cols:
            lon = float(loc[r, 0, c]); lat = float(loc[r, 1, c])
            if np.isfinite(lon) and np.isfinite(lat) and -90 <= lat <= 90:
                gcps.append(GroundControlPoint(row=r + 0.5, col=c + 0.5, x=lon, y=lat, z=0.0))
    print(f"[CH1-fix] {len(gcps)} GCPs (lon {min(g.x for g in gcps):.2f}..{max(g.x for g in gcps):.2f}, "
          f"lat {min(g.y for g in gcps):.2f}..{max(g.y for g in gcps):.2f})")

    # write GCPs into the sensor GeoTIFF (metadata only)
    with rasterio.open(sensor_tif, "r+") as ds:
        ds.gcps = (gcps, CRS.from_proj4(MOON_GEOG))

    for f in glob.glob(out_tif[:-4] + ".*") + glob.glob(out_img[:-4] + ".*"):
        try: os.remove(f)
        except OSError: pass

    print("[CH1-fix] gdalwarp -tps -> Equirectangular Moon 2000 ...")
    _run([GDALWARP, "-tps", "-t_srs", MOON_EQC, "-r", "bilinear",
          "-srcnodata", "-999", "-dstnodata", "-999", "-tr", str(P.RES), str(P.RES),
          "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES", "-overwrite", sensor_tif, out_tif])

    _run([GDALTRANS, "-of", "ENVI", "-co", "INTERLEAVE=BSQ", out_tif, out_img])
    P._finalize_hdr(out_img, wl)

    with rasterio.open(out_tif) as ds:
        print(f"[CH1-fix] DONE {ds.width} x {ds.height} bands {ds.count} "
              f"bounds {[round(b,0) for b in ds.bounds]}")


def _run(cmd):
    print("  $", os.path.basename(cmd[0]), " ".join(cmd[1:])[:110])
    p = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
    if p.returncode != 0:
        print(p.stdout[-1500:]); print(p.stderr[-1500:]); raise RuntimeError("cmd failed")


if __name__ == "__main__":
    main()
