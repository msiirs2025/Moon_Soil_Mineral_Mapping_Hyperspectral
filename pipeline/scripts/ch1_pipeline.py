"""
CH1 (Chandrayaan-1 M3) selenoreferencing.

Input : M3 L2 reflectance (rfl, 85 bands, BIL) + LOC backplane (lon/lat/radius).
Method: rigorous per-pixel geolocation-array warp (gdalwarp -geoloc) using the LOC
        lon/lat of every pixel -> Equirectangular Moon 2000 (sphere R=1737.4 km).
Note  : this is a near-polar orbital strip (55 N -> ~87 S). The final ~2500 lines
        cross immediately beside the south pole where longitude is undefined/degenerate
        for an equirectangular grid; those lines (lat < -80) are excluded from the eqc
        product (eqc is singular near the pole). They remain in the sensor-geometry data;
        a south-polar-stereographic tail can be produced on request.
"""
import os, sys, subprocess, re
import numpy as np
import rasterio
from rasterio.transform import Affine
sys.path.insert(0, os.path.dirname(__file__))
from common import (OUT, MOON_GEOG, MOON_EQC, QGIS_BIN, GDALWARP, GDALTRANS)

M3_DIR = r"D:\new project moon\ch1\cartOrder (1)\cartorder"
RFL = os.path.join(M3_DIR, "m3g20090731t045352_v01_rfl.img")
RFL_HDR = os.path.join(M3_DIR, "m3g20090731t045352_v01_rfl.hdr")
LOC = os.path.join(M3_DIR, "m3g20090731t045352_v03_loc.img")
L, B, S = 29939, 85, 304          # lines, bands, samples
NODATA = -999.0
LAT_CLIP = -80.0                   # exclude near-pole tail from eqc product
RES = 140.0                       # m/pixel (M3 global mode)

C1 = os.path.join(OUT, "ch1", "seleno")
ENV = dict(os.environ, GDAL_DATA=os.path.join(os.path.dirname(QGIS_BIN), "share", "gdal"),
           PROJ_LIB=os.path.join(os.path.dirname(QGIS_BIN), "share", "proj"),
           PATH=QGIS_BIN + os.pathsep + os.environ.get("PATH", ""))


def parse_wavelengths():
    txt = open(RFL_HDR).read()
    m = re.search(r"wavelength\s*=\s*\{([^}]*)\}", txt)
    return [float(x) for x in m.group(1).replace("\n", " ").split(",") if x.strip()]


def run(cmd):
    print("  $", os.path.basename(cmd[0]), " ".join(cmd[1:])[:120])
    p = subprocess.run(cmd, capture_output=True, text=True, env=ENV)
    if p.returncode != 0:
        print(p.stdout[-1500:]); print(p.stderr[-1500:]); raise RuntimeError("cmd failed")
    return p


def main():
    os.makedirs(C1, exist_ok=True)
    wl = parse_wavelengths()
    print(f"[CH1] {len(wl)} bands, {wl[0]}-{wl[-1]} nm")

    rfl = np.memmap(RFL, dtype='<f4', mode='r', shape=(L, B, S))       # BIL: line,band,sample
    loc = np.memmap(LOC, dtype='<f8', mode='r', shape=(L, 3, S))
    latc = loc[:, 1, S // 2]
    cut = int(np.argmax(latc < LAT_CLIP)) if np.any(latc < LAT_CLIP) else L
    print(f"[CH1] clipping near-pole tail at line {cut}/{L} (lat<{LAT_CLIP}); "
          f"lon range kept {loc[:cut,0,:].min():.2f}..{loc[:cut,0,:].max():.2f}")

    # 1) sensor-geometry rasters (clipped): reflectance (85b), lon, lat
    sensor_tif = os.path.join(C1, "m3_sensor_rfl.tif")
    lon_tif = os.path.join(C1, "m3_lon.tif")
    lat_tif = os.path.join(C1, "m3_lat.tif")
    prof = dict(driver="GTiff", height=cut, width=S, dtype="float32", crs=None,
                transform=Affine.identity(), compress="lzw", BIGTIFF="IF_SAFER")
    if not (os.path.exists(sensor_tif) and os.path.exists(lon_tif) and os.path.exists(lat_tif)):
        print("[CH1] writing sensor reflectance + geolocation rasters ...")
        with rasterio.open(sensor_tif, "w", count=B, nodata=NODATA, **prof) as d:
            for b in range(B):
                d.write(np.ascontiguousarray(rfl[:cut, b, :]), b + 1)
                d.set_band_description(b + 1, f"{wl[b]:.2f} nm")
        with rasterio.open(lon_tif, "w", count=1, dtype="float64",
                           **{k: v for k, v in prof.items() if k != "dtype"}) as d:
            d.write(np.ascontiguousarray(loc[:cut, 0, :].astype('float64')), 1)
        with rasterio.open(lat_tif, "w", count=1, dtype="float64",
                           **{k: v for k, v in prof.items() if k != "dtype"}) as d:
            d.write(np.ascontiguousarray(loc[:cut, 1, :].astype('float64')), 1)
    else:
        print("[CH1] sensor/lon/lat rasters already present -> reuse")

    # 2) VRT with GEOLOCATION metadata pointing at lon/lat rasters
    vrt = os.path.join(C1, "m3_geoloc.vrt")
    run([GDALTRANS, "-of", "VRT", sensor_tif, vrt])
    _inject_geoloc(vrt, lon_tif, lat_tif)

    # 3) warp with per-pixel geolocation -> Moon eqc
    out_tif = os.path.join(C1, "m3g20090731t045352_seleno_rfl.tif")
    print("[CH1] gdalwarp -geoloc -> Equirectangular Moon 2000 ...")
    run([GDALWARP, "-geoloc", "-t_srs", MOON_EQC, "-r", "bilinear",
         "-srcnodata", "-999", "-dstnodata", "-999", "-tr", str(RES), str(RES),
         "-ovr", "NONE", "-co", "COMPRESS=LZW", "-co", "BIGTIFF=YES",
         "-overwrite", vrt, out_tif])

    # 4) ENVI copy (.img/.hdr) + wavelengths + data ignore
    out_img = os.path.join(C1, "m3g20090731t045352_seleno_rfl.img")
    run([GDALTRANS, "-of", "ENVI", "-co", "INTERLEAVE=BSQ", out_tif, out_img])
    _finalize_hdr(out_img, wl)

    with rasterio.open(out_tif) as ds:
        print(f"[CH1] DONE  {ds.width} x {ds.height}  bands {ds.count}  "
              f"crs eqc  bounds {[round(b,0) for b in ds.bounds]}")
    print("   ->", out_tif, "\n   ->", out_img)


def _inject_geoloc(vrt, lon_tif, lat_tif):
    srs = rasterio.crs.CRS.from_proj4(MOON_GEOG).to_wkt().replace('"', "&quot;")
    lon_p = lon_tif.replace("\\", "/"); lat_p = lat_tif.replace("\\", "/")
    block = ('  <Metadata domain="GEOLOCATION">\n'
             f'    <MDI key="SRS">{srs}</MDI>\n'
             f'    <MDI key="X_DATASET">{lon_p}</MDI>\n    <MDI key="X_BAND">1</MDI>\n'
             f'    <MDI key="Y_DATASET">{lat_p}</MDI>\n    <MDI key="Y_BAND">1</MDI>\n'
             '    <MDI key="PIXEL_OFFSET">0</MDI>\n    <MDI key="PIXEL_STEP">1</MDI>\n'
             '    <MDI key="LINE_OFFSET">0</MDI>\n    <MDI key="LINE_STEP">1</MDI>\n'
             '    <MDI key="GEOREFERENCING_CONVENTION">PIXEL_CENTER</MDI>\n'
             '  </Metadata>\n')
    t = open(vrt).read()
    idx = t.index(">", t.index("<VRTDataset")) + 1      # after the opening tag
    t = t[:idx] + "\n" + block + t[idx:]
    open(vrt, "w").write(t)


def _finalize_hdr(out_img, wl):
    hdr = out_img[:-4] + ".hdr" if os.path.exists(out_img[:-4] + ".hdr") else out_img + ".hdr"
    if not os.path.exists(hdr):
        hdr = out_img.replace(".img", ".hdr")
    t = open(hdr).read()
    if "data ignore value" not in t:
        t = t.rstrip() + "\ndata ignore value = -999.0\n"
    t = t.rstrip() + "\nwavelength units = Nanometers\nwavelength = {\n"
    v = ["%.4f" % x for x in wl]
    for i in range(0, len(v), 6):
        t += " " + ", ".join(v[i:i+6]) + (",\n" if i + 6 < len(v) else "}\n")
    open(hdr, "w").write(t)


if __name__ == "__main__":
    main()
