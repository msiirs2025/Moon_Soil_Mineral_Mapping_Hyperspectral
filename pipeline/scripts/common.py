"""
Common I/O + constants for the Chandrayaan CH1(M3) / CH2(IIRS) processing pipeline.

Design goals:
- Every step writes an ENVI cube (.img raw BSQ float32 + .hdr, fully ENVI-visible with
  wavelengths) AND a multiband GeoTIFF.
- No dependency on osgeo python bindings (not installed): ENVI is written by hand,
  GeoTIFF/warp via rasterio (bundled GDAL 3.10) and the QGIS gdalwarp/gdal_translate exes.
"""
import os, struct, subprocess
import numpy as np

# ---------------------------------------------------------------- paths / tools
BASE      = r"D:\new project moon"
OUT       = os.path.join(BASE, "pipeline", "outputs")
REF_DIR   = os.path.join(BASE, "pipeline", "reference", "CH2IIRS")
SOLAR_TXT = os.path.join(REF_DIR, "Solar flux.txt")

QGIS_BIN  = r"C:\Program Files\QGIS 3.40.12\bin"
GDALWARP  = os.path.join(QGIS_BIN, "gdalwarp.exe")
GDALTRANS = os.path.join(QGIS_BIN, "gdal_translate.exe")

# ---------------------------------------------------------------- CH2 IIRS meta
IIRS_QUB = os.path.join(BASE, r"ch2\ch2\ch2\ch2_iir_nri_20250804T0253259678_d_img_d18",
                        r"data\raw\20250804\ch2_iir_nri_20250804T0253259678_d_img_d18.qub")
IIRS_XML = IIRS_QUB[:-4] + ".xml"
IIRS_SPM = os.path.join(BASE, r"ch2\ch2\ch2\ch2_iir_nri_20250804T0253259678_d_img_d18",
                        r"miscellaneous\calibrated\20250804\ch2_iir_nri_20250804T0253259678_d_img_d18.spm")
IIRS_S, IIRS_L, IIRS_B = 250, 13569, 256      # samples, lines, bands
IIRS_LINE_EXP = 0.05306                        # s per line (53.06 ms)
RAD_SCALE = 0.01                               # stored integer -> radiance (CH2IIRS convention)

# Lunar body (IAU Moon 2000, sphere)
MOON_R = 1737400.0
MOON_GEOG = f"+proj=longlat +a={MOON_R} +b={MOON_R} +no_defs"
MOON_EQC  = f"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a={MOON_R} +b={MOON_R} +units=m +no_defs"

# Physical constants (exactly as CH2IIRS reference, for fidelity)
C_LIGHT, H_PLANCK, K_BOLTZ = 3e8, 6.626e-34, 1.38e-23

# IIRS scene corner geographic coords from PDS4 label (lon, lat)
IIRS_CORNERS = {   # (line, sample) -> (lon_deg, lat_deg)
    "UL": (0,          0,          27.280194, 30.874388),
    "UR": (0,          IIRS_S - 1, 26.399155, 30.900110),
    "LL": (IIRS_L - 1, 0,          25.928290, -5.080639),
    "LR": (IIRS_L - 1, IIRS_S - 1, 25.145961, -5.057851),
}

_DTYPE_ENVI = {np.dtype('float32'): 4, np.dtype('float64'): 5,
               np.dtype('uint16'): 12, np.dtype('int16'): 2, np.dtype('uint8'): 1}


def load_solar_flux():
    """Return (wavelength_nm[256], solar_flux[256]) from the CH2IIRS Solar flux.txt."""
    d = np.loadtxt(SOLAR_TXT)
    return d[:, 0].copy(), d[:, 1].copy()


def write_envi(path_img, cube, wavelengths=None, band_names=None,
               wl_units="Nanometers", data_ignore=None, description="",
               map_info=None):
    """Write float32 BSQ ENVI cube (bands, lines, samples) + .hdr.
    map_info: optional list already formatted for the 'map info' ENVI keyword."""
    cube = np.ascontiguousarray(cube.astype('float32'))
    B, L, S = cube.shape
    os.makedirs(os.path.dirname(path_img), exist_ok=True)
    cube.tofile(path_img)                       # BSQ = band, line, sample order
    _write_hdr(path_img + ".hdr", S, L, B, 4, wavelengths, band_names,
               wl_units, data_ignore, description, map_info)


def _write_hdr(hdr, S, L, B, dtype, wavelengths, band_names,
               wl_units, data_ignore, description, map_info):
    with open(hdr, "w") as f:
        f.write("ENVI\n")
        f.write("description = {\n  %s }\n" % description)
        f.write(f"samples = {S}\nlines   = {L}\nbands   = {B}\n")
        f.write("header offset = 0\nfile type = ENVI Standard\n")
        f.write(f"data type = {dtype}\ninterleave = bsq\n")
        f.write("sensor type = Unknown\nbyte order = 0\n")
        if data_ignore is not None:
            f.write(f"data ignore value = {data_ignore}\n")
        if map_info:
            f.write("map info = {%s}\n" % ", ".join(str(x) for x in map_info))
        if band_names:
            f.write("band names = {\n " + ",\n ".join(band_names) + "}\n")
        if wavelengths is not None:
            f.write(f"wavelength units = {wl_units}\n")
            f.write("wavelength = {\n")
            vals = ["%.4f" % w for w in wavelengths]
            for i in range(0, len(vals), 6):
                f.write(" " + ", ".join(vals[i:i + 6]))
                f.write(",\n" if i + 6 < len(vals) else "}\n")


def write_geotiff(path_tif, cube, wavelengths=None, crs=None, transform=None,
                  gcps=None, gcps_crs=None, nodata=None):
    """Write a multiband GeoTIFF (bands, lines, samples) via rasterio."""
    import rasterio
    from rasterio.transform import Affine
    cube = np.asarray(cube, dtype='float32')
    B, L, S = cube.shape
    os.makedirs(os.path.dirname(path_tif), exist_ok=True)
    prof = dict(driver="GTiff", height=L, width=S, count=B, dtype="float32",
                compress="lzw", interleave="band", BIGTIFF="IF_SAFER")
    if transform is not None:
        prof["transform"] = transform
        if crs is not None:
            prof["crs"] = crs
    else:
        prof["transform"] = Affine.identity()
    if nodata is not None:
        prof["nodata"] = nodata
    with rasterio.open(path_tif, "w", **prof) as dst:
        for b in range(B):
            dst.write(cube[b], b + 1)
        if gcps is not None:
            dst.gcps = (gcps, gcps_crs)
        if wavelengths is not None:
            for b in range(B):
                dst.set_band_description(b + 1, f"{wavelengths[b]:.3f} nm")
    return path_tif


def run(cmd):
    print("  $", " ".join(os.path.basename(c) if i == 0 else c
                           for i, c in enumerate(cmd)))
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        print(p.stdout[-2000:]); print(p.stderr[-2000:])
        raise RuntimeError(f"command failed ({p.returncode})")
    return p
