"""
CH2 (Chandrayaan-2 IIRS) processing pipeline.

Input  : IIRS L1 radiance cube (provided .qub, integer, radiance = DN * 0.01)
Steps  : 1 radiometric calibration (radiance)  -> step1_radiance
         .  incidence-angle image (_angle)      -> aux
         2 thermal removal (I/F reflectance)    -> step2_thermal   (247 bands)
         3 photometric correction (Lommel-Seeliger) -> step3_photometric
   (selenoreferencing + destriping are in ch2_seleno_destripe.py)

Replicates the ISRO CH2IIRS reference (Verma, Chauhan & Chauhan 2022, Icarus 383,115075)
for the thermal / reflectance math, exactly.
"""
import os, re, sys
import numpy as np
from datetime import datetime
from scipy.interpolate import interp1d
sys.path.insert(0, os.path.dirname(__file__))
from common import (IIRS_QUB, IIRS_XML, IIRS_SPM, IIRS_S, IIRS_L, IIRS_B, IIRS_LINE_EXP,
                    RAD_SCALE, OUT, C_LIGHT, H_PLANCK, K_BOLTZ, load_solar_flux,
                    write_envi, write_geotiff)

C2 = os.path.join(OUT, "ch2")

# 247-band empirical spectral correction (CH2IIRS Corrector.py, verbatim)
COEFF = np.array([1.007951631,1.003136498,1.012137177,0.983221552,0.982301846,0.980570808,0.974865433,0.97640194,0.992019454,1.007100398,1.004167929,0.986444903,1.005016461,1.013553818,1.009490688,1.003308071,1.011969658,1.0141375,1.006193945,1.015513554,1.020218579,1.016652203,1.023330972,1.020211072,1.015149832,1.00860599,1.000767149,0.994761146,1.067891991,1.044841918,0.972851837,0.9762951,0.979046208,0.97102683,0.982455709,0.987001059,0.988436954,0.979160136,0.9986847,0.997752521,0.998447812,1.016854934,0.997900542,0.983831379,0.969376425,0.979670198,0.973878228,0.996289037,1.014184498,1.01477613,1.001481496,0.975778394,0.997684295,0.984763646,0.983395789,0.98823384,0.988012215,0.952374697,0.937581808,1.013480796,1.014063354,1.022099753,1.03204087,1.035484491,1.041175268,1.036634478,1.034566162,1.03809732,1.046558655,1.021576737,1.028588801,1.007904219,0.972531419,0.956629793,0.955287872,0.979925337,1.002466257,1.014420813,1.008463966,1.00003099,1.008471522,0.99282371,0.99980302,0.987547313,0.984068212,0.992388063,0.988139868,1.012650312,1.052139241,1.012472306,0.965437548,0.939243701,0.967441753,0.963544372,0.990952989,0.973117569,1.013248635,1.020203633,1.029157186,1.020485688,1.020579985,1.006238104,1.014835128,1.006275303,1.010898129,1.00045314,1.008028897,1.012421421,1.027533403,1.031189861,1.029572772,0.998242,0.988472694,0.964883574,0.965888213,0.949632716,0.959468989,0.953967103,0.97185828,0.975099523,0.99858193,0.998573727,1.021937907,1.019142949,1.038073792,1.031876174,1.041953094,1.029414494,1.031461328,1.008333598,1.0158236,0.997626672,0.997772942,0.951297807,0.98848544,0.983055255,0.973287938,1.020004362,1.007683124,1.002357289,1.009943767,0.967795252,0.999227534,0.97305723,0.984516543,0.996874775,1.05340332,1.037099582,1.015145219,1.02252679,0.987330432,0.896931034,0.854160633,1.141899078,1.108682789,1.023030659,0.939513645,0.914993685,0.994857247,1.264567511,0.759581027,1.008047196,1.012150301,0.962843991,1.019548525,0.989332932,1.049050719,0.980318607,1.018721793,0.965286873,1.027176253,0.963709485,1.010598964,0.971110269,1.076940686,1.020255436,0.929940419,0.996612177,0.958811048,1.059205186,1.002448358,0.973559634,1.002937935,1.013683258,1.0474072,1.000925241,0.985500261,0.982660401,0.938587583,1.046900212,0.952743318,1.026649622,1.115011575,0.887761742,1.005459078,0.95483364,0.967275555,1.13852458,0.986343909,0.965651948,1.101161236,0.988116073,0.853138055,0.974375472,1.177442497,0.927739407,0.932486397,0.962300997,0.985181241,0.829206405,1.493934047,1.059590538,0.71322689,0.875132631,1.121072252,0.938541303,1.139668914,0.897165185,1.117762208,0.961511703,0.528956345,1.711570255,1.20749058,1.024243553,0.713165588,1.06875679,0.570841199,0.932178556,2.546499369,2.664596812,1.224978705,0.567434034,0.726464718,0.745085856,1.21801428,1.140071493,0.828417838,1.063662997,0.953993832,1.084346492,0.826805177,1.166156321,0.873581953,0.906603941,1.403382254,0.898905495,1.006518779])

BSTART, BEND = 7, 254           # reference [7:-2] on 256 bands -> indices 7..253 (247 bands)


def _iirs_wavelengths_all():
    wl, _ = load_solar_flux()
    return wl                                   # 256 calibrated band centres (nm)


# ---------------------------------------------------------------- STEP 1
def step1_radiance():
    """Ingest L1 integer cube -> physical radiance (mW/cm^2/sr/um) = DN*0.01, 256 bands."""
    print("[CH2 step1] radiometric calibration (radiance)")
    wl = _iirs_wavelengths_all()
    src = np.memmap(IIRS_QUB, dtype='<u2', mode='r', shape=(IIRS_B, IIRS_L, IIRS_S))
    out_img = os.path.join(C2, "step1_radiance",
                           "ch2_iirs_step1_radiance.img")
    os.makedirs(os.path.dirname(out_img), exist_ok=True)
    # stream write ENVI BSQ band-by-band
    with open(out_img, "wb") as f:
        for b in range(IIRS_B):
            (src[b].astype('float32') * RAD_SCALE).tofile(f)
    from common import _write_hdr
    _write_hdr(out_img + ".hdr", IIRS_S, IIRS_L, IIRS_B, 4, wl, None,
               "Nanometers", None,
               "IIRS L1 radiance (DN x 0.01), mW/cm^2/sr/um; radiometric calibration step", None)
    # GeoTIFF (sensor geometry, no CRS)
    import rasterio
    from rasterio.transform import Affine
    tif = out_img.replace(".img", ".tif")
    prof = dict(driver="GTiff", height=IIRS_L, width=IIRS_S, count=IIRS_B,
                dtype="float32", compress="lzw", interleave="band",
                BIGTIFF="IF_SAFER", transform=Affine.identity())
    with rasterio.open(tif, "w", **prof) as dst:
        for b in range(IIRS_B):
            dst.write(src[b].astype('float32') * RAD_SCALE, b + 1)
            dst.set_band_description(b + 1, f"{wl[b]:.3f} nm")
    print("   ->", out_img, "\n   ->", tif)


# ---------------------------------------------------------------- ANGLE
def gen_angle():
    """Per-line solar incidence-angle image from .spm (CH2IIRS Incidence.py, verbatim math)."""
    print("[CH2 aux] incidence-angle image from .spm")
    dom = open(IIRS_XML, encoding="utf-8", errors="ignore").read()
    s = re.search(r"<start_date_time>(.*?)</start_date_time>", dom).group(1)
    s = re.sub(r'(\.\d{6})\d+(Z)', r'\1\2', s)
    start_time = datetime.strptime(s, '%Y-%m-%dT%H:%M:%S.%fZ').timestamp()
    tm, ang = [], []
    with open(IIRS_SPM) as inp:
        for line in inp:
            sp = line.strip().split()
            if len(sp) < 19:
                continue
            elev = abs(float(sp[-1]))                     # sun elevation
            yr = sp[2][3:]                                # block-len+year merged
            t = datetime.strptime(str([yr] + sp[3:9]),
                                  "['%Y', '%m', '%d', '%H', '%M', '%S', '%f']").timestamp()
            tm.append(t); ang.append(elev)
    f = interp1d(tm, ang, kind='linear', bounds_error=False, fill_value='extrapolate')
    inc = np.zeros((IIRS_L, IIRS_S), dtype='float32')
    for i in range(IIRS_L):
        inc[i, :] = 90.0 - f(start_time + IIRS_LINE_EXP * i)
    out = os.path.join(C2, "aux", "ch2_iirs_angle.img")
    write_envi(out, inc[None, :, :], description="Solar incidence angle (deg) per line")
    np.save(os.path.join(C2, "aux", "angle.npy"), inc)
    print(f"   incidence range {inc.min():.2f}..{inc.max():.2f} deg  -> {out}")
    return inc


# ---------------------------------------------------------------- STEP 2
def step2_thermal(chunk=1500):
    """Thermal emission removal + I/F reflectance (CH2IIRS Corrector.py, verbatim)."""
    print("[CH2 step2] thermal removal -> reflectance (247 bands)")
    wl_all, fx_all = load_solar_flux()
    wl = wl_all[BSTART:BEND]                       # 247 nm
    fx = np.round(fx_all[BSTART:BEND], 4) * 10.0   # 247 solar flux (x10)
    lmb = wl * 1e-9                                # m
    tb_lo = int(np.argmin(np.abs(wl - 4500)))
    tb_hi = int(np.argmin(np.abs(wl - 4874)))
    e = 0.95
    coeff = COEFF.reshape(-1, 1, 1)
    fx3 = fx.reshape(-1, 1, 1)
    lmb3 = lmb.reshape(-1, 1, 1)
    lmbS = lmb[tb_lo:tb_hi].reshape(-1, 1, 1)
    NB = BEND - BSTART

    rad_src = np.memmap(os.path.join(C2, "step1_radiance", "ch2_iirs_step1_radiance.img"),
                        dtype='<f4', mode='r', shape=(IIRS_B, IIRS_L, IIRS_S))
    out = np.zeros((NB, IIRS_L, IIRS_S), dtype='float32')
    temp = np.zeros((IIRS_L, IIRS_S), dtype='float32')

    for y0 in range(0, IIRS_L, chunk):
        y1 = min(y0 + chunk, IIRS_L)
        rad = rad_src[BSTART:BEND, y0:y1, :].astype('float64')      # (247,ny,S) physical radiance
        rad = np.clip(rad, 1e-9, None)
        ra = rad[tb_lo:tb_hi]
        q = np.log(e * 2 * H_PLANCK * C_LIGHT**2 * 1e-6 / (ra * lmbS**5) + 1)
        Tband = H_PLANCK * C_LIGHT / (lmbS * K_BOLTZ * q)
        mT = np.nanmean(Tband, axis=0)                              # (ny,S)
        temp[y0:y1] = mT.astype('float32')
        B = 1e-6 * (2 * H_PLANCK * C_LIGHT**2 / lmb3**5) / \
            (np.exp(H_PLANCK * C_LIGHT / (lmb3 * K_BOLTZ * mT[None])) - 1)   # thermal radiance
        ref = 3.14 * (rad - e * B) / fx3
        ref[ref < 0] = 0
        ref /= coeff
        out[:, y0:y1, :] = _smooth_spectral(ref).astype('float32')
        print(f"   lines {y0}-{y1}  meanT={np.nanmean(mT):.1f}K")

    # write temperature (1 band) + reflectance (247 bands)
    write_envi(os.path.join(C2, "aux", "ch2_iirs_Temp.img"), temp[None],
               description="Brightness temperature (K), Planck inversion 4500-4874nm, e=0.95")
    _write_step("step2_thermal", "ch2_iirs_step2_reflectance_thermcorr", out, wl,
                "Thermally-corrected I/F reflectance (thermal removal), 247 bands 847-4993 nm")
    print("   done step2")


def _smooth_spectral(x):
    """3-point spectral moving average, CH2IIRS convention (endpoints special-cased)."""
    n = x.shape[0]
    z1 = np.zeros((1,) + x.shape[1:]); z2 = np.zeros((2,) + x.shape[1:])
    a1 = np.concatenate([x, z2], 0)
    a2 = np.concatenate([z1, x, z1], 0)
    a3 = np.concatenate([z2, x], 0)
    mid = ((a1 + a2 + a3) / 3.0)[2:-2]
    out = np.empty_like(x)
    out[0] = (x[1] + x[2]) / 2.0
    out[1:n-1] = mid
    out[n-1] = (x[-1] + x[-2]) / 2.0
    return out


# ---------------------------------------------------------------- STEP 3
def step3_photometric(inc):
    """Lommel-Seeliger photometric correction to standard geometry i=30,e=0,g=30.

    D(i,e)=cos i/(cos i+cos e); with nadir emission e~0 -> D(i,0)=cos i/(cos i+1).
    Correct to standard: ref_std = ref * D(30,0)/D(i,0) = ref*(cos i+1)/cos i * 0.46410.
    (The empirical phase-curve term of Verma et al. 2023 needs ISRO's non-public
     polynomial fits; the analytic Lommel-Seeliger disk normalisation is applied here.)
    """
    print("[CH2 step3] photometric correction (Lommel-Seeliger -> i=30,e=0,g=30)")
    wl_all, _ = load_solar_flux()
    wl = wl_all[BSTART:BEND]
    src = np.memmap(os.path.join(C2, "step2_thermal", "ch2_iirs_step2_reflectance_thermcorr.img"),
                    dtype='<f4', mode='r', shape=(BEND - BSTART, IIRS_L, IIRS_S))
    ang = np.clip(inc, None, 80.0)
    cosi = np.cos(np.deg2rad(ang))
    D30 = np.cos(np.deg2rad(30.0)) / (np.cos(np.deg2rad(30.0)) + 1.0)   # 0.46410
    factor = ((cosi + 1.0) / cosi) * D30                               # (L,S)
    out = np.empty((BEND - BSTART, IIRS_L, IIRS_S), dtype='float32')
    for b in range(BEND - BSTART):
        out[b] = (src[b] * factor).astype('float32')
    _write_step("step3_photometric", "ch2_iirs_step3_reflectance_photom", out, wl,
                "Photometrically corrected reflectance (Lommel-Seeliger, i=30 e=0 g=30)")
    print("   done step3")


def _write_step(folder, name, cube, wl, desc):
    img = os.path.join(C2, folder, name + ".img")
    write_envi(img, cube, wavelengths=wl, description=desc)
    write_geotiff(img.replace(".img", ".tif"), cube, wavelengths=wl)
    print("   ->", img)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("1", "all"):
        step1_radiance()
    if what in ("angle", "all"):
        gen_angle()
    if what in ("2", "all"):
        step2_thermal()
    if what in ("3", "all"):
        inc = np.load(os.path.join(C2, "aux", "angle.npy"))
        step3_photometric(inc)
