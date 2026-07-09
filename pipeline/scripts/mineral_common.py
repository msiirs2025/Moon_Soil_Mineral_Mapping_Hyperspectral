"""
Shared helpers for the RELAB-referenced mineral analysis of the ENVI ROI subsets
(CH1 = M3, CH2 = IIRS).

Design (matches the rest of the pipeline):
- Read the extensionless ENVI subset binaries directly by memmap using their .hdr.
- All comparison against the RELAB laboratory library is CAPPED AT 2600 nm (the library's
  upper limit): CH1 461-2600 nm, CH2 830-2600 nm. Bands beyond 2600 nm are never matched
  or plotted.
- No osgeo; GeoTIFF via rasterio, plots via matplotlib.
"""
import os, re, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from common import OUT, MOON_R, MOON_EQC

# ---------------------------------------------------------------- paths
SUBSET_DIR = os.path.join(OUT, "Subset")
SPECLIB    = os.path.join(OUT, "spectral_library")
ANALYSIS   = os.path.join(OUT, "mineral_analysis")

# RELAB reference file per mineral (user-supplied relab_map)
RELAB_MAP = {
    "olivine":       "LR-CMP-212.txt",
    "clinopyroxene": "LR-CMP-219.txt",
    "orthopyroxene": "DL-CMP-024-A.txt",
    "plagioclase":   "LR-CMP-217.txt",
    "ilmenite":      "LR-CMP-218.txt",
    "mg_spinel":     "SP-AHT-044.txt",
}
MINERAL_ORDER = ["olivine", "clinopyroxene", "orthopyroxene",
                 "plagioclase", "ilmenite", "mg_spinel"]

# distinct plotting colour per mineral (Fig-3 style boxes / traces)
MINERAL_COLOR = {
    "olivine":       "#e41a1c",   # red
    "clinopyroxene": "#377eb8",   # blue
    "orthopyroxene": "#4daf4a",   # green
    "plagioclase":   "#984ea3",   # purple
    "ilmenite":      "#ff7f00",   # orange
    "mg_spinel":     "#a6761d",   # brown
}

# ROI subsets (extensionless ENVI binary + <path>.hdr)
# fill_lo: values <= fill_lo are treated as invalid (ENVI resize stores M3's -999 as
# ~-998.9999, so an exact test misses it; reflectance is always > 0 for both sensors).
SUBSETS = {
    "ch1": dict(path=os.path.join(SUBSET_DIR, "ch1_subset", "ch1_subset"),
                label="CH1 / M3",   nodata=-999.0, fill_lo=-900.0),
    "ch2": dict(path=os.path.join(SUBSET_DIR, "ch2_subset", "ch2_subset"),
                label="CH2 / IIRS", nodata=0.0,    fill_lo=0.0),
}

FILL_BAND_FRAC = 0.5      # drop a band from comparison if >50% of its pixels are fill

CMP_MAX_NM = 2600.0                       # RELAB upper limit; cap all comparison here

# Fixed, identical axes for every "Spectral Profile" graph (M3 & IIRS, mineral & crater)
# so the two sensors' plots start from a common origin and are equally spaced / comparable.
PROFILE_XLIM    = (500.0, 2600.0)         # nm  (common x origin for M3 & IIRS)
PROFILE_XTICK   = 250.0                   # nm between x ticks
PROFILE_YTICK   = 0.005                   # reflectance between y ticks
PROFILE_FIGSIZE = (9.0, 9.0)              # taller panel -> curves separate, don't overlap


def profile_ylim(curves, pad_frac=0.06, tick=PROFILE_YTICK):
    """Zoom the y-axis to the plotted data with a little padding, snapped to tick
    multiples for clean tick labels.  `curves` = iterable of 1-D reflectance arrays."""
    v = np.concatenate([np.asarray(c, dtype="float64").ravel() for c in curves])
    v = v[np.isfinite(v)]
    lo, hi = float(v.min()), float(v.max())
    pad = (hi - lo) * pad_frac
    lo = np.floor((lo - pad) / tick) * tick
    hi = np.ceil((hi + pad) / tick) * tick
    return lo, hi

# mafic-diagnostic FCC targets (R, G, B) in nm; nearest band is chosen per sensor
FCC_TARGETS = {"ch1": (1978.0, 990.0, 750.0),
               "ch2": (1976.0, 999.0, 830.0)}

# diagnostic absorption centres to annotate on the plots (nm)
ABS_LINES = (1000.0, 2000.0)


# ---------------------------------------------------------------- ENVI header
def parse_hdr(hdr):
    """Return dict(samples, lines, bands, wl[np], map_info[list|None]) from an ENVI .hdr."""
    txt = open(hdr).read()

    def _int(key):
        return int(re.search(rf"{key}\s*=\s*(\d+)", txt).group(1))

    S, L, B = _int("samples"), _int("lines"), _int("bands")
    m = re.search(r"wavelength\s*=\s*\{([^}]*)\}", txt)
    wl = np.array([float(x) for x in m.group(1).replace("\n", " ").split(",")
                   if x.strip()], dtype="float64")
    mi = re.search(r"map info\s*=\s*\{([^}]*)\}", txt)
    map_info = [x.strip() for x in mi.group(1).split(",")] if mi else None
    return dict(samples=S, lines=L, bands=B, wl=wl, map_info=map_info)


def read_subset(name):
    """Memmap a subset -> (cube[bands,lines,samples], hdr_dict, nodata)."""
    info = SUBSETS[name]
    h = parse_hdr(info["path"] + ".hdr")
    cube = np.memmap(info["path"], dtype="<f4", mode="r",
                     shape=(h["bands"], h["lines"], h["samples"]))
    return cube, h, info["nodata"]


def overlap_idx(wl):
    """Band indices with wavelength <= 2600 nm (the RELAB comparison window)."""
    return np.where(wl <= CMP_MAX_NM)[0]


def comparison_idx(cube, h, fill_lo):
    """Bands used for matching: wavelength <= 2600 nm AND not mostly-fill.
    (e.g. M3 bands 1-2 @461/500 nm are 100% fill in the subset and are dropped.)"""
    oi = overlap_idx(h["wl"])
    keep = []
    for b in oi:
        frac = np.mean(np.asarray(cube[b], dtype="float32") <= fill_lo)
        if frac < FILL_BAND_FRAC:
            keep.append(int(b))
    return np.array(keep)


def valid_mask(cube, fill_lo, idx):
    """Pixels finite and > fill_lo across all bands in idx -> (lines,samples) bool."""
    X = np.asarray(cube[idx], dtype="float32")
    return np.all(np.isfinite(X) & (X > fill_lo), axis=0)


def _hull_continuum(wl, y):
    """Upper convex-hull continuum of a single spectrum (the standard planetary
    continuum: the taut band stretched over the reflectance peaks)."""
    n = len(y)
    idx = [0]
    for i in range(1, n):
        while len(idx) >= 2:
            k, j = idx[-2], idx[-1]
            cross = (wl[j] - wl[k]) * (y[i] - y[k]) - (y[j] - y[k]) * (wl[i] - wl[k])
            if cross >= 0:            # middle point sits below the hull -> drop it
                idx.pop()
            else:
                break
        idx.append(i)
    return np.interp(wl, wl[idx], y[idx])


def continuum_remove_1d(wl, y):
    """Continuum-removed spectrum (hull quotient), 1.0 at the continuum, dips at bands."""
    c = _hull_continuum(wl, np.asarray(y, dtype="float64"))
    return np.asarray(y, dtype="float64") / np.where(c < 1e-9, 1e-9, c)


def continuum_remove_cube(X, wl, mask):
    """Hull continuum removal for every valid pixel of a (bands,lines,samples) cube.
    Invalid pixels are left at 1.0 (no feature)."""
    n, L, S = X.shape
    Xf = np.asarray(X, dtype="float64").reshape(n, -1)
    out = np.ones((n, L * S), dtype="float32")
    for p in np.where(mask.reshape(-1))[0]:
        out[:, p] = continuum_remove_1d(wl, Xf[:, p])
    return out.reshape(n, L, S)


def nearest_band(wl, target_nm):
    """0-based band index whose wavelength is closest to target_nm."""
    return int(np.argmin(np.abs(wl - target_nm)))


# ---------------------------------------------------------------- RELAB library
def load_relab(mineral):
    """Return (wavelength_nm, reflectance) for a mineral from its RELAB .txt."""
    path = os.path.join(SPECLIB, RELAB_MAP[mineral])
    d = np.loadtxt(path, comments="#")
    return d[:, 0].copy(), d[:, 1].copy()


def resample_ref(ref_wl, ref_rf, sensor_wl):
    """Linear-interpolate a RELAB spectrum onto sensor_wl (assumed within ref range)."""
    return np.interp(sensor_wl, ref_wl, ref_rf)


# ---------------------------------------------------------------- SAM
def sff_fit(Xcr, ref_cr, mask):
    """Spectral Feature Fitting score: how well each pixel's continuum-removed absorption
    features match the reference's, in position AND relative depth.

    Fits ref_cr to each pixel by least squares (offset + scale) and scores by the fit
    correlation r (= sign(scale)*sqrt(R^2)); r->1 means the same bands at the same places
    with proportional depths. Unlike SAM (whole-vector angle), this responds only to the
    absorption features, which is the diagnostic information for mineral identification.

    Returns (misfit, r): misfit = 1 - r (lower is better; invalid pixels +inf), and r map.
    """
    n, L, S = Xcr.shape
    A = np.asarray(Xcr, dtype="float64").reshape(n, -1)
    r = np.asarray(ref_cr, dtype="float64")
    Am = A - A.mean(0, keepdims=True)
    rm = r - r.mean()
    num = (Am * rm[:, None]).sum(0)
    den = np.sqrt((Am ** 2).sum(0)) * np.sqrt((rm ** 2).sum()) + 1e-12
    corr = (num / den).reshape(L, S)
    misfit = 1.0 - corr
    misfit[~mask] = np.inf
    return misfit, corr


def pick_points(ang, n=4, min_sep=3):
    """Pick the n lowest-angle pixels, each >= min_sep px from the already-picked ones.
    Returns list of (line, sample)."""
    L, S = ang.shape
    order = np.argsort(ang, axis=None)
    pts = []
    for flat in order:
        a = ang.flat[flat]
        if not np.isfinite(a):
            break
        l, s = divmod(int(flat), S)
        if all((l - pl) ** 2 + (s - ps) ** 2 >= min_sep ** 2 for pl, ps in pts):
            pts.append((l, s))
            if len(pts) == n:
                break
    return pts


# ---------------------------------------------------------------- geo helpers
def geo_transform(h):
    """rasterio Affine + CRS from an ENVI 'map info' (eqc Moon 2000)."""
    from rasterio.transform import Affine
    from rasterio.crs import CRS
    mi = h["map_info"]
    rpx, rpy = float(mi[1]), float(mi[2])           # 1-based reference pixel
    mx, my = float(mi[3]), float(mi[4])             # its map coord (upper-left corner)
    px, py = float(mi[5]), float(mi[6])             # pixel size
    ulx = mx - (rpx - 1.0) * px
    uly = my + (rpy - 1.0) * py
    return Affine(px, 0, ulx, 0, -py, uly), CRS.from_proj4(MOON_EQC)


def pixel_lonlat(h, line, sample):
    """Pixel (line,sample) centre -> (easting, northing, lon_deg, lat_deg) on Moon eqc."""
    tr, _ = geo_transform(h)
    e, n = tr * (sample + 0.5, line + 0.5)          # eqc metres
    lon = np.degrees(e / MOON_R)
    lat = np.degrees(n / MOON_R)
    return e, n, lon, lat


def lonlat_to_pixel(h, lon, lat):
    """Inverse of pixel_lonlat: (lon_deg,lat_deg) -> (line, sample) index on this grid.
    Returns the pixel whose footprint contains the point (floor of the fractional index)."""
    tr, _ = geo_transform(h)
    e = np.radians(lon) * MOON_R
    n = np.radians(lat) * MOON_R
    colf, rowf = (~tr) * (e, n)                      # continuous pixel coords
    return int(np.floor(rowf)), int(np.floor(colf))


def subset_extent(h):
    """Map extent of a subset in eqc metres -> (e_min, e_max, n_min, n_max)."""
    tr, _ = geo_transform(h)
    ulx, uly, px, py = tr.c, tr.f, tr.a, -tr.e
    return (ulx, ulx + px * h["samples"], uly - py * h["lines"], uly)


def dual_valid_overlap(master_h, master_valid, other_h, other_valid):
    """Boolean (Lm,Sm) mask of master pixels that (a) fall inside the OTHER sensor's
    footprint and (b) map to a valid OTHER pixel AND are themselves valid.  Guarantees a
    selected master point has a real, sampleable counterpart in the other sensor."""
    trm, _ = geo_transform(master_h)
    ulx, uly, px, py = trm.c, trm.f, trm.a, -trm.e
    Lm, Sm = master_h["lines"], master_h["samples"]
    cols = np.arange(Sm); rows = np.arange(Lm)
    E = ulx + px * (cols + 0.5)                      # (Sm,)
    N = uly - py * (rows + 0.5)                       # (Lm,)
    EE, NN = np.meshgrid(E, N)                        # (Lm,Sm) eqc metres

    oe0, oe1, on0, on1 = subset_extent(other_h)
    inside = (EE >= oe0) & (EE <= oe1) & (NN >= on0) & (NN <= on1)

    tro, _ = geo_transform(other_h)
    oulx, ouly, opx, opy = tro.c, tro.f, tro.a, -tro.e
    Lo, So = other_h["lines"], other_h["samples"]
    ocol = np.floor((EE - oulx) / opx).astype(int)
    orow = np.floor((ouly - NN) / opy).astype(int)
    ocol = np.clip(ocol, 0, So - 1)
    orow = np.clip(orow, 0, Lo - 1)
    other_ok = other_valid[orow, ocol]
    return master_valid & inside & other_ok


# ---------------------------------------------------------------- diagnostic band depths
def band_depth(Xcr, swl, lo, hi):
    """Deepest absorption in the [lo,hi] nm window from a continuum-removed cube
    (Xcr ~ 1 at the continuum, dips in bands).  Returns depth map (lines,samples),
    depth = 1 - min(hull-quotient) over the window; 0 where the window has no bands."""
    idx = np.where((swl >= lo) & (swl <= hi))[0]
    if idx.size == 0:
        return np.zeros(Xcr.shape[1:], dtype="float32")
    return 1.0 - np.min(np.asarray(Xcr[idx], dtype="float32"), axis=0)


def _norm01(a, mask):
    """Robust 2-98% normalisation of a[mask] to ~[0,1] (for cross-parameter combination)."""
    v = a[mask]
    if v.size == 0:
        return np.zeros_like(a)
    lo, hi = np.percentile(v, [2, 98])
    return np.clip((a - lo) / (hi - lo + 1e-12), 0.0, 1.0)


def diagnostic_scores(X, Xcr, swl, mask):
    """Physically-diagnostic, library-INDEPENDENT selection score per mineral, computed
    directly from the observed cube.  Each score is maximised at pixels showing that
    mineral's characteristic absorptions (Adams 1974; Cloutis 1986; Pieters lunar work):

      olivine        strong broad 1 um, NO 2 um            ->  d1 - d2
      clinopyroxene  1 um + LONG-wavelength 2 um (~2.3 um)  ->  d1 + d2long - d2short
      orthopyroxene  1 um + SHORT-wavelength 2 um (~1.9 um) ->  d1 + d2short - d2long
      plagioclase    weak 1.25 um, weak mafic bands, bright ->  d125 - 0.5(d1+d2)
      ilmenite       opaque: dark + spectrally featureless  ->  -(d1+d2+d125) - albedo
      mg_spinel      strong 2 um, weak/absent 1 um          ->  d2 - d1
    """
    d1   = band_depth(Xcr, swl,  900.0, 1150.0)      # 1 um mafic band
    d125 = band_depth(Xcr, swl, 1200.0, 1350.0)      # plagioclase 1.25 um
    d2s  = band_depth(Xcr, swl, 1800.0, 2050.0)      # opx 2 um (short)
    d2l  = band_depth(Xcr, swl, 2150.0, 2450.0)      # cpx 2 um (long)
    d2   = np.maximum(d2s, d2l)                       # any 2 um band
    ab   = _norm01(X[nearest_band(swl, 1600.0)], mask)   # continuum albedo proxy
    return {
        "olivine":       d1 - d2,
        "clinopyroxene": d1 + d2l - d2s,
        "orthopyroxene": d1 + d2s - d2l,
        "plagioclase":   d125 - 0.5 * (d1 + d2),
        "ilmenite":      -(d1 + d2 + d125) - ab,
        "mg_spinel":     d2 - d1,
    }


def pick_top(score, mask, n=4, min_sep=3):
    """Pick the n highest-score pixels within mask, each >= min_sep px apart.
    Returns list of (line, sample)."""
    s = np.where(mask, np.asarray(score, dtype="float64"), -np.inf)
    order = np.argsort(s, axis=None)[::-1]
    S = s.shape[1]
    pts = []
    for flat in order:
        if not np.isfinite(s.flat[flat]):
            break
        l, d = divmod(int(flat), S)
        if all((l - pl) ** 2 + (d - ps) ** 2 >= min_sep ** 2 for pl, ps in pts):
            pts.append((l, d))
            if len(pts) == n:
                break
    return pts


# ---------------------------------------------------------------- ENVI spectral library
def write_speclib(path_sli, wl, spectra, names):
    """Write an ENVI spectral library: spectra is (n_spectra, n_bands) float32."""
    spectra = np.ascontiguousarray(np.asarray(spectra, dtype="<f4"))
    n, nb = spectra.shape
    os.makedirs(os.path.dirname(path_sli), exist_ok=True)
    spectra.tofile(path_sli)
    with open(path_sli + ".hdr", "w") as f:
        f.write("ENVI\n")
        f.write("description = { Mineral analysis extracted + RELAB reference spectra "
                "(<=2600 nm) }\n")
        f.write(f"samples = {nb}\nlines = {n}\nbands = 1\n")
        f.write("header offset = 0\nfile type = ENVI Spectral Library\n")
        f.write("data type = 4\ninterleave = bsq\nsensor type = Unknown\nbyte order = 0\n")
        f.write("wavelength units = Nanometers\n")
        f.write("spectra names = {\n " + ",\n ".join(names) + "}\n")
        vals = ["%.4f" % w for w in wl]
        f.write("wavelength = {\n")
        for i in range(0, len(vals), 6):
            f.write(" " + ", ".join(vals[i:i + 6]))
            f.write(",\n" if i + 6 < len(vals) else "}\n")
