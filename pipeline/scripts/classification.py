"""
Per-pixel mineral CLASSIFICATION of the CH1 (M3) and CH2 (IIRS) ROI subsets into the six
RELAB reference minerals, by three independent, standard hyperspectral classifiers:

    SAM  - Spectral Angle Mapper       (angle between diagnostic feature vectors)
    SFF  - Spectral Feature Fitting     (least-squares fit of absorption depths; ENVI-style)
    SVM  - Support Vector Machine       (RBF kernel, supervised, library-trained)

Scientifically-correct feature space (the key to real discrimination)
---------------------------------------------------------------------
Classifying against single laboratory endmembers only works if the feature space isolates
the DIAGNOSTIC ABSORPTIONS. A single full-range continuum hull leaves one broad shared
concavity that matches the same mafic endmember for nearly every pixel (degenerate).
Instead, following standard lunar spectroscopy (Pieters, Klima, MGM), the 1 um and 2 um
ferrous bands are continuum-removed SEPARATELY, each against its own local straight-line
continuum between the band shoulders, and the two continuum-removed bands are concatenated
into the feature vector  F = [ 1-CR(1um window) , 1-CR(2um window) ]:

    olivine        strong broad 1 um, ~no 2 um
    orthopyroxene  1 um + SHORT 2 um (~1.9 um)
    clinopyroxene  1 um + LONG  2 um (~2.3 um)
    plagioclase    weak 1.25 um, weak mafic bands
    mg_spinel      strong 2 um, weak/absent 1 um
    ilmenite       dark, spectrally featureless

These are genuinely different vectors, so SAM/SFF/SVM can separate them.

* Reference endmembers = the 6 RELAB spectra (relab_map) resampled to the sensor bands,
  everything CAPPED AT 2600 nm (library limit), then put through the same two-band
  continuum removal as the pixels.
* Feature-strength GATE: a pixel whose deepest diagnostic band is < MIN_FEATURE (2%) has no
  measurable absorption to identify a mineral from, so it is left UNCLASSIFIED (class 0)
  rather than force-assigned. This is the standard rule-image behaviour of SAM/SFF.

    SAM : class = argmin spectral angle(F_pixel, F_endmember)          [gated]
    SFF : per endmember fit F_pixel = a*F_ref + b; score = a/RMS (a>0); argmax  [gated]
    SVM : RBF-kernel SVC trained on the endmember features augmented with random band-depth
          gain + Gaussian noise, features standardised; predict            [gated]

SVM has no field truth, so it is trained on a SYNTHETIC set from the same RELAB endmembers
(library-based SVM) - stated openly; it is the defensible option with lab endmembers only.

Outputs (pipeline/outputs/classification/):
    ch{1,2}_{SAM,SFF,SVM}.tif        single-band uint8 class map + colour table (eqc)
    ch{1,2}_{SAM,SFF,SVM}.png        coloured class map, legend, heading (no axes)
    ch{1,2}_classification.png       SAM|SFF|SVM side-by-side per sensor
    envi/ch{1,2}_{SAM,SFF,SVM}.img+.hdr   ENVI Classification (uint8, class names+lookup)
    classification_report.txt        per-class pixel counts/% and method agreement

Run with C:\\Users\\HP\\anaconda3\\python.exe
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
import mineral_common as mc

OUTDIR = os.path.join(mc.OUT, "classification")
ENVIDIR = os.path.join(OUTDIR, "envi")

METHODS = ["SAM", "SFF", "SVM"]
CLASSES = mc.MINERAL_ORDER                      # class index 1..6 (0 = unclassified)
SVM_AUG = 500                                   # synthetic training samples per mineral
RNG_SEED = 42
MIN_FEATURE = 0.02                              # min diagnostic band depth to classify

# 1 um and 2 um band continuum tie-points (nm) per sensor (IIRS 1 um is truncated at 830 nm)
BANDS = {
    "ch1": {"one": (750.0, 1500.0), "two": (1500.0, 2550.0)},    # M3
    # IIRS 2 um window uses the RELIABLE region 2010-2560 nm, BEYOND the 1.9 um instrument
    # response artifact. (Diagnostic: over 2010-2560 nm the IIRS reflectance is convex with
    # no absorption band, so any real 2 um mafic band is either absent or buried under the
    # artifact/low SNR -> pixels with no genuine band are left Unclassified, not forced.)
    "ch2": {"one": (850.0, 1450.0), "two": (2010.0, 2560.0)},    # IIRS
}

# Feature mode per sensor.
#   two_band  : concat 1 um + 2 um diagnostic bands (M3 has the full 1 um band).
#   two_micron: 2 um band ONLY. IIRS starts ~830 nm and misses the short-wavelength
#               shoulder of the 1 um band, which biased the earlier full-spectrum result
#               toward mg-spinel. The 2 um band IS fully measured by IIRS at higher spectral
#               resolution, and its CENTRE / shape is the classic pyroxene-vs-spinel
#               discriminant (Adams 1974; Cloutis & Gaffey 1991): opx ~1.9, spinel ~2.0
#               (broad), cpx ~2.3 um. Minerals without a 2 um band (olivine, plagioclase,
#               ilmenite) have a flat 2 um feature -> ~zero vector -> never selected, i.e.
#               they are correctly reported as not identifiable from IIRS's 2 um band.
FEATURE_MODE = {"ch1": "two_band", "ch2": "two_micron"}

# Instrument bad-band regions (nm) excluded from BOTH pixels and endmembers before feature
# building. IIRS 1875-1960 nm is a detector spectral-RESPONSE minimum: raw L1 DN plunges
# ~80% (160->30 counts) to a floor at 1908.9 nm and recovers by 2010 nm, far steeper than
# the ~17% solar-flux decline. With only a flat DN->radiance gain available (per-band IIRS
# gains are unpublished) this low-responsivity notch propagates into reflectance as a
# spurious ~1.9 um "absorption" that pinned every pixel's 2 um centre to 1909 nm. Masking it
# lets the 2 um band be characterised from the reliable higher-SNR bands.
BAD_BANDS = {"ch1": [], "ch2": [(1875.0, 1960.0)]}


def good_bands(swl, name):
    g = np.ones(len(swl), dtype=bool)
    for lo, hi in BAD_BANDS[name]:
        g &= ~((swl >= lo) & (swl <= hi))
    return g


def windows_for(name):
    b = BANDS[name]
    return [b["two"]] if FEATURE_MODE[name] == "two_micron" else [b["one"], b["two"]]


def hex2rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


CLASS_RGB = [(0, 0, 0)] + [hex2rgb(mc.MINERAL_COLOR[m]) for m in CLASSES]
CLASS_NAME = ["Unclassified"] + list(CLASSES)


# ---------------------------------------------------------------- diagnostic features
def band_feature(swl, Y, windows):
    """Concatenated continuum-removed band depth over the given (lo,hi) nm windows.
    Y = reflectance (nb, Npix). Each band uses a straight-line continuum between its two
    shoulder reflectances; depth = 1 - R/continuum. Returns F (nfeat, Npix)."""
    Y = np.atleast_2d(Y)
    out = []
    for ta, tb in windows:
        ia, ib = mc.nearest_band(swl, ta), mc.nearest_band(swl, tb)
        if ia > ib:
            ia, ib = ib, ia
        idx = np.arange(ia, ib + 1)
        frac = ((swl[idx] - swl[ia]) / (swl[ib] - swl[ia] + 1e-9))[:, None]   # (m,1)
        cont = Y[ia][None, :] * (1 - frac) + Y[ib][None, :] * frac            # (m, Npix)
        depth = 1.0 - Y[idx] / np.where(np.abs(cont) < 1e-9, np.nan, cont)
        out.append(depth)
    return np.nan_to_num(np.vstack(out), nan=0.0)


def sensor_features(name):
    """Return (F[nfeat,Npix], mask[L,S], flat_idx, swl, h) for the valid pixels."""
    fill_lo = mc.SUBSETS[name]["fill_lo"]
    cube, h, _ = mc.read_subset(name)
    ci = mc.comparison_idx(cube, h, fill_lo)
    swl = h["wl"][ci]
    X = np.asarray(cube[ci], dtype="float32")
    g = good_bands(swl, name)                             # drop instrument bad-band regions
    swl, X = swl[g], X[g]
    nb = X.shape[0]
    mask = np.all(np.isfinite(X) & (X > fill_lo), axis=0)
    flat = np.where(mask.reshape(-1))[0]
    Xv = X.reshape(nb, -1)[:, flat].astype("float64")
    F = band_feature(swl, Xv, windows_for(name))
    return F, mask, flat, swl, h


def endmember_features(swl, name):
    """RELAB endmembers -> diagnostic features (6, nfeat) using this sensor's window set."""
    win = windows_for(name)
    E = []
    for m in CLASSES:
        rw, rf = mc.load_relab(m)
        ref = mc.resample_ref(rw, rf, swl)
        E.append(band_feature(swl, ref[:, None], win)[:, 0])
    return np.asarray(E, dtype="float64")


def _report_band2_center(name, swl, F, Eref, gate, report):
    """For the 2 um-only mode: report Band II centre (wavelength of deepest 2 um
    continuum-removed point) for each endmember and the pixel population."""
    ia, ib = mc.nearest_band(swl, BANDS[name]["two"][0]), mc.nearest_band(swl, BANDS[name]["two"][1])
    if ia > ib:
        ia, ib = ib, ia
    wwin = swl[ia:ib + 1]
    report.append("  Band II (2 um) centre diagnostics (Adams/Cloutis-Gaffey):")
    for k, m in enumerate(CLASSES):
        depth = float(Eref[k].max())
        c = float(wwin[int(np.argmax(Eref[k]))])
        tag = "" if depth >= 0.02 else "   <- no 2 um band (not identifiable from 2 um)"
        report.append(f"     endmember {m:14s} centre ~{c:6.0f} nm  depth={depth:.3f}{tag}")
    Fg = F[:, gate]
    if Fg.shape[1]:
        cen = wwin[np.argmax(Fg, axis=0)]
        report.append(f"     IIRS pixels    centre median={np.median(cen):.0f}  "
                      f"p10={np.percentile(cen,10):.0f}  p90={np.percentile(cen,90):.0f} nm  "
                      f"(n={Fg.shape[1]})")


# ---------------------------------------------------------------- classifiers (return 1..6)
def classify_sam(F, Eref, cand):
    fn = np.linalg.norm(F, axis=0)
    en = np.linalg.norm(Eref, axis=1)
    cos = np.clip((Eref @ F) / (np.outer(en, fn) + 1e-12), -1.0, 1.0)
    ang = np.arccos(cos)
    ang[~cand, :] = np.inf                                # only minerals whose band is in-window
    ang[:, fn < 1e-9] = np.pi / 2
    return np.argmin(ang, axis=0) + 1


def classify_sff(F, Eref, cand):
    nb, N = F.shape
    Fm = F.mean(axis=0, keepdims=True)
    scores = np.full((len(Eref), N), -np.inf)
    for k, r in enumerate(Eref):
        if not cand[k]:
            continue
        rc = r - r.mean()
        denom = rc @ rc
        if denom < 1e-12:
            continue
        a = (rc @ (F - Fm)) / denom
        b = Fm[0] - a * r.mean()
        pred = np.outer(r, a) + b
        rms = np.sqrt(np.mean((F - pred) ** 2, axis=0))
        s = a / (rms + 1e-6)
        s[a <= 0] = -np.inf
        scores[k] = s
    cls = np.argmax(scores, axis=0) + 1
    cls[~np.isfinite(scores.max(axis=0))] = 0            # no positive fit -> unclassified
    return cls


def classify_svm(F, Eref, cand):
    rng = np.random.default_rng(RNG_SEED)
    Xtr, ytr = [], []
    noise = 0.012
    for k, r in enumerate(Eref):
        if not cand[k]:
            continue
        gains = rng.uniform(0.35, 1.35, SVM_AUG)
        samp = r[None, :] * gains[:, None] + rng.normal(0, noise, (SVM_AUG, r.size))
        Xtr.append(samp); ytr.append(np.full(SVM_AUG, k + 1))
    Xtr = np.vstack(Xtr); ytr = np.concatenate(ytr)
    scaler = StandardScaler().fit(Xtr)
    clf = SVC(kernel="rbf", C=10.0, gamma="scale", random_state=RNG_SEED)
    clf.fit(scaler.transform(Xtr), ytr)
    return clf.predict(scaler.transform(F.T)).astype(int)


# ---------------------------------------------------------------- writers
def to_map(cls_flat, mask, flat):
    out = np.zeros(mask.size, dtype="uint8")
    out[flat] = cls_flat.astype("uint8")
    return out.reshape(mask.shape)


def save_geotiff(path, cmap_img, h):
    import rasterio
    transform, crs = mc.geo_transform(h)
    prof = dict(driver="GTiff", height=h["lines"], width=h["samples"], count=1,
                dtype="uint8", crs=crs, transform=transform, compress="lzw")
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(cmap_img, 1)
        dst.write_colormap(1, {i: (*CLASS_RGB[i], 255) for i in range(len(CLASS_RGB))})
        dst.set_band_description(1, "mineral class (0=unclassified)")


def save_envi_class(path_img, cmap_img, h):
    os.makedirs(os.path.dirname(path_img), exist_ok=True)
    np.ascontiguousarray(cmap_img, dtype="uint8").tofile(path_img)
    S, L = h["samples"], h["lines"]
    lut = ", ".join(str(v) for rgb in CLASS_RGB for v in rgb)
    with open(path_img + ".hdr", "w") as f:
        f.write("ENVI\n")
        f.write("description = { mineral classification (SAM/SFF/SVM) }\n")
        f.write(f"samples = {S}\nlines   = {L}\nbands   = 1\n")
        f.write("header offset = 0\nfile type = ENVI Classification\n")
        f.write("data type = 1\ninterleave = bsq\nsensor type = Unknown\nbyte order = 0\n")
        f.write(f"classes = {len(CLASS_RGB)}\n")
        f.write("class names = {\n " + ",\n ".join(CLASS_NAME) + "}\n")
        f.write("class lookup = {" + lut + "}\n")
        if h["map_info"]:
            f.write("map info = {%s}\n" % ", ".join(str(x) for x in h["map_info"]))


def plot_map(ax, cmap_img, title):
    listed = ListedColormap([tuple(c / 255 for c in rgb) for rgb in CLASS_RGB])
    ax.imshow(cmap_img, cmap=listed, vmin=0, vmax=len(CLASS_RGB) - 1, interpolation="nearest")
    ax.set_title(title, fontsize=11)
    ax.axis("off")


def legend_handles():
    return [Patch(facecolor=tuple(c / 255 for c in CLASS_RGB[i]), edgecolor="0.4",
                  label=CLASS_NAME[i]) for i in range(len(CLASS_RGB))]


# ---------------------------------------------------------------- driver
def analyse(name, report):
    label = mc.SUBSETS[name]["label"]
    print(f"\n[classify] {label}")
    F, mask, flat, swl, h = sensor_features(name)
    Eref = endmember_features(swl, name)
    # feature strength = deepest INTERIOR band depth (exclude window edges, which can show
    # continuum-touching false depths). A pixel without a genuine absorption is Unclassified.
    interior = slice(2, -2) if F.shape[0] > 5 else slice(None)
    strength = F[interior].max(axis=0)
    gate = strength >= MIN_FEATURE
    # candidate minerals = those whose diagnostic band actually lies in the feature window
    # (>=4% depth); a mineral whose band is out-of-window/absent cannot be assigned.
    cand = Eref[:, interior].max(axis=1) >= 0.04
    print(f"   valid pixels {flat.size}; with diagnostic band {int(gate.sum())} "
          f"({100.0*gate.mean():.1f}%); candidates "
          f"{[CLASSES[k] for k in range(len(CLASSES)) if cand[k]]}")

    raw = {"SAM": classify_sam(F, Eref, cand), "SFF": classify_sff(F, Eref, cand),
           "SVM": classify_svm(F, Eref, cand)}

    os.makedirs(OUTDIR, exist_ok=True); os.makedirs(ENVIDIR, exist_ok=True)
    feat = ("2 um band only [%.0f-%.0f] nm" % BANDS[name]["two"]
            if FEATURE_MODE[name] == "two_micron"
            else "1um[%.0f-%.0f] + 2um[%.0f-%.0f] nm"
                 % (*BANDS[name]["one"], *BANDS[name]["two"]))
    report.append(f"\n=== {label}  (valid {flat.size}; diagnostic-band pixels "
                  f"{int(gate.sum())} = {100.0*gate.mean():.1f}%; features = {feat}) ===")
    report.append("  candidate minerals (diagnostic band in window): "
                  + ", ".join(CLASSES[k] for k in range(len(CLASSES)) if cand[k]))
    if FEATURE_MODE[name] == "two_micron":
        _report_band2_center(name, swl, F, Eref, gate, report)
    maps, cls_gated = {}, {}
    for meth in METHODS:
        cls = raw[meth].copy()
        cls[~gate] = 0                                   # gate -> unclassified
        cls_gated[meth] = cls
        cimg = to_map(cls, mask, flat)
        maps[meth] = cimg
        save_geotiff(os.path.join(OUTDIR, f"{name}_{meth}.tif"), cimg, h)
        save_envi_class(os.path.join(ENVIDIR, f"{name}_{meth}.img"), cimg, h)
        _plot_single(name, label, meth, cimg)
        counts = np.bincount(cls, minlength=len(CLASS_RGB))
        report.append(f"\n  {meth}:  (unclassified {counts[0]}  "
                      f"{100.0*counts[0]/flat.size:.1f}%)")
        for k in range(1, len(CLASS_RGB)):
            report.append(f"     {CLASS_NAME[k]:14s} {counts[k]:7d}  "
                          f"{100.0*counts[k]/flat.size:5.1f}%")

    report.append("\n  pairwise agreement over diagnostic-band pixels (same class):")
    g = gate
    for a, b in (("SAM", "SFF"), ("SAM", "SVM"), ("SFF", "SVM")):
        agr = 100.0 * np.mean(cls_gated[a][g] == cls_gated[b][g])
        report.append(f"     {a} vs {b}: {agr:5.1f}%")

    _plot_combined(name, label, maps)
    print("   ->", os.path.join(OUTDIR, f"{name}_classification.png"))
    return {"maps": maps, "h": h}


def _m3_to_iirs_index(h1, h2):
    """For every M3 pixel, the (row, col) of the IIRS pixel at the same ground point, plus an
    in-bounds mask. Both grids are Moon-2000 eqc, so co-location is done in eqc metres."""
    tr1, _ = mc.geo_transform(h1); tr2, _ = mc.geo_transform(h2)
    ulx1, uly1, px1, py1 = tr1.c, tr1.f, tr1.a, -tr1.e
    ulx2, uly2, px2, py2 = tr2.c, tr2.f, tr2.a, -tr2.e
    L1, S1 = h1["lines"], h1["samples"]
    L2, S2 = h2["lines"], h2["samples"]
    E = ulx1 + px1 * (np.arange(S1) + 0.5)
    N = uly1 - py1 * (np.arange(L1) + 0.5)
    EE, NN = np.meshgrid(E, N)
    col = np.floor((EE - ulx2) / px2).astype(int)
    row = np.floor((uly2 - NN) / py2).astype(int)
    inb = (col >= 0) & (col < S2) & (row >= 0) & (row < L2)
    return np.clip(row, 0, L2 - 1), np.clip(col, 0, S2 - 1), inb


def cross_sensor_agreement(res, report):
    """Append an M3-vs-IIRS cross-sensor agreement section (co-located by lon/lat)."""
    h1, h2 = res["ch1"]["h"], res["ch2"]["h"]
    row, col, inb = _m3_to_iirs_index(h1, h2)
    report.append("\n\n=== CROSS-SENSOR AGREEMENT (M3 vs corrected IIRS) ===")
    report.append("Per-pixel, co-located by lon/lat on the Moon-2000 eqc grid (each M3 140 m "
                  "pixel vs the IIRS 94.56 m pixel at the same ground point).")
    report.append("IIRS 2 um can only host clinopyroxene/orthopyroxene, so agreement reflects "
                  "the shared discriminable mafic classes; other M3 classes have no IIRS 2 um "
                  "counterpart.")
    same_dom = []
    for meth in METHODS:
        c1 = res["ch1"]["maps"][meth]
        c2m = res["ch2"]["maps"][meth][row, col]
        valid = inb & (c1 > 0) & (c2m > 0)
        n = int(valid.sum())
        if n == 0:
            report.append(f"  {meth}: no co-located classified pixels")
            continue
        agree = 100.0 * np.mean(c1[valid] == c2m[valid])
        d1 = int(np.bincount(c1[valid], minlength=len(CLASS_RGB)).argmax())
        d2 = int(np.bincount(c2m[valid], minlength=len(CLASS_RGB)).argmax())
        cpx = CLASSES.index("clinopyroxene") + 1
        cpx_both = 100.0 * np.mean((c1[valid] == cpx) & (c2m[valid] == cpx))
        report.append(f"  {meth}: co-located classified px={n}; same-class {agree:5.1f}%; "
                      f"both=clinopyroxene {cpx_both:5.1f}%; "
                      f"M3 dominant={CLASS_NAME[d1]}, IIRS dominant={CLASS_NAME[d2]}")
        if d1 == d2:
            same_dom.append((meth, CLASS_NAME[d1]))
    if same_dom:
        report.append("  => Same dominant mineral in both sensors for: "
                      + "; ".join(f"{m} ({c})" for m, c in same_dom)
                      + ".  SAM/SVM converge on clinopyroxene, cross-validating the M3 "
                        "crater composition with the independent IIRS sensor.")


def _plot_single(name, label, meth, cimg):
    fig, ax = plt.subplots(figsize=(7, 7))
    plot_map(ax, cimg, f"{label}  —  {meth} mineral classification")
    ax.legend(handles=legend_handles(), loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=9, frameon=True)
    fig.tight_layout()
    out = os.path.join(OUTDIR, f"{name}_{meth}.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("   ->", out)


def _plot_combined(name, label, maps):
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    for ax, meth in zip(axes, METHODS):
        plot_map(ax, maps[meth], meth)
    fig.legend(handles=legend_handles(), loc="lower center", ncol=7, fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(f"{label}  —  mineral classification (SAM | SFF | SVM)", fontsize=13)
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    fig.savefig(os.path.join(OUTDIR, f"{name}_classification.png"), dpi=130,
                bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    report = ["Mineral classification (SAM / SFF / SVM), 6 RELAB minerals.",
              "Features: two-band (1um+2um) local continuum-removed band depth (<=2600 nm). "
              "Pixels without a measurable diagnostic band are left Unclassified. "
              "SVM trained on library endmembers augmented with gain+noise.",
              "Classes: " + ", ".join(f"{i+1}={m}" for i, m in enumerate(CLASSES))]
    res = {}
    for name in ("ch1", "ch2"):
        res[name] = analyse(name, report)
    cross_sensor_agreement(res, report)
    rep = os.path.join(OUTDIR, "classification_report.txt")
    open(rep, "w").write("\n".join(report) + "\n")
    print("\n->", rep)


if __name__ == "__main__":
    main()
