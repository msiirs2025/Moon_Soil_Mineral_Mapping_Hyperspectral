"""
Band-Shaping-Algorithm band-parameter images (BR, BT, BS, BC) for the CH1 (M3) and
CH2 (IIRS) ROI subsets, reproducing Table 2 of the reference paper
"Evaluation of Chandrayaan-2 IIRS Data for Lunar Surface Exploration" (Singh & Kumar, IIRS).

Table 2 (B_x = reflectance at x nm; nearest band used per sensor):

    Parameter            IIRS                         M3
    Band Ratio    (BR)   B2026 / B1015                B2018 / B1009
    Band Tilt     (BT)   B914  / B1015                B910  / B1009
    Band Strength (BS)   B1015 / B847                 B1009 / B750
    Band Curvature(BC)   B847/B914 + B1015/B914       B750/B910 + B1009/B910

These four parameters characterise the 1 um / 2 um ferrous absorptions and, per the paper,
are combined into an enhanced FCC to discriminate mafic minerals.

Outputs (pipeline/outputs/band_parameters/):
    ch{1,2}_{BR,BT,BS,BC}.tif      georeferenced float32 parameter image (Moon-2000 eqc)
    ch{1,2}_{BR,BT,BS,BC}.png      quick-look (2-98% stretch, colourbar)
    ch{1,2}_band_parameters.png    combined 2x2 panel per sensor
    ch{1,2}_BSA_fcc.png            enhanced FCC  R=BR, G=BT, B=BS  (PNG only; channels are
                                   the single-band BR/BT/BS GeoTIFFs -> no RGB GeoTIFF)
    envi/ch{1,2}_{BR,BT,BS,BC}.img+.hdr   ENVI single-band float32 (georeferenced)
    band_parameters_report.txt     bands actually used per parameter

Run with the interpreter named by PYTHON_EXE in the repo-root .env file.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mineral_common as mc
import common as cm

OUTDIR = os.path.join(mc.OUT, "band_parameters")
ENVIDIR = os.path.join(OUTDIR, "envi")            # ENVI .img + .hdr live here

# target wavelengths (nm) per sensor, from Table 2
TARGETS = {
    "ch1": {"cont": 750.0, "short": 910.0, "one": 1009.0, "two": 2018.0},   # M3
    "ch2": {"cont": 847.0, "short": 914.0, "one": 1015.0, "two": 2026.0},   # IIRS
}

# parameter order, display label, and formula string per sensor (filled at runtime)
PARAM_ORDER = ["BR", "BT", "BS", "BC"]
PARAM_LABEL = {"BR": "Band Ratio (BR)", "BT": "Band Tilt (BT)",
               "BS": "Band Strength (BS)", "BC": "Band Curvature (BC)"}
PARAM_CMAP = {p: "gray" for p in PARAM_ORDER}     # single grey colour for every parameter


def compute_params(R):
    """R = dict of reflectance arrays keyed cont/short/one/two -> dict of 4 parameter maps."""
    eps = 1e-6
    def d(a, b):
        return a / np.where(np.abs(b) < eps, np.nan, b)
    return {
        "BR": d(R["two"], R["one"]),                        # 2um / 1um
        "BT": d(R["short"], R["one"]),                      # short shoulder / 1um
        "BS": d(R["one"], R["cont"]),                       # 1um / continuum
        "BC": d(R["cont"], R["short"]) + d(R["one"], R["short"]),
    }


def stretch(a, mask, lo=2, hi=98):
    p1, p2 = np.percentile(a[mask], [lo, hi])
    return np.clip((a - p1) / (p2 - p1 + 1e-12), 0.0, 1.0)


def save_geotiff(path, arr, h):
    import rasterio
    transform, crs = mc.geo_transform(h)
    prof = dict(driver="GTiff", height=h["lines"], width=h["samples"], count=1,
                dtype="float32", crs=crs, transform=transform, nodata=np.nan,
                compress="lzw")
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(arr.astype("float32"), 1)


def analyse(name, report):
    label = mc.SUBSETS[name]["label"]
    sensor = "M3" if name == "ch1" else "IIRS"
    fill_lo = mc.SUBSETS[name]["fill_lo"]
    print(f"\n[band-param] {label}")
    cube, h, _ = mc.read_subset(name)
    wl = h["wl"]
    tgt = TARGETS[name]

    bidx = {k: mc.nearest_band(wl, v) for k, v in tgt.items()}
    R = {k: np.asarray(cube[b], dtype="float32") for k, b in bidx.items()}
    used = np.array(list(bidx.values()))
    mask = mc.valid_mask(cube, fill_lo, used)

    report.append(f"\n=== {label} ===  (nearest band to each Table-2 wavelength)")
    for k in ("cont", "short", "one", "two"):
        report.append(f"   {k:6s} target {tgt[k]:7.1f} nm -> band {bidx[k]+1:3d} "
                      f"= {wl[bidx[k]]:7.1f} nm")

    params = compute_params(R)
    # formula strings (actual nearest-band wavelengths, rounded)
    w = {k: round(wl[bidx[k]]) for k in bidx}
    formulas = {
        "BR": f"B{w['two']} / B{w['one']}",
        "BT": f"B{w['short']} / B{w['one']}",
        "BS": f"B{w['one']} / B{w['cont']}",
        "BC": f"B{w['cont']}/B{w['short']} + B{w['one']}/B{w['short']}",
    }

    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(ENVIDIR, exist_ok=True)
    for p in PARAM_ORDER:
        arr = np.where(mask, params[p], np.nan)
        save_geotiff(os.path.join(OUTDIR, f"{name}_{p}.tif"), arr, h)
        # ENVI .img + .hdr (float32 BSQ, georeferenced via the subset map info)
        cm.write_envi(os.path.join(ENVIDIR, f"{name}_{p}.img"),
                      arr[None, :, :], band_names=[f"{p} {formulas[p]}"],
                      data_ignore="nan", map_info=h["map_info"],
                      description=f"{label} {PARAM_LABEL[p]} = {formulas[p]}")
        _plot_single(name, sensor, p, arr, mask, formulas[p])
        vals = arr[mask & np.isfinite(arr)]
        report.append(f"   {PARAM_LABEL[p]:20s} = {formulas[p]:34s} "
                      f"range [{np.nanpercentile(vals,2):.3f}, {np.nanpercentile(vals,98):.3f}]")

    _plot_panel(name, label, sensor, params, mask, formulas)
    _enhanced_fcc(name, label, params, mask, h)


def _plot_single(name, sensor, p, arr, mask, formula):
    fig, ax = plt.subplots(figsize=(6.2, 6.6))
    disp = np.where(mask, arr, np.nan)
    v1, v2 = np.nanpercentile(disp[mask], [2, 98])
    ax.imshow(disp, cmap=PARAM_CMAP[p], vmin=v1, vmax=v2)
    ax.set_title(f"{sensor}  —  {PARAM_LABEL[p]}\n{formula}", fontsize=11)
    ax.axis("off")
    fig.tight_layout()
    out = os.path.join(OUTDIR, f"{name}_{p}.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print("   ->", out)


def _plot_panel(name, label, sensor, params, mask, formulas):
    fig, axes = plt.subplots(2, 2, figsize=(11, 11))
    for ax, p in zip(axes.ravel(), PARAM_ORDER):
        disp = np.where(mask, params[p], np.nan)
        v1, v2 = np.nanpercentile(disp[mask], [2, 98])
        ax.imshow(disp, cmap=PARAM_CMAP[p], vmin=v1, vmax=v2)
        ax.set_title(f"{PARAM_LABEL[p]}\n{formulas[p]}", fontsize=10)
        ax.axis("off")
    fig.suptitle(f"{label}  —  Band-Shaping-Algorithm parameters (paper Table 2)",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(OUTDIR, f"{name}_band_parameters.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print("   ->", out)


def _enhanced_fcc(name, label, params, mask, h):
    """Enhanced false-colour composite from the band parameters (R=BR, G=BT, B=BS),
    each 2-98% stretched -- the paper's mineral-discrimination map. PNG quick-look only;
    the three channels are already exported as single-band BR/BT/BS GeoTIFFs, so no RGB
    GeoTIFF is written (every .tif stays single-band)."""
    rgb = np.zeros((h["lines"], h["samples"], 3), dtype="float32")
    for i, p in enumerate(("BR", "BT", "BS")):
        rgb[..., i] = stretch(params[p], mask)
    rgb[~mask] = 0.0
    rgb8 = (rgb * 255).astype("uint8")
    plt.imsave(os.path.join(OUTDIR, f"{name}_BSA_fcc.png"), rgb8)
    print(f"   -> {name}_BSA_fcc.png  (R=BR, G=BT, B=BS; channels are the BR/BT/BS GeoTIFFs)")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    report = ["Band-Shaping-Algorithm parameters (reference paper Table 2)",
              "BR=Band Ratio, BT=Band Tilt, BS=Band Strength, BC=Band Curvature.",
              "Enhanced FCC channel map: R=BR, G=BT, B=BS."]
    for name in ("ch1", "ch2"):
        analyse(name, report)
    rep = os.path.join(OUTDIR, "band_parameters_report.txt")
    open(rep, "w").write("\n".join(report) + "\n")
    print("\n->", rep)


if __name__ == "__main__":
    main()
