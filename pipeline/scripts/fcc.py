"""
Mafic-diagnostic False-Colour Composite (FCC) for the CH1 (M3) and CH2 (IIRS) ROI subsets.

Band assignment (nearest band to each target, chosen per sensor):
    R ~ 2.0 um  (pyroxene 2 um band)
    G ~ 1.0 um  (mafic 1 um band)
    B ~ shortest / continuum
Each channel gets an independent 2-98% linear stretch over valid pixels.

Outputs (pipeline/outputs/mineral_analysis/):
    <name>_fcc.tif   georeferenced 3-band uint8 RGB (Moon 2000 eqc)
    <name>_fcc.png   quick-look

compute_rgb() is reused by mineral_spectra.py to draw the sample points on the FCC.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mineral_common as mc


def _stretch(band, mask, lo=2, hi=98):
    """2-98% linear stretch to [0,1] using only valid pixels."""
    p1, p2 = np.percentile(band[mask], [lo, hi])
    return np.clip((band - p1) / (p2 - p1 + 1e-12), 0.0, 1.0)


def compute_rgb(name):
    """Return (rgb_uint8[L,S,3], mask[L,S], (bR,bG,bB), hdr) for a subset."""
    cube, h, nod = mc.read_subset(name)
    wl = h["wl"]
    fill_lo = mc.SUBSETS[name]["fill_lo"]
    tR, tG, tB = mc.FCC_TARGETS[name]
    bR, bG, bB = (mc.nearest_band(wl, t) for t in (tR, tG, tB))
    mask = mc.valid_mask(cube, fill_lo, np.array([bR, bG, bB]))
    rgb = np.zeros((h["lines"], h["samples"], 3), dtype="float32")
    for i, b in enumerate((bR, bG, bB)):
        rgb[..., i] = _stretch(np.asarray(cube[b], dtype="float32"), mask)
    rgb[~mask] = 0.0
    return (rgb * 255).astype("uint8"), mask, (bR, bG, bB), h


def build(name):
    print(f"[FCC] {mc.SUBSETS[name]['label']}")
    rgb, mask, (bR, bG, bB), h = compute_rgb(name)
    wl = h["wl"]
    print(f"   R band {bR+1:3d} = {wl[bR]:7.1f} nm")
    print(f"   G band {bG+1:3d} = {wl[bG]:7.1f} nm")
    print(f"   B band {bB+1:3d} = {wl[bB]:7.1f} nm")

    os.makedirs(mc.ANALYSIS, exist_ok=True)
    png = os.path.join(mc.ANALYSIS, f"{name}_fcc.png")
    tif = os.path.join(mc.ANALYSIS, f"{name}_fcc.tif")

    # PNG quick-look
    plt.imsave(png, rgb)

    # georeferenced RGB GeoTIFF
    import rasterio
    transform, crs = mc.geo_transform(h)
    prof = dict(driver="GTiff", height=h["lines"], width=h["samples"], count=3,
                dtype="uint8", crs=crs, transform=transform, photometric="RGB",
                compress="lzw")
    with rasterio.open(tif, "w", **prof) as dst:
        for i in range(3):
            dst.write(rgb[..., i], i + 1)
        dst.set_band_description(1, f"R {wl[bR]:.1f} nm")
        dst.set_band_description(2, f"G {wl[bG]:.1f} nm")
        dst.set_band_description(3, f"B {wl[bB]:.1f} nm")
    print("   ->", png, "\n   ->", tif)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for name in (["ch1", "ch2"] if which == "all" else [which]):
        build(name)
