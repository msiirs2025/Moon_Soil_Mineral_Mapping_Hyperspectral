"""
RELAB-referenced mineral spectral signatures for the CH1 (M3) and CH2 (IIRS) ROI subsets,
following Figure 3 of the reference paper (FCC + sampled-point spectra).

Point selection (scientifically-correct, library-INDEPENDENT):
  * Sample points are chosen by each mineral's diagnostic ABSORPTION BAND DEPTHS measured
    directly from the observed cube (1 um / 2 um position & strength, 1.25 um for
    plagioclase, dark+featureless for ilmenite) -- NOT by matching the RELAB library.
    The RELAB spectra are used ONLY for the plotted shape comparison, never for selection.
  * Selection is done ONCE on M3 (CH1) -- it carries the full 1 um+2 um mafic diagnostic
    range with cleaner bands -- restricted to the M3 n IIRS geographic overlap and to
    pixels that are valid in BOTH sensors.
  * The chosen pixels are converted to lon/lat (Moon-2000 eqc) and the SAME GROUND POINTS
    are then sampled in both M3 and IIRS.  Identical location in the two datasets.

For each sensor x each of 6 minerals:
  1. Sample the 4 shared points -> observed spectrum (<= 2600 nm, the RELAB limit).
  2. Plot the 4 observed spectra + the RELAB reference (scaled to the observed continuum),
     both as reflectance and continuum-removed.
  3. Mark the points on the FCC.

Outputs (pipeline/outputs/mineral_analysis/):
  <name>_<mineral>.png            per-mineral comparison plot        (6 x 2)
  <name>_points_on_fcc.png        all 24 points on the FCC           (2)
  extracted_spectra.csv           tidy long-format spectra           (1)
  extracted_spectra_<name>.sli    ENVI spectral library per sensor   (2)
  match_report.txt                per-mineral selection summary      (1)

Run with the interpreter named by PYTHON_EXE in the repo-root .env file.
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import mineral_common as mc
from fcc import compute_rgb

N_POINTS = 4
MIN_SEP = 3


def load_sensor(name):
    """Read a subset and derive everything used for selection/sampling."""
    fill_lo = mc.SUBSETS[name]["fill_lo"]
    cube, h, _ = mc.read_subset(name)
    ci = mc.comparison_idx(cube, h, fill_lo)              # bands <= 2600 nm, not fill
    swl = h["wl"][ci]
    X = np.asarray(cube[ci], dtype="float32")            # (n, L, S) reflectance
    mask = np.all(np.isfinite(X) & (X > fill_lo), axis=0)
    return dict(name=name, h=h, fill_lo=fill_lo, ci=ci, swl=swl, X=X, mask=mask)


def select_points_on_m3(m3, other):
    """Diagnostic-band-depth point selection on M3, restricted to the M3 n IIRS overlap
    (and dual-valid).  Returns mineral -> list of dicts(line,sample,lon,lat,score)."""
    print("[select] diagnostic band-depth selection on CH1/M3 (shared with IIRS)")
    print("   continuum removal (hull) ...")
    Xcr = mc.continuum_remove_cube(m3["X"], m3["swl"], m3["mask"])
    sel_mask = mc.dual_valid_overlap(m3["h"], m3["mask"], other["h"], other["mask"])
    print(f"   candidate pixels (valid in both, inside overlap): {int(sel_mask.sum())}")
    scores = mc.diagnostic_scores(m3["X"], Xcr, m3["swl"], sel_mask)

    shared = {}
    for mineral in mc.MINERAL_ORDER:
        pts = mc.pick_top(scores[mineral], sel_mask, n=N_POINTS, min_sep=MIN_SEP)
        recs = []
        for (l, s) in pts:
            _, _, lon, lat = mc.pixel_lonlat(m3["h"], l, s)
            recs.append(dict(line=l, sample=s, lon=lon, lat=lat,
                             score=float(scores[mineral][l, s])))
        shared[mineral] = recs
        print(f"   {mineral:14s} pts(M3 s,l) {[(r['sample'], r['line']) for r in recs]}"
              f"  score {np.mean([r['score'] for r in recs]):+.3f}")
    return shared


def analyse_sensor(sd, shared, csv_rows, report_lines):
    """Sample the shared ground points in one sensor and build its plots/spectra."""
    name, h, swl, X, mask = sd["name"], sd["h"], sd["swl"], sd["X"], sd["mask"]
    label = mc.SUBSETS[name]["label"]
    print(f"\n[spectra] {label}  (sampling shared points)")
    Xcr = mc.continuum_remove_cube(X, swl, mask)

    rgb, _, _, _ = compute_rgb(name)
    overlay_pts, profile = {}, {}
    lib_spectra, lib_names = [], []

    report_lines.append(f"\n=== {label}  (comparison window {swl[0]:.0f}-{swl[-1]:.0f} nm, "
                        f"{len(swl)} bands; points defined on M3 by diagnostic band depths, "
                        f"sampled here at the identical lon/lat) ===")

    for mineral in mc.MINERAL_ORDER:
        recs = shared[mineral]
        obs, obs_cr, pix, m3score = [], [], [], []
        report_lines.append(f"\n {mineral}  ({mc.RELAB_MAP[mineral]})  "
                            f"mean M3 diag score = "
                            f"{np.mean([r['score'] for r in recs]):+.3f}")
        for k, r in enumerate(recs):
            l, s = mc.lonlat_to_pixel(h, r["lon"], r["lat"])
            l = int(np.clip(l, 0, h["lines"] - 1))
            s = int(np.clip(s, 0, h["samples"] - 1))
            pix.append((l, s))
            sp = X[:, l, s]
            obs.append(sp)
            obs_cr.append(Xcr[:, l, s])
            m3score.append(r["score"])
            report_lines.append(f"   P{k+1}: lon={r['lon']:8.4f} lat={r['lat']:8.4f}  "
                                f"-> sample={s:3d} line={l:3d}  M3 score={r['score']:+.3f}")
            for w, rf in zip(swl, sp):
                csv_rows.append([name, label, mineral, k + 1, s, l,
                                 f"{r['lon']:.5f}", f"{r['lat']:.5f}", f"{r['score']:.4f}",
                                 f"{w:.2f}", f"{rf:.6f}"])
            lib_spectra.append(sp)
            lib_names.append(f"{name}_{mineral}_P{k+1}")
        overlay_pts[mineral] = pix

        # reference scaled to the observed continuum (shape comparison only)
        mean_obs = np.mean(obs, axis=0)
        profile[mineral] = mean_obs
        ref_wl, ref_rf = mc.load_relab(mineral)
        ref = mc.resample_ref(ref_wl, ref_rf, swl)
        ref_cr = mc.continuum_remove_1d(swl, ref)
        scale = np.median(mean_obs) / (np.median(ref) + 1e-12)
        lib_spectra.append(ref * scale)
        lib_names.append(f"{name}_{mineral}_RELAB_ref_scaled")

        _plot_mineral(name, label, mineral, swl, obs, obs_cr, m3score, ref * scale, ref_cr)
        print(f"   {mineral:14s} pts(s,l) {[(s, l) for (l, s) in pix]}")

    _plot_overlay(name, label, rgb, overlay_pts)
    _plot_profile(name, profile, swl)

    sli = os.path.join(mc.ANALYSIS, f"extracted_spectra_{name}.sli")
    mc.write_speclib(sli, swl, np.array(lib_spectra), lib_names)
    print("   ->", sli)


def _plot_mineral(name, label, mineral, swl, obs, obs_cr, score, ref_scaled, ref_cr):
    color = mc.MINERAL_COLOR[mineral]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # left: reflectance (observed + RELAB reference scaled to observed continuum)
    for k, sp in enumerate(obs):
        ax1.plot(swl, sp, lw=1.2, color=color, alpha=0.45 + 0.14 * k, label=f"P{k+1}")
    ax1.plot(swl, ref_scaled, "k--", lw=1.8, label="RELAB ref (scaled)")
    ax1.set_ylabel("Reflectance")
    ax1.set_title("Reflectance")

    # right: continuum-removed (diagnostic absorption form)
    for k, sp in enumerate(obs_cr):
        ax2.plot(swl, sp, lw=1.2, color=color, alpha=0.45 + 0.14 * k, label=f"P{k+1}")
    ax2.plot(swl, ref_cr, "k--", lw=1.8, label="RELAB ref")
    ax2.set_ylabel("Continuum-removed reflectance")
    ax2.set_title("Continuum-removed (diagnostic bands)")

    for ax in (ax1, ax2):
        for xln in mc.ABS_LINES:
            if swl[0] <= xln <= swl[-1]:
                ax.axvline(xln, color="grey", ls=":", lw=0.8)
        ax.set_xlim(swl[0], mc.CMP_MAX_NM)
        ax.set_xlabel("Wavelength (nm)")
        ax.legend(fontsize=8, ncol=2)
        ax.grid(alpha=0.25)

    fig.suptitle(f"{label}  —  {mineral}   (4 band-depth-selected points, same lon/lat "
                 f"in M3 & IIRS; mean M3 diag score {np.mean(score):+.3f})", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(mc.ANALYSIS, f"{name}_{mineral}.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)


def _plot_profile(name, profile, swl):
    """Paper Figure-3 style 'Spectral Profile' graph: one reflectance curve per mineral
    (mean of its 4 shared points), all overlaid on a single axis with a legend box."""
    sensor = "M3" if name == "ch1" else "IIRS"
    fig, ax = plt.subplots(figsize=mc.PROFILE_FIGSIZE)
    for mineral in mc.MINERAL_ORDER:
        ax.plot(swl, profile[mineral], lw=1.6, color=mc.MINERAL_COLOR[mineral],
                label=mineral)
    ax.set_xlim(*mc.PROFILE_XLIM)
    ax.set_ylim(*mc.profile_ylim([profile[m] for m in mc.MINERAL_ORDER]))
    ax.xaxis.set_major_locator(MultipleLocator(mc.PROFILE_XTICK))
    ax.yaxis.set_major_locator(MultipleLocator(mc.PROFILE_YTICK))
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance")
    ax.set_title(f"Spectral Profile ({sensor})")
    leg = ax.legend(loc="upper left", fontsize=9, frameon=True, framealpha=0.9)
    for txt, mineral in zip(leg.get_texts(), mc.MINERAL_ORDER):
        txt.set_color(mc.MINERAL_COLOR[mineral])          # coloured labels, as in the paper
    fig.tight_layout()
    out = os.path.join(mc.ANALYSIS, f"{name}_mineral_spectral_profile.png")
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print("   ->", out)


def _plot_overlay(name, label, rgb, overlay_pts):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(rgb)
    for mineral, pts in overlay_pts.items():
        ys = [l for (l, s) in pts]
        xs = [s for (l, s) in pts]
        ax.scatter(xs, ys, s=70, marker="s", facecolors="none",
                   edgecolors=mc.MINERAL_COLOR[mineral], linewidths=1.6, label=mineral)
    ax.set_title(f"{label}  —  FCC with band-depth-selected mineral points\n"
                 f"(same ground locations in M3 & IIRS)")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.8)
    ax.set_xlabel("sample"); ax.set_ylabel("line")
    fig.tight_layout()
    out = os.path.join(mc.ANALYSIS, f"{name}_points_on_fcc.png")
    fig.savefig(out, dpi=130)
    plt.close(fig)


def main():
    os.makedirs(mc.ANALYSIS, exist_ok=True)
    m3 = load_sensor("ch1")
    iirs = load_sensor("ch2")
    shared = select_points_on_m3(m3, iirs)

    csv_rows, report_lines = [], []
    report_lines.append("Mineral point-selection report")
    report_lines.append("Selection: per-mineral DIAGNOSTIC BAND DEPTHS measured on M3 "
                        "(library-independent), restricted to the M3 n IIRS overlap and "
                        "dual-valid pixels; the SAME lon/lat are sampled in both sensors.")
    report_lines.append("All comparison / plotting capped at 2600 nm (RELAB library limit).")

    for sd in (m3, iirs):
        analyse_sensor(sd, shared, csv_rows, report_lines)

    csv_path = os.path.join(mc.ANALYSIS, "extracted_spectra.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sensor", "label", "mineral", "point", "sample", "line",
                    "lon_deg", "lat_deg", "m3_diag_score", "wavelength_nm", "reflectance"])
        w.writerows(csv_rows)
    print("\n->", csv_path, f"({len(csv_rows)} rows)")

    rep = os.path.join(mc.ANALYSIS, "match_report.txt")
    open(rep, "w").write("\n".join(report_lines) + "\n")
    print("->", rep)


if __name__ == "__main__":
    main()
