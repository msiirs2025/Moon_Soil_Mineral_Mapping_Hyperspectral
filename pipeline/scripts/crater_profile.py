"""
Crater geomorphic-zone spectral profiles (paper Figure-3 methodology) for the CH1 (M3)
and CH2 (IIRS) ROI subsets.

Six points are sampled, one per geomorphic zone, and their reflectance spectra are
overlaid on a single "Spectral Profile" graph per sensor (as in the reference paper Fig 3):

    center        -> crater centre
    interior      -> crater floor, between centre and wall
    inner_edge    -> inner edge of the crater (inner wall)
    rim           -> crater edge / rim crest (boundary)
    outer_edge    -> outer edge of the crater (outer wall / proximal ejecta)
    outside_wall  -> outside the crater walls (background terrain)

Method (revised):
- The zone is set by RADIUS (fraction of the rim radius R); each zone is placed at a
  DIFFERENT AZIMUTH so the six points are DISTRIBUTED AROUND THE CRATER (different sectors),
  NOT strung along one straight radial line.
- Points are chosen ONCE on M3 (centre = image centre, rim R = outermost bright radial
  peak), restricted to pixels valid in BOTH sensors (M3 n IIRS overlap), then converted to
  lon/lat.  The SAME GROUND POINTS are sampled in both M3 and IIRS, so the two
  *_crater_points_on_fcc.png overlays coincide on the ground.
- Each point is a 3x3 average (valid pixels) for stability.  Spectra shown <= 2600 nm.

Outputs (pipeline/outputs/mineral_analysis/):
    ch{1,2}_crater_profile.png          Spectral Profile graph (6 zone curves)
    ch{1,2}_crater_points_on_fcc.png    the 6 points on the FCC (same ground pts both)
    crater_profile_points.csv           tidy long-format spectra
    crater_profile_report.txt           centre, rim, per-zone azimuth & coords

Run with the interpreter named by PYTHON_EXE in the repo-root .env file.
"""
import os, sys, csv
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.ndimage import uniform_filter1d
import mineral_common as mc
from fcc import compute_rgb

# (key, label, radius fraction of rim R, colour)
ZONES = [
    ("center",       "Center",          0.00, "#000000"),
    ("interior",     "Crater interior", 0.45, "#d62728"),
    ("inner_edge",   "Inner edge",      0.75, "#ff7f0e"),
    ("rim",          "Rim (boundary)",  1.00, "#2ca02c"),
    ("outer_edge",   "Outer edge",      1.25, "#1f77b4"),
    ("outside_wall", "Outside wall",    1.60, "#9467bd"),
]


def load_sensor(name):
    fill_lo = mc.SUBSETS[name]["fill_lo"]
    cube, h, _ = mc.read_subset(name)
    ci = mc.comparison_idx(cube, h, fill_lo)
    swl = h["wl"][ci]
    X = np.asarray(cube[ci], dtype="float32")
    mask = np.all(np.isfinite(X) & (X > fill_lo), axis=0)
    return dict(name=name, h=h, ci=ci, swl=swl, X=X, mask=mask)


def detect_center_rim(M, mask):
    """Crater centre (= image centre) and rim radius R (outermost prominent radial peak)."""
    L, S = M.shape
    cy, cx = L / 2.0, S / 2.0
    yy, xx = np.mgrid[0:L, 0:S]
    rr = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    rmax = int(min(cy, cx, L - cy, S - cx))
    prof = np.array([np.nanmean(np.where((rr >= r) & (rr < r + 1) & mask, M, np.nan))
                     for r in range(rmax)])
    profs = uniform_filter1d(np.nan_to_num(prof, nan=np.nanmin(prof)), 5)
    lo = int(0.25 * rmax)
    peaks = [lo + i for i in range(1, len(profs) - lo - 1)
             if profs[lo + i] > profs[lo + i - 1] and profs[lo + i] >= profs[lo + i + 1]]
    if peaks:
        base = np.nanmin(profs[lo:]); top = np.nanmax(profs[lo:])
        strong = [p for p in peaks if profs[p] - base >= 0.25 * (top - base + 1e-9)]
        rim = max(strong) if strong else max(peaks)
    else:
        rim = int(0.6 * rmax)
    return (cy, cx), rim, rmax


def place_zones(center, rim, sel_mask):
    """Place each zone at its radius but a DIFFERENT azimuth (distributed around the crater),
    every point valid in both sensors (sel_mask).  A global rotation offset is optimised for
    border clearance; each non-centre zone gets a distinct target sector (evenly spaced), with
    a small local nudge to land on a valid pixel while staying in its sector.
    Returns list of (key, label, color, frac, azimuth_deg, line, sample)."""
    cy, cx = center
    L, S = sel_mask.shape
    non_center = [z for z in ZONES if z[2] > 0]
    targets = np.linspace(0, 360, len(non_center), endpoint=False)   # even sectors

    def valid(l, s):
        return 0 <= l < L and 0 <= s < S and sel_mask[l, s]

    def try_offset(off):
        placed = {}
        for (key, zlabel, frac, color), tgt in zip(non_center, targets):
            hit = None
            for d in range(0, 60, 2):                        # stay near the sector
                for sgn in ((1, -1) if d else (1,)):
                    deg = (off + tgt + sgn * d) % 360
                    th = np.deg2rad(deg)
                    l = int(round(cy + frac * rim * np.sin(th)))
                    s = int(round(cx + frac * rim * np.cos(th)))
                    if valid(l, s):
                        hit = (deg, l, s); break
                if hit:
                    break
            if hit is None:
                return None
            placed[key] = hit
        # score = min border clearance across placed points (larger = safer)
        clear = min(min(l, L - 1 - l, s, S - 1 - s) for _, l, s in placed.values())
        return placed, clear

    best = None
    for off in range(0, 360, 5):
        r = try_offset(off)
        if r and (best is None or r[1] > best[1]):
            best = r
    if best is None:
        raise RuntimeError("could not place all crater zones on dual-valid pixels")
    placed = best[0]

    # centre point: image centre, nudged to the nearest dual-valid pixel if needed
    cl, cs = int(round(cy)), int(round(cx))
    if not valid(cl, cs):
        for rad in range(1, rim):
            cand = [(cl + dl, cs + ds) for dl in range(-rad, rad + 1)
                    for ds in range(-rad, rad + 1) if max(abs(dl), abs(ds)) == rad]
            cand = [(l, s) for l, s in cand if valid(l, s)]
            if cand:
                cl, cs = min(cand, key=lambda p: (p[0] - cy) ** 2 + (p[1] - cx) ** 2)
                break

    out = []
    for key, zlabel, frac, color in ZONES:
        if frac == 0:
            out.append((key, zlabel, color, frac, float("nan"), cl, cs))
        else:
            deg, l, s = placed[key]
            out.append((key, zlabel, color, frac, deg, l, s))
    return out


def sample_3x3(X, l, s, mask):
    """Mean of valid pixels in the 3x3 window around (l,s) -> spectrum (n,)."""
    L, S = mask.shape
    l0, l1 = max(0, l - 1), min(L, l + 2)
    s0, s1 = max(0, s - 1), min(S, s + 2)
    win = X[:, l0:l1, s0:s1]
    wm = mask[l0:l1, s0:s1]
    if wm.any():
        return win[:, wm].mean(1)
    return X[:, l, s]


def analyse_sensor(sd, zone_pts_lonlat, csv_rows, report):
    """Sample the shared crater ground points in one sensor and build its plots."""
    name, h, swl, X, mask = sd["name"], sd["h"], sd["swl"], sd["X"], sd["mask"]
    label = mc.SUBSETS[name]["label"]
    sensor = "M3" if name == "ch1" else "IIRS"
    print(f"\n[crater] {label}  (sampling shared points)")
    report.append(f"\n=== {label} ===  (spectra <= 2600 nm, {len(swl)} bands, 3x3 mean; "
                  f"points defined on M3, sampled here at identical lon/lat)")

    profile, points = {}, []
    for key, zlabel, color, frac, deg, lon, lat in zone_pts_lonlat:
        l, s = mc.lonlat_to_pixel(h, lon, lat)
        l = int(np.clip(l, 0, h["lines"] - 1))
        s = int(np.clip(s, 0, h["samples"] - 1))
        sp = sample_3x3(X, l, s, mask)
        profile[key] = (zlabel, color, sp)
        points.append((key, zlabel, color, l, s))
        az = "  --" if np.isnan(deg) else f"{deg:5.1f}"
        report.append(f"   {zlabel:16s} r={frac:.2f}R az={az} deg  "
                      f"lon={lon:8.4f} lat={lat:8.4f} -> sample={s:3d} line={l:3d}")
        for w, r in zip(swl, sp):
            csv_rows.append([name, label, key, zlabel, s, l, f"{deg:.1f}",
                             f"{lon:.5f}", f"{lat:.5f}", f"{w:.2f}", f"{r:.6f}"])

    _plot_profile(name, sensor, swl, profile)
    _plot_points(name, label, compute_rgb(name)[0], points)


def _plot_profile(name, sensor, swl, profile):
    fig, ax = plt.subplots(figsize=mc.PROFILE_FIGSIZE)
    for key, zlabel, frac, color in ZONES:
        _, col, sp = profile[key]
        ax.plot(swl, sp, lw=1.6, color=col, label=zlabel)
    ax.set_xlim(*mc.PROFILE_XLIM)
    ax.set_ylim(*mc.profile_ylim([profile[k][2] for k in profile]))
    ax.xaxis.set_major_locator(MultipleLocator(mc.PROFILE_XTICK))
    ax.yaxis.set_major_locator(MultipleLocator(mc.PROFILE_YTICK))
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Reflectance")
    ax.set_title(f"Spectral Profile ({sensor})")
    leg = ax.legend(loc="upper left", fontsize=9, frameon=True, framealpha=0.9)
    for txt, (key, zlabel, frac, color) in zip(leg.get_texts(), ZONES):
        txt.set_color(color)
    fig.tight_layout()
    out = os.path.join(mc.ANALYSIS, f"{name}_crater_profile.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("   ->", out)


def _plot_points(name, label, rgb, points):
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.imshow(rgb)
    # centre + rim circle for geomorphic context (no straight transect line)
    ck = next(p for p in points if p[0] == "center")
    cl, cs = ck[3], ck[4]
    rim_zone = next(p for p in points if p[0] == "rim")
    R = np.hypot(rim_zone[3] - cl, rim_zone[4] - cs)
    ax.add_patch(plt.Circle((cs, cl), R, fill=False, color="white", lw=0.8,
                            alpha=0.5, ls="--", zorder=1))
    for key, zlabel, color, l, s in points:
        ax.scatter([s], [l], s=90, marker="o", facecolors="none",
                   edgecolors=color, linewidths=2.0, zorder=3)
        ax.annotate(zlabel, (s, l), textcoords="offset points", xytext=(6, 4),
                    fontsize=8, color=color)
    ax.set_title(f"{label}  —  crater zone points (distributed around the crater;\n"
                 f"same ground locations in M3 & IIRS)")
    ax.set_xlabel("sample"); ax.set_ylabel("line")
    fig.tight_layout()
    out = os.path.join(mc.ANALYSIS, f"{name}_crater_points_on_fcc.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("   ->", out)


def main():
    os.makedirs(mc.ANALYSIS, exist_ok=True)
    m3 = load_sensor("ch1")
    iirs = load_sensor("ch2")

    # select once on M3, dual-valid with IIRS
    sel_mask = mc.dual_valid_overlap(m3["h"], m3["mask"], iirs["h"], iirs["mask"])
    M = np.where(m3["mask"], m3["X"].mean(0), np.nan)
    center, rim, rmax = detect_center_rim(M, m3["mask"])
    zone_pts = place_zones(center, rim, sel_mask)
    print(f"[crater] center=({center[0]:.0f},{center[1]:.0f})  rim R={rim}px "
          f"(M3); zones distributed around crater")

    # convert M3 pixels -> lon/lat (the shared ground points)
    zone_pts_lonlat = []
    for key, zlabel, color, frac, deg, l, s in zone_pts:
        _, _, lon, lat = mc.pixel_lonlat(m3["h"], l, s)
        zone_pts_lonlat.append((key, zlabel, color, frac, deg, lon, lat))

    csv_rows, report = [], ["Crater geomorphic-zone spectral profiles (paper Fig-3 method)"]
    report.append("Points selected ONCE on M3 (centre=image centre, rim R=outermost bright "
                  "radial peak), one per zone at a DIFFERENT azimuth (distributed around the "
                  "crater), restricted to M3 n IIRS dual-valid pixels; identical lon/lat "
                  "sampled in both sensors.")
    report.append(f"M3 centre=(line {center[0]:.0f}, sample {center[1]:.0f})  rim R={rim}px")

    for sd in (m3, iirs):
        analyse_sensor(sd, zone_pts_lonlat, csv_rows, report)

    csv_path = os.path.join(mc.ANALYSIS, "crater_profile_points.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sensor", "label", "zone", "zone_label", "sample", "line",
                    "azimuth_deg", "lon_deg", "lat_deg", "wavelength_nm", "reflectance"])
        w.writerows(csv_rows)
    rep = os.path.join(mc.ANALYSIS, "crater_profile_report.txt")
    open(rep, "w").write("\n".join(report) + "\n")
    print("\n->", csv_path, f"({len(csv_rows)} rows)\n->", rep)


if __name__ == "__main__":
    main()
