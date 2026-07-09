# Chandrayaan lunar hyperspectral processing pipeline

Two channels processed independently and delivered as **ENVI (`.img`+`.hdr`) and GeoTIFF,
all bands, for every step**.

- **CH2 = Chandrayaan‑2 IIRS**, treated as **Level‑1 radiance** (`radiance = DN × 0.01`,
  per the ISRO CH2IIRS convention): radiometric calibration → thermal removal →
  photometric correction → selenoreferencing → destriping.
- **CH1 = Chandrayaan‑1 M3**, Level‑2 reflectance: **selenoreferencing only**.

All outputs live under `pipeline/outputs/`. Scripts under `pipeline/scripts/`.
Run with `C:\Users\HP\anaconda3\python.exe`.

---

## CH2 — IIRS (250 samples × 13569 lines × 256 bands, BSQ uint16)

Method reproduces the ISRO reference **CH2IIRS** QGIS plugin (Verma, Chauhan & Chauhan,
2022, *Icarus* 383, 115075) — the same tool that produced the `_Temp`/`_angle`/`_inc_corrRef`
files already in the data folder. Reproduction verified against those files
(incidence angle: exact; temperature: mean |Δ| ≈ 0.1 K; reflectance: sub‑percent).

| Step | Script | Output folder | Bands |
|---|---|---|---|
| 1 Radiometric calibration (radiance = DN×0.01, mW/cm²/sr/µm) | `ch2_pipeline.py 1` | `step1_radiance` | 256 |
| aux Incidence angle (from `.spm`, 90−sun_elev, per line) | `ch2_pipeline.py angle` | `aux` | 1 |
| 2 Thermal removal → I/F reflectance | `ch2_pipeline.py 2` | `step2_thermal` | 247 |
| 3 Photometric correction (Lommel‑Seeliger → i=30,e=0,g=30) | `ch2_pipeline.py 3` | `step3_photometric` | 247 |
| 4 Selenoreferencing (Moon2000 eqc, 94.56 m) | `ch2_seleno_destripe.py 4` | `step4_seleno` | 247 |
| 5 Destriping (per‑detector) + selenoreference | `ch2_seleno_destripe.py 5` | `step5_destripe` | 247 |

**Thermal removal / reflectance** (exact reference math, per pixel):
- Brightness temperature by inverting Planck over 4500–4874 nm with emissivity ε=0.95:
  `T = hc / (λ k · ln(ε·2hc²·1e‑6 / (L·λ⁵) + 1))`, averaged over those bands.
- Thermal radiance `B(λ,T)` (Planck) removed and converted to I/F reflectance:
  `ref = π·(L − ε·B) / E_sun(λ)`, then divided by the 247‑band empirical spectral
  correction and 3‑point spectrally smoothed. `E_sun` from CH2IIRS `Solar flux.txt`.
- Bands 1–7 and 255–256 are dropped (reference convention) → **247 bands, 847–4993 nm**.

**Photometric correction:** Lommel‑Seeliger disk function `D(i,e)=cos i/(cos i+cos e)`,
nadir emission (e≈0), normalised to standard geometry:
`ref_std = ref · (cos i + 1)/cos i · D(30°,0)`  with `D(30°,0)=0.4641`.
(The empirical phase‑curve term of Verma et al. 2023 needs ISRO's non‑public polynomial
fits; the analytic disk normalisation is applied instead.)

**Selenoreferencing:** 4 PDS4 corner coordinates → GCP warp (bilinear) to
**Equirectangular Moon 2000** (sphere R=1737.4 km). This is a first‑order georeference;
sub‑pixel accuracy would need a SPICE‑derived per‑pixel geometry backplane.

**Destriping:** per‑detector (per cross‑track sample column) moment matching, done in
**sensor geometry** (the correct domain — stripes are detector‑aligned) then georeferenced.
The requested order lists destriping last; because striping must be removed before warping,
`step4` is the georeferenced photometric product and `step5` is the destriped product
(sensor + georeferenced) as the definitive final cube.

---

## CH1 — M3 (304 samples × 29939 lines × 85 bands, BIL float32)

Selenoreferencing only (data is already L2 reflectance, 460–3000 nm).

- **Method:** rigorous per‑pixel **geolocation‑array** warp (`gdalwarp -geoloc`) using the
  LOC backplane longitude/latitude of *every* pixel → Equirectangular Moon 2000, 140 m/pixel.
- **Polar tail:** this is a near‑polar orbital strip (55°N → ~87°S). The last ~2500 lines
  cross immediately beside the south pole, where longitude is degenerate for an
  equirectangular grid; lines with lat < −80° are **excluded** from the eqc product (eqc is
  singular near the pole). They remain in the sensor‑geometry data; a south‑polar‑
  stereographic tail can be generated on request.
- Output: `pipeline/outputs/ch1/seleno/m3g20090731t045352_seleno_rfl.{tif,img}`.

---

## Mineral analysis of the ROI subsets (RELAB-referenced)

Post-processing of the ENVI ROI subsets (`outputs/Subset/ch1_subset`, `ch2_subset`),
reproducing Figure 3 of the reference paper (FCC + sampled-point spectra vs a spectral
library). Scripts: `fcc.py`, `mineral_spectra.py`, `mineral_common.py`.

| Step | Script | Output (`outputs/mineral_analysis/`) |
|---|---|---|
| Mafic FCC (R≈2.0 µm, G≈1.0 µm, B=continuum) | `fcc.py` | `ch{1,2}_fcc.{tif,png}` |
| SAM match vs RELAB, 4 points/mineral, plots | `mineral_spectra.py` | `ch{1,2}_<mineral>.png`, `ch{1,2}_points_on_fcc.png`, `extracted_spectra.csv`, `extracted_spectra_ch{1,2}.sli`, `match_report.txt` |
| Paper Fig‑3 style "Spectral Profile" (6 minerals overlaid, mean of 4 points) | `mineral_spectra.py` | `ch{1,2}_mineral_spectral_profile.png` |
| Paper Fig‑3 crater radial transect (6 geomorphic zones) | `crater_profile.py` | `ch{1,2}_crater_profile.png`, `ch{1,2}_crater_points_on_fcc.png`, `crater_profile_points.csv`, `crater_profile_report.txt` |

- **6 minerals** (RELAB library in `outputs/spectral_library/`): olivine, clino/ortho‑pyroxene,
  plagioclase, ilmenite, Mg‑spinel.
- **Matching:** Spectral Angle Mapper on **continuum‑removed** spectra (keys on the diagnostic
  1 µm / 2 µm absorptions, not brightness); 4 best‑matching pixels ≥3 px apart per mineral.
- **Capped at 2600 nm** (RELAB limit): M3 uses 541–2577 nm (bands 1–2 @461/500 nm are 100 % fill
  and dropped); IIRS uses 830–2600 nm (the thermal tail has no library reference).
- Caveat: SAM matches spectral **shape**, not abundance — best‑match points are candidates, not
  confirmed mineralogy; lunar spectra are mixtures. M3 spectra are cleaner than IIRS here.
- **Crater transect** (`crater_profile.py`): the paper's Fig‑3 geomorphic sampling. Crater centre =
  ROI centre; rim radius R = outermost peak of the radial mean‑reflectance profile (M3 R=61 px,
  IIRS R=94 px). Six points along one auto‑chosen azimuth at 0/0.45/0.75/1.0/1.25/1.6 R
  (centre, interior, inner edge, rim, outer edge, outside wall), each a 3×3 mean, overlaid as a
  single "Spectral Profile" graph per sensor.

---

## Opening in ENVI
Every `.img` has a matching `.hdr` with band count, data type, interleave, wavelengths
(and `map info` for georeferenced steps). Open the `.hdr` (or the `.img`) directly in ENVI.
GeoTIFFs open in ENVI/QGIS/ArcGIS directly.

## Provenance / caveats
- CH2 absolute radiance rests on the `DN×0.01` scaling being ISRO's intended L1 scaling
  (the file's PDS4 label reads `processing_level = Raw`). Spectral shapes, band ratios and
  reflectance are valid regardless; only absolute magnitude depends on it.
- Reference: https://github.com/prabhakaralok/CH2IIRS ·
  Verma, Chauhan & Chauhan 2022, *Icarus* 383, 115075.
