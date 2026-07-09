# Roadmap & Status

Status as of 2026‑07‑08.

## Current state — complete and verified

Both channels are fully processed end‑to‑end; all products exist on disk under
`pipeline/outputs/` (~87 GB) as ENVI + GeoTIFF.

### CH2 — Chandrayaan‑2 IIRS  ✅
- [x] **Step 1** Radiometric calibration — L1 integer cube → radiance (`DN × 0.01`, mW/cm²/sr/µm), 256 bands.
- [x] **Aux** Incidence‑angle image from `.spm` (90 − sun_elev, per line). *Verified exact vs ISRO `_angle`.*
- [x] **Step 2** Thermal removal → I/F reflectance, 247 bands (847–4993 nm). *Temp mean |Δ| ≈ 0.1 K, reflectance sub‑% vs `_Temp` / `_inc_corrRef`.*
- [x] **Step 3** Photometric correction — Lommel‑Seeliger disk normalisation to i=30°, e=0°, g=30°.
- [x] **Step 4** Selenoreferencing — 4 PDS4 corner GCPs → Moon 2000 eqc, 94.56 m.
- [x] **Step 5** Destriping (per‑detector moment matching, sensor domain) + selenoreference → **definitive final cube**.

### CH1 — Chandrayaan‑1 M3  ✅
- [x] Selenoreferencing of L2 reflectance (85 bands, 460–3000 nm) via **GCP + thin‑plate‑spline** warp → Moon 2000 eqc, 140 m.

## Key decisions (why the pipeline looks the way it does)

1. **Algorithm fidelity over reinvention.** The CH2 thermal/reflectance math is copied
   verbatim from the ISRO CH2IIRS plugin (cloned in `pipeline/reference/CH2IIRS/`),
   including the 247‑band empirical spectral‑correction coefficients and the exact band
   trim `[7:-2]`. This is what makes cross‑validation against the ISRO files meaningful.

2. **Destripe before warp (step order).** The requested order lists destriping last, but
   striping is a per‑detector‑column sensor artefact and must be removed in *sensor
   geometry* before georeferencing. So step 4 is the georeferenced photometric product and
   step 5 is destripe(sensor)→selenoref, delivered as the final cube. Both are kept.

3. **CH1 warp: GCP+TPS, not geoloc‑array.** The first attempt (`ch1_pipeline.py`,
   `gdalwarp -geoloc` over the full LOC backplane) collapsed longitude. The working
   approach (`ch1_fix_warp.py`) lays a dense GCP grid from LOC and warps with `-tps`.

4. **CH1 polar tail clipped.** This is a near‑polar strip (55°N → ~87°S). Lines with
   lat < −80° are excluded from the equirectangular product (eqc is singular at the pole).
   Data still exists in sensor geometry.

## Known caveats / limitations

- **CH2 absolute radiance rests on an assumption.** `DN × 0.01` is the CH2IIRS L1 scaling,
  but the file's PDS4 label reads `processing_level = Raw`. Spectral shapes, band ratios,
  and reflectance are valid regardless; only absolute magnitude depends on it.
  ISRO's absolute DN→radiance gain coefficients are not published anywhere.
- **Photometric correction is analytic only.** The empirical phase‑curve polynomial of
  Verma et al. 2023 (`PhaseF/*.npy`) is not public; only the Lommel‑Seeliger disk
  normalisation is applied.
- **Selenoreferencing is first‑order.** 4‑corner (CH2) / GCP‑grid (CH1) georeferencing.
  Sub‑pixel accuracy would need a SPICE‑derived per‑pixel geometry backplane.

## Possible next steps (on request — not started)

- [ ] South‑polar stereographic tail for the clipped CH1 lat < −80° segment.
- [ ] SPICE per‑pixel geometry backplane for CH2 to replace 4‑corner georeferencing.
- [ ] Apply the empirical phase‑curve photometric term if ISRO's polynomial becomes available.
- [ ] Mineral / band‑ratio science products (e.g. 1 µm / 2 µm mafic absorptions, OH/H₂O 3 µm).
- [ ] Co‑registration of CH1 and CH2 products onto a common grid for joint analysis.
- [ ] Trim/compress or archive the 87 GB `outputs/` (currently regenerable from scripts).
