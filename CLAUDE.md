# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Tkinter app (`app.py`; the thickness calculator is Python standard library only) that estimates oxide/hydroxide/nitride overlayer thickness from XPS peak areas using the Strohmeier equation (B. R. Strohmeier, *Surf. Interface Anal.* 1990, 15, 51; HarwellXPS Tech Note #1: https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf). `README.md` covers user-facing usage.

## Commands

```bash
python app.py                                                  # launch the Tk GUI (blocks until closed)
python -c "import py_compile; py_compile.compile('app.py', doraise=True)"   # syntax check without opening a window
python -m unittest discover -s tests -v                        # run all tests (no display; numpy tests skip if numpy is missing)
pip install numpy matplotlib                                   # optional: enables the "Simulated spectrum" tab
python -m unittest tests.test_regression.BiesingerAlSpreadsheet   # a single test class
```

There is no linter or formatter. The calculation core is `strohmeier_thickness()` (a pure function, importable without starting the GUI thanks to the `__main__` guard), fed by `volume_density()`; everything else is Tkinter. To check a change to the maths, call those two directly rather than driving the GUI. `tests/test_regression.py` pins the app to published numbers (see below) — if it fails after a preset edit, the preset drifted from its source paper.

## Things that aren't obvious from a skim

- The single- and mixed-layer equations from the tech note are one code path: overlayer N, λ and RSF are intensity-weighted averages over the active components. With one component and RSF = 1 it reduces exactly to Strohmeier's equation. RSFs cancel for same-core-level metal/oxide pairs, so they default to 1.0.
- N = atoms_per_formula × density / molar_mass (e.g. Cu₂O uses 2 atoms); only ratios matter, so units are arbitrary but must be consistent.
- Take-off angle θ is measured from the **surface** (90° = normal emission), per the tech note. The validity limit reported is ≈ 3× the overlayer IMFP.
- Several presets are taken from the original papers, and `tests/test_regression.py` pins them to published numbers (a failing test after a preset edit means the preset drifted from the paper):
  - **Al / Al₂O₃**: Strohmeier 1990 — IMFP 2.8/2.6 nm (Al Kα) and 2.4/2.2 nm (Mg Kα), θ from the **surface**, Nm/Nox = 1.5. The paper rounds the log coefficient to 1.4 (`d = 28 ln(1.4 Io/Im + 1)`); the app uses the exact value (≈1.393), so expect ≤0.4 Å differences. The 3.40 g/cm³ Al₂O₃ density is back-calculated to give Nm/Nox = 1.5 (the paper quotes ~1.3–1.5), not a measured density. Biesinger's Al spreadsheet uses the same inputs and is also reproduced to 0.1 Å.
  - **Cu oxides/hydroxide**: Biesinger 2017 — densities 8.96/6.00/6.31/3.37 g/cm³, IMFPs 1.54/1.77/1.76/2.19 nm. These IMFPs are for **Cu LMM Auger electrons (KE ≈ 918 eV)**, not Cu 2p3/2 (KE ≈ 554 eV), so the presets are keyed `... (LMM)`. The parallel `... (2p3/2)` presets (1.01/1.23/1.22/1.51 nm) are **computed** with `tpp2m_imfp_nm()` (Nv 11/28/17/25, band gaps 2.1/1.25/1.85 eV), not published; the same function reproduces Biesinger's three oxide IMFPs at 918 eV to 0.2%, but his Cu(0) 1.54 nm is ~6% above TPP-2M at 918 eV (it equals TPP-2M at ~1000 eV), so the 2p3/2 metal value may sit slightly low. Never mix the two sets — thickness scales ~linearly with the IMFPs (the 63/13/25 example gives 3.2 nm with LMM vs 2.2 nm with 2p3/2). His ~3.2 nm worked example (63% Cu₂O, 13% Cu(OH)₂, 25% metal) is reproduced only with intensity-weighted N and λ, which is what the app does (averaging the Nm/N ratios instead gives 3.4 nm). The HarwellXPS tech note's Cu numbers (8.93 g/cm³, CuO 6.39) differ by ~1–2%.
  - **W / WO₃**: Carlson & McGuire 1972/73 — 19.3 and 6.8 g/cm³ (WO₃ at 6.8, not the bulk 7.16), IMFP 1.28/2.63 nm, 45° take-off. Their eq. (5) is the same equation as Strohmeier's; the test round-trips it.
  - Strohmeier and Carlson & McGuire were read in full; of Biesinger 2017 only the thickness section and Table 6 were checked. The tech note's own 4.4 nm Al example comes from the tech note's spectrum, not from Strohmeier's paper, so it cannot be reproduced.
- Si, SiO₂, Si₃N₄, Ti, TiO₂: densities and band gaps come from the user's Thermo Avantage `Densities.txt` (Si 2.33/1.12 eV, SiO₂ 2.196/8.9, Si₃N₄ 3.44/4.65, Ti 4.54, TiO₂ 4.26/3.05); IMFPs are **computed** with `tpp2m_imfp_nm()` (Nv 4/16/32/4/16) at Al Kα, not published. The same file lists Al 2.6989 and Al₂O₃ 3.97 g/cm³ (sapphire), which would give Nm/Nox ≈ 1.29, but the Al presets deliberately stay on Strohmeier's 1.5 (Al₂O₃ 3.40 back-calculated). Ti/TiO₂ binding energies are typical literature values, not from the file.
- Band gaps live in `BAND_GAPS_EV` (value, source) and in each preset's `band_gap_ev` (0 = metal, None = unknown); only `tpp2m_imfp_nm()` uses them. **Policy set by the user: never overwrite an existing entry** — the Avantage values win over any other list (e.g. their pasted list had Si 1.11, Ge 0.66, GaAs 1.43, SiO₂ 9.0, all ignored; only InP 1.27, GaP 2.25, CdTe 1.44, PbS 0.37 were new, tagged unsourced). `tests/test_regression.py::BandGapTable` enforces this.
- Extra pairs (Mg/MgO, Ge/GeO₂, Zn/ZnO, Zr/ZrO₂, Hf/HfO₂, Sn/SnO₂, SiC→SiO₂, GaAs & GaP→Ga₂O₃, CdTe→TeO₂, plus InP/PbS below) were added at the user's request. Densities and gaps use the Avantage file where it has them (Mg, Ge, Sn, SiC, GaAs, ZnO, ZrO₂, HfO₂); the rest (Zn, Zr, Hf metals; MgO, GeO₂, SnO₂, Ga₂O₃, TeO₂; GaP, CdTe) are **my handbook values, not from the user's files**, and the oxide choices are judgement calls (GeO₂ 3.64 amorphous, TeO₂ 6.04 paratellurite). Valence-electron counts (`TPP2M_INPUTS`) count s/p (+ unfilled d) electrons with filled d shells as core (Cu stays 11 to match Biesinger). Sanity check: Nm/Nox equals the Pilling–Bedworth ratio for 1:1 oxides and matches published values.
- InP (In 3d5/2 → In₂O₃) and PbS (S 2p → PbSO₄) were added once the user supplied their densities (4.81, 7.60 g/cm³). In₂O₃ and PbSO₄ densities are handbook values and their band gaps (2.9, 4.5 eV) are rough estimates (`ROUGH` in `BAND_GAPS_EV`). The In 3d5/2 shift (~0.3–0.5 eV) is barely resolvable and InPO₄ (the other native-oxide component) is not modelled; PbS→PbSO₄ is clean because the S 2p sulfate shift is ~7.5 eV (Pb 4f cannot separate PbO well, so PbO is not modelled).
- Deliberately NOT added: ZnS (no clean oxide route), diamond/graphite (no meaningful oxide overlayer), As₂O₃ for GaAs (no band gap/density from the user), CdO for CdTe (Cd 3d can't resolve it).
- SiON is a rough placeholder, not sourced; the GUI makes all presets editable. Binding energies don't enter the thickness calculation, but they position the peaks in the simulated spectrum.

## Simulated spectrum tab (`spectrum_sim.py`)

Optional (numpy + matplotlib; `app.py` guards the import and sets `HAVE_PLOT`, and without them the tab just shows an install hint). `spectrum_sim.py` is pure numpy with no Tk, so it is unit-tested directly (`tests/test_spectrum_sim.py`); `app.py` only holds the Tk/matplotlib glue (`_build_spectrum_tab`, `_refresh_spectrum`, `_draw_spectrum`, fed by `self.last_calc` which `_on_calculate` fills).

- It is an **estimate, not a fit**, and every figure says so (`DISCLAIMER`). `oxide_area_fractions()` is the forward form of the Strohmeier equation and reuses the `N_ov_avg` / `lam_ov_avg` / `rsf_ov_avg` that `strohmeier_thickness()` returns, so simulating at the calculated d reproduces the entered area ratio exactly (tested). Intensities are relative (areas normalised to 100%).
- `CORE_LEVEL_MODEL` (spin-orbit splitting, branching ratio, default FWHM, Lorentzian mix) holds **typical literature values from general knowledge, not from the user's files**; FWHM is editable in the tab. The preset BE is taken as the lower-BE j-component (2p3/2, 3d5/2, 4f7/2) with the partner at +splitting. `model_for()` returns None for Auger presets (`Cu LMM`, `Zn LMM`), which are deliberately not simulated; a test enforces that every preset with a BE has a model and every Auger preset has none.
- Not modelled: shake-up satellites (CuO, Cu(OH)₂), asymmetric metal tails, plasmon/loss features, Auger lines. Background is a self-consistent Shirley step (15% of peak height). The default window clips a little Lorentzian tail (≤ ~3% of area), which the tests allow for.
- A user spectrum is loaded as 2-column text/CSV (BE, intensity; headers skipped) and min–max scaled to the simulation height, so it is a shape comparison, not a quantitative overlay.

## Git

GitHub rejects pushes whose commits carry a private email address (GH007). This repo's local git config uses the noreply identity `171147279+aktiv8@users.noreply.github.com`; keep it that way, and don't change the global git config.
