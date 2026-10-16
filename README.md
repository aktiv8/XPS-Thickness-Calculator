# XPS Overlayer Thickness Calculator

A small Tkinter GUI implementing the **Strohmeier equation** for estimating
oxide/hydroxide/nitride overlayer thickness from a single XPS spectrum, as
described in HarwellXPS Technical Note #1:
https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf
(original method: B. R. Strohmeier, *Surf. Interface Anal.*, 1990, **15**, 51).

## Running it

The thickness calculator needs no third-party packages (only the Python
standard library / Tkinter, which ships with standard Windows Python installs).
The optional **Simulated spectrum** pane needs numpy and matplotlib:

```bash
pip install numpy matplotlib   # optional
python app.py
```

## What it does

- Pick a **substrate** (Cu, Al, Si, or custom) and up to 4 **overlayer**
  components (CuO, Cu2O, Cu(OH)2, Al2O3, SiO2, Si3N4, SiON, or custom) —
  tick the ones present and enter their fitted peak areas.
- Density, molar mass, atoms/formula unit, typical binding energy and IMFP
  are pre-filled from a small reference table and are all editable.
- Enter the photoelectron take-off angle and the substrate ("metal") peak
  area, then click **Calculate thickness**.
- Handles both single-component (metal vs. one oxide) and mixed-layer
  (e.g. oxide + hydroxide present together) cases automatically, using
  intensity-weighted average density/IMFP for the mixed case, matching the
  two forms of the equation given in the tech note.
- An optional RSF (relative sensitivity factor) field on each row
  generalises the equation for cases mixing different core levels/elements;
  leave RSF = 1.0 for the standard same-element case, where it cancels.
- Flags whether the result is within the method's validity limit
  (~3x the overlayer IMFP).
- Results can be exported to CSV via **Export history to CSV...**.
- **One window:** inputs on the left, the calculation (big blue
  **Calculate thickness** button, a large thickness readout and a green/amber/red
  validity badge) on the right, and the simulated spectrum underneath. Drag the
  divider between them to trade height. On short screens (laptops) it switches
  to a compact layout; pressing Enter in any input box also calculates.
- **Simulated spectrum (estimate)** pane: after a calculation it draws what that
  thickness would look like for the substrate and oxide core-level peaks
  (spin-orbit doublets, Shirley-type background) next to a to-scale layer
  diagram, clearly labelled as an estimate. Use **Load my spectrum...**
  (2-column text/CSV of binding energy, intensity) to overlay your own data,
  the thickness box/slider for what-if values, the BE shift for charge
  correction and the FWHM boxes to match your resolution. The figure can be
  saved as PNG from its toolbar. Peak splitting/FWHM are typical defaults.
- **What the simulation leaves out is stated on the figure:** for simplicity it
  does not draw shake-up satellites, plasmon/loss features, Auger lines or
  asymmetric metal tails. For Cu 2p3/2 with CuO or Cu(OH)₂ the Cu(II) shake-up
  satellite region is shaded and labelled "not simulated" (the satellites carry
  roughly 35–39 % of the Cu(II) area, so the drawn Cu(II) main peak is
  overstated). Auger presets (e.g. Cu LMM) are not simulated.

See the **Method / References** tab inside the app for the full equations,
derivation notes, and important caveats about the reference density/IMFP
defaults (in particular: IMFP is source/instrument-dependent and film
density can differ from bulk values — verify both for rigorous work).

Several presets come straight from the original papers, and the tests pin them
to the published numbers:

| System | Source | Inputs |
|---|---|---|
| Al / Al₂O₃ (Al Kα and Mg Kα) | Strohmeier, *Surf. Interface Anal.* 1990, **15**, 51 | IMFP 2.8 / 2.6 nm (Al Kα), 2.4 / 2.2 nm (Mg Kα), Nm/Nox = 1.5 |
| Cu / Cu₂O / CuO / Cu(OH)₂ | Biesinger, *Surf. Interface Anal.* 2017, **49**, 1325 | 8.96 / 6.00 / 6.31 / 3.37 g/cm³; IMFP 1.54 / 1.77 / 1.76 / 2.19 nm |
| W / WO₃ | Carlson & McGuire, *J. Electron Spectrosc.* 1972/73, **1**, 161 | 19.3 / 6.8 g/cm³; IMFP 1.28 / 2.63 nm; 45° take-off |

The Al preset also reproduces M. C. Biesinger's Al oxide/hydroxide spreadsheet to
0.1 Å. The Cu IMFPs are for Cu LMM Auger electrons (KE ≈ 918 eV), as in the
paper, so use the `(LMM)` presets with Cu LMM areas. Separate `(2p3/2)` presets
(IMFP 1.01 / 1.23 / 1.22 / 1.51 nm for Cu / Cu₂O / CuO / Cu(OH)₂ at KE ≈ 554 eV) are
**computed** with TPP-2M (`tpp2m_imfp_nm()` in `app.py`), not published values; the same
code reproduces Biesinger's oxide IMFPs at 918 eV to 0.2 %. Don't mix the two sets —
thickness scales almost linearly with the IMFPs.

Si, SiO₂, Si₃N₄, Ti and TiO₂ take density and band gap from the Thermo Avantage
`Densities.txt` (2.33/2.196/3.44/4.54/4.26 g/cm³), with IMFPs **computed** by TPP-2M
(3.09/3.89/3.00/2.19/2.17 nm at Al Kα). Note the 3.44 g/cm³ Si₃N₄ density is high for
thin amorphous films. SiON is still a rough placeholder.

More pairs: Mg/MgO, Ge/GeO₂, Zn/ZnO (Zn LMM line), Zr/ZrO₂, Hf/HfO₂, Sn/SnO₂, SiC (→ SiO₂), GaAs and GaP
(Ga 3d → Ga₂O₃), CdTe (Te 3d5/2 → TeO₂), InP (In 3d5/2 → In₂O₃; small shift) and PbS (S 2p → PbSO₄). Densities/band gaps come from the Avantage file where it has
them, otherwise standard handbook values (not from the user's files); IMFPs are computed with TPP-2M.
Their Nm/Nox ratios agree with the Pilling–Bedworth ratios (Mg 0.79, Zr 1.55, Hf 1.62, Sn 1.34, Zn 1.58).
Not included: ZnS, diamond and graphite (see CLAUDE.md for why).

`BAND_GAPS_EV` in `app.py` is a reference table of band gaps (used only for TPP-2M IMFP
estimates): the Avantage file's values, Biesinger's Cu oxide gaps, and four user-supplied
additions (InP 1.27, GaP 2.25, CdTe 1.44, PbS 0.37 eV). Existing entries are never overwritten.

## Tests

```bash
python -m unittest discover -s tests -v
```

No display needed; the numpy-based spectrum tests skip themselves if numpy is missing. `tests/test_regression.py` checks published numbers from
the three papers above and Biesinger's Al spreadsheet (65/35 → 35.8 Å,
23/77 → 9.7 Å, 85/15 → 61.2 Å), plus a few properties of the equation.

## Files

- `app.py` — the application (calculator, presets, GUI).
- `spectrum_sim.py` — simulated-spectrum maths (numpy only, no GUI).
- `tests/test_regression.py` — regression tests against published numbers.
- `tests/test_spectrum_sim.py` — tests for the spectrum simulation.
