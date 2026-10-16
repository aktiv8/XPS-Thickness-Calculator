# XPS Overlayer Thickness Calculator

A small Tkinter GUI implementing the **Strohmeier equation** for estimating
oxide/hydroxide/nitride overlayer thickness from a single XPS spectrum, as
described in HarwellXPS Technical Note #1:
https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf
(original method: B. R. Strohmeier, *Surf. Interface Anal.*, 1990, **15**, 51).

## Running it

No third-party packages required (uses only the Python standard library /
Tkinter, which ships with standard Windows Python installs):

```bash
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

See the **Method / References** tab inside the app for the full equations,
derivation notes, and important caveats about the reference density/IMFP
defaults (in particular: IMFP is source/instrument-dependent and film
density can differ from bulk values — verify both for rigorous work).

## Files

- `app.py` — the application (single file).
