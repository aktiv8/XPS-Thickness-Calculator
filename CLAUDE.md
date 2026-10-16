# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Single-file Tkinter app (`app.py`, Python standard library only, no dependencies) that estimates oxide/hydroxide/nitride overlayer thickness from XPS peak areas using the Strohmeier equation (B. R. Strohmeier, *Surf. Interface Anal.* 1990, 15, 51; HarwellXPS Tech Note #1: https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf). `README.md` covers user-facing usage.

## Commands

```bash
python app.py                                                  # launch the Tk GUI (blocks until closed)
python -c "import py_compile; py_compile.compile('app.py', doraise=True)"   # syntax check without opening a window
```

There are no tests, linter or formatter. The calculation core is `strohmeier_thickness()` (a pure function, importable without starting the GUI thanks to the `__main__` guard), fed by `volume_density()`; everything else is Tkinter. To check a change to the maths, call those two directly rather than driving the GUI.

## Things that aren't obvious from a skim

- The single- and mixed-layer equations from the tech note are one code path: overlayer N, λ and RSF are intensity-weighted averages over the active components. With one component and RSF = 1 it reduces exactly to Strohmeier's equation. RSFs cancel for same-core-level metal/oxide pairs, so they default to 1.0.
- N = atoms_per_formula × density / molar_mass (e.g. Cu₂O uses 2 atoms); only ratios matter, so units are arbitrary but must be consistent.
- Take-off angle θ is measured from the **surface** (90° = normal emission), per the tech note. The validity limit reported is ≈ 3× the overlayer IMFP.
- Preset densities, IMFPs and binding energies in `SUBSTRATES` / `OVERLAYERS` are representative defaults, not authoritative; the GUI makes them editable. Binding energies are informational only and don't enter the calculation. The Cu₂O and Cu(OH)₂ Nm/Nox ratios differ ~2% from the tech note because of different reference atomic weights and densities.

## Git

GitHub rejects pushes whose commits carry a private email address (GH007). This repo's local git config uses the noreply identity `171147279+aktiv8@users.noreply.github.com`; keep it that way, and don't change the global git config.
