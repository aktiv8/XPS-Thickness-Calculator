"""
XPS Overlayer Thickness Calculator
===================================

Implements the Strohmeier uniform-overlayer equation for estimating the
thickness of a thin oxide / hydroxide / nitride layer on a metal (or other
substrate) from XPS peak areas, as described in:

    B. R. Strohmeier, Surf. Interface Anal., 1990, 15, 51
    HarwellXPS Technical Note #1 - "Oxide Thickness Determination
    Using the Strohmeier Equation"
    https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf

Core (single overlayer) equation:

    d = lambda_ox * sin(theta) * ln[ (N_m * lambda_m * I_ox) / (N_ox * lambda_ox * I_m) + 1 ]

Mixed overlayer (e.g. oxide + hydroxide) equation:

    d = lambda_ave * sin(theta) * ln[ (N_m * lambda_m * (I_oxide + I_hydroxide))
                                        / (N_ave * lambda_ave * I_m) + 1 ]

where lambda_ave / N_ave are the intensity-weighted averages of the
individual overlayer components' IMFPs and volume densities.

This tool generalises both equations with an optional RSF (relative
sensitivity factor) term for the substrate and each overlayer component.
For the classic Strohmeier case (metal vs. its own oxide, same core level)
all RSFs cancel and should be left at 1.0. RSFs matter if you intend to mix
intensities from different core levels / elements (e.g. an overlayer signal
read from N 1s while the substrate signal is read from Si 2p) - the general
form used here is:

    d = lambda_ov * sin(theta) * ln[ (RSF_sub * N_sub * lambda_sub * I_ov)
                                       / (RSF_ov * N_ov * lambda_ov * I_sub) + 1 ]

Only the intensity (peak area) values, take-off angle, and RSFs are
measurement-specific. Densities, molar masses and IMFPs are pre-filled with
representative literature/TPP-2M-type values for common systems, but ARE
EDITABLE - you should verify/replace them (e.g. via NIST SRD 82, or your own
reference measurements) for your specific instrument, X-ray source and
material phase (amorphous vs. crystalline density strongly affects results).

The thickness calculator needs no third-party packages (standard library + tkinter only). The optional
"Simulated spectrum" pane (see spectrum_sim.py) needs numpy and matplotlib; without them it shows an install hint.
"""

import csv
import math
import time
import tkinter as tk
from dataclasses import dataclass, field
from typing import Optional
from tkinter import ttk, messagebox, filedialog
import os
import textwrap

# Optional: the simulated-spectrum pane needs numpy + matplotlib. The thickness calculator itself does not.
try:
    import numpy as np  # noqa: F401  (used by spectrum_sim and the plotting code below)
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    import spectrum_sim as sim
    HAVE_PLOT = True
except Exception:
    HAVE_PLOT = False


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

@dataclass
class Compound:
    label: str
    density_g_cm3: float
    molar_mass_g_mol: float
    atoms_per_formula: float   # number of "substrate element" atoms per formula unit
    core_level: str
    be_ev: float                # typical binding energy of this component (eV)
    imfp_nm: float               # typical inelastic mean free path (nm)
    rsf: float = 1.0
    band_gap_ev: Optional[float] = None   # eV; 0 = metal, None = unknown. Only TPP-2M uses it.
    note: str = ""


# Band gaps (eV), name -> (value, source). Only the TPP-2M IMFP estimate uses these, not the thickness equation.
# Existing entries are never overwritten by later lists: where a material appears in more than one source the
# first source here wins (Avantage, then Biesinger), and the user-supplied list only adds materials not already present.
AVANTAGE = "Thermo Avantage Densities.txt"
BIESINGER = "Biesinger, Surf. Interface Anal. 2017, 49, 1325"
USER_LIST = "user-supplied list (unsourced)"
HANDBOOK = "standard literature value, not from the user's files"
ROUGH = "rough estimate, not from the user's files"
BAND_GAPS_EV = {
    "SiO2": (8.9, AVANTAGE), "Sn (gray)": (0.08, AVANTAGE), "HfO2": (5.7, AVANTAGE), "ZrO2": (6.0, AVANTAGE),
    "Si": (1.12, AVANTAGE), "TiO2": (3.05, AVANTAGE), "Si3N4": (4.65, AVANTAGE), "Al2O3": (8.7, AVANTAGE),
    "GaAs": (1.35, AVANTAGE), "SiC": (2.86, AVANTAGE), "ZnO": (3.2, AVANTAGE), "ZnS": (3.6, AVANTAGE),
    "Diamond": (5.5, AVANTAGE), "Graphite": (0.0, AVANTAGE), "Ge": (0.67, AVANTAGE),
    "Cu2O": (2.1, BIESINGER), "CuO": (1.25, BIESINGER), "Cu(OH)2": (1.85, BIESINGER),
    "MgO": (7.8, HANDBOOK), "GeO2": (5.6, HANDBOOK), "SnO2": (3.6, HANDBOOK), "Ga2O3": (4.8, HANDBOOK),
    "TeO2": (3.7, HANDBOOK), "In2O3": (2.9, ROUGH), "PbSO4": (4.5, ROUGH),
    "InP": (1.27, USER_LIST), "GaP": (2.25, USER_LIST), "CdTe": (1.44, USER_LIST), "PbS": (0.37, USER_LIST),
}

# Substrate ("metal") reference compounds
SUBSTRATES = {
    "Cu (metal, LMM)": Compound(
        "Cu (metal, LMM)", density_g_cm3=8.96, molar_mass_g_mol=63.546,
        atoms_per_formula=1, core_level="Cu LMM", be_ev=0.0, imfp_nm=1.54,
        note="Biesinger, Surf. Interface Anal. 2017, 49, 1325. IMFP 1.54 nm for Cu L3M4,5M4,5 Auger electrons, KE ~918 eV "
             "(Cu LMM KE 918.61 eV; Cu 2p3/2 BE 932.63 eV). Use with areas from the Cu LMM fit; for Cu 2p3/2 areas "
             "pick the '2p3/2' presets instead (do not mix the two sets).",
    ),
    "Cu (metal, 2p3/2)": Compound(
        "Cu (metal, 2p3/2)", density_g_cm3=8.96, molar_mass_g_mol=63.546,
        atoms_per_formula=1, core_level="Cu 2p3/2", be_ev=932.63, imfp_nm=1.01,
        note="IMFP 1.01 nm COMPUTED with TPP-2M (Nv = 11) at KE 554 eV (Al Ka - Cu 2p3/2); not a published value. "
             "For comparison, Biesinger's 1.54 nm (LMM set) is ~6% above TPP-2M at 918 eV, so this metal value may be "
             "slightly low relative to a NIST SRD 71 element calculation - edit if you have one.",
    ),
    "W (metal)": Compound(
        "W (metal)", density_g_cm3=19.3, molar_mass_g_mol=183.84,
        atoms_per_formula=1, core_level="W 4f7/2", be_ev=31.4, imfp_nm=1.28,
        note="Carlson & McGuire, J. Electron Spectrosc. 1972/73, 1, 161: IMFP 12.8 A at KE ~1455 eV (Al Ka). Their "
             "geometry was 45 deg take-off. BE is a typical literature value, not from the paper.",
    ),
    "Al (metal)": Compound(
        "Al (metal)", density_g_cm3=2.70, molar_mass_g_mol=26.98,
        atoms_per_formula=1, core_level="Al 2p", be_ev=72.9, imfp_nm=2.6,
        note="Al Ka source, KE ~1413 eV. IMFP 2.6 nm = Strohmeier, Surf. Interface Anal. 1990, 15, 51 (from Tanuma, "
             "Powell & Penn 1988); also used in Biesinger's Al calculator.",
    ),
    "Al (metal, Mg Ka)": Compound(
        "Al (metal, Mg Ka)", density_g_cm3=2.70, molar_mass_g_mol=26.98,
        atoms_per_formula=1, core_level="Al 2p", be_ev=72.9, imfp_nm=2.2,
        note="Mg Ka source, KE ~1180 eV. IMFP 2.2 nm = Strohmeier (1990). Pair with the 'Al2O3 (Mg Ka)' overlayer.",
    ),
    "Si (metal/elemental)": Compound(
        "Si (elemental)", density_g_cm3=2.33, molar_mass_g_mol=28.09,
        atoms_per_formula=1, core_level="Si 2p3/2", be_ev=99.3, imfp_nm=3.09,
        note="Density and band gap (1.12 eV) from Thermo Avantage Densities.txt. IMFP 3.09 nm COMPUTED with TPP-2M "
             "(Nv = 4) at KE ~1387 eV (Al Ka); not a published value.",
    ),
    "Ti (metal)": Compound(
        "Ti (metal)", density_g_cm3=4.54, molar_mass_g_mol=47.867,
        atoms_per_formula=1, core_level="Ti 2p3/2", be_ev=453.9, imfp_nm=2.19,
        note="Density from Thermo Avantage Densities.txt. IMFP 2.19 nm COMPUTED with TPP-2M (Nv = 4) at KE ~1033 eV "
             "(Al Ka); not a published value. BE is a typical literature value, not from the file.",
    ),
    "Mg (metal)": Compound(
        "Mg (metal)", density_g_cm3=1.7, molar_mass_g_mol=24.305,
        atoms_per_formula=1, core_level="Mg 2p", be_ev=49.7, imfp_nm=3.85,
        note="Density 1.7 g/cm3 (Thermo Avantage Densities.txt); metal (no band gap). IMFP 3.85 nm COMPUTED with TPP-2M "
             "(Nv = 2) at KE ~1437 eV (Al Ka); not a published value. BE is a typical literature value, not from your "
             "files. Mg 2p oxide shift is only ~1 eV; the Mg KLL Auger line separates metal and oxide better.",
    ),
    "Ge (elemental)": Compound(
        "Ge (elemental)", density_g_cm3=5.323, molar_mass_g_mol=72.63,
        atoms_per_formula=1, core_level="Ge 3d", be_ev=29.4, imfp_nm=2.99,
        note="Density 5.323 g/cm3 (Thermo Avantage Densities.txt); band gap 0.67 eV (Thermo Avantage Densities.txt). IMFP "
             "2.99 nm COMPUTED with TPP-2M (Nv = 4) at KE ~1457 eV (Al Ka); not a published value. BE is a typical "
             "literature value, not from your files.",
    ),
    "Zn (metal, LMM)": Compound(
        "Zn (metal, LMM)", density_g_cm3=7.14, molar_mass_g_mol=65.38,
        atoms_per_formula=1, core_level="Zn LMM", be_ev=0, imfp_nm=2.39,
        note="Density 7.14 g/cm3 (standard handbook value, not from your files); metal (no band gap). IMFP 2.39 nm "
             "COMPUTED with TPP-2M (Nv = 2) at KE ~992 eV (Al Ka); not a published value. Zn 2p3/2 cannot resolve metal "
             "from ZnO (<0.3 eV), so use the Zn L3M4,5M4,5 line (KE ~992 eV metal, ~988 eV oxide); BE left blank.",
    ),
    "Zr (metal)": Compound(
        "Zr (metal)", density_g_cm3=6.52, molar_mass_g_mol=91.224,
        atoms_per_formula=1, core_level="Zr 3d5/2", be_ev=178.8, imfp_nm=2.72,
        note="Density 6.52 g/cm3 (standard handbook value, not from your files); metal (no band gap). IMFP 2.72 nm "
             "COMPUTED with TPP-2M (Nv = 4) at KE ~1308 eV (Al Ka); not a published value. BE is a typical literature "
             "value, not from your files.",
    ),
    "Hf (metal)": Compound(
        "Hf (metal)", density_g_cm3=13.31, molar_mass_g_mol=178.49,
        atoms_per_formula=1, core_level="Hf 4f7/2", be_ev=14.2, imfp_nm=2.79,
        note="Density 13.31 g/cm3 (standard handbook value, not from your files); metal (no band gap). IMFP 2.79 nm "
             "COMPUTED with TPP-2M (Nv = 4) at KE ~1472 eV (Al Ka); not a published value. BE is a typical literature "
             "value, not from your files.",
    ),
    "Sn (metal)": Compound(
        "Sn (metal)", density_g_cm3=7.31, molar_mass_g_mol=118.71,
        atoms_per_formula=1, core_level="Sn 3d5/2", be_ev=484.9, imfp_nm=2.31,
        note="Density 7.31 g/cm3 (Thermo Avantage Densities.txt); metal (no band gap). IMFP 2.31 nm COMPUTED with TPP-2M "
             "(Nv = 4) at KE ~1002 eV (Al Ka); not a published value. BE is a typical literature value, not from your "
             "files. Density is the Avantage 'Sn white' entry.",
    ),
    "SiC": Compound(
        "SiC", density_g_cm3=3.21, molar_mass_g_mol=40.097,
        atoms_per_formula=1, core_level="Si 2p", be_ev=100.4, imfp_nm=2.90,
        note="Density 3.21 g/cm3 (Thermo Avantage Densities.txt); band gap 2.86 eV (Thermo Avantage Densities.txt). IMFP "
             "2.90 nm COMPUTED with TPP-2M (Nv = 8) at KE ~1386 eV (Al Ka); not a published value. BE is a typical "
             "literature value, not from your files. 1 Si atom per formula unit; oxide is the SiO2 overlayer preset (Si 2p "
             "at ~103.4 eV).",
    ),
    "GaAs (Ga 3d)": Compound(
        "GaAs (Ga 3d)", density_g_cm3=5.31, molar_mass_g_mol=144.645,
        atoms_per_formula=1, core_level="Ga 3d", be_ev=19.3, imfp_nm=3.02,
        note="Density 5.31 g/cm3 (Thermo Avantage Densities.txt); band gap 1.35 eV (Thermo Avantage Densities.txt). IMFP "
             "3.02 nm COMPUTED with TPP-2M (Nv = 8) at KE ~1467 eV (Al Ka); not a published value. BE is a typical "
             "literature value, not from your files. 1 Ga atom per formula unit; pair with Ga2O3. Native oxide also "
             "contains As oxides, which this preset does not model.",
    ),
    "GaP (Ga 3d)": Compound(
        "GaP (Ga 3d)", density_g_cm3=4.14, molar_mass_g_mol=100.697,
        atoms_per_formula=1, core_level="Ga 3d", be_ev=19.2, imfp_nm=3.04,
        note="Density 4.14 g/cm3 (standard handbook value, not from your files); band gap 2.25 eV (user-supplied list "
             "(unsourced)). IMFP 3.04 nm COMPUTED with TPP-2M (Nv = 8) at KE ~1467 eV (Al Ka); not a published value. BE "
             "is a typical literature value, not from your files. 1 Ga atom per formula unit; pair with Ga2O3. Native "
             "oxide is often partly GaPO4, which this preset does not model.",
    ),
    "CdTe (Te 3d5/2)": Compound(
        "CdTe (Te 3d5/2)", density_g_cm3=5.85, molar_mass_g_mol=240.014,
        atoms_per_formula=1, core_level="Te 3d5/2", be_ev=572.1, imfp_nm=2.37,
        note="Density 5.85 g/cm3 (standard handbook value, not from your files); band gap 1.44 eV (user-supplied list "
             "(unsourced)). IMFP 2.37 nm COMPUTED with TPP-2M (Nv = 8) at KE ~914 eV (Al Ka); not a published value. BE is "
             "a typical literature value, not from your files. 1 Te atom per formula unit; pair with TeO2 (Te 3d5/2 shifts "
             "~4 eV). CdO cannot be resolved from CdTe by Cd 3d and is not modelled.",
    ),
    "InP (In 3d5/2)": Compound(
        "InP (In 3d5/2)", density_g_cm3=4.81, molar_mass_g_mol=145.792,
        atoms_per_formula=1, core_level="In 3d5/2", be_ev=444.6, imfp_nm=2.40,
        note="Density 4.81 g/cm3 (your message); band gap 1.27 eV (user-supplied list (unsourced)). IMFP 2.40 nm COMPUTED "
             "with TPP-2M (Nv = 8) at KE ~1042 eV (Al Ka); not a published value. BE is a typical literature value, not "
             "from your files. 1 In atom per formula unit; pair with In2O3. The In 3d5/2 oxide shift is only ~0.3-0.5 eV, "
             "so the fit is hard; the In MNN Auger line separates them better. Native oxide also contains InPO4, which "
             "this preset does not model.",
    ),
    "PbS (S 2p)": Compound(
        "PbS (S 2p)", density_g_cm3=7.6, molar_mass_g_mol=239.26,
        atoms_per_formula=1, core_level="S 2p3/2", be_ev=161, imfp_nm=2.63,
        note="Density 7.6 g/cm3 (your message); band gap 0.37 eV (user-supplied list (unsourced)). IMFP 2.63 nm COMPUTED "
             "with TPP-2M (Nv = 10) at KE ~1326 eV (Al Ka); not a published value. BE is a typical literature value, not "
             "from your files. 1 S atom per formula unit; pair with PbSO4 (S 2p3/2 sulfate shift ~7.5 eV, easily "
             "resolved). PbO is not modelled (Pb 4f cannot separate it well).",
    ),
    "Custom / manual entry": Compound(
        "Custom", density_g_cm3=1.0, molar_mass_g_mol=1.0,
        atoms_per_formula=1, core_level="", be_ev=0.0, imfp_nm=1.0,
        note="Enter your own values.",
    ),
}

# Overlayer reference compounds (oxides / hydroxides / nitrides)
OVERLAYERS = {
    "None": None,
    "CuO (LMM)": Compound(
        "CuO (LMM)", density_g_cm3=6.31, molar_mass_g_mol=79.55,
        atoms_per_formula=1, core_level="Cu LMM", be_ev=0.0, imfp_nm=1.76,
        note="Biesinger 2017: 6.31 g/cm3, IMFP 1.76 nm at Cu LMM KE ~918 eV (Cu LMM KE 917.57 eV; Cu 2p3/2 BE 933.76 eV); "
             "Nm/Nox = 1.78.",
    ),
    "Cu2O (LMM)": Compound(
        "Cu2O (LMM)", density_g_cm3=6.00, molar_mass_g_mol=143.09,
        atoms_per_formula=2, core_level="Cu LMM", be_ev=0.0, imfp_nm=1.77,
        note="Biesinger 2017: 6.00 g/cm3, IMFP 1.77 nm at Cu LMM KE ~918 eV (Cu LMM KE 916.99 eV; Cu 2p3/2 BE 932.18 eV); "
             "Nm/Nox = 1.68 (2 Cu atoms per formula unit). Cu 2p3/2 BE overlaps Cu0, so use the Cu LMM line to separate them.",
    ),
    "Cu(OH)2 (LMM)": Compound(
        "Cu(OH)2 (LMM)", density_g_cm3=3.37, molar_mass_g_mol=97.56,
        atoms_per_formula=1, core_level="Cu LMM", be_ev=0.0, imfp_nm=2.19,
        note="Biesinger 2017: 3.37 g/cm3, IMFP 2.19 nm at Cu LMM KE ~918 eV (Cu LMM KE 916.25 eV; Cu 2p3/2 BE 934.67 eV); "
             "Nm/Nox = 4.08.",
    ),
    "CuO (2p3/2)": Compound(
        "CuO (2p3/2)", density_g_cm3=6.31, molar_mass_g_mol=79.55,
        atoms_per_formula=1, core_level="Cu 2p3/2", be_ev=933.76, imfp_nm=1.22,
        note="IMFP 1.22 nm COMPUTED with TPP-2M (Nv = 17, band gap 1.25 eV) at KE 554 eV (Al Ka); not a published value. "
             "Same method and inputs reproduce Biesinger's 918 eV value to 0.1%.",
    ),
    "Cu2O (2p3/2)": Compound(
        "Cu2O (2p3/2)", density_g_cm3=6.00, molar_mass_g_mol=143.09,
        atoms_per_formula=2, core_level="Cu 2p3/2", be_ev=932.18, imfp_nm=1.23,
        note="IMFP 1.23 nm COMPUTED with TPP-2M (Nv = 28, band gap 2.1 eV) at KE 554 eV (Al Ka); not a published value. "
             "Cu 2p3/2 BE overlaps Cu0, so a Cu2O/Cu0 split from 2p3/2 alone is unreliable - prefer the Cu LMM line.",
    ),
    "Cu(OH)2 (2p3/2)": Compound(
        "Cu(OH)2 (2p3/2)", density_g_cm3=3.37, molar_mass_g_mol=97.56,
        atoms_per_formula=1, core_level="Cu 2p3/2", be_ev=934.67, imfp_nm=1.51,
        note="IMFP 1.51 nm COMPUTED with TPP-2M (Nv = 25, band gap 1.85 eV) at KE 554 eV (Al Ka); not a published value.",
    ),
    "WO3": Compound(
        "WO3", density_g_cm3=6.8, molar_mass_g_mol=231.84,
        atoms_per_formula=1, core_level="W 4f7/2", be_ev=35.8, imfp_nm=2.63,
        note="Carlson & McGuire 1972/73: anodic WO3 taken as 6.8 g/cm3 (bulk 7.16), IMFP 26.3 A, Nm/Nox = 3.58 with the "
             "W preset; use theta = 45 deg to match their geometry. BE is a typical literature value, not from the paper.",
    ),
    "Al2O3": Compound(
        "Al2O3", density_g_cm3=3.40, molar_mass_g_mol=101.96,
        atoms_per_formula=2, core_level="Al 2p", be_ev=74.4, imfp_nm=2.8,
        note="Strohmeier (1990): IMFP 2.8 nm (Al Ka) and Nm/Nox = 1.5 (he quotes ~1.3-1.5 from gamma-Al2O3 densities and "
             "uses 1.5). The 3.40 g/cm3 density is back-calculated to give that ratio, not a measured density "
             "(native/anodic films ~3.0-3.2, sapphire 3.95). Reproduces Biesinger's Al calculator - edit if your film differs.",
    ),
    "Al2O3 (Mg Ka)": Compound(
        "Al2O3 (Mg Ka)", density_g_cm3=3.40, molar_mass_g_mol=101.96,
        atoms_per_formula=2, core_level="Al 2p", be_ev=74.4, imfp_nm=2.4,
        note="Mg Ka source (KE ~1180 eV): IMFP 2.4 nm (Strohmeier 1990), Nm/Nox = 1.5. Pair with 'Al (metal, Mg Ka)'.",
    ),
    "SiO2": Compound(
        "SiO2", density_g_cm3=2.196, molar_mass_g_mol=60.08,
        atoms_per_formula=1, core_level="Si 2p3/2", be_ev=103.4, imfp_nm=3.89,
        note="Density 2.196 g/cm3 and band gap 8.9 eV from Thermo Avantage Densities.txt (amorphous; quartz is 2.65). "
             "IMFP 3.89 nm COMPUTED with TPP-2M (Nv = 16) at KE ~1387 eV (Al Ka); not a published value.",
    ),
    "Si3N4": Compound(
        "Si3N4", density_g_cm3=3.44, molar_mass_g_mol=140.28,
        atoms_per_formula=3, core_level="Si 2p3/2", be_ev=101.8, imfp_nm=3.00,
        note="Density 3.44 g/cm3 and band gap 4.65 eV from Thermo Avantage Densities.txt; thin amorphous/LPCVD films "
             "are often lower (~2.9-3.1), which raises Nm/Nox. IMFP 3.00 nm COMPUTED with TPP-2M (Nv = 32) at KE ~1387 eV "
             "(Al Ka); not a published value.",
    ),
    "TiO2": Compound(
        "TiO2", density_g_cm3=4.26, molar_mass_g_mol=79.866,
        atoms_per_formula=1, core_level="Ti 2p3/2", be_ev=458.8, imfp_nm=2.17,
        note="Density 4.26 g/cm3 and band gap 3.05 eV from Thermo Avantage Densities.txt. IMFP 2.17 nm COMPUTED with "
             "TPP-2M (Nv = 16) at KE ~1033 eV (Al Ka); not a published value. BE is a typical literature value, not from the file.",
    ),
    "MgO": Compound(
        "MgO", density_g_cm3=3.58, molar_mass_g_mol=40.304,
        atoms_per_formula=1, core_level="Mg 2p", be_ev=50.8, imfp_nm=3.22,
        note="Density 3.58 g/cm3 (standard handbook value, not from your files); band gap 7.8 eV (standard handbook value, "
             "not from your files). IMFP 3.22 nm COMPUTED with TPP-2M (Nv = 8) at KE ~1436 eV (Al Ka); not a published "
             "value. BE is a typical literature value, not from your files. Mg 2p oxide shift is only ~1 eV; the Mg KLL "
             "Auger line separates metal and oxide better.",
    ),
    "GeO2": Compound(
        "GeO2", density_g_cm3=3.64, molar_mass_g_mol=104.628,
        atoms_per_formula=1, core_level="Ge 3d", be_ev=32.9, imfp_nm=3.03,
        note="Density 3.64 g/cm3 (standard handbook value, not from your files); band gap 5.6 eV (standard handbook value, "
             "not from your files). IMFP 3.03 nm COMPUTED with TPP-2M (Nv = 16) at KE ~1454 eV (Al Ka); not a published "
             "value. BE is a typical literature value, not from your files. 3.64 g/cm3 is amorphous GeO2 (hexagonal 4.23, "
             "rutile 6.24) - edit for your film.",
    ),
    "ZnO (LMM)": Compound(
        "ZnO (LMM)", density_g_cm3=5.61, molar_mass_g_mol=81.379,
        atoms_per_formula=1, core_level="Zn LMM", be_ev=0, imfp_nm=1.95,
        note="Density 5.61 g/cm3 (Thermo Avantage Densities.txt); band gap 3.2 eV (Thermo Avantage Densities.txt). IMFP "
             "1.95 nm COMPUTED with TPP-2M (Nv = 8) at KE ~988 eV (Al Ka); not a published value. Use with areas from the "
             "Zn L3M4,5M4,5 line (KE ~988 eV oxide); BE left blank.",
    ),
    "ZrO2": Compound(
        "ZrO2", density_g_cm3=5.68, molar_mass_g_mol=123.222,
        atoms_per_formula=1, core_level="Zr 3d5/2", be_ev=182.3, imfp_nm=2.42,
        note="Density 5.68 g/cm3 (Thermo Avantage Densities.txt); band gap 6 eV (Thermo Avantage Densities.txt). IMFP 2.42 "
             "nm COMPUTED with TPP-2M (Nv = 16) at KE ~1304 eV (Al Ka); not a published value. BE is a typical literature "
             "value, not from your files. Band gap is the Avantage average (5-7 eV depending on phase).",
    ),
    "HfO2": Compound(
        "HfO2", density_g_cm3=9.68, molar_mass_g_mol=210.488,
        atoms_per_formula=1, core_level="Hf 4f7/2", be_ev=16.8, imfp_nm=2.26,
        note="Density 9.68 g/cm3 (Thermo Avantage Densities.txt); band gap 5.7 eV (Thermo Avantage Densities.txt). IMFP "
             "2.26 nm COMPUTED with TPP-2M (Nv = 16) at KE ~1470 eV (Al Ka); not a published value. BE is a typical "
             "literature value, not from your files. Density is the crystalline value; thin amorphous HfO2 films are "
             "lower.",
    ),
    "SnO2": Compound(
        "SnO2", density_g_cm3=6.95, molar_mass_g_mol=150.708,
        atoms_per_formula=1, core_level="Sn 3d5/2", be_ev=486.6, imfp_nm=1.80,
        note="Density 6.95 g/cm3 (standard handbook value, not from your files); band gap 3.6 eV (standard handbook value, "
             "not from your files). IMFP 1.80 nm COMPUTED with TPP-2M (Nv = 16) at KE ~1000 eV (Al Ka); not a published "
             "value. BE is a typical literature value, not from your files. SnO (Sn 3d5/2 ~486.0 eV) is not modelled.",
    ),
    "Ga2O3": Compound(
        "Ga2O3", density_g_cm3=5.88, molar_mass_g_mol=187.443,
        atoms_per_formula=2, core_level="Ga 3d", be_ev=20.7, imfp_nm=2.56,
        note="Density 5.88 g/cm3 (standard handbook value, not from your files); band gap 4.8 eV (standard handbook value, "
             "not from your files). IMFP 2.56 nm COMPUTED with TPP-2M (Nv = 24) at KE ~1466 eV (Al Ka); not a published "
             "value. BE is a typical literature value, not from your files. 2 Ga atoms per formula unit; beta-Ga2O3 "
             "density.",
    ),
    "TeO2": Compound(
        "TeO2", density_g_cm3=6.04, molar_mass_g_mol=159.598,
        atoms_per_formula=1, core_level="Te 3d5/2", be_ev=576, imfp_nm=1.76,
        note="Density 6.04 g/cm3 (standard handbook value, not from your files); band gap 3.7 eV (standard handbook value, "
             "not from your files). IMFP 1.76 nm COMPUTED with TPP-2M (Nv = 18) at KE ~911 eV (Al Ka); not a published "
             "value. BE is a typical literature value, not from your files. 6.04 g/cm3 is tetragonal paratellurite "
             "(orthorhombic tellurite is 5.67) - edit for your film.",
    ),
    "In2O3": Compound(
        "In2O3", density_g_cm3=7.18, molar_mass_g_mol=277.633,
        atoms_per_formula=2, core_level="In 3d5/2", be_ev=444.9, imfp_nm=1.88,
        note="Density 7.18 g/cm3 (standard handbook value, not from your files); band gap 2.9 eV (rough estimate, not from "
             "your files). IMFP 1.88 nm COMPUTED with TPP-2M (Nv = 24) at KE ~1042 eV (Al Ka); not a published value. BE "
             "is a typical literature value, not from your files. 2 In atoms per formula unit. See the InP note on the "
             "small In 3d5/2 shift.",
    ),
    "PbSO4": Compound(
        "PbSO4", density_g_cm3=6.29, molar_mass_g_mol=303.256,
        atoms_per_formula=1, core_level="S 2p3/2", be_ev=168.6, imfp_nm=2.32,
        note="Density 6.29 g/cm3 (standard handbook value, not from your files); band gap 4.5 eV (rough estimate, not from "
             "your files). IMFP 2.32 nm COMPUTED with TPP-2M (Nv = 34) at KE ~1318 eV (Al Ka); not a published value. BE "
             "is a typical literature value, not from your files. 6.29 g/cm3 is anglesite. The 4.5 eV band gap is a ROUGH "
             "estimate (IMFP is only mildly sensitive to it).",
    ),
    "SiON (oxynitride, approx.)": Compound(
        "SiON", density_g_cm3=2.50, molar_mass_g_mol=87.09,
        atoms_per_formula=1, core_level="Si 2p3/2", be_ev=102.3, imfp_nm=3.30,
        note="Rough intermediate values only - oxynitride stoichiometry varies; edit for your film.",
    ),
    "Custom / manual entry": Compound(
        "Custom", density_g_cm3=1.0, molar_mass_g_mol=1.0,
        atoms_per_formula=1, core_level="", be_ev=0.0, imfp_nm=1.0,
        note="Enter your own values.",
    ),
}


# Band gap per preset (eV): metals are 0, None = not known (the thickness equation itself doesn't use it).
_PRESET_BAND_GAP = {
    "InP (In 3d5/2)": BAND_GAPS_EV["InP"][0],
    "PbS (S 2p)": BAND_GAPS_EV["PbS"][0],
    "In2O3": BAND_GAPS_EV["In2O3"][0],
    "PbSO4": BAND_GAPS_EV["PbSO4"][0],
    "Mg (metal)": 0.0,
    "Ge (elemental)": BAND_GAPS_EV["Ge"][0],
    "Zn (metal, LMM)": 0.0,
    "Zr (metal)": 0.0,
    "Hf (metal)": 0.0,
    "Sn (metal)": 0.0,
    "SiC": BAND_GAPS_EV["SiC"][0],
    "GaAs (Ga 3d)": BAND_GAPS_EV["GaAs"][0],
    "GaP (Ga 3d)": BAND_GAPS_EV["GaP"][0],
    "CdTe (Te 3d5/2)": BAND_GAPS_EV["CdTe"][0],
    "MgO": BAND_GAPS_EV["MgO"][0],
    "GeO2": BAND_GAPS_EV["GeO2"][0],
    "ZnO (LMM)": BAND_GAPS_EV["ZnO"][0],
    "ZrO2": BAND_GAPS_EV["ZrO2"][0],
    "HfO2": BAND_GAPS_EV["HfO2"][0],
    "SnO2": BAND_GAPS_EV["SnO2"][0],
    "Ga2O3": BAND_GAPS_EV["Ga2O3"][0],
    "TeO2": BAND_GAPS_EV["TeO2"][0],
    "Cu (metal, LMM)": 0.0, "Cu (metal, 2p3/2)": 0.0, "W (metal)": 0.0, "Al (metal)": 0.0, "Al (metal, Mg Ka)": 0.0,
    "Ti (metal)": 0.0, "Si (metal/elemental)": BAND_GAPS_EV["Si"][0],
    "CuO (LMM)": BAND_GAPS_EV["CuO"][0], "CuO (2p3/2)": BAND_GAPS_EV["CuO"][0],
    "Cu2O (LMM)": BAND_GAPS_EV["Cu2O"][0], "Cu2O (2p3/2)": BAND_GAPS_EV["Cu2O"][0],
    "Cu(OH)2 (LMM)": BAND_GAPS_EV["Cu(OH)2"][0], "Cu(OH)2 (2p3/2)": BAND_GAPS_EV["Cu(OH)2"][0],
    "Al2O3": BAND_GAPS_EV["Al2O3"][0], "Al2O3 (Mg Ka)": BAND_GAPS_EV["Al2O3"][0],
    "SiO2": BAND_GAPS_EV["SiO2"][0], "Si3N4": BAND_GAPS_EV["Si3N4"][0], "TiO2": BAND_GAPS_EV["TiO2"][0],
}
for _presets in (SUBSTRATES, OVERLAYERS):
    for _key, _compound in _presets.items():
        if _compound is not None and _key in _PRESET_BAND_GAP:
            _compound.band_gap_ev = _PRESET_BAND_GAP[_key]


# Inputs behind every COMPUTED (TPP-2M) IMFP preset: name -> (valence electrons per formula unit, KE in eV).
# Density, molar mass and band gap come from the preset itself, so each imfp_nm can be regenerated and tested.
TPP2M_INPUTS = {
    "InP (In 3d5/2)": (8, 1042.0),
    "PbS (S 2p)": (10, 1325.6),
    "In2O3": (24, 1041.7),
    "PbSO4": (34, 1318.0),
    "Mg (metal)": (2, 1436.9),
    "Ge (elemental)": (4, 1457.2),
    "Zn (metal, LMM)": (2, 992.0),
    "Zr (metal)": (4, 1307.8),
    "Hf (metal)": (4, 1472.4),
    "Sn (metal)": (4, 1001.7),
    "SiC": (8, 1386.2),
    "GaAs (Ga 3d)": (8, 1467.3),
    "GaP (Ga 3d)": (8, 1467.4),
    "CdTe (Te 3d5/2)": (8, 914.5),
    "MgO": (8, 1435.8),
    "GeO2": (16, 1453.7),
    "ZnO (LMM)": (8, 988.0),
    "ZrO2": (16, 1304.3),
    "HfO2": (16, 1469.8),
    "SnO2": (16, 1000.0),
    "Ga2O3": (24, 1465.9),
    "TeO2": (18, 910.6),
    "Cu (metal, 2p3/2)": (11, 554.0),
    "CuO (2p3/2)": (17, 554.0),
    "Cu2O (2p3/2)": (28, 554.0),
    "Cu(OH)2 (2p3/2)": (25, 554.0),
    "Si (metal/elemental)": (4, 1387.3),
    "SiO2": (16, 1387.3),
    "Si3N4": (32, 1387.3),
    "Ti (metal)": (4, 1032.7),
    "TiO2": (16, 1032.7),
}


def volume_density(density_g_cm3: float, molar_mass_g_mol: float, atoms_per_formula: float) -> float:
    """Nx per tech note: atoms_per_formula * density / molar_mass (mol cm^-3 equivalent)."""
    if molar_mass_g_mol <= 0:
        raise ValueError("Molar mass must be > 0")
    return atoms_per_formula * density_g_cm3 / molar_mass_g_mol


def tpp2m_imfp_nm(ke_ev: float, density_g_cm3: float, molar_mass_g_mol: float,
                  valence_electrons: float, band_gap_ev: float = 0.0) -> float:
    """
    TPP-2M predictive IMFP (Tanuma, Powell & Penn, Surf. Interface Anal. 1994, 21, 165), in nm.
    valence_electrons = Nv per formula unit (e.g. Cu 11, CuO 17, Cu2O 28, Cu(OH)2 25); band gap 0 for metals.
    Used to derive the Cu 2p3/2 presets; not wired into the GUI - edit the IMFP fields by hand.
    """
    u = valence_electrons * density_g_cm3 / molar_mass_g_mol
    e_p = 28.8 * math.sqrt(u)
    beta = -0.10 + 0.944 / math.sqrt(e_p ** 2 + band_gap_ev ** 2) + 0.069 * density_g_cm3 ** 0.1
    gamma = 0.191 * density_g_cm3 ** -0.5
    c = 1.97 - 0.91 * u
    d = 53.4 - 20.8 * u
    lam_angstrom = ke_ev / (e_p ** 2 * (beta * math.log(gamma * ke_ev) - c / ke_ev + d / ke_ev ** 2))
    return lam_angstrom / 10.0


@dataclass
class OverlayerInput:
    label: str
    N: float
    lam: float
    rsf: float
    area: float


def strohmeier_thickness(theta_deg: float, N_sub: float, lam_sub: float, rsf_sub: float,
                          I_sub: float, components: list) -> dict:
    """
    General Strohmeier-type calculation.

    components: list of OverlayerInput, each with its own N, lambda, RSF and
    peak area. Reduces exactly to the single-component Strohmeier equation
    when len(components) == 1 and RSFs == 1, and to the tech note's
    oxide+hydroxide mixture equation when there are 2+ components with RSFs == 1.
    """
    active = [c for c in components if c.area > 0]
    if not active:
        raise ValueError("Enter at least one overlayer peak area > 0.")
    if I_sub <= 0:
        raise ValueError("Substrate (metal) peak area must be > 0.")
    if not (0 < theta_deg <= 90):
        raise ValueError("Take-off angle must be between 0 and 90 degrees.")

    sum_I = sum(c.area for c in active)
    N_ov_avg = sum(c.N * c.area for c in active) / sum_I
    lam_ov_avg = sum(c.lam * c.area for c in active) / sum_I
    rsf_ov_avg = sum(c.rsf * c.area for c in active) / sum_I

    theta = math.radians(theta_deg)
    numerator = rsf_sub * N_sub * lam_sub * sum_I
    denominator = rsf_ov_avg * N_ov_avg * lam_ov_avg * I_sub
    ratio = numerator / denominator
    d_nm = lam_ov_avg * math.sin(theta) * math.log(ratio + 1.0)

    return {
        "d_nm": d_nm,
        "d_angstrom": d_nm * 10.0,
        "N_ov_avg": N_ov_avg,
        "lam_ov_avg": lam_ov_avg,
        "rsf_ov_avg": rsf_ov_avg,
        "N_sub": N_sub,
        "N_ratio_sub_over_ov": N_sub / N_ov_avg if N_ov_avg else float("nan"),
        "validity_limit_nm": 3.0 * lam_ov_avg,
        "within_validity": d_nm <= 3.0 * lam_ov_avg,
        "sum_I_overlayer": sum_I,
    }


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

N_OVERLAYER_ROWS = 4

FIELD_WIDTH = 8

# Palette (one place, so the look is easy to tweak)
BG = "#eef1f5"
CARD = "#ffffff"
INK = "#1f2937"
MUTED = "#6b7280"
ACCENT = "#1f6feb"
ACCENT_HOVER = "#1a5fd0"
ACCENT_PRESSED = "#174fb0"
BORDER = "#8d96a3"
BTN = "#e5e7eb"
BTN_HOVER = "#d6dbe2"
BTN_PRESSED = "#c3cad3"
BADGE = {   # kind -> (background, foreground)
    "idle": ("#e5e7eb", "#4b5563"),
    "ok": ("#d1fae5", "#065f46"),
    "warn": ("#fef3c7", "#92400e"),
    "bad": ("#fee2e2", "#991b1b"),
}


def setup_style(root, compact=False):
    """Switch to the 'clam' theme and define the button/card/entry styles used by the app."""
    style = ttk.Style(root)
    try:
        style.theme_use("clam")
    except tk.TclError:
        pass
    base = 8 if compact else 9
    font = ("Segoe UI", base)
    root.configure(bg=BG)

    style.configure(".", background=BG, foreground=INK, font=font, bordercolor=BORDER, focuscolor=ACCENT)
    style.configure("TFrame", background=BG)
    style.configure("Card.TFrame", background=CARD)
    style.configure("TLabel", background=BG, foreground=INK)
    style.configure("Card.TLabel", background=CARD, foreground=INK)
    style.configure("Muted.TLabel", background=BG, foreground=MUTED)
    style.configure("CardMuted.TLabel", background=CARD, foreground=MUTED)
    style.configure("Header.TLabel", background=BG, foreground=INK, font=("Segoe UI", base + 4, "bold"))
    style.configure("ColHead.TLabel", background=CARD, foreground=MUTED, font=("Segoe UI", base - 1, "bold"))
    style.configure("RowHead.TLabel", background=CARD, foreground=INK, font=("Segoe UI", base, "bold"))
    style.configure("Card.TCheckbutton", background=CARD, foreground=INK)
    style.map("Card.TCheckbutton", background=[("active", CARD)])

    style.configure("Card.TLabelframe", background=CARD, bordercolor=BORDER, relief="solid", borderwidth=1)
    style.configure("Card.TLabelframe.Label", background=CARD, foreground=ACCENT, font=("Segoe UI", base + 1, "bold"))
    style.configure("TSeparator", background="#d1d5db")

    # Buttons: a clearly bordered, raised secondary style and a bold blue primary style, with hover/pressed states.
    style.configure("TButton", background=BTN, foreground=INK, bordercolor=BORDER, lightcolor="#ffffff",
                    darkcolor="#c1c7cf", relief="raised", borderwidth=1, padding=(10, 5), focuscolor=BTN, font=font)
    style.map("TButton",
              background=[("disabled", "#f3f4f6"), ("pressed", BTN_PRESSED), ("active", BTN_HOVER)],
              foreground=[("disabled", "#9ca3af")],
              relief=[("pressed", "sunken"), ("!pressed", "raised")])
    style.configure("Primary.TButton", background=ACCENT, foreground="#ffffff", bordercolor=ACCENT_PRESSED,
                    lightcolor="#4a8df0", darkcolor=ACCENT_PRESSED, relief="raised", borderwidth=1,
                    padding=(18, 9), focuscolor=ACCENT, font=("Segoe UI", base + 2, "bold"))
    style.map("Primary.TButton",
              background=[("disabled", "#9db9e8"), ("pressed", ACCENT_PRESSED), ("active", ACCENT_HOVER)],
              foreground=[("disabled", "#ffffff")],
              relief=[("pressed", "sunken"), ("!pressed", "raised")])

    # Input fields
    style.configure("TEntry", fieldbackground="#ffffff", bordercolor=BORDER, lightcolor=BORDER, darkcolor=BORDER,
                    padding=3)
    style.map("TEntry", bordercolor=[("focus", ACCENT)], lightcolor=[("focus", ACCENT)],
              darkcolor=[("focus", ACCENT)])
    style.configure("TCombobox", fieldbackground="#ffffff", bordercolor=BORDER, lightcolor=BORDER, darkcolor=BORDER,
                    arrowcolor=INK, padding=3)
    style.map("TCombobox", fieldbackground=[("readonly", "#ffffff")], selectbackground=[("readonly", "#ffffff")],
              selectforeground=[("readonly", INK)], bordercolor=[("focus", ACCENT)])

    style.configure("TNotebook", background=BG, borderwidth=0, tabmargins=(2, 4, 2, 0))
    style.configure("TNotebook.Tab", background="#dfe4ea", foreground=MUTED, padding=(18, 7),
                    font=("Segoe UI", base + 1, "bold"))
    style.map("TNotebook.Tab", background=[("selected", CARD), ("active", "#e9edf2")],
              foreground=[("selected", ACCENT)], padding=[("selected", (18, 9))])

    style.configure("Horizontal.TScale", background=CARD, troughcolor="#d1d5db")
    style.configure("TPanedwindow", background=BG)
    style.configure("Sash", sashthickness=10, gripcount=16, handlesize=14)


class Tooltip:
    """Small hover tooltip: appears after a short delay, disappears on leave/click."""

    def __init__(self, widget, text, delay=450):
        self.widget, self.text, self.delay = widget, text, delay
        self._tip = None
        self._after = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._cancel()
        self._after = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after is not None:
            self.widget.after_cancel(self._after)
            self._after = None

    def _show(self):
        if self._tip is not None or not self.text:
            return
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")
        tk.Label(self._tip, text=self.text, justify="left", background="#111827", foreground="#ffffff",
                 padx=8, pady=4, font=("Segoe UI", 8), wraplength=320).pack()

    def _hide(self, _event=None):
        self._cancel()
        if self._tip is not None:
            self._tip.destroy()
            self._tip = None


def _to_float(text):
    """Float from an entry's text, or None if blank/invalid."""
    try:
        return float(str(text).strip())
    except ValueError:
        return None


class CompoundRow:
    """One editable row (substrate or overlayer) bound to a Compound preset dropdown."""

    def __init__(self, parent, row_idx, choices, default_key, is_substrate=False, pady=2):
        self.is_substrate = is_substrate
        self.active_var = tk.BooleanVar(value=is_substrate)  # substrate always active
        self.choice_var = tk.StringVar(value=default_key)

        c = 0
        if not is_substrate:
            self.active_chk = ttk.Checkbutton(parent, variable=self.active_var, style="Card.TCheckbutton")
            self.active_chk.grid(row=row_idx, column=c, padx=(2, 4), pady=pady)
            Tooltip(self.active_chk, "Tick to include this component in the (mixed) overlayer")
        c += 1

        self.combo = ttk.Combobox(parent, textvariable=self.choice_var, values=list(choices.keys()),
                                   width=21, state="readonly")
        self.combo.grid(row=row_idx, column=c, padx=2, pady=pady, sticky="w")
        c += 1

        self.density_var = tk.StringVar()
        self.molar_var = tk.StringVar()
        self.atoms_var = tk.StringVar()
        self.be_var = tk.StringVar()
        self.core_var = tk.StringVar()
        self.imfp_var = tk.StringVar()
        self.rsf_var = tk.StringVar(value="1")
        self.area_var = tk.StringVar()

        entries = [
            (self.density_var, FIELD_WIDTH),
            (self.molar_var, FIELD_WIDTH),
            (self.atoms_var, 5),
            (self.core_var, 9),
            (self.be_var, 7),
            (self.imfp_var, 7),
            (self.rsf_var, 5),
            (self.area_var, FIELD_WIDTH),
        ]
        self.entry_widgets = []
        for var, width in entries:
            w = ttk.Entry(parent, textvariable=var, width=width)
            w.grid(row=row_idx, column=c, padx=2, pady=pady)
            self.entry_widgets.append(w)
            c += 1

        self.combo.bind("<<ComboboxSelected>>", lambda e: self._apply_preset(choices))
        self._apply_preset(choices)

    def _apply_preset(self, choices):
        key = self.choice_var.get()
        comp = choices.get(key)
        if comp is None:
            # "None" (overlayer disabled) - leave fields blank/untouched
            if key == "None":
                self.active_var.set(False)
            return
        if not self.is_substrate:
            self.active_var.set(True)
        self.density_var.set(f"{comp.density_g_cm3:g}")
        self.molar_var.set(f"{comp.molar_mass_g_mol:g}")
        self.atoms_var.set(f"{comp.atoms_per_formula:g}")
        self.be_var.set(f"{comp.be_ev:g}" if comp.be_ev else "")
        self.core_var.set(comp.core_level)
        self.imfp_var.set(f"{comp.imfp_nm:g}")
        self.rsf_var.set(f"{comp.rsf:g}")

    def read(self):
        """Return (label, N, lambda, rsf, area) or raise ValueError with a helpful message."""
        label = self.choice_var.get()
        try:
            density = float(self.density_var.get())
            molar = float(self.molar_var.get())
            atoms = float(self.atoms_var.get())
            lam = float(self.imfp_var.get())
            rsf = float(self.rsf_var.get()) if self.rsf_var.get().strip() else 1.0
            area = float(self.area_var.get()) if self.area_var.get().strip() else 0.0
        except ValueError:
            raise ValueError(f"Row '{label}': please enter numeric values for density, "
                              f"molar mass, atoms/formula, IMFP, RSF and peak area.")
        N = volume_density(density, molar, atoms)
        return label, N, lam, rsf, area


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("XPS Overlayer Thickness Calculator (Strohmeier method)")
        width, height = self._initial_geometry()
        self.compact = height < 880       # laptops: smaller fonts, tighter rows, no hint paragraph
        setup_style(self, compact=self.compact)
        x = max(0, (self.winfo_screenwidth() - width) // 2)
        y = max(0, (self.winfo_screenheight() - height) // 2 - 20)
        self.geometry(f"{width}x{height}+{x}+{y}")
        self.minsize(min(1100, width), min(680, height))

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=8)

        self.calc_tab = ttk.Frame(nb)
        self.help_tab = ttk.Frame(nb)
        nb.add(self.calc_tab, text="Calculator & simulation")
        nb.add(self.help_tab, text="Method / References")

        self.history = []  # list of dict rows for CSV export
        self.last_calc = None   # inputs + result of the last successful calculation (feeds the spectrum pane)
        self.sim_data = None    # (x, y, name) of a user-loaded spectrum

        self._build_calc_tab()
        self._build_help_tab()
        self.after(150, self._place_sash)

    def _initial_geometry(self):
        """Window size from the detected screen: up to 1500 x 980, otherwise ~95% of the screen."""
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        return min(int(sw * 0.95), 1500), min(int(sh * 0.93), 980)

    def _place_sash(self):
        """Give the top pane exactly the height its widgets need (capped), the rest goes to the spectrum."""
        try:
            self.update_idletasks()
            wanted = self.top_pane.winfo_reqheight() + 6
            cap = int(self.paned.winfo_height() * (0.50 if self.compact else 0.58))
            self.paned.sashpos(0, max(160, min(wanted, cap)))
        except tk.TclError:
            pass

    # ------------------------------------------------------------------
    def _build_calc_tab(self):
        if not self.compact:      # on short screens the window title already names the app; keep the height for the plot
            header = ttk.Frame(self.calc_tab)
            header.pack(fill="x", padx=10, pady=(8, 2))
            ttk.Label(header, text="Overlayer thickness from XPS peak areas", style="Header.TLabel").pack(side="left")
            ttk.Label(header, text="Strohmeier equation (1990)", style="Muted.TLabel").pack(side="left", padx=14,
                                                                                               pady=(6, 0))

        self.paned = ttk.PanedWindow(self.calc_tab, orient="vertical")
        self.paned.pack(fill="both", expand=True, padx=8, pady=(4, 8))
        self.top_pane = ttk.Frame(self.paned)
        self.bottom_pane = ttk.Frame(self.paned)
        self.paned.add(self.top_pane, weight=0)
        self.paned.add(self.bottom_pane, weight=1)

        self.top_pane.columnconfigure(0, weight=3, minsize=760)
        self.top_pane.columnconfigure(1, weight=2, minsize=380)
        self.top_pane.rowconfigure(0, weight=1)
        self._build_inputs_card(self.top_pane).grid(row=0, column=0, sticky="nsew", padx=(0, 8), pady=(0, 4))
        self._build_results_card(self.top_pane).grid(row=0, column=1, sticky="nsew", pady=(0, 4))
        self._build_spectrum_pane(self.bottom_pane)

    # --- inputs (left) -------------------------------------------------
    def _build_inputs_card(self, parent):
        card = ttk.LabelFrame(parent, text="  Inputs  ", style="Card.TLabelframe", padding=(10, 4, 10, 8))

        card.columnconfigure(0, weight=1)
        pady = 1 if self.compact else 2
        if not self.compact:
            hint = ttk.Label(card, style="CardMuted.TLabel", justify="left",
                             text="Tick the overlayer component(s) present (e.g. oxide + hydroxide), enter their fitted "
                                  "peak areas, then Calculate. Every preset value is editable - verify them for your "
                                  "material and instrument.")
            hint.grid(row=0, column=0, sticky="ew", pady=(0, 6))
            hint.bind("<Configure>", lambda e: hint.configure(wraplength=max(200, e.width - 4)))

        table = ttk.Frame(card, style="Card.TFrame")
        table.grid(row=2, column=0, sticky="w")

        angle = ttk.Frame(table if self.compact else card, style="Card.TFrame")
        if self.compact:    # fold the angle into the Substrate row: it saves a whole row of height on short screens
            angle.grid(row=1, column=2, columnspan=8, sticky="e", pady=(4, 0))
        else:
            angle.grid(row=1, column=0, sticky="w", pady=(0, 6))
        ttk.Label(angle, style="Card.TLabel",
                  text=("Take-off angle θ (deg from surface; 90 = normal):" if self.compact
                        else "Take-off angle θ (deg from the surface; 90 = normal emission):")).pack(side="left")
        self.theta_var = tk.StringVar(value="90")
        theta_entry = ttk.Entry(angle, textvariable=self.theta_var, width=6)
        theta_entry.pack(side="left", padx=6)
        theta_entry.bind("<Return>", lambda e: self._on_calculate())

        headers =["", "Compound", "Density\n(g/cm³)", "Molar mass\n(g/mol)", "Atoms/\nformula",
                   "Core level", "BE\n(eV)", "IMFP\n(nm)", "RSF", "Peak area\n(I)"]
        for i, h in enumerate(headers):
            ttk.Label(table, text=h, style="ColHead.TLabel", justify="center").grid(row=0, column=i, padx=2)

        ttk.Label(table, text="Substrate", style="RowHead.TLabel").grid(row=1, column=0, columnspan=2, sticky="w",
                                                                         pady=(4, 0))
        self.substrate_row = CompoundRow(table, 2, SUBSTRATES, "Cu (metal, LMM)", is_substrate=True, pady=pady)

        if not self.compact:
            ttk.Separator(table, orient="horizontal").grid(row=3, column=0, columnspan=10, sticky="ew", pady=6)
        ttk.Label(table, text="Overlayer component(s)  (tick to include)", style="RowHead.TLabel").grid(
            row=4, column=0, columnspan=6, sticky="w", pady=(4, 0) if self.compact else 0)

        self.overlayer_rows = []
        defaults = ["CuO (LMM)", "Cu2O (LMM)", "Cu(OH)2 (LMM)", "None"]
        for i in range(N_OVERLAYER_ROWS):
            row = CompoundRow(table, 5 + i, OVERLAYERS, defaults[i] if i < len(defaults) else "None",
                               is_substrate=False, pady=pady)
            self.overlayer_rows.append(row)

        for row in [self.substrate_row] + self.overlayer_rows:
            for w in row.entry_widgets:
                w.bind("<Return>", lambda e: self._on_calculate())
        return card

    # --- calculation (right) ---------------------------------------------
    def _build_results_card(self, parent):
        card = ttk.LabelFrame(parent, text="  Calculation  ", style="Card.TLabelframe", padding=10)
        card.columnconfigure(0, weight=1)
        card.rowconfigure(3, weight=1)

        self.calc_btn = ttk.Button(card, text="Calculate thickness", style="Primary.TButton",
                                   command=self._on_calculate)
        Tooltip(self.calc_btn, "Calculate the overlayer thickness from the peak areas entered on the left "
                               "(pressing Enter in any input box does the same)")

        readout = ttk.Frame(card, style="Card.TFrame")
        if self.compact:    # button and readout share a row to save height on short screens
            self.calc_btn.grid(row=0, column=0, sticky="ew")
            readout.grid(row=0, column=1, sticky="e", padx=(12, 0))
        else:
            self.calc_btn.grid(row=0, column=0, columnspan=2, sticky="ew")
            readout.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0))
        self.readout_var = tk.StringVar(value="-.-- nm")
        self.readout_sub_var = tk.StringVar(value="Overlayer thickness will appear here")
        self.readout_lbl = tk.Label(readout, textvariable=self.readout_var,
                                    font=("Segoe UI", 30 if not self.compact else 24, "bold"), fg="#b8bec7", bg=CARD,
                                    anchor="w")
        self.readout_lbl.pack(anchor="w")
        tk.Label(readout, textvariable=self.readout_sub_var, font=("Segoe UI", 9), fg=MUTED, bg=CARD,
                 anchor="w").pack(anchor="w")

        self.badge = tk.Label(card, text="No calculation yet", font=("Segoe UI", 9, "bold"), padx=8, pady=4,
                              anchor="w", justify="left", wraplength=360)
        self.badge.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(6, 6) if self.compact else (8, 8))
        self._set_badge("idle", "No calculation yet")

        text_frame = ttk.Frame(card, style="Card.TFrame")
        text_frame.grid(row=3, column=0, columnspan=2, sticky="nsew")
        text_frame.columnconfigure(0, weight=1)
        text_frame.rowconfigure(0, weight=1)
        self.result_text = tk.Text(text_frame, height=5, width=44, wrap="word", state="disabled",
                                   font=("Consolas", 9 if not self.compact else 8), relief="flat", bg="#f8fafc",
                                   highlightthickness=1, highlightbackground="#d1d5db", padx=6, pady=4)
        scroll = ttk.Scrollbar(text_frame, orient="vertical", command=self.result_text.yview)
        self.result_text.configure(yscrollcommand=scroll.set)
        self.result_text.grid(row=0, column=0, sticky="nsew")
        scroll.grid(row=0, column=1, sticky="ns")

        buttons = ttk.Frame(card, style="Card.TFrame")
        buttons.grid(row=4, column=0, columnspan=2, sticky="ew", pady=(6, 0) if self.compact else (8, 0))
        buttons.columnconfigure((0, 1), weight=1, uniform="b")
        clear_btn = ttk.Button(buttons, text="Clear peak areas", command=self._clear_areas)
        clear_btn.grid(row=0, column=0, sticky="ew", padx=(0, 4))
        export_btn = ttk.Button(buttons, text="Export history to CSV...", command=self._export_csv)
        export_btn.grid(row=0, column=1, sticky="ew", padx=(4, 0))
        Tooltip(clear_btn, "Empty all the peak-area boxes (the presets stay)")
        Tooltip(export_btn, "Save every calculation from this session to a CSV file")
        return card

    def _set_badge(self, kind, text):
        bg, fg = BADGE[kind]
        self.badge.configure(text=text, bg=bg, fg=fg)

    def _clear_areas(self):
        self.substrate_row.area_var.set("")
        for row in self.overlayer_rows:
            row.area_var.set("")

    def _on_calculate(self):
        try:
            theta = float(self.theta_var.get())
            sub_label, N_sub, lam_sub, rsf_sub, I_sub = self.substrate_row.read()

            components = []
            comp_meta = []   # (binding energy text, core level text) per component, for the spectrum pane
            for row in self.overlayer_rows:
                if not row.active_var.get() or row.choice_var.get() == "None":
                    continue
                label, N, lam, rsf, area = row.read()
                if area <= 0:
                    continue
                components.append(OverlayerInput(label, N, lam, rsf, area))
                comp_meta.append((row.be_var.get(), row.core_var.get()))

            if not components:
                raise ValueError("Tick at least one overlayer row as active and enter its peak area.")

            res = strohmeier_thickness(theta, N_sub, lam_sub, rsf_sub, I_sub, components)
        except ValueError as e:
            messagebox.showerror("Input error", str(e))
            return
        except ZeroDivisionError:
            messagebox.showerror("Input error", "Division by zero - check densities/molar masses/areas.")
            return

        comp_names = ", ".join(f"{c.label} (I={c.area:g})" for c in components)
        limit_txt = "within" if res["within_validity"] else "EXCEEDS"
        lines = []      # the result comes first: the box is small and the inputs are visible on the left anyway
        lines.append(f"==> Overlayer thickness d = {res['d_nm']:.3f} nm  ({res['d_angstrom']:.2f} Angstrom)")
        lines.append(f"Validity: reliable to ~3x lambda_ov = {res['validity_limit_nm']:.2f} nm; "
                     f"d is {limit_txt} this limit.")
        if not res["within_validity"]:
            lines.append("  -> May be unreliable; consider angle-resolved XPS, depth profiling, ellipsometry or TEM.")
        lines.append("-" * 44)
        lines.append(f"Substrate : {sub_label}")
        lines.append(f"  N_sub = {N_sub:.5f}  lambda = {lam_sub:g} nm  I = {I_sub:g}")
        lines.append(f"Overlayer : {comp_names}")
        lines.append(f"  N_ov (weighted)      = {res['N_ov_avg']:.5f}")
        lines.append(f"  lambda_ov (weighted) = {res['lam_ov_avg']:.4f} nm")
        lines.append(f"  N_sub / N_ov         = {res['N_ratio_sub_over_ov']:.4f}")
        lines.append(f"Take-off angle theta = {theta:g} deg")

        self.result_text.configure(state="normal")
        self.result_text.delete("1.0", "end")
        self.result_text.insert("1.0", "\n".join(lines))
        self.result_text.configure(state="disabled")

        limit = res["validity_limit_nm"]
        self.readout_var.set(f"{res['d_nm']:.2f} nm")
        self.readout_lbl.configure(fg=ACCENT)
        self.readout_sub_var.set(f"{res['d_angstrom']:.1f} Å     (θ = {theta:g}°)")
        if res["d_nm"] > limit:
            self._set_badge("bad", f"Exceeds the ~3λ limit ({limit:.1f} nm): result may be unreliable")
        elif res["d_nm"] > 0.8 * limit:
            self._set_badge("warn", f"Close to the ~3λ limit ({limit:.1f} nm): treat with caution")
        else:
            self._set_badge("ok", f"Within the validity limit (~3λ = {limit:.1f} nm)")

        self.history.append({
            "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "substrate": sub_label,
            "I_sub": I_sub,
            "overlayer_components": comp_names,
            "theta_deg": theta,
            "thickness_nm": round(res["d_nm"], 4),
            "thickness_angstrom": round(res["d_angstrom"], 3),
            "validity_limit_nm": round(res["validity_limit_nm"], 3),
            "within_validity": res["within_validity"],
        })

        self.last_calc = {
            "theta": theta, "res": res,
            "sub": {"label": sub_label, "N": N_sub, "lam": lam_sub, "rsf": rsf_sub,
                    "be": _to_float(self.substrate_row.be_var.get()), "core": self.substrate_row.core_var.get()},
            "comps": [{"label": c.label, "area": c.area, "be": _to_float(be), "core": core}
                      for c, (be, core) in zip(components, comp_meta)],
        }
        self._refresh_spectrum(reset=True)

    def _export_csv(self):
        if not self.history:
            messagebox.showinfo("No data", "Run at least one calculation first.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv",
                                             filetypes=[("CSV files", "*.csv")],
                                             initialfile="xps_thickness_results.csv")
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(self.history[0].keys()))
            writer.writeheader()
            writer.writerows(self.history)
        messagebox.showinfo("Exported", f"Saved {len(self.history)} result(s) to:\n{path}")

    # ------------------------------------------------------------------
    # Simulated spectrum pane (needs numpy + matplotlib)
    # ------------------------------------------------------------------
    OXIDE_COLOURS = ["#d62728", "#ff7f0e", "#2ca02c", "#9467bd"]
    METAL_COLOUR = "#1f77b4"

    def _build_spectrum_pane(self, parent):
        pane = ttk.LabelFrame(parent, text="  Simulated spectrum (estimate)  ", style="Card.TLabelframe",
                              padding=(10, 4, 10, 6))
        pane.pack(fill="both", expand=True)
        self.spec_pane = pane
        if not HAVE_PLOT:
            ttk.Label(pane, style="Card.TLabel", justify="left",
                      text="The simulated spectrum needs numpy and matplotlib:\n\n"
                           "    pip install numpy matplotlib\n\n"
                           "then restart the app. The thickness calculator works without them.").pack(anchor="nw", pady=10)
            return

        self.sim_d_var = tk.StringVar()
        self.sim_shift_var = tk.StringVar(value="0")
        self.sim_fwhm_m_var = tk.StringVar()
        self.sim_fwhm_c_var = tk.StringVar()
        self.sim_status = tk.StringVar(value="Run a calculation to see an estimated spectrum here.")
        self._sim_busy = False

        ctl = ttk.Frame(pane, style="Card.TFrame")
        ctl.pack(fill="x")
        ttk.Label(ctl, text="Thickness (nm):", style="Card.TLabel").pack(side="left")
        d_entry = ttk.Entry(ctl, textvariable=self.sim_d_var, width=7)
        d_entry.pack(side="left", padx=(4, 6))
        self.sim_scale = ttk.Scale(ctl, from_=0.0, to=10.0, length=130, command=self._on_slider)
        self.sim_scale.pack(side="left", padx=(0, 6))
        reset_btn = ttk.Button(ctl, text="Reset", command=lambda: self._refresh_spectrum(reset=True))
        reset_btn.pack(side="left", padx=(0, 14))

        entries = [d_entry]
        ttk.Label(ctl, text="BE shift (eV):", style="Card.TLabel").pack(side="left", padx=(0, 2))
        shift_entry = ttk.Entry(ctl, textvariable=self.sim_shift_var, width=5)
        shift_entry.pack(side="left", padx=(0, 12))
        ttk.Label(ctl, text="FWHM (eV):  metal", style="Card.TLabel").pack(side="left", padx=(0, 2))
        fm_entry = ttk.Entry(ctl, textvariable=self.sim_fwhm_m_var, width=5)
        fm_entry.pack(side="left", padx=(0, 6))
        ttk.Label(ctl, text="oxide", style="Card.TLabel").pack(side="left", padx=(0, 2))
        fc_entry = ttk.Entry(ctl, textvariable=self.sim_fwhm_c_var, width=5)
        fc_entry.pack(side="left", padx=(0, 14))
        entries += [shift_entry, fm_entry, fc_entry]
        load_btn = ttk.Button(ctl, text="Load my spectrum...", command=self._load_user_spectrum)
        load_btn.pack(side="left")
        clear_btn = ttk.Button(ctl, text="Clear data", command=self._clear_user_spectrum)
        clear_btn.pack(side="left", padx=6)
        for widget in entries:
            widget.bind("<Return>", lambda e: self._draw_spectrum())
            widget.bind("<FocusOut>", lambda e: self._draw_spectrum())
        Tooltip(d_entry, "What-if thickness: type a value or drag the slider")
        Tooltip(self.sim_scale, "Drag to see the spectrum for other thicknesses (what-if)")
        Tooltip(reset_btn, "Set the what-if thickness back to the calculated value")
        Tooltip(shift_entry, "Shift the simulated peaks, e.g. to match your charge correction")
        Tooltip(fm_entry, "Width (FWHM) of each metal peak - match your instrument resolution")
        Tooltip(fc_entry, "Width (FWHM) of each oxide/compound peak")
        Tooltip(load_btn, "Overlay your own measured spectrum (2-column text/CSV: binding energy, intensity)")
        Tooltip(clear_btn, "Remove the overlaid spectrum")

        info = ttk.Frame(pane, style="Card.TFrame")
        info.pack(fill="x", pady=(3, 0))
        self.fig = Figure(figsize=(9.5, 4.4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=pane)
        # Toolbar first (pinned right), then the status text takes the remaining width. The label asks for a tiny
        # width so its unwrapped text can never squeeze the toolbar; wraplength follows the width it is given.
        toolbar = NavigationToolbar2Tk(self.canvas, info, pack_toolbar=False)
        toolbar.set_message = lambda s: None          # a coordinate read-out means nothing on a relative-intensity plot
        for child in [toolbar] + list(toolbar.winfo_children()):
            try:
                child.configure(bg=CARD)
            except tk.TclError:
                pass
        toolbar.configure(width=250)                  # matplotlib sizes the toolbar frame to the whole figure width
        toolbar.pack(side="right")
        status = ttk.Label(info, textvariable=self.sim_status, style="CardMuted.TLabel", justify="left", width=10)
        status.pack(side="left", fill="x", expand=True)
        status.bind("<Configure>", lambda e: status.configure(wraplength=max(200, e.width - 8)))
        widget = self.canvas.get_tk_widget()
        widget.pack(fill="both", expand=True, pady=(2, 0))
        self._resize_job = None
        self._last_canvas_size = (0, 0)
        # add="+" keeps matplotlib's own <Configure> handler (it resizes the figure); a plain bind would replace it.
        widget.bind("<Configure>", self._on_canvas_resize, add="+")
        self._draw_placeholder("Run a calculation to see the estimated spectrum.")

    def _on_canvas_resize(self, event):
        """Redraw (debounced) when the pane is resized so the pixel-based margins stay right."""
        size = (event.width, event.height)
        if size == self._last_canvas_size or self.last_calc is None:
            self._last_canvas_size = size
            return
        self._last_canvas_size = size
        if self._resize_job is not None:
            self.after_cancel(self._resize_job)
        self._resize_job = self.after(150, self._draw_spectrum)

    def _draw_placeholder(self, message):
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True, fontsize=11, transform=ax.transAxes)
        self.canvas.draw_idle()

    def _on_slider(self, value):
        if self._sim_busy:
            return
        self.sim_d_var.set(f"{float(value):.3f}")
        self._draw_spectrum()

    def _simulation_problem(self):
        """Return (model, None) if the last calculation can be drawn, else (None, explanation)."""
        lc = self.last_calc
        model = sim.model_for(lc["sub"]["core"])
        if model is None:
            return None, ("No spectrum can be simulated for '%s'.\nOnly core-level presets with a binding energy are "
                          "supported; Auger (LMM) presets and custom rows without a recognised core level are not."
                          % (lc["sub"]["core"] or "blank core level"))
        for item in [lc["sub"]] + lc["comps"]:
            if not item["be"]:
                return None, f"'{item['label']}' has no binding energy - fill in the BE column to simulate."
        for comp in lc["comps"]:
            if sim.model_for(comp["core"]) is not model:
                return None, (f"Substrate ({lc['sub']['core']}) and overlayer '{comp['label']}' ({comp['core']}) use "
                              "different core levels, so they can't be drawn in one window.")
        return model, None

    def _refresh_spectrum(self, reset=False):
        if not HAVE_PLOT or self.last_calc is None:
            return
        model, problem = self._simulation_problem()
        if problem:
            self.sim_status.set(problem.replace("\n", " "))
            self._draw_placeholder(problem)
            return
        if reset:
            res = self.last_calc["res"]
            d_max = max(3.0 * res["lam_ov_avg"], 1.5 * res["d_nm"], 1.0)
            self._sim_busy = True
            self.sim_scale.configure(to=d_max)
            self.sim_scale.set(res["d_nm"])
            self._sim_busy = False
            self.sim_d_var.set(f"{res['d_nm']:.3f}")
            self.sim_fwhm_m_var.set(f"{model.fwhm_metal:g}")
            self.sim_fwhm_c_var.set(f"{model.fwhm_compound:g}")
        self._draw_spectrum()

    def _load_user_spectrum(self):
        path = filedialog.askopenfilename(title="Load a 2-column spectrum (binding energy, intensity)",
                                          filetypes=[("Text / CSV", "*.txt *.csv *.dat *.asc"), ("All files", "*.*")])
        if not path:
            return
        try:
            x, y = sim.load_spectrum_file(path)
        except (OSError, ValueError) as e:
            messagebox.showerror("Could not load spectrum", str(e))
            return
        self.sim_data = (x, y, os.path.basename(path))
        self._draw_spectrum()

    def _clear_user_spectrum(self):
        self.sim_data = None
        self._draw_spectrum()

    def _draw_spectrum(self):
        if not HAVE_PLOT or self.last_calc is None:
            return
        model, problem = self._simulation_problem()
        if problem:
            self._draw_placeholder(problem)
            return
        lc = self.last_calc
        res = lc["res"]
        try:
            d_nm = float(self.sim_d_var.get())
            shift = float(self.sim_shift_var.get() or 0.0)
            fwhm_m = float(self.sim_fwhm_m_var.get()) if self.sim_fwhm_m_var.get().strip() else None
            fwhm_c = float(self.sim_fwhm_c_var.get()) if self.sim_fwhm_c_var.get().strip() else None
            if d_nm < 0 or (fwhm_m is not None and fwhm_m <= 0) or (fwhm_c is not None and fwhm_c <= 0):
                raise ValueError
        except ValueError:
            self.sim_status.set("Check the controls: thickness must be >= 0, FWHM > 0, BE shift a number.")
            return

        areas = [c["area"] for c in lc["comps"]]
        sub_frac, ov_fracs, _ = sim.oxide_area_fractions(
            d_nm, lc["theta"], lc["sub"]["N"], lc["sub"]["lam"], lc["sub"]["rsf"],
            res["N_ov_avg"], res["lam_ov_avg"], res["rsf_ov_avg"], areas)
        comps = [sim.SimComponent(lc["sub"]["label"], lc["sub"]["be"], sub_frac, "metal")]
        comps += [sim.SimComponent(c["label"], c["be"], f, "compound") for c, f in zip(lc["comps"], ov_fracs)]

        # What the simulation leaves out is stated on the figure; satellite regions are shaded and kept in the window.
        note, regions = sim.ignored_features(model, comps, shift)
        lo, hi = sim.default_window(model, comps, fwhm_m, fwhm_c, shift)
        for region in regions:
            lo, hi = min(lo, region.be_lo - 1.0), max(hi, region.be_hi + 1.0)
        simulation = sim.simulate(model, comps, fwhm_metal=fwhm_m, fwhm_compound=fwhm_c, be_shift=shift,
                                  window=(lo, hi))

        what_if = abs(d_nm - res["d_nm"]) > 5e-4
        status = ((f"WHAT-IF thickness {d_nm:.3f} nm (calculated: {res['d_nm']:.3f} nm). " if what_if
                   else f"Calculated thickness {d_nm:.3f} nm. ")
                  + ("" if self.compact else  # the legend already lists the areas; keep the line short on laptops
                     "Areas: " + ", ".join(f"{c.label} {100 * c.fraction:.0f}%" for c in comps) + ". ")
                  + f"Take-off angle {lc['theta']:g} deg. "
                  + ("Cu(II) satellites are ignored for simplicity (see the note under the plot)." if regions
                     else "Satellites/loss features are not drawn (see the note under the plot)."))

        fig = self.fig
        fig.clear()
        w_px, h_px = [max(float(v), 120.0) for v in (fig.get_size_inches() * fig.dpi)]   # only guards divide-by-zero
        footer = textwrap.fill("Note: " + note, int((w_px - 24) / 5.3))
        n_footer = footer.count("\n") + 1
        bottom_px = (44 + 11.5 * n_footer) if h_px < 330 else (56 + 12.5 * n_footer)   # tick labels + x label + note
        top_px = 28 if h_px < 330 else 30                                              # the ESTIMATE title
        gs = fig.add_gridspec(1, 2, width_ratios=[4.4, 1.0], wspace=0.30, left=58 / w_px, right=1 - 14 / w_px,
                              top=1 - top_px / h_px, bottom=bottom_px / h_px)
        ax = fig.add_subplot(gs[0])
        axd = fig.add_subplot(gs[1])
        fig.suptitle(sim.DISCLAIMER, color="#b00020", fontsize=11, fontweight="bold", y=1 - 5 / h_px, va="top")

        x, bg = simulation.be, simulation.background
        for i, (label, kind, frac, y) in enumerate(simulation.components):
            col = self.METAL_COLOUR if i == 0 else self.OXIDE_COLOURS[(i - 1) % len(self.OXIDE_COLOURS)]
            ax.fill_between(x, bg, bg + y, color=col, alpha=0.35, linewidth=0)
            ax.plot(x, bg + y, color=col, lw=1.0, label=f"{label} ({100 * frac:.0f}%)")
        ax.plot(x, bg, "--", color="0.45", lw=1.0, label="Shirley-type background")
        ax.plot(x, simulation.total, color="k", lw=1.6, label="Simulated total")
        for region in regions:
            ax.axvspan(region.be_lo, region.be_hi, facecolor="none", edgecolor="#9a9a9a", hatch="///",
                       linewidth=0.0, alpha=0.55, label=region.label)

        if self.sim_data is not None:
            dx, dy, name = self.sim_data
            m = (dx >= x.min()) & (dx <= x.max())
            if m.sum() >= 3:
                lo_y, hi_y = float(dy[m].min()), float(dy[m].max())
                if hi_y > lo_y:
                    span = simulation.total.max() - simulation.total.min()
                    scaled = (dy[m] - lo_y) / (hi_y - lo_y) * span + simulation.total.min()
                    ax.plot(dx[m], scaled, ".", ms=3.5, color="#555555", label=f"{name} (scaled to simulation)")
            else:
                status += f" Loaded data '{name}' does not overlap this BE window."
        self.sim_status.set(status)

        ax.invert_xaxis()
        ax.set_xlabel("Binding energy (eV)")
        ax.set_ylabel("Relative intensity (a.u.)")
        ax.set_yticks([])
        ax.text(0.5, 0.5, "ESTIMATE", transform=ax.transAxes, fontsize=64, color="red", alpha=0.07,
                ha="center", va="center", rotation=25, fontweight="bold")
        ax.legend(loc="upper left", fontsize=8 if h_px >= 330 else 7, frameon=False, ncol=1 if h_px >= 330 else 2,
                  title=f"{model.shell} region (areas normalised to 100%)", title_fontsize=8, alignment="left")
        fig.text(10 / w_px, 5 / h_px, footer, fontsize=7.5, style="italic", color="#444444", ha="left", va="bottom")

        # --- to-scale layer diagram -------------------------------------------------
        sub_h = 3.0 * lc["sub"]["lam"]
        top = max(0.3 * (d_nm + sub_h), 0.5) * (1.0 + 0.5 * (len(lc["comps"]) - 1))   # room for one label line per component
        axd.add_patch(Rectangle((0, -sub_h), 1, sub_h, facecolor=self.METAL_COLOUR, alpha=0.45, edgecolor="k"))
        axd.add_patch(Rectangle((0, 0), 1, max(d_nm, 1e-9), facecolor=self.OXIDE_COLOURS[0], alpha=0.6, edgecolor="k"))
        axd.set_xlim(0, 1)
        axd.set_ylim(-sub_h, d_nm + top)
        axd.set_xticks([])
        axd.set_ylabel("height (nm)")
        axd.text(0.5, -sub_h / 2, f"{lc['sub']['label']}\n(shown to 3λ = {sub_h:.1f} nm)",
                 ha="center", va="center", fontsize=8)
        names = "\n".join(c["label"] for c in lc["comps"])      # one per line: the diagram axes are narrow
        axd.text(0.5, d_nm + top / 2, f"{names}\nd = {d_nm:.2f} nm" + (" (what-if)" if what_if else ""),
                 ha="center", va="center", fontsize=7.5, fontweight="bold")
        axd.set_title("Layer (to scale)", fontsize=9)
        for side in ("top", "right"):
            axd.spines[side].set_visible(False)
        self.canvas.draw_idle()

    # ------------------------------------------------------------------
    def _build_help_tab(self):
        frame = ttk.Frame(self.help_tab)
        frame.pack(fill="both", expand=True, padx=6, pady=6)
        frame.columnconfigure(0, weight=1)
        frame.rowconfigure(0, weight=1)
        text = tk.Text(frame, wrap="word", font=("Segoe UI", 10), padx=14, pady=12, relief="flat", bg=CARD,
                       highlightthickness=1, highlightbackground="#d1d5db")
        scroll = ttk.Scrollbar(frame, orient="vertical", command=text.yview)
        text.configure(yscrollcommand=scroll.set)
        text.grid(row=0, column=0, sticky="nsew")
        scroll.grid(row=0, column=1, sticky="ns")
        text.insert("1.0", HELP_TEXT)
        text.configure(state="disabled")


HELP_TEXT = """XPS OVERLAYER THICKNESS CALCULATOR - METHOD NOTES

Source method
--------------
B. R. Strohmeier, Surface and Interface Analysis, 1990, 15, 51.
HarwellXPS Technical Note #1, "Oxide Thickness Determination Using the
Strohmeier Equation":
https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf

Model
-----
A uniform overlayer of oxide/hydroxide/nitride sits on top of a
semi-infinite substrate (e.g. the parent metal). Curve-fitting the
high-resolution core-level spectrum gives the percentage peak areas of the
substrate component (I_sub, e.g. metallic Cu0) and of the overlayer
component(s) (I_ov, e.g. Cu2+ in CuO). The method needs only ONE spectrum
per sample and no depth profiling, but is only valid for thin films -
approximately up to 3x the inelastic mean free path (IMFP) of the relevant
photoelectrons in the overlayer.

Single-component equation (Strohmeier, 1990):

    d = lambda_ox * sin(theta) * ln[ (N_m * lambda_m * I_ox)
                                       / (N_ox * lambda_ox * I_m) + 1 ]

Mixed overlayer (e.g. oxide + hydroxide present together), per tech note:

    d = lambda_ave * sin(theta) * ln[ (N_m * lambda_m * (I_oxide + I_hydroxide))
                                        / (N_ave * lambda_ave * I_m) + 1 ]

    lambda_ave and N_ave are the intensity-weighted averages of the
    individual overlayer components' IMFP and volume density. This tool
    computes those automatically from the ticked/active overlayer rows.

Generalised form used in this tool (adds optional RSF terms):

    d = lambda_ov * sin(theta) * ln[ (RSF_sub * N_sub * lambda_sub * I_ov)
                                       / (RSF_ov * N_ov * lambda_ov * I_sub) + 1 ]

    For the classic same-element case (e.g. Cu0 vs Cu2+, both from Cu 2p3/2)
    the RSF is identical for substrate and overlayer and cancels - leave
    both RSF fields at 1.0. Only change the RSF fields if you are combining
    intensities measured from genuinely different core levels/elements
    (e.g. normalising an N 1s overlayer signal against a Si 2p substrate
    signal), in which case enter each region's standard relative
    sensitivity factor (instrument- or library-specific, e.g. Scofield/
    Wagner values as tabulated by your XPS software).

Volume density N
-----------------
N = (atoms of the substrate element per formula unit) * density / molar mass

e.g. for Cu2O: N = (2 x 6.00) / 143.09. This matches the tech note's
worked calculation for the Cu/CuO/Cu2O/Cu(OH)2 system (only the
ratio Nm/Nox enters the equation, so absolute units are not important as
long as they are used consistently, which this tool does automatically).

Take-off angle
--------------
theta is measured from the SAMPLE SURFACE (grazing = small theta, normal
emission = 90 degrees), matching the Strohmeier/tech-note convention.
Some instruments/software quote the angle from the surface NORMAL instead
- if so, use theta_this_tool = 90 - theta_your_software.

On the pre-filled reference values
-----------------------------------
Several presets come straight from the original papers and are pinned by
tests/test_regression.py:
  - Al / Al2O3 (Al Ka and Mg Ka): Strohmeier, Surf. Interface Anal. 1990,
    15, 51. IMFP 2.8 nm oxide / 2.6 nm metal (Al Ka) and 2.4 / 2.2 nm
    (Mg Ka), Nm/Nox = 1.5. Also used in M. C. Biesinger's Al oxide/hydroxide
    spreadsheet. The 3.40 g/cm3 Al2O3 density is back-calculated to give
    Nm/Nox = 1.5 (Strohmeier quotes ~1.3-1.5 from gamma-Al2O3 densities).
  - Cu / Cu2O / CuO / Cu(OH)2: Biesinger, Surf. Interface Anal. 2017, 49,
    1325 (densities 8.96, 6.00, 6.31, 3.37 g/cm3; IMFPs 1.54, 1.77, 1.76,
    2.19 nm). These IMFPs are for Cu LMM Auger electrons (KE ~918 eV), so
    use the '(LMM)' presets with areas from the Cu LMM fit. Separate
    '(2p3/2)' presets (IMFP 1.01 / 1.23 / 1.22 / 1.51 nm for Cu / Cu2O / CuO /
    Cu(OH)2, KE ~554 eV with Al Ka) are COMPUTED here with TPP-2M, not
    published values. Do not mix the two sets: thickness scales almost
    linearly with the IMFPs, so the wrong set changes the answer by ~30%.
  - W / WO3: Carlson & McGuire, J. Electron Spectrosc. 1972/73, 1, 161
    (19.3 and 6.8 g/cm3, IMFP 1.28 and 2.63 nm, 45 deg take-off).

Si, SiO2, Si3N4, Ti and TiO2 use densities and band gaps from the Thermo
Avantage Densities.txt, with IMFPs COMPUTED here by TPP-2M (not published
values). SiON is a rough placeholder.

Also included: Mg/MgO, Ge/GeO2, Zn/ZnO (use the Zn LMM line), Zr/ZrO2, Hf/HfO2,
Sn/SnO2, SiC (with the SiO2 overlayer), GaAs and GaP (Ga 3d, with Ga2O3) and
CdTe (Te 3d5/2, with TeO2), InP (In 3d5/2, with In2O3) and PbS (S 2p, with
PbSO4). Densities and band gaps come from the Avantage file
where it has them, otherwise from standard handbook values (not from your files);
IMFPs are COMPUTED with TPP-2M. Check Nm/Nox against the Pilling-Bedworth ratio
if in doubt, and edit anything that doesn't match your film. All presets are
ALL EDITABLE. In particular:
  - IMFP is the most sensitive/uncertain input and depends on the exact
    kinetic energy (i.e. your X-ray source) and the true film density.
    For rigorous work, obtain IMFPs from the NIST Electron Inelastic-Mean-
    Free-Path Database (SRD 82) or a TPP-2M calculator for your specific
    photon energy and material, rather than relying on the defaults here.
  - Film density can differ substantially from bulk/crystalline values
    for thin, amorphous, native or anodic films (e.g. native Al2O3 vs.
    sapphire; thermal SiO2 vs. quartz) - adjust accordingly.
  - Binding energy fields are for peak identification/reference only and
    do not enter the thickness calculation directly.

Simulated spectrum pane (estimate)
-----------------------------------
Below the inputs, the "Simulated spectrum" pane draws what the calculated
thickness would look like for the substrate and oxide core-level peaks
(spin-orbit doublets, Shirley-type background) next to a to-scale layer
diagram, so you can compare it by eye with your own data. Use "Load my
spectrum..." (2-column text/CSV: binding energy, intensity) to overlay it; it
is scaled to the height of the simulation. Use the thickness box/slider for
what-if values, the BE shift for charge correction, and the FWHM boxes to
match your resolution. Drag the divider above the pane to trade height
between the inputs and the plot.

It is an ESTIMATE, not a fit: splitting, branching ratio and FWHM are typical
literature defaults, and intensities are relative (areas normalised to 100%).
For simplicity the simulation does NOT draw shake-up/shake-off satellites
(e.g. the Cu(II) satellites of CuO and Cu(OH)2 in Cu 2p3/2, which carry roughly
35-39% of the Cu(II) area), plasmon/loss features, Auger lines or asymmetric
metal tails. A note under the plot says so, and for Cu(II) the missing satellite
region is shaded. Auger presets (Cu LMM, Zn LMM) are not simulated. The pane
needs numpy and matplotlib (pip install numpy matplotlib); the calculator
itself does not.

Validity / limitations
------------------------
- Assumes a uniform, laterally-continuous overlayer of constant
  composition (no depth-dependent stoichiometry).
- Reliable up to roughly 3x the overlayer IMFP; beyond that the method
  becomes insensitive to further thickness increases and should not be
  trusted (the tool flags this automatically).
- Assumes photoelectron energies for substrate and overlayer peaks are
  similar enough that a single effective attenuation length is
  reasonable per component (true for chemically-shifted peaks of the
  same core level; use the RSF-generalised form with care otherwise).

References
----------
[1] B. R. Strohmeier, Surf. Interface Anal., 1990, 15, 51.
[2] M. C. Biesinger, Surf. Interface Anal. (in press), DOI: 10.1002/sia.6239.
[3] M. C. Biesinger, B. P. Payne, L. W. M. Lau, A. Gerson, R. S. C. Smart,
    Surf. Interface Anal., 2009, 41, 324.
"""


if __name__ == "__main__":
    app = App()
    app.mainloop()
