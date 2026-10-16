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

No third-party dependencies - pure standard library (tkinter).
"""

import csv
import math
import time
import tkinter as tk
from dataclasses import dataclass, field
from tkinter import ttk, messagebox, filedialog


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
    note: str = ""


# Substrate ("metal") reference compounds
SUBSTRATES = {
    "Cu (metal)": Compound(
        "Cu (metal)", density_g_cm3=8.96, molar_mass_g_mol=63.546,
        atoms_per_formula=1, core_level="Cu 2p3/2", be_ev=932.63, imfp_nm=1.54,
        note="Biesinger, Surf. Interface Anal. 2017, 49, 1325. IMFP 1.54 nm (NIST SRD 71, TPP-2M) is for Cu LMM Auger "
             "electrons, KE ~918 eV (Cu L3M4,5M4,5 KE 918.61 eV). If your areas come from Cu 2p3/2 (KE ~554 eV with "
             "Al Ka) the IMFPs are smaller - replace all Cu/oxide IMFPs consistently.",
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
        atoms_per_formula=1, core_level="Si 2p3/2", be_ev=99.3, imfp_nm=2.87,
        note="Al Ka source, KE ~1387 eV. IMFP is a TPP-2M-type estimate - verify.",
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
    "CuO": Compound(
        "CuO", density_g_cm3=6.31, molar_mass_g_mol=79.55,
        atoms_per_formula=1, core_level="Cu 2p3/2", be_ev=933.76, imfp_nm=1.76,
        note="Biesinger 2017: 6.31 g/cm3, IMFP 1.76 nm at Cu LMM KE ~918 eV (Cu LMM KE 917.57 eV); Nm/Nox = 1.78.",
    ),
    "Cu2O": Compound(
        "Cu2O", density_g_cm3=6.00, molar_mass_g_mol=143.09,
        atoms_per_formula=2, core_level="Cu 2p3/2", be_ev=932.18, imfp_nm=1.77,
        note="Biesinger 2017: 6.00 g/cm3, IMFP 1.77 nm at Cu LMM KE ~918 eV (Cu LMM KE 916.99 eV); Nm/Nox = 1.68 "
             "(2 Cu atoms per formula unit). Cu 2p3/2 BE overlaps Cu0, so use the Cu LMM Auger line to separate them.",
    ),
    "Cu(OH)2": Compound(
        "Cu(OH)2", density_g_cm3=3.37, molar_mass_g_mol=97.56,
        atoms_per_formula=1, core_level="Cu 2p3/2", be_ev=934.67, imfp_nm=2.19,
        note="Biesinger 2017: 3.37 g/cm3, IMFP 2.19 nm at Cu LMM KE ~918 eV (Cu LMM KE 916.25 eV); Nm/Nox = 4.08.",
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
        "SiO2", density_g_cm3=2.20, molar_mass_g_mol=60.08,
        atoms_per_formula=1, core_level="Si 2p3/2", be_ev=103.4, imfp_nm=3.68,
        note="Amorphous thermal-oxide density (2.2); quartz is 2.65 - edit if relevant.",
    ),
    "Si3N4": Compound(
        "Si3N4", density_g_cm3=3.10, molar_mass_g_mol=140.28,
        atoms_per_formula=3, core_level="Si 2p3/2", be_ev=101.8, imfp_nm=3.00,
        note="Amorphous LPCVD-type silicon nitride density (~2.9-3.1); edit for your film.",
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


def volume_density(density_g_cm3: float, molar_mass_g_mol: float, atoms_per_formula: float) -> float:
    """Nx per tech note: atoms_per_formula * density / molar_mass (mol cm^-3 equivalent)."""
    if molar_mass_g_mol <= 0:
        raise ValueError("Molar mass must be > 0")
    return atoms_per_formula * density_g_cm3 / molar_mass_g_mol


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

FIELD_WIDTH = 9


class CompoundRow:
    """One editable row (substrate or overlayer) bound to a Compound preset dropdown."""

    def __init__(self, parent, row_idx, choices, default_key, is_substrate=False):
        self.is_substrate = is_substrate
        self.active_var = tk.BooleanVar(value=is_substrate)  # substrate always active
        self.choice_var = tk.StringVar(value=default_key)

        c = 0
        if not is_substrate:
            self.active_chk = ttk.Checkbutton(parent, variable=self.active_var)
            self.active_chk.grid(row=row_idx, column=c, padx=2)
        c += 1

        self.combo = ttk.Combobox(parent, textvariable=self.choice_var, values=list(choices.keys()),
                                   width=22, state="readonly")
        self.combo.grid(row=row_idx, column=c, padx=2, sticky="w")
        c += 1

        self.density_var = tk.StringVar()
        self.molar_var = tk.StringVar()
        self.atoms_var = tk.StringVar()
        self.be_var = tk.StringVar()
        self.core_var = tk.StringVar()
        self.imfp_var = tk.StringVar()
        self.rsf_var = tk.StringVar(value="1.0")
        self.area_var = tk.StringVar()

        entries = [
            (self.density_var, FIELD_WIDTH),
            (self.molar_var, FIELD_WIDTH),
            (self.atoms_var, 6),
            (self.core_var, 10),
            (self.be_var, 8),
            (self.imfp_var, FIELD_WIDTH),
            (self.rsf_var, 6),
            (self.area_var, FIELD_WIDTH),
        ]
        self.entry_widgets = []
        for var, width in entries:
            w = ttk.Entry(parent, textvariable=var, width=width)
            w.grid(row=row_idx, column=c, padx=2)
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
        self.geometry("1180x760")
        self.minsize(1040, 640)

        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True, padx=8, pady=8)

        self.calc_tab = ttk.Frame(nb)
        self.help_tab = ttk.Frame(nb)
        nb.add(self.calc_tab, text="Calculator")
        nb.add(self.help_tab, text="Method / References")

        self._build_calc_tab()
        self._build_help_tab()

        self.history = []  # list of dict rows for CSV export

    # ------------------------------------------------------------------
    def _build_calc_tab(self):
        top = ttk.Frame(self.calc_tab)
        top.pack(fill="x", padx=6, pady=(6, 0))

        ttk.Label(
            top,
            text="Uniform-overlayer thickness from XPS peak areas (Strohmeier, 1990 / HarwellXPS Tech Note #1)",
            font=("Segoe UI", 11, "bold"),
        ).pack(anchor="w")
        ttk.Label(
            top,
            text="Fill in / edit the substrate row and up to 4 overlayer rows (tick 'active' to include a "
                 "component in a mixed layer, e.g. oxide + hydroxide). Densities, molar masses, binding "
                 "energies and IMFPs are pre-filled reference defaults - verify/edit for your material and instrument.",
            wraplength=1100, justify="left",
        ).pack(anchor="w", pady=(2, 8))

        angle_frame = ttk.Frame(self.calc_tab)
        angle_frame.pack(fill="x", padx=6)
        ttk.Label(angle_frame, text="Photoelectron take-off angle theta (deg from surface, 90 = normal emission):").pack(side="left")
        self.theta_var = tk.StringVar(value="90")
        ttk.Entry(angle_frame, textvariable=self.theta_var, width=8).pack(side="left", padx=6)

        table = ttk.Frame(self.calc_tab)
        table.pack(fill="x", padx=6, pady=8)

        headers_sub = ["", "Compound", "Density\n(g/cm3)", "Molar mass\n(g/mol)", "Atoms/\nformula",
                       "Core level", "BE\n(eV)", "IMFP\n(nm)", "RSF", "Peak area\n(I)"]
        # substrate has no "active" checkbox column content but keep alignment
        for i, h in enumerate(headers_sub):
            ttk.Label(table, text=h, font=("Segoe UI", 8, "bold"), justify="center").grid(row=0, column=i, padx=2)

        ttk.Label(table, text="Substrate", font=("Segoe UI", 9, "bold")).grid(row=1, column=0, sticky="w")
        self.substrate_row = CompoundRow(table, 2, SUBSTRATES, "Cu (metal)", is_substrate=True)

        ttk.Separator(table, orient="horizontal").grid(row=3, column=0, columnspan=10, sticky="ew", pady=6)
        ttk.Label(table, text="Overlayer component(s)", font=("Segoe UI", 9, "bold")).grid(row=4, column=0, columnspan=3, sticky="w")

        self.overlayer_rows = []
        defaults = ["CuO", "Cu2O", "Cu(OH)2", "None"]
        for i in range(N_OVERLAYER_ROWS):
            row = CompoundRow(table, 5 + i, OVERLAYERS, defaults[i] if i < len(defaults) else "None",
                               is_substrate=False)
            self.overlayer_rows.append(row)

        btn_frame = ttk.Frame(self.calc_tab)
        btn_frame.pack(fill="x", padx=6, pady=(4, 4))
        ttk.Button(btn_frame, text="Calculate thickness", command=self._on_calculate).pack(side="left")
        ttk.Button(btn_frame, text="Export history to CSV...", command=self._export_csv).pack(side="left", padx=8)
        ttk.Button(btn_frame, text="Clear peak areas", command=self._clear_areas).pack(side="left")

        result_frame = ttk.LabelFrame(self.calc_tab, text="Result")
        result_frame.pack(fill="both", expand=True, padx=6, pady=6)
        self.result_text = tk.Text(result_frame, height=10, wrap="word", state="disabled",
                                    font=("Consolas", 10))
        self.result_text.pack(fill="both", expand=True, padx=4, pady=4)

    def _clear_areas(self):
        self.substrate_row.area_var.set("")
        for row in self.overlayer_rows:
            row.area_var.set("")

    def _on_calculate(self):
        try:
            theta = float(self.theta_var.get())
            sub_label, N_sub, lam_sub, rsf_sub, I_sub = self.substrate_row.read()

            components = []
            for row in self.overlayer_rows:
                if not row.active_var.get() or row.choice_var.get() == "None":
                    continue
                label, N, lam, rsf, area = row.read()
                if area <= 0:
                    continue
                components.append(OverlayerInput(label, N, lam, rsf, area))

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
        lines = []
        lines.append(f"Substrate: {sub_label}   N_sub = {N_sub:.5f}   lambda_sub = {lam_sub:g} nm   I_sub = {I_sub:g}")
        lines.append(f"Overlayer component(s): {comp_names}")
        lines.append(f"Take-off angle theta = {theta:g} deg")
        lines.append("-" * 70)
        lines.append(f"Intensity-weighted N_overlayer  = {res['N_ov_avg']:.5f}")
        lines.append(f"Intensity-weighted lambda_overlayer = {res['lam_ov_avg']:.4f} nm")
        lines.append(f"N_sub / N_overlayer ratio        = {res['N_ratio_sub_over_ov']:.4f}")
        lines.append("")
        lines.append(f"==> Overlayer thickness d = {res['d_nm']:.3f} nm  ({res['d_angstrom']:.2f} Angstrom)")
        lines.append("")
        limit_txt = "within" if res["within_validity"] else "EXCEEDS"
        lines.append(f"Validity check: method is reliable to ~3x lambda_overlayer = "
                     f"{res['validity_limit_nm']:.2f} nm. Calculated thickness is {limit_txt} this limit.")
        if not res["within_validity"]:
            lines.append("  -> Result may be unreliable; consider angle-resolved XPS, sputter-depth "
                         "profiling, or an alternative technique (ellipsometry, TEM) for thicker films.")

        self.result_text.configure(state="normal")
        self.result_text.delete("1.0", "end")
        self.result_text.insert("1.0", "\n".join(lines))
        self.result_text.configure(state="disabled")

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
    def _build_help_tab(self):
        text = tk.Text(self.help_tab, wrap="word", font=("Segoe UI", 10), padx=10, pady=10)
        text.pack(fill="both", expand=True)
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
    2.19 nm). These IMFPs are for Cu LMM Auger electrons (KE ~918 eV); if
    your areas come from Cu 2p3/2 (KE ~554 eV) use smaller 2p3/2 IMFPs.
  - W / WO3: Carlson & McGuire, J. Electron Spectrosc. 1972/73, 1, 161
    (19.3 and 6.8 g/cm3, IMFP 1.28 and 2.63 nm, 45 deg take-off).

For the other systems (Si/SiO2/Si3N4/SiON), densities, molar masses, typical
binding energies and IMFP values supplied as defaults are representative
literature / TPP-2M-style estimates. All presets are
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
