"""
Regression tests for the Strohmeier calculation.

Run from the repo root (stdlib only, no display needed):

    python -m unittest discover -s tests -v

References pinned here (all published numbers, not app output):
  - Strohmeier, Surf. Interface Anal. 1990, 15, 51 (Al/Al2O3, Al Ka and Mg Ka)
  - Carlson & McGuire, J. Electron Spectrosc. 1972/73, 1, 161 (W/WO3, 45 deg)
  - Biesinger, Surf. Interface Anal. 2017, 49, 1325 (Cu oxides/hydroxide)
  - Biesinger's "Aluminum Oxide/Hydroxide Thickness Calculator"
spreadsheet (Surface Science Western, 2008), which uses Strohmeier's original
Al values: 90 deg take-off, IMFP_ox = 28 A, IMFP_m = 26 A, Nm/Nox = 1.5.
The sheet's three worked examples (Io, Im as percentage areas -> d in A) are
reproduced below. The sheet rounds d to 0.1 A, hence the tolerance.
"""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import app as A  # noqa: E402  (importing app does not open a window)

# (Io %, Im %, d in Angstrom from the spreadsheet)
BIESINGER_AL_ROWS = [
    (65, 35, 35.8),
    (23, 77, 9.7),
    (85, 15, 61.2),
]
TOL_ANGSTROM = 0.1


def _thickness(theta, n_sub, lam_sub, i_sub, n_ov, lam_ov, i_ov, rsf_sub=1.0, rsf_ov=1.0):
    return A.strohmeier_thickness(
        theta, n_sub, lam_sub, rsf_sub, i_sub,
        [A.OverlayerInput("ov", n_ov, lam_ov, rsf_ov, i_ov)],
    )


class BiesingerAlSpreadsheet(unittest.TestCase):
    def test_explicit_inputs(self):
        # Nm = 1.5, Nox = 1.0, lambda in nm (2.6 / 2.8), theta = 90
        for io, im, ref in BIESINGER_AL_ROWS:
            d = _thickness(90, 1.5, 2.6, im, 1.0, 2.8, io)["d_angstrom"]
            self.assertAlmostEqual(d, ref, delta=TOL_ANGSTROM, msg=f"Io={io} Im={im}")

    def test_built_in_al_presets(self):
        sub, ox = A.SUBSTRATES["Al (metal)"], A.OVERLAYERS["Al2O3"]
        n_sub = A.volume_density(sub.density_g_cm3, sub.molar_mass_g_mol, sub.atoms_per_formula)
        n_ox = A.volume_density(ox.density_g_cm3, ox.molar_mass_g_mol, ox.atoms_per_formula)
        self.assertAlmostEqual(n_sub / n_ox, 1.5, delta=0.005)
        self.assertEqual((sub.imfp_nm, ox.imfp_nm), (2.6, 2.8))
        for io, im, ref in BIESINGER_AL_ROWS:
            d = _thickness(90, n_sub, sub.imfp_nm, im, n_ox, ox.imfp_nm, io)["d_angstrom"]
            self.assertAlmostEqual(d, ref, delta=TOL_ANGSTROM, msg=f"Io={io} Im={im}")


def _preset_N(sub_key, ov_key):
    sub, ov = A.SUBSTRATES[sub_key], A.OVERLAYERS[ov_key]
    n_sub = A.volume_density(sub.density_g_cm3, sub.molar_mass_g_mol, sub.atoms_per_formula)
    n_ov = A.volume_density(ov.density_g_cm3, ov.molar_mass_g_mol, ov.atoms_per_formula)
    return sub, ov, n_sub, n_ov


class Strohmeier1990(unittest.TestCase):
    """Numbers quoted in the paper's text. It rounds Nm*lam_m/(No*lam_o) to 1.4, so allow ~1 A."""

    def test_al_ka_presets_ratio_3p30_gives_48A(self):   # 99.99% Al sheet, Al Ka: Io/Im = 3.30 -> 48 A
        sub, ov, n_sub, n_ov = _preset_N("Al (metal)", "Al2O3")
        d = _thickness(90, n_sub, sub.imfp_nm, 100, n_ov, ov.imfp_nm, 330)["d_angstrom"]
        self.assertAlmostEqual(d, 48, delta=1.0)

    def test_mg_ka_presets_ratio_4p97_gives_50A(self):   # same sheet, Mg Ka: Io/Im = 4.97 -> 50 A
        sub, ov, n_sub, n_ov = _preset_N("Al (metal, Mg Ka)", "Al2O3 (Mg Ka)")
        self.assertEqual((sub.imfp_nm, ov.imfp_nm), (2.2, 2.4))
        d = _thickness(90, n_sub, sub.imfp_nm, 100, n_ov, ov.imfp_nm, 497)["d_angstrom"]
        self.assertAlmostEqual(d, 50, delta=1.0)

    def test_max_measurable_ratio_18_gives_78A_mg_and_91A_al(self):
        mg = _preset_N("Al (metal, Mg Ka)", "Al2O3 (Mg Ka)")
        al = _preset_N("Al (metal)", "Al2O3")
        d_mg = _thickness(90, mg[2], mg[0].imfp_nm, 100, mg[3], mg[1].imfp_nm, 1800)["d_angstrom"]
        d_al = _thickness(90, al[2], al[0].imfp_nm, 100, al[3], al[1].imfp_nm, 1800)["d_angstrom"]
        self.assertAlmostEqual(d_mg, 78, delta=1.0)
        self.assertAlmostEqual(d_al, 91, delta=1.0)


class CarlsonMcGuire1972Tungsten(unittest.TestCase):
    def test_w_presets_reproduce_paper_constants(self):
        sub, ov, n_sub, n_ov = _preset_N("W (metal)", "WO3")
        self.assertAlmostEqual(n_sub / n_ov, 3.58, delta=0.01)                       # paper: D2/D1 = 3.58
        k = n_sub / n_ov * sub.imfp_nm / ov.imfp_nm
        self.assertAlmostEqual(k, 1.74, delta=0.01)                                 # paper: D2*s1/(D1*s2) = 1.74

    def test_equation_5_round_trip_at_45_degrees(self):
        # Paper eq. (5): Im/Iox = K * exp(-x/lam) / (1 - exp(-x/lam)), x = sqrt(2) * d, lam = 26.3 A, K = 1.74.
        # Feed that ratio to the app (theta = 45 deg) and it must return d. Thicknesses from the paper's Table 1.
        sub, ov, n_sub, n_ov = _preset_N("W (metal)", "WO3")
        for d_true in (4.2, 14.9, 30.6, 63.3, 125.0):
            s = math.sqrt(2) * d_true / 26.3
            r_m_over_ox = 1.74 * math.exp(-s) / (1 - math.exp(-s))
            d = _thickness(45, n_sub, sub.imfp_nm, r_m_over_ox, n_ov, ov.imfp_nm, 1.0)["d_angstrom"]
            self.assertAlmostEqual(d, d_true, delta=0.1, msg=f"d={d_true}")


class Biesinger2017Copper(unittest.TestCase):
    def test_density_ratios_match_paper(self):
        s = A.SUBSTRATES["Cu (metal, LMM)"]
        n_m = A.volume_density(s.density_g_cm3, s.molar_mass_g_mol, 1)
        for key, ref in (("Cu2O", 1.68), ("CuO", 1.78), ("Cu(OH)2", 4.08)):
            o = A.OVERLAYERS[key + " (LMM)"]
            n_o = A.volume_density(o.density_g_cm3, o.molar_mass_g_mol, o.atoms_per_formula)
            self.assertAlmostEqual(n_m / n_o, ref, delta=0.01, msg=key)

    def test_imfps_match_paper(self):
        self.assertEqual(A.SUBSTRATES["Cu (metal, LMM)"].imfp_nm, 1.54)
        for key, ref in (("Cu2O", 1.77), ("CuO", 1.76), ("Cu(OH)2", 2.19)):
            self.assertEqual(A.OVERLAYERS[key + " (LMM)"].imfp_nm, ref, key)

    def test_worked_example_63_Cu2O_13_CuOH2_25_metal_is_about_3p2_nm(self):
        sub = A.SUBSTRATES["Cu (metal, LMM)"]
        n_m = A.volume_density(sub.density_g_cm3, sub.molar_mass_g_mol, 1)
        comps = []
        for key, area in (("Cu2O", 63), ("Cu(OH)2", 13)):
            o = A.OVERLAYERS[key + " (LMM)"]
            comps.append(A.OverlayerInput(key, A.volume_density(o.density_g_cm3, o.molar_mass_g_mol,
                                                                  o.atoms_per_formula), o.imfp_nm, 1.0, area))
        d = A.strohmeier_thickness(90, n_m, sub.imfp_nm, 1.0, 25, comps)["d_nm"]
        self.assertAlmostEqual(d, 3.2, delta=0.05)   # paper: "~3.2 nm"


class CuTpp2mImfps(unittest.TestCase):
    """The Cu 2p3/2 presets are computed (not published) with TPP-2M; Biesinger's 918 eV values validate the code."""
    # name: (density, molar mass, Nv, band gap eV, Biesinger 918 eV IMFP nm or None)
    CU = {"Cu2O": (6.00, 143.09, 28, 2.1, 1.77), "CuO": (6.31, 79.55, 17, 1.25, 1.76),
          "Cu(OH)2": (3.37, 97.56, 25, 1.85, 2.19)}

    def test_tpp2m_reproduces_biesinger_oxide_imfps_at_918_eV(self):
        for name, (rho, m, nv, eg, ref) in self.CU.items():
            self.assertAlmostEqual(A.tpp2m_imfp_nm(918, rho, m, nv, eg), ref, delta=0.01, msg=name)

    def test_2p32_presets_match_tpp2m_at_554_eV(self):
        for name, (rho, m, nv, eg, _) in self.CU.items():
            self.assertAlmostEqual(A.OVERLAYERS[name + " (2p3/2)"].imfp_nm,
                                   A.tpp2m_imfp_nm(554, rho, m, nv, eg), delta=0.005, msg=name)
        self.assertAlmostEqual(A.SUBSTRATES["Cu (metal, 2p3/2)"].imfp_nm,
                               A.tpp2m_imfp_nm(554, 8.96, 63.546, 11), delta=0.005)

    def test_2p32_and_lmm_sets_share_densities(self):
        for name in ("Cu2O", "CuO", "Cu(OH)2"):
            a, b = A.OVERLAYERS[name + " (LMM)"], A.OVERLAYERS[name + " (2p3/2)"]
            self.assertEqual((a.density_g_cm3, a.molar_mass_g_mol, a.atoms_per_formula),
                             (b.density_g_cm3, b.molar_mass_g_mol, b.atoms_per_formula))


class EquationProperties(unittest.TestCase):
    def test_mixed_layer_reduces_to_single_when_components_identical(self):
        one = _thickness(90, 1.5, 2.6, 35, 1.0, 2.8, 65)["d_nm"]
        two = A.strohmeier_thickness(
            90, 1.5, 2.6, 1.0, 35,
            [A.OverlayerInput("a", 1.0, 2.8, 1.0, 40), A.OverlayerInput("b", 1.0, 2.8, 1.0, 25)],
        )["d_nm"]
        self.assertAlmostEqual(one, two, places=12)

    def test_takeoff_angle_scales_with_sin_theta(self):
        d90 = _thickness(90, 1.5, 2.6, 35, 1.0, 2.8, 65)
        # Same areas at 30 deg: only the sin(theta) prefactor changes for the same ratio term
        d30 = _thickness(30, 1.5, 2.6, 35, 1.0, 2.8, 65)
        self.assertAlmostEqual(d30["d_nm"], d90["d_nm"] * math.sin(math.radians(30)), places=12)

    def test_rsf_cancels_when_equal(self):
        base = _thickness(90, 1.5, 2.6, 35, 1.0, 2.8, 65)["d_nm"]
        scaled = _thickness(90, 1.5, 2.6, 35, 1.0, 2.8, 65, rsf_sub=4.2, rsf_ov=4.2)["d_nm"]
        self.assertAlmostEqual(base, scaled, places=12)

    def test_volume_density_counts_atoms_per_formula(self):
        # Cu2O: [(6.14 x 2) / 143.09] per the HarwellXPS tech note
        self.assertAlmostEqual(A.volume_density(6.14, 143.09, 2), 0.0858, places=3)

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            _thickness(90, 1.5, 2.6, 0, 1.0, 2.8, 65)       # no substrate signal
        with self.assertRaises(ValueError):
            _thickness(0, 1.5, 2.6, 35, 1.0, 2.8, 65)        # angle out of range
        with self.assertRaises(ValueError):
            A.strohmeier_thickness(90, 1.5, 2.6, 1.0, 35, [])  # no overlayer


if __name__ == "__main__":
    unittest.main()
