"""
Regression tests for the Strohmeier calculation.

Run from the repo root (stdlib only, no display needed):

    python -m unittest discover -s tests -v

Reference: M. C. Biesinger's "Aluminum Oxide/Hydroxide Thickness Calculator"
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
