"""
Tests for the simulated-spectrum module (numpy only, no display needed).

    python -m unittest discover -s tests -v
"""
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

try:
    import numpy as np
    import spectrum_sim as S
    HAVE_NUMPY = True
except ImportError:          # the calculator itself does not need numpy
    HAVE_NUMPY = False

import app as A  # noqa: E402


def _trapz(y, x):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


@unittest.skipUnless(HAVE_NUMPY, "numpy not installed")
class ForwardModel(unittest.TestCase):
    def _case(self, comps):
        """comps: [(label, N, lam, rsf, area)] with the Al-like substrate N=1.5, lam=2.6 nm."""
        res = A.strohmeier_thickness(
            90, 1.5, 2.6, 1.0, 35,
            [A.OverlayerInput(l, n, lam, r, a) for l, n, lam, r, a in comps])
        sub, fracs, _ = S.oxide_area_fractions(
            res["d_nm"], 90, 1.5, 2.6, 1.0, res["N_ov_avg"], res["lam_ov_avg"], res["rsf_ov_avg"],
            [c[4] for c in comps])
        return res, sub, fracs

    def test_round_trip_single_layer(self):
        _, sub, fracs = self._case([("Al2O3", 1.0, 2.8, 1.0, 65)])
        self.assertAlmostEqual(sub, 0.35, places=9)
        self.assertAlmostEqual(fracs[0], 0.65, places=9)

    def test_round_trip_mixed_layer(self):
        _, sub, fracs = self._case([("a", 1.0, 2.8, 1.0, 40), ("b", 0.4, 3.4, 1.0, 20)])
        self.assertAlmostEqual(sub + sum(fracs), 1.0, places=12)
        # entered metal fraction is 35 / (35 + 60)
        self.assertAlmostEqual(sub, 35 / 95.0, places=9)
        self.assertAlmostEqual(fracs[0] / fracs[1], 2.0, places=9)    # split follows the entered areas

    def test_round_trip_with_unequal_rsf(self):
        res = A.strohmeier_thickness(90, 1.5, 2.6, 2.0, 35, [A.OverlayerInput("x", 1.0, 2.8, 5.0, 65)])
        sub, fracs, _ = S.oxide_area_fractions(res["d_nm"], 90, 1.5, 2.6, 2.0, 1.0, 2.8, 5.0, [65])
        self.assertAlmostEqual(sub, 35 / 100.0, places=9)

    def test_round_trip_off_normal_angle(self):
        res = A.strohmeier_thickness(45, 1.5, 2.6, 1.0, 35, [A.OverlayerInput("x", 1.0, 2.8, 1.0, 65)])
        sub, _, _ = S.oxide_area_fractions(res["d_nm"], 45, 1.5, 2.6, 1.0, 1.0, 2.8, 1.0, [65])
        self.assertAlmostEqual(sub, 0.35, places=9)

    def test_limits_and_monotonic(self):
        kw = dict(theta_deg=90, n_sub=1.5, lam_sub=2.6, rsf_sub=1.0, n_ov=1.0, lam_ov=2.8, rsf_ov=1.0,
                  component_areas=[1.0])
        self.assertEqual(S.oxide_area_fractions(0.0, **kw)[1][0], 0.0)
        prev = -1.0
        for d in (0.5, 1, 2, 4, 8, 16, 100, 1e6):
            frac = S.oxide_area_fractions(d, **kw)[1][0]
            self.assertGreater(frac, prev)
            prev = frac
        self.assertAlmostEqual(prev, 1.0, places=9)

    def test_bad_inputs(self):
        with self.assertRaises(ValueError):
            S.oxide_area_fractions(-1, 90, 1.5, 2.6, 1.0, 1.0, 2.8, 1.0, [1.0])
        with self.assertRaises(ValueError):
            S.oxide_area_fractions(1, 90, 1.5, 2.6, 1.0, 1.0, 2.8, 1.0, [0.0])


@unittest.skipUnless(HAVE_NUMPY, "numpy not installed")
class LineShapes(unittest.TestCase):
    def test_area_is_one(self):
        x = np.linspace(-200, 200, 400001)
        for mix in (0.0, 0.3, 1.0):
            area = _trapz(S.pseudo_voigt(x, 0.0, 1.0, mix), x)
            self.assertAlmostEqual(area, 1.0, delta=0.01 if mix else 1e-6, msg=f"mix={mix}")

    def test_fwhm_of_pure_gaussian_and_lorentzian(self):
        x = np.linspace(-5, 5, 200001)
        for mix in (0.0, 1.0):
            y = S.pseudo_voigt(x, 0.0, 1.3, mix)
            above = x[y >= y.max() / 2.0]
            self.assertAlmostEqual(above.max() - above.min(), 1.3, delta=0.001, msg=f"mix={mix}")

    def test_shirley_rises_monotonically_from_base_by_step(self):
        x = np.linspace(70, 78, 801)
        peak = S.pseudo_voigt(x, 74.0, 1.0)
        bg = S.shirley_background(peak, step=0.3, base=0.05)
        self.assertAlmostEqual(bg[0], 0.05, places=12)
        self.assertAlmostEqual(bg[-1], 0.35, places=12)
        self.assertTrue(np.all(np.diff(bg) >= -1e-15))


@unittest.skipUnless(HAVE_NUMPY, "numpy not installed")
class Simulate(unittest.TestCase):
    def _al(self, fractions=(0.35, 0.65), **kw):
        m = S.model_for("Al 2p")
        comps = [S.SimComponent("Al", 72.9, fractions[0], "metal"), S.SimComponent("Al2O3", 74.4, fractions[1], "compound")]
        return S.simulate(m, comps, **kw)

    def test_component_areas_match_fractions(self):
        # A very wide window captures the Lorentzian tails, so areas must equal the fractions almost exactly...
        wide = self._al(margin=80.0, point_spacing=0.05)
        for (label, kind, frac, y) in wide.components:
            self.assertAlmostEqual(_trapz(y, wide.be), frac, delta=0.005, msg=label)
        # ...while the default display window clips a little tail (a couple of % at most), as real fit windows do.
        sim = self._al()
        for (label, kind, frac, y) in sim.components:
            self.assertAlmostEqual(_trapz(y, sim.be), frac, delta=0.03, msg=label)
            self.assertLess(_trapz(y, sim.be), frac + 1e-9, msg=label)

    def test_doublet_ratio_and_splitting(self):
        m = S.model_for("Al 2p")
        sim = S.simulate(m, [S.SimComponent("Al", 72.9, 1.0, "metal")], fwhm_metal=0.15)   # narrow: resolved
        y = sim.components[0][3]
        main = sim.be[np.argmax(y)]
        self.assertAlmostEqual(main, 72.9, delta=0.02)
        # the higher-BE partner is a local maximum 0.44 eV above, with half the area (so ~half the height)
        win = (sim.be > 73.2) & (sim.be < 73.6)
        partner = sim.be[win][np.argmax(y[win])]
        self.assertAlmostEqual(partner, 72.9 + 0.44, delta=0.02)
        self.assertAlmostEqual(y[win].max() / y.max(), 0.5, delta=0.03)

    def test_single_peak_level_has_no_partner(self):
        m = S.model_for("Cu 2p3/2")
        self.assertIsNone(m.splitting_ev)
        sim = S.simulate(m, [S.SimComponent("Cu", 932.63, 1.0, "metal")])
        self.assertAlmostEqual(sim.be[np.argmax(sim.components[0][3])], 932.63, delta=0.02)

    def test_be_shift_moves_everything(self):
        a, b = self._al(), self._al(be_shift=-1.5)
        self.assertAlmostEqual(a.be[np.argmax(a.peaks)] - 1.5, b.be[np.argmax(b.peaks)], delta=0.03)

    def test_background_and_total(self):
        sim = self._al()
        self.assertTrue(np.all(np.diff(sim.background) >= -1e-12))
        self.assertTrue(np.all(sim.total >= sim.background - 1e-12))
        np.testing.assert_allclose(sim.total, sim.peaks + sim.background)

    def test_no_components(self):
        with self.assertRaises(ValueError):
            S.simulate(S.model_for("Al 2p"), [])


@unittest.skipUnless(HAVE_NUMPY, "numpy not installed")
class ModelLookup(unittest.TestCase):
    def test_j_suffix_is_stripped(self):
        self.assertIs(S.model_for("Si 2p3/2"), S.CORE_LEVEL_MODEL["Si 2p"])
        self.assertIs(S.model_for("W 4f7/2"), S.CORE_LEVEL_MODEL["W 4f"])
        self.assertIs(S.model_for("Ge 3d"), S.CORE_LEVEL_MODEL["Ge 3d"])

    def test_auger_blank_and_unknown_are_unsupported(self):
        for s in ("Cu LMM", "Zn LMM", "", "Xx 9z", "nonsense"):
            self.assertIsNone(S.model_for(s), s)

    def test_every_preset_with_a_binding_energy_can_be_simulated(self):
        for presets in (A.SUBSTRATES, A.OVERLAYERS):
            for key, c in presets.items():
                if c is None or key.startswith("Custom"):
                    continue
                if c.be_ev:
                    self.assertIsNotNone(S.model_for(c.core_level), f"{key} ({c.core_level})")
                else:   # Auger presets carry no BE and must not claim a simulable core level
                    self.assertIsNone(S.model_for(c.core_level), f"{key} ({c.core_level})")

    def test_model_values_are_sane(self):
        for m in S.CORE_LEVEL_MODEL.values():
            self.assertGreater(m.fwhm_metal, 0)
            self.assertGreater(m.fwhm_compound, 0)
            self.assertTrue(0 <= m.lorentzian_mix <= 1)
            if m.splitting_ev is not None:
                self.assertGreater(m.splitting_ev, 0)
                self.assertGreater(m.area_ratio, 1.0)


@unittest.skipUnless(HAVE_NUMPY, "numpy not installed")
class LoadSpectrum(unittest.TestCase):
    def _load(self, text):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False, encoding="utf-8") as fh:
            fh.write(text)
            path = fh.name
        try:
            return S.load_spectrum_file(path)
        finally:
            os.unlink(path)

    def test_delimiters_headers_and_descending_order(self):
        rows = [(78.0, 10), (77.0, 20), (76.0, 30), (75.0, 40), (74.0, 50), (73.0, 60)]
        for delim in (",", "\t", " ", ";", "  "):
            text = "# my spectrum\nBE,counts\n" + "\n".join(f"{x}{delim}{y}" for x, y in rows) + "\n"
            x, y = self._load(text)
            self.assertEqual(list(x), sorted(r[0] for r in rows), repr(delim))
            self.assertEqual(y[0], 60)      # sorted ascending, so 73.0 first

    def test_extra_columns_are_ignored(self):
        x, y = self._load("\n".join(f"{70 + i} {i} 99 99" for i in range(8)))
        self.assertEqual(len(x), 8)
        self.assertEqual(list(y), list(range(8)))

    def test_too_few_rows(self):
        with self.assertRaises(ValueError):
            self._load("1 2\n3 4\n")


if __name__ == "__main__":
    unittest.main()
