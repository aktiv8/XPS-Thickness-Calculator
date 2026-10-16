"""
Simulated core-level XPS spectrum for a substrate + overlayer, driven by a layer thickness.

THIS IS AN ESTIMATE, NOT A FIT. It turns the forward form of the Strohmeier equation into relative peak areas,
draws each component as a (spin-orbit split) pseudo-Voigt peak, adds a Shirley-type background, and returns
curves a GUI can plot next to the user's measured spectrum.

Known simplifications (also stated in the app and README):
  - splitting / branching ratio / FWHM values are typical literature defaults from general knowledge, not from
    the user's files; FWHM is editable in the app because it depends on the instrument and the chemistry
  - no shake-up satellites (e.g. CuO, Cu(OH)2), no asymmetric metal tails, no plasmon/loss features,
    no Auger lines, relative intensities only
  - the preset binding energy is taken as the lower-BE (j = l - 1/2) component of a doublet (2p3/2, 3d5/2, 4f7/2)

Needs numpy only (no Tk, no matplotlib), so it is unit-testable on its own.
"""
import math
import re
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

DISCLAIMER = "ESTIMATE - simulated from the calculated thickness, not a fit"

_FWHM_TO_SIGMA = 1.0 / (2.0 * math.sqrt(2.0 * math.log(2.0)))


@dataclass(frozen=True)
class CoreLevelModel:
    shell: str
    splitting_ev: Optional[float]      # None = single peak
    area_ratio: Optional[float]        # low-BE : high-BE area ratio (2 for p, 1.5 for d, 4/3 for f)
    fwhm_metal: float                  # eV, of one component of the doublet
    fwhm_compound: float
    lorentzian_mix: float = 0.3        # fraction of Lorentzian in the pseudo-Voigt


# Typical literature values (general knowledge, not from the user's files). Edit if you have better ones.
CORE_LEVEL_MODEL = {m.shell: m for m in (
    CoreLevelModel("Al 2p", 0.44, 2.0, 0.55, 1.35),
    CoreLevelModel("Mg 2p", 0.28, 2.0, 0.60, 1.30),
    CoreLevelModel("Si 2p", 0.61, 2.0, 0.50, 1.25),
    CoreLevelModel("Ge 3d", 0.58, 1.5, 0.50, 1.20),
    CoreLevelModel("Ga 3d", 0.46, 1.5, 0.60, 1.20),
    CoreLevelModel("In 3d", 7.55, 1.5, 0.80, 1.20),
    CoreLevelModel("Ti 2p", 5.90, 2.0, 0.90, 1.10),
    CoreLevelModel("W 4f", 2.15, 4.0 / 3.0, 0.50, 1.10),
    CoreLevelModel("Zr 3d", 2.40, 1.5, 0.70, 1.30),
    CoreLevelModel("Hf 4f", 1.69, 4.0 / 3.0, 0.60, 1.30),
    CoreLevelModel("Sn 3d", 8.41, 1.5, 0.90, 1.30),
    CoreLevelModel("Te 3d", 10.39, 1.5, 1.00, 1.40),
    CoreLevelModel("S 2p", 1.18, 2.0, 0.80, 1.30),
    # Cu 2p: the presets are Cu 2p3/2 only. Shake-up satellites of CuO / Cu(OH)2 are NOT modelled.
    CoreLevelModel("Cu 2p", None, None, 0.90, 2.00),
)}

_LEVEL_RE = re.compile(r"^\s*([A-Z][a-z]?)\s*(\d[spdf])")


def model_for(core_level: str) -> Optional[CoreLevelModel]:
    """'Si 2p3/2' -> Si 2p model; Auger labels ('Cu LMM'), blanks and unknown levels -> None."""
    m = _LEVEL_RE.match(core_level or "")
    if not m:
        return None
    return CORE_LEVEL_MODEL.get(f"{m.group(1)} {m.group(2)}")


# ---------------------------------------------------------------------------
# Physics: forward form of the Strohmeier equation
# ---------------------------------------------------------------------------

def oxide_area_fractions(d_nm: float, theta_deg: float, n_sub: float, lam_sub: float, rsf_sub: float,
                         n_ov: float, lam_ov: float, rsf_ov: float,
                         component_areas: Sequence[float]) -> Tuple[float, List[float], float]:
    """
    Invert d = lam_ov sin(theta) ln[(RSF_sub N_sub lam_sub I_ov)/(RSF_ov N_ov lam_ov I_sub) + 1]:

        I_ov / I_sub = (RSF_ov N_ov lam_ov)/(RSF_sub N_sub lam_sub) * (exp(d/(lam_ov sin theta)) - 1)

    Returns (substrate fraction, per-component overlayer fractions, I_ov/I_sub). Fractions sum to 1 (area basis).
    The overlayer share is split between components in proportion to `component_areas` (the user's entered areas).
    """
    if d_nm < 0:
        raise ValueError("Thickness must be >= 0.")
    total = float(sum(component_areas))
    if total <= 0:
        raise ValueError("Need at least one overlayer area > 0.")
    arg = d_nm / (lam_ov * math.sin(math.radians(theta_deg)))
    ratio = (rsf_ov * n_ov * lam_ov) / (rsf_sub * n_sub * lam_sub) * math.expm1(min(arg, 700.0))
    substrate = 1.0 / (1.0 + ratio)
    overlayer = ratio / (1.0 + ratio)
    return substrate, [overlayer * a / total for a in component_areas], ratio


# ---------------------------------------------------------------------------
# Line shapes and background
# ---------------------------------------------------------------------------

def pseudo_voigt(x: np.ndarray, centre: float, fwhm: float, lorentzian_mix: float = 0.3) -> np.ndarray:
    """Area-normalised sum of a Gaussian and a Lorentzian of the same FWHM (mix = Lorentzian fraction)."""
    if fwhm <= 0:
        raise ValueError("FWHM must be > 0.")
    sigma = fwhm * _FWHM_TO_SIGMA
    gamma = fwhm / 2.0
    dx = x - centre
    gauss = np.exp(-0.5 * (dx / sigma) ** 2) / (sigma * math.sqrt(2.0 * math.pi))
    lorentz = (gamma / math.pi) / (dx ** 2 + gamma ** 2)
    return (1.0 - lorentzian_mix) * gauss + lorentzian_mix * lorentz


def shirley_background(peak: np.ndarray, step: float, base: float = 0.0) -> np.ndarray:
    """
    Self-consistent Shirley background for a known peak shape on an ascending-BE grid: the background rises from
    `base` (low-BE end) by `step` across the peak in proportion to the peak area accumulated at lower BE.
    """
    cum = np.concatenate(([0.0], np.cumsum(0.5 * (peak[1:] + peak[:-1]))))
    total = cum[-1]
    frac = cum / total if total > 0 else np.zeros_like(cum)
    return base + step * frac


# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------

@dataclass
class SimComponent:
    label: str
    be_ev: float
    fraction: float          # share of the total peak area (0-1)
    kind: str                # "metal" or "compound"


@dataclass
class Simulation:
    be: np.ndarray                   # ascending BE grid (eV)
    components: List[Tuple[str, str, float, np.ndarray]]   # (label, kind, fraction, curve without background)
    peaks: np.ndarray                # sum of components, no background
    background: np.ndarray
    total: np.ndarray


def simulate(model: CoreLevelModel, components: Sequence[SimComponent], fwhm_metal: Optional[float] = None,
             fwhm_compound: Optional[float] = None, be_shift: float = 0.0, step_fraction: float = 0.15,
             base_fraction: float = 0.02, point_spacing: float = 0.02, margin: float = 4.0) -> Simulation:
    """Relative-intensity spectrum with the peak areas normalised to 1 in total; height is a.u."""
    if not components:
        raise ValueError("No components to simulate.")
    fw = {"metal": fwhm_metal or model.fwhm_metal, "compound": fwhm_compound or model.fwhm_compound}
    split = model.splitting_ev or 0.0
    centres = [c.be_ev + be_shift for c in components]
    lo = min(centres) - margin - 2.0 * max(fw.values())
    hi = max(centres) + split + margin + 2.0 * max(fw.values())
    n = int(round((hi - lo) / point_spacing)) + 1
    x = np.linspace(lo, hi, n)

    curves = []
    for comp, centre in zip(components, centres):
        w = fw[comp.kind]
        if model.splitting_ev is None:
            y = comp.fraction * pseudo_voigt(x, centre, w, model.lorentzian_mix)
        else:
            low = model.area_ratio / (1.0 + model.area_ratio)
            y = comp.fraction * (low * pseudo_voigt(x, centre, w, model.lorentzian_mix)
                                 + (1.0 - low) * pseudo_voigt(x, centre + model.splitting_ev, w,
                                                              model.lorentzian_mix))
        curves.append((comp.label, comp.kind, comp.fraction, y))

    peaks = np.sum([c[3] for c in curves], axis=0)
    height = float(peaks.max())
    background = shirley_background(peaks, step=step_fraction * height, base=base_fraction * height)
    return Simulation(be=x, components=curves, peaks=peaks, background=background, total=peaks + background)


# ---------------------------------------------------------------------------
# Loading the user's measured spectrum
# ---------------------------------------------------------------------------

_SPLIT_RE = re.compile(r"[,\t; ]+")


def load_spectrum_file(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read a 2-column text/CSV file (BE, intensity). Header/comment lines are skipped; returns ascending BE."""
    xs, ys = [], []
    with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
        for line in fh:
            parts = [p for p in _SPLIT_RE.split(line.strip()) if p]
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
            except ValueError:
                continue
    if len(xs) < 5:
        raise ValueError("Could not find at least 5 numeric 'BE, intensity' rows in that file.")
    x, y = np.array(xs), np.array(ys)
    order = np.argsort(x)
    return x[order], y[order]
