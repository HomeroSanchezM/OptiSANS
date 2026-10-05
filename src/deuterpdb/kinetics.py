"""Exchange probability models for labile hydrogens.

`simple` mode: p = d2o_percent / 100.
`kinetic` mode: two populations (accessible / buried), each with
p(t) = (d2o/100) * (1 - exp(-k*t)).
t = 0 gives 0 (no exchange); t -> inf gives d2o/100 (equilibrium at the D2O
concentration).
Modeled after DeuteRater's fraction new = a - a * exp(-k * t) with a = 1.
"""

from __future__ import annotations

import math


def simple_probability(d2o_percent: float) -> float:
    """Exchange probability in the simple mode."""
    return d2o_percent / 100.0


def kinetic_probability(d2o_percent: float, rate: float, time_minutes: float) -> float:
    """Fraction of a population exchanged after time_minutes at rate (1/min).

    The D2O level sets the equilibrium asymptote: p(t) = (d2o/100) * (1 - exp(-k*t)).
    t = 0 gives 0 (no exchange); t -> inf gives d2o/100 (equilibrium at the D2O
    concentration). A slower rate lags behind a faster one at any finite time.
    """
    if time_minutes <= 0.0:
        return 0.0
    return (d2o_percent / 100.0) * (1.0 - math.exp(-rate * time_minutes))
