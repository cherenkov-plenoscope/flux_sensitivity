import numpy as np
from . import integral_sensitivity
from . import critical_rate


def _find_intersection_two_lines(b1, m1, b2, m2):
    """
    Find the intersection of two affine functions:
    f1(x) = m1 * x + b1, and f2(x) = m2 * x + b2 in x.

    Parameters
    ----------
    m1 : float
        Slope of f1(x).
    b1 : float
        Support of f1(x).
    m2 : float
        Slope of f2(x).
    b2 : float
        Support of f2(x).

    Returns
    -------
    x : float
        Intersection of f1(x) and f2(x).
     """
    return (b2 - b1) / (m1 - m2)


def _estimate_tangent_of_consecutive_power_laws(A_ns, G_ns):
    """
    Estimate the curve described by the intersection-points of two
    consecutive power-laws in a list of N power-laws [f_0(x), ..., f_N(x)].

    f_0(x) = A_0 * x ** G_0,
           .
           .
    f_N(x) = A_N * x ** G_N

    Parameters
    ----------
    A_ns : list of floats
        N power-laws normalizations: A_ns = [A_0, A_1, A_2, ... , A_N]
    G_ns : list of floats
        N power-laws exponents: G_ns = [G_0, G_1, G_2, ... , G_N]

    Returns
    -------
    (x, y) : (array of floats, array of floats)
        List of N-1 intersections for N power-laws.
    """
    assert len(A_ns) == len(G_ns)
    num = len(A_ns)
    assert num >= 2

    log10_A_ns = np.log10(np.array(A_ns))
    G_ns = np.array(G_ns)

    x = []
    y = []
    for i in range(num - 1):
        log10_x = _find_intersection_two_lines(
            b1=log10_A_ns[i], m1=G_ns[i], b2=log10_A_ns[i + 1], m2=G_ns[i + 1],
        )
        log10_y = log10_A_ns[i] + G_ns[i] * (log10_x)
        x.append(10 ** log10_x)
        y.append(10 ** log10_y)
    return (np.array(x), np.array(y))


def estimate_tangent_of_consecutive_power_laws(
    flux_densities, spectral_indices
):
    return _estimate_tangent_of_consecutive_power_laws(
        A_ns=flux_densities, G_ns=spectral_indices,
    )


def estimate_integral_spectral_exclusion_zone(
    signal_effective_area_m2,
    energy_bin_edges_GeV,
    critical_signal_rate_per_s,
    power_law_spectral_indices=np.linspace(start=-5, stop=-0.5, num=137),
    power_law_pivot_energy_GeV=1.0,
):
    """
    Estimates a curve in the space of differential flux and energy.
    All power-laws below this curve can not be detected.

    Parameters
    ----------
    signal_effective_area_m2 : list of N floats
        The effective area where signal is collected in the on-region.
    energy_bin_edges_GeV : list of (N+1) floats
        The edges of the energy-bins used for the effective area.
    critical_signal_rate_per_s : float
        The minimal rate of signal in the on-region $R_S$ required to
        claim a detection.
    power_law_spectral_indices : list of floats
    """
    assert critical_signal_rate_per_s > 0.0

    power_law_spectral_indices = np.array(power_law_spectral_indices)
    assert len(power_law_spectral_indices) >= 2
    assert np.all(np.gradient(power_law_spectral_indices) > 0.0)

    assert power_law_pivot_energy_GeV > 0.0

    power_law_flux_densities = integral_sensitivity.estimate_flux_densities_of_critical_power_laws(
        signal_effective_area_m2=signal_effective_area_m2,
        energy_bin_edges_GeV=energy_bin_edges_GeV,
        critical_signal_rate_per_s=critical_signal_rate_per_s,
        power_law_spectral_indices=power_law_spectral_indices,
        power_law_pivot_energy_GeV=power_law_pivot_energy_GeV,
    )

    (
        energy_GeV,
        diff_flux_per_m2_per_GeV_per_s,
    ) = estimate_tangent_of_consecutive_power_laws(
        flux_densities=power_law_flux_densities,
        spectral_indices=power_law_spectral_indices,
    )
    return energy_GeV, diff_flux_per_m2_per_GeV_per_s
