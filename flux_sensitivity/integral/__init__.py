import numpy as np
import binning_utils
from . import critical_rate


def power_law_spectrum(energy, flux_density, spectral_index, pivot_energy):
    """
    power-law:
    f(energy) = flux_density * (energy/pivot_energy) ** spectral_index
    """
    return flux_density * (energy / pivot_energy) ** (spectral_index)


def estimate_signal_rate_for_power_law_per_s(
    signal_area_m2,
    energy_bin_edges_GeV,
    power_law_flux_density_per_m2_per_GeV_per_s,
    power_law_spectral_index,
    power_law_pivot_energy_GeV,
):
    """
    Estimate the rate of signal-counts in the on-region $R_S / s^{-1}$
    given the emission's energy-spectrum is a power-law:

    f(energy) = flux_density * (energy/pivot_energy) ** spectral_index

    Parameters
    ----------
    signal_area_m2 : array of N floats / m^2
        The signal's effective area in the on-region.
    energy_bin_edges_GeV : array of (N+1) floats / GeV
        The edges of the energy-bins used for the signal's effective area.
    power_law_flux_density_per_m2_per_GeV_per_s : float / m^{-2} s^{-1} (GeV)^{-1}
        The power-law's flux-density.
    power_law_spectral_index : float / 1
        The power-law's spectral-index.
    power_law_pivot_energy_GeV : float / GeV
        The power-law's pivot-energy.

    Returns
    -------
    The signal-rate $R_S$ : float
    """
    assert np.all(signal_area_m2 >= 0.0)
    assert len(energy_bin_edges_GeV) == len(signal_area_m2) + 1
    assert np.all(energy_bin_edges_GeV > 0.0)
    assert np.all(np.gradient(energy_bin_edges_GeV) > 0.0)

    assert power_law_flux_density_per_m2_per_GeV_per_s > 0.0
    assert power_law_pivot_energy_GeV > 0.0

    energy_bin_centers_GeV = binning_utils.centers(
        bin_edges=energy_bin_edges_GeV,
    )
    energy_bin_width_GeV = binning_utils.widths(
        bin_edges=energy_bin_edges_GeV,
    )

    differential_flux_per_m2_per_s_per_GeV = power_law_spectrum(
        energy=energy_bin_centers_GeV,
        flux_density=power_law_flux_density_per_m2_per_GeV_per_s,
        spectral_index=power_law_spectral_index,
        pivot_energy=power_law_pivot_energy_GeV,
    )

    differential_signal_rate_per_s_per_GeV = (
        differential_flux_per_m2_per_s_per_GeV * signal_area_m2
    )

    signal_rate_per_s = np.sum(
        differential_signal_rate_per_s_per_GeV * energy_bin_width_GeV
    )
    return signal_rate_per_s


def _relative_ratio(a, b):
    return np.abs(a - b) / (0.5 * (a + b))


def estimate_flux_densities_of_critical_power_laws(
    signal_area_m2,
    energy_bin_edges_GeV,
    critical_signal_rate_per_s,
    power_law_spectral_indices,
    power_law_pivot_energy_GeV=1.0,
    margin=1e-2,
    upper_flux_density_per_m2_per_GeV_per_s=1e6,
    max_num_iterations=10000,
):
    """
    Estimates the flux-density of the critical power-laws which are still
    detectable.

    power-law:
    f(energy) = flux_density * (energy/pivot_energy) ** spectral_index

    Parameters
    ----------
    signal_area_m2 : array of N floats
        The signal's effective area in the on-region.
    energy_bin_edges_GeV : list of (N+1) floats
        The edges of the energy-bins used for the effective area.
    critical_signal_rate_per_s : float
        The critical rate of signal in the on-region which is required
        to claim a detection.
    power_law_spectral_indices : list of floats
        A list of the power-law's spectral-indices.
    power_law_pivot_energy_GeV : float
        Same for all power-laws.
    margin : float
        Stopping-criteria for iterative search of flux-density. When the
        relative ratio of the critical_signal_rate_per_s and the acual
        signal_rate_per_s caused by a particular power-law is below this
        margin, the search is complete.
    upper_flux_density_per_m2_per_GeV_per_s : float
        Starting point for iterative search of flux-density. This should be
        larger than the expected flux.
    max_num_iterations : int
        Stop the iterative search of flux-density in any case after
        this many iterations.

    Returns
    -------
    power_law_flux_densities : list of floats
    """
    assert critical_signal_rate_per_s > 0.0

    assert power_law_pivot_energy_GeV > 0.0
    assert margin > 0.0
    assert upper_flux_density_per_m2_per_GeV_per_s > 0.0
    assert max_num_iterations > 0

    power_law_flux_densities = []

    for i, power_law_spectral_index in enumerate(power_law_spectral_indices):

        flux_dens = float(upper_flux_density_per_m2_per_GeV_per_s)

        iteration = 0
        while True:
            assert iteration < max_num_iterations

            signal_rate_per_s = estimate_signal_rate_for_power_law_per_s(
                signal_area_m2=signal_area_m2,
                energy_bin_edges_GeV=energy_bin_edges_GeV,
                power_law_flux_density_per_m2_per_GeV_per_s=flux_dens,
                power_law_spectral_index=power_law_spectral_index,
                power_law_pivot_energy_GeV=power_law_pivot_energy_GeV,
            )

            ratio = _relative_ratio(
                signal_rate_per_s, critical_signal_rate_per_s
            )

            if ratio < margin:
                break

            rr = ratio / 3
            if signal_rate_per_s > critical_signal_rate_per_s:
                flux_dens *= 1 - rr
            else:
                flux_dens *= 1 + rr

            iteration += 1

        power_law_flux_densities.append(float(flux_dens))
    return power_law_flux_densities


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


def estimate_spectral_exclusion_zone(
    signal_area_m2,
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
    signal_area_m2 : list of N floats
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

    power_law_flux_densities = estimate_flux_densities_of_critical_power_laws(
        signal_area_m2=signal_area_m2,
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
