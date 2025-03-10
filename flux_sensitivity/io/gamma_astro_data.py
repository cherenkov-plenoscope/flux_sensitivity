import astropy
from astropy.io import fits
import binning_utils as bu
import numpy as np


def average_field_of_view_radial(X, theta_bin_edges_deg, roi_opening_deg):
    avg_shape = X.shape[1:]
    X_avg = np.zeros(avg_shape)
    total_theta_solid_angle_deg2 = 0.0
    for t in range(len(theta_bin_edges_deg) - 1):
        theta_start_deg = theta_bin_edges_deg[t]
        theta_stop_deg = theta_bin_edges_deg[t + 1]
        theta_solid_angle_deg2 = np.pi * (
            theta_stop_deg**2 - theta_start_deg**2
        )
        theta_deg = 0.5 * (theta_start_deg + theta_stop_deg)
        if theta_deg <= roi_opening_deg:
            total_theta_solid_angle_deg2 += theta_solid_angle_deg2
            X_avg += theta_solid_angle_deg2 * X[t, :]

    X_avg /= total_theta_solid_angle_deg2
    return X_avg


def average_field_of_view_grid(
    X, detx_bin_edges_deg, dety_bin_edges_deg, roi_opening_deg
):
    avg_shape = X.shape[2:]
    X_avg = np.zeros(avg_shape)

    N_detx = len(detx_bin_edges_deg) - 1
    N_dety = len(dety_bin_edges_deg) - 1

    total_solid_angle = 0

    for x in range(N_detx):
        for y in range(N_dety):
            detx_start_deg = detx_bin_edges_deg[x]
            detx_stop_deg = detx_bin_edges_deg[x + 1]
            dety_start_deg = dety_bin_edges_deg[x]
            dety_stop_deg = dety_bin_edges_deg[x + 1]

            detx_deg = 0.5 * (detx_start_deg + detx_stop_deg)
            dety_deg = 0.5 * (dety_start_deg + dety_stop_deg)

            det_solid_angle = (detx_stop_deg - detx_start_deg) * (
                dety_stop_deg - dety_start_deg
            )
            if np.hypot(detx_deg, dety_deg) < roi_opening_deg:
                total_solid_angle += det_solid_angle
                X_avg += det_solid_angle * X[x, y, :]

    X_avg /= total_solid_angle
    return X_avg


def find_common_energy_bin_edges(components, num_bins_per_decade):
    multiple_edges = [
        components[k]["energy_bin_edges_GeV"] for k in components
    ]

    E_start = bu.max_lowest_edge(multiple_edges=multiple_edges)
    E_stop = bu.min_highest_edge(multiple_edges=multiple_edges)

    (start_decade, start_bin) = bu.power10.find_upper_decade_and_bin(
        x=E_start, num_bins_per_decade=num_bins_per_decade
    )
    (stop_decade, stop_bin) = bu.power10.find_lower_decade_and_bin(
        x=E_stop, num_bins_per_decade=num_bins_per_decade
    )

    return bu.power10.space(
        start_decade=start_decade,
        start_bin=start_bin,
        stop_decade=stop_decade,
        stop_bin=stop_bin,
        num_bins_per_decade=num_bins_per_decade,
    )


def _read_energy_dispersion(hdu):
    hdu
    return {
        "energy_bin_edges_GeV": bu.merge_low_high_edges(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "Mu_bin_edges": bu.merge_low_high_edges(
            low=hdu.data["MIGRA_LO"][0, :],
            high=hdu.data["MIGRA_HI"][0, :],
        ),
        "theta_bin_edges_deg": bu.merge_low_high_edges(
            low=hdu.data["THETA_LO"][0, :],
            high=hdu.data["THETA_HI"][0, :],
        ),
        "dPdMu_vs_theta": hdu.data["MATRIX  "][0],
    }


def _read_effective_area(hdu):
    return {
        "energy_bin_edges_GeV": bu.merge_low_high_edges(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "theta_bin_edges_deg": bu.merge_low_high_edges(
            low=hdu.data["THETA_LO"][0, :],
            high=hdu.data["THETA_HI"][0, :],
        ),
        "area_vs_theta_m2": hdu.data["EFFAREA"][0],
    }


def _read_background(hdu):
    per_MeV_to_per_GeV = 1e3
    return {
        "energy_bin_edges_GeV": bu.merge_low_high_edges(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "detx_bin_edges_deg": bu.merge_low_high_edges(
            low=hdu.data["DETX_LO"][0, :],
            high=hdu.data["DETX_HI"][0, :],
        ),
        "dety_bin_edges_deg": bu.merge_low_high_edges(
            low=hdu.data["DETY_LO"][0, :],
            high=hdu.data["DETY_HI"][0, :],
        ),
        "background_vs_detx_vs_dety_per_s_per_sr_per_GeV": per_MeV_to_per_GeV
        * hdu.data["BKG"][0],
    }


def _read_point_spread_function(hdu):
    return {
        "energy_bin_edges_GeV": bu.merge_low_high_edges(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "theta_bin_edges_deg": bu.merge_low_high_edges(
            low=hdu.data["THETA_LO"][0, :],
            high=hdu.data["THETA_HI"][0, :],
        ),
        "sigma_vs_theta_deg": hdu.data["SIGMA_1"][0],
    }


def read_instrument_response_function(path):
    irf = {}
    with astropy.io.fits.open(path) as f:
        irf["effective_area"] = _read_effective_area(hdu=f["EFFECTIVE AREA"])
        irf["energy_dispersion"] = _read_energy_dispersion(
            hdu=f["ENERGY DISPERSION"]
        )
        irf["background"] = _read_background(hdu=f["BACKGROUND"])
        irf["point_spread_function"] = _read_point_spread_function(
            hdu=f["POINT SPREAD FUNCTION"]
        )

    return irf


def average_instrument_response_over_field_of_view(irf, roi_opening_deg):
    assert roi_opening_deg > 0.0

    irf["effective_area"]["area_m2"] = average_field_of_view_radial(
        X=irf["effective_area"]["area_vs_theta_m2"],
        theta_bin_edges_deg=irf["effective_area"]["theta_bin_edges_deg"],
        roi_opening_deg=roi_opening_deg,
    )

    irf["energy_dispersion"]["dPdMu"] = average_field_of_view_radial(
        X=irf["energy_dispersion"]["dPdMu_vs_theta"],
        theta_bin_edges_deg=irf["energy_dispersion"]["theta_bin_edges_deg"],
        roi_opening_deg=roi_opening_deg,
    )

    irf["background"][
        "background_per_s_per_sr_per_GeV"
    ] = average_field_of_view_grid(
        X=irf["background"]["background_vs_detx_vs_dety_per_s_per_sr_per_GeV"],
        detx_bin_edges_deg=irf["background"]["detx_bin_edges_deg"],
        dety_bin_edges_deg=irf["background"]["dety_bin_edges_deg"],
        roi_opening_deg=roi_opening_deg,
    )

    irf["point_spread_function"]["sigma_deg"] = average_field_of_view_radial(
        X=irf["point_spread_function"]["sigma_vs_theta_deg"],
        theta_bin_edges_deg=irf["point_spread_function"][
            "theta_bin_edges_deg"
        ],
        roi_opening_deg=roi_opening_deg,
    )

    return irf


def integrate_dPdMu_to_get_probability_reco_given_true(
    dPdMu,
    dPdMu_energy_bin_edges,
    dPdMu_Mu_bin_edges,
    energy_bin_edges,
):
    """
    Convert a dispersion-matrux dPdMu to a conditional probability-matrix.

    dP/dMu is the probability-density-function and must fulfill:

    Int_0^\\inf \\frac{dP}{d\\mu} d\\mu = 1

    with \\mu being the ratio:

    \\mu = \\frac{E_{reco}}{E_{true}}.

    The conditional probability R to get E_{reco} given E_{true} is:

    R(i, j) = (\\Delta E_{true})^{-1} Int_\\Delta E_{true} R(i, E_{true}) dE_{true}

    with:

    R(i, E_true) = Int_\\mu(\\Delta E_{reco}) dP/d\\mu(E_true, \\mu) d\\mu

    Parameters
    ----------
    dPdMu : array(N, M)
        The dispersion-matrix.
    dPdMu_Mu_bin_edges : array(N)
        The true energy bin-edges used in dPdMu.
    dPdMu_energy_bin_edges : array(M)
        The margin's \\mu bin-edges used in dPdMu.
        This range sould be in awide range such as 0.2 < \\mu < 5.
    energy_bin_edges : array(J)
        The targeted energy-bin-edges of the conditional-probability-matrix
        to be computed.

    Returns
    -------
    conditional-probability-matrix : array(J, J)
        The conditional probability to measure E_{reco} given E_{true}.
        Using energy_bin_edges.
        ax0: true
        ax1: reco
    """
    assert bu.is_strictly_monotonic_increasing(dPdMu_Mu_bin_edges)
    assert bu.is_strictly_monotonic_increasing(dPdMu_energy_bin_edges)
    assert bu.is_strictly_monotonic_increasing(energy_bin_edges)

    assert np.all(dPdMu_Mu_bin_edges >= 0.0)
    assert np.all(dPdMu_energy_bin_edges >= 0.0)
    assert np.all(energy_bin_edges >= 0.0)

    assert np.all(dPdMu >= 0.0)

    dPdMu_Mu_bin_centers = bu.centers(dPdMu_Mu_bin_edges)
    Mu_bin_widths = bu.widths(dPdMu_Mu_bin_edges)
    dPdMu_energy_bin_centers = bu.centers(dPdMu_energy_bin_edges)
    dPdMu_energy_bin_widths = bu.widths(dPdMu_energy_bin_edges)
    energy_bin_centers = bu.centers(energy_bin_edges)

    N = len(dPdMu_Mu_bin_edges) - 1
    M = len(dPdMu_energy_bin_edges) - 1
    J = len(energy_bin_edges) - 1

    R = np.nan * np.ones(shape=(J, J))

    for r in range(J):
        for t in range(J):
            Etrue_center = energy_bin_centers[t]
            Ereco_start = energy_bin_edges[r]
            Ereco_stop = energy_bin_edges[r + 1]
            assert Ereco_stop > Ereco_start

            Mu_start = Ereco_start / Etrue_center
            Mu_stop = Ereco_stop / Etrue_center

            assert Mu_start < Mu_stop
            Mu_mask = np.logical_and(
                dPdMu_Mu_bin_centers >= Mu_start,
                dPdMu_Mu_bin_centers < Mu_stop,
            )

            # integrate over Mu
            # -----------------
            R_I_Etrue = np.zeros(dPdMu.shape[0])
            for u in range(N):
                R_I_Etrue += dPdMu[u, :] * Mu_mask[u] * Mu_bin_widths[u]

            assert R_I_Etrue.shape[0] == dPdMu.shape[0]

            Etrue_start = energy_bin_edges[t]
            Etrue_stop = energy_bin_edges[t + 1]
            assert Etrue_stop > Etrue_start
            Etrue_mask = np.logical_and(
                dPdMu_energy_bin_centers >= Etrue_start,
                dPdMu_energy_bin_centers < Etrue_stop,
            )
            delta_Etrue = Etrue_stop - Etrue_start

            # integrate over Etrue
            # --------------------
            R[t, r] = 0.0
            for et in range(M):
                R[t, r] += (
                    R_I_Etrue[et]
                    * Etrue_mask[et]
                    * dPdMu_energy_bin_widths[et]
                )
            R[t, r] /= delta_Etrue

    # normalize
    # ---------
    # ax0 -> true
    # ax1 -> reco
    probability_reco_given_true = np.nan * np.ones(shape=(J, J))

    for t in range(J):
        norm = np.sum(R[t, :])
        for r in range(J):
            probability_reco_given_true[t, r] = R[t, r] / norm

    # check
    for t in range(J):
        if np.sum(probability_reco_given_true[t, :]) > 0.0:
            assert 0.99 < np.sum(probability_reco_given_true[t, :]) < 1.01

    return probability_reco_given_true


def integrate_background_rate_in_onregion(
    background_per_s_per_sr_per_GeV,
    point_spread_function_sigma_deg,
    energy_bin_edges_GeV,
):
    psf_sigma_rad = np.deg2rad(point_spread_function_sigma_deg)
    energy_bin_widths_GeV = bu.widths(energy_bin_edges_GeV)

    background_rate_per_s = np.zeros(len(energy_bin_edges_GeV) - 1)
    for e in range(len(energy_bin_edges_GeV) - 1):
        psf_solid_angle_sr = np.pi * psf_sigma_rad[e] ** 2
        background_rate_per_s[e] = (
            background_per_s_per_sr_per_GeV[e]
            * energy_bin_widths_GeV[e]
            * psf_solid_angle_sr
        )
    return background_rate_per_s
