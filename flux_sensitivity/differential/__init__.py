import numpy as np
import os
import propagate_uncertainties
import binning_utils
from .. import critical_rate
from . import scenarios


def estimate_differential_sensitivity(
    energy_bin_edges_GeV,
    signal_area_in_scenario_m2,
    critical_signal_rate_in_scenario_per_s,
):
    """
    Estimates the differential flux-sensitivity in N energy-bin for a
    certain scenario.

    Parameters
    ----------
    energy_bin_edges_GeV : array of (N+1) floats / GeV
        Edges of the energy-bins.
    signal_area_in_scenario_m2 : array of N floats / m^{2}
        Signal's effective area for each energy-bin w.r.t. to the applied
        scenario.
    critical_signal_rate_in_scenario_per_s : array of N floats / s^{-1}
        Signal's minimal rate to claim a detection for each energy-bin w.r.t.
        the scenario, the same scenario as applied in
        signal_area_in_scenario_m2.

    Returns
    -------
    differential_flux_vs_energy : array of N floats / m^{-2} s^{-1} (GeV)^{-1}
        The minimal differential flux required to claim a detection in a
        certain scenario.
    """
    num_energy_bins = len(energy_bin_edges_GeV) - 1
    assert num_energy_bins >= 1, "Need at least two bin-edges."

    assert np.all(energy_bin_edges_GeV > 0.0), "Energy must be positive."
    assert np.all(
        np.gradient(energy_bin_edges_GeV) > 0.0
    ), "Energy bin-edges must be increasing."

    assert num_energy_bins == len(signal_area_in_scenario_m2)
    assert num_energy_bins == len(critical_signal_rate_in_scenario_per_s)

    # work
    # ----
    dVdE_per_s_per_m2_per_GeV = np.nan * np.ones(num_energy_bins)

    for ebin in range(num_energy_bins):
        dE_GeV = energy_bin_edges_GeV[ebin + 1] - energy_bin_edges_GeV[ebin]
        assert dE_GeV > 0.0
        if signal_area_in_scenario_m2[ebin] > 0:
            dV_per_s_per_m2 = (
                critical_signal_rate_in_scenario_per_s[ebin]
                / signal_area_in_scenario_m2[ebin]
            )
        else:
            dV_per_s_per_m2 = np.nan

        dVdE_per_s_per_m2_per_GeV[ebin] = dV_per_s_per_m2 / dE_GeV
    return dVdE_per_s_per_m2_per_GeV


SCENARIOS = {
    "perfect_energy": {"energy_axes_label": "",},
    "broad_spectrum": {"energy_axes_label": "reco. ",},
    "line_spectrum": {"energy_axes_label": "reco. ",},
    "bell_spectrum": {"energy_axes_label": "",},
}


def make_energy_confusion_matrices_for_signal_and_background(
    probability_reco_given_true,
    probability_reco_given_true_abs_unc,
    scenario_key="broad_spectrum",
):
    shape = probability_reco_given_true.shape
    assert probability_reco_given_true_abs_unc.shape == shape

    if scenario_key == "perfect_energy":
        G = np.eye(N=shape[0])
        G_au = np.zeros(shape=shape)  # zero uncertainty

        B = np.eye(N=shape[0])
        B_au = np.zeros(shape=shape)

    elif scenario_key == "broad_spectrum":
        G = np.array(probability_reco_given_true)
        G_au = np.array(probability_reco_given_true_abs_unc)  # adopt as is

        B = np.eye(N=shape[0])
        B_au = np.zeros(shape=shape)

    elif scenario_key == "line_spectrum":
        # only the diagonal
        eye = np.eye(N=shape[0])
        G = eye * np.diag(probability_reco_given_true)
        G_au = eye * np.diag(probability_reco_given_true_abs_unc)

        B = np.eye(N=shape[0])
        B_au = np.zeros(shape=shape)

    elif scenario_key == "bell_spectrum":
        containment = 0.68
        G = containment * np.eye(N=shape[0])  # true energy for gammas
        G_au = np.zeros(shape=shape)  # zero uncertainty

        B = scenarios.bell_spectrum.init_B(
            probability_reco_given_true=probability_reco_given_true,
            containment=containment,
        )
        B_au = np.zeros(shape=shape)

    else:
        raise KeyError("Unknown scenario_key: '{:s}'".format(scenario_key))

    return {
        "G_matrix": G,
        "G_matrix_au": G_au,
        "B_matrix": B,
        "B_matrix_au": B_au,
        "energy_axes_label": SCENARIOS[scenario_key]["energy_axes_label"],
    }


def estimate_critical_signal_rate_vs_energy(
    background_rate_onregion_per_s,
    onregion_over_offregion_ratio,
    observation_time_s,
    instrument_systematic_uncertainty_relative,
    detection_threshold_std,
    estimator_statistics,
):
    critical_signal_rate_per_s = np.nan * np.ones(
        shape=background_rate_onregion_per_s.shape
    )

    for ebin in range(len(background_rate_onregion_per_s)):
        if background_rate_onregion_per_s[ebin] > 0.0:
            critical_signal_rate_per_s[
                ebin
            ] = critical_rate.estimate_critical_signal_rate(
                background_rate_onregion_per_s=background_rate_onregion_per_s[
                    ebin
                ],
                onregion_over_offregion_ratio=onregion_over_offregion_ratio,
                observation_time_s=observation_time_s,
                instrument_systematic_uncertainty_relative=instrument_systematic_uncertainty_relative,
                detection_threshold_std=detection_threshold_std,
                estimator_statistics=estimator_statistics,
            )
        else:
            critical_signal_rate_per_s[ebin] = float("nan")
    return critical_signal_rate_per_s


def apply_scenario_to_signal_effective_area(
    signal_area_m2, signal_area_m2_au, scenario_G_matrix, scenario_G_matrix_au,
):
    """
    Apply scenario (matrix G) to the signal's effective area.

    Parameters
    ----------
    signal_area_m2 : array of N floats / m^{2}
        Signal's effective area for each energy-bin in true energy.
    signal_area_m2_au : array of N floats / m^{2}
        Absolute uncertainty of signal_area_m2.
    scenario_G_matrix : array (N x N) / 1
        The scenario's matrix 'G' to confuse the signal's effective area.
    scenario_G_matrix_au : array (N x N) / 1
        Absolute uncertainty of scenario_G_matrix.

    Returns
    -------
    scenario_signal_area_m2, scenario_signal_area_m2_au
    """
    Atrue = signal_area_m2
    Atrue_au = signal_area_m2_au
    G = scenario_G_matrix
    G_au = scenario_G_matrix_au

    assert len(Atrue) == len(Atrue_au)
    assert np.all(Atrue >= 0)
    assert np.all(Atrue_au >= 0)
    N = len(Atrue)

    assert G.shape == G_au.shape
    assert G.shape[0] == G.shape[1]
    assert np.all(G >= 0)
    assert np.all(G_au >= 0)
    assert G.shape[0] == N

    Ascenario = np.zeros(N)
    Ascenario_au = np.zeros(N)

    for ereco in range(N):
        tmp = np.zeros(N)
        tmp_au = np.zeros(N)

        for etrue in range(N):
            tmp[etrue], tmp_au[etrue] = propagate_uncertainties.prod(
                x=[G[etrue, ereco], Atrue[etrue],],
                x_au=[G_au[etrue, ereco], Atrue_au[etrue],],
            )

        Ascenario[ereco], Ascenario_au[ereco] = propagate_uncertainties.sum(
            x=tmp, x_au=tmp_au
        )
    return Ascenario, Ascenario_au


def apply_scenario_to_background_rate(
    rate_in_reco_energy_per_s,
    rate_in_reco_energy_per_s_au,
    scenario_B_matrix,
    scenario_B_matrix_au,
):
    """
    Apply scenario (matrix B) to the background's rate in reco. energy.

    Parameters
    ----------
    rate_in_reco_energy_per_s : array of N floats / m^{2}
        Background's rate for each energy-bin in reco energy.
    rate_in_reco_energy_per_s_au : array of N floats / m^{2}
        Absolute uncertainty of rate_in_reco_energy_per_s.
    scenario_B_matrix : array (N x N) / 1
        The scenario's matrix 'B' to mask the energy-bins where
        background is collected in.
    scenario_B_matrix_au : array (N x N) / 1
        Absolute uncertainty of scenario_G_matrix.

    Returns
    -------
    rate_scenario_per_s, rate_scenario_per_s_au
        The scenario's rate of background in reco. energy.
    """
    Rreco = rate_in_reco_energy_per_s
    Rreco_au = rate_in_reco_energy_per_s_au
    B = scenario_B_matrix
    B_au = scenario_B_matrix_au

    assert len(Rreco) == len(Rreco_au)
    assert np.all(Rreco >= 0)
    assert np.all(Rreco_au >= 0)
    N = len(Rreco)

    assert B.shape == B_au.shape
    assert B.shape[0] == B.shape[1]
    assert np.all(B >= 0)
    assert np.all(B_au >= 0)
    assert B.shape[0] == N

    Rscenario = np.zeros(N)
    Rscenario_au = np.zeros(N)

    for ereco in range(N):
        tmp = np.zeros(N)
        tmp_au = np.zeros(N)

        for etrue in range(N):
            tmp[etrue], tmp_au[etrue] = propagate_uncertainties.prod(
                x=[B[ereco, etrue], Rreco[etrue],],
                x_au=[B_au[ereco, etrue], Rreco_au[etrue],],
            )

        Rscenario[ereco], Rscenario_au[ereco] = propagate_uncertainties.sum(
            x=tmp, x_au=tmp_au
        )
    return Rscenario, Rscenario_au


def assert_energy_reco_given_true_ax0true_ax1reco_is_normalized(
    energy_reco_given_true_ax0true_ax1reco, margin=1e-2
):
    M = energy_reco_given_true_ax0true_ax1reco
    assert M.shape[0] == M.shape[1]
    num_energy_bins = M.shape[0]

    for etrue in range(num_energy_bins):
        check = np.sum(M[etrue, :])
        if check > 0:
            assert (
                (1.0 - margin) < check < (1.0 + margin)
            ), "Expected sum(P(reco|true)) = 1, but it's {:f}".format(check)


def estimate_rate_in_reco_energy(
    energy_bin_edges_GeV,
    acceptance_m2_sr,
    acceptance_m2_sr_au,
    differential_flux_per_m2_per_sr_per_s_per_GeV,
    differential_flux_per_m2_per_sr_per_s_per_GeV_au,
    energy_reco_given_true_ax0true_ax1reco,
    energy_reco_given_true_ax0true_ax1reco_au,
):
    num_energy_bins = len(energy_bin_edges_GeV) - 1
    dE_GeV = binning_utils.widths(bin_edges=energy_bin_edges_GeV)
    dE_GeV_au = np.zeros(num_energy_bins)  # abs. uncertainty is zero

    rate_in_reco_energy_per_s = np.zeros(num_energy_bins)
    rate_in_reco_energy_per_s_au = np.zeros(num_energy_bins)

    for ereco in range(num_energy_bins):
        _summands = np.zeros(num_energy_bins)
        _summands_au = np.zeros(num_energy_bins)
        for etrue in range(num_energy_bins):
            _multiplicands = [
                differential_flux_per_m2_per_sr_per_s_per_GeV[etrue],
                energy_reco_given_true_ax0true_ax1reco[etrue, ereco],
                acceptance_m2_sr[etrue],
                dE_GeV[etrue],
            ]
            _multiplicands_au = [
                differential_flux_per_m2_per_sr_per_s_per_GeV_au[etrue],
                energy_reco_given_true_ax0true_ax1reco_au[etrue, ereco],
                acceptance_m2_sr_au[etrue],
                dE_GeV_au[etrue],
            ]
            (
                _summands[etrue],
                _summands_au[etrue],
            ) = propagate_uncertainties.prod(
                x=_multiplicands, x_au=_multiplicands_au
            )
        (
            rate_in_reco_energy_per_s[ereco],
            rate_in_reco_energy_per_s_au[ereco],
        ) = propagate_uncertainties.sum(x=_summands, x_au=_summands_au)

    return rate_in_reco_energy_per_s, rate_in_reco_energy_per_s_au


def estimate_rate_in_true_energy(
    energy_bin_edges_GeV,
    acceptance_m2_sr,
    acceptance_m2_sr_au,
    differential_flux_per_m2_per_sr_per_s_per_GeV,
    differential_flux_per_m2_per_sr_per_s_per_GeV_au,
):
    num_energy_bins = len(energy_bin_edges_GeV) - 1
    dE_GeV = binning_utils.widths(bin_edges=energy_bin_edges_GeV)
    dE_GeV_au = np.zeros(num_energy_bins)  # abs. uncertainty is zero

    rate_in_true_energy_per_s = np.zeros(num_energy_bins)
    rate_in_true_energy_per_s_au = np.zeros(num_energy_bins)

    for etrue in range(num_energy_bins):
        (
            rate_in_true_energy_per_s[etrue],
            rate_in_true_energy_per_s_au[etrue],
        ) = propagate_uncertainties.prod(
            x=[
                differential_flux_per_m2_per_sr_per_s_per_GeV[etrue],
                acceptance_m2_sr[etrue],
                dE_GeV[etrue],
            ],
            x_au=[
                differential_flux_per_m2_per_sr_per_s_per_GeV_au[etrue],
                acceptance_m2_sr_au[etrue],
                dE_GeV_au[etrue],
            ],
        )

    return rate_in_true_energy_per_s, rate_in_true_energy_per_s_au


def assert_integral_rates_are_similar_in_reco_and_true_energy(
    rate_in_reco_energy_per_s, rate_in_true_energy_per_s, margin=0.3,
):
    """
    Integral rate over all energy-bins must not change (much) under
    energy migration.
    """
    integral_Rtrue = np.sum(rate_in_true_energy_per_s[:])
    integral_Rreco = np.sum(rate_in_reco_energy_per_s[:])
    assert (1.0 - margin) < integral_Rtrue / integral_Rreco < (1.0 + margin), (
        "Expected integral rates in both true and reco. energy to be "
        "roughly the same. But they are not."
    )
