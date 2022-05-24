import numpy as np
import os
import propagate_uncertainties as pru
from . import integral_sensitivity
from . import critical_rate


def estimate_differential_sensitivity(
    energy_bin_edges_GeV, signal_effective_area_m2, critical_signal_rate_per_s,
):
    """
    Estimates the differential flux-sensitivity in N energy-bin.

    Parameters
    ----------
    energy_bin_edges_GeV : array of (N+1) floats / GeV
        The edges of the energy-bins.
    signal_effective_area_m2 : array of N floats / m^{2}
        The signal's effective area for each energy-bin.
    critical_signal_rate_per_s : array of N floats / s^{-1}
        The signal's minimal rate to claim a detection for each energy-bin.

    Returns
    -------
    differential_flux_vs_energy : array of N floats / m^{-2} s^{-1} (GeV)^{-1}
        The minimal differential flux required to claim a detection in an
        energy-bin.
    """
    num_energy_bins = len(energy_bin_edges_GeV) - 1
    assert num_energy_bins >= 1, "Need at least two bin-edges."

    assert np.all(energy_bin_edges_GeV > 0.0), "Energy must be positive."
    assert np.all(
        np.gradient(energy_bin_edges_GeV) > 0.0
    ), "Energy bin-edges must be increasing."

    assert num_energy_bins == len(signal_effective_area_m2)
    assert num_energy_bins == len(critical_signal_rate_per_s)

    # work
    # ----
    dVdE_per_s_per_m2_per_GeV = np.nan * np.ones(num_energy_bins)

    for ebin in range(num_energy_bins):
        dE_GeV = energy_bin_edges_GeV[ebin + 1] - energy_bin_edges_GeV[ebin]
        assert dE_GeV > 0.0
        if signal_effective_area_m2[ebin] > 0:
            dV_per_s_per_m2 = (
                critical_signal_rate_per_s[ebin]
                / signal_effective_area_m2[ebin]
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

        B = make_mask_for_energy_confusion_matrix_for_bell_spectrum(
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
    expected_background_rate_in_onregion_per_s,
    onregion_over_offregion_ratio,
    observation_time_s,
    instrument_systematic_uncertainty_relative,
    detection_threshold_std,
    estimator_statistics,
):
    critical_signal_rate_per_s = np.nan * np.ones(
        shape=expected_background_rate_in_onregion_per_s.shape
    )

    for ebin in range(len(expected_background_rate_in_onregion_per_s)):
        if expected_background_rate_in_onregion_per_s[ebin] > 0.0:
            critical_signal_rate_per_s[
                ebin
            ] = critical_rate.estimate_critical_signal_rate(
                expected_background_rate_in_onregion_per_s=expected_background_rate_in_onregion_per_s[
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


def next_containment_and_weight(
    accumulated_containment, bin_containment, target_containment,
):
    assert 0 <= accumulated_containment <= 1
    assert 0 <= bin_containment <= 1
    assert 0 < target_containment <= 1

    missing_containment = target_containment - accumulated_containment
    assert missing_containment > 0

    if bin_containment > 0:
        weight = np.min([missing_containment / bin_containment, 1])
    else:
        weight = 0

    if weight == 1:
        return accumulated_containment + bin_containment, 1
    else:
        return target_containment, weight


def make_mask_for_energy_confusion_matrix_for_bell_spectrum(
    probability_reco_given_true, containment=0.68
):
    # ax0 -> true
    # ax1 -> reco
    num_bins = probability_reco_given_true.shape[0]
    M = probability_reco_given_true
    mask = np.zeros(shape=(num_bins, num_bins))

    # estimate containment regions:
    for etrue in range(num_bins):
        if np.sum(M[etrue, :]) > 0.0:
            assert 0.99 < np.sum(M[etrue, :]) < 1.01

            accumulated_containment = 0.0
            reco_best = np.argmax(M[etrue, :])

            accumulated_containment, weight = next_containment_and_weight(
                accumulated_containment=accumulated_containment,
                bin_containment=M[etrue, reco_best],
                target_containment=containment,
            )

            mask[etrue, reco_best] = weight
            start = reco_best - 1
            stop = reco_best + 1
            i = 0
            while accumulated_containment < containment:
                if start > 0:
                    accumulated_containment, w = next_containment_and_weight(
                        accumulated_containment=accumulated_containment,
                        bin_containment=M[etrue, start],
                        target_containment=containment,
                    )
                    mask[etrue, start] = w
                    start -= 1
                if accumulated_containment == containment:
                    break

                if stop + 1 < num_bins:
                    accumulated_containment, w = next_containment_and_weight(
                        accumulated_containment=accumulated_containment,
                        bin_containment=M[etrue, stop],
                        target_containment=containment,
                    )
                    mask[etrue, stop] = w
                    stop += 1
                if accumulated_containment == containment:
                    break

                if start == 0 and stop + 1 == num_bins:
                    break

                i += 1
                assert i < 2 * num_bins
    return mask


def make_area_in_reco_energy(
    area, area_au, G_matrix, G_matrix_au,
):
    A = area
    A_au = area_au
    G = G_matrix
    G_au = G_matrix_au
    assert len(A) == len(A_au)
    assert np.all(A >= 0)
    assert np.all(A_au >= 0)

    assert G.shape == G_au.shape
    assert G.shape[0] == G.shape[1]
    assert np.all(G >= 0)
    assert np.all(G_au >= 0)

    assert len(A) == G.shape[0]

    num_bins = len(A)
    A_out = np.zeros(num_bins)
    A_out_au = np.zeros(num_bins)

    for er in range(num_bins):
        tmp = np.zeros(num_bins)
        tmp_au = np.zeros(num_bins)

        for et in range(num_bins):
            tmp[et], tmp_au[et] = pru.prod(
                x=[G[et, er], A[et],], x_au=[G_au[et, er], A_au[et],],
            )

        A_out[er], A_out_au[er] = pru.sum(x=tmp, x_au=tmp_au)
    return A_out, A_out_au


def integrate_rates_in_reco_energy_with_mask(
    Rreco, Rreco_au, integration_mask, integration_mask_au
):
    assert len(Rreco) == len(Rreco_au)
    assert np.all(Rreco >= 0)
    assert np.all(Rreco_au >= 0)
    num_energy_bins = len(Rreco)

    assert integration_mask.shape == integration_mask_au.shape
    assert np.all(integration_mask >= 0)
    assert np.all(integration_mask_au >= 0)

    assert integration_mask.shape[0] == integration_mask.shape[1]
    assert integration_mask.shape[0] == num_energy_bins

    imask = integration_mask
    imask_au = integration_mask_au

    Rreco_total = np.zeros(num_energy_bins)
    Rreco_total_au = np.zeros(num_energy_bins)

    for ereco in range(num_energy_bins):
        tmp_sum = np.zeros(num_energy_bins)
        tmp_sum_au = np.zeros(num_energy_bins)
        for etrue in range(num_energy_bins):
            tmp_sum[etrue], tmp_sum_au[etrue] = pru.prod(
                x=[imask[ereco, etrue], Rreco[etrue],],
                x_au=[imask_au[ereco, etrue], Rreco_au[etrue],],
            )
        Rreco_total[ereco], Rreco_total_au[ereco] = pru.sum(
            x=tmp_sum, x_au=tmp_sum_au
        )
    return Rreco_total, Rreco_total_au
