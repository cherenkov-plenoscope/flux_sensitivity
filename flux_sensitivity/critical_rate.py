"""
Variables
---------
T_obs, T_\\text{obs}
    Effective observation-time.

N_on, N_\\text{on}
    On source counts.

N_off, N_\\text{off}
    On source counts.

N_S, N_\\text{S}
    Signal counts.
    N_S = N_on - alpha * N_off

hatN_B, \\hat{N}_B
    Expected num. background-counts in on-region. Included in N_on.
    hatN_B = alpha * N_off

S, S
    Significance.

"""
import numpy as np
import lima1983analysis


def estimate_critical_signal_rate(
    expected_background_rate_in_onregion_per_s,
    onregion_over_offregion_ratio,
    observation_time_s,
    instrument_systematic_uncertainty_relative,
    detection_threshold_std,
    estimator_statistics="LiMaEq17",
    mixer_statistics_systematics="hypot",
):
    """
    "Algorithm C"
    Estimate the critical rate of signal counts $R_S$ in the on-rerion
    required to claim a detection.

    Parameters
    ----------
    expected_background_rate_in_onregion_per_s : float / s^{-1}
        Expected rate of background counts in the on-region.
        This is $\\hat{R}_B = \\hat{N}_B / T_\\text{obs}$.
    onregion_over_offregion_ratio : float / 1
        Ratio of on- over off-region.
        This is $\\alpha$
    observation_time_s : float / s
        Effective time of observation, $T_\\text{obs}$.
    instrument_systematic_uncertainty_relative : float / 1
        The instrument's relative systematic uncertaintiy.
    detection_threshold_std : float / 1
        Significance in stdandard deviations.
    estimator_statistics : str
        Statistical estimator for required num. signal-counts in onregion.
        ["sqrt", "LiMaEq9", "LiMaEq17"]

    Returns
    -------
    R_S : float
        The signal's minimal rate in the on-region required to claim a
        detection.
    """
    hatR_B = float(expected_background_rate_in_onregion_per_s)
    alpha = float(onregion_over_offregion_ratio)
    T_obs = float(observation_time_s)
    S = float(detection_threshold_std)
    U_sys_rel_unc = float(instrument_systematic_uncertainty_relative)

    assert hatR_B >= 0.0
    assert T_obs > 0.0
    assert alpha > 0.0
    assert U_sys_rel_unc >= 0.0
    assert S > 0.0

    hatN_B = hatR_B * T_obs
    N_off = hatN_B / alpha

    if estimator_statistics == "sqrt":
        _N_off_std = np.sqrt(N_off)
        _N_on_std = _N_off_std * alpha
        N_S_stat = S * _N_on_std
    elif estimator_statistics == "LiMaEq9":
        N_S_stat = lima1983analysis.estimate_N_s_eq9(
            N_off=N_off, alpha=alpha, S=S,
        )
    elif estimator_statistics == "LiMaEq17":
        N_S_stat = lima1983analysis.estimate_N_s_eq17(
            N_off=N_off,
            alpha=alpha,
            S=S,
            margin=1e-4,
            max_num_iterations=10 * 1000,
        )
    else:
        raise KeyError(
            "Unknown estimator for statistics: '{:s}'".format(
                estimator_statistics
            )
        )

    N_S_sys = S * U_sys_rel_unc * hatN_B

    if mixer_statistics_systematics == "max":
        N_S = np.max([N_S_stat, N_S_sys])
    elif mixer_statistics_systematics == "hypot":
        N_S = np.hypot(N_S_stat, N_S_sys)
    else:
        raise KeyError(
            "Unknown mixing of statistics and systematics: '{:s}'".format(
                mixer_statistics_systematics
            )
        )

    R_S = N_S / T_obs
    return R_S
