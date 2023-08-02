import lima1983analysis
import propagate_uncertainties as pru


def estimate_critical_signal_rate(
    background_rate_onregion_in_scenario_per_s,
    background_rate_onregion_in_scenario_per_s_au,
    onregion_over_offregion_ratio,
    observation_time_s,
    instrument_systematic_uncertainty_relative,
    detection_threshold_std,
    estimator_statistics="LiMaEq17",
    combiner_statistics_systematics="hypot",
):
    """
    "Algorithm C"
    Estimate the critical rate of signal counts $R_S$ in the on-rerion
    required to claim a detection.

    Parameters
    ----------
    background_rate_onregion_in_scenario_per_s : float / s^{-1}
        Expected rate of background counts in the on-region for a certain
        scenario.
        This is $\\hat{R}_B = \\hat{N}_B / T_\\text{obs}$.
    background_rate_onregion_in_scenario_per_s_au : float / s^{-1}
        Absulute uncertainty.
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
    combiner_statistics_systematics : str
        ["max", "hypot"]

    Returns
    -------
    R_S, R_S_au : (float, float)
        The signal's minimal rate in the on-region required to claim a
        detection.
    """
    hatR_B = float(background_rate_onregion_in_scenario_per_s)
    hatR_B_au = float(background_rate_onregion_in_scenario_per_s_au)
    alpha = float(onregion_over_offregion_ratio)
    T_obs = float(observation_time_s)
    S = float(detection_threshold_std)
    U_sys_rel_unc = float(instrument_systematic_uncertainty_relative)

    assert hatR_B >= 0.0
    assert hatR_B_au >= 0.0
    assert T_obs > 0.0
    assert alpha > 0.0
    assert U_sys_rel_unc >= 0.0
    assert S > 0.0

    hatN_B, hatN_B_au = pru.multiply(x=hatR_B, x_au=hatR_B_au, y=T_obs, y_au=0)
    N_off, N_off_au = pru.divide(x=hatN_B, x_au=hatN_B_au, y=alpha, y_au=0)

    if estimator_statistics == "sqrt":
        N_S_stat, N_S_stat_au = estimator_sqrt(
            N_off=N_off, N_off_au=N_off_au, alpha=alpha, S=S,
        )
    elif estimator_statistics == "LiMaEq9":
        N_S_stat, N_S_stat_au = estimator_lima1983analysis_N_s_eq9(
            N_off=N_off, N_off_au=N_off_au, alpha=alpha, S=S,
        )
    elif estimator_statistics == "LiMaEq17":
        N_S_stat, N_S_stat_au = estimator_lima1983analysis_N_s_eq17(
            N_off=N_off,
            N_off_au=N_off_au,
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

    # estimate the required signal-rate to satisfy the systematic uncertainty
    N_s_sys, N_s_sys_au = pru.multiply(
        x=(S * U_sys_rel_unc), x_au=0.0, y=hatN_B, y_au=hatN_B_au
    )

    # combine statistics and systematics
    if combiner_statistics_systematics == "max":
        N_s, N_s_au = pru.max(
            x=[N_S_stat, N_s_sys], x_au=[N_S_stat_au, N_s_sys_au]
        )
    elif combiner_statistics_systematics == "hypot":
        N_s, N_s_au = pru.hypot(
            x=N_S_stat, x_au=N_S_stat_au, y=N_s_sys, y_au=N_s_sys_au
        )
    else:
        raise KeyError(
            "Unknown mixing of statistics and systematics: '{:s}'".format(
                combiner_statistics_systematics
            )
        )

    R_s, R_s_au = pru.divide(x=N_s, x_au=N_s_au, y=T_obs, y_au=0.0)
    return R_s, R_s_au


def estimator_sqrt(N_off, N_off_au, alpha, S):
    N_off_std, N_off_std_au = pru.sqrt(x=N_off, x_au=N_off_au)
    N_on_std, N_on_std_au = pru.multiply(
        x=N_off_std, x_au=N_off_std_au, y=alpha, y_au=0.0
    )
    return pru.multiply(x=S, x_au=0.0, y=N_on_std, y_au=N_on_std_au)


def estimator_lima1983analysis_N_s_eq9(
    N_off, N_off_au, alpha, S, numeric_derivative_epsilon_relative=1e-2
):
    lima_eq9 = lima1983analysis.estimate_N_s_eq9
    N_off_plus = N_off * (1.0 + numeric_derivative_epsilon_relative)
    N_off_minus = N_off * (1.0 - numeric_derivative_epsilon_relative)
    N_s_plus = lima_eq9(N_off=N_off_plus, alpha=alpha, S=S)
    N_s_minus = lima_eq9(N_off=N_off_minus, alpha=alpha, S=S)
    dN_s_dN_off = derivative(
        f_plus=N_s_plus,
        f_minus=N_s_minus,
        x_plus=N_off_plus,
        x_minus=N_off_minus,
    )
    N_s = lima_eq9(N_off=N_off, alpha=alpha, S=S)
    N_s_au = pru.auN(dfdx=[dN_s_dN_off], x_au=[N_off_au])
    return N_s, N_s_au


def estimator_lima1983analysis_N_s_eq17(
    N_off,
    N_off_au,
    alpha,
    S,
    margin=1e-4,
    max_num_iterations=10 * 1000,
    numeric_derivative_epsilon_relative=1e-2,
):
    lima_eq17 = lima1983analysis.estimate_N_s_eq17
    N_off_plus = N_off * (1.0 + numeric_derivative_epsilon_relative)
    N_off_minus = N_off * (1.0 - numeric_derivative_epsilon_relative)
    N_s_plus = lima_eq17(
        N_off=N_off_plus,
        alpha=alpha,
        S=S,
        margin=margin,
        max_num_iterations=max_num_iterations,
    )
    N_s_minus = lima_eq17(
        N_off=N_off_minus,
        alpha=alpha,
        S=S,
        margin=margin,
        max_num_iterations=max_num_iterations,
    )
    dN_s_dN_off = derivative(
        f_plus=N_s_plus,
        f_minus=N_s_minus,
        x_plus=N_off_plus,
        x_minus=N_off_minus,
    )
    N_s = lima_eq17(
        N_off=N_off,
        alpha=alpha,
        S=S,
        margin=margin,
        max_num_iterations=max_num_iterations,
    )
    N_s_au = pru.auN(dfdx=[dN_s_dN_off], x_au=[N_off_au])
    return N_s, N_s_au


def derivative(f_plus, f_minus, x_plus, x_minus):
    df = f_plus - f_minus
    dx = x_plus - x_minus
    return df / dx
