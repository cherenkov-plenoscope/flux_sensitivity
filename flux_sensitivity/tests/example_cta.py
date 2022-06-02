import flux_sensitivity as fs
import binning_utils
import pkg_resources
import glob
import os


RESOURCE_DIR = pkg_resources.resource_filename(
    "flux_sensitivity", os.path.join("tests", "resources",)
)

NUM_BINS_PER_DECADE = 5
ROI_OPENING_DEG = 3.25

DETECTION_THRESHOLD_STD = 5.0
ESTIMATOR_STATISTICS = "LiMaEq17"
SYSTEMATIC_UNCERTAINTY = 1e-2
OBSERVATION_TIME = 60 * 30
ON_OVER_OFF_RATIO = 4


cta_irf_path = os.path.join(
    RESOURCE_DIR,
    "cta",
    "Prod5-South-20deg-AverageAz-14MSTs37SSTs.1800s-v0.1.fits.gz",
)

irf = fs.io.gamma_astro_data.read_instrument_response_function(
    path=cta_irf_path
)


irf = fs.io.gamma_astro_data.average_instrument_response_over_field_of_view(
    irf=irf, roi_opening_deg=ROI_OPENING_DEG
)

energy_bin_edges_GeV = fs.io.gamma_astro_data.find_common_energy_bin_edges(
    components=irf, num_bins_per_decade=NUM_BINS_PER_DECADE,
)

probability_reco_given_true = fs.io.gamma_astro_data.integrate_dPdMu_to_get_probability_reco_given_true(
    dPdMu=irf["energy_dispersion"]["dPdMu"],
    dPdMu_energy_bin_edges=irf["energy_dispersion"]["energy_bin_edges_GeV"],
    dPdMu_Mu_bin_edges=irf["energy_dispersion"]["Mu_bin_edges"],
    energy_bin_edges=energy_bin_edges_GeV,
)
probability_reco_given_true_au = np.zeros(probability_reco_given_true.shape)

signal_area_m2 = np.interp(
    x=binning_utils.centers(energy_bin_edges_GeV),
    xp=binning_utils.centers(irf["effective_area"]["energy_bin_edges_GeV"]),
    fp=irf["effective_area"]["area_m2"],
)
signal_area_m2_au = np.zeros(signal_area_m2.shape)

background_per_s_per_sr_per_GeV = np.interp(
    x=binning_utils.centers(energy_bin_edges_GeV),
    xp=binning_utils.centers(irf["background"]["energy_bin_edges_GeV"]),
    fp=irf["background"]["background_per_s_per_sr_per_GeV"],
)

point_spread_function_sigma_deg = np.interp(
    x=binning_utils.centers(energy_bin_edges_GeV),
    xp=binning_utils.centers(
        irf["point_spread_function"]["energy_bin_edges_GeV"]
    ),
    fp=irf["point_spread_function"]["sigma_deg"],
)

background_rate_per_s = fs.io.gamma_astro_data.integrate_background_rate_in_onregion(
    background_per_s_per_sr_per_GeV=background_per_s_per_sr_per_GeV,
    point_spread_function_sigma_deg=point_spread_function_sigma_deg,
    energy_bin_edges_GeV=energy_bin_edges_GeV,
)
background_rate_per_s_au = np.zeros(background_rate_per_s.shape)


for dk in fs.differential.SCENARIOS:
    print(dk)

    scenario = fs.differential.make_energy_confusion_matrices_for_signal_and_background(
        probability_reco_given_true=probability_reco_given_true,
        probability_reco_given_true_abs_unc=probability_reco_given_true_au,
        scenario_key=dk,
    )

    (
        A_gamma_scenario,
        A_gamma_scenario_au,
    ) = fs.differential.apply_scenario_to_signal_effective_area(
        signal_area_m2=signal_area_m2,
        signal_area_m2_au=signal_area_m2_au,
        scenario_G_matrix=scenario["G_matrix"],
        scenario_G_matrix_au=scenario["G_matrix_au"],
    )

    (
        R_background_scenario,
        R_background_scenario_au,
    ) = fs.differential.apply_scenario_to_background_rate(
        rate_in_reco_energy_per_s=background_rate_per_s,
        rate_in_reco_energy_per_s_au=background_rate_per_s_au,
        scenario_B_matrix=scenario["B_matrix"],
        scenario_B_matrix_au=scenario["B_matrix_au"],
    )

    (
        R_gamma_scenario,
        R_gamma_scenario_au,
    ) = fs.differential.estimate_critical_signal_rate_vs_energy(
        background_rate_onregion_in_scenario_per_s=R_background_scenario,
        background_rate_onregion_in_scenario_per_s_au=R_background_scenario_au,
        onregion_over_offregion_ratio=ON_OVER_OFF_RATIO,
        observation_time_s=OBSERVATION_TIME,
        instrument_systematic_uncertainty_relative=SYSTEMATIC_UNCERTAINTY,
        detection_threshold_std=DETECTION_THRESHOLD_STD,
        estimator_statistics=ESTIMATOR_STATISTICS,
    )

    (dVdE, dVdE_au,) = fs.differential.estimate_differential_sensitivity(
        energy_bin_edges_GeV=energy_bin_edges_GeV,
        signal_area_in_scenario_m2=A_gamma_scenario,
        signal_area_in_scenario_m2_au=A_gamma_scenario_au,
        critical_signal_rate_in_scenario_per_s=R_gamma_scenario,
        critical_signal_rate_in_scenario_per_s_au=R_gamma_scenario_au,
    )
