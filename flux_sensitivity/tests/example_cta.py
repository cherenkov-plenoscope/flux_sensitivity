import flux_sensitivity as fs
import binning_utils
import pkg_resources
import glob
import os
import json_numpy


RESOURCE_DIR = pkg_resources.resource_filename(
    "flux_sensitivity", os.path.join("tests", "resources",)
)
IRF_NAME = "Prod5-South-20deg-AverageAz-14MSTs37SSTs.1800s-v0.1"

out_dir = os.path.join(".", IRF_NAME)

CONFIG = {}
CONFIG["num_bins_per_decade"] = 5
CONFIG["roi_opening_deg"] = 3.25

CONFIG["detection_threshold_std"] = 5.0
CONFIG["estimator_statistics"] = "LiMaEq17"
CONFIG["systematic_uncertainty_relative"] = 1e-2
CONFIG["observation_time_s"] = 60 * 30
CONFIG["on_over_off_ratio"] = 4

cta_irf_path = os.path.join(RESOURCE_DIR, "cta", IRF_NAME + ".fits.gz")

irf = fs.io.gamma_astro_data.read_instrument_response_function(
    path=cta_irf_path
)

irf = fs.io.gamma_astro_data.average_instrument_response_over_field_of_view(
    irf=irf, roi_opening_deg=CONFIG["roi_opening_deg"]
)

energy_bin_edges_GeV = fs.io.gamma_astro_data.find_common_energy_bin_edges(
    components=irf, num_bins_per_decade=CONFIG["num_bins_per_decade"],
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


os.makedirs(out_dir, exist_ok=True)
blk = {}
blk["energy_bin_edges_GeV"] = energy_bin_edges_GeV
blk["probability_reco_given_true"] = probability_reco_given_true
blk["probability_reco_given_true_au"] = probability_reco_given_true_au
blk["signal_area_m2"] = signal_area_m2
blk["signal_area_m2_au"] = signal_area_m2_au
blk["background_rate_per_s"] = background_rate_per_s
blk["background_rate_per_s_au"] = background_rate_per_s_au
json_numpy.write(os.path.join(out_dir, "irf.json"), blk)
json_numpy.write(os.path.join(out_dir, "config.json"), CONFIG)

scenario_dir = os.path.join(out_dir, "scenarios")
os.makedirs(scenario_dir, exist_ok=True)

for scenario_key in fs.differential.SCENARIOS:
    scenario = fs.differential.init_scenario_matrices_for_signal_and_background(
        probability_reco_given_true=probability_reco_given_true,
        probability_reco_given_true_au=probability_reco_given_true_au,
        scenario_key=scenario_key,
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
        onregion_over_offregion_ratio=CONFIG["on_over_off_ratio"],
        observation_time_s=CONFIG["observation_time_s"],
        instrument_systematic_uncertainty_relative=CONFIG["systematic_uncertainty_relative"],
        detection_threshold_std=CONFIG["detection_threshold_std"],
        estimator_statistics=CONFIG["estimator_statistics"],
    )

    (dVdE, dVdE_au,) = fs.differential.estimate_differential_sensitivity(
        energy_bin_edges_GeV=energy_bin_edges_GeV,
        signal_area_in_scenario_m2=A_gamma_scenario,
        signal_area_in_scenario_m2_au=A_gamma_scenario_au,
        critical_signal_rate_in_scenario_per_s=R_gamma_scenario,
        critical_signal_rate_in_scenario_per_s_au=R_gamma_scenario_au,
    )

    out = {}
    out["scenario"] = scenario

    out["A_gamma_scenario"] = A_gamma_scenario
    out["A_gamma_scenario_au"] = A_gamma_scenario_au

    out["R_background_scenario"] = R_background_scenario
    out["R_background_scenario_au"] = R_background_scenario_au

    out["dVdE"] = dVdE
    out["dVdE_au"] = dVdE_au

    path = os.path.join(scenario_dir, "{:s}.json".format(scenario_key))
    json_numpy.write(path, out)
