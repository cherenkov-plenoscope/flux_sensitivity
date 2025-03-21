#!/usr/bin/python
import argparse
import flux_sensitivity as fs
import binning_utils
import os
import json_utils
import numpy as np


parser = argparse.ArgumentParser(
    prog="estimate_differential_sensitivity",
    description=(
        "Reads an instrument's response function from a "
        "gamma-astro-data-fits-file and estimates the instrument's "
        "differential sensitivity"
    ),
)

parser.add_argument(
    "irf_file",
    metavar="IRF_PATH",
    type=str,
    help=(
        "gamma-astro-data-fits-file containing the instrument's "
        "response function."
    ),
)
parser.add_argument(
    "out_dir",
    metavar="OUT_DIR",
    type=str,
    help="writes the output here.",
)

parser.add_argument(
    "--num_bins_per_decade",
    metavar="NUM",
    type=int,
    help="Number of bins per decade",
    default=5,
)
parser.add_argument(
    "--roi_opening_deg",
    metavar="ANGLE",
    type=float,
    help="Opening angle of region-of-interest in field-of-view.",
    default=3.25,
)
parser.add_argument(
    "--detection_threshold_std",
    metavar="STDDEV",
    type=float,
    help=(
        "Standard-deviation a signal has to reach above "
        "the background's noise."
    ),
    default=5.0,
)
parser.add_argument(
    "--systematic_uncertainty_relative",
    metavar="SYS",
    type=float,
    help="Relative systematic uncertainty of the instrument.",
    default=1e-2,
)
parser.add_argument(
    "--observation_time_s",
    metavar="TOBS",
    type=float,
    help="The observation-time in s.",
    default=1800,
)
parser.add_argument(
    "--on_over_off_ratio",
    metavar="ALPHA",
    type=float,
    help="Ratio of on-region over off-region.",
    default=0.2,
)
parser.add_argument(
    "--estimator_statistics",
    metavar="METHOD",
    type=str,
    help=(
        "The name of the estimator for the statistics. "
        "['sqrt', 'LiMaEq9', 'LiMaEq17']"
    ),
    default="LiMaEq17",
)


args = parser.parse_args()
irf_path = args.irf_file
IRF_NAME = os.path.basename(irf_path)
out_dir = args.out_dir

CONFIG = {}
CONFIG["name"] = IRF_NAME
CONFIG["num_bins_per_decade"] = args.num_bins_per_decade
CONFIG["roi_opening_deg"] = args.roi_opening_deg
CONFIG["detection_threshold_std"] = args.detection_threshold_std
CONFIG["estimator_statistics"] = args.estimator_statistics
CONFIG["systematic_uncertainty_relative"] = (
    args.systematic_uncertainty_relative
)
CONFIG["observation_time_s"] = args.observation_time_s
CONFIG["on_over_off_ratio"] = args.on_over_off_ratio


irf = fs.io.gamma_astro_data.read_instrument_response_function(path=irf_path)

irf = fs.io.gamma_astro_data.average_instrument_response_over_field_of_view(
    irf=irf, roi_opening_deg=CONFIG["roi_opening_deg"]
)

energy_bin_edges_GeV = fs.io.gamma_astro_data.find_common_energy_bin_edges(
    components=irf,
    num_bins_per_decade=CONFIG["num_bins_per_decade"],
)

probability_reco_given_true = (
    fs.io.gamma_astro_data.integrate_dPdMu_to_get_probability_reco_given_true(
        dPdMu=irf["energy_dispersion"]["dPdMu"],
        dPdMu_energy_bin_edges=irf["energy_dispersion"][
            "energy_bin_edges_GeV"
        ],
        dPdMu_Mu_bin_edges=irf["energy_dispersion"]["Mu_bin_edges"],
        energy_bin_edges=energy_bin_edges_GeV,
    )
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

background_rate_onregion_per_s = (
    fs.io.gamma_astro_data.integrate_background_rate_in_onregion(
        background_per_s_per_sr_per_GeV=background_per_s_per_sr_per_GeV,
        point_spread_function_sigma_deg=point_spread_function_sigma_deg,
        energy_bin_edges_GeV=energy_bin_edges_GeV,
    )
)
background_rate_onregion_per_s_au = np.zeros(
    background_rate_onregion_per_s.shape
)


os.makedirs(out_dir, exist_ok=True)
blk = {}
blk["energy_bin_edges_GeV"] = energy_bin_edges_GeV
blk["probability_reco_given_true"] = probability_reco_given_true
blk["probability_reco_given_true_au"] = probability_reco_given_true_au
blk["signal_area_m2"] = signal_area_m2
blk["signal_area_m2_au"] = signal_area_m2_au
blk["background_rate_onregion_per_s"] = background_rate_onregion_per_s
blk["background_rate_onregion_per_s_au"] = background_rate_onregion_per_s_au
json_utils.write(os.path.join(out_dir, "irf.json"), blk)
json_utils.write(os.path.join(out_dir, "config.json"), CONFIG)

scenario_dir = os.path.join(out_dir, "scenarios")
os.makedirs(scenario_dir, exist_ok=True)

for scenario_key in fs.differential.SCENARIOS:
    scenario = (
        fs.differential.init_scenario_matrices_for_signal_and_background(
            probability_reco_given_true=probability_reco_given_true,
            probability_reco_given_true_au=probability_reco_given_true_au,
            scenario_key=scenario_key,
        )
    )

    (
        signal_area_in_scenario_m2,
        signal_area_in_scenario_m2_au,
    ) = fs.differential.apply_scenario_to_signal_effective_area(
        signal_area_m2=signal_area_m2,
        signal_area_m2_au=signal_area_m2_au,
        scenario_G_matrix=scenario["G_matrix"],
        scenario_G_matrix_au=scenario["G_matrix_au"],
    )

    (
        background_rate_onregion_in_scenario_per_s,
        background_rate_onregion_in_scenario_per_s_au,
    ) = fs.differential.apply_scenario_to_background_rate(
        rate_in_reco_energy_per_s=background_rate_onregion_per_s,
        rate_in_reco_energy_per_s_au=background_rate_onregion_per_s_au,
        scenario_B_matrix=scenario["B_matrix"],
        scenario_B_matrix_au=scenario["B_matrix_au"],
    )

    (
        critical_signal_rate_in_scenario_per_s,
        critical_signal_rate_in_scenario_per_s_au,
    ) = fs.differential.estimate_critical_signal_rate_vs_energy(
        background_rate_onregion_in_scenario_per_s=background_rate_onregion_in_scenario_per_s,
        background_rate_onregion_in_scenario_per_s_au=background_rate_onregion_in_scenario_per_s_au,
        onregion_over_offregion_ratio=CONFIG["on_over_off_ratio"],
        observation_time_s=CONFIG["observation_time_s"],
        instrument_systematic_uncertainty_relative=CONFIG[
            "systematic_uncertainty_relative"
        ],
        detection_threshold_std=CONFIG["detection_threshold_std"],
        estimator_statistics=CONFIG["estimator_statistics"],
    )

    (
        dVdE_per_m2_per_GeV_per_s,
        dVdE_per_m2_per_GeV_per_s_au,
    ) = fs.differential.estimate_differential_sensitivity(
        energy_bin_edges_GeV=energy_bin_edges_GeV,
        signal_area_in_scenario_m2=signal_area_in_scenario_m2,
        signal_area_in_scenario_m2_au=signal_area_in_scenario_m2_au,
        critical_signal_rate_in_scenario_per_s=critical_signal_rate_in_scenario_per_s,
        critical_signal_rate_in_scenario_per_s_au=critical_signal_rate_in_scenario_per_s_au,
    )

    out = {}
    out["scenario"] = scenario

    out["signal_area_in_scenario_m2"] = signal_area_in_scenario_m2
    out["signal_area_in_scenario_m2_au"] = signal_area_in_scenario_m2_au

    out["background_rate_onregion_in_scenario_per_s"] = (
        background_rate_onregion_in_scenario_per_s
    )
    out["background_rate_onregion_in_scenario_per_s_au"] = (
        background_rate_onregion_in_scenario_per_s_au
    )

    out["dVdE_per_m2_per_GeV_per_s"] = dVdE_per_m2_per_GeV_per_s
    out["dVdE_per_m2_per_GeV_per_s_au"] = dVdE_per_m2_per_GeV_per_s_au

    path = os.path.join(scenario_dir, "{:s}.json".format(scenario_key))
    json_utils.write(path, out)
