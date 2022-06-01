import flux_sensitivity as fs
import pkg_resources
import glob
import os


RESOURCE_DIR = pkg_resources.resource_filename(
    "flux_sensitivity", os.path.join("tests", "resources",)
)

cta_irf_path = os.path.join(
    RESOURCE_DIR,
    "cta",
    "Prod5-South-20deg-AverageAz-14MSTs37SSTs.1800s-v0.1.fits.gz",
)

print("READ CTA")

ROI_OPENING_DEG = 3.25

irf = fs.io.gamma_astro_data.read_instrument_response_function(
    path=cta_irf_path
)
irf = fs.io.gamma_astro_data.average_instrument_response_over_field_of_view(
    irf=irf, roi_opening_deg=ROI_OPENING_DEG
)

energy_bin_edges_GeV = fs.io.gamma_astro_data.find_common_energy_bin_edges(
    components=irf
)

probability_reco_given_true = fs.io.gamma_astro_data.integrate_dPdMu_to_get_probability_reco_given_true(
    dPdMu=irf["energy_dispersion"]["dPdMu"],
    dPdMu_energy_bin_edges=irf["energy_dispersion"]["energy_bin_edges_GeV"],
    dPdMu_Mu_bin_edges=irf["energy_dispersion"]["Mu_bin_edges"],
    energy_bin_edges=energy_bin_edges_GeV,
)
