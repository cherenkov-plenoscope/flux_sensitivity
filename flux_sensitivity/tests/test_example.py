import flux_sensitivity
import pkg_resources
import os
import pandas
import astropy
from astropy.io import fits


RESOURCE_DIR = pkg_resources.resource_filename(
    "flux_sensitivity", os.path.join("tests", "resources",)
)


def read_csv(fname):
    arr = pandas.read_csv(
        os.path.join(RESOURCE_DIR, "portal", fname)
    ).to_numpy()
    if len(arr.shape) == 2:
        if arr.shape[1] == 1:
            arr = arr.reshape((arr.shape[0],))
    return arr


energy_bin_edges_GeV = read_csv("energy_bin_edges_GeV.csv")

A_gamma_m2 = read_csv("effective_area_gamma_m2.csv")
A_gamma_m2_au = read_csv("effective_area_gamma_m2_absolute_uncertainty.csv")

Q_m2_sr = {}
Q_m2_sr_au = {}
for pk in ["proton", "electron", "helium"]:
    Q_m2_sr[pk] = read_csv("effective_acceptance_{:s}_m2_sr.csv".format(pk))
    Q_m2_sr_au[pk] = read_csv(
        "effective_acceptance_{:s}_m2_sr_absolute_uncertainty.csv".format(pk),
    )


M = {}
M_au = {}
for pk in ["gamma", "proton", "electron", "helium"]:
    M[pk] = read_csv(
        "energy_conditional_probability_reco_given_true_ax0true_ax1reco_{:s}.csv".format(
            pk
        ),
    )

    M_au[pk] = read_csv(
        "energy_conditional_probability_reco_given_true_ax0true_ax1reco_{:s}_au.csv".format(
            pk
        ),
    )


dFdE_per_m2_per_sr_per_GeV_per_s = {}
for pk in ["proton", "electron", "helium"]:
    dFdE_per_m2_per_sr_per_GeV_per_s[pk] = read_csv(
        "differential_flux_{:s}_per_m2_per_sr_per_GeV_per_s.csv".format(pk),
    )


def test_resources():
    print(dFdE_per_m2_per_sr_per_GeV_per_s)

    for dk in flux_sensitivity.differential.SCENARIOS:
        pass

    print("READ CTA")

    f = astropy.io.fits.open(
        os.path.join(
            RESOURCE_DIR,
            "cta",
            "Prod5-South-20deg-AverageAz-14MSTs37SSTs.1800s-v0.1.fits.gz",
        )
    )
    f.readall()
    f.info()
    f.close()
