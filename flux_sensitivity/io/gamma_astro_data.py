import astropy
from astropy.io import fits
import scipy
from scipy.interpolate import interp2d
import binning_utils
import numpy as np


def relative_deviation(a, b):
    if a == 0 and b == 0:
        return 0
    return np.abs(a - b) / np.abs(0.5 * (a + b))


def merge_bin_edges_low_high(low, high):
    assert len(low) == len(high)
    N = len(low)
    bin_edges = np.zeros(N + 1)
    for n in range(N):
        bin_edges[n] = low[n]
    bin_edges[N] = high[N - 1]

    for n in range(N):
        assert relative_deviation(a=bin_edges[n + 1], b=high[n]) < 1e-2

    return bin_edges


def average_field_of_view_radial(X, theta_bin_edges_deg, roi_opening_deg):
    avg_shape = X.shape[1:]
    X_avg = np.zeros(avg_shape)
    total_theta_solid_angle_deg2 = 0.0
    for t in range(len(theta_bin_edges_deg) - 1):
        theta_start_deg = theta_bin_edges_deg[t]
        theta_stop_deg = theta_bin_edges_deg[t + 1]
        theta_solid_angle_deg2 = np.pi * (
            theta_stop_deg ** 2 - theta_start_deg ** 2
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


def find_common_energy_bin_edges(components):
    N_min = 1e99
    N_min_key = ""
    for key in components:
        N_key = len(components[key]["energy_bin_edges_GeV"]) - 1
        if N_key < N_min:
            N_min = N_key
            N_min_key = str(key)

    E_min = 0
    E_max = 1e99
    for key in components:
        E_start = components[key]["energy_bin_edges_GeV"][0]
        E_stop = components[key]["energy_bin_edges_GeV"][-1]
        if E_start > E_min:
            E_min = E_start
        if E_stop < E_max:
            E_max = E_stop

    bin_edges = components[N_min_key]["energy_bin_edges_GeV"]
    bins_fine = np.logical_and(bin_edges >= E_min, bin_edges <= E_max)
    num_bins = np.sum(bins_fine)

    energy_bin_edges_GeV = np.geomspace(E_min, E_max, num_bins + 1)
    return energy_bin_edges_GeV


def bin_centers(edges):
    return 0.5 * (edges[0:-1] + edges[1:])


def log10interp2d(x, y, fp, xp, yp):
    def ll(x):
        return x

    mm_f = scipy.interpolate.interp2d(
        x=ll(xp), y=ll(yp), z=fp, kind="linear"
    )
    return mm_f(ll(x), ll(y))



def integrate_energy_dispertsion(dPdMu, energy_bin_edges_GeV, target_bin_edges_GeV):
    """
    \\mu = E_{reco} / E_{true}

    Int_0^\\infty \\frac{dP}{dMu} dMu = 1
    """
    N = dPdMu.shape[0]
    assert N == dPdMu.shape[1]
    assert N == len(energy_bin_edges_GeV) - 1

    Nt = len(target_bin_edges_GeV) - 1

    energy_bin_centers_GeV = bin_centers(bin_edges=energy_bin_edges_GeV)

    R = np.zeros(shape=(Nt, Nt))
    for itEtrue in range(Nt):
        for jtEreco in range(Nt):

            tEreco_start = target_bin_edges_GeV[jtEreco]
            tEreco_stop = target_bin_edges_GeV[jtEreco + 1]

            mask = np.logical_and(
                energy_bin_centers_GeV >= tEreco_start,
                energy_bin_centers_GeV < tEreco_stop
            )

            R_I_Etrue = dPdMu[mask, iEtrue]

    return R



def _read_energy_dispersion(hdu):
    hdu
    return {
        "energy_bin_edges_GeV": merge_bin_edges_low_high(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "Mu_bin_edges" : merge_bin_edges_low_high(
            low=hdu.data["MIGRA_LO"][0, :],
            high=hdu.data["MIGRA_HI"][0, :],
        ),
        "theta_bin_edges_deg": merge_bin_edges_low_high(
            low=hdu.data["THETA_LO"][0, :], high=hdu.data["THETA_HI"][0, :],
        ),
        "dPdMu_vs_theta": hdu.data["MATRIX  "][0],
    }


def _read_effective_area(hdu):
    return {
        "energy_bin_edges_GeV": merge_bin_edges_low_high(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "theta_bin_edges_deg": merge_bin_edges_low_high(
            low=hdu.data["THETA_LO"][0, :], high=hdu.data["THETA_HI"][0, :],
        ),
        "area_vs_theta_m2": hdu.data["EFFAREA"][0],
    }


def _read_background(hdu):
    return {
        "energy_bin_edges_GeV": merge_bin_edges_low_high(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "detx_bin_edges_deg": merge_bin_edges_low_high(
            low=hdu.data["DETX_LO"][0, :], high=hdu.data["DETX_HI"][0, :],
        ),
        "dety_bin_edges_deg": merge_bin_edges_low_high(
            low=hdu.data["DETY_LO"][0, :], high=hdu.data["DETY_HI"][0, :],
        ),
        "background_vs_detx_vs_dety_per_s_per_sr_per_MeV": hdu.data["BKG"][0],
    }


def _read_point_spread_function(hdu):
    return {
        "energy_bin_edges_GeV": merge_bin_edges_low_high(
            low=hdu.data["ENERG_LO"][0, :] * 1e3,
            high=hdu.data["ENERG_HI"][0, :] * 1e3,
        ),
        "theta_bin_edges_deg": merge_bin_edges_low_high(
            low=hdu.data["THETA_LO"][0, :], high=hdu.data["THETA_HI"][0, :],
        ),
        "sigma_vs_theta_deg": hdu.data["SIGMA_1"][0],
    }


def read_instrument_response_function(path):
    irf = {}
    with astropy.io.fits.open(path) as f:
        irf["effective_area"] = _read_effective_area(hdu=f["EFFECTIVE AREA"])
        irf["energy_dispersion"] =  _read_energy_dispersion(hdu=f["ENERGY DISPERSION"])
        irf["background"] = _read_background(hdu=f["BACKGROUND"])
        irf["point_spread_function"] = _read_point_spread_function(hdu=f["POINT SPREAD FUNCTION"])

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

    irf["background"]["background_per_s_per_sr_per_MeV"] = average_field_of_view_grid(
        X=irf["background"]["background_vs_detx_vs_dety_per_s_per_sr_per_MeV"],
        detx_bin_edges_deg=irf["background"]["detx_bin_edges_deg"],
        dety_bin_edges_deg=irf["background"]["dety_bin_edges_deg"],
        roi_opening_deg=roi_opening_deg,
    )

    irf["point_spread_function"]["sigma_deg"] = average_field_of_view_radial(
        X=irf["point_spread_function"]["sigma_vs_theta_deg"],
        theta_bin_edges_deg=irf["point_spread_function"]["theta_bin_edges_deg"],
        roi_opening_deg=roi_opening_deg,
    )

    return irf


def integrate_dPdMu_to_get_probability_reco_given_true(
    dPdMu,
    dPdMu_energy_bin_edges,
    dPdMu_Mu_bin_edges,
    energy_bin_edges,
):
    dPdMu_Mu_bin_centers = bin_centers(edges=dPdMu_Mu_bin_edges)
    Mu_bin_widths = binning_utils.widths(bin_edges=dPdMu_Mu_bin_edges)
    dPdMu_energy_bin_centers = bin_centers(edges=dPdMu_energy_bin_edges)
    energy_bin_centers = bin_centers(edges=energy_bin_edges)

    N = len(dPdMu_energy_bin_edges) - 1
    Nt = len(energy_bin_edges) - 1

    R = np.nan = np.ones(shape=(Nt, Nt))

    for r in range(Nt):
        for t in range(Nt):

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
            for u in range(len(dPdMu_Mu_bin_centers)):
                R_I_Etrue += dPdMu[u, :] * Mu_mask[u] * Mu_bin_widths[u]

            assert R_I_Etrue.shape[0] == dPdMu.shape[0]

            Etrue_start = energy_bin_edges[t]
            Etrue_stop = energy_bin_edges[t + 1]
            assert Etrue_stop > Etrue_start
            Etrue_mask = np.logical_and(
                dPdMu_energy_bin_centers >= Etrue_start,
                dPdMu_energy_bin_centers < Etrue_stop
            )
            delta_Etrue = Etrue_stop - Etrue_start

            # integrate over Etrue
            # --------------------
            R[t, r] = np.sum(R_I_Etrue[Etrue_mask]) / delta_Etrue

    # normalize
    # ---------

    # ax0 -> true
    # ax1 -> reco
    probability_reco_given_true = np.nan = np.ones(shape=(Nt, Nt))

    for t in range(Nt):
        norm = np.sum(R[t, :])
        for r in range(Nt):
            probability_reco_given_true[t, r] = R[t, r] / norm

    # check
    for t in range(Nt):
        if np.sum(probability_reco_given_true[t, :]) > 0.0:
            assert 0.99 < np.sum(probability_reco_given_true[t, :]) < 1.01

    return probability_reco_given_true



def my_instrument_response_function(path, roi_opening_deg):
    irf = read_instrument_response_function(path=path)
    irf = average_instrument_response_over_field_of_view(irf=irf, roi_opening_deg=roi_opening_deg)

    energy_bin_edges_GeV = find_common_energy_bin_edges(components=irf)

    signal_area_m2 = np.interp(
        x=bin_centers(edges=energy_bin_edges_GeV),
        xp=bin_centers(edges=A["energy_bin_edges_GeV"]),
        fp=A["area_m2"],
    )

    bg_per_s_per_sr_per_GeV = np.interp(
        x=bin_centers(edges=energy_bin_edges_GeV),
        xp=bin_centers(edges=B["energy_bin_edges_GeV"]),
        fp=B["bg_per_s_per_sr_per_MeV"] * 1e3,
    )

    psf_sigma_deg = np.interp(
        x=bin_centers(edges=energy_bin_edges_GeV),
        xp=bin_centers(edges=P["energy_bin_edges_GeV"]),
        fp=P["sigma_deg"],
    )

    energy_confusion = log10interp2d(
        x=bin_centers(edges=energy_bin_edges_GeV),
        y=bin_centers(edges=energy_bin_edges_GeV),
        fp=C["confusion"],
        xp=bin_centers(edges=C["energy_bin_edges_GeV"]),
        yp=bin_centers(edges=C["energy_bin_edges_GeV"]),
    )

    # back ground rate
    bg_rate_per_s_per = np.zeros(len(energy_bin_edges_GeV) - 1)
    for e in range(len(energy_bin_edges_GeV) - 1):
        E_start_GeV = energy_bin_edges_GeV[e]
        E_stop_GeV = energy_bin_edges_GeV[e + 1]
        dE_GeV = E_stop_GeV - E_start_GeV

        psf_sigma_rad = np.deg2rad(psf_sigma_deg[e])
        psf_solid_angle_sr = np.pi * psf_sigma_rad ** 2
        bg_rate_per_s_per[e] = (
            bg_per_s_per_sr_per_GeV[e] * dE_GeV * psf_solid_angle_sr
        )

    irf = {}
    irf["energy_bin_edges_GeV"] = energy_bin_edges_GeV
    irf["signal_area_m2"] = signal_area_m2
    irf["background_per_s_per_sr_per_GeV"] = bg_per_s_per_sr_per_GeV
    irf["background_rate_per_s"] = bg_rate_per_s_per
    irf["psf_sigma_deg"] = psf_sigma_deg
    irf["energy_confusion"] = energy_confusion
    return irf
