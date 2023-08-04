#!/usr/bin/python
import argparse
import numpy as np
import plenoirf
import flux_sensitivity
import os
import sebastians_matplotlib_addons as seb
import spectral_energy_distribution_units as sed
from plenoirf.analysis import spectral_energy_distribution as sed_styles
import json_utils
import binning_utils
import cosmic_fluxes

parser = argparse.ArgumentParser(
    prog="plot_differential_sensitivity", description=("Make plots"),
)
parser.add_argument(
    "input_dir",
    metavar="INPUT_DIR",
    type=str,
    help="Path of directory with the IRF and diff. sensitivity.",
)
args = parser.parse_args()


AA = json_utils.tree.read(args.input_dir)
out_dir = os.path.join(args.input_dir, "plot")
os.makedirs(out_dir, exist_ok=True)

CONFIG = AA["config"]
IRF = AA["irf"]

energy_bin = binning_utils.Binning(bin_edges=IRF["energy_bin_edges_GeV"])
x_lim_GeV = energy_bin["decade_limits"]
y_lim_per_m2_per_s_per_GeV = np.array([1e-4, 1e-14])

y_lim_background_rate = np.array([1e-5, 1e0])
y_lim_gamma_area = np.array([1e0, 1e7])

observation_time = CONFIG["observation_time_s"]
observation_time_str = "{:.0f}s".format(observation_time)
internal_sed_style = sed_styles.PLENOIRF_SED_STYLE
COSMIC_RAYS = {"total": "black"}

# diff sens
# ---------
SED_STYLES = {
    "portal": sed_styles.PLENOIRF_SED_STYLE,
    "science": sed_styles.SCIENCE_SED_STYLE,
    "fermi": sed_styles.FERMI_SED_STYLE,
    "cta": sed_styles.CHERENKOV_TELESCOPE_ARRAY_SED_STYLE,
}

crab_flux = cosmic_fluxes.read_crab_nebula_flux_from_resources()
fermi = plenoirf.other_instruments.fermi_lat
cta = plenoirf.other_instruments.cherenkov_telescope_array_south


def fig_add_scenario_marker(
    fig, scenario_key, pos_relative=(0.0, 0.0), size=0.66
):
    ax = seb.add_axes(fig=fig, span=[0, 0, 1, 1], style=seb.AXES_BLANK)

    text_color = "black"
    if scenario_key in seb.matplotlib.colors.cnames:
        # use color
        x, y = pos_relative

        c = seb.matplotlib.colors.cnames[scenario_key]
        rgb = np.array([int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)])
        argb = 255 - rgb
        alum = np.sum(argb) / 3

        if alum < 128:
            argb = [0, 0, 0]
        else:
            argb = [255, 255, 255]

        text_color = "#{:02X}{:02X}{:02X}".format(argb[0], argb[1], argb[2])

        dx = 0.2 * size
        dy = 0.1 * size

        ax.fill_between(
            x=[x, x + dx],
            y1=[y, y],
            y2=[y + dy, y + dy],
            facecolor=scenario_key,
            transform=ax.transAxes,
        )

    # plain text
    ax.text(
        x=pos_relative[0] + dx * 0.1,
        y=pos_relative[1] + dy * 0.3,
        s=scenario_key,
        color=text_color,
        transform=ax.transAxes,
        fontsize=12 * size,
    )


# plot
# ----


# irf
# ===

fig = seb.figure(plenoirf.summary.figure.FIGURE_STYLE)
ax = seb.add_axes(fig=fig, span=plenoirf.summary.figure.AX_SPAN)
seb.ax_add_histogram(
    ax=ax,
    bin_edges=energy_bin["edges"],
    bincounts=AA["irf"]["background_rate_onregion_per_s"],
    linestyle="-",
    linecolor="black",
    linealpha=1.0,
    bincounts_upper=AA["irf"]["background_rate_onregion_per_s"]
    + AA["irf"]["background_rate_onregion_per_s_au"],
    bincounts_lower=AA["irf"]["background_rate_onregion_per_s"]
    - AA["irf"]["background_rate_onregion_per_s_au"],
    face_color="black",
    face_alpha=0.2,
    label=None,
    draw_bin_walls=False,
)
ax.set_ylabel("rate / s$^{-1}$")
ax.set_xlabel("reco. energy / GeV")
ax.set_ylim(y_lim_background_rate)
ax.loglog()
fig.savefig(os.path.join(out_dir, "irf_background_rate_onregion.jpg"))
seb.close(fig)


fig = seb.figure(plenoirf.summary.figure.FIGURE_STYLE)
ax = seb.add_axes(fig=fig, span=plenoirf.summary.figure.AX_SPAN)
ax.plot(
    energy_bin["centers"], IRF["signal_area_m2"], "+k",
)
ax.set_ylabel("area / m$^{2}$")
ax.set_xlabel("energy / GeV")
ax.set_ylim(y_lim_gamma_area)
ax.loglog()
fig.savefig(os.path.join(out_dir, "irf_signal_area.jpg"))
seb.close(fig)


fig = seb.figure(seb.FIGURE_1_1)
ax_c = seb.add_axes(fig=fig, span=[0.16, 0.16, 0.7, 0.7])
ax_cb = seb.add_axes(fig=fig, span=[0.88, 0.16, 0.02, 0.7])
_pcm_confusion = ax_c.pcolormesh(
    energy_bin["edges"],
    energy_bin["edges"],
    np.transpose(AA["irf"]["probability_reco_given_true"]),
    cmap="Greys",
    norm=seb.plt_colors.PowerNorm(gamma=0.5),
    vmin=0,
    vmax=1,
)
ax_c.grid(color="k", linestyle="-", linewidth=0.66, alpha=0.1)
seb.plt.colorbar(_pcm_confusion, cax=ax_cb, extend="max")
ax_c.set_aspect("equal")
ax_c.set_ylabel("reco. energy / GeV")
ax_c.loglog()
ax_c.set_xlabel("energy / GeV")
fig.savefig(os.path.join(out_dir, "irf_probability_reco_given_true.jpg"))
seb.close(fig)

# scenarios
# =========

for dk in flux_sensitivity.differential.SCENARIOS:
    elabel = flux_sensitivity.differential.SCENARIOS[dk]["energy_axes_label"]
    SCENARIO = AA["scenarios"][dk]

    fig = seb.figure(plenoirf.summary.figure.FIGURE_STYLE)
    fig_add_scenario_marker(fig=fig, scenario_key=dk)
    ax = seb.add_axes(fig=fig, span=plenoirf.summary.figure.AX_SPAN)
    for ck in COSMIC_RAYS:
        ck_Rt = SCENARIO["background_rate_onregion_in_scenario_per_s"]
        ck_Rt_au = SCENARIO["background_rate_onregion_in_scenario_per_s_au"]
        seb.ax_add_histogram(
            ax=ax,
            bin_edges=energy_bin["edges"],
            bincounts=ck_Rt,
            linestyle="-",
            linecolor=COSMIC_RAYS[ck],
            linealpha=1.0,
            bincounts_upper=ck_Rt + ck_Rt_au,
            bincounts_lower=ck_Rt - ck_Rt_au,
            face_color=COSMIC_RAYS[ck],
            face_alpha=0.2,
            label=None,
            draw_bin_walls=False,
        )
    ax.set_ylabel("rate / s$^{-1}$")
    ax.set_xlabel("reco. energy / GeV")
    ax.set_ylim(y_lim_background_rate)
    ax.loglog()
    fig.savefig(
        os.path.join(
            out_dir, "background_rate_in_scenario_{:s}.jpg".format(dk)
        )
    )
    seb.close(fig)

    fig = seb.figure(plenoirf.summary.figure.FIGURE_STYLE)
    fig_add_scenario_marker(fig=fig, scenario_key=dk)
    ax = seb.add_axes(fig=fig, span=plenoirf.summary.figure.AX_SPAN)
    seb.ax_add_histogram(
        ax=ax,
        bin_edges=energy_bin["edges"],
        bincounts=SCENARIO["signal_area_in_scenario_m2"],
        linestyle="-",
        linecolor="k",
        linealpha=1.0,
        bincounts_upper=SCENARIO["signal_area_in_scenario_m2"]
        + SCENARIO["signal_area_in_scenario_m2_au"],
        bincounts_lower=SCENARIO["signal_area_in_scenario_m2"]
        - SCENARIO["signal_area_in_scenario_m2_au"],
        face_color="k",
        face_alpha=0.2,
        label="Reco. energy, according to G-matrix for this scenario",
        draw_bin_walls=False,
    )
    ax.plot(
        energy_bin["centers"],
        IRF["signal_area_m2"],
        "+k",
        label="True energy.",
    )
    ax.set_ylabel("area / m$^{2}$")
    ax.set_xlabel(elabel + "energy / GeV")
    ax.set_ylim(y_lim_gamma_area)
    ax.loglog()
    ax.legend(loc="best", fontsize=6)
    fig.savefig(
        os.path.join(out_dir, "signal_area_in_scenario_{:s}.jpg".format(dk))
    )
    seb.close(fig)

    # G_matrix
    # ---------------------------
    G_matrix = SCENARIO["scenario"]["G_matrix"]
    fig = seb.figure(seb.FIGURE_1_1)
    fig_add_scenario_marker(fig=fig, scenario_key=dk)
    ax_c = seb.add_axes(fig=fig, span=[0.16, 0.16, 0.7, 0.7])
    ax_cb = seb.add_axes(fig=fig, span=[0.88, 0.16, 0.02, 0.7])
    _pcm_confusion = ax_c.pcolormesh(
        energy_bin["edges"],
        energy_bin["edges"],
        np.transpose(G_matrix),
        cmap="Greys",
        norm=seb.plt_colors.PowerNorm(gamma=0.5),
        vmin=0,
        vmax=1,
    )
    ax_c.grid(color="k", linestyle="-", linewidth=0.66, alpha=0.1)
    seb.plt.colorbar(_pcm_confusion, cax=ax_cb, extend="max")
    ax_c.set_aspect("equal")
    ax_c.set_ylabel("reco. energy / GeV")
    ax_c.loglog()
    ax_c.set_xlabel("energy / GeV")
    fig.savefig(
        os.path.join(out_dir, "G_matrix_in_scenario_{:s}.jpg".format(dk))
    )
    seb.close(fig)

    # B_matrix
    # --------
    B_matrix = SCENARIO["scenario"]["B_matrix"]
    fig = seb.figure(seb.FIGURE_1_1)
    fig_add_scenario_marker(fig=fig, scenario_key=dk)
    ax_c = seb.add_axes(fig=fig, span=[0.16, 0.16, 0.7, 0.7])
    ax_cb = seb.add_axes(fig=fig, span=[0.88, 0.16, 0.02, 0.7])
    _pcm_confusion = ax_c.pcolormesh(
        energy_bin["edges"],
        energy_bin["edges"],
        np.transpose(B_matrix),
        cmap="Greys",
        norm=seb.plt_colors.PowerNorm(gamma=0.5),
        vmin=0,
        vmax=1,
    )
    ax_c.grid(color="k", linestyle="-", linewidth=0.66, alpha=0.1)
    seb.plt.colorbar(_pcm_confusion, cax=ax_cb, extend="max")
    ax_c.set_aspect("equal")
    ax_c.set_ylabel("reco. energy / GeV")
    ax_c.loglog()
    ax_c.set_xlabel("energy / GeV")
    fig.savefig(
        os.path.join(out_dir, "B_matrix_in_scenario_{:s}.jpg".format(dk))
    )
    seb.close(fig)

    components = []

    # Crab reference fluxes
    # ---------------------
    for i in range(4):
        com = {}
        scale_factor = np.power(10.0, (-1) * i)
        com["energy"] = [np.array(crab_flux["energy"]["values"])]
        com["differential_flux"] = [
            scale_factor * np.array(crab_flux["differential_flux"]["values"])
        ]
        com[
            "label"
        ] = None  # "{:1.1e} Crab".format(scale_factor) if i == 0 else None
        com["color"] = "k"
        com["alpha"] = 0.25 / (1.0 + i)
        com["linestyle"] = "--"
        components.append(com.copy())

    # Fermi-LAT diff
    # --------------
    fermi_diff = fermi.differential_sensitivity(l=0, b=90)
    com = {}
    com["energy"] = [np.array(fermi_diff["energy"]["values"])]
    com["differential_flux"] = [
        np.array(fermi_diff["differential_flux"]["values"])
    ]
    com["label"] = fermi.LABEL + ", 10y, (l=0, b=90)"
    com["color"] = fermi.COLOR
    com["alpha"] = 1.0
    com["linestyle"] = "-"
    components.append(com)

    # CTA South 1800s
    # ---------------
    cta_diff = cta.differential_sensitivity(observation_time=observation_time)
    com = {}
    com["energy"] = [np.array(cta_diff["energy"]["values"])]
    com["differential_flux"] = [
        np.array(cta_diff["differential_flux"]["values"])
    ]
    com["label"] = cta.LABEL + ", " + observation_time_str
    com["color"] = cta.COLOR
    com["alpha"] = 1.0
    com["linestyle"] = "-"
    components.append(com)

    # My IRF
    # ------
    com = {}
    com["energy"] = []
    com["differential_flux"] = []
    com["differential_flux_au"] = []

    for ii in range(energy_bin["num"]):
        com["energy"].append(
            [energy_bin["edges"][ii], energy_bin["edges"][ii + 1]]
        )
        _dFdE_sens = SCENARIO["dVdE_per_m2_per_GeV_per_s"][ii]
        com["differential_flux"].append([_dFdE_sens, _dFdE_sens])

        _dFdE_sens_au = SCENARIO["dVdE_per_m2_per_GeV_per_s_au"][ii]
        com["differential_flux_au"].append([_dFdE_sens_au, _dFdE_sens_au])

    com["label"] = CONFIG["name"]
    com["color"] = "black"
    com["alpha"] = 1.0
    com["linestyle"] = "-"
    components.append(com)

    for sedk in SED_STYLES:
        sed_style = SED_STYLES[sedk]
        sed_style_dirname = "sed_style_" + sedk
        os.makedirs(os.path.join(out_dir, sed_style_dirname), exist_ok=True)

        fig = seb.figure(plenoirf.summary.figure.FIGURE_STYLE)
        fig_add_scenario_marker(fig=fig, scenario_key=dk)
        ax = seb.add_axes(fig=fig, span=plenoirf.summary.figure.AX_SPAN)

        for com in components:

            for ii in range(len(com["energy"])):
                _energy, _dFdE = sed.convert_units_with_style(
                    x=com["energy"][ii],
                    y=com["differential_flux"][ii],
                    input_style=internal_sed_style,
                    target_style=sed_style,
                )

                ax.plot(
                    _energy,
                    _dFdE,
                    label=com["label"] if ii == 0 else None,
                    color=com["color"],
                    alpha=com["alpha"],
                    linestyle=com["linestyle"],
                )

                if "differential_flux_au" in com:
                    _, _dFdE_au = sed.convert_units_with_style(
                        x=com["energy"][ii],
                        y=com["differential_flux_au"][ii],
                        input_style=internal_sed_style,
                        target_style=sed_style,
                    )
                    _d_ru = _dFdE_au / _dFdE
                    _dFdE_lu = _dFdE - _dFdE_au
                    _dFdE_uu = _dFdE + _dFdE_au
                    ax.fill_between(
                        x=_energy,
                        y1=_dFdE_lu,
                        y2=_dFdE_uu,
                        label=None,
                        color=com["color"],
                        alpha=com["alpha"] * 0.15,
                        linewidth=0.0,
                    )

        _x_lim, _y_lim = sed.convert_units_with_style(
            x=x_lim_GeV,
            y=y_lim_per_m2_per_s_per_GeV,
            input_style=internal_sed_style,
            target_style=sed_style,
        )

        ax.set_xlim(np.sort(_x_lim))
        ax.set_ylim(np.sort(_y_lim))
        ax.loglog()
        ax.legend(loc="best", fontsize=6)
        etype = flux_sensitivity.differential.SCENARIOS[dk][
            "energy_axes_label"
        ]
        ax.set_xlabel(
            etype + " " + sed_style["x_label"] + " /" + sed_style["x_unit"]
        )
        ax.set_ylabel(sed_style["y_label"] + " /\n " + sed_style["y_unit"])
        fig.savefig(
            os.path.join(
                out_dir,
                sed_style_dirname,
                "differential_sensitivity_{:s}.jpg".format(dk),
            )
        )
        seb.close(fig)
