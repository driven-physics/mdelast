"""Plotting helpers for simulation time series and thermal averages.

``PlotSeries`` writes diagnostic figures into a simulation's figure directory.
Plot generation is optional and is not required for running the MD engine.
"""

import os
from os import path
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

from md_driver_module import *

import logging

if len(logging.root.handlers) == 0:
    logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("Plot")


def get_splitted_axes(delta=0.0):
    """Create two adjacent axes for a broken x-axis plot."""

    gs1 = gridspec.GridSpec(1, 2)
    gs1.update(wspace=delta, hspace=0.0)  # set the spacing between axes.
    ax1 = plt.subplot(gs1[0], facecolor="w")
    ax2 = plt.subplot(gs1[1], facecolor="w")
    return ax1, ax2


def plot_splitted_x(
    *,
    x1,
    x2,
    x3,
    x4,
    x=[],
    y=[],
    labels=[],
    colors=[],
    lws=[],
    lss=[],
    x_twin=[],
    y_twin=[],
    labels_twin=[],
    colors_twin=[],
    lws_twin=[],
    lss_twin=[],
    xscale_left="linear",
    xscale_right="linear",
    yscale="linear",
    yscale_twin="linear",
    xlabel=None,
    ylabel=None,
    ylabel_twin=None,
    ymin=None,
    ymax=None,
    ymin_twin=None,
    ymax_twin=None,
    delta=0.03,
    d=0.01,
    plot_separation=True,
):
    """
    Broken axis example, where the x-axis will have a portion cut out.
    https://stackoverflow.com/questions/32185411/break-in-x-axis-of-matplotlib
    """

    assert {len(x)} == {len(var) for var in [y, labels, colors, lws, lss]}
    assert {len(x_twin)} == {
        len(var) for var in [y_twin, labels_twin, colors_twin, lws_twin, lss_twin]
    }

    lines = []
    ax1, ax2 = get_splitted_axes(delta)

    # plot the same data on both axes
    for xi, yi, label, c, lw, ls in zip(x, y, labels, colors, lws, lss):
        lines += ax1.plot(xi, yi, label=label, c=c, lw=lw, ls=ls)
        ax2.plot(xi, yi, c=c, lw=lw, ls=ls)

    if plot_separation:
        ax1.axvline(x=x2, c="grey", ls=":", lw=3.0)
        ax2.axvline(x=x3, c="grey", ls=":", lw=3.0)

    # plot the "/" separations
    kwargs = dict(transform=ax1.transAxes, color="grey", clip_on=False)
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((-d, +d), (-d, +d), **kwargs)

    ax1.set_xscale(xscale_left)
    ax2.set_xscale(xscale_right)
    ax1.set_yscale(yscale)
    ax2.set_yscale(yscale)

    # labels
    ax1.set_xlabel(xlabel)
    ax2.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    # twin right y axes
    twin = len(x_twin) > 0
    if twin:
        ax1_2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        ax2_2 = ax2.twinx()  # instantiate a second axes that shares the same x-axis

        ax1_2.tick_params(left=False, right=False, labelleft=False, labelright=False)
        ax2_2.tick_params(left=False, right=True, labelleft=False, labelright=True)
        ax1_2.spines["right"].set_visible(False)
        ax1_2.spines["left"].set_visible(False)
        ax2_2.spines["left"].set_visible(False)
        ax2_2.spines["right"].set_visible(False)
        ax1_2.set_yticklabels([])

        for xi, yi, label, c, lw, ls in zip(
            x_twin, y_twin, labels_twin, colors_twin, lws_twin, lss_twin
        ):
            lines += ax1_2.plot(xi, yi, label=label, c=c, lw=lw, ls=ls)
            ax2_2.plot(xi, yi, c=c, lw=lw, ls=ls)

        ax1_2.set_yscale(yscale_twin)
        ax2_2.set_yscale(yscale_twin)

        ax2_2.set_ylabel(ylabel_twin)

        if ymin_twin:
            ax1_2.set_ylim(ymin=ymin_twin)
            ax2_2.set_ylim(ymin=ymin_twin)
        if ymax_twin:
            ax1_2.set_ylim(ymax=ymax_twin)
            ax2_2.set_ylim(ymax=ymax_twin)

    else:
        ax1_2, ax2_2 = None, None

    # hide the spines between ax and ax2
    ax1.tick_params(left=False, right=False, labelleft=True, labelright=False)
    ax2.tick_params(left=False, right=False, labelleft=False, labelright=False)
    ax1.spines["right"].set_visible(False)
    ax2.spines["left"].set_visible(False)
    ax2.set_yticklabels([])

    # limits
    ax1.set_xlim(x1, x2)
    ax2.set_xlim(x3, x4)

    if ymin:
        ax1.set_ylim(ymin=ymin)
        ax2.set_ylim(ymin=ymin)
    if ymax:
        ax1.set_ylim(ymax=ymax)
        ax2.set_ylim(ymax=ymax)

    # plot grid
    ax1.grid()
    ax2.grid()

    # legend
    labels = [line.get_label() for line in lines]
    ax_legend = ax1_2 if twin else ax1
    ax_legend.legend(lines, labels, loc="best", frameon=True, framealpha=0.85).set_draggable(True)

    return ax1, ax2, ax1_2, ax2_2


class PlotSeries:
    """Generate plots from a :class:`series.Series` simulation record.

    The class retains the historical ``ensamble`` spelling used by the driver.
    Figures are written under ``dir_figs`` and are displayed only when
    ``show_graphs`` is true.
    """

    def __init__(
        self,
        *,
        cycles_thermalization,
        ensamble,
        temperature_thermostat,
        pressure_barostat,
        dir_figs,
        show_graphs,
        atom_names=None,
        fontsize=12,
        eam=False,
    ):
        self.cycles_thermalization = cycles_thermalization
        self.ensamble = ensamble
        self.temperature_thermostat = temperature_thermostat
        self.pressure_barostat = pressure_barostat
        self.dir_figs = dir_figs
        self.show_graphs = show_graphs
        self.atom_names = atom_names if atom_names is not None else []
        self.eam = eam
        os.makedirs(self.dir_figs, exist_ok=True)
        plt.rcParams.update({"font.size": fontsize})

    def get_x_params(self, series, plot_steps):
        if plot_steps:
            # plotting steps instead of times
            times = 1 + np.array(list(range(len(series.time))))
            suptitle = f"MD steps={len(series.time)},  "
            xlabel = "MD steps"
        else:
            times = np.array(series.time) / Units.ps
            suptitle = f"t={times[-1]:g}ps,  "
            xlabel = "$t \\ \\ \\mathrm{{[ps]}}$"

        if self.ensamble in {"NVT", "NPT"}:
            suptitle += f" T={self.temperature_thermostat:g}K"

        tmax = times[-1]
        return times, times[self.cycles_thermalization :], suptitle, xlabel, tmax

    def temperature(self, *, series, plot_steps=True):
        times, times_thermal, suptitle, xlabel, tmax = self.get_x_params(series, plot_steps)

        # ------ temperature
        label_kinetic = (
            f"$\\bar{{T}}_{{k}}={series.temperature_kin.mean:.2f}$ K,"
            f" $\\overline{{\\delta T_{{k}}}}/\\overline{{T}}_{{k}}={series.temperature_kin.std_mean:.3f}$"
        )
        data = 1.0 / np.array(series.temperature_conf_inverse.data)[self.cycles_thermalization :]
        mean = np.average(data)
        std = np.sqrt(np.average((data - mean) ** 2))
        label_conf = (
            f"$\\bar{{T}}_{{c}}={mean:.2f}$ K,"
            f" $\\overline{{\\delta T_{{c}}}}/\\overline{{T}}_{{c}}={std/abs(mean):.3f}$"
        )

        plot_splitted_x(
            x1=1.0,
            x2=tmax / 10.0,
            x3=tmax / 10.0,
            x4=tmax,
            x=[times, times_thermal],
            y=[series.temperature_kin.data, series.temperature_kin.mean_cum],
            labels=["$T_{{K}}$", "$\\overline{{T}}_{{K}}$"],
            colors=["r", "k"],
            lws=[1.0, 2.0],
            lss=["-", "--"],
            x_twin=[times_thermal],
            y_twin=[np.array(series.temperature_kin.std_cum)],
            labels_twin=["$\\delta T_{{k}}$"],
            colors_twin=["b"],
            lws_twin=[2.0],
            lss_twin=["-"],
            xscale_left="log",
            xlabel=xlabel,
            ylabel="T [K]",
            ylabel_twin="$\\delta T_{{k}}$ [K]",
            delta=0.03,
        )

        plt.savefig(path.join(self.dir_figs, f"temps_splitted.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"temps_splitted.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        plot_splitted_x(
            x1=1.0,
            x2=tmax / 10.0,
            x3=tmax / 10.0,
            x4=tmax,
            x=[times, times_thermal],
            y=[
                np.array(series.temperature_kin.data),
                np.array(series.temperature_kin.mean_cum),
            ],
            labels=["$T_{{k}}$", "$\\overline{{T_{{k}}}}$"],
            colors=["r", "k"],
            lws=[1.0, 2],
            lss=["-", "--"],
            x_twin=[times_thermal, times_thermal[1:]],
            y_twin=[
                np.abs(
                    np.array(series.temperature_kin.mean_cum) / series.temperature_kin.mean - 1.0
                ),
                (
                    np.array(series.temperature_kin.std_cum)
                    / np.abs(series.temperature_kin.mean_cum)
                )[1:],
            ],
            labels_twin=["error", "$\\delta T_{{k}} / \\overline{{T_{{k}}}}$"],
            colors_twin=["g", "b"],
            lws_twin=[1.0, 2.0],
            lss_twin=[":", "-"],
            yscale_twin="log",
            ymin_twin=1e-6,
            ymax_twin=1.0,
            xscale_left="log",
            xlabel=xlabel,
            ylabel="$T_K$ [K]",
            delta=0.03,
        )

        plt.savefig(path.join(self.dir_figs, f"temperature_err_splitted.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"temperature_err_splitted.pdf"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"temp_kinetic_splitted.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"temp_kinetic_splitted.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        tc = 1.0 / np.asarray(series.temperature_conf_inverse.data)
        tc_mean_cum = 1.0 / np.asarray(series.temperature_conf_inverse.mean_cum)
        tc_mean = 1.0 / series.temperature_conf_inverse.mean
        tc_thermal = tc[self.cycles_thermalization :]
        if tc_thermal.size:
            tc_cum_sum = np.cumsum(tc_thermal)
            tc_cum_sq = np.cumsum(tc_thermal**2)
            tc_counts = np.arange(1, tc_thermal.size + 1)
            tc_std_cum = np.sqrt(
                np.maximum(tc_cum_sq / tc_counts - (tc_cum_sum / tc_counts) ** 2, 0.0)
            )
        else:
            tc_std_cum = np.array([])

        plot_splitted_x(
            x1=1.0,
            x2=tmax / 10.0,
            x3=tmax / 10.0,
            x4=tmax,
            x=[times, times_thermal],
            y=[tc, tc_mean_cum],
            labels=["$T_{{c}}$", "$\\overline{{T_{{c}}}}$"],
            colors=["r", "k"],
            lws=[1.0, 2],
            lss=["-", "--"],
            x_twin=[times_thermal, times_thermal[1:]],
            y_twin=[np.abs(tc_mean_cum / tc_mean - 1.0), (tc_std_cum / np.abs(tc_mean_cum))[1:]],
            labels_twin=["error", "$\\delta T_{{c}} / \\overline{{T_{{c}}}}$"],
            colors_twin=["g", "b"],
            lws_twin=[1.0, 2.0],
            lss_twin=[":", "-"],
            yscale_twin="log",
            ymin_twin=1e-6,
            ymax_twin=1.0,
            xscale_left="log",
            xlabel=xlabel,
            ylabel="$T_c$ [K]",
            delta=0.03,
        )

        plt.savefig(path.join(self.dir_figs, f"temp_config_splitted.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"temp_config_splitted.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        for xscale in ("linear", "log"):
            if self.ensamble in {"NVT", "NPT"}:
                tk_tnh_err = abs(series.temperature_kin.mean / self.temperature_thermostat - 1.0)
                plt.suptitle(f"  $|\\overline{{T_{{k}}}}/T_{{Thermostat}}-1|=${tk_tnh_err:g}")
                plt.axhline(
                    y=self.temperature_thermostat,
                    c="k",
                    ls="--",
                    lw=2,
                    label=f"$T_{{NH}}={self.temperature_thermostat:g}$K",
                )

            times_plot = np.asarray(times)
            times_thermal_plot = np.asarray(times_thermal)

            tk = np.asarray(series.temperature_kin.data)
            tk_mean_cum = np.asarray(series.temperature_kin.mean_cum)
            tc = 1.0 / np.asarray(series.temperature_conf_inverse.data)
            tc_mean_cum = 1.0 / np.asarray(series.temperature_conf_inverse.mean_cum)

            if xscale == "log":
                # log-x cannot show t=0; drop non-positive points
                mask = times_plot > 0
                times_plot = times_plot[mask]
                tk = tk[mask]
                tc = tc[mask]

                mask_th = times_thermal_plot > 0
                times_thermal_plot = times_thermal_plot[mask_th]
                tk_mean_cum = tk_mean_cum[mask_th]
                tc_mean_cum = tc_mean_cum[mask_th]

            plt.plot(times_plot, tk, lw=0.5, c="b", ls="-")
            plt.plot(
                times_thermal_plot,
                tk_mean_cum,
                c="b",
                lw=2,
                ls="--",
                label=r"$\overline{{T}}_{{K}}(t)$",
            )
            plt.axhline(y=series.temperature_kin.mean, c="b", ls="--", label=label_kinetic)

            plt.plot(times_plot, tc, lw=0.5, c="r", ls="-")
            plt.plot(
                times_thermal_plot, tc_mean_cum, c="r", ls="--", label=r"$\overline{{T}}_{{c}}(t)$"
            )
            plt.axhline(
                y=1.0 / series.temperature_conf_inverse.mean, c="r", ls="--", label=label_conf
            )

            plt.xlabel(xlabel)
            plt.ylabel("T [K]")
            plt.xscale(xscale)
            plt.legend(framealpha=0.95).set_draggable(True)
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.grid()
            plt.savefig(path.join(self.dir_figs, f"temps_{xscale}.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"temps_{xscale}.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

        # temperature gaussians
        plt.hist(
            series.temperature_kin.data[self.cycles_thermalization :],
            bins=100,
            density=True,
            facecolor="b",
            alpha=0.5,
            label=label_kinetic,
        )
        plt.hist(
            1.0 / np.array(series.temperature_conf_inverse.data)[self.cycles_thermalization :],
            bins=100,
            density=True,
            facecolor="r",
            alpha=0.5,
            label=label_conf,
        )
        plt.plot(
            series.temperature_kin.data_gauss,
            series.temperature_kin.gauss,
            color="k",
            ls="--",
            lw=3,
        )

        x = np.linspace(mean - 4 * std, mean + 4 * std, 1000)
        gaussian = (
            lambda x, mu, sigma: 1.0
            / (2.0 * np.pi * sigma**2.0) ** 0.5
            * np.exp(-((x - mu) ** 2.0) / (2.0 * sigma**2))
        )

        plt.plot(x, gaussian(x, mean, std), color="b", ls="--", lw=3)

        if self.ensamble in {"NVT", "NPT"}:
            plt.axvline(
                x=self.temperature_thermostat,
                c="lime",
                ls="--",
                lw=3,
                label=f"$T_{{Nose-Hoover}}={self.temperature_thermostat:g}$K",
            )

        plt.grid()
        plt.legend(framealpha=0.95, fontsize=10)
        plt.xlabel("$T \\ \\mathrm{{[K]}}$")
        plt.ylabel("$f(T) \\ \\mathrm{{[1/K]}}$")
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, "temps_hist.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, "temps_hist.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        return

    def pressure(self, *, series, plot_steps=True):
        times, times_thermal, suptitle, xlabel, tmax = self.get_x_params(series, plot_steps)

        label = f"$\\overline{{P}}={series.pressure.mean/ Units.GPa:g}$ GPa"
        for xscale in {"linear"}:  # {"log", "linear"}:
            plt.plot(times, np.array(series.pressure.data) / Units.GPa, lw=1, c="b", label="$P(t)$")
            plt.plot(
                times_thermal,
                np.array(series.pressure.mean_cum) / Units.GPa,
                ls="--",
                c="lime",
                lw=2,
                label=r"$\overline{{P}}(t)$",
            )

            plt.axhline(
                y=series.pressure.mean / Units.GPa,
                c="b",
                ls="--",
                lw=2,
                label=f"$\\overline{{P}}={series.pressure.mean/ Units.GPa:g}$ GPa",
            )
            if self.ensamble == "NPT":
                plt.axhline(
                    y=self.pressure_barostat / Units.GPa,
                    c="lime",
                    ls="--",
                    lw=3,
                    label=f"$P_{{Barostat}}={self.pressure_barostat/Units.GPa:g}$ GPa",
                )
                pressure_baro_err = abs(series.pressure.mean / self.pressure_barostat - 1.0)
                plt.suptitle(f"$|\\overline{{P}}/P_{{Barostat}}-1|=${pressure_baro_err:g}")

            plt.xscale(xscale)
            plt.xlabel(xlabel)
            plt.ylabel("$P \\ \\  \\mathrm{{[GPa]}}$")
            plt.legend(framealpha=0.95)

            plt.grid()
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"pressure_{xscale}.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"pressure_{xscale}.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

        # pressure gaussians
        label_pressure = (
            f"$\\bar{{P}}={series.pressure.mean/ Units.GPa:g}GPa , "
            f"\\ \\overline{{\\delta P}}/|\\overline{{P}}|={series.pressure.std_mean:g}$"
        )
        plt.hist(
            np.array(series.pressure.data[self.cycles_thermalization :]) / Units.GPa,
            bins=100,
            density=True,
            facecolor="b",
            alpha=0.5,
        )
        plt.plot(
            np.array(series.pressure.data_gauss) / Units.GPa,
            series.pressure.gauss * Units.GPa,
            color="k",
            ls="--",
            lw=3,
            label=label_pressure,
        )

        plt.grid()
        plt.legend(framealpha=0.95)
        plt.xlabel("$P \\ \\mathrm{{[GPa]}}$")
        plt.ylabel("$f(P) \\ \\mathrm{{[1/GPa]}}$")
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, "pressure_hist.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, "pressure_hist.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        # --- instantaneous pressure tensor (Voigt: xx, yy, zz, yz, xz, xy)
        if hasattr(series, "stress_tensor") and series.stress_tensor is not None:
            cname = ["xx", "yy", "zz", "yz", "xz", "xy"]
            colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
            plot_splitted_x(
                x1=1.0,
                x2=tmax / 10.0,
                x3=tmax / 10.0,
                x4=tmax,
                x=[times] * 6,
                y=[series.stress_tensor.get_time_array()[i, :] / Units.GPa for i in range(6)],
                labels=cname,
                colors=colors,
                lws=[0.6] * 6,
                lss=["--"] * 6,
                xscale_left="log",
                xlabel=xlabel,
                ylabel="$P_{{\\alpha\\beta}}(t)$ [GPa]",
                delta=0.03,
            )

            plt.savefig(
                path.join(self.dir_figs, f"pressure_tensor_splitted.png"), bbox_inches="tight"
            )
            plt.savefig(
                path.join(self.dir_figs, f"pressure_tensor_splitted.pdf"), bbox_inches="tight"
            )
            if self.show_graphs:
                plt.show()
            plt.close()

            colors_1 = ["r", "k", "b", "g", "y", "m"]
            for i in range(6):
                plt.hist(
                    series.stress_tensor.get_time_array()[i, self.cycles_thermalization :]
                    / Units.GPa,
                    bins=50,
                    density=True,
                    color=colors_1[i],
                    alpha=0.3,
                )
                plt.plot(
                    np.array(series.stress_tensor.data[i].data_gauss) / Units.GPa,
                    series.stress_tensor.data[i].gauss * Units.GPa,
                    ls="--",
                    lw=1.5,
                    color=colors_1[i],
                    label=cname[i],
                )

            plt.grid()
            plt.legend(framealpha=0.85).set_draggable(True)
            plt.xlabel("$P_{{\\alpha \\beta}} \\ \\mathrm{{[GPa]}}$")
            plt.ylabel("$f(P_{{\\alpha \\beta}}) \\ \\mathrm{{[1/GPa]}}$")
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"pressure_tensor_hist.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"pressure_tensor_hist.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

    def elasticity(self, *, series, plot_steps=True):
        if self.ensamble not in {"NVT", "NPT"}:
            return
        times, times_thermal, suptitle, xlabel, tmax = self.get_x_params(series, plot_steps)

        # bulk modulus
        if self.ensamble == "NPT":
            plt.plot(
                times_thermal,
                np.abs(
                    np.array(series.bulk_modulus.data) / series.thermal_results.bulk_modulus - 1.0
                ),
                lw=1,
                c="b",
            )
            plt.xlabel(xlabel)
            plt.ylabel("$| \\overline{{K}}(t)/\\overline{{K}}-1|$")
            plt.grid()
            plt.yscale("log")
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"err_bulk_modulus.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"err_bulk_modulus.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

            plt.plot(
                times_thermal[len(times_thermal) // 2 :],
                np.array(series.bulk_modulus.data[len(times_thermal) // 2 :]) / Units.GPa,
                lw=1,
                c="b",
            )
            plt.xlabel(xlabel)
            plt.ylabel("$K \\ \\mathrm{{[GPa]}}$")
            plt.grid()
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"bulk_modulus_ylin.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"bulk_modulus_ylin.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

            plt.plot(
                times_thermal - times_thermal[0],
                np.array(series.bulk_modulus.data) / Units.GPa,
                lw=1,
                c="b",
            )
            plt.xlabel(xlabel)
            plt.yscale("log")
            plt.xscale("log")
            plt.ylabel("$K \\ \\mathrm{{[GPa]}}$")
            plt.grid()
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"bulk_modulus_ylog.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"bulk_modulus_ylog.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

        cname = ["xx", "yy", "zz", "yz", "xz", "xy"]
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
        # ---- born stress tensor
        plot_splitted_x(
            x1=1.0,
            x2=tmax / 10.0,
            x3=tmax / 10.0,
            x4=tmax,
            x=[times] * 6,
            y=[series.born_stress_tensor.get_time_array()[i, :] / Units.GPa for i in range(6)],
            labels=cname,
            colors=colors,
            lws=[0.6] * 6,
            lss=["--"] * 6,
            xscale_left="log",
            xlabel=xlabel,
            ylabel="$\\sigma^{{B}}_{{\\alpha\\beta}}(t)$ [GPa]",
            delta=0.03,
        )

        plt.savefig(
            path.join(self.dir_figs, f"born_stress_tensor_splitted.png"), bbox_inches="tight"
        )
        plt.savefig(
            path.join(self.dir_figs, f"born_stress_tensor_splitted.pdf"), bbox_inches="tight"
        )
        if self.show_graphs:
            plt.show()
        plt.close()

        # born stress gaussians
        colors_1 = ["r", "k", "b", "g", "y", "m"]
        for i in range(6):
            plt.hist(
                series.born_stress_tensor.get_time_array()[i, self.cycles_thermalization :]
                / Units.GPa,
                bins=50,
                density=True,
                color=colors_1[i],
                alpha=0.3,
            )
            plt.plot(
                np.array(series.born_stress_tensor.data[i].data_gauss) / Units.GPa,
                series.born_stress_tensor.data[i].gauss * Units.GPa,
                ls="--",
                lw=1.5,
                color=colors_1[i],
                label=cname[i],
            )

        plt.grid()
        plt.legend(framealpha=0.85).set_draggable(True)
        plt.xlabel("$\\sigma^{{B}}_{{\\alpha \\beta}} \\ \\mathrm{{[GPa]}}$")
        plt.ylabel("$f( \\sigma^{{B}}_{{\\alpha \\beta}}) \\ \\mathrm{{[1/GPa]}}$")
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, f"born_stress_tensor_hist.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"born_stress_tensor_hist.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        # TOTAL ELASTIC CONSTANTS

        plt.plot(
            np.array(series.elasticity_tensor.data[0][0].mean_cum) / Units.GPa, label="$C_{{11}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor.data[1][1].mean_cum) / Units.GPa, label="$C_{{22}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor.data[2][2].mean_cum) / Units.GPa, label="$C_{{33}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor_11.mean_cum) / Units.GPa,
            label="$\\frac{{1}}{{3}}(C_{{11}}+C_{{22}}+C_{{33}})$",
        )
        plt.ylabel("$\\overline{{C}}(t)$")
        plt.grid()
        plt.legend()
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, f"elasticity_tensor_11.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"elasticity_tensor_11.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        plt.plot(
            np.array(series.elasticity_tensor.data[0][1].mean_cum) / Units.GPa, label="$C_{{12}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor.data[0][2].mean_cum) / Units.GPa, label="$C_{{13}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor.data[1][2].mean_cum) / Units.GPa, label="$C_{{23}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor_12.mean_cum) / Units.GPa,
            label="$\\frac{{1}}{{3}}(C_{{12}}+C_{{13}}+C_{{23}})$",
        )
        plt.ylabel("$\\overline{{C}}(t)$")
        plt.grid()
        plt.legend()
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, f"elasticity_tensor_12.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"elasticity_tensor_12.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        plt.plot(
            np.array(series.elasticity_tensor.data[3][3].mean_cum) / Units.GPa, label="$C_{{44}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor.data[4][4].mean_cum) / Units.GPa, label="$C_{{55}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor.data[5][5].mean_cum) / Units.GPa, label="$C_{{66}}$"
        )
        plt.plot(
            np.array(series.elasticity_tensor_44.mean_cum) / Units.GPa,
            label="$\\frac{{1}}{{3}}(C_{{44}}+C_{{55}}+C_{{66}})$",
        )
        plt.ylabel("$\\overline{{C}}(t)$")
        plt.grid()
        plt.legend()
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, f"elasticity_tensor_44.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"elasticity_tensor_44.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        for c_type in ["11", "12", "44"]:
            C = getattr(series, f"elasticity_tensor_{c_type}")
            plot_splitted_x(
                x1=1.0,
                x2=tmax / 10.0,
                x3=tmax / 10.0,
                x4=tmax,
                x=[times, times_thermal],
                y=[
                    np.array(C.data) / Units.GPa,
                    np.array(C.mean_cum) / Units.GPa,
                ],
                labels=[f"$C_{{{c_type}}}$", f"$\\overline{{C_{{{c_type}}}}}$"],
                colors=["r", "k"],
                lws=[1.0, 2],
                lss=["-", "--"],
                x_twin=[times_thermal, times_thermal[1:]],
                y_twin=[
                    np.abs(np.array(C.mean_cum) / C.mean - 1.0),
                    (np.array(C.std_cum) / np.abs(C.mean_cum))[1:],
                ],
                labels_twin=["error", f"$\\delta C_{{{c_type}}} / \\overline{{C_{{{c_type}}}}}$"],
                colors_twin=["g", "b"],
                lws_twin=[1.0, 2.0],
                lss_twin=[":", "-"],
                yscale_twin="log",
                ymin_twin=1e-8,
                ymax_twin=1e-1,
                xscale_left="log",
                xlabel=xlabel,
                ylabel=f"$C_{{{c_type}}} \\ \\mathrm{{[GPa]}}$",
                delta=0.03,
            )
            plt.savefig(
                path.join(self.dir_figs, f"elasticity_tensor_{c_type}_splitted.png"),
                bbox_inches="tight",
            )
            plt.savefig(
                path.join(self.dir_figs, f"elasticity_tensor_{c_type}_splitted.pdf"),
                bbox_inches="tight",
            )
            if self.show_graphs:
                plt.show()
            plt.close()

        # BORN ELASTIC CONSTANTS
        for c_type in ["11", "12", "44"]:
            C = getattr(series, f"born_elasticity_tensor_{c_type}")
            plot_splitted_x(
                x1=1.0,
                x2=tmax / 10.0,
                x3=tmax / 10.0,
                x4=tmax,
                x=[times, times_thermal],
                y=[
                    np.array(C.data) / Units.GPa,
                    np.array(C.mean_cum) / Units.GPa,
                ],
                labels=[f"$C^{{B}}_{{{c_type}}}$", f"$\\overline{{C^{{B}}_{{{c_type}}}}}$"],
                colors=["r", "k"],
                lws=[1.0, 2],
                lss=["-", "--"],
                x_twin=[times_thermal, times_thermal[1:]],
                y_twin=[
                    np.abs(np.array(C.mean_cum) / C.mean - 1.0),
                    (np.array(C.std_cum) / np.abs(C.mean_cum))[1:],
                ],
                labels_twin=["error", f"$\\delta C_{{{c_type}}} / \\overline{{C_{{{c_type}}}}}$"],
                colors_twin=["g", "b"],
                lws_twin=[1.0, 2.0],
                lss_twin=[":", "-"],
                yscale_twin="log",
                ymin_twin=1e-8,
                ymax_twin=1e-1,
                xscale_left="log",
                xlabel=xlabel,
                ylabel=f"$C^{{B}}_{{{c_type}}} \\ \\mathrm{{[GPa]}}$",
                delta=0.03,
            )
            plt.savefig(
                path.join(self.dir_figs, f"born_elasticity_tensor_{c_type}_splitted.png"),
                bbox_inches="tight",
            )
            plt.savefig(
                path.join(self.dir_figs, f"born_elasticity_tensor_{c_type}_splitted.pdf"),
                bbox_inches="tight",
            )
            if self.show_graphs:
                plt.show()
            plt.close()

        # born gaussians
        for c_type in ["11", "12", "44"]:
            C = getattr(series, f"born_elasticity_tensor_{c_type}")
            label = (
                f"$\\overline{{C^{{B}}_{{{c_type}}}}} = {C.mean/Units.GPa:g}GPa , "
                f"\\ \\overline{{\\delta C^{{B}}}}/|\\overline{{C^{{B}}}}|="
                f"{C.std_mean:g}$"
            )

            plt.hist(
                np.array(C.data[self.cycles_thermalization :]) / Units.GPa,
                bins=100,
                density=True,
                alpha=0.5,
                facecolor="b",
            )
            plt.plot(
                np.array(C.data_gauss) / Units.GPa, C.gauss * Units.GPa, ls="--", lw=2.0, color="k"
            )
            plt.grid()
            plt.xlabel(f"$C^{{B}}_{{{c_type}}} \\ \\mathrm{{[GPa]}}$")
            plt.ylabel(f"$f(C^{{B}}_{{{c_type}}}) \\ \\mathrm{{[1/GPa]}}$")
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(
                path.join(self.dir_figs, f"born_elasticity_tensor_diag_{c_type}_hist.png"),
                bbox_inches="tight",
            )
            plt.savefig(
                path.join(self.dir_figs, f"born_elasticity_tensor_diag_{c_type}_hist.pdf"),
                bbox_inches="tight",
            )
            if self.show_graphs:
                plt.show()
            plt.close()

        # elasticity stress fluctuations term

    def volume_fluctuations(self, *, series, plot_steps=True):
        if self.ensamble != "NPT":
            return

        times, times_thermal, suptitle, xlabel, tmax = self.get_x_params(series, plot_steps)

        # Volume fluctuations
        label = (
            f"$\\overline{{V/V_{{0}}}}={series.V_V0.mean:g}, "
            f"\\ \\overline{{\\delta V}}/|\\overline{{V}}|={series.V_V0.std_mean:g}$"
        )

        for xscale in {"log", "linear"}:
            plt.plot(times, series.V_V0.data, lw=1, c="b")
            plt.plot(
                times_thermal, np.array(series.V_V0.mean_cum), c="r", ls="--", lw=2, label=label
            )

            plt.xlabel(xlabel)
            plt.xscale(xscale)
            plt.ylabel("$V(t)/V_{{0}}$")
            plt.legend()
            plt.grid()
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"vv0_{xscale}.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"vv0_{xscale}.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

        # volume gaussians
        plt.hist(series.V_V0.data[self.cycles_thermalization :], bins=100, density=True)
        plt.plot(series.V_V0.data_gauss, series.V_V0.gauss, ls="--", lw=3, label=label)

        plt.grid()
        plt.legend(framealpha=0.95)
        plt.xlabel("$V/V_{{0}}$")
        plt.ylabel("$f(V/V_{{0}})$")
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.suptitle(suptitle)
        plt.savefig(path.join(self.dir_figs, "vv0_hist.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, "vv0_hist.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        # Length

        plot_splitted_x(
            x1=1.0,
            x2=tmax / 10.0,
            x3=tmax / 10.0,
            x4=tmax,
            x=[times, times_thermal],
            y=[
                np.array(series.L_L0.data),
                np.array(series.L_L0.mean_cum),
            ],
            labels=["$L/L_{{0}}$", "$\\overline{{L}}/L_{{0}}$"],
            colors=["r", "k"],
            lws=[1.0, 2],
            lss=["-", "--"],
            x_twin=[times_thermal, times_thermal[1:]],
            y_twin=[
                np.abs(np.array(series.L_L0.mean_cum) / series.L_L0.mean - 1.0),
                (np.array(series.L_L0.std_cum) / np.abs(series.L_L0.mean_cum))[1:],
            ],
            labels_twin=["error", "$\\delta L / \\overline{{L}}$"],
            colors_twin=["g", "b"],
            lws_twin=[1.0, 2.0],
            lss_twin=[":", "-"],
            yscale_twin="log",
            ymin_twin=1e-8,
            ymax_twin=1e-2,
            xscale_left="log",
            xlabel=xlabel,
            ylabel="$L/L_{{0}}$",
            delta=0.03,
        )

        plt.savefig(path.join(self.dir_figs, f"ll0_err_splitted.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"ll0_err_splitted.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        label = (
            f"$\\overline{{L/L_{{0}}}}={series.L_L0.mean:g}, "
            f"\\ \\overline{{\\delta L}}/|\\overline{{L}}|={series.L_L0.std_mean:g}$"
        )
        for xscale in {"log", "linear"}:
            plt.plot(times, series.L_L0.data, lw=1, c="b")
            plt.plot(
                times_thermal, np.array(series.L_L0.mean_cum), c="r", ls="--", lw=2, label=label
            )
            plt.xlabel(xlabel)
            plt.xscale(xscale)
            plt.ylabel("$L(t)/L_{0}$")
            plt.legend()
            plt.grid()
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"ll0_{xscale}.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"ll0_{xscale}.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

        # Length uncertainty

        # Absolute length

        # length gaussians
        plt.hist(
            series.L_L0.data[self.cycles_thermalization :],
            bins=100,
            density=True,
            facecolor="b",
            alpha=0.5,
        )
        plt.plot(series.L_L0.data_gauss, series.L_L0.gauss, ls="--", lw=3, color="k", label=label)

        plt.grid()
        plt.xlabel("$L/L_{{0}}$")
        plt.ylabel("$f(L/L_{{0}})$")
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, "ll0_hist.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, "ll0_hist.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        # thermal expnsion coeff
        plot_splitted_x(
            x1=1.0,
            x2=tmax / 10.0,
            x3=tmax / 10.0,
            x4=tmax,
            x=[times, times_thermal],
            y=[
                np.array(series.volume_thermal_expansion.data),
                np.array(series.volume_thermal_expansion.mean_cum),
            ],
            labels=["$\\alpha_{P}$", "$\\overline{{\\alpha_{P}}}$"],
            colors=["r", "k"],
            lws=[1.0, 2],
            lss=["-", "--"],
            x_twin=[times_thermal, times_thermal[1:]],
            y_twin=[
                np.abs(
                    np.array(series.volume_thermal_expansion.mean_cum)
                    / series.volume_thermal_expansion.mean
                    - 1.0
                ),
                (
                    np.array(series.volume_thermal_expansion.std_cum)
                    / np.abs(series.volume_thermal_expansion.mean_cum)
                )[1:],
            ],
            labels_twin=["error", "$\\delta \\alpha_{P}(T)$"],
            colors_twin=["g", "b"],
            lws_twin=[1.0, 2.0],
            lss_twin=[":", "-"],
            yscale_twin="log",
            xscale_left="log",
            xlabel=xlabel,
            ylabel="$\\alpha_{P}$",
            delta=0.03,
        )

        plt.savefig(
            path.join(self.dir_figs, f"volume_thermal_expansion_splitted.png"), bbox_inches="tight"
        )
        plt.savefig(
            path.join(self.dir_figs, f"volume_thermal_expansion_splitted.pdf"), bbox_inches="tight"
        )
        if self.show_graphs:
            plt.show()
        plt.close()

        label = f"$\\overline{{\\alpha_{{P}}}}={series.volume_thermal_expansion.mean:g}$"
        for xscale in {"log", "linear"}:
            plt.plot(times, series.volume_thermal_expansion.data, lw=1, c="b")
            plt.plot(
                times_thermal,
                np.array(series.volume_thermal_expansion.mean_cum),
                c="r",
                ls="--",
                lw=2,
                label=label,
            )

            plt.xlabel(xlabel)
            plt.xscale(xscale)
            plt.ylabel("$\\alpha_{P}$")
            plt.legend()
            plt.grid()
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, f"alphap_{xscale}.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"alphap_{xscale}.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

        # compressibility
        plt.plot(times_thermal, np.array(series.compressibility.data) * Units.GPa, lw=1, c="b")
        plt.xlabel(xlabel)
        plt.ylabel("$\\beta_{{T}} \\ \\mathrm{{[1/GPa]}}$")
        plt.grid()
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, f"compressibility_ylin.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, f"compressibility_ylin.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

    def conservation(self, *, series, plot_steps=True):
        times, times_thermal, suptitle, xlabel, tmax = self.get_x_params(series, plot_steps)

        for xscale in {"log", "linear"}:
            plt.plot(times, np.array(series.etot_cons_err.data) / 1e-5, "-")

            plt.grid()
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.xlabel(xlabel)
            plt.xscale(xscale)
            plt.ylabel("$|E_{cons}/E_{cons}(0)-1|  \\ \\ [10^{{-5}}]$ ")
            plt.savefig(path.join(self.dir_figs, f"conservation_{xscale}.png"), bbox_inches="tight")
            plt.savefig(path.join(self.dir_figs, f"conservation_{xscale}.pdf"), bbox_inches="tight")
            if self.show_graphs:
                plt.show()
            plt.close()

    def energies_cv(self, *, series, plot_steps=True):
        times, times_thermal, suptitle, xlabel, tmax = self.get_x_params(series, plot_steps)
        times_arr = np.asarray(times)

        epot = np.array(series.epot_atom.data) / Units.ev
        ekin = np.array(series.ekin_atom.data) / Units.ev
        etot = np.array(series.etot_atom.data) / Units.ev

        y_series = [epot, ekin, etot]
        labels = [
            f"$\\overline{{E}}_{{p}}={series.epot_atom.mean/Units.ev:g}eV, \\ \\overline{{\\delta E}}_{{p}}={series.epot_atom.std/Units.ev:g}$eV",
            f"$\\overline{{E}}_{{k}}={series.ekin_atom.mean/Units.ev:g}eV, \\ \\overline{{\\delta E}}_{{k}}={series.ekin_atom.std/Units.ev:g}$eV",
            f"$\\overline{{E}}_{{tot}}={series.etot_atom.mean/Units.ev:g}eV, \\ \\overline{{\\delta E}}_{{tot}}={series.etot_atom.std/Units.ev:g}$eV",
        ]
        colors = ["r", "b", "g"]
        lss = ["-", "-", "-"]
        mean_hlines = [
            (series.epot_atom.mean / Units.ev, "r"),
            (series.ekin_atom.mean / Units.ev, "b"),
        ]

        if self.eam:
            epot_eam = np.array(series.epot_eam_atom.data) / Units.ev
            epot_pair = epot - epot_eam
            y_series = [epot, epot_pair, epot_eam, ekin, etot]
            labels = [
                f"$\\overline{{E}}_{{p}}={series.epot_atom.mean/Units.ev:g}eV, \\ \\overline{{\\delta E}}_{{p}}={series.epot_atom.std/Units.ev:g}$eV",
                f"$\\overline{{E}}_{{pair}}={np.mean(epot_pair):g}eV, \\ \\overline{{\\delta E}}_{{pair}}={np.std(epot_pair):g}$eV",
                f"$\\overline{{E}}_{{eam}}={np.mean(epot_eam):g}eV, \\ \\overline{{\\delta E}}_{{eam}}={np.std(epot_eam):g}$eV",
                f"$\\overline{{E}}_{{k}}={series.ekin_atom.mean/Units.ev:g}eV, \\ \\overline{{\\delta E}}_{{k}}={series.ekin_atom.std/Units.ev:g}$eV",
                f"$\\overline{{E}}_{{tot}}={series.etot_atom.mean/Units.ev:g}eV, \\ \\overline{{\\delta E}}_{{tot}}={series.etot_atom.std/Units.ev:g}$eV",
            ]
            colors = ["r", "orange", "purple", "b", "g"]
            lss = ["-", "--", "--", "-", "-"]
            mean_hlines = [
                (series.epot_atom.mean / Units.ev, "r"),
                (series.ekin_atom.mean / Units.ev, "b"),
            ]

        if self.ensamble == "NVE":
            title = (
                f"$|\\overline{{\\delta E}}_{{pot}}/\\overline{{\\delta E}}_{{kin}}-1|="
                f"{series.thermal_results.epot_ekin_std_err:g}, "
                f"$|c_{{V,E_{{kin}}}}/c_{{V,E_{{pot}}}}-1|={series.thermal_results.cv_epot_ekin_err:g}$\n"
                f"$c_{{V,E_{{kin}}}}={series.thermal_results.cv_ekin*Units.Navogadro/Units.joule:g}$ [J/(mol K)], "
                f"($c_{{V,id}}={series.thermal_results.cv_id*Units.Navogadro/Units.joule:g}$)"
            )
        elif self.ensamble == "NVT":
            title = (
                f"$|\\overline{{\\delta E}}_{{kin}}^{{2}}/1.5 N (k_{{B}} T)^{{2}}-1|="
                f"{series.thermal_results.ekin_atom_std_err:g}, "
                f"$|c_{{V,E_{{tot}}}}/c_{{V,E_{{pot}}}}-1|={series.thermal_results.cv_etot_epot_err:g}$\n"
                f"$c_{{V,E_{{tot}}}}={series.thermal_results.cv_etot*Units.Navogadro/Units.joule:g}$ [J/(mol K)], "
                f"($c_{{V,id}}={series.thermal_results.cv_id*Units.Navogadro/Units.joule:g}$)"
            )
        elif self.ensamble == "NPT":
            title = ""
        else:
            logger.fatal("bad ensamble")
            quit()

        if not self.eam:
            title += (
                f"$|\\overline{{E}}_{{tot}}/\\overline{{E}}_{{tot, g(r)}}-1|="
                f"{series.thermal_results.etot_atom_g_err:g}"
            )

        def _plot_energies_on_ax(
            ax, times_plot, y_plot, *, draw_mean_hlines=False, use_legend=True
        ):
            for yi, label, c, ls in zip(y_plot, labels, colors, lss):
                ax.plot(times_plot, yi, c=c, ls=ls, lw=0.8, label=label)
            if draw_mean_hlines:
                for ymean, c in mean_hlines:
                    ax.axhline(ymean, c=c, ls="--")
            ax.grid(True, alpha=0.3)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("E/atom [eV]")
            if use_legend:
                ax.legend(framealpha=0.95, fontsize=8)

        # --- linear (default)
        fig, ax = plt.subplots()
        _plot_energies_on_ax(ax, times_arr, y_series, draw_mean_hlines=True)
        if title:
            fig.suptitle(title, fontsize=10.0)
        fig.autofmt_xdate()
        ax.autoscale(enable=True, axis="both", tight=True)
        fig.savefig(path.join(self.dir_figs, "energies.png"), bbox_inches="tight")
        fig.savefig(path.join(self.dir_figs, "energies.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close(fig)

        # --- splitted x-axis (early / late zoom)
        plot_splitted_x(
            x1=1.0,
            x2=tmax / 10.0,
            x3=tmax / 10.0,
            x4=tmax,
            x=[times_arr] * len(y_series),
            y=y_series,
            labels=labels,
            colors=colors,
            lws=[0.8] * len(y_series),
            lss=lss,
            xscale_left="log",
            xlabel=xlabel,
            ylabel="E/atom [eV]",
            delta=0.03,
        )
        plt.savefig(path.join(self.dir_figs, "energies_splitted.png"), bbox_inches="tight")
        plt.savefig(path.join(self.dir_figs, "energies_splitted.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close()

        # --- log x-axis
        mask = times_arr > 0
        times_log = times_arr[mask]
        y_log = [yi[mask] for yi in y_series]
        fig, ax = plt.subplots()
        _plot_energies_on_ax(ax, times_log, y_log)
        ax.set_xscale("log")
        if title:
            fig.suptitle(title, fontsize=10.0)
        ax.autoscale(enable=True, axis="both", tight=True)
        fig.savefig(path.join(self.dir_figs, "energies_log.png"), bbox_inches="tight")
        fig.savefig(path.join(self.dir_figs, "energies_log.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close(fig)

        # --- delta from average (post-thermalization)
        ct = self.cycles_thermalization
        epot_th = epot[ct:]
        ekin_th = ekin[ct:]
        etot_th = etot[ct:]

        delta_series = [epot - np.mean(epot_th), ekin - np.mean(ekin_th), etot - np.mean(etot_th)]
        delta_labels = [
            f"$E_{{p}} - \\overline{{E}}_{{p}}$,  $\\overline{{E}}_{{p}}={np.mean(epot_th):.3f}$ eV",
            f"$E_{{k}} - \\overline{{E}}_{{k}}$,  $\\overline{{E}}_{{k}}={np.mean(ekin_th):.3f}$ eV",
            f"$E_{{tot}} - \\overline{{E}}_{{tot}}$,  $\\overline{{E}}_{{tot}}={np.mean(etot_th):.3f}$ eV",
        ]
        delta_colors = ["r", "b", "g"]
        delta_lss = ["-", "-", "-"]

        if self.eam:
            epot_eam_th = epot_eam[ct:]
            epot_pair_th = epot_pair[ct:]
            delta_series = [
                epot - np.mean(epot_th),
                epot_pair - np.mean(epot_pair_th),
                epot_eam - np.mean(epot_eam_th),
                ekin - np.mean(ekin_th),
                etot - np.mean(etot_th),
            ]
            delta_labels = [
                f"$E_{{p}} - \\overline{{E}}_{{p}}$,  $\\overline{{E}}_{{p}}={np.mean(epot_th):.3f}$ eV",
                f"$E_{{pair}} - \\overline{{E}}_{{pair}}$,  $\\overline{{E}}_{{pair}}={np.mean(epot_pair_th):.3f}$ eV",
                f"$E_{{eam}} - \\overline{{E}}_{{eam}}$,  $\\overline{{E}}_{{eam}}={np.mean(epot_eam_th):.3f}$ eV",
                f"$E_{{k}} - \\overline{{E}}_{{k}}$,  $\\overline{{E}}_{{k}}={np.mean(ekin_th):.3f}$ eV",
                f"$E_{{tot}} - \\overline{{E}}_{{tot}}$,  $\\overline{{E}}_{{tot}}={np.mean(etot_th):.3f}$ eV",
            ]
            delta_colors = ["r", "orange", "purple", "b", "g"]
            delta_lss = ["-", "--", "--", "-", "-"]

        mask_pos = times_arr > 0
        times_log = times_arr[mask_pos]
        delta_log = [yi[mask_pos] for yi in delta_series]

        fig, ax = plt.subplots()
        for yi, label, c, ls in zip(delta_log, delta_labels, delta_colors, delta_lss):
            ax.plot(times_log, yi, c=c, ls=ls, lw=0.8, label=label)
        ax.axhline(0, c="k", ls=":", lw=0.5)
        ax.set_xscale("log")
        ax.grid(True, alpha=0.3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("$\\Delta E$/atom [eV]")
        ax.legend(framealpha=0.95, fontsize=11)
        ax.autoscale(enable=True, axis="both", tight=True)
        fig.savefig(path.join(self.dir_figs, "energies_delta.png"), bbox_inches="tight")
        fig.savefig(path.join(self.dir_figs, "energies_delta.pdf"), bbox_inches="tight")
        if self.show_graphs:
            plt.show()
        plt.close(fig)

        if self.ensamble == "NVE":
            plt.plot(ekin, epot, "o")
            plt.grid()
            plt.legend()
            plt.xlabel("ekin/atom [eV]")
            plt.ylabel("epot/atom [eV]")
            plt.title(f"slope fit={series.thermal_results.epot_ekin_slope_fit:g}")
            plt.autoscale(enable=True, axis="both", tight=True)
            plt.savefig(path.join(self.dir_figs, "epot_ekin.png"))
            plt.savefig(path.join(self.dir_figs, "epot_ekin.pdf"))
            if self.show_graphs:
                plt.show()
            plt.close()

        plt.plot(series.thermal_results.rg / Units.angstrom, series.thermal_results.gr, "-")
        plt.xlabel("r [A]")

        plt.axhline(1.0, ls="--")
        plt.grid()

        plt.ylabel("$g(r)$")
        plt.autoscale(enable=True, axis="both", tight=True)

        title += f" T={series.temperature_kin.mean:g}K P={series.pressure.mean/Units.GPa:g}GPa"

        if not self.eam:
            title += f", $|P/P_{{g(r)}}-1|$={series.thermal_results.pressure_g_err:g}"
        plt.suptitle(title, fontsize=10)
        plt.savefig(path.join(self.dir_figs, "gr.png"))
        plt.savefig(path.join(self.dir_figs, "gr.pdf"))
        if self.show_graphs:
            plt.show()
        plt.close()

    def partial_pair_correlation(self, *, series):
        if not hasattr(series.thermal_results, "gr_partial"):
            return
        if len(self.atom_names) < 2:
            return

        rg = series.thermal_results.rg / Units.angstrom
        gr_partial = series.thermal_results.gr_partial
        idx = 0
        for b in range(len(self.atom_names)):
            for a in range(b + 1):
                label = f"$g_{{\\mathrm{{{self.atom_names[a]},{self.atom_names[b]}}}}}(r)$"
                plt.plot(rg, gr_partial[idx], label=label)
                idx += 1

        plt.xlabel("r [A]")
        plt.axhline(1.0, ls="--")
        plt.grid()
        plt.ylabel("$g_{a,b}(r)$")
        plt.legend(loc="best", frameon=True, framealpha=0.85)
        plt.autoscale(enable=True, axis="both", tight=True)
        plt.savefig(path.join(self.dir_figs, "gr_partial.png"))
        plt.savefig(path.join(self.dir_figs, "gr_partial.pdf"))
        if self.show_graphs:
            plt.show()
        plt.close()
