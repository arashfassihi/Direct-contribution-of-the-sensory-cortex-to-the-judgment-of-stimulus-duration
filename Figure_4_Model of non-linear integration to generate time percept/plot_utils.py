#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 14:09:10 2019

@author: lucianopaz
"""

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from utils import to_tuple


__all__ = ["crude_pair_plot", "plot_ppc", "unstack_plot_data"]


def get_flat_points(trace, model=None, varnames=None):
    try:
        raw_vnames = trace.varnames
    except AttributeError:
        raw_vnames = list(trace.keys())
    if varnames is not None:
        if model is not None:
            model_varnames = [model.name_of(v) for v in varnames]
            raw_vnames = [
                v
                for v in raw_vnames
                if v in varnames + model_varnames
                or model.name_of(v) in varnames + model_varnames
            ]
        else:
            raw_vnames = [v for v in raw_vnames if v in varnames]
    try:
        test_point = trace.point(0)
        iterator = trace.points()
    except AttributeError:
        test_point = {
            k: vv for k, vv in zip(trace.keys(), [v[0] for v in trace.values()])
        }

        def get_iterator():
            keys = list(trace.keys())
            L = len(trace[keys[0]])
            for i in range(L):
                yield {k: np.asarray(trace[k][i]) for k in keys}

        iterator = get_iterator()
    vsizes = {vname: test_point[vname].size for vname in raw_vnames}
    vnames = []
    points = []
    for i, raw_point in enumerate(iterator):
        point = []
        for vname in raw_vnames:
            vsize = vsizes[vname]
            vals = raw_point[vname].flatten()
            for var_ind, val in enumerate(vals):
                point.append(val)
                if i == 0:
                    if model is None:
                        name = vname
                    else:
                        name = model.name_of(vname)
                    if vsize == 1:
                        vnames.append(name)
                    else:
                        vnames.append(name + "_{}".format(var_ind))
        points.append(point)
    points = np.array(points)
    return points, vnames


def crude_pair_plot(trace, model=None, varnames=None):
    points, vnames = get_flat_points(trace, model=model, varnames=varnames)
    try:
        divs = trace.get_sampler_stats("diverging")
    except Exception:
        divs = np.zeros(len(points))
    fig, axs = plt.subplots(
        nrows=len(vnames), ncols=len(vnames), figsize=(8 + len(vnames), 8 + len(vnames))
    )
    for i, row_vname in enumerate(vnames):
        for j, col_vname in enumerate(vnames):
            ax = axs[i, j]
            if j > i:
                ax.axis("off")
                continue
            elif j == i:
                sns.kdeplot(points[:, i], ax=ax)
                if i < len(vnames) - 1:
                    ax.tick_params(axis="x", labelbottom=False)
            else:
                inds0 = np.logical_not(divs)
                if any(inds0):
                    ax.scatter(
                        x=points[inds0, j],
                        y=points[inds0, i],
                        c="b",
                        marker="o",
                        alpha=0.2,
                    )
                if any(divs):
                    ax.scatter(
                        x=points[divs, j],
                        y=points[divs, i],
                        c="r",
                        marker="o",
                        alpha=0.5,
                    )
                if i < len(vnames) - 1:
                    ax.tick_params(axis="x", labelbottom=False)
                if j > 0:
                    ax.tick_params(axis="y", labelleft=False)
            if j == 0:
                ax.set_ylabel(row_vname)
            if i == len(vnames) - 1:
                ax.set_xlabel(col_vname)
    fig.tight_layout()
    return fig, axs


def unstack_plot_data(t, n, psth, ppc=None, group_labels=None, return_firing_rate=True):
    if group_labels is not None:
        if group_labels.ndim == 2:
            temp = np.concatenate([group_labels, t[None, :]], axis=0)
        else:
            temp = np.array([group_labels, t])
        labels, melt_labels = np.unique(temp, return_inverse=True, axis=1)
        if ppc is not None:
            n, psth, ppc = sum_psths(melt_labels, n, psth, ppc)
        else:
            n, psth = sum_psths(melt_labels, n, psth)
        labels, t = labels
        ulabels = np.unique(labels)
        n_ = []
        psth_ = []
        ppc_ = []
        t_ = []
        for label in ulabels:
            inds = labels == label
            t_.append(t[inds])
            n_.append(n[..., inds])
            psth_.append(psth[..., inds])
            if ppc is not None:
                ppc_.append(ppc[..., inds])
    else:
        slice_ends = np.nonzero(np.diff(t) < 0)[0]
        if slice_ends.size > 0:
            slice_ends += 1
            slicers = (
                [slice(None, slice_ends[0])]
                + [slice(s0, s1) for s0, s1 in zip(slice_ends[:-1], slice_ends[1:])]
                + [slice(slice_ends[-1], None)]
            )
        else:
            slicers = [slice(None)]
        ulabels = list(range(len(slicers)))
        n_ = []
        psth_ = []
        ppc_ = []
        t_ = []
        for slicer in slicers:
            t_.append(t[slicer])
            n_.append(n[slicer])
            psth_.append(psth[slicer])
            if ppc is not None:
                ppc_.append(ppc[(slice(None),) + slicer])
    if return_firing_rate:
        if ppc is not None:
            for i, (t, n, psth, ppc) in enumerate(zip(t_, n_, psth_, ppc_)):
                dt = np.diff(t)
                dt = np.concatenate([np.array([dt[0]]), dt])
                psth_[i] = psth / n / dt * 1e3
                ppc_[i] = ppc / n / dt * 1e3
        else:
            for i, (t, n, psth) in enumerate(zip(t_, n_, psth_)):
                dt = np.diff(t)
                dt = np.concatenate([np.array([dt[0]]), dt])
                psth_[i] = psth / n / dt * 1e3
    if ppc is None:
        ppc_ = None
    return t_, n_, psth_, ppc_, ulabels


def sum_psths(group_labels, *args):
    ulabels = np.unique(group_labels)
    outputs = [[] for _ in args]
    for label in ulabels:
        inds = group_labels == label
        for i, arg in enumerate(args):
            outputs[i].append(np.sum(arg[..., inds], axis=-1))
    for i, out in enumerate(outputs):
        outputs[i] = np.moveaxis(np.array(out), 0, -1)
    return tuple(outputs)


def plot_ppc(
    t,
    n,
    psth,
    ppc,
    stacked_df,
    groupby=["T2", "light"],
    labelby=["light"],
    axesby=["T2"],
    cmap="winter_r",
    ppc_step=100,
    plot_firing_rate=True,
    plot_ppc_stats=False,
    figure=None,
    gridspec_kw=None,
    plot_legend=True,
    *args,
    **kwargs
):
    groupby = to_tuple(groupby)
    labelby = to_tuple(labelby)
    axesby = to_tuple(axesby)

    ugroupby, group_labels = np.unique(
        stacked_df.loc[:, groupby], axis=0, return_inverse=True
    )
    t, n, psth, ppc, unstacked_labels = unstack_plot_data(
        t=t,
        n=n,
        psth=psth,
        ppc=ppc,
        group_labels=group_labels,
        return_firing_rate=plot_firing_rate,
    )
    group_labels = np.array([ugroupby[l] for l in unstacked_labels.astype(np.int)])
    if plot_firing_rate:
        ylabel = "Firing rate [Hz]"
    else:
        ylabel = "PSTH"

    labelby_axis = []
    for label in labelby:
        if label not in groupby:
            raise ValueError(
                "Label {} was not used to unstack the data and "
                "cannot be used in labelby.".format(label)
            )
        labelby_axis.append(groupby.index(label))
    labelby_values = ugroupby[:, labelby_axis]
    ulabels, temp_inds = np.unique(labelby_values, axis=0, return_inverse=True)
    ulabels = [
        ", ".join(
            [
                "{} = {:1.2f}".format(label_name, label_value)
                for label_name, label_value in zip(labelby, ulabel)
            ]
        )
        for ulabel in ulabels
    ]
    group_to_labels_map = {}
    for old, new in enumerate(temp_inds):
        group_to_labels_map[old] = ulabels[new]
    colors = {
        label: plt.get_cmap(cmap)(x)
        for label, x in zip(ulabels, np.linspace(0, 1, len(ulabels)))
    }

    axesby_axis = []
    for label in axesby:
        if label not in groupby:
            raise ValueError(
                "Label {} was not used to unstack the data and "
                "cannot be used in axesby.".format(label)
            )
        axesby_axis.append(groupby.index(label))
    if len(axesby_axis) > 0:
        axesby_values = ugroupby[:, axesby_axis]
        uaxes, temp_inds = np.unique(axesby_values, axis=0, return_inverse=True)
    else:
        uaxes = [0]
        temp_inds = np.zeros(ugroupby.shape[0], dtype=np.int)
    group_to_axes_map = {}
    for old, new in enumerate(temp_inds):
        group_to_axes_map[old] = new

    if figure is None:
        figure, axs = plt.subplots(
            len(uaxes), 1, figsize=(8, 12), sharex=True, gridspec_kw=gridspec_kw
        )
        gs = None
        if len(uaxes) == 1:
            axs = [axs]
    else:
        if gridspec_kw is None:
            gridspec_kw = {}
        gs = plt.GridSpec(len(uaxes), 1, figure=figure, **gridspec_kw)
        ax0 = figure.add_subplot(gs[0])
        axs = [ax0]
        axs.extend(
            [figure.add_subplot(gs[i], sharex=ax0) for i in range(1, len(uaxes))]
        )
        axs = np.array(axs)
    if ppc is not None:
        iterator = enumerate(zip(t, psth, ppc))
    else:
        iterator = enumerate(zip(t, psth))
    for group_key, it in iterator:
        if ppc is not None:
            t_, psth_, ppc_ = it
        else:
            t_, psth_ = it
            ppc_ = None
        ax = axs[group_to_axes_map[group_key]]
        label = group_to_labels_map[group_key]
        color = colors[label]
        _plot_ppc(
            t=t_,
            psth=psth_,
            ppc=ppc_,
            label=label,
            ax=ax,
            ppc_step=ppc_step,
            color=color,
            plot_ppc_stats=plot_ppc_stats,
            *args,
            **kwargs
        )
        ax.set_ylabel(ylabel)
    if plot_legend:
        axs[0].legend(loc="best")
    axs[-1].set_xlabel("time [ms]")
    return figure, axs


def _plot_ppc(
    t,
    psth,
    ppc,
    label,
    ax,
    ppc_step=100,
    color=None,
    plot_ppc_stats=False,
    *args,
    **kwargs
):
    l, = ax.step(t, psth, linestyle="-", label=label, color=color, *args, **kwargs)
    if ppc is not None:
        if plot_ppc_stats:
            mppc = np.mean(ppc, axis=0)
            sppc = np.std(ppc, axis=0)
            ax.plot(t, mppc, "--", color=l.get_color(), *args, **kwargs)
            ax.fill_between(t, mppc + sppc, mppc - sppc, alpha=0.3, color=l.get_color())
        else:
            ax.step(
                t,
                ppc[::ppc_step].T,
                "--",
                alpha=0.1,
                color=l.get_color(),
                *args,
                **kwargs
            )
