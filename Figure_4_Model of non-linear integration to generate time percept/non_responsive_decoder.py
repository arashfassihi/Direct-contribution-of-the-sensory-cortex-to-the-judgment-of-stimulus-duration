#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 21:29:25 2019

@author: lucianopaz
"""

import numpy as np
import pymc3 as pm
from theano import tensor as tt
from compress_pickle import dump, load
from data_interface import load_fittable_data, DT
from utils import alpha_function
import sys
import os
import re
import pandas as pd
import matplotlib as mt
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_ppc
from tqdm import tqdm
import warnings

mt.use("pdf")
mt.style.use("seaborn")


def get_data(nid=None, light_on_shift=DT, **kwargs):
    if nid is not None:
        query = "nid == {}".format(nid)
    else:
        query = None
    return load_fittable_data(light_on_shift=light_on_shift, query=query, **kwargs)


def model_factory(df, kind="non_responsive", dt=DT, light_cond="all",
                  onset_jitter=None):
    if light_cond == "all":
        pass
    else:
        df = df.query("light == {}".format(light_cond))
    t = df["t"].values
    offset_t = df["offset_t"].values
    psth = df["psth"].values
    n = df["ntrials"].values
    nid = df["nid"].values
    light = df["light"].values

    unid, nid_inds = np.unique(nid, return_inverse=True)
    ucat, inds = np.unique(np.array([light, nid]), return_inverse=True, axis=1)

    with pm.Model(kind) as model:
        if onset_jitter is None:
            onset_jitter = pm.Normal('onset_jitter', 0, 1, shape=ucat.shape[1])
        mech_stim_on = tt.and_(
            t >= onset_jitter[inds],
            offset_t < onset_jitter[inds]
        )
        mech_stim_off = tt.or_(
            t < onset_jitter[inds],
            offset_t >= onset_jitter[inds]
        )
        background = pm.HalfNormal("background", sigma=1.0, shape=len(unid))
        if kind != "non_responsive":
            steady = pm.HalfNormal("steady", sigma=1.0, shape=ucat.shape[1])
            mu_t = tt.zeros(len(t))
            mu_t = tt.inc_subtensor(
                mu_t[mech_stim_off], background[nid_inds][mech_stim_off]
            )
            mu_t = tt.inc_subtensor(mu_t[mech_stim_on], steady[inds][mech_stim_on])
            mu_t = mu_t * n# * 1e-3 * dt
        else:
            mu_t = background[nid_inds] * n# * 1e-3 * dt
        pm.Poisson("psth", mu=mu_t, observed=psth)
    return model


def compare_models(df, light_cond="all", return_details=False, *args, **kwargs):
    with model_factory(df, kind="non_responsive", light_cond=light_cond) as nr_model:
        nr_trace = pm.sample(2000, tune=1000, *args, **kwargs)
    with model_factory(df, kind="responsive", light_cond=light_cond) as r_model:
        r_trace = pm.sample(2000, tune=1000, *args, **kwargs)
    model_dict = {nr_model: nr_trace, r_model: r_trace}
    if return_details:
        return pm.compare(model_dict), nr_model, r_model, nr_trace, r_trace
    else:
        return pm.compare(model_dict)


def main(nid, fit_path="non_responsive", figsaver=None):
    df = get_data(nid=nid)
    comp_all = compare_models(df, light_cond="all", return_details=True)
    if figsaver is not None:
        fig = plot_prediction(df, nid, *comp_all[1:],
                              significant=analyze(comp_all[0]))
        figsaver.savefig(fig, bbox_inches="tight")
#    comp_light = compare_models(df, light_cond=1)
#    comp_mech = compare_models(df, light_cond=0)
    f = os.path.join(fit_path, "nid_{}.gz".format(nid))
    dump(comp_all[0], f)
    return comp_all
#    dump((comp_all[0], comp_light, comp_mech), f)
#    return comp_all, comp_light, comp_mech


def fit_jitter(df, light_cond="all"):
    jitter_range = np.arange(-3, 4) * 10
    logps = np.zeros((len(jitter_range), len(jitter_range)))
    for i, loff_jitter in enumerate(jitter_range):
        for j, lon_jitter in enumerate(jitter_range):
            onset_jitter = np.array([loff_jitter, lon_jitter])
            reg = -np.sum((onset_jitter/5)**2)
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                with model_factory(df, kind="responsive", light_cond=light_cond,
                                   onset_jitter=onset_jitter) as model:
                    params = pm.find_MAP(progressbar=False)
                logps[i, j] = model.logp(params) + reg
    idx = np.argmax(logps)
    idx0, idx1 = np.divmod(idx, len(jitter_range))
    jitter = [jitter_range[idx0], jitter_range[idx1]]
    return jitter, logps


def main_jitter(nid, fit_path="non_responsive", figsaver=None):
    f = os.path.join(fit_path, "nid_{}.gz".format(nid))
    comp_all = load(f)
    if isinstance(comp_all, tuple):
        comp_all = comp_all[0]
    if analyze(comp_all):
        df = get_data(nid=nid)
        jitter = fit_jitter(df)
    else:
        jitter = np.zeros(2), None
    f = os.path.join(fit_path, "jitter_nid_{}.gz".format(nid))
    dump(jitter, f)
    return jitter


def plot_prediction(
    df,
    nid,
    non_responsive_model,
    responsive_model,
    non_responsive_trace,
    responsive_trace,
    significant,
):
    with non_responsive_model:
        nr_ppc = pm.sample_posterior_predictive(non_responsive_trace)
    with responsive_model:
        r_ppc = pm.sample_posterior_predictive(responsive_trace)
    t = df.t.values
    n = df.ntrials.values
    psth = df.psth.values
    fig, axs = plot_ppc(
        t,
        n,
        psth,
        nr_ppc[non_responsive_model.name_for("psth")],
        df,
        gridspec_kw={"left": 0, "right": 0.45, "bottom": 0, "top": 0.95},
        plot_ppc_stats=True,
    )
    axs[0].set_title("Non responsive model")
    fig, axs = plot_ppc(
        t,
        n,
        psth,
        r_ppc[responsive_model.name_for("psth")],
        df,
        gridspec_kw={"left": 0.55, "right": 1, "bottom": 0, "top": 0.95},
        figure=fig,
        plot_ppc_stats=True,
    )
    axs[0].set_title("Responsive model")
    fig.suptitle("Neuron {}".format(nid))
    return fig


def significant_responsive(comp):
    return (
        (comp.loc["responsive", "dWAIC"] + 2 * comp["dSE"].max()) <
        comp.loc["non_responsive", "dWAIC"]
    )


def analyze(comp):
    if isinstance(comp, (tuple, list)):
        return [significant_responsive(c) for c in comp]
    else:
        return significant_responsive(comp)


def analyze_all():
    fit_path = "non_responsive"
    nids = []
    resp = []
    columns = None
    for fname in os.listdir(fit_path):
        try:
            nid = int(re.search("[0-9]+", fname).group())
        except AttributeError:
            continue
        comp = load(os.path.join(fit_path, fname))
        an = analyze(comp)
        if columns is None or isinstance(an, (tuple, list)):
            columns = ["all", "light", "mech"]
        else:
            columns = ["all"]
        resp.append(an)
        nids.append(nid)
    return pd.DataFrame(np.array(resp), columns=columns, index=nids).sort_index()


if __name__ == "__main__":
    if len(sys.argv) < 3:
        offset = 0
        njobs = 1
    else:
        offset, njobs = sys.argv[1:3]
    jitters = []
    for nid in tqdm(range(int(offset), 228, int(njobs))):
        jitters.append(main_jitter(nid))
    dump(jitters, "neuron_jitters.pkl")
#    with PdfPages(os.path.join("non_responsive", "fits2.pdf")) as figsaver:
#        for nid in range(int(offset), 228, int(njobs)):
#            main(nid, figsaver=figsaver)
#    out = analyze_all()
#    dump(out, "responsive_neurons2.pkl")
