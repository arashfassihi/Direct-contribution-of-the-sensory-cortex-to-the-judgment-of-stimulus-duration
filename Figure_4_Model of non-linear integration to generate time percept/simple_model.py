#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:59:37 2019

@author: lucianopaz
"""

import sys
import os
import numpy as np
import matplotlib as mt

mt.use("pdf")
from matplotlib import pyplot as plt
import pandas as pd
from data_interface import load_fittable_data, DT
from theano import tensor as tt
import pymc3 as pm
from pymc3.backends.report import WarningType
import tqdm
from compress_pickle import dump, load
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_ppc
import argparse
import seaborn as sns
from sklearn.mixture import BayesianGaussianMixture


dt = DT


def plot_psth(df):
    gb = df.groupby(["T2", "light", "t"])
    d = gb["ntrials", "psth"].sum()
    y = d.psth / d.ntrials
    y.name = "spike_prob"
    if "ppc" not in df:
        plot_df = y.to_frame().reset_index()
        plot_ppc = False
    else:

        def mean_ppc(df):
            return np.array([v / n for n, v in zip(df.ntrials, df.ppc)]).mean()

        def std_ppc(df):
            return np.array([v / n for n, v in zip(df.ntrials, df.ppc)]).std()

        pred = gb["ntrials", "ppc"].apply(mean_ppc)
        err = gb["ntrials", "ppc"].apply(std_ppc)
        pred.name = "prediction"
        err.name = "pred_err"
        plot_df = pd.concat([y, pred, err], axis=1).reset_index()
        plot_ppc = True
    uT2 = sorted(plot_df.T2.unique())
    f, axs = plt.subplots(
        len(plot_df.T2.unique()), 1, sharex=True, sharey=True, figsize=(8, 9)
    )
    colors = ["b", "g"]
    for (T2, light), df in plot_df.groupby(["T2", "light"]):
        ax = axs[uT2.index(T2)]
        color = colors[light]
        ls = "--" if plot_ppc else "-"
        ax.plot(df.t, df.spike_prob, color=color, linestyle=ls)
        if plot_ppc:
            ax.plot(df.t, df.prediction, color=color)
            ax.fill_between(
                df.t,
                df.prediction + df.pred_err,
                df.prediction - df.pred_err,
                color=color,
                alpha=0.3,
            )


def get_data_arrays(df):
    t = df["t"].values
    offset_t = df["offset_t"].values
    psth = df["psth"].values
    n = df["ntrials"].values
    nid = df["nid"].values
    light = df["light"].values

    unid, nid_inds = np.unique(nid, return_inverse=True)
    return t, offset_t, psth, n, nid, light, unid, nid_inds


def get_training_model(df):
    t, offset_t, psth, n, nid, light, unid, nid_inds = get_data_arrays(df)
    with pm.Model() as model:
        nu = pm.Gamma("nu", mu=30, sigma=5, shape=len(unid))
        lam = pm.Gamma("lam", mu=2.9, sigma=0.4, shape=len(unid))
        #        nu = pm.Gamma("nu", mu=40, sigma=10, shape=len(unid))
        #        lam = pm.Gamma("lam", mu=20, sigma=5, shape=len(unid))
        tau = pm.Gamma("tau", mu=50, sigma=10, shape=len(unid))
        mech_onset = pm.Normal("mech_onset", 0, sigma=10, shape=len(unid))
        opto_onset = pm.Normal("opto_onset", 0, sigma=10, shape=len(unid))
        tau_opto = pm.Gamma("tau_opto", mu=50, sigma=10, shape=len(unid))
        tau_inh = pm.Gamma("tau_inh", mu=30, sigma=10, shape=len(unid))
        I0 = pm.Normal("I0", mu=0, sigma=1.0, shape=len(unid))
        I1 = pm.Normal("I1", mu=0.0, sigma=1.0, shape=len(unid))
        I1p = pm.Normal("I1p", mu=0.0, sigma=1.0, shape=len(unid))
        I2 = pm.Normal("I2", mu=0.0, sigma=1.0, shape=len(unid))
        I2p = pm.Normal("I2p", mu=0.0, sigma=1.0, shape=len(unid))
        I3p = pm.Normal("I3p", mu=0.0, sigma=1.0, shape=len(unid))
        mech_w = (
            tt.nnet.sigmoid(1 * (t - mech_onset[nid_inds]))
            * tt.nnet.sigmoid(-1 * (offset_t - mech_onset[nid_inds]))
        ) * (1 - light)
        opto_w = (
            tt.nnet.sigmoid(1 * (t - opto_onset[nid_inds]))
            * tt.nnet.sigmoid(-1 * (offset_t - opto_onset[nid_inds]))
        ) * light
        opto_w2 = tt.nnet.sigmoid(1 * (offset_t - opto_onset[nid_inds])) * light
        input_current = (
            I0[nid_inds]
            + mech_w
            * (
                I1[nid_inds]
                + (I1p[nid_inds] - I1[nid_inds])
                * tt.exp(-(t - mech_onset[nid_inds]) / tau[nid_inds])
            )
            + opto_w
            * (
                I1[nid_inds]
                + (I1p[nid_inds] - I1[nid_inds])
                * tt.exp(-(t - mech_onset[nid_inds]) / tau[nid_inds])
                + I2[nid_inds]
                + (I2p[nid_inds] - I2[nid_inds])
                * tt.exp(-(t - opto_onset[nid_inds]) / tau_opto[nid_inds])
            )
            + opto_w2
            * I3p[nid_inds]
            * tt.exp(-(offset_t - opto_onset[nid_inds]) / tau_inh[nid_inds])
        )
        io_func = (
            1 - tt.nnet.sigmoid(-10 * input_current) ** (1 / nu[nid_inds])
        ) * lam[nid_inds]
        mu_t = (io_func + 1e-4) * n
        pm.Poisson("psth", mu=mu_t, observed=psth)
    return model


def tuning_curve_model(I, nu_mu=40.0, nu_sigma=10.0, lam_mu=20.0, lam_sigma=5.0):
    with pm.Model() as model:
        nu = pm.Gamma("nu", mu=nu_mu, sigma=nu_sigma)
        lam = pm.Gamma("lam", mu=lam_mu, sigma=lam_sigma)
        pm.Deterministic(
            "F", (1e-4 + (1 - tt.nnet.sigmoid(-10 * I) ** (1 / nu)) * lam) * 1e3 / DT
        )

    return model


def train_model(model):
    with model:
        trace = pm.sample(draws=2000, tune=1000, target_accept=0.85)
    return trace


def test_model(model, trace):
    with model:
        ppc = pm.sample_posterior_predictive(trace)
    return ppc


def store_predictions(worker_id=0, n_workers=1):
    stacked_df = load_fittable_data(light_on_shift=0, unmasked_offset=200)
    unid = np.array(sorted(stacked_df.nid.unique()))
    for nid in tqdm.tqdm(unid[worker_id::n_workers]):
        df = stacked_df.query("nid == {}".format(nid))
        model = get_training_model(df)
        trace = load(os.path.join("simple_model_fits", "{}.gz".format(nid)))
        with model:
            ppc = pm.sample_posterior_predictive(trace, progressbar=False)
            dump(ppc, os.path.join("simple_model_fits", "ppc_{}.gz").format(nid))


def welford_variance_update(obs, prev, axis=None):
    x, w_obs = obs
    mean_prev, M2_prev,  w_prev, w2_prev = prev
    w = np.sum(w_obs + w_prev, axis=axis, keepdims=True)
    w2 = np.sum(w_obs**2 + w2_prev, axis=axis, keepdims=True)
    mean = (
        mean_prev +
        np.sum(w_obs * (x - mean_prev), axis=axis, keepdims=True) / w
    )
    M2 = (
        M2_prev +
        np.sum(w_obs * (x - mean_prev) * (x - mean), axis=axis, keepdims=True)
    )
    return mean, M2, w, w2


def welford_variance_finalize(mean, M2, w, w2, axis=None):
    return (
        np.sum(mean, axis=axis),
        np.sum(M2 / w, axis=axis),
        np.sum(M2 / (w - 1), axis=axis),
        np.sum(M2 / (w - w2 / w), axis=axis)
    )


def predicted_mean_var(*args, **kwargs):
    stacked_df = load_fittable_data(light_on_shift=0, unmasked_offset=200)
    groupby = ["light", "T2", "t"]
    gb = stacked_df.groupby(groupby)
    columns = ["ntrials", "_mean", "_var"]
    index = pd.MultiIndex.from_tuples(gb.indices.keys(), names=groupby)

    unid = np.array(sorted(stacked_df.nid.unique()))
    output = pd.DataFrame(
        data=np.zeros((len(index), len(columns))),
        index=index,
        columns=columns,
    )
    for nid in tqdm.tqdm(unid):
        df = stacked_df.query("nid == {}".format(nid))
        ppc = load(os.path.join("simple_model_fits", "ppc_{}.gz".format(nid)))
        df = df.assign(ppc=[v for v in ppc["psth"].T])
        for ind, subdf in df.groupby(groupby):
            w = subdf.ntrials.values
            x = np.array([v for v in subdf.ppc.values]).T * 100
            output.loc[ind, "ntrials"] += w
            output.loc[ind, "_mean"] += np.mean(x, axis=0)
    output.loc[:, "_mean"] = output.loc[:, "_mean"] / output.loc[:, "ntrials"]
    for nid in tqdm.tqdm(unid):
        df = stacked_df.query("nid == {}".format(nid))
        ppc = load(os.path.join("simple_model_fits", "ppc_{}.gz".format(nid)))
        df = df.assign(ppc=[v for v in ppc["psth"].T])
        for ind, subdf in df.groupby(groupby):
            w = subdf.ntrials.values
            x = np.array([v for v in subdf.ppc.values]).T * 100
            output.loc[ind, "_var"] += np.sum(
                (x - w * output.loc[ind, "_mean"])**2,
                axis=0,
            )
    output.loc[:, "_var"] = output.loc[:, "_var"] / output.loc[:, "ntrials"] / len(x)
    output["_std"] = np.sqrt(output["_var"])
    dump(output, os.path.join("simple_model_fits", "population_ppc.gz".format(nid)))


def plot_pred(df, ppc):
    temp_df = df.assign(
        ppc=np.mean(ppc["psth"], axis=0), eppc=np.std(ppc["psth"], axis=0)
    )
    plot_psth(temp_df)


def fit(worker_id=0, n_workers=1):
    stacked_df = load_fittable_data(light_on_shift=0, unmasked_offset=200)
    unid = np.array(sorted(stacked_df.nid.unique()))
    for nid in tqdm.tqdm(unid[worker_id::n_workers]):
        df = stacked_df.query("nid == {}".format(nid))
        model = get_training_model(df)
        try:
            trace = train_model(model)
        except Exception as e:
            dump(e, "simple_model_fits/{}.gz".format(nid))
            continue
        dump(trace, "simple_model_fits/{}.gz".format(nid))


def plot_trace(worker_id=0, n_workers=1):
    fnames = sorted(
        [f for f in os.listdir("simple_model_fits") if f.endswith(".gz")],
        key=lambda x: int(x.split(".")[0]),
    )
    saver = PdfPages(os.path.join("simple_model_fits", "traces.pdf"))
    with saver:
        for fname in tqdm.tqdm(fnames[worker_id::n_workers]):
            trace = load(os.path.join("simple_model_fits", fname))
            axs = pm.traceplot(trace)
            fig = axs[0, 0].figure
            fig.suptitle("Neuron {}".format(fname.split(".")[0]), fontsize=16)
            saver.savefig(fig)
            plt.close(fig.number)


def mean_parameters(worker_id=0, n_workers=1):
    fnames = sorted(
        [f for f in os.listdir("simple_model_fits") if f.endswith(".gz")],
        key=lambda x: int(x.split(".")[0]),
    )
    data = []
    for fname in tqdm.tqdm(fnames[worker_id::n_workers]):
        nid = int(fname.split(".")[0])
        trace = load(os.path.join("simple_model_fits", fname))
        if not trace.report.ok:
            if not all(
                (w[0] == WarningType.BAD_ACCEPTANCE for w in trace.report._warnings)
            ):
                continue
        df = pm.trace_to_dataframe(trace)
        df = (
            df.agg([np.mean, np.std])
            .rename(columns=lambda x: x.split("__0")[0])
            .unstack()
            .to_frame()
            .T.assign(nid=nid)
            .set_index("nid")
        )
        data.append(df)
    df = pd.concat(data)
    dump(df, os.path.join("simple_model_fits", "parameters.pkl"))
    return df


def analyze_parameters(worker_id=0, n_workers=1):
    df = load(os.path.join("simple_model_fits", "parameters.pkl"))
    # Raw pairplot
    sns.pairplot(data=df.xs("mean", axis=1, level=1))
    plt.savefig(os.path.join("simple_model_fits", "all_parameters.pdf"))
    plt.close()
    # Cluster and then pairplot
    gmm = BayesianGaussianMixture(
        n_components=4,
        covariance_type="full",
        n_init=5,
        weight_concentration_prior_type="dirichlet_process",
        weight_concentration_prior=1e-6,
    )
    train_df = df.xs("mean", axis=1, level=1).drop(
        columns=["mech_onset", "opto_onset", "I0", "lam", "tau", "tau_inh", "tau_opto"]
    )
    label = pd.Series(data=gmm.fit_predict(train_df), index=train_df.index)
    plt.bar(range(len(gmm.weights_)), gmm.weights_)
    plt.ylabel("Class weights")
    plt.savefig(os.path.join("simple_model_fits", "parameter_clustering_weights.pdf"))
    plt.close()
    # All parameters according to clusters
    viz_df1 = df.xs("mean", axis=1, level=1).assign(label=label)
    sns.pairplot(viz_df1, hue="label")
    plt.savefig(os.path.join("simple_model_fits", "parameter_full_clustering.pdf"))
    plt.close()
    # Main parameters according to clusters
    viz_df2 = train_df.assign(label=label)
    sns.pairplot(
        viz_df2,
        vars=[c for c in viz_df2.columns if c != "label"],
        diag_kind="kde",
        hue="label",
    )
    plt.savefig(os.path.join("simple_model_fits", "parameter_reduced_clustering.pdf"))
    plt.close()


def analyze_significance(worker_id=0, n_workers=1):
    fnames = sorted(
        [f for f in os.listdir("simple_model_fits") if f.endswith(".gz")],
        key=lambda x: int(x.split(".")[0]),
    )
    pvals = []
    for fname in tqdm.tqdm(fnames[worker_id::n_workers]):
        trace = load(os.path.join("simple_model_fits", fname))
        pval = [
            int(fname.split(".")[0]),
            np.mean(trace["I1"] > 0),
            np.mean(trace["I1p"] > 0),
            np.mean(trace["I2"] > 0),
            np.mean(trace["I2p"] > 0),
            np.mean(trace["I3p"] > 0),
        ]
        pvals.append(pval)
    df = pd.DataFrame(
        data=np.array(pvals), columns=["nid", "I1", "I1p", "I2", "I2p", "I3p"]
    )
    df = df.assign(
        mech=(
            (df.I1 >= 0.995) | (df.I1 <= 0.005) | (df.I1p >= 0.995) | (df.I1p <= 0.005)
        )
    )
    df = df.assign(
        opto=(
            (df.I2 >= 0.995)
            | (df.I2 <= 0.005)
            | (df.I2p >= 0.995)
            | (df.I2p <= 0.005)
            | (df.I3p >= 0.995)
            | (df.I3p <= 0.005)
        )
    )
    df = df.assign(resp=df.mech | df.opto)
    dump(df, os.path.join("simple_model_fits", "significant.pkl"))


def plot_significance(*args, **kwargs):
    df = load(os.path.join("simple_model_fits", "significant.pkl"))
    pie_df = pd.concat(
        [
            (df.mech & df.opto),
            (np.logical_not(df.mech) & df.opto),
            (df.mech & np.logical_not(df.opto)),
            (np.logical_not(df.mech) & np.logical_not(df.opto)),
        ],
        axis=1,
    )
    pie_df = pie_df.rename(
        columns={
            k: v
            for k, v in enumerate(
                [
                    "Mech and Opto",
                    "Not mech but Opto",
                    "Mech and not Opto",
                    "neither Mech nor Opto",
                ]
            )
        }
    )
    f = plt.figure(figsize=(10, 10))
    _, _, autotexts = plt.pie(
        pie_df.mean().values,
        autopct="%1.1f%%",
        textprops={"fontsize": 19},
        explode=[0, 0, 0, 0.1],
        labels=list(pie_df.columns.values),
        labeldistance=1.1,
        shadow=True,
    )
    for a in autotexts:
        a.set_color("w")
    plt.savefig(os.path.join("simple_model_fits", "significant_response.pdf"))
    return f


def plot_predictions(worker_id=0, n_workers=1):
    stacked_df = load_fittable_data(light_on_shift=0, unmasked_offset=200)
    unid = np.array(sorted(stacked_df.nid.unique()))
    saver = PdfPages(os.path.join("simple_model_fits", "predictions.pdf"))
    with saver:
        for nid in tqdm.tqdm(unid[worker_id::n_workers]):
            df = stacked_df.query("nid == {}".format(nid))
            model = get_training_model(df)
            try:
                ppc = load(os.path.join("simple_model_fits", "ppc_{}.gz".format(nid)))
            except IOError:
                trace = load(os.path.join("simple_model_fits", "{}.gz".format(nid)))
                with model:
                    ppc = pm.sample_posterior_predictive(trace, progressbar=False)
                    dump(ppc, os.path.join("simple_model_fits", "ppc_{}.gz".format(nid)))
                    t, offset_t, psth, n = get_data_arrays(df)[:4]
                    fig = plot_ppc(t, n, psth, ppc["psth"], df, plot_ppc_stats=True)[0]
                    fig.suptitle("Neuron {}".format(nid), fontsize=16)
                    saver.savefig(fig)
                    plt.close(fig.number)


def _test():
    stacked_df = load_fittable_data(
        light_on_shift=0,
        unmasked_offset=200,
        version="_v2"
    )
    nid = 4
    ndf = stacked_df.query("nid == {}".format(nid))
    model = get_training_model(ndf)
    trace = load(os.path.join("simple_model_fits", "{}.gz".format(nid)))
    with model:
        ppc = pm.sample_posterior_predictive(trace, progressbar=True)
        t, offset_t, psth, n = get_data_arrays(ndf)[:4]
        fig, axs = plot_ppc(
            t,
            n,
            psth,
            ppc["psth"],
            ndf,
            plot_ppc_stats=True,
            axesby=[],
            groupby=["light"],
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit and analyze simple models for vS1 activity"
    )
    parser.add_argument(
        "--fit", help="Fit the model", dest="tasks", action="append_const", const="fit"
    )
    parser.add_argument(
        "--plot_trace",
        help="Plot the traces of the model",
        dest="tasks",
        action="append_const",
        const="plot_trace",
    )
    parser.add_argument(
        "--plot_predictions",
        help="Plot the predictions of the model",
        dest="tasks",
        action="append_const",
        const="plot_predictions",
    )
    parser.add_argument(
        "--mean_parameters",
        help="Compute mean of the fitted parameters for neurons",
        dest="tasks",
        action="append_const",
        const="mean_parameters",
    )
    parser.add_argument(
        "--analyze_parameters",
        help="Analyze mean parameters of neurons",
        dest="tasks",
        action="append_const",
        const="analyze_parameters",
    )
    parser.add_argument(
        "--analyze_significance",
        help="Analyze significant response of neurons",
        dest="tasks",
        action="append_const",
        const="analyze_significance",
    )
    parser.add_argument(
        "--plot_significance",
        help="Plot pie chart of significantly responsive of neurons",
        dest="tasks",
        action="append_const",
        const="plot_significance",
    )
    parser.add_argument(
        "--store_predictions",
        help="Store the posterior predictive traces for each neuron",
        dest="tasks",
        action="append_const",
        const="store_predictions",
    )
    parser.add_argument(
        "--predicted_mean_var",
        help="Computed the predicted population's mean and variance",
        dest="tasks",
        action="append_const",
        const="predicted_mean_var",
    )
    parser.add_argument(
        "--worker_id",
        nargs=1,
        default=[0],
        type=int,
        help="Worker id of a parallel workers pool",
        dest="worker_id",
    )
    parser.add_argument(
        "--n_workers",
        nargs=1,
        default=[1],
        type=int,
        help="Total size of parallel pool of workers",
        dest="n_workers",
    )
    args = vars(parser.parse_args(sys.argv[1:]))
    worker_id = args["worker_id"][0]
    n_workers = args["n_workers"][0]
    tasks = args.get("tasks", None)
    if tasks is None:
        tasks = ["fit"]
    print(tasks)
    for task in tasks:
        eval("{}(worker_id, n_workers)".format(task))
