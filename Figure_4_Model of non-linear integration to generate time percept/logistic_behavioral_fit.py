#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 11:52:39 2019

@author: lucianopaz
"""

import os
import numpy as np
import pymc3 as pm
import arviz as az
from theano import tensor as tt
from scipy.special import expit
from scipy.special import logit as nlogit
from pymc3.math import sigmoid, logit
from data_interface import load_behavioral_data
from matplotlib import pyplot as plt
from compress_pickle import load, dump


def numpy_prob(x, beta, pse, lapse_rate, lapse_bias):
    return lapse_rate * lapse_bias + (1 - lapse_rate) * expit(beta * (x - pse))


def prob(x, beta, pse, lapse_rate, lapse_bias):
    return lapse_rate * lapse_bias + (1 - lapse_rate) * sigmoid(beta * (x - pse))


def numpy_inv_prob(y, beta, pse, lapse_rate, lapse_bias):
    return nlogit((y - lapse_rate * lapse_bias) / (1 - lapse_rate)) / beta + pse


def inv_prob(y, beta, pse, lapse_rate, lapse_bias):
    return logit((y - lapse_rate * lapse_bias) / (1 - lapse_rate)) / beta + pse


def model_factory(ntd, ntrials=1, choices=None, track_p=False):
    if choices.ndim < 2:
        n = None
    else:
        n = choices.shape[1]
    with pm.Model() as model:
        beta = pm.Normal("beta", 10, 1, shape=n)
        pse = pm.Normal("pse", 0, 10, shape=n)
        lapse_rate = pm.Uniform("lapse_rate", 0, 1, shape=n)
        lapse_bias = pm.Uniform("lapse_bias", 0, 1, shape=n)
        p = prob(ntd, beta, pse, lapse_rate, lapse_bias)
        if track_p:
            p = pm.Deterministic("p", p)
        pm.Binomial("choices", ntrials, p, observed=choices, shape=n)
    return model


def fit(ntd, ntrials, choices, draws=2000, tune=1000, **kwargs):
    with model_factory(ntd, ntrials, choices) as model:
        pchoice = pm.Deterministic(
            "pchoice",
            100 * tt.sum(model.choices.distribution.p * ntrials, axis=0) /
            np.sum(ntrials, axis=0)
        )
        bpse = inv_prob(
            0.5, model.beta, model.pse, model.lapse_rate, model.lapse_bias
        )
        pm.Deterministic(
            "contralateral_inflection_point_shift",
            model.pse[0] - model.pse[1]
        )
        pm.Deterministic(
            "ipsilateral_inflection_point_shift",
            model.pse[2] - model.pse[3]
        )
        pm.Deterministic(
            "inflection_point_shift_difference",
            model.contralateral_inflection_point_shift -
            model.ipsilateral_inflection_point_shift
        )
        pm.Deterministic(
            "contralateral_pse_shift",
            bpse[0] - bpse[1]
        )
        pm.Deterministic(
            "ipsilateral_pse_shift",
            bpse[2] - bpse[3]
        )
        pm.Deterministic(
            "pse_shift_difference",
            model.contralateral_pse_shift -
            model.ipsilateral_pse_shift
        )
        pm.Deterministic(
            "contralateral_behavioral_bias",
            pchoice[1] - pchoice[0]
        )
        pm.Deterministic(
            "ipsilateral_behavioral_bias",
            pchoice[3] - pchoice[2]
        )
        pm.Deterministic(
            "behavioral_bias_difference",
            model.contralateral_behavioral_bias -
            model.ipsilateral_behavioral_bias
        )
        trace = pm.sample(draws=draws, tune=tune, **kwargs)
    return model, trace


def store_behavioral_fit():
    behavioral_data = load_behavioral_data()
    T = behavioral_data.index.values.astype("float")
    ntd = (T - 334) / (T + 334)
    ntrials = np.concatenate(
        (
            behavioral_data.xs("contra", axis=1).values[:, :2],
            behavioral_data.xs("ipsi", axis=1).values[:, :2]
        ),
        axis=1
    )
    choices = np.concatenate(
        (
            behavioral_data.xs("contra", axis=1).values[:, 2:],
            behavioral_data.xs("ipsi", axis=1).values[:, 2:]
        ),
        axis=1
    )
    model, trace = fit(ntd[..., None], ntrials, choices)
    dump(
        {"trace": pm.trace_to_dataframe(trace)},
        os.path.join("simple_model_fits", "logistic_behavioral_fit.pkl")
    )


def load_behavioral_fit():
    try:
        trace = load(
            os.path.join("simple_model_fits", "logistic_behavioral_fit.pkl")
        )
    except Exception:
        store_behavioral_fit()
        trace = load(
            os.path.join("simple_model_fits", "logistic_behavioral_fit.pkl")
        )
    return trace


if __name__ == "__main__":
    behavioral_data = load_behavioral_data()
    T = behavioral_data.index.values.astype("float")
    ntd = (T - 334) / (T + 334)
    ntrials = np.concatenate(
        (
            behavioral_data.xs("contra", axis=1).values[:, :2],
            behavioral_data.xs("ipsi", axis=1).values[:, :2]
        ),
        axis=1
    )
    choices = np.concatenate(
        (
            behavioral_data.xs("contra", axis=1).values[:, 2:],
            behavioral_data.xs("ipsi", axis=1).values[:, 2:]
        ),
        axis=1
    )
    model, trace = fit(ntd[..., None], ntrials, choices)
    pm.plot_posterior(
        trace,
        var_names=[
            "contralateral_inflection_point_shift",
            "ipsilateral_inflection_point_shift",
            "inflection_point_shift_difference",
            "contralateral_pse_shift",
            "ipsilateral_pse_shift",
            "pse_shift_difference",
            "contralateral_behavioral_bias",
            "ipsilateral_behavioral_bias",
            "behavioral_bias_difference",
        ]
    )

    plot_ntd = np.linspace(-0.4, 0.4, 1000)
    plot_p = numpy_prob(
        plot_ntd[..., None, None],
        trace["beta"],
        trace["pse"],
        trace["lapse_rate"],
        trace["lapse_bias"],
    )

    blue = np.array([62, 105, 178]) / 255
    green = np.array([106, 189, 69]) / 255
    y = choices / ntrials
    fig, axs = plt.subplots(1, 2, figsize=(10, 7), sharey=True, sharex=True)
    axs[0].plot(plot_ntd, np.mean(plot_p[..., 0], axis=1), color=green, linewidth=2)
    axs[0].plot(plot_ntd, np.mean(plot_p[..., 1], axis=1), color=blue, linewidth=2)
    axs[0].plot(plot_ntd, plot_p[:, ::50, 0], color=green, alpha=0.01)
    axs[0].plot(plot_ntd, plot_p[:, ::50, 1], color=blue, alpha=0.01)
    axs[0].plot(ntd, y[:, 0], "o", color=green)
    axs[0].plot(ntd, y[:, 1], "o", color=blue)
    axs[0].set_xlabel("NTD")
    axs[0].set_ylabel("P(T2 > T1)")
    axs[0].set_title("Cpntralateral")

    axs[1].plot(plot_ntd, np.mean(plot_p[..., 2], axis=1), color=green, linewidth=2)
    axs[1].plot(plot_ntd, np.mean(plot_p[..., 3], axis=1), color=blue, linewidth=2)
    axs[1].plot(plot_ntd, plot_p[:, ::50, 2], color=green, alpha=0.01)
    axs[1].plot(plot_ntd, plot_p[:, ::50, 3], color=blue, alpha=0.01)
    axs[1].plot(ntd, y[:, 2], "o", color=green)
    axs[1].plot(ntd, y[:, 3], "o", color=blue)
    axs[1].set_xlabel("NTD")
    axs[1].set_title("Ipsilateral")

    plt.figure()
    plt.plot(trace["pse"][:, 0], trace["pchoice"][:, 0], "o", color=green, alpha=0.1)
    plt.plot(trace["pse"][:, 1], trace["pchoice"][:, 1], "o", color=blue, alpha=0.1)