#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:16:51 2019

@author: lucianopaz
"""

from collections import namedtuple
import numpy as np
import pymc3 as pm
import theano
from theano import tensor as tt
from theano.tensor.signal.conv import conv2d
from pymc3.math import erf
from pymc3.theanof import take_along_axis
from data_interface import load_behavioral_data
from compress_pickle import load, dump
from tqdm import tqdm
import os
from behavioral_fit import Population, setup_drives, plot_predictions, shifts_plot


__all__ = [
    "Population",
    "leaky_integrator",
    "predicted_psychometric",
    "model_factory",
    "fit",
    "store_fits",
]


def convolve(in1, kernel, mode="full", **kwargs):
    if in1.ndim == 1:
        in1 = in1[None, :]
    if kernel.ndim == 1:
        kernel = kernel[None, :]
    return conv2d(in1, kernel, border_mode=mode)


def leaky_integrator(t, drive, tau=600.0):
    kernel = tt.exp(-t / tau) / tau
    if kernel.ndim == 1:
        kernel = kernel[None, :]

    LI = convolve(drive, kernel)[..., : kernel.shape[-1]]
    return LI


def predicted_psychometric(
    drives,
    T_ind,
    ref_drives,
    ref_ind,
    T_ref_ind,
    tau=6.4,
    background_mean=0.0,
    background_var=0.0,
    lapse_rate=0.1,
    lapse_bias=0.5,
    t=None,
):
    # Rescale parameters to their true domain
    if t is None:
        t = np.arange(drives.shape[-1], dtype=theano.config.floatX)

    # Precompute some factors
    back_drive = background_mean * tau * (1 - tt.exp(-t / tau))
    back_drive_var = background_var * 0.5 * tau * (1 - tt.exp(-2 * t / tau))

    # Get the reference against which to compare
    ref = leaky_integrator(t, ref_drives, tau) + back_drive
    ref_var = leaky_integrator(t, ref_drives, 0.5 * tau) + back_drive_var

    # Compute the drives
    li = leaky_integrator(t, drives, tau) + back_drive
    li_var = leaky_integrator(t, drives, 0.5 * tau) + back_drive_var

    # Compute dprime and psychometric
    if T_ind.ndim < li.ndim:
        li_T = take_along_axis(li, T_ind[..., None], axis=-1)[..., 0]
        li_var_T = take_along_axis(li_var, T_ind[..., None], axis=-1)[..., 0]
    else:
        li_T = take_along_axis(li, T_ind, axis=-1)
        li_var_T = take_along_axis(li_var, T_ind, axis=-1)
    dprime = (li_T - ref[ref_ind, T_ref_ind]) / tt.sqrt(
        2 * (li_var_T + ref_var[ref_ind, T_ref_ind])
    )
    pure_psychometric = 0.5 + 0.5 * erf(dprime)
    p = lapse_bias * lapse_rate + (1 - lapse_rate) * pure_psychometric
    return p


def model_factory(
    ntrials, choices, drives, T_ind, ref_drives, ref_ind, T_ref_ind, t=None
):
    with pm.Model() as model:
        tau = pm.Gamma("tau", mu=400, sd=200)
        background_mean = pm.Normal("background_mean", mu=10, sd=20)
        background_var = pm.Gamma("background_var", mu=400, sd=200)
        lapse_rate = pm.Beta("lapse_rate", alpha=2, beta=10)
        lapse_bias = pm.Beta("lapse_bias", alpha=2, beta=2)
        p = predicted_psychometric(
            drives=drives,
            T_ind=T_ind,
            ref_drives=ref_drives,
            ref_ind=ref_ind,
            T_ref_ind=T_ref_ind,
            tau=tau,
            background_mean=1000 * background_mean,
            background_var=1e5 * background_var,
            lapse_rate=lapse_rate,
            lapse_bias=lapse_bias,
            t=t,
        )
        pm.Binomial("choices", ntrials, p, observed=choices)
    return model


def fit(parameters_df, behavioral_data, population_kwargs=None, drives_options=None):
    # Extract key data
    T = behavioral_data.index.values.astype("float")
    ntrials = behavioral_data.xs("contra", axis=1).values[:, :2].T
    choices = behavioral_data.xs("contra", axis=1).values[:, 2:].T
    t = np.arange(np.max(T) + 1, dtype="float")

    # Build population
    if population_kwargs is None:
        population_kwargs = {}
    pop = Population(parameters_df=parameters_df, **population_kwargs)

    # Setupe population drive to LI
    if drives_options is None:
        drives_options = {}
    n_samples = drives_options.pop("n_samples", 5000)
    try:
        drives, contra_parameters, ipsi_parameters = load(
            os.path.join("simple_model_fits", "behavioral_drives.gz")
        )
    except Exception:
        drives, contra_parameters, ipsi_parameters = setup_drives(
            pop,
            T,
            t,
            kinds=["no_light", "contralateral_light"],
            n_samples=n_samples,
            **drives_options,
        )
        dump(
            (drives, contra_parameters, ipsi_parameters),
            os.path.join("simple_model_fits", "behavioral_drives.gz"),
        )
    ref_drives = drives[:, [1, 0], 3]
    T_ind = np.repeat(T[None, :], drives.shape[1], axis=0).astype("int")
    T_ref_ind = T[3].astype("int")
    ref_ind = np.arange(drives.shape[1], dtype="int")[:, None]
    T_ind, T_ref_ind, ref_ind = np.broadcast_arrays(T_ind, T_ref_ind, ref_ind)

    results = []
    Results = namedtuple("results", ["true_parameters"])
    TrueParameters = namedtuple(
        "true_parameters",
        [
            "tau",
            "background_mean",
            "background_var",
            "lapse_rate",
            "lapse_bias",
        ],
    )
    for drive, ref_drive in tqdm(zip(drives, ref_drives), total=len(drives)):
        with model_factory(
            ntrials, choices, drive, T_ind, ref_drive, ref_ind, T_ref_ind
        ):
            trace = pm.sample(draws=1, tune=1000, cores=1, chains=1)
            point = trace.point(0)
            # Using namedtuples to get a compatible result format with
            # behavioral_fit.py
            result = Results(
                TrueParameters(
                    point["tau"],
                    point["background_mean"],
                    point["background_var"],
                    point["lapse_rate"],
                    point["lapse_bias"]
                )
            )
            results.append(result)
    return results, contra_parameters, ipsi_parameters


def store_fits():
    try:
        parameters_df = load(os.path.join("simple_model_fits", "parameters.pkl"))
    except Exception:
        from simple_model import mean_parameters

        parameters_df = mean_parameters()
    behavioral_data = load_behavioral_data()
    res, contra_parameters, ipsi_parameters = fit(
        parameters_df,
        behavioral_data,
        population_kwargs=None,
        drives_options={"reps": 100},
    )
    dump(
        {
            "results": res,
            "contralateral_parameters": contra_parameters,
            "ipsilateral_parameters": ipsi_parameters,
        },
        os.path.join("simple_model_fits", "behavioral_parameters_mcmc.pkl"),
    )


if __name__ == "__main__":
    store_fits()
    plot_predictions(model_type="mcmc")
    shifts_plot(model_type="mcmc")
