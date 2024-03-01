#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 10:34:13 2019

@author: lucianopaz
"""

from autograd import numpy as np
import pandas as pd
from autograd.scipy.special import erf, expit, logit, gammaln
from autograd.scipy.signal import convolve
from autograd import grad
from scipy.optimize import minimize, basinhopping
from sklearn.mixture import BayesianGaussianMixture
from tqdm import tqdm
from data_interface import load_behavioral_data
from compress_pickle import load, dump
import os
from matplotlib import pyplot as plt
import matplotlib.colors as mc


__all__ = [
    "Population",
    "leaky_integrator",
    "rescale_parameters",
    "parameters_to_unbound_domain",
    "predicted_psychometric",
    "loss_function",
    "parameter_regularization",
    "regularized_loss_function",
    "setup_drives",
]


def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import colorsys

    try:
        c = mc.cnames[color]
    except Exception:
        c = color
    c = mc.to_rgb(c)
    c = colorsys.rgb_to_hls(*c)
    return colorsys.hls_to_rgb(c[0], max(0, min(1, 1 - amount * (1 - c[1]))), c[2])


def take_along_axis(arr, indices, axis):
    if axis < 0:
        axis = arr.ndim + axis
    if axis < 0 or axis >= arr.ndim:
        raise ValueError("axis is out of number of dimensions")
    arr_shape = arr.shape
    return arr[_make_along_axis_idx(arr_shape, indices, axis)]


def _make_along_axis_idx(arr_shape, indices, axis):
    # compute dimensions to iterate over
    if len(arr_shape) != indices.ndim:
        raise ValueError("`indices` and `arr` must have the same number of dimensions")
    shape_ones = (1,) * indices.ndim
    dest_dims = list(range(axis)) + [None] + list(range(axis + 1, indices.ndim))

    # build a fancy index, consisting of orthogonal aranges, with the
    # requested index inserted at the right location
    fancy_index = []
    for dim, n in zip(dest_dims, arr_shape):
        if dim is None:
            fancy_index.append(indices)
        else:
            ind_shape = shape_ones[:dim] + (-1,) + shape_ones[dim + 1 :]
            fancy_index.append(np.arange(n).reshape(ind_shape))

    return tuple(fancy_index)


def _richards(x, nu, lam):
    return lam * (1 - 1 / (1 + np.exp(10 * x)) ** (1 / nu)) + 1e-4


def binomial_logp(k, n, p):
    return (
        gammaln(n)
        - gammaln(k)
        - gammaln(n - k)
        + k * np.log(p)
        + (n - k) * np.log1p(-p)
    )


class Population:
    def __init__(self, parameters_df, n_components=5, gmm_on_parameters=None):
        gmm = BayesianGaussianMixture(
            n_components=n_components,
            covariance_type="full",
            n_init=10,
            weight_concentration_prior_type="dirichlet_process",
            weight_concentration_prior=1e-6,
        )
        self.parameters = parameters_df.xs("mean", axis=1, level=1)
        self.set_gmm_parameters(gmm_on_parameters)
        self.gmm = gmm.fit(self.gmm_parameters)

    def set_gmm_parameters(self, gmm_on_parameters=None):
        if gmm_on_parameters is None:
            gmm_on_parameters = ["nu", "lam", "I0", "I1", "I1p", "I2", "I2p", "I3p"]
        self.gmm_parameters = self.parameters[gmm_on_parameters]
        self.gmm_columns = gmm_on_parameters
        self.median_columns = [c for c in self.parameters if c not in gmm_on_parameters]
        self.medians = self.parameters[self.median_columns].median(axis=0)
        if "mech_onset" in self.medians:
            self.medians["mech_onset"] = 0
        if "mech_offset" in self.medians:
            self.medians["mech_offset"] = 0

    def sample(self, n_samples=1):
        samples, label = self.gmm.sample(n_samples)
        output = pd.DataFrame(data=samples, columns=self.gmm_columns)
        output = output.assign(
            **{k: v * np.ones(n_samples) for k, v in self.medians.items()}
        )
        return output.assign(label=label)

    def response_rate(self, t, T, mech=True, light=False, samples=None, reduce=True):
        if samples is None:
            samples = self.sample()[0]
        onset = 0.0
        nu = samples["nu"].values
        lam = samples["lam"].values
        I0 = samples["I0"].values
        I1 = samples["I1"].values
        I1p = samples["I1p"].values
        I2 = samples["I2"].values
        I2p = samples["I2p"].values
        I3p = samples["I3p"].values
        tau = samples["tau"].values
        tau_opto = samples["tau_opto"].values
        tau_inh = samples["tau_inh"].values

        offset_t = t - T
        w = expit(1 * (t - onset)) * expit(-1 * (offset_t - onset))
        if light:
            w2 = expit(1 * (offset_t - onset))
        else:
            w2 = 0.0
        if mech and light:
            input_current = I0 + (
                w
                * (
                    I1
                    + (I1p - I1) * np.exp(-(t - onset) / tau)
                    + I2
                    + (I2p - I2) * np.exp(-(t - onset) / tau_opto)
                )
                + w2 * I3p * np.exp(-(offset_t - onset) / tau_inh)
            )
        elif mech:
            input_current = I0 + (w * (I1 + (I1p - I1) * np.exp(-(t - onset) / tau)))
        elif light:
            input_current = I0 + (
                w * (I2 + (I2p - I2) * np.exp(-(t - onset) / tau_opto))
                + w2 * I3p * np.exp(-(offset_t - onset) / tau_inh)
            )
        else:
            I0, t = np.broadcast_arrays(I0, t)
            input_current = I0
        rate = _richards(input_current, nu, lam)
        if reduce:
            rate = np.sum(rate, axis=-1)
        return rate

    def no_light_drive(
        self, t, T, contralateral_parameters=None, ipsilateral_parameters=None
    ):
        if contralateral_parameters is None:
            contralateral_parameters = self.sample(n_samples=1000)
        if ipsilateral_parameters is None:
            ipsilateral_parameters = self.sample(n_samples=1000)
        contra_rate = self.response_rate(
            t=t[..., None],
            T=T[..., None],
            mech=True,
            light=False,
            samples=contralateral_parameters,
            reduce=True,
        )
        ipsi_rate = self.response_rate(
            t=t[..., None],
            T=T[..., None],
            mech=False,
            light=False,
            samples=ipsilateral_parameters,
            reduce=True,
        )
        return contra_rate + ipsi_rate

    def contralateral_light_drive(
        self, t, T, contralateral_parameters=None, ipsilateral_parameters=None
    ):
        if contralateral_parameters is None:
            contralateral_parameters = self.sample(n_samples=1000)
        if ipsilateral_parameters is None:
            ipsilateral_parameters = self.sample(n_samples=1000)
        contra_rate = self.response_rate(
            t=t[..., None],
            T=T[..., None],
            mech=True,
            light=True,
            samples=contralateral_parameters,
            reduce=True,
        )
        ipsi_rate = self.response_rate(
            t=t[..., None],
            T=T[..., None],
            mech=False,
            light=False,
            samples=ipsilateral_parameters,
            reduce=True,
        )
        return contra_rate + ipsi_rate

    def ipsilateral_light_drive(
        self, t, T, contralateral_parameters=None, ipsilateral_parameters=None
    ):
        if contralateral_parameters is None:
            contralateral_parameters = self.sample(n_samples=1000)
        if ipsilateral_parameters is None:
            ipsilateral_parameters = self.sample(n_samples=1000)
        contra_rate = self.response_rate(
            t=t[..., None],
            T=T[..., None],
            mech=False,
            light=True,
            samples=contralateral_parameters,
            reduce=True,
        )
        ipsi_rate = self.response_rate(
            t=t[..., None],
            T=T[..., None],
            mech=True,
            light=False,
            samples=ipsilateral_parameters,
            reduce=True,
        )
        return contra_rate + ipsi_rate


def leaky_integrator(t, drive, tau=600.0):
    kernel = np.exp(-t / tau) / tau
    if kernel.ndim < drive.ndim:
        kernel = kernel[(np.newaxis,) * (drive.ndim - kernel.ndim) + (slice(None),)]
    LI = convolve(drive, kernel)[..., : kernel.shape[-1]]
    return LI


def rescale_parameters(
    tau_log, background_mean, background_var_log, lapse_rate_logodds, lapse_bias_logodds
):
    return (
        np.exp(tau_log),
        background_mean,
        np.exp(background_var_log),
        expit(lapse_rate_logodds),
        expit(lapse_bias_logodds),
    )


def parameters_to_unbound_domain(
    tau, background_mean, background_var, lapse_rate, lapse_bias
):
    return (
        np.log(tau),
        background_mean,
        np.log(background_var),
        logit(lapse_rate),
        logit(lapse_bias),
    )


def predicted_psychometric(
    drives,
    T_ind,
    ref_drives,
    ref_ind,
    T_ref_ind,
    tau_log=6.4,
    background_mean=0,
    background_var_log=-np.inf,
    lapse_rate_logodds=-np.inf,
    lapse_bias_logodds=0.0,
    t=None,
):
    # Rescale parameters to their true domain
    tau, background_mean, background_var, lapse_rate, lapse_bias = rescale_parameters(
        tau_log,
        background_mean,
        background_var_log,
        lapse_rate_logodds,
        lapse_bias_logodds,
    )
    if t is None:
        t = np.arange(drives.shape[-1], dtype=np.float)

    # Precompute some factors
    back_drive = background_mean * tau * (1 - np.exp(-t / tau))
    back_drive_var = background_var * 0.5 * tau * (1 - np.exp(-2 * t / tau))

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
    if ref.ndim > 2:
        advanced_ind = np.arange(len(ref)).reshape((-1, 1, 1))
        dprime = (li_T - ref[advanced_ind, ref_ind, T_ref_ind]) / np.sqrt(
            2 * (li_var_T + ref_var[advanced_ind, ref_ind, T_ref_ind])
        )
    else:
        dprime = (li_T - ref[ref_ind, T_ref_ind]) / np.sqrt(
            2 * (li_var_T + ref_var[ref_ind, T_ref_ind])
        )
    pure_psychometric = 0.5 + 0.5 * erf(dprime)
    p = lapse_bias * lapse_rate + (1 - lapse_rate) * pure_psychometric
    return p


def loss_function(
    parameter_array, ntrials, choices, drives, T_ind, ref_drives, ref_ind, T_ref_ind, t
):
    tau_log, background_mean, background_var_log, lapse_rate_logodds, lapse_bias_logodds = (
        parameter_array
    )
    p = predicted_psychometric(
        drives,
        T_ind,
        ref_drives,
        ref_ind,
        T_ref_ind,
        tau_log=tau_log,
        background_mean=background_mean,
        background_var_log=background_var_log,
        lapse_rate_logodds=lapse_rate_logodds,
        lapse_bias_logodds=lapse_bias_logodds,
        t=t,
    )
    return -np.sum(binomial_logp(choices, ntrials, p))


def parameter_regularization(parameter_array, w, loc):
    p = np.array(rescale_parameters(*parameter_array))
    return (
        np.sum((w[:3] * (p[:3] - loc[:3])) ** 2) +  # Gaussian
        np.sum(
            (loc[3:] - 1) * np.log(p[3:]) +  # Beta
            (w[3:] - 1) * np.log(1 - p[3:]) -
            (gammaln(loc[3:]) + gammaln(w[3:]) - gammaln(loc[3:] + w[3:]))
        )
    )


def regularized_loss_function(
    parameter_array,
    w,
    loc,
    ntrials,
    choices,
    drives,
    T_ind,
    ref_drives,
    ref_ind,
    T_ref_ind,
    t,
):
    return loss_function(
        parameter_array,
        ntrials,
        choices,
        drives,
        T_ind,
        ref_drives,
        ref_ind,
        T_ref_ind,
        t,
    ) + parameter_regularization(parameter_array, w, loc)


def setup_drives(
    population,
    T,
    t,
    n_samples,
    reps=200,
    kinds=["no_light", "contralateral_light"],
    contralateral_parameters=None,
    ipsilateral_parameters=None,
):
    _kind_map = {
        "no_light": "no_light_drive",
        "contralateral_light": "contralateral_light_drive",
        "ipsilateral_light": "ipsilateral_light_drive",
    }
    drives = np.empty((reps, len(kinds), len(T), len(t)), dtype="float")
    cp = []
    ip = []
    for i, rep in enumerate(tqdm(range(reps), desc="Setting up population drives")):
        if contralateral_parameters is None:
            contra_parameters = population.sample(n_samples=n_samples)
            cp.append(contra_parameters)
        else:
            contra_parameters = contralateral_parameters[i]
        if ipsilateral_parameters is None:
            ipsi_parameters = population.sample(n_samples=n_samples)
            ip.append(ipsi_parameters)
        else:
            ipsi_parameters = ipsilateral_parameters[i]
        for i, kind in enumerate(kinds):
            method = _kind_map[kind]
            drive = getattr(population, method)(
                t,
                T[..., None],
                contralateral_parameters=contra_parameters,
                ipsilateral_parameters=ipsi_parameters,
            )
            drives[rep, i] = drive
    return drives, cp, ip


def fit(
    parameters_df,
    behavioral_data,
    method="L-BFGS-B",
    regularize_parameters=True,
    population_kwargs=None,
    drives_options=None,
    minimize_options=None,
):
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
    if drives_options is not None:
        reps = drives_options.get("reps", None)
    else:
        reps = None
    if reps is not None and len(drives) > reps:
        drives = drives[:reps]
        contra_parameters = contra_parameters[:reps]
        ipsi_parameters = ipsi_parameters[:reps]
    ref_drives = drives[:, [1, 0], 3]
    T_ind = np.repeat(T[None, :], drives.shape[1], axis=0).astype("int")
    T_ref_ind = T[3].astype("int")
    ref_ind = np.arange(drives.shape[1], dtype="int")[:, None]
    T_ind, T_ref_ind, ref_ind, _ = np.broadcast_arrays(
        T_ind, T_ref_ind, ref_ind, drives[..., 0]
    )

    if regularize_parameters:
        fun = regularized_loss_function
        w = np.array([1e-2, 0, 0, 1, 1])
        loc = np.array([6e3, 0, 1500, 1, 1])
        args = (
            w,
            loc,
            ntrials,
            choices,
            drives,
            T_ind,
            ref_drives,
            ref_ind,
            T_ref_ind,
            t,
        )
    else:
        fun = loss_function
        args = (ntrials, choices, drives, T_ind, ref_drives, ref_ind, T_ref_ind, t)
#    x0 = np.array(
#        parameters_to_unbound_domain(600.0, 7, 1000, 0.15, 0.5)
#    )
    x0 = np.array(
        parameters_to_unbound_domain(6e3, 7, 1000, 0.15, 0.5)
    )
    jac = grad(fun, 0)
    results = []
    if minimize_options is None:
        minimize_options = dict()
    minimize_options.setdefault("disp", True)
    res = minimize(
        fun, x0=x0, args=args, method=method, jac=jac, options=minimize_options
    )
    res.true_parameters = np.array(rescale_parameters(*res.x))
    results = res
    return results, contra_parameters, ipsi_parameters


class SimpleLogistic:
    @staticmethod
    def rescale(beta, pse, lapse_rate_logodds, lapse_bias_logodds):
        return beta, pse, expit(lapse_rate_logodds), expit(lapse_bias_logodds)

    @staticmethod
    def invert(beta, pse, lapse_rate, lapse_bias):
        return beta, pse, logit(lapse_rate), logit(lapse_bias)

    @staticmethod
    def pred(params, x):
        beta, pse, lapse_rate, lapse_bias = SimpleLogistic.rescale(*params)
        p = lapse_rate * lapse_bias + (1 - lapse_rate) * expit(beta * (x - pse))
        return p

    @staticmethod
    def loss(params, x, choices, ntrials):
        p = SimpleLogistic.pred(params, x)
        return -np.sum(binomial_logp(choices, ntrials, p))

    def fit(
        self, behavioral_data, side="contra", method="L-BFGS-B", minimize_options=None
    ):
        # Extract key data
        T = behavioral_data.index.values.astype("float")
        ntd = (T - 334) / (T + 334)
        ntrials = behavioral_data.xs(side, axis=1).values[:, :2].T
        choices = behavioral_data.xs(side, axis=1).values[:, 2:].T

        x0 = np.array(self.invert(10, 0.0, 0.05, 0.5))
        fun = self.loss
        jac = grad(fun, 0)
        results = []
        for i in range(len(ntrials)):
            args = (ntd, choices[i], ntrials[i])
            res = minimize(
                fun, x0=x0, args=args, method=method, jac=jac, options=minimize_options
            )
            res.true_parameters = np.array(self.rescale(*res.x))
            results.append(res)
        return results

    def posterior_predictive(self, fit_results, ntd):
        return np.array([self.pred(res.x, ntd) for res in fit_results])


def store_fits(
    method="L-BFGS-B",
    regularize_parameters=True,
    population_kwargs=None,
    drives_options=None,
    minimize_options=None,
):
    try:
        parameters_df = load(os.path.join("simple_model_fits", "parameters.pkl"))
    except Exception:
        from simple_model import mean_parameters

        parameters_df = mean_parameters()
    behavioral_data = load_behavioral_data()
    res, contra_parameters, ipsi_parameters = fit(
        parameters_df,
        behavioral_data,
        method=method,
        regularize_parameters=regularize_parameters,
        population_kwargs=population_kwargs,
        drives_options=drives_options,
        minimize_options=minimize_options,
    )
    dump(
        {
            "results": res,
            "contralateral_parameters": contra_parameters,
            "ipsilateral_parameters": ipsi_parameters,
        },
        os.path.join("simple_model_fits", "behavioral_parameters.pkl"),
    )


def load_predictions(override=False, model_type="mle"):
    extension = "" if model_type == "mle" else "_mcmc"
    if not override:
        try:
            return load(
                os.path.join(
                    "simple_model_fits", "predicted_psychometric{}.pkl".format(extension)
                )
            )
        except Exception:
            if model_type != "mle":
                raise
    elif model_type != "mle":
        raise TypeError(
            "{}.load_predictions cannot compute the predicted psychometrics "
            "for model_types different from 'mle'. Please run load_predictions "
            "from the appropriate module, such as logistic_behavioral_fit or "
            "behavioral_fit_mcmc.".format(str(__file__)[:-3])
        )
    parameters_df = load(os.path.join("simple_model_fits", "parameters.pkl"))
    pop = Population(parameters_df=parameters_df)
    temp = load(
        os.path.join(
            "simple_model_fits", "behavioral_parameters{}.pkl".format(extension)
        )
    )
    results = temp["results"]
    contra_parameters = temp["contralateral_parameters"]
    ipsi_parameters = temp["ipsilateral_parameters"]
    T = np.array([780])
    t = np.arange(T[0], dtype="float")
    drives = setup_drives(
        pop,
        T,
        t,
        n_samples=5000,
        reps=len(contra_parameters),
        kinds=["no_light", "contralateral_light", "ipsilateral_light"],
        contralateral_parameters=contra_parameters,
        ipsilateral_parameters=ipsi_parameters,
    )[0]
    ref_drives = drives[:, :, 0]
    drives = drives[:, [0, 1, 0, 2], 0]
    T_ind = np.arange(143, 780, dtype="int")
    T_ref_ind = 334
    ref_ind = np.array([1, 0, 2, 0])[:, None]

    if not isinstance(results, list):
        T_ind, T_ref_ind, ref_ind, _ = np.broadcast_arrays(
            T_ind, T_ref_ind, ref_ind, np.empty(drives.shape[:-1] + (1,))
        )
        psychometric = predicted_psychometric(
            drives, T_ind, ref_drives, ref_ind, T_ref_ind, *results.x, t=t
        )
        Ts = T_ind[0, 0]
    else:
        T_ind, T_ref_ind, ref_ind = np.broadcast_arrays(
            T_ind, T_ref_ind, ref_ind
        )
        Ts = T_ind[0]
        psychometric = np.empty((len(drives), 4, T_ind.shape[-1]))
        for i, (drive, ref_drive, res) in tqdm(
            enumerate(zip(drives, ref_drives, results)), total=len(drives)
        ):
            p = predicted_psychometric(
                drive, T_ind, ref_drive, ref_ind, T_ref_ind, *res.x, t=t
            )
            psychometric[i] = p
    dump(
        (Ts, psychometric),
        os.path.join(
            "simple_model_fits", "predicted_psychometric{}.pkl".format(extension)
        ),
    )
    return Ts, psychometric


def plot_predictions(override=False, model_type="mle"):
    T, psychometric = load_predictions(override=override, model_type=model_type)
    # We change the lapse rates of the ipsilateral condition
    # First we remove the lapses inferred on the contralateral side
    if model_type == "mle":
        temp = load(os.path.join("simple_model_fits", "behavioral_parameters.pkl"))
    elif model_type == "mcmc":
        temp = load(os.path.join("simple_model_fits", "behavioral_parameters_mcmc.pkl"))
    results = temp["results"]
    cl, cb = results.true_parameters[-2:]
    psychometric[:, 2:] = (psychometric[:, 2:] - cl * cb)/(1 - cl)
    # Now we use the lapse rates that were inferred separately for the ipsi side
    il = 0.44
    ib = 0.48
    psychometric[:, 2:] = il * ib + (1 - il) * psychometric[:, 2:]
    psychometric *= 100

    behavioral_data = load_behavioral_data()
    contra = behavioral_data.xs("contra", level=0, axis=1)
    ipsi = behavioral_data.xs("ipsi", level=0, axis=1)
    subj_ntd = ((contra.index - 334) / (contra.index + 334)).values

    blue = np.array([62, 105, 178]) / 255
    green = np.array([106, 189, 69]) / 255
    ntd = (T - 334) / (T + 334)
    with plt.style.context(["default", "leaky_neurons_figure_style"]):
        f, axs = plt.subplots(1, 2, sharex=True, sharey=False, figsize=(12, 7))
        axs[0].axvline(0, color="gray")
        axs[0].axhline(50, color="gray")
        axs[1].axvline(0, color="gray")
        axs[1].axhline(50, color="gray")
        axs[0].plot(ntd, psychometric[:, 0].T, color=green, alpha=0.05, linewidth=1)
        axs[0].plot(ntd, psychometric[:, 1].T, color=blue, alpha=0.05, linewidth=1)
#        axs[0].plot(ntd, np.mean(psychometric[:, 0], axis=0), color=green)
#        axs[0].plot(ntd, np.mean(psychometric[:, 1], axis=0), color=blue)
        axs[0].plot(
            subj_ntd,
            100 * contra["nchoice2_light1"] / contra["ntrials_light1"],
            "o",
            color=green,
        )
        axs[0].plot(
            subj_ntd,
            100 * contra["nchoice2_light2"] / contra["ntrials_light2"],
            "o",
            color=blue,
        )
        axs[1].plot(ntd, psychometric[:, 2].T, color=green, alpha=0.05, linewidth=1)
        axs[1].plot(ntd, psychometric[:, 3].T, color=blue, alpha=0.05, linewidth=1)
#        axs[1].plot(ntd, np.mean(psychometric[:, 2], axis=0), color=green)
#        axs[1].plot(ntd, np.mean(psychometric[:, 3], axis=0), color=blue)
        axs[1].plot(
            subj_ntd,
            100 * ipsi["nchoice2_light1"] / ipsi["ntrials_light1"],
            "o",
            color=green,
        )
        axs[1].plot(
            subj_ntd,
            100 * ipsi["nchoice2_light2"] / ipsi["ntrials_light2"],
            "o",
            color=blue,
        )
        #        axs[0].set_title("Contralateral barrel\ncortex activation")
        axs[0].set_yticks(np.arange(0, 101, 20))
        axs[0].set_xticks([-0.2, 0, 0.2])
        axs[0].set_ylabel("percent called T2>T1")
        axs[0].set_xlabel("Duration difference")
        axs[0].set_ylim([0, 100])
        axs[0].set_xlim([-0.36, 0.36])

        #        axs[1].set_title("Ipsilateral barrel\ncortex activation", color="red")
        axs[1].set_ylabel("percent called T2>T1")
        axs[1].set_xlabel("Duration difference")
        axs[1].set_ylim([0, 100])
        axs[1].set_xlim([-0.37, 0.37])
        plt.tight_layout()
        extension = "" if model_type == "mle" else "_mcmc"
        plt.savefig(
            os.path.join(
                "simple_model_fits", "predicted_psychometric{}.pdf".format(extension)
            ),
            bbox_inches="tight",
        )
        plt.savefig(
            os.path.join(
                "simple_model_fits", "predicted_psychometric{}.svg".format(extension)
            ),
            bbox_inches="tight",
        )


def _pse(x, p):
    return np.interp(0.5, p, x)


pse = np.vectorize(_pse, signature="(i),(i)->()")


def predicted_shifts(model_type="mle"):
    behavioral_data = load_behavioral_data()
    ntrials = np.concatenate(
        (
            behavioral_data.xs("contra", axis=1).values[:, :2],
            behavioral_data.xs("ipsi", axis=1).values[:, :2],
        ),
        axis=1,
    )
    T_exp = behavioral_data.index.values.astype(int)
    if model_type == "mle":
        temp = load(os.path.join("simple_model_fits", "behavioral_parameters.pkl"))
    elif model_type == "mcmc":
        temp = load(os.path.join("simple_model_fits", "behavioral_parameters_mcmc.pkl"))
    results = temp["results"]
    if not isinstance(results, list):
        lapses = results.true_parameters[-2:]
        lapse_rate = lapses[0]
        lapse_bias = lapses[1]
    else:
        lapses = np.array([r.true_parameters[-2:] for r in results])
        lapse_rate = lapses[:, 0][..., None, None]
        lapse_bias = lapses[:, 1][..., None, None]
    T, psychometric = load_predictions()
    T_inds = T_exp - T[0]
    ntd = (T - 334) / (T + 334)
    PSE = pse(ntd, psychometric)
    inflection_point = pse(ntd, (psychometric - lapse_rate * lapse_bias) / (1 - lapse_rate))
    behavioral_bias = 100 * (
        np.sum(
            take_along_axis(psychometric, T_inds[None, None], axis=-1) * ntrials.T,
            axis=-1,
        )
        / np.sum(ntrials, axis=0)
    )

    PSE = PSE[:, ::2] - PSE[:, 1::2]
    inflection_point = inflection_point[:, ::2] - inflection_point[:, 1::2]
    behavioral_bias = behavioral_bias[:, 1::2] - behavioral_bias[:, ::2]
    return pd.DataFrame(
        data=np.concatenate((inflection_point, PSE, behavioral_bias), axis=1),
        columns=[
            "contralateral_inflection_point_shift",
            "ipsilateral_inflection_point_shift",
            "contralateral_pse_shift",
            "ipsilateral_pse_shift",
            "contralateral_behavioral_bias",
            "ipsilateral_behavioral_bias",
        ],
    )


def shifts_plot(model_type="mle"):
    from logistic_behavioral_fit import load_behavioral_fit

    extension = "" if model_type == "mle" else "_{}".format(model_type)
    model_shifts = predicted_shifts(model_type=model_type)
    data_shifts = load_behavioral_fit()["trace"]
    np.random.seed(123456789)
    data_shifts = (
        data_shifts.drop(
            columns=[c for c in data_shifts.columns if c not in model_shifts]
        )
        .sample(len(model_shifts))
        .reset_index(drop=True)
    )
    dfs = []
    for column in model_shifts.columns:
        temp = column.split("_")
        side = temp[0].capitalize()
        shift_type = " ".join(temp[1:]).capitalize()
        model_col = model_shifts[column].to_frame().rename(columns=lambda x: "shift")
        data_col = data_shifts[column].to_frame().rename(columns=lambda x: "shift")
        dfs.append(model_col.assign(side=side, shift_type=shift_type, kind="model"))
        dfs.append(data_col.assign(side=side, shift_type=shift_type, kind="data"))
    df = pd.concat(dfs, ignore_index=True)
    base_colors = ["black", "red"]
    centers = np.arange(2)
    scat_markers = ["o", "P"]
    mean_markers = ["o", "P"]
    for shift_type in df.shift_type.unique():
        fig, ax = plt.subplots(1, 1, figsize=(3, 7))
        sides = ["Contralateral", "Ipsilateral"]
        for center, side, base_color in zip(centers, sides, base_colors):
            for displacement, (kind, scat_marker, marker) in enumerate(
                zip(["model", "data"], scat_markers, mean_markers)
            ):
                if displacement == 0:
                    color = lighten_color(base_color, 0.65)
                else:
                    color = base_color
                if center == centers[0]:
                    label = kind
                else:
                    label = None
                df2 = df.query(
                    "shift_type == @shift_type and "
                    "side == @side and "
                    "kind == @kind"
                )
                x = center + 0.4 * (displacement - 0.5)
                jitter = 0.04
                y = df2["shift"].values
                ax.plot(
                    x + jitter * np.random.randn(len(y)),
                    y,
                    linestyle="",
                    marker=scat_marker,
                    color=color,
                    alpha=0.1,
                )
                ax.errorbar(
                    [x],
                    [np.mean(y)],
                    [np.std(y)],
                    linestyle="",
                    marker=marker,
                    markersize=12,
                    color=color,
                    label=label,
                    capsize=3,
                )
        ax.legend()
        ax.set_xticks(centers)
        ax.set_xticklabels(sides)
        ax.set_ylabel(shift_type)
        for color, xticklabel in zip (base_colors, ax.get_xticklabels()):
            xticklabel.set_color(color)
        plt.savefig(
            os.path.join(
                "simple_model_fits",
                "predicted_{}{}{}.pdf".format(
                    shift_type.lower().replace(" ", "_"),
                    "" if shift_type.endswith("shift") else "_shift",
                    extension
                )
            ),
            bbox_inches="tight",
        )


if __name__ == "__main__":
    store_fits(regularize_parameters=True)
    plot_predictions(override=True)
    shifts_plot()
