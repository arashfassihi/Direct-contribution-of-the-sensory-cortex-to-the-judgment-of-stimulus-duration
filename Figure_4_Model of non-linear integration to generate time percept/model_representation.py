#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:06:15 2019

@author: lucianopaz
"""

import numpy as np
from theano import tensor as tt
import pymc3 as pm
from pymc3.distributions import transforms as tr
from data_interface import load_fittable_data, DT
from utils import alpha_function, beta_function


dt = DT
stacked_df = load_fittable_data(query="nid < 10", light_on_shift=dt)
mfr = stacked_df.groupby("nid").apply(
    lambda x: np.mean((x["psth"] * x["ntrials"])[(x["t"] >= 0) & (x["offset_t"] <= 0)])
    / dt
)
mre = stacked_df.groupby("nid")["mech_resp"].apply(lambda x: np.unique(x)[0])
ore = stacked_df.groupby("nid")["opto_resp"].apply(lambda x: np.unique(x)[0])
valid = (mfr > 0.5) | mre | ore
valid_nid = list(valid.index[valid])
stacked_df = stacked_df.query("nid in @valid_nid")

t = stacked_df["t"].values
offset_t = stacked_df["offset_t"].values
psth = stacked_df["psth"].values
n = stacked_df["ntrials"].values
nid = stacked_df["nid"].values
light = stacked_df["light"].values

mech_stim_on = np.logical_and(t >= 0, offset_t < 0)
opto_stim_on = np.logical_and(np.logical_and(t >= -10, offset_t < -10), light)
opto_stim_offset = np.logical_and(offset_t >= -10, light)

unid, nid_inds = np.unique(nid, return_inverse=True)
ulight, light_inds = np.unique(light, return_inverse=True)


def mixture_builder(
    name_root,
    n=100,
    npop=2,
    alpha=None,
    p_n=None,
    mu_pop=None,
    alpha_kwargs=None,
    mu_kwargs=None,
    set_default_mu_order=True,
    mixture_kwargs=None,
):
    if alpha_kwargs is None:
        alpha_kwargs = {}
    alpha_kwargs.setdefault("alpha", 1.0)
    alpha_kwargs.setdefault("beta", 1.0)

    if mu_kwargs is None:
        mu_kwargs = {}
    mu_kwargs.setdefault("mu", 0.0)
    mu_kwargs.setdefault("sigma", 0.5)
    if set_default_mu_order:
        mu_kwargs.setdefault("transform", tr.ordered)
        mu_kwargs.setdefault("testval", np.linspace(0, 1, npop))

    if mixture_kwargs is None:
        mixture_kwargs = {}
    mixture_kwargs.setdefault("sigma", 1.0)

    if alpha is None:
        alpha = pm.Gamma("{}_pop_w".format(name_root), shape=npop, **alpha_kwargs)
    if p_n is None:
        p_n = pm.Dirichlet("{}_n_prob".format(name_root), alpha, shape=(n, npop))
    if mu_pop is None:
        mu_pop = pm.Normal("{}_mu".format(name_root), shape=npop, **mu_kwargs)
    x_i = pm.NormalMixture(
        name_root, w=p_n, mu=mu_pop, shape=n, comp_shape=(n, npop), **mixture_kwargs
    )
    return alpha, p_n, mu_pop, x_i


with pm.Model() as model:
    mixed_pop_model = model
    npop = 2
    lam_pop = pm.Normal("lam_pop", mu=3.4, sigma=1.0)
    lam = pm.Gamma("lam", mu=tt.exp(lam_pop), sigma=5, shape=len(unid))

    alpha1, p_n1, _, I1 = mixture_builder("I1", len(unid))
    I1p = mixture_builder("I1p", len(unid), npop, alpha=alpha1, p_n=p_n1)[-1]
    alpha2, p_n2, _, I2 = mixture_builder("I2", len(unid), npop)
    I2p = mixture_builder("I2p", len(unid), npop, alpha=alpha2, p_n=p_n2)[-1]
    I3p = mixture_builder("I3p", len(unid), npop, alpha=alpha2, p_n=p_n2)[-1]
    I0 = pm.Normal("I0", mu=0, sigma=1.0, shape=len(unid))

    tau = pm.Gamma("tau", mu=50, sigma=1)
    tau_opto = pm.Gamma("tau_opto", mu=20, sigma=1)
    tau_inh = pm.Gamma("tau_inh", mu=20, sigma=1)
    light_lag = np.zeros(len(unid))
    opto_stim_on = tt.switch(
        tt.and_(
            tt.and_(
                tt.ge(t, light_lag[nid_inds]), tt.lt(offset_t, light_lag[nid_inds])
            ),
            light == 1,
        ),
        1.0,
        0.0,
    )
    opto_stim_offset = tt.switch(
        tt.and_(tt.ge(offset_t, light_lag[nid_inds]), light == 1), 1.0, 0.0
    )
    inds = nid_inds
    input_current = (
        I0[nid_inds]
        + mech_stim_on * (I1[inds] + (I1p - I1)[inds] * alpha_function(t, tau))
        + light
        * opto_stim_on
        * (
            I2[inds]
            + (I2p - I2)[inds] * beta_function(t - light_lag[nid_inds], 10, tau_opto)
        )
        + opto_stim_offset
        * (I3p[inds] * alpha_function(offset_t - light_lag[nid_inds], tau_inh))
    )
    io_func = tt.nnet.sigmoid(input_current) * lam[nid_inds]
    mu_t = (io_func + 1) * n * 1e-3 * dt
    obs_psth = pm.Poisson("psth", mu=mu_t, observed=psth)

gviz = pm.model_to_graphviz(model)
