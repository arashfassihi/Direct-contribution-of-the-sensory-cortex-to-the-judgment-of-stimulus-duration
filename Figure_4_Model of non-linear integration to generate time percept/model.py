#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 15:23:16 2018

@author: lucianopaz
"""

import os
import re
import numpy as np
from theano import tensor as tt
import pymc3 as pm
from pymc3.distributions import transforms as tr
from matplotlib import pyplot as plt
from compress_pickle import dump, load
from data_interface import load_fittable_data, DT, merge_psths
from utils import alpha_function, beta_function
from plot_utils import plot_ppc, crude_pair_plot
import pandas as pd


dt = DT
stacked_df = load_fittable_data(
    query='nid < 10',
    light_on_shift=0,
    filter_non_responsive=True,
)
jitters = []
nid = []
for f in [fs for fs in os.listdir('non_responsive') if fs.startswith('jitter_')]:
    jitters.append(load(os.path.join('non_responsive', f))[0])
    nid.append(int(re.search('[0-9]+', f).group()))
jitters = np.array(jitters)
jitters = pd.DataFrame(data=jitters,
                       index=nid,
                       columns=['light_off', 'light_on']).sort_index()
for nid, row in jitters.iterrows():
    for light, jitter in enumerate(row):
        inds = np.logical_and(stacked_df.nid == nid,
                              stacked_df.light== light)
        stacked_df.loc[inds, ['t', 'offset_t']] -= jitter
#area_1 = [46, 139]
#area_2 = [139, 228]
#inds1 = stacked_df['nid'] >= area_1[0] & stacked_df['nid'] < area_1[1]
#inds2 = stacked_df['nid'] >= area_2[0] & stacked_df['nid'] < area_2[1]
#stacked_df['t'][inds1] -= 20
#stacked_df['offset_t'][inds1] -= 20
#stacked_df['t'][inds2] -= 10
#stacked_df['offset_t'][inds2] -= 10
mfr = stacked_df.groupby('nid').apply(
    lambda x: np.mean((x['psth'] *
                       x['ntrials'])[(x['t'] >=0) &
                                     (x['offset_t'] <= 0)]) / dt
)
mre = stacked_df.groupby('nid')['mech_resp'].apply(lambda x: np.unique(x)[0])
ore = stacked_df.groupby('nid')['opto_resp'].apply(lambda x: np.unique(x)[0])
valid = (mfr > 0.5) | mre | ore
valid_nid = list(valid.index[valid])
stacked_df = stacked_df.query('nid in @valid_nid')

#stacked_df = load_fittable_data(light_on_shift=dt)
#stacked_df = merge_psths(stacked_df).reset_index()
#stacked_df = stacked_df.assign(nid=-1)
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

# with pm.Model('constant') as const_model:
#    background = pm.HalfCauchy('background', beta=10.)
#    mu = pm.HalfNormal('mu', sigma=100.,
#                       shape=psth.shape[1] // len(unid))
#    mu_t = (mech_stim_on[:, :2] * mu * n[:, :2] * 1e-3 * dt +
#            background)
#    obs_psth = pm.Poisson('psth', mu=mu_t,
#                          observed=psth[:, :2])
#    const_trace = PyMC3Sampler()(2000, tune=1000)
#    pm.traceplot(const_trace)
#    crude_pair_plot(const_trace, const_model)
#    ppc = pm.sample_posterior_predictive(const_trace)
#    plt.figure()
#    plot_ppc(t, psth, ppc[const_model.name_for('psth')], 500)

# with pm.Model('const_link') as model:
#    const_link_model = model
#    lam = pm.Gamma('lam', mu=30, sigma=5)
#    background = pm.Normal('background', mu=0, sigma=1.)
#    mu = pm.Normal('mu', mu=0., sigma=1.,
#                   shape=len(unid) * len(ulight))
#    input_current = background + mech_stim_on * mu[nid_inds +
#                                                   light_inds * len(unid)]
#    io_func = tt.nnet.sigmoid(input_current) * lam
#    mu_t = (io_func + 1) * n * 1e-3 * dt
#    obs_psth = pm.Poisson('psth', mu=mu_t,
#                          observed=psth)
##    tuning_trace = sampler.tune(tune=1000)
#    const_link_trace = trace = pm.sample(draws=2000, tune=1000)
#    pm.traceplot(trace)
#    crude_pair_plot(trace, model)
##    arviz.plot_pair(trace, divergences=True)
#    const_link_ppc = ppc = pm.sample_posterior_predictive(trace)
#    plt.figure()
#    plot_ppc(t, psth, ppc[model.name_for('psth')], 500)

# with pm.Model('link_t') as model:
#    link_t_model = model
#    lam = pm.Gamma('lam', mu=30, sigma=5)
#    background = pm.Normal('background', mu=0, sigma=1.)
#    base = pm.Normal('base', mu=0., sigma=1.,
#                     shape=len(unid) * len(ulight))
#    tau = pm.Gamma('tau', mu=50, sigma=10,
#                   shape=len(unid) * len(ulight))
#    peak = pm.Normal('mu', mu=0., sigma=1.,
#                     shape=len(unid) * len(ulight))
#    postinhpeak = pm.Normal('postinhpeak', mu=-0.2, sigma=1.,
#                            shape=len(unid))
#    postinhbase = pm.Normal('postinhbase', mu=0., sigma=0.1,
#                            shape=len(unid))
#    postinhtau = pm.Gamma('postinhtau', mu=20, sigma=10,
#                          shape=len(unid))
#    inds = nid_inds + light_inds * len(unid)
#    input_current = (background +
#                     (base[inds] +
#                      (peak[inds] - base[inds]) *
#                      tt.exp(-t / tau[inds])) * mech_stim_on +
#                     mech_stim_offset * (postinhbase[nid_inds] +
#                      (postinhpeak[nid_inds] - postinhbase[nid_inds]) *
#                      tt.exp(-offset_t / postinhtau[nid_inds]))
#                     )
#    io_func = tt.nnet.sigmoid(input_current) * lam
#    mu_t = (io_func + 1) * n * 1e-3 * dt
#    obs_psth = pm.Poisson('psth', mu=mu_t,
#                          observed=psth)
##    tuning_trace = sampler.tune(tune=1000)
#    link_t_trace = trace = pm.sample(draws=2000, tune=1000)
##    link_t_trace = trace = pm.sample(draws=10, tune=10)
#    pm.traceplot(trace)
##    crude_pair_plot(trace, model)
##    arviz.plot_pair(trace, divergences=True)
#    link_t_ppc = ppc = pm.sample_posterior_predictive(trace)
##    link_t_ppc = ppc = pm.sample_posterior_predictive(trace, samples=10)
#    fig, axs = plot_ppc(t=t,
#                        n=n,
#                        psth=psth,
#                        ppc=ppc[model.name_for('psth')],
#                        stacked_df=stacked_df,
#                        ppc_step=100,
#                        groupby=['T2', 'light'],
#                        labelby=['light'],
#                        axesby=['T2'],
#                        plot_ppc_stats=True)
#    fig.tight_layout()

# with pm.Model() as model:
#    link_t_model = model
#    lam = pm.Gamma('lam', mu=30, sigma=5)
#    background = pm.Normal('I0', mu=0, sigma=1., shape=len(unid))
#    base = pm.Normal('I1', mu=0., sigma=1.,
#                     shape=len(unid) * len(ulight))
#    tau = pm.Gamma('tau', mu=50, sigma=10,
#                   shape=len(unid) * len(ulight))
#    peak = pm.Normal('I1p', mu=0., sigma=1.,
#                     shape=len(unid) * len(ulight))
#    postinhpeak = pm.Normal('I2p', mu=0., sigma=1.,
#                            shape=len(unid))
##    postinhtau = 20.
#    postinhtau = pm.Gamma('postinhtau', mu=20, sigma=10,
#                          shape=len(unid))
#    inds = nid_inds + light_inds * len(unid)
#    input_current = (background[nid_inds]
#                     + (base[inds] +
#                        (peak[inds] - base[inds]) *
#                        tt.exp(-t / tau[inds])) * mech_stim_on
#                     + mech_stim_offset *
#                     (-postinhpeak[nid_inds]) *
#                     alpha_function(offset_t, postinhtau[nid_inds])
#                     )
#    io_func = tt.nnet.sigmoid(input_current) * lam
#    mu_t = (io_func + 1) * n * 1e-3 * dt
#    obs_psth = pm.Poisson('psth', mu=mu_t,
#                          observed=psth)
#    link_t_trace = trace = pm.sample(draws=2000, tune=1000)
##    link_t_trace = trace = pm.sample(draws=10, tune=10)
#    pm.traceplot(trace)
#    plt.tight_layout()
#    plt.savefig('parameter traces.pdf', bbox_inches='tight')
##    crude_pair_plot(trace, model)
##    arviz.plot_pair(trace, divergences=True)
#    link_t_ppc = ppc = pm.sample_posterior_predictive(trace)
##    link_t_ppc = ppc = pm.sample_posterior_predictive(trace, samples=10)
#    fig, axs = plot_ppc(t=t,
#                        n=n,
#                        psth=psth,
#                        ppc=ppc[model.name_for('psth')],
#                        stacked_df=stacked_df,
#                        ppc_step=100,
#                        groupby=['T2', 'light'],
#                        labelby=['light'],
#                        axesby=['T2'],
#                        plot_ppc_stats=True)
#    fig.tight_layout()
#    plt.savefig('Model fit.pdf', bbox_inches='tight')
# plt.show()

#with pm.Model() as model:
#    opto_I_model = model
#    lam = pm.Gamma("lam", mu=30, sigma=5)
#    I0 = pm.Normal("I0", mu=0, sigma=1.0, shape=len(unid))
#    I1 = pm.Normal("I1", mu=0.0, sigma=1.0, shape=len(unid))
#    tau = pm.Gamma("tau", mu=50, sigma=10, shape=len(unid))
#    I1p = pm.Normal("I1p", mu=0.0, sigma=1.0, shape=len(unid))
#    I2 = pm.Normal("I2", mu=0.0, sigma=1.0, shape=len(unid))
#    I2p = pm.Normal("I2p", mu=1.0, sigma=1.0, shape=len(unid))
#    I3p = pm.Normal("I3p", mu=-1.0, sigma=1.0, shape=len(unid))
#    tau_opto = pm.Gamma("tau_opto", mu=20, sigma=10, shape=len(unid))
#    tau_inh = pm.Gamma("tau_inh", mu=20, sigma=10, shape=len(unid))
#    light_lag = pm.Normal("light_lag", mu=-5, sd=1, shape=len(unid))
#    opto_stim_on = tt.switch(
#        tt.and_(
#            tt.and_(
#                tt.ge(t, light_lag[nid_inds]), tt.lt(offset_t, light_lag[nid_inds])
#            ),
#            light == 1,
#        ),
#        1.0,
#        0.0,
#    )
#    opto_stim_offset = tt.switch(
#        tt.and_(tt.ge(offset_t, light_lag[nid_inds]), light == 1), 1.0, 0.0
#    )
#    inds = nid_inds
#    input_current = (
#        I0[nid_inds]
#        + mech_stim_on * (I1[inds] + (I1p - I1)[inds] * alpha_function(t, tau[inds]))
#        + light
#        * opto_stim_on
#        * (
#            I2[inds]
#            + (I2p - I2)[inds]
#            * beta_function(t - light_lag[nid_inds], 10, tau_opto[inds])
#        )
#        + opto_stim_offset
#        * (I3p[inds] * alpha_function(offset_t - light_lag[nid_inds], tau_inh[inds]))
#    )
#    io_func = tt.nnet.sigmoid(input_current) * lam
#    mu_t = (io_func + 1) * n * 1e-3 * dt
#    obs_psth = pm.Poisson("psth", mu=mu_t, observed=psth)
#    opto_I_trace = trace = pm.sample(draws=2000, tune=1000)
#    #    opto_I_trace = trace = pm.sample(draws=10, tune=10)
#    pm.traceplot(trace)
#    plt.tight_layout()
#    plt.savefig("parameter traces.pdf", bbox_inches="tight")
#    #    crude_pair_plot(trace, model)
#    #    arviz.plot_pair(trace, divergences=True)
#    opto_I_ppc = ppc = pm.sample_posterior_predictive(trace)
#    #    opto_I_ppc = ppc = pm.sample_posterior_predictive(trace, samples=10)
#    fig, axs = plot_ppc(
#        t=t,
#        n=n,
#        psth=psth,
#        ppc=ppc[model.name_for("psth")],
#        stacked_df=stacked_df,
#        ppc_step=100,
#        groupby=["T2", "light"],
#        labelby=["light"],
#        axesby=["T2"],
#        plot_ppc_stats=True,
#    )
#    fig.tight_layout()
#    plt.savefig("Model fit.pdf", bbox_inches="tight")
#plt.show()


#with pm.Model() as model:
#    pop_model = model
#    npop = 3
#    lam = pm.Gamma("lam", mu=30, sigma=5)
#    w = pm.Dirichlet('w', a=np.ones(npop))
#    I1_pop = pm.Normal('I1_pop', mu=0, sd=1, shape=npop,
#                       transform=tr.ordered,
#                       testval=np.linspace(-1, 1, npop))
#    I1 = pm.NormalMixture('I1',
#                          w=w, mu=I1_pop, sigma=1., shape=len(unid))
#    I1p_pop = pm.Normal('I1p_pop', mu=0, sd=1, shape=npop,
#                        transform=tr.ordered,
#                        testval=np.linspace(-1, 1, npop))
#    I1p = pm.NormalMixture('I1p',
#                           w=w, mu=I1p_pop, sigma=1., shape=len(unid))
#    I2_pop = pm.Normal('I2_pop', mu=0, sd=1, shape=npop,
#                       transform=tr.ordered,
#                       testval=np.linspace(-1, 1, npop))
#    I2 = pm.NormalMixture('I2',
#                          w=w, mu=I2_pop, sigma=1., shape=len(unid))
#    I2p_pop = pm.Normal('I2p_pop', mu=0, sd=1, shape=npop,
#                        transform=tr.ordered,
#                        testval=np.linspace(-1, 1, npop))
#    I2p = pm.NormalMixture('I2p',
#                           w=w, mu=I2p_pop, sigma=1., shape=len(unid))
#    I0 = pm.Normal("I0", mu=0, sigma=1.0, shape=len(unid))
#    tau = pm.Gamma("tau", mu=50, sigma=10, shape=len(unid))
#    I3p = pm.Normal("I3p", mu=-1.0, sigma=1.0, shape=len(unid))
#    tau_opto = pm.Gamma("tau_opto", mu=20, sigma=10, shape=len(unid))
#    tau_inh = pm.Gamma("tau_inh", mu=20, sigma=10, shape=len(unid))
#    light_lag = pm.Normal("light_lag", mu=-5, sd=1, shape=len(unid))
#    opto_stim_on = tt.switch(
#        tt.and_(
#            tt.and_(
#                tt.ge(t, light_lag[nid_inds]), tt.lt(offset_t, light_lag[nid_inds])
#            ),
#            light == 1,
#        ),
#        1.0,
#        0.0,
#    )
#    opto_stim_offset = tt.switch(
#        tt.and_(tt.ge(offset_t, light_lag[nid_inds]), light == 1), 1.0, 0.0
#    )
#    inds = nid_inds
#    input_current = (
#        I0[nid_inds]
#        + mech_stim_on * (I1[inds] + (I1p - I1)[inds] * alpha_function(t, tau[inds]))
#        + light
#        * opto_stim_on
#        * (
#            I2[inds]
#            + (I2p - I2)[inds]
#            * beta_function(t - light_lag[nid_inds], 10, tau_opto[inds])
#        )
#        + opto_stim_offset
#        * (I3p[inds] * alpha_function(offset_t - light_lag[nid_inds], tau_inh[inds]))
#    )
#    io_func = tt.nnet.sigmoid(input_current) * lam
#    mu_t = (io_func + 1) * n * 1e-3 * dt
#    obs_psth = pm.Poisson("psth", mu=mu_t, observed=psth)
#    pop_trace = trace = pm.sample(draws=2000, tune=1000)
#    pm.traceplot(trace)
#    plt.tight_layout()
#    plt.savefig("parameter traces.pdf", bbox_inches="tight")
#    crude_pair_plot(trace, model, varnames=['I1_pop',
#                                            'I2_pop',
#                                            'I2p_pop'])
#    pop_ppc = ppc = pm.sample_posterior_predictive(trace)
#    fig, axs = plot_ppc(
#        t=t,
#        n=n,
#        psth=psth,
#        ppc=ppc[model.name_for("psth")],
#        stacked_df=stacked_df,
#        ppc_step=100,
#        groupby=["T2", "light"],
#        labelby=["light"],
#        axesby=["T2"],
#        plot_ppc_stats=True,
#    )
#    fig.tight_layout()
#    plt.savefig("Model fit.pdf", bbox_inches="tight")
#plt.show()


#with pm.Model() as model:
#    unmixed_pop_model = model
#    npop = 2
#    lam_pop = pm.Normal('lam_pop', mu=3.4, sigma=1.)
#    lam = pm.Gamma("lam", mu=tt.exp(lam_pop), sigma=5, shape=len(unid))
#
#    w1_pop = pm.Gamma('w1_pop', alpha=np.ones(npop), beta=np.ones(npop),
#                      shape=npop)
#    I1_pop = pm.Normal('I1_pop', mu=0, sd=0.5, shape=npop,
#                       transform=tr.ordered,
#                       testval=np.linspace(0, 1, npop))
#    w1 = pm.Dirichlet('w1', a=w1_pop, shape=(len(unid), npop))
#    I1 = pm.NormalMixture('I1', w=w1, mu=I1_pop, sd=1.,
#                          shape=len(unid), comp_shape=npop)
#
#    I1p_pop = pm.Normal('I1p_pop', mu=0, sd=0.5, shape=npop,
#                        transform=tr.ordered,
#                        testval=np.linspace(0, 1, npop))
#    I1p = pm.NormalMixture('I1p', w=w1, mu=I1p_pop, sd=1.,
#                           shape=len(unid), comp_shape=npop)
#    I1p = pm.Normal('I1p', mu=0., sd=1.,
#                    shape=len(unid))
#
#    w2 = pm.Dirichlet('w2', a=np.ones(npop))
#    I2_pop = pm.Normal('I2_pop', mu=0, sd=0.5, shape=npop,
#                       transform=tr.ordered,
#                       testval=np.linspace(0, 1, npop))
#    I2 = pm.NormalMixture('I2', w=w1, mu=I2_pop, sd=1.,
#                          shape=len(unid), comp_shape=npop)
#    I2 = pm.Normal('I2', mu=0., sd=1.,
#                   shape=len(unid))

#    I2p_pop = pm.Normal('I2p_pop', mu=0, sd=0.5, shape=npop,
#                        transform=tr.ordered,
#                        testval=np.linspace(0, 1, npop))
#    I2p = pm.NormalMixture('I2p', w=w1, mu=I2p_pop, sd=1.,
#                           shape=len(unid), comp_shape=npop)
#    I2p = pm.Normal('I2p', mu=0., sd=1.,
#                    shape=len(unid))
#
#    I3p_pop = pm.Normal('I3p_pop', mu=0, sd=1.)
#    I3p_pop_tau = pm.Gamma('I3p_pop_tau', alpha=1., beta=2.)
#    I3p_offset = pm.Normal('I3p_offset', mu=0, sd=1, shape=len(unid))
#    I3p = pm.Deterministic('I3p', I3p_pop + I3p_offset / tt.sqrt(I3p_pop_tau))
#
#    I0 = pm.Normal("I0", mu=0, sigma=1.0, shape=len(unid))
#
#    tau = pm.Gamma("tau", mu=50, sigma=1)
#    tau_opto = pm.Gamma("tau_opto", mu=20, sigma=1)
#    tau_inh = pm.Gamma("tau_inh", mu=20, sigma=1)
#    light_lag = np.zeros(len(unid))
#    opto_stim_on = tt.switch(
#        tt.and_(
#            tt.and_(
#                tt.ge(t, light_lag[nid_inds]), tt.lt(offset_t, light_lag[nid_inds])
#            ),
#            light == 1,
#        ),
#        1.0,
#        0.0,
#    )
#    opto_stim_offset = tt.switch(
#        tt.and_(tt.ge(offset_t, light_lag[nid_inds]), light == 1), 1.0, 0.0
#    )
#    inds = nid_inds
#    input_current = (
#        I0[nid_inds]
#        + mech_stim_on * (I1[inds] + (I1p - I1)[inds] * alpha_function(t, tau))
#        + light
#        * opto_stim_on
#        * (
#            I2[inds]
#            + (I2p - I2)[inds]
#            * beta_function(t - light_lag[nid_inds], 10, tau_opto)
#        )
#        + opto_stim_offset
#        * (I3p[inds] * alpha_function(offset_t - light_lag[nid_inds], tau_inh))
#    )
#    io_func = tt.nnet.sigmoid(input_current) * lam[nid_inds]
#    mu_t = (io_func + 1) * n * 1e-3 * dt
#    obs_psth = pm.Poisson("psth", mu=mu_t, observed=psth)
#    unmixed_pop_trace = trace = pm.sample(draws=2000, tune=1000)
#    dump(trace, "trace.gz")
#
#    pm.traceplot(trace)
#    plt.tight_layout()
#    plt.savefig("parameter traces.pdf", bbox_inches="tight")
#    crude_pair_plot(trace, model, varnames=['w1',
#                                            'w2',
#                                            'lam',
#                                            'I1_pop',
#                                            'I2_pop',
#                                            'I2p_pop',
#                                            'I3p_pop'])
#    unmixed_pop_ppc = ppc = pm.sample_posterior_predictive(trace)
#    fig, axs = plot_ppc(
#        t=t,
#        n=n,
#        psth=psth,
#        ppc=ppc[model.name_for("psth")],
#        stacked_df=stacked_df,
#        ppc_step=100,
#        groupby=["T2", "light"],
#        labelby=["light"],
#        axesby=["T2"],
#        plot_ppc_stats=True,
#    )
#    fig.tight_layout()
#    plt.savefig("Model fit.pdf", bbox_inches="tight")
#plt.show()


def mixture_builder(name_root,
                    n=100,
                    npop=2,
                    alpha=None,
                    p_n=None,
                    pop_mu=None,
                    alpha_kwargs=None,
                    mu_kwargs=None,
                    set_default_mu_order=True,
                    mixture_kwargs=None,
                    ):
    if alpha_kwargs is None:
        alpha_kwargs = {}
    alpha_kwargs.setdefault('alpha', 1.)
    alpha_kwargs.setdefault('beta', 1.)

    if mu_kwargs is None:
        mu_kwargs = {}
    mu_kwargs.setdefault('mu', 0.)
    mu_kwargs.setdefault('sigma', 0.5)
    if set_default_mu_order:
        if npop > 1:
            mu_kwargs.setdefault('transform', tr.ordered)
            mu_kwargs.setdefault('testval', np.linspace(0, 1, npop))

    if mixture_kwargs is None:
        mixture_kwargs = {}
    mixture_kwargs.setdefault('sigma', 1.)

    if alpha is None:
        if npop > 1:
            alpha = pm.Gamma('{}_pop_w'.format(name_root), shape=npop, **alpha_kwargs)
    if p_n is None:
        if npop > 1:
            p_n = pm.Dirichlet('{}_n_prob'.format(name_root), alpha, shape=(n, npop))
    if pop_mu is None:
        pop_mu = pm.Normal('{}_pop_mu'.format(name_root), shape=npop, **mu_kwargs)
    if npop > 1:
        mu = pm.NormalMixture('{}_mu'.format(name_root), w=p_n, mu=pop_mu, shape=n,
                               comp_shape=(1, npop), **mixture_kwargs)
    else:
        mu = pop_mu
    offset = pm.Normal('{}_offset'.format(name_root), mu=0, sigma=1, shape=n)
    x_i = pm.Deterministic(name_root, mu + offset)
    return alpha, p_n, pop_mu, mu, offset, x_i


with pm.Model() as model:
    mixed_pop_model = model
    n_mech_pop = 1
    n_opto_pop = 2
    lam_pop = pm.Normal('lam_pop', mu=3.4, sigma=1.)
    lam = pm.Gamma("lam", mu=tt.exp(lam_pop), sigma=5, shape=len(unid))

    alpha1, p_n1, _, _, _, I1 = mixture_builder('I1', len(unid), n_mech_pop)
    I1p = mixture_builder('I1p', len(unid), n_mech_pop, alpha=alpha1, p_n=p_n1)[-1]
    alpha2, p_n2, _, _, _, I2 = mixture_builder('I2', len(unid), n_opto_pop)
    I2p = mixture_builder('I2p', len(unid), n_opto_pop, alpha=alpha2, p_n=p_n2)[-1]
    I3p = mixture_builder('I3p', len(unid), n_opto_pop, alpha=alpha2, p_n=p_n2)[-1]
    I0 = mixture_builder("I0", len(unid), 1)[-1]

    tau = pm.Gamma("tau", mu=50, sigma=1)
    tau_opto = pm.Gamma("tau_opto", mu=20, sigma=1)
    tau_inh = pm.Gamma("tau_inh", mu=20, sigma=1)
    light_lag = pm.Normal("light_lag", 0., 1., shape=len(unid))
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
            + (I2p - I2)[inds]
            * beta_function(t - light_lag[nid_inds], 10, tau_opto)
        )
        + opto_stim_offset
        * (I3p[inds] * alpha_function(offset_t - light_lag[nid_inds], tau_inh))
    )
    io_func = tt.nnet.sigmoid(input_current) * lam[nid_inds]
    mu_t = (io_func + 1e-4) * n# * 1e-3 * dt
    obs_psth = pm.Poisson("psth", mu=mu_t, observed=psth)
    fitter = pm.FullRankADVI()
    approx = fitter.fit(n=int(2e6))
    trace = approx.sample(draws=2000)
    mixed_pop_trace = trace = pm.sample(draws=2000, tune=1000,
                                        init='advi', target_accept=0.9)
    dump(trace, "trace.gz")

    pm.traceplot(trace)
    plt.tight_layout()
    plt.savefig("parameter traces.pdf", bbox_inches="tight")
    crude_pair_plot(trace, model, varnames=['I1_pop_w',
                                            'I2_pop_w',
                                            'lam',])
#                                            'I1_mu',
#                                            'I2_mu',
#                                            'I2p_mu',
#                                            'I3p_mu'])
    plt.savefig("parameter pairplot.pdf", bbox_inches="tight")
    mixed_pop_ppc = ppc = pm.sample_posterior_predictive(trace)
    fig, axs = plot_ppc(
        t=t,
        n=n,
        psth=psth,
        ppc=ppc[model.name_for("psth")],
        stacked_df=stacked_df,
        ppc_step=100,
        groupby=["T2", "light"],
        labelby=["light"],
        axesby=["T2"],
        plot_ppc_stats=True,
    )
    fig.tight_layout()
    plt.savefig("Model fit.pdf", bbox_inches="tight")
plt.show()


#with pm.Model() as model:
#    pop_model = model
#    npop = 2
#    lam_pop = pm.Normal('lam_pop', mu=3.4, sigma=1.)
#    lam = pm.Gamma("lam", mu=tt.exp(lam_pop), sigma=5, shape=len(unid))
##    lam_ = pm.Gamma("lam", mu=30, sigma=5)
##    lam = tt.stack([lam_] * len(unid))
#    # Mechanical responsive mixture
#    w1 = pm.Dirichlet('w1', a=np.ones(npop))
#    I1_pop = pm.Normal('I1_pop', mu=0, sd=0.5, shape=npop,
#                       transform=tr.ordered,
#                       testval=np.linspace(0, 1, npop))
#    I1 = pm.NormalMixture('I1', w=w1, mu=I1_pop, sd=1., shape=len(unid))
#    I1p_pop = pm.Normal('I1p_pop', mu=0, sd=0.5, shape=npop,
#                        transform=tr.ordered,
#                        testval=np.linspace(0, 1, npop))
#    I1p = pm.NormalMixture('I1p', w=w1, mu=I1p_pop, sd=1., shape=len(unid))
#
#    # Opto responsive mixture
##    w2 = pm.Dirichlet('w2', a=np.ones(npop))
##    I2_pop = pm.Normal('I2_pop', mu=0, sd=1, shape=npop,
##                       transform=tr.ordered,
##                       testval=np.linspace(-1, 1, npop))
##    I2 = pm.NormalMixture('I2',
##                          w=w2, mu=I2_pop, sigma=1., shape=len(unid))
##    # I2p's ordering should be established by I2's ordering through w2
##    I2p_pop = pm.Normal('I2p_pop', mu=0, sd=1, shape=npop,
##                        transform=tr.ordered,
##                        testval=np.linspace(-1, 1, npop))
##    I2p = pm.NormalMixture('I2p',
##                           w=w2, mu=I2p_pop, sigma=1., shape=len(unid))
##    I3p_pop = pm.Normal('I3p_pop', mu=0, sd=1, shape=npop)
##    I3p = pm.NormalMixture('I3p',
##                           w=w2, mu=I3p_pop, sigma=1., shape=len(unid))
#    I2 = pm.Normal('I2', mu=0, sigma=1., shape=len(unid))
#    I2p = pm.Normal('I2p', mu=0., sigma=1., shape=len(unid))
#    I3p = pm.Normal('I3p', mu=0., sigma=1., shape=len(unid))
#    I0 = pm.Normal("I0", mu=0, sigma=1.0, shape=len(unid))
#    tau = pm.Gamma("tau", mu=50, sigma=10)
#    tau_opto = pm.Gamma("tau_opto", mu=20, sigma=10)
#    tau_inh = pm.Gamma("tau_inh", mu=20, sigma=10)
##    light_lag = pm.Normal("light_lag", mu=-5, sd=1, shape=len(unid))
#    light_lag = np.zeros(len(unid))
#    opto_stim_on = tt.switch(
#        tt.and_(
#            tt.and_(
#                tt.ge(t, light_lag[nid_inds]), tt.lt(offset_t, light_lag[nid_inds])
#            ),
#            light == 1,
#        ),
#        1.0,
#        0.0,
#    )
#    opto_stim_offset = tt.switch(
#        tt.and_(tt.ge(offset_t, light_lag[nid_inds]), light == 1), 1.0, 0.0
#    )
#    inds = nid_inds
#    input_current = (
#        I0[nid_inds]
#        + mech_stim_on * (I1[inds] + (I1p - I1)[inds] * alpha_function(t, tau))
#        + light
#        * opto_stim_on
#        * (
#            I2[inds]
#            + (I2p - I2)[inds]
#            * beta_function(t - light_lag[nid_inds], 10, tau_opto)
#        )
#        + opto_stim_offset
#        * (I3p[inds] * alpha_function(offset_t - light_lag[nid_inds], tau_inh))
#    )
#    io_func = tt.nnet.sigmoid(input_current) * lam[nid_inds]
#    mu_t = (io_func + 1) * n * 1e-3 * dt
#    obs_psth = pm.Poisson("psth", mu=mu_t, observed=psth)
#    pop_trace = trace = pm.sample(draws=2000, tune=1000)
#    pm.traceplot(trace)
#    plt.tight_layout()
#    plt.savefig("parameter traces.pdf", bbox_inches="tight")
#    crude_pair_plot(trace, model, varnames=['w1',
#                                            'w2',
#                                            'lam',
#                                            'I1_pop',
#                                            'I2_pop',
#                                            'I2p_pop',
#                                            'I3p_pop'])
#    pop_ppc = ppc = pm.sample_posterior_predictive(trace)
#    fig, axs = plot_ppc(
#        t=t,
#        n=n,
#        psth=psth,
#        ppc=ppc[model.name_for("psth")],
#        stacked_df=stacked_df,
#        ppc_step=100,
#        groupby=["T2", "light"],
#        labelby=["light"],
#        axesby=["T2"],
#        plot_ppc_stats=True,
#    )
#    fig.tight_layout()
#    plt.savefig("Model fit.pdf", bbox_inches="tight")
#plt.show()
#
#dump(trace, "trace.gz")
