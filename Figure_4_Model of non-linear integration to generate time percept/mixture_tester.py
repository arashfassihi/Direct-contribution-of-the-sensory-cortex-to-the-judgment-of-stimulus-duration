#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 12:01:53 2019

@author: lucianopaz
"""

import numpy as np
import pymc3 as pm
from pymc3.distributions import transforms as tr
from plot_utils import crude_pair_plot
from matplotlib import pyplot as plt
import seaborn as sns


def NormalMixture(name, w, mu, sd, *args, **kwargs):
    try:
        npop = mu.distribution.shape[0]
    except AttributeError:
        try:
            npop = mu.shape[0].eval()
        except Exception:
            npop = len(mu)
    comp_dists = [pm.Normal.dist(mu[i], sd=sd) for i in range(npop)]
    return pm.Mixture(name, w=w,
                      comp_dists=comp_dists, *args, **kwargs)


def factory(observed=None, n=5, npop=2, point=None):
    if point is None:
        point = {}
    with pm.Model() as model:
        # Mechanical responsive mixture
        if 'w1' in point:
            w1 = point.get('w1')
        else:
            w1 = pm.Dirichlet('w1', a=np.ones(npop))
        if 'I1_pop' in point:
            I1_pop= point.get('I1_pop')
        else:
            I1_pop= pm.Normal('I1_pop', mu=0, sd=0.2, shape=npop,
                              transform=tr.ordered,
                              testval=np.linspace(0, 1, npop))
        I1_pop = np.arange(npop)
#        I1_offset = pm.Normal('I1_offset', mu=0, sd=1, shape=n)
        I1_offset = np.arange(n)
        if observed is None:
            
            I1 = pm.NormalMixture('I1', w=w1, mu=I1_offset[..., None] + I1_pop,
                                  sd=0.2,
                                  comp_shape=(n, npop),
                                  shape=n)
        else:
            I1 = pm.NormalMixture('I1', w=w1, mu=I1_offset[..., None] + I1_pop,
                                  sd=0.2,
                                  comp_shape=(n, npop),
                                  shape=n,
                                  observed=observed)
    return model


p = 0.4
point = {'w1': np.array([p, 1 - p]),
         'I1_pop_mu': np.array([-0.5, 0.9])}

with factory(point=point):
    prior = pm.sample_prior_predictive(2000)
    crude_pair_plot(prior, varnames=['w1',
                                     'I1_pop',
                                     'I1'])


with factory(observed=prior['I1']) as model:
    trace = pm.sample()
    pm.traceplot(trace)
    crude_pair_plot(trace, varnames=['w1',
                                     'I1_pop',
                                     'I1'])
    ppc = pm.sample_posterior_predictive(trace, 100)
    fig, axs = plt.subplots(prior['I1'].shape[-1], 1,
                            sharex=True,
                            figsize=(6, 10))
    for ax, pri, pos in zip(axs, prior['I1'].T, np.squeeze(ppc['I1']).T):
        for pos_ in pos.T:
            ax.hist(pos_, 100, color='b', alpha=0.05,
                    histtype='step', density=True)
        ax.hist(pri, 100, color='r', linewidth=2,
                histtype='step', density=True)
    fig.tight_layout()


#import numpy as np
#from theano import tensor as tt
#import pymc3 as pm
#from matplotlib import pyplot as plt
#
#
#npop = 5
#nd = 3
#with pm.Model():
#    mus = np.tile(np.arange(npop)[None, :], (nd, 1))
#    m1 = pm.NormalMixture('m1',
#                          w=np.ones(npop) / npop,
#                          mu=mus,
#                          sigma=1e-5,
#                          comp_shape=(nd, npop),
#                          shape=nd)
#    m2 = pm.NormalMixture('m2',
#                          w=np.ones((nd, npop)) / npop,
#                          mu=mus,
#                          sigma=1e-5,
#                          comp_shape=(nd, npop),
#                          shape=nd)
#
#m1_val = m1.random(size=100)
#m2_val = m2.random(size=100)
#
#print('m1 rows come from the same compenent? {}'.
#      format(all(np.all(np.diff(m1_val) < 1e-3, axis=-1))))
#print('m2 rows come from the same compenent? {}'.
#      format(all(np.all(np.diff(m2_val) < 1e-3, axis=-1))))
#
#plt.figure()
#plt.subplot(121)
#plt.plot(np.all(np.diff(m1_val) < 1e-3, axis=-1))
#plt.subplot(122)
#plt.plot(np.all(np.diff(m2_val) < 1e-3, axis=-1))
#
#
#npop = 5
#nd = 3
#with pm.Model():
#    mus = tt.as_tensor_variable(np.tile(np.arange(npop)[None, :], (nd, 1)))
#    z1 = pm.Categorical('z1', p=np.ones(npop) / npop)
#    m1 = pm.Normal('m1', mu=mus[..., z1], sigma=1e-5, shape=nd)
#    z2 = pm.Categorical('z2', p=np.ones(npop) / npop, shape=nd)
#    mu2 = tt.as_tensor_variable([mus[i, z2[i]] for i in range(nd)])
#    m2 = pm.Normal('m2', mu=mu2, sigma=1e-5, shape=nd)
#
#m1_val = m1.random(size=100)
#m2_val = m2.random(size=100)
#
#print('m1 rows come from the same compenent? {}'.
#      format(all(np.all(np.diff(m1_val) < 1e-3, axis=-1))))
#print('m2 rows come from the same compenent? {}'.
#      format(all(np.all(np.diff(m2_val) < 1e-3, axis=-1))))
#
#plt.figure()
#plt.subplot(121)
#plt.plot(np.all(np.diff(m1_val) < 1e-3, axis=-1))
#plt.subplot(122)
#plt.plot(np.all(np.diff(m2_val) < 1e-3, axis=-1))