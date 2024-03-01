#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 09:51:31 2019

@author: lucianopaz
"""
import numpy as np
from matplotlib import pyplot as plt
from data_interface import load_fittable_data, DT, merge_psths
from utils import _np_alpha as alpha_function
from utils import _np_beta as beta_function
from scipy.special import expit


dt = DT

stacked_df = load_fittable_data(light_on_shift=dt)
stacked_df = merge_psths(stacked_df).reset_index()
stacked_df = stacked_df.assign(nid=-1)
t = stacked_df["t"].values + dt
offset_t = stacked_df["offset_t"].values + dt
psth = stacked_df["psth"].values
n = stacked_df["ntrials"].values
nid = stacked_df["nid"].values
light = stacked_df["light"].values

mech_stim_on = np.logical_and(t >= 0, offset_t < 0)
opto_stim_on = np.logical_and(np.logical_and(t >= -10, offset_t < -10), light)
opto_stim_offset = np.logical_and(offset_t >= -10, light)

unid, nid_inds = np.unique(nid, return_inverse=True)
ulight, light_inds = np.unique(light, return_inverse=True)

inds = nid_inds
I0 = np.zeros(len(unid))
I1 = np.ones(len(unid))
I1p = 2 * np.ones(len(unid))
I2 = np.ones(len(unid))
I2p = 2 * np.ones(len(unid))
I3 = np.ones(len(unid))
I3p = 2 * np.ones(len(unid))
lam = 40.0
tau = 50.0 * np.ones(len(unid))
tau_opto = 50.0 * np.ones(len(unid))
tau_inh = 50.0 * np.ones(len(unid))
light_lag = -10. * np.ones(len(unid))

mech_current = mech_stim_on * (
    I1[inds] + (I1p - I1)[inds] * alpha_function(t, tau[inds])
)
opto_current = (
    light
    * opto_stim_on
    * (
        I2[inds]
        + (I2p - I2)[inds] * beta_function(t - light_lag[nid_inds], 10, tau_opto[inds])
    )
)
inh_current = opto_stim_offset * (
    -I3p[inds] * alpha_function(offset_t - light_lag[nid_inds], tau_inh[inds])
)

input_current = I0[nid_inds] + mech_current + opto_current + inh_current
io_func = lam * expit(input_current)
mu_t = (io_func + 1) * n * 1e-3 * dt

for T2 in stacked_df.T2.unique():
    fig, axs = plt.subplots(4, 1, figsize=(10, 10), sharex=True)
    _inds = stacked_df.T2 == T2
    axs[0].plot(t[_inds], mech_current[_inds], '.')
    axs[1].plot(t[_inds], opto_current[_inds], '.')
    axs[2].plot(t[_inds], inh_current[_inds], '.')
    axs[3].plot(t[_inds], mu_t[_inds], '.')
