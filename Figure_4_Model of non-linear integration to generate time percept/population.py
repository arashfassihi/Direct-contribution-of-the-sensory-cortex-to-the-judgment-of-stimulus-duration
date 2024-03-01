#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 07:41:52 2019

@author: lucianopaz
"""

import numpy as np
import pandas as pd
from scipy.special import erf, expit
from sklearn.mixture import BayesianGaussianMixture
from numba import njit, float64, vectorize

convolve = np.vectorize(np.convolve, signature="(i),(j)->(k)")


def richards(x, nu, lam):
    return lam * (1 - 1 / (1 + np.exp(10 * x)) ** (1 / nu))


@vectorize(
    [
        float64(
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
            float64,
        )
    ],
    target="parallel",
)
def response_rate(
    t, T, mech, light, nu, lam, I0, I1, I1p, I2, I2p, I3p, tau, tau_opto, tau_inh, onset
):
    t = t - onset
    offset_t = t - T
    if mech:
        w = expit(t) * expit(-offset_t)
    if light:
        w2 = expit(offset_t)
    if mech and light:
        input_current = I0 + (
            w
            * (
                I1
                + (I1p - I1) * np.exp(-t / tau)
                + I2
                + (I2p - I2) * np.exp(-t / tau_opto)
            )
            + w2 * I3p * np.exp(-offset_t / tau_inh)
        )
    elif mech:
        input_current = I0 + (w * (I1 + (I1p - I1) * np.exp(-t / tau)))
    elif light:
        input_current = I0 + (
            w * (I2 + (I2p - I2) * np.exp(-t / tau_opto))
            + w2 * I3p * np.exp(-offset_t / tau_inh)
        )
    else:
        input_current = I0
    return lam * (1 - 1 / (1 + np.exp(10 * input_current)) ** (1 / nu))


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
            gmm_on_parameters = ["nu", "lam", "I1", "I1p", "I2", "I2p", "I3p"]
        self.gmm_parameters = self.parameters[gmm_on_parameters]
        self.gmm_columns = gmm_on_parameters
        self.median_columns = [c for c in self.parameters if c not in gmm_on_parameters]
        self.medians = self.parameters[self.median_columns].median(axis=0)

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
        rate = richards(input_current, nu, lam)
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

    def only_contralateral_light_drive(
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
            mech=False,
            light=False,
            samples=ipsilateral_parameters,
            reduce=True,
        )
        return contra_rate + ipsi_rate

    def only_ipsilateral_light_drive(
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
            light=False,
            samples=contralateral_parameters,
            reduce=True,
        )
        ipsi_rate = self.response_rate(
            t=t[..., None],
            T=T[..., None],
            mech=False,
            light=True,
            samples=ipsilateral_parameters,
            reduce=True,
        )
        return contra_rate + ipsi_rate

    def leaky_integrator(self, t, drive, tau=600.0):
        kernel = np.exp(-t / tau)
        LI = convolve(drive, kernel)[..., : kernel.shape[-1]]
        return LI

    def psychometric(
        self,
        Ts,
        n_samples=1000,
        reps=100,
        tau=600.0,
        background_mean=None,
        background_var=None,
    ):
        if background_mean is None:
            background_mean = -0.1
        if background_var is None:
            background_var = 1.0
        background_mean *= n_samples
        background_var *= n_samples
        T = np.asarray(Ts[-1])
        t = np.linspace(0, 694, 695)
        drives = []
        for rep in range(reps):
            contralateral_parameters = self.sample(n_samples=n_samples)
            ipsilateral_parameters = self.sample(n_samples=n_samples)
            no_light_drive = self.no_light_drive(
                t,
                T,
                contralateral_parameters=contralateral_parameters,
                ipsilateral_parameters=ipsilateral_parameters,
            )
            contra_light_drive = self.contralateral_light_drive(
                t,
                T,
                contralateral_parameters=contralateral_parameters,
                ipsilateral_parameters=ipsilateral_parameters,
            )
            ipsi_light_drive = self.ipsilateral_light_drive(
                t,
                T,
                contralateral_parameters=contralateral_parameters,
                ipsilateral_parameters=ipsilateral_parameters,
            )
            drives.append([no_light_drive, contra_light_drive, ipsi_light_drive])
        drives = np.moveaxis(np.array(drives), 0, -2)
        LI = []
        LI_var = []
        for drive in drives:
            li = self.leaky_integrator(t, drive, tau) + background_mean * tau * (
                1 - np.exp(-t / tau)
            )
            li_var = self.leaky_integrator(
                t, drive, 0.5 * tau
            ) + background_var * 0.5 * tau * (1 - np.exp(-2 * t / tau))
            LI.append(li)
            LI_var.append(li_var)
        LI = np.array(LI)
        LI_var = np.array(LI_var)
        T_ref_ind = np.argmin(np.abs(t - Ts[len(Ts) // 2]))
        ref = LI[..., T_ref_ind]
        ref_vars = LI_var[..., T_ref_ind]
        dprime = (LI - ref[:, None, :, None]) / np.sqrt(
            2 * (LI_var + ref_vars[:, None, :, None])
        )
        psychometric = 0.5 + 0.5 * erf(dprime)
        return drives, LI, LI_var, dprime, psychometric
