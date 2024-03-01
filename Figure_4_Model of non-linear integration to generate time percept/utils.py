#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 14:10:32 2019

@author: lucianopaz
"""

import numpy as np
from theano import tensor as tt


__all__ = ["to_tuple", "alpha_function", "beta_function"]


def to_tuple(x):
    try:
        return tuple(x)
    except TypeError:
        return (x,)


def _np_alpha(t, tau, peak=1.0, tmin=None, tmax=None):
    alpha_ = peak * t / tau * np.exp(1 - t / tau)
    if tmin is None and tmax is None:
        return alpha_
    elif tmin is None:
        return np.where(t < tmax, alpha_, 0.0)
    elif tmax is None:
        return np.where(t > tmin, alpha_, 0.0)
    else:
        return np.where(np.logical_and(t > tmin, t < tmax), alpha_, 0.0)


def _np__beta(t, tau_rise, tau_fall, peak=1.0, tmin=None, tmax=None):
    tpeak = tau_fall * tau_rise / (tau_fall - tau_rise) * np.log(tau_fall / tau_rise)
    norm = np.exp(-tpeak / tau_fall) - np.exp(-tpeak / tau_rise)
    beta_ = peak * (np.exp(-t / tau_fall) - np.exp(-t / tau_rise)) / norm
    if tmin is None and tmax is None:
        return beta_
    elif tmin is None:
        return np.where(t < tmax, beta_, 0.0)
    elif tmax is None:
        return np.where(t > tmin, beta_, 0.0)
    else:
        return np.where(np.logical_and(t > tmin, t < tmax), beta_, 0.0)


def _np_beta(t, tau_rise, tau_fall, peak=1.0, tol=1e-3, tmin=None, tmax=None):
    return np.where(
        np.abs(tau_rise - tau_fall) < tol,
        _np_alpha(t, tau_rise, peak, tmin, tmax),
        np.where(
            tau_fall > tau_rise,
            _np__beta(t, tau_rise, tau_fall, peak, tmin, tmax),
            _np__beta(t, tau_fall, tau_rise, peak, tmin, tmax),
        ),
    )


def alpha_function(t, tau, peak=1.0, tmin=None, tmax=None):
    alpha_ = peak * t / tau * tt.exp(1 - t / tau)
    if tmin is None and tmax is None:
        return alpha_
    elif tmin is None:
        return tt.switch(tt.lt(t, tmax), alpha_, 0.0)
    elif tmax is None:
        return tt.switch(tt.gt(t, tmin), alpha_, 0.0)
    else:
        return tt.switch(tt.and_(tt.gt(t, tmin), tt.lt(t, tmax)), alpha_, 0.0)


def _beta_function(t, tau_rise, tau_fall, peak=1.0, tmin=None, tmax=None):
    tpeak = tau_fall * tau_rise / (tau_fall - tau_rise) * tt.log(tau_fall / tau_rise)
    norm = tt.exp(-tpeak / tau_fall) - tt.exp(-tpeak / tau_rise)
    beta_ = peak * (tt.exp(-t / tau_fall) - tt.exp(-t / tau_rise)) / norm
    if tmin is None and tmax is None:
        return beta_
    elif tmin is None:
        return tt.switch(tt.lt(t, tmax), beta_, 0.0)
    elif tmax is None:
        return tt.switch(tt.gt(t, tmin), beta_, 0.0)
    else:
        return tt.switch(tt.and_(tt.gt(t, tmin), tt.lt(t, tmax)), beta_, 0.0)


def beta_function(t, tau_rise, tau_fall, peak=1.0, tol=1e-3, tmin=None, tmax=None):
    return tt.switch(
        tt.lt(abs(tau_rise - tau_fall), tol),
        alpha_function(t, tau_rise, peak, tmin, tmax),
        tt.switch(
            tau_fall > tau_rise,
            _beta_function(t, tau_rise, tau_fall, peak, tmin, tmax),
            _beta_function(t, tau_fall, tau_rise, peak, tmin, tmax),
        ),
    )
