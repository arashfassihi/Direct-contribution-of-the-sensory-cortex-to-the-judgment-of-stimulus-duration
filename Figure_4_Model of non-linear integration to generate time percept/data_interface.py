#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 14:07:39 2019

@author: lucianopaz
"""

import os
from collections import defaultdict
import numpy as np
from scipy.io import loadmat
import pandas as pd
from compress_pickle import load


__all__ = [
    "load_raw_data",
    "set_time_vector",
    "mask_psth",
    "get_stacked_stim_course",
    "load_fittable_data",
    "merge_psths",
]


DT = 10.0
PRE = 400.0
LIGHT_ON_SHIFT = 20.0
DELAY = 500.0
ACTION = 500.0


def load_raw_data(path=".", version=""):
    raw_psth = loadmat(os.path.join(path, "PSTH{}.mat".format(version)))["PSTH"]
    raw_resp = loadmat(os.path.join(path, "resp{}.mat".format(version)))["resp"]

    light_map = {0: "Optooff", 1: "Optoon"}
    vibration_responsive = np.squeeze(raw_resp["Vibration"][0][0])
    opto_responsive = np.squeeze(raw_resp["Otpo"][0][0])
    is_single_unit = np.squeeze(raw_resp["IsSingleunit"][0][0])
    n_suffix = "_NumberofTrials"
    T1 = 334.0
    T2s = np.array([161.0, 205.0, 264.0, 334.0, 423.0, 545.0, 694.0])
    data = []
    for light, label in light_map.items():
        frate = raw_psth[label][0][0]
        ntrials = raw_psth[label + n_suffix][0][0]
        for c, T2 in enumerate(T2s):
            temp1 = frate[0][c]
            temp2 = ntrials[0][c]
            for nid, vals in enumerate(
                zip(temp1, temp2, vibration_responsive, opto_responsive, is_single_unit)
            ):
                rate, n, vresp, oreps, sing = vals
                n = n[0]
                if n > 0:
                    psth = np.round(rate * n * 1e-3)
                    if version == "_v2":
                        psth = psth[60:]
                    data.append([nid, vresp, oreps, sing, T1, T2, light, n, psth])
    columns = [
        "nid",
        "mech_resp",
        "opto_resp",
        "single_unit",
        "T1",
        "T2",
        "light",
        "ntrials",
        "psth",
    ]
    return pd.DataFrame(data=data, columns=columns)


def load_behavioral_data(path=".", version=""):
    raw = loadmat(os.path.join(path, "behavior{}.mat".format(version)))
    contra = raw["contra"]
    ipsi = raw["ipsi"]
    index = pd.Index(contra[:, 0].astype(np.float), name="T2")
    columns = pd.MultiIndex.from_product(
        [
            ["contra", "ipsi"],
            ["ntrials_light1", "ntrials_light2", "nchoice2_light1", "nchoice2_light2"],
        ],
        names=["hemisphere", "data_header"],
    )
    return pd.DataFrame(
        data=np.concatenate([contra[:, 1:], ipsi[:, 1:]], axis=1).astype(np.int),
        index=index,
        columns=columns,
    )


def set_time_vector(series, dt=DT, pre=PRE, light_on_shift=LIGHT_ON_SHIFT):
    if series.light == 1:
        shift = light_on_shift
    else:
        shift = 0.0
    t = np.arange(len(series.psth), dtype=np.float) * dt - pre + shift
    d = series.to_dict()
    d["t"] = t
    return pd.Series(d)


def mask_psth(series, unmasked_onset=-np.inf, unmasked_offset=np.inf):
    mask = np.logical_or(
        series["t"] > series["T2"] + unmasked_offset, series["t"] < unmasked_onset
    )
    masked_psth = np.ma.array(series["psth"], mask=mask)
    d = series.to_dict()
    d["psth"] = masked_psth
    return pd.Series(d)


def get_stacked_stim_course(df, dt=DT):
    column_names = df.columns
    values = defaultdict(list)
    for ind, row in df.iterrows():
        psth = row["psth"]
        if isinstance(psth, np.ma.masked_array):
            inds = np.logical_not(psth.mask)
        else:
            inds = None
        row_values = []
        for name in column_names:
            v = row[name]
            if isinstance(v, (np.ndarray, np.ma.masked_array)):
                if inds is not None:
                    v = v[inds]
            v = np.asarray(v)
            row_values.append(v)
        broadcasted_row = np.broadcast_arrays(*row_values)
        T2_ind = list(column_names).index("T2")
        for name, value in zip(column_names, broadcasted_row):
            values[name].append(value)
            if name == "t":
                offset_t = value - np.ceil(broadcasted_row[T2_ind] / dt) * dt
                values["offset_t"].append(offset_t)
    return pd.DataFrame.from_dict({k: np.concatenate(v) for k, v in values.items()})


def load_fittable_data(
    path=".",
    dt=DT,
    pre=PRE,
    light_on_shift=LIGHT_ON_SHIFT,
    unmasked_onset=-300.0,
    unmasked_offset=270.0,
    query=None,
    filter_non_responsive=False,
    version="_v2",
):
    df = (
        load_raw_data(version=version)
        .apply(set_time_vector, axis=1, dt=dt, pre=pre, light_on_shift=light_on_shift)
        .apply(
            mask_psth,
            axis=1,
            unmasked_onset=unmasked_onset,
            unmasked_offset=unmasked_offset,
        )
    )
    if query is not None:
        df = df.query(query)
    if filter_non_responsive:
        responsive = load("responsive_neurons.pkl")
        valid = responsive[responsive.values].index.values
        df = df.query("nid in @valid")
    return get_stacked_stim_course(df, dt=dt)


def merge_psths(df, groupby=None):
    if groupby is None:
        groupby = [c for c in df.columns if c not in ["nid", "ntrials", "psth"]]
    groupby = list(groupby)
    out = defaultdict(lambda x: [[], []])
    for g, gdf in df.groupby(groupby):
        psth = gdf["psth"].values
        if len(psth) > 1:
            psth = sum(psth)
        else:
            psth = psth[0]
        ntrials = gdf["ntrials"].values
        if len(ntrials) > 1:
            ntrials = sum(ntrials)
        else:
            ntrials = ntrials[0]
        out[g] = [ntrials, psth]
    out = pd.DataFrame.from_dict(data=out, orient="index", columns=["ntrials", "psth"])
    return out.set_index(pd.MultiIndex.from_tuples(out.index.values, names=groupby))
