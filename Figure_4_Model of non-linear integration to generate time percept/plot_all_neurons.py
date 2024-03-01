#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 15:00:19 2019

@author: lucianopaz
"""


import numpy as np
import matplotlib as mt
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from data_interface import load_fittable_data, DT
from plot_utils import plot_ppc
from matplotlib.colors import LinearSegmentedColormap
from tqdm import tqdm
#mt.use('pdf')


a = 0.2
lim = (3 / 8 - a) / (3 / 4 - a)

cdict1 = {'red': ((0.0, 0., 0.),
                  (1.0, 0., 0.)),
          'green': ((0.0, 0.8, 0.8),
                    (1.0, 0.36, 0.36)),
          'blue': ((0.0, 0.0, 0.0),
                   (1.0, 1.0, 1.0)),
          }
bg_cmap = LinearSegmentedColormap('BG', cdict1)


dt = DT
stacked_df = load_fittable_data(
    light_on_shift=0,
)
#area_1 = [46, 139]
#area_2 = [139, 228]
#inds1 = np.logical_and(stacked_df['nid'].values >= area_1[0],
#                       stacked_df['nid'].values < area_1[1])
#inds2 = np.logical_and(stacked_df['nid'].values >= area_2[0],
#                       stacked_df['nid'].values < area_2[1])
#stacked_df['t'][inds1] -= 20
#stacked_df['offset_t'][inds1] -= 20
#stacked_df['t'][inds2] -= 10
#stacked_df['offset_t'][inds2] -= 10

unid = sorted(stacked_df.nid.unique())

Ts = [stacked_df.T2.min(), 545.]
stacked_df = stacked_df.query('T2 in @Ts')

t = stacked_df["t"].values
psth = stacked_df["psth"].values
n = stacked_df["ntrials"].values

fig, axs = plot_ppc(
    t=t,
    n=n,
    psth=psth,
    ppc=None,
    stacked_df=stacked_df,
    ppc_step=100,
    groupby=["T2", "light"],
    labelby=["light"],
    axesby=["T2"],
    plot_ppc_stats=True,
    marker='',
    cmap=bg_cmap,
    linewidth=2,
)
fig.suptitle('All neurons')
plt.savefig('All neurons.svg', bbox_inches='tight')
#plt.close(fig)


#with PdfPages('All neurons.pdf') as f:
#    t = stacked_df["t"].values
#    psth = stacked_df["psth"].values
#    n = stacked_df["ntrials"].values
#
#    fig, axs = plot_ppc(
#        t=t,
#        n=n,
#        psth=psth,
#        ppc=None,
#        stacked_df=stacked_df,
#        ppc_step=100,
#        groupby=["T2", "light"],
#        labelby=["light"],
#        axesby=["T2"],
#        plot_ppc_stats=True,
#        marker='o',
#    )
#    fig.suptitle('All neurons')
#    f.savefig(fig)
#    plt.close(fig)
#
#    for nid in tqdm(unid):
#        df = stacked_df.query('nid == {}'.format(nid))
#        t = df["t"].values
#        psth = df["psth"].values
#        n = df["ntrials"].values
#
#        fig, axs = plot_ppc(
#            t=t,
#            n=n,
#            psth=psth,
#            ppc=None,
#            stacked_df=df,
#            ppc_step=100,
#            groupby=["T2", "light"],
#            labelby=["light"],
#            axesby=["T2"],
#            plot_ppc_stats=True,
#            marker='o',
#        )
#        fig.suptitle('Neuron {}'.format(nid))
#        f.savefig(fig)
#        plt.close(fig)

