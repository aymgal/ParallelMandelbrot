import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
import os
from matplotlib.gridspec import GridSpec

import paramplots as prp
import utils as u

# some global plot parameters
ms = 9 # marker size
lw = 1.8 # line width

def plot_timetosol(datadir):
    datafile = 'srl_N_x_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (Nx, _, _, _, _, _, time_100, _) \
    = u.get_data(datapath, sort_with='Nx')

    datafile = 'srl_N_x_200_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (Nx, _, _, _, _, _, time_200, _) \
    = u.get_data(datapath, sort_with='Nx')

    datafile = 'srl_N_x_300_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (Nx, _, _, _, _, _, time_300, _) \
    = u.get_data(datapath, sort_with='Nx')

    fig = plt.figure(figsize=(6, 5))
    gs = GridSpec(1, 1, width_ratios=[1,], height_ratios=[1,])
    gs.update(wspace=0.05, hspace=0.05)
    ax1 = plt.subplot(gs[0])
    ax1.loglog(Nx, time_100, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, label=r"$n_{\rm{max}}=100$")
    ax1.loglog(Nx, time_200, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, label=r"$n_{\rm{max}}=200$")
    ax1.loglog(Nx, time_300, marker='v', ms=ms, lw=lw, ls='--', 
               basex=2, label=r"$n_{\rm{max}}=300$")
    ax1.set_ylim(0.8e-1, 1e1)
    ax1.xaxis.set_ticks([1024, 1448, 2048, 2896, 4096, 5793])
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.set_xlabel(r"$N$")
    ax1.set_ylabel(r"$T\ \rm{[s]}$")
    ax1.legend(loc='upper left')
    return fig