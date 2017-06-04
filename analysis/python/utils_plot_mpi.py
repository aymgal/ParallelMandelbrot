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

def plot_weak_scaling(datadir):
    # -> data with division on 100 rows
    datafile = 'mpi_weak_x_100_100_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc, _, time_100, _) \
    = u.get_data(datapath, sort_with='n_proc')

    # -> data with division on 250 rows
    datafile = 'mpi_weak_x_100_250_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (Nx, _, _, _, n_proc, _, time_250, _) \
    = u.get_data(datapath, sort_with='n_proc')

    # -> data with division on 500 rows
    datafile = 'mpi_weak_x_100_500_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_500, _) \
    = u.get_data(datapath, sort_with='n_proc')

    # -> data with division in max number of rows (= N)
    datafile = 'mpi_weak_x_100_x_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_varrow, _) \
    = u.get_data(datapath, sort_with='n_proc')

    # -> equivalent data from serial code
    datafile = 'srl_N_x_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath, sort_with='Nx')

    # speedups
    S_100    = time_srl / time_100
    S_250    = time_srl / time_250
    S_500    = time_srl / time_500
    S_varrow = time_srl / time_varrow

    # efficiencies
    E_100    = S_100 / n_proc
    E_250    = S_250 / n_proc
    E_500    = S_500 / n_proc
    E_varrow = S_varrow / n_proc
    E_varrow_corr = S_varrow / (n_proc-1)

    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(2, 2, width_ratios=[0.5, 1], height_ratios=[1, 1])
    gs.update(wspace=0.4, hspace=0.15)

    ax1 = plt.subplot(gs[1, 0])
    ax1.loglog(n_proc, time_100, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, label=r"$n_{\rm{row}}=100$")
    ax1.loglog(n_proc, time_250, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, label=r"$n_{\rm{row}}=250$")
    ax1.loglog(n_proc, time_varrow, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, color='#D62728',
               label=r"$n_{\rm{row}}=N$")
    ax1.set_xlim(1.5, 80)
    ax1.xaxis.set_ticks([2, 4, 8, 16, 32, 64])
    # ax1.yaxis.set_ticks([7e-2, 8e-2, 9e-2, 1e-1, 1.2e-1], minor=False)
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    # ax1.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_base10))
    ax1.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", 
                   fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{run\ time},\ $" + r"$T\ \rm{[s]}$", 
                   fontsize=prp.labelsize_tex)

    ax2 = plt.subplot(gs[0, 1])
    ax2.loglog(n_proc, S_100, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, label=r"$n_{\rm{row}}=100$")
    ax2.loglog(n_proc, S_250, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, label=r"$n_{\rm{row}}=250$")
    ax2.loglog(n_proc, S_varrow, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, color='#D62728',
               label=r"$n_{\rm{row}}=N$")
    ax2.loglog(n_proc, u.identity(n_proc), lw=2, ls=':', color='black',
               basex=2, basey=2, zorder=2)
    ax2.text(0.08, 0.4, r"$\rm{ideal}$",
             verticalalignment='top', horizontalalignment='left', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.set_ylim(-10, 90)
    ax2.yaxis.set_ticks([1, 4, 16, 64])
    ax2.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.set_ylabel(r"$\rm{speedup},\ $"+r"$S_p$", fontsize=prp.labelsize_tex)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax3 = plt.subplot(gs[1, 1], sharex=ax2)
    ax3.semilogx(n_proc, E_100, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=7, label=r"$n_{\rm{row}}=100$")
    ax3.semilogx(n_proc, E_250, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=8, label=r"$n_{\rm{row}}=250$")
    ax3.semilogx(n_proc, E_varrow, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, color='#D62728',
               label=r"$n_{\rm{row}}=N$")
    ax3.semilogx(n_proc, E_varrow_corr, marker='x', ms=ms, lw=lw, ls='--', 
               color='#666666', basex=2, zorder=4, 
               label=r"$n_{\rm{row}}=N,\ \rm{norm.\ to}\ $"+r"$p=p_{\rm{w}}$")
    ax3.semilogx(n_proc, u.gustafson_law(n_proc, 0.05)/n_proc, lw=2, ls=':', 
               color='#888888', basex=2, zorder=100)
    ax3.annotate(r"$\rm{Gustafson},\ \beta=0.05$", xy=(40, 0.94), xytext=(10, 1.38),
                 arrowprops=dict(color='#888888', 
                                 arrowstyle='->'),
                 color='#888888', fontsize=prp.labelsize_tex)
    ax3.semilogx(n_proc, u.identity(n_proc)/n_proc, lw=2, ls=':', 
                 basex=2, color='black', zorder=2)
    ax3.text(0.08, 0.94, r"$\rm{ideal}$",
             verticalalignment='top', horizontalalignment='left', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax3.transAxes)
    ax3.yaxis.set_ticks([0.25, 0.5, 0.75, 1.0])
    ax3.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax3.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_float_p2))
    ax3.set_xlim(1.5, 80)
    ax3.set_ylim(0.4, 1.16)
    ax3.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", fontsize=prp.labelsize_tex)
    ax3.set_ylabel(r"$\rm{efficiency},\ $"+r"$E_p$", fontsize=prp.labelsize_tex)

    ax3.legend(loc='lower center', bbox_to_anchor=(-0.66, 1.1),
               fontsize=prp.labelsize_tex)

    return fig


def plot_strong_scaling(datadir):
    datafile = 'mpi_strong_10000_100_x_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc, _, time_varrows, time_err) \
    = u.get_data(datapath, sort_with='n_proc')

    datafile = 'mpi_strong_10000_100_100_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc, _, time_100, time_err) \
    = u.get_data(datapath, sort_with='n_proc')

    datafile = 'mpi_strong_10000_100_250_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc, _, time_250, time_err) \
    = u.get_data(datapath, sort_with='n_proc')

    datafile = 'mpi_strong_10000_100_1000_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_10000, _) \
    = u.get_data(datapath, sort_with='n_proc')

    datafile = 'srl_N_10000_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath)
    time_srl = time_srl[0] # time_srl contains only one value in this case

    # speedups
    S_varrows = time_srl / time_varrows
    S_100 = time_srl / time_100
    S_250 = time_srl / time_250
    S_10000 = time_srl / time_10000

    fig = plt.figure(figsize=(12, 3))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.2], height_ratios=[1,])
    gs.update(wspace=0.45, hspace=0.05)

    ax1 = plt.subplot(gs[0])
    ax1.loglog(n_proc, time_100, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=7, 
               label=r"$n_{\rm{row}}=100$")
    ax1.loglog(n_proc, time_250, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=8, 
               label=r"$n_{\rm{row}}=250$")
    ax1.loglog(n_proc, time_10000, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=10, color='#D62728',
               label=r"$n_{\rm{row}}=N$")
    ax1.loglog(n_proc, time_varrows, marker='v', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=5, color='#2CA02C',
               label=r"$n_{\rm{row}}=p-1$")
    ax1.yaxis.set_ticks([0.3, 1, 3, 10])
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_base10))
    ax1.text(0.95, 0.95, r"$N=10\,000$",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax1.transAxes)
    ax1.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", 
                   fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{run\ time},\ $" + r"$T\ \rm{[s]}$", 
                   fontsize=prp.labelsize_tex)

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.loglog(n_proc, S_100, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=7, 
               label=r"$n_{\rm{row}}=100$")
    ax2.loglog(n_proc, S_250, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=8, 
               label=r"$n_{\rm{row}}=250$")
    ax2.loglog(n_proc, S_10000, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, color='#D62728',
               label=r"$n_{\rm{row}}=N$")
    ax2.loglog(n_proc, S_varrows, marker='v', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=5, color='#2CA02C',
               label=r"$n_{\rm{row}}=p_{\rm{w}}$")
    ax2.loglog(n_proc, u.amdahl_law(n_proc, 0.05), lw=2, ls=':', 
               color='#888888', basex=2, basey=2, zorder=2)
    ax2.annotate(r"$\rm{Amdahl},\ \alpha=0.05$", xy=(20, 10), xytext=(6, 1.2),
                 arrowprops=dict(color='#888888', 
                                 arrowstyle='->'),
                 color='#888888', fontsize=prp.labelsize_tex)
    ax2.loglog(n_proc, u.identity(n_proc), lw=2, ls=':', color='black',
               basex=2, basey=2, zorder=2)
    ax2.text(0.8, 0.92, r"$\rm{ideal}$",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.xaxis.set_ticks([2, 4, 8, 16, 32, 64])
    ax2.yaxis.set_ticks([1, 2, 4, 8, 16, 32, 64])
    ax2.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.set_xlim(1.6, 76)
    ax2.set_ylim(0.8, 76)
    ax2.text(0.05, 0.92, r"$N=10\,000$",
             verticalalignment='top', horizontalalignment='left', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", 
                   fontsize=prp.labelsize_tex)
    ax2.set_ylabel(r"$\rm{speedup},\ $"+r"$S_p$", fontsize=prp.labelsize_tex)

    ax2.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), 
               fontsize=prp.labelsize_tex)

    return fig


