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
    # -> data with division on N rows (best case), 4 threads / proc
    datafile = 'hyb_weak_x_100_x_x_4.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc_4, _, time_4, _) \
    = u.get_data(datapath, sort_with='n_proc')

    # -> data with division on N rows (best case), 8 threads / proc
    datafile = 'hyb_weak_x_100_x_x_8.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc_8, _, time_8, _) \
    = u.get_data(datapath, sort_with='n_proc')

    # -> equivalent data from serial code
    datafile = 'srl_N_x_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath, sort_with='Nx')
    time_srl_red_for4 = time_srl[2:] # because 4thr/rank ran only from size 2048 to 5793
    time_srl_red_for8 = time_srl[3:] # because 8thr/rank ran only from size 2896 to 5793

    # speedups
    S_4  = time_srl_red_for4 / time_4
    S_8  = time_srl_red_for8 / time_8
    # efficiencies
    E_4 = S_4 / n_proc_4
    E_8 = S_8 / n_proc_8
    E_4_corr = S_4 / (n_proc_4-1)
    E_8_corr = S_8 / (n_proc_8-1)

    # Comparison datas :
    # -> best mpi data
    datafile = 'mpi_weak_x_100_x_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc_mpi, _, time_mpi, _) \
    = u.get_data(datapath, sort_with='n_proc')
    S_mpi = time_srl / time_mpi
    E_mpi = S_mpi / n_proc_mpi
    E_mpi_corr = S_mpi / (n_proc_mpi-1)


    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(2, 2, width_ratios=[0.5, 1], height_ratios=[1, 1])
    gs.update(wspace=0.4, hspace=0.15)

    ax1 = plt.subplot(gs[1, 0])
    ax1.loglog(n_proc_4, time_4, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, label=r"$t=4$")
    ax1.loglog(n_proc_8, time_8, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, label=r"$t=8$")
    ax1.loglog(n_proc_mpi, time_mpi, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, label=r"$\rm{pure\ MPI\ }$"+r"$(t=1)$")
    ax1.set_xlim(1.5, 20)
    ax1.xaxis.set_ticks([2, 4, 8, 16])
    # ax1.yaxis.set_ticks([7e-2, 8e-2, 9e-2, 1e-1, 1.2e-1], minor=False)
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_base10))
    ax1.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{run\ time},\ $" + r"$T\ \rm{[s]}$", 
                   fontsize=prp.labelsize_tex)

    ax2 = plt.subplot(gs[0, 1])
    ax2.loglog(n_proc_4, S_4, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, label=r"$t=4$")
    ax2.loglog(n_proc_8, S_8, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, label=r"$t=8$")
    ax2.loglog(n_proc_mpi, S_mpi, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, 
               label=r"$\rm{pure\ MPI\ }$"+r"$(t=1)$")
    ax2.loglog(n_proc_mpi, u.identity(n_proc_mpi), lw=2, ls=':', color='black',
               basex=2, basey=2, zorder=2)
    ax2.text(0.08, 0.56, r"$\rm{ideal}$",
             verticalalignment='top', horizontalalignment='left', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.set_ylim(-10, 30)
    ax2.yaxis.set_ticks([1, 4, 16])
    ax2.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.set_ylabel(r"$\rm{speedup},\ $"+r"$S_p$", fontsize=prp.labelsize_tex)
    plt.setp(ax2.get_xticklabels(), visible=False)

    ax3 = plt.subplot(gs[1, 1], sharex=ax2)
    ax3.semilogx(n_proc_4, E_4, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=7, label=r"$t=4$")
    ax3.semilogx(n_proc_4, E_4_corr, marker='o', ms=ms, lw=lw, ls=':', 
               basex=2, zorder=3, color='#666666')
    ax3.semilogx(n_proc_8, E_8, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=7, label=r"$t=8$")
    ax3.semilogx(n_proc_8, E_8_corr, marker='s', ms=ms, lw=lw, ls=':', 
               basex=2, zorder=5, color='#AAAAAA')
    ax3.semilogx(n_proc_mpi, E_mpi, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=7, label=r"$\rm{pure\ MPI\ }$"+r"$(t=1)$")
    ax3.semilogx(n_proc_mpi, E_mpi_corr, marker='d', ms=ms, lw=lw, ls=':', 
               basex=2, zorder=3, 
               label=r"$\rm{(all\ dotted)\ same,\ norm.\ to\ }$"+r"$p_{\rm{w}}$",
               color='#AAAAAA')
    ax3.semilogx(n_proc_mpi, u.identity(n_proc_mpi)/n_proc_mpi, lw=2, ls=':', 
               basex=2, color='black', zorder=2)
    ax3.text(0.8, 0.95, r"$\rm{ideal}$",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax3.transAxes)
    ax3.yaxis.set_ticks([0.25, 0.5, 0.75, 1.0])
    ax3.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax3.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_float_p2))
    ax3.set_xlim(1.5, 20)
    ax3.set_ylim(0.25, 1.2)
    ax3.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", fontsize=prp.labelsize_tex)
    ax3.set_ylabel(r"$\rm{efficiency},\ $"+r"$E_p$", fontsize=prp.labelsize_tex)

    ax3.legend(loc='lower center', bbox_to_anchor=(-0.63, 1.1),
               fontsize=prp.labelsize_tex)

    return fig


def plot_strong_scaling(datadir):
    # -> data with division on 5000 rows (best case), 4 threads / proc
    datafile = 'hyb_strong_5000_100_5000_x_4.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc_4, _, time_4, time_err) \
    = u.get_data(datapath, sort_with='n_proc')

    # -> data with division on 5000 rows (best case), 8 threads / proc
    datafile = 'hyb_strong_5000_100_5000_x_8.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc_8, _, time_8, time_err) \
    = u.get_data(datapath, sort_with='n_proc')

    # Comparison data:
    # -> pure-MPI, with division on 5000 rows
    datafile = 'mpi_strong_5000_100_5000_x_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, n_proc_mpi, _, time_mpi, _) \
    = u.get_data(datapath, sort_with='n_proc')

    datafile = 'srl_N_5000_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath)
    time_srl = time_srl[0] # time_srl contains only one value in this case

    # speedups
    S_4 = time_srl  / time_4
    S_8 = time_srl  / time_8
    S_mpi = time_srl / time_mpi

    fig = plt.figure(figsize=(12, 3))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.2], height_ratios=[1,])
    gs.update(wspace=0.45, hspace=0.05)

    ax1 = plt.subplot(gs[0])
    ax1.loglog(n_proc_4, time_4, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=8, 
               label=r"$t=4$")
    ax1.loglog(n_proc_8, time_8, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=7, 
               label=r"$t=8$")
    ax1.loglog(n_proc_mpi, time_mpi, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=10, 
               label=r"$\rm{pure\ MPI\ }$"+r"$(t=1)$")
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_base10))
    ax1.set_ylim(0.1, 5)
    ax1.text(0.95, 0.93, r"$n_{\rm{row}}=N=5\,000$",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax1.transAxes)
    ax1.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{run\ time},\ $"+r"$T\ $"+r"$\rm{[s]}$", 
                   fontsize=prp.labelsize_tex)

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.loglog(n_proc_4, S_4, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=8, 
               label=r"$t=4$")
    ax2.loglog(n_proc_8, S_8, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=7, 
               label=r"$t=8$")
    ax2.loglog(n_proc_mpi, S_mpi, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, 
               label=r"$\rm{pure\ MPI\ }$"+r"$(t=1)$")
    ax2.loglog(n_proc_4, u.identity(n_proc_4), lw=2, ls=':', color='black',
               basex=2, basey=2, zorder=2)
    ax2.text(0.05, 0.93, r"$n_{\rm{row}}=N=5\,000$",
             verticalalignment='top', horizontalalignment='left', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.text(0.25, 0.45, r"$\rm{ideal}$",
             verticalalignment='bottom', horizontalalignment='right', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.set_xlim(1.5, 20)
    ax2.set_ylim(-1, 28)
    ax2.xaxis.set_ticks([2, 4, 8, 16])
    ax2.yaxis.set_ticks([1, 2, 4, 8, 16])
    ax2.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", fontsize=prp.labelsize_tex)
    ax2.set_ylabel(r"$\rm{speedup},\ $"+r"$S_p$", fontsize=prp.labelsize_tex)

    ax2.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
               fontsize=prp.labelsize_tex)

    return fig


def plot_strong_scaling_threads(datadir):
    datafile = 'hyb_strong_1000_100_1_2_x.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, n_thread, time_1000, time_err) \
    = u.get_data(datapath, sort_with='n_thread')

    datafile = 'hyb_strong_5000_100_1_2_x.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, n_thread, time_5000, time_err) \
    = u.get_data(datapath, sort_with='n_thread')

    datafile = 'hyb_strong_10000_100_1_2_x.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_10000, time_err) \
    = u.get_data(datapath, sort_with='n_thread')

    datafile = 'srl_N_1000_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath)
    time_srl_1000 = time_srl[0] # time_srl contains only one value in this case

    datafile = 'srl_N_5000_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath)
    time_srl_5000 = time_srl[0] # time_srl contains only one value in this case

    datafile = 'srl_N_10000_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath)
    time_srl_10000 = time_srl[0] # time_srl contains only one value in this case

    # speedups
    S_1000  = time_srl_1000  / time_1000
    S_5000  = time_srl_5000  / time_5000
    S_10000 = time_srl_10000 / time_10000

    fig = plt.figure(figsize=(12, 3))
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.2], height_ratios=[1,])
    gs.update(wspace=0.45, hspace=0.05)

    ax1 = plt.subplot(gs[0])
    ax1.loglog(n_thread, time_1000, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=10, 
               label=r"$N=1\,000$")
    ax1.loglog(n_thread, time_5000, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=10, 
               label=r"$N=5\,000$")
    ax1.loglog(n_thread, time_10000, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=10, zorder=10, 
               label=r"$N=10\,000$")
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_base10))
    ax1.text(0.95, 0.95, r"$p=2\ (p_{\rm{w}}=1)$",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax1.transAxes)
    ax1.set_xlabel(r"$\rm{number\ of\ threads},\ $"+r"$t$", fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{run\ time},\ $" + r"$T\ \rm{[s]}$", 
                   fontsize=prp.labelsize_tex)

    ax2 = plt.subplot(gs[1], sharex=ax1)
    ax2.loglog(n_thread, S_1000, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, 
               label=r"$N=1\,000$")
    ax2.loglog(n_thread, S_5000, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, 
               label=r"$N=5\,000$")
    ax2.loglog(n_thread, S_10000, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, basey=2, zorder=10, 
               label=r"$N=10\,000$")
    ax2.loglog(n_thread, u.amdahl_law(n_thread, 0.05), lw=2, ls=':', 
               color='#888888', basex=2, basey=2, zorder=2)
    ax2.annotate(r"$\rm{Amdahl},\ \alpha=0.05$", xy=(13, 8), xytext=(2.4, 1.2),
                 arrowprops=dict(color='#888888', 
                                 arrowstyle='->'),
                 color='#888888', fontsize=prp.labelsize_tex)
    ax2.loglog(n_thread, u.identity(n_thread), lw=2, ls=':', color='black',
               basex=2, basey=2, zorder=2)
    ax2.text(0.05, 0.95, r"$p=2\ (p_{\rm{w}}=1)$",
             verticalalignment='top', horizontalalignment='left', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.text(0.65, 0.92, r"$\rm{ideal}$",
             verticalalignment='top', horizontalalignment='left', 
             color='black', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.set_xlim(0.8, 20)
    ax2.set_ylim(0.8, 20)
    ax2.xaxis.set_ticks([1, 2, 4, 8, 16])
    ax2.yaxis.set_ticks([1, 2, 4, 8, 16])
    ax2.set_xlabel(r"$\rm{number\ of\ threads},\ $"+r"$t$", fontsize=prp.labelsize_tex)
    ax2.set_ylabel(r"$\rm{speedup},\ $"+r"$S_t$", fontsize=prp.labelsize_tex)

    ax2.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
               fontsize=prp.labelsize_tex)

    return fig


