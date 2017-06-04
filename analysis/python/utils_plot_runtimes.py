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


def plot_all_runtimes(datadir):
    # -> serial, without I/O
    datafile = 'srl_N_x_100_1_1_1.dat'
    datapath = os.path.join(datadir, datafile)
    (Nx, _, _, _, _, _, time_srl, _) \
    = u.get_data(datapath, sort_with='Nx')
    Nx = Nx[:-1] # remove size 5793, but not others
    time_srl = time_srl[:-1]

    # -> serial, with I/O
    datafile = 'srl_N_x_100_1_1_1_IO.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_srl_IO, _) \
    = u.get_data(datapath, sort_with='Nx')

    # -> best MPI, without I/O, on 16 processors (15 workers)
    datafile = 'mpi_N_x_100_x_16_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_mpi16, _) \
    = u.get_data(datapath, sort_with='Nx')

    # -> best MPI, with I/O, on 16 processors (15 workers)
    datafile = 'mpi_N_x_100_x_16_1_IO.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_mpi16_IO, _) \
    = u.get_data(datapath, sort_with='Nx')

    # -> best MPI, without I/O, on 64 processors (63 workers)
    datafile = 'mpi_N_x_100_x_64_1.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_mpi64, _) \
    = u.get_data(datapath, sort_with='Nx')

    # -> best MPI, with I/O, on 64 processors (63 workers)
    datafile = 'mpi_N_x_100_x_64_1_IO.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_mpi64_IO, _) \
    = u.get_data(datapath, sort_with='Nx')

    # -> hybrid, without I/O, on 16 processors (15 workers) and 4 threads/proc
    datafile = 'hyb_N_x_100_x_16_4.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_hyb, _) \
    = u.get_data(datapath, sort_with='Nx')

    # -> hybrid, with I/O, on 16 processors (15 workers) and 4 threads/proc
    datafile = 'hyb_N_x_100_x_16_4_IO.dat'
    datapath = os.path.join(datadir, datafile)
    (_, _, _, _, _, _, time_hyb_IO, _) \
    = u.get_data(datapath, sort_with='Nx')


    fig = plt.figure(figsize=(12, 4))
    gs = GridSpec(1, 2, width_ratios=[1, 0.3], height_ratios=[1,])
    gs.update(wspace=0.2, hspace=0.05)

    ax1 = plt.subplot(gs[0])
    ax1.loglog(Nx, time_srl, marker='P', ms=ms, lw=lw, ls=':', 
               basex=2, zorder=10, color='#1F77B4',
               label=r"$\rm{serial}$")
    ax1.loglog(Nx, time_srl_IO, marker='P', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, color='#1F77B4', dashes=(5, 1.5),
               label=r"$\rm{serial,\ I/O}$")
    ax1.loglog(Nx, time_mpi16, marker='o', ms=ms, lw=lw, ls=':', 
               basex=2, zorder=10, color='#FF7F0E',
               label=r"$\rm{MPI},\ $"+r"$p=16$")
    ax1.loglog(Nx, time_mpi16_IO, marker='o', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, color='#FF7F0E', dashes=(5, 1.5),
               label=r"$\rm{MPI},\ $"+r"$p=16,\ \rm{I/O}$")
    ax1.loglog(Nx, time_mpi64, marker='s', ms=ms, lw=lw, ls=':', 
               basex=2, zorder=10, color='#2CA02C',
               label=r"$\rm{MPI},\ $"+r"$p=64$")
    ax1.loglog(Nx, time_mpi64_IO, marker='s', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, color='#2CA02C', dashes=(5, 1.5),
               label=r"$\rm{MPI},\ $"+r"$p=64,\ \rm{I/O}$")
    ax1.loglog(Nx, time_hyb, marker='d', ms=ms, lw=lw, ls=':', 
               basex=2, zorder=10, color='#D62728',
               label=r"$\rm{hybrid},\ $"+r"$p=16,\ t=4$")
    ax1.loglog(Nx, time_hyb_IO, marker='d', ms=ms, lw=lw, ls='--', 
               basex=2, zorder=10, color='#D62728', dashes=(5, 1.5),
               label=r"$\rm{hybrid},\ $"+r"$p=16,\ t=4,\ \rm{I/O}$")
    ax1.xaxis.set_ticks([1024, 1448, 2048, 2896, 4096])
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.set_xlabel(r"$\rm{image\ size,}\ $"+r"$N\ \rm{[px]}$", fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{run\ time},\ $"+r"$T\ \rm{[s]}$", 
                   fontsize=prp.labelsize_tex)

    ax1.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
               fontsize=prp.labelsize_tex)

    return fig
