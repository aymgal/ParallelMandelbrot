import argparse as argp
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import ScalarFormatter, FuncFormatter

import paramplots as prp

reportdir = '../../report'
datadir = '../cluster/outputs'
figdir  = 'figures'
save = True # save plots
save_in_report = True # if False save figures in the directory of the report

def parser_arguments():
    parser = argp.ArgumentParser(prog="Plot builder",
                description="Plot builder from mandelbrot program outputs",
                formatter_class=argp.RawTextHelpFormatter)
    help_m = """    choice of program version
    0 for serial
    1 for OpenMP
    2 for MPI
    3 for hybrid
    """
    # parser.add_argument('-v', '--version', dest='version', type=int,
    #                     metavar='\b', default=0, help=help_m)
    parser.add_argument('-t', '--timetosol', dest='timetosol',
                        help="plot of time to solution vs. number of procs", 
                        action='store_true')
    parser.add_argument('-s', '--speedup', dest='speedup',
                        help="plot of speedup vs. number of procs", 
                        action='store_true')
    parser.add_argument('-e', '--efficiency', dest='efficiency',
                        help="plot of efficiency vs. number of procs", 
                        action='store_true')
    parser.add_argument('-ws', '--weak-scaling', dest='wscaling',
                        help="plot of weak scaling vs. number of procs", 
                        action='store_true')
    parser.add_argument('-sl', '--speedup-laws', dest='laws', 
                        help="plot Amdahl's law and Gustafson's law",
                        action='store_true')
    return parser.parse_args()

def integer_ticks(value, index):
    # return r'${:d}$'.format(int(value))
    return '{:d}'.format(int(value))


def identity(p):
    return p

def amdahl_law(p, alpha):
    return 1.0 / (alpha + (1.0 - alpha)/p)

def gustafson_law(p, beta):
    return p - beta * (p - 1.0)

def average_every_n(array, n, std_too=False):
    """array should be 1D, 
    and every successive n values supposed to represent the same quantity

    stackoverflow.com/(...)/averaging-over-every-n-elements-of-a-numpy-array
    """
    if std_too:
        return [np.mean(array.reshape(-1, n), axis=1),
                np.std(array.reshape(-1, n), axis=1) ]
    else:
        return np.mean(array.reshape(-1, n), axis=1)


if __name__ == '__main__':

    pargs = parser_arguments()

    # if pargs.version == 0:
    #     datafile = 'srl_'
    # if pargs.version == 1:
    #     datafile = 'omp_'
    # elif pargs.version == 2:
    #     datafile = 'mpi_'
    # elif pargs.version == 3:
    #     datafile = 'hyb_'

    # if pargs.wscaling:
    #     datafile += 'weak_x_100_100_x_0.dat'

    # datapath = os.path.join(datadir, datafile)
    # datalen = sum(1 for line in open(datapath)) # number of lines in file

    # Nx, Ny, n_iter, n_row, \
    # n_proc, n_thread, time = np.loadtxt(datapath, delimiter=' ', unpack=True)

    # # since every run is done 10 times, we must average values
    # Nx       = average_every_n(Nx, 10)
    # Ny       = average_every_n(Ny, 10)
    # n_iter   = average_every_n(n_iter, 10)
    # n_row    = average_every_n(n_row, 10)
    # n_proc   = average_every_n(n_proc, 10)
    # n_thread = average_every_n(n_thread, 10)
    # time     = average_every_n(time, 10)

    # fig = plt.figure(figsize=(6, 5))
    # gs = GridSpec(1, 1, width_ratios=[1,], height_ratios=[1,])
    # gs.update(wspace=0.05, hspace=0.05) # set the spacing between axes
    # ax1 = plt.subplot(gs[0])
    # ax1.loglog(n_proc, time)
    # plt.show()

    # ---------
    if pargs.laws:

        figname = 'speedup_laws'

        p_rge = np.logspace(0, 5, 1e2, base=2)

        fig = plt.figure(figsize=(11, 4))
        gs = GridSpec(1, 2, width_ratios=[1, 1], height_ratios=[1,])
        gs.update(wspace=0.3, hspace=0.05)

        ax1 = plt.subplot(gs[0])
        ax1.loglog(p_rge, identity(p_rge), lw=1, ls=':', color='black')
        ax1.text(0.85, 0.94, "ideal",
                 verticalalignment='top', horizontalalignment='right', 
                 color='black', fontsize=prp.globalfsize, 
                 bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
                 transform=ax1.transAxes)
        ax1.loglog(p_rge, amdahl_law(p_rge, 0.05), ls='-', basex=2, basey=2,
                   label=r"$\alpha=0.05$")
        ax1.loglog(p_rge, amdahl_law(p_rge, 0.35), ls='--', basex=2, basey=2,
                   label=r"$\alpha=0.35$")
        # ax1.xaxis.set_major_formatter(ScalarFormatter())
        ax1.xaxis.set_major_formatter(FuncFormatter(integer_ticks))
        ax1.yaxis.set_major_formatter(FuncFormatter(integer_ticks))
        ax1.set_xlim(1, p_rge.max())
        ax1.set_ylim(0, p_rge.max())
        ax1.set_xlabel(r"$p$", fontsize=prp.labelsize_tex)
        ax1.set_ylabel(r"$S_p$", fontsize=prp.labelsize_tex)
        ax1.set_title(r"Amdahl's law")
        ax1.legend(loc='upper left', frameon=False, handletextpad=0.1)

        ax2 = plt.subplot(gs[1])
        ax2.loglog(p_rge, identity(p_rge), lw=1, ls=':', color='black')
        ax2.text(0.85, 0.94, "ideal",
                 verticalalignment='top', horizontalalignment='right', 
                 color='black', fontsize=prp.globalfsize, 
                 bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
                 transform=ax2.transAxes)
        ax2.loglog(p_rge, gustafson_law(p_rge, 0.05), ls='-', basex=2, basey=2,
                   label=r"$\beta=0.05$")
        ax2.loglog(p_rge, gustafson_law(p_rge, 0.35), ls='--', basex=2, basey=2,
                   label=r"$\beta=0.35$")
        ax2.xaxis.set_major_formatter(FuncFormatter(integer_ticks))
        ax2.yaxis.set_major_formatter(FuncFormatter(integer_ticks))
        ax2.set_xlim(1, p_rge.max())
        ax2.set_ylim(0, p_rge.max())
        ax2.set_xlabel(r"$p$", fontsize=prp.labelsize_tex)
        ax2.set_ylabel(r"$S_p$", fontsize=prp.labelsize_tex)
        ax2.set_title("Gustafson's law")
        ax2.legend(loc='upper left', frameon=False, handletextpad=0.1)

        plt.show()
        if save:
            if save_in_report:
                figdir = os.path.join(reportdir, figdir)
            prp.savefig(fig, figdir, figname)

