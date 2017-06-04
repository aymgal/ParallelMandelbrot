import argparse as argp
import math
import matplotlib.pyplot as plt
import numpy as np
import os

import paramplots as prp
import utils_plot_srl as psrl
import utils_plot_mpi as pmpi
import utils_plot_hyb as phyb
import utils_plot_runtimes as prun
import utils_plot_laws as plaws

reportdir = '../../report'
datadir = '../cluster/outputs'
figdir  = 'figures'
save = False # save plots
save_in_report = True # if False save figures in the directory of the report


if __name__ == '__main__':

    parser = argp.ArgumentParser(prog="Plot builder",
                description="Plot builder from mandelbrot program outputs",
                formatter_class=argp.RawTextHelpFormatter)
    help_m = """choice of program version
    0 for serial (nothing specific to this yet)
    1 for OpenMP (nothing specific to this yet)
    2 for MPI
    3 for hybrid"""
    help_runtime = """plot run time comparison of serial, MPI 
and hybrid versions, each case with and without I/O vs problem size N"""
    parser.add_argument('-c', '--code-version', dest='code_version', type=int,
                        metavar='\b', default=0, help=help_m)
    parser.add_argument('-w', '--weak-scaling', dest='wscaling',
                        help="plot weak scaling test", 
                        action='store_true')
    parser.add_argument('-s', '--strong-scaling', dest='sscaling',
                        help="plot strong scaling test", 
                        action='store_true')
    parser.add_argument('-t', '--s-scaling-thr', dest='sscaling_thr',
                        help="plot strong scaling test (nb of threads varying)", 
                        action='store_true')
    parser.add_argument('--serial-N', dest='serial_N',
                        help="plot serial run time vs N", 
                        action='store_true')
    parser.add_argument('--all-runtime', dest='runtime_comp',
                        help=help_runtime, 
                        action='store_true')
    parser.add_argument('--laws', dest='laws', 
                        help="plot Amdahl's law and Gustafson's law",
                        action='store_true')
    parser.add_argument('--save-fig', dest='save', 
                        help="save figure in the directory set in this script",
                        action='store_true')

    pargs = parser.parse_args()

# --------------------------------------------------------------------------- #
    fig = None

    # WEAK SCALING
    if ((pargs.code_version == 2) and pargs.wscaling):
        figname = 'mpi_weak_scaling'
        fig = pmpi.plot_weak_scaling(datadir)

    if ((pargs.code_version == 3) and pargs.wscaling):
        figname = 'hyb_weak_scaling'
        fig = phyb.plot_weak_scaling(datadir)

# --------------------------------------------------------------------------- #
    
    # STRONG SCALING
    if pargs.code_version == 2 and pargs.sscaling:
        figname = 'mpi_strong_scaling'
        fig = pmpi.plot_strong_scaling(datadir)

    if pargs.code_version == 3 and pargs.sscaling_thr:
        figname = 'hyb_strong_scaling_thr'
        fig = phyb.plot_strong_scaling_threads(datadir)

    if pargs.code_version == 3 and pargs.sscaling:
        figname = 'hyb_strong_scaling'
        fig = phyb.plot_strong_scaling(datadir)

# --------------------------------------------------------------------------- #
    
    # TIME-TO-SOLUTION vs N
    if pargs.code_version == 0 and pargs.serial_N:
        figname = 'srl_timetosol'
        fig = psrl.plot_timetosol(datadir)

    if pargs.runtime_comp:
        figname = 'runtime_comp'
        fig = prun.plot_all_runtimes(datadir)

# --------------------------------------------------------------------------- #
    
    # SPEEDUP LAWS
    if pargs.laws:
        figname = 'speedup_laws'
        fig = plaws.plot_scaling_laws()

# --------------------------------------------------------------------------- #
    if fig:
        plt.show()

        if pargs.save:
            if save_in_report: figdir = os.path.join(reportdir, figdir)
            prp.savefig(fig, figdir, figname)


