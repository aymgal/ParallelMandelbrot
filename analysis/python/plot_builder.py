import argparse as argp
import glob
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

datadir = '../datafiles/'

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
    parser.add_argument('-v', '--version', dest='version', type=int,
                        metavar='\b', default=0, help=help_m)
    parser.add_argument('-t', '--timetosol', dest='eff',
                        help="plot of time to solution vs. number of procs", 
                        action='store_true')
    parser.add_argument('-s', '--speedup', dest='sup',
                        help="plot of speedup vs. number of procs", 
                        action='store_true')
    parser.add_argument('-e', '--efficiency', dest='eff',
                        help="plot of efficiency vs. number of procs", 
                        action='store_true')
    return parser.parse_args()


if __name__ == '__main__':

    parser_args = parser_arguments()

    datafiles = glob.glob(datadir+'*.dat')

    # N = np.zeros()
    # for datafile in datafiles:
    #     print 
