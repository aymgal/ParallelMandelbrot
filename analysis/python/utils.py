import numpy as np


# ---------------------------------------

def ticks_int(value, index):
    return r'${:d}$'.format(int(value))

def ticks_float_p2(value, index):
    return r'${:.2f}$'.format(value)

def ticks_base10(value, index):
    """
    This function decompose value in base*10^{exp} and return a latex string.
    If 0<=value<99: return the value as it is.
    if 0.1<value<0: returns as it is rounded to the first decimal
    otherwise returns $base*10^{exp}$
    I've designed the function to be use with values for which the decomposition
    returns integers.
    From http://stackoverflow.com/questions/19239297/matplotlib-bad-ticks-labels-for-loglog-twin-axis
    """
    exp = np.floor(np.log10(value))
    base = value/10**exp
    if exp == 0 or exp == 1:
        return '${0:d}$'.format(int(value))
    if exp == -1:
        return '${0:.1f}$'.format(value)
    else:
        return '${0:d}\\cdot10^{{{1:d}}}$'.format(int(base), int(exp))

# ---------------------------------------

def identity(p):
    return p

def amdahl_law(p, alpha):
    return 1.0 / (alpha + (1.0 - alpha)/p)

def gustafson_law(p, beta):
    return p - beta * (p - 1.0)

# ---------------------------------------

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

def get_data(datapath, sort_with=None):
    Nx, Ny, n_iter, n_row, \
    n_proc, n_thread, time = np.loadtxt(datapath, delimiter=' ', unpack=True)

    # because output from runs on cluster 
    # are not necessarily in the right order (increasing number of proc)
    if sort_with:
        if sort_with == 'n_proc':
            sort_idx = np.argsort(n_proc)
        elif sort_with == 'n_thread':
            sort_idx = np.argsort(n_thread)
        elif sort_with == 'n_iter':
            sort_idx = np.argsort(n_iter)
        elif sort_with == 'Nx':
            sort_idx = np.argsort(Nx)
    else:
        sort_idx = np.arange(Nx.size)

    # since every run is done 10 times, we must average values
    Nx       = average_every_n(Nx[sort_idx],       n=10)
    Ny       = average_every_n(Ny[sort_idx],       n=10)
    n_iter   = average_every_n(n_iter[sort_idx],   n=10)
    n_row    = average_every_n(n_row[sort_idx],    n=10)
    n_proc   = average_every_n(n_proc[sort_idx],   n=10)
    n_thread = average_every_n(n_thread[sort_idx], n=10)
    time, time_std = average_every_n(time[sort_idx], n=10, std_too=True)

    return (Nx, Ny, n_iter, n_row, n_proc, n_thread, time, time_std)


