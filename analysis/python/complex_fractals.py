# Simply to take quick glance at the Mandelbrodt set project
# (animated plot with varying n_max)

from __future__ import print_function

import math as ma
import matplotlib.animation as manim
import matplotlib.cm as cm
import matplotlib.colors as col
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from matplotlib.gridspec import GridSpec

import paramplots as prp
import utils as u

fractal_mode = 'Mandelbrodt'
# fractal_mode = 'Julia'

# julia_C = -0.7 + 0.27015j
julia_C = -0.8j
# julia_C = 1.0 - (1.0+ma.sqrt(5.0))/2.0
# julia_C = -0.7269+0.1889j

N_points = 512 # resolution of the image
fixed_n_max = 6e1

# color_map = cm.jet
temp_jet = cm.jet(np.linspace(0, 1, 256))
# color_map.set_under('black')
# color_map.set_over('white')
temp_black = cm.Greys(np.linspace(1, 1, 1))
colors = np.vstack((temp_black, temp_jet))
color_map = col.LinearSegmentedColormap.from_list('colormap', colors)

show_plot = True
save_movie = False
save_plot = False

def iter_polyquad(C, z_init=None, n_max=40, z_thresh=2, mode='std'):
    if z_init is None:
        z = 0.0
    else:
        z = z_init

    if mode == 'dist':
        dz = 0.0

    for n in xrange(1, n_max+1):
        z_next = z**2 + C
        if mode == 'dist':
            dz = 2.0 * z * dz + 1.0 # derivative
        z = z_next
        if abs(z) > z_thresh:
            return n

    if mode == 'dist':
        dist = ma.log(abs(z**2))*abs(z)/abs(dz)
        if dist < 1.0e-4: return -1

    return -1


if __name__ == '__main__':

    figsize = (6, 7)

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(2, 1, height_ratios=[0.15, 4])
    gs.update(wspace=0.05, hspace=0.3) # set the spacing between axes

    ax1 = plt.subplot(gs[1])
    ax1.set_xlabel(r"$\rm{Re}$"+r"$(c)$", fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{Im}$"+r"$(c)$", fontsize=prp.labelsize_tex)

    # set the range of max number of iterations per pixel
    n_max_rge = np.logspace(1, 1.8, 6, dtype=int)       # animated
    # n_max_rge = np.array([fixed_n_max,], dtype=int)   # plain

    print("iterations :", n_max_rge.size)

    X = np.linspace(-2.0, 1.0, N_points)
    Y = np.linspace(-1.5, 1.5, N_points)
    Z = np.zeros((Y.size, X.size))

    ims = []
    for i, n_max in enumerate(n_max_rge):
        print("|", sep='', end='') # progress bar
        for iy, y in enumerate(Y):

            for ix, x in enumerate(X):

                if fractal_mode == 'Julia':
                    Z[iy, ix] = iter_polyquad(julia_C, x + 1j*y, n_max=n_max)

                elif fractal_mode == 'Mandelbrodt':
                    Z[iy, ix] = iter_polyquad(x + 1j*y, n_max=n_max)

        im = ax1.imshow(Z, cmap = color_map,
                       norm=col.Normalize(clip=False),
                       interpolation='none',
                       extent = (X.min(), X.max(), Y.min(), Y.max()))
        ims.append([im])
    ims = ims[1:] # looks better on colored plot

    if n_max_rge.size > 1 and not save_plot:
        # enables animation
        ani = manim.ArtistAnimation(fig, ims, interval=500, blit=True, 
                                    repeat=True, repeat_delay=800)
    else:
        pass

    if n_max_rge.size == 1:
        ax2 = plt.subplot(gs[0])
        cb  = fig.colorbar(im, cax=ax2, 
                           orientation='horizontal', 
                           ticks=[0, int(fixed_n_max/2), fixed_n_max])
        ax2.set_title(r"$\rm{number\ of\ iterations}$", 
                      fontsize=prp.labelsize_tex)

    ax1.text(0.04, 0.96, r"$N={}$".format(N_points),
             verticalalignment='top', horizontalalignment='left', 
             color='white', fontsize=prp.labelsize_tex, 
             bbox={'facecolor': 'black', 'alpha': 0, 'pad': 6},
             transform=ax1.transAxes)
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_float_p1))
    ax1.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_float_p1))

    if fractal_mode == 'Julia':
        figname = 'fractals_Julia_{}+{}i_{}_{}.mp4'.format(julia_C.real, 
                                                           julia_C.imag,
                                                           n_max_rge.size, 
                                                           N_points)
    elif fractal_mode == 'Mandelbrodt':
        figname = 'fractals_Mandelbrodt_{}_{}.mp4'.format(n_max_rge.size,
                                                          N_points)
    if save_movie:
        ani.save(figname, write=write)

    if show_plot:
        plt.show()


