import matplotlib.pyplot as plt
import matplotlib.ticker as tick
import numpy as np
from matplotlib.gridspec import GridSpec

import paramplots as prp
import utils as u


def plot_scaling_laws():

    p_rge = np.logspace(0, 5, 1e2, base=2)

    fig = plt.figure(figsize=(13, 4))
    gs = GridSpec(2, 2, width_ratios=[0.4, 0.6], height_ratios=[1, 1])
    gs.update(wspace=0.3, hspace=0.1)

    ax1 = plt.subplot(gs[:, 0])
    ax1.loglog(p_rge, u.identity(p_rge), lw=2, ls=':', color='black')
    ax1.text(0.85, 0.94, "ideal",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.globalfsize, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax1.transAxes)
    ax1.loglog(p_rge, u.amdahl_law(p_rge, 0.05), ls='-', basex=2, basey=2,
               label=r"$\alpha=0.05$")
    ax1.loglog(p_rge, u.amdahl_law(p_rge, 0.35), ls='--', basex=2, basey=2,
               label=r"$\alpha=0.35$")
    ax1.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax1.set_ylim(0, p_rge.max())
    ax1.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", fontsize=prp.labelsize_tex)
    ax1.set_ylabel(r"$\rm{speedup},\ $"+r"$S_p$", fontsize=prp.labelsize_tex)
    ax1.set_title(r"$\rm{Amdahl's\ law}$", fontsize=prp.labelsize_tex)
    ax1.legend(loc='upper left', frameon=False, handletextpad=0.1)

    ax2 = plt.subplot(gs[0, 1], sharex=ax1)
    ax2.loglog(p_rge, u.identity(p_rge), lw=2, ls=':', color='black')
    ax2.text(0.8, 0.94, "ideal",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.globalfsize, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax2.transAxes)
    ax2.loglog(p_rge, u.gustafson_law(p_rge, 0.05), ls='-', basex=2, basey=2,
               label=r"$\beta=0.05$")
    ax2.loglog(p_rge, u.gustafson_law(p_rge, 0.35), ls='--', basex=2, basey=2,
               label=r"$\beta=0.35$")
    # ax2.yaxis.set_ticks([1, 4, 16, 64])
    ax2.yaxis.set_ticks([2, 8, 32])
    ax2.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax2.set_ylim(0, p_rge.max())
    ax2.set_ylabel(r"$\rm{speedup},\ $"+r"$S_p$", fontsize=prp.labelsize_tex)
    ax2.set_title(r"$\rm{Gustafson's\ law}$", fontsize=prp.labelsize_tex)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.legend(loc='upper left', frameon=False, handletextpad=0.1)

    ax3 = plt.subplot(gs[1, 1], sharex=ax2)
    ax3.semilogx(p_rge, u.identity(p_rge)/p_rge, lw=2, ls=':', color='black')
    ax3.text(0.94, 0.94, "ideal",
             verticalalignment='top', horizontalalignment='right', 
             color='black', fontsize=prp.globalfsize, 
             bbox={'facecolor': 'white', 'alpha': 0, 'pad': 6},
             transform=ax3.transAxes)
    ax3.semilogx(p_rge, u.gustafson_law(p_rge, 0.05)/p_rge, ls='-', 
               basex=2, label=r"$\beta(N)=0.05$")
    ax3.semilogx(p_rge, u.gustafson_law(p_rge, 0.35)/p_rge, ls='--', 
               basex=2, label=r"$\beta(N)=0.35$")
    ax3.xaxis.set_ticks([1, 2, 4, 8, 16, 32, 64])
    ax3.yaxis.set_ticks([0.5, 0.75, 1])
    ax3.xaxis.set_major_formatter(tick.FuncFormatter(u.ticks_int))
    ax3.yaxis.set_major_formatter(tick.FuncFormatter(u.ticks_float_p2))
    ax3.set_xlim(1, p_rge.max())
    ax3.set_ylim(0.5, 1.15)
    ax3.set_xlabel(r"$\rm{number\ of\ processors},\ $"+r"$p$", 
                   fontsize=prp.labelsize_tex)
    ax3.set_ylabel(r"$\rm{efficiency},\ E_p$", fontsize=prp.labelsize_tex)

    return fig


