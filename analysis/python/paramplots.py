import os
import numpy as np
import matplotlib as mpl

globalfsize = 14

figsize       = (8, 6)
linewidth     = 2
fontsize      = globalfsize
legendsize    = globalfsize
labelsize     = globalfsize
labelsize_tex = 18
titlesize     = globalfsize
ticksize      = globalfsize
fmt = 'pdf'

# LaTeX symbols shortcuts
mfH2 = r"$f_{\mathrm{H}_2}$"

mpl.rcParams['axes.labelpad']        = 5
mpl.rcParams['axes.labelsize']       = labelsize
mpl.rcParams['axes.linewidth']       = 1
mpl.rcParams['axes.titlesize']       = titlesize
mpl.rcParams['figure.figsize']       = figsize
mpl.rcParams['font.family']          = 'serif'
# mpl.rcParams['font.serif']           = ['Palatino']
mpl.rcParams['font.size']            = fontsize
mpl.rcParams['grid.color']           = '#8A8A8A'
mpl.rcParams['grid.linewidth']       = 0.5
mpl.rcParams['legend.fontsize']      = legendsize
mpl.rcParams['legend.numpoints']     = 3
mpl.rcParams['legend.scatterpoints'] = 3
mpl.rcParams['legend.markerscale']   = 2.5
mpl.rcParams['lines.color']          = 'b'
mpl.rcParams['lines.linewidth']      = linewidth
mpl.rcParams['patch.linewidth']      = 1
mpl.rcParams['savefig.bbox']         = 'tight'
mpl.rcParams['savefig.dpi']          = 300
mpl.rcParams['savefig.format']       = fmt
mpl.rcParams['savefig.pad_inches']   = 0.05
mpl.rcParams['text.usetex']          = False # buuuug
mpl.rcParams['xtick.labelsize']      = ticksize
mpl.rcParams['xtick.major.pad']      = 10
mpl.rcParams['xtick.major.size']     = 6
mpl.rcParams['xtick.major.width']    = 1
mpl.rcParams['xtick.minor.pad']      = 10
mpl.rcParams['xtick.minor.size']     = 4
mpl.rcParams['xtick.minor.width']    = 1
mpl.rcParams['ytick.labelsize']      = ticksize
mpl.rcParams['ytick.major.pad']      = 10
mpl.rcParams['ytick.major.size']     = 6
mpl.rcParams['ytick.major.width']    = 1
mpl.rcParams['ytick.minor.pad']      = 10
mpl.rcParams['ytick.minor.size']     = 4
mpl.rcParams['ytick.minor.width']    = 1

#### FIXES FOR 2.0
# mpl.style.use('classic') # reset all back to 1.x version
#** math mode
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm']      = 'serif'
#** tick placements
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top']       = True
mpl.rcParams['ytick.right']     = True
#** legend styling
mpl.rcParams['legend.fancybox']      = False
mpl.rcParams['legend.loc']           = 'upper right'
mpl.rcParams['legend.numpoints']     = 2
mpl.rcParams['legend.fontsize']      = 'large'
mpl.rcParams['legend.framealpha']    = None
mpl.rcParams['legend.scatterpoints'] = 3
mpl.rcParams['legend.edgecolor']     = 'inherit'
#** figure display
mpl.rcParams['figure.dpi'] = 80

#
def savefig(figobj, figdir, figname):
	savedir = figdir
	if not os.path.exists(savedir):
	    os.makedirs(savedir)
	savepath = os.path.join(savedir, figname)
	print "Saving figure to {}.{} ...".format(savepath, fmt)
	figobj.savefig(savepath)
	print "Figure saved."

def figsize(scale):
    fig_width_pt = 483.697                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width, fig_height]
    return fig_size


