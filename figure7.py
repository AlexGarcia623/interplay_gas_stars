import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,ListedColormap

from plot_functions import get_Z_Mstar_SFR, getScatter, skew, ztoSnaps, sSFRcut
from helpers import estimateSkew, getMedians
from getAlpha import switch_sim

SAVEDIR = './Figures (pdf)/'

fig = plt.figure(figsize=(18,9))
widths=[1,1,1]
heights=[1,4]
gs = gridspec.GridSpec(ncols = 3, nrows = 2, width_ratios = widths, height_ratios = heights, wspace=0.00, hspace=0.00,
                       left=0.1, right=0.9, bottom=0.1, top=0.9)

ALL_SIMS = ['TNG',"ORIGINAL",'EAGLE']

ax_big_left = None
labelsOff   = False

# 0 -> z=0, 1 -> z=1, 2 -> z=2
snap_idx = 0 # same as redshift

for index, sim in enumerate(ALL_SIMS):

    DATA     = './Data/' + sim + '/'

    ax       = fig.add_subplot(gs[1,index])
    ax_histx = fig.add_subplot(gs[0,index],sharex=ax)

    snapshots, snap2z, BLUE_DIR = switch_sim(sim)

    if (index > 0):
        labelsOff = True
    
    skew(DATA,ax_histx,ax,sim,snapshots[snap_idx], snap2z, labelsOff)

    if index == 1:
        ax.text( 0.15,0.035,r'${\rm one\!-\!to\!-\!one}$', transform=ax.transAxes,color='gray',fontsize=14 )
    
    ax_histx.axis('off')
    
plt.tight_layout()

plt.savefig(SAVEDIR + 'Figure7.pdf', bbox_inches='tight')