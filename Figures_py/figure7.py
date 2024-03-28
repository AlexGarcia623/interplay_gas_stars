'''
This file is used to create Figure 7 of "Interplay of Stellar
and Gas-Phase Metallicities: Unveiling Insights for Stellar 
Feedback Modeling with Illustris, IllustrisTNG, and EAGLE"

Paper: https://ui.adsabs.harvard.edu/abs/2024MNRAS.tmp..787G/abstract

Code written by: Alex Garcia, 2023-24
'''
# Standard Imports
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm,ListedColormap
# Import from this library
import sys, os
sys.path.append(os.path.dirname(os.getcwd()))

from interplay_gas_stars.plot_functions import (
    get_Z_Mstar_SFR, getScatter, skew, ztoSnaps, sSFRcut
)
from interplay_gas_stars.helpers import (
    estimateSkew, getMedians
)
from interplay_gas_stars.getAlpha import (
    switch_sim
)

SAVEDIR = './Figures (pdf)/' # Where to save files

fig = plt.figure(figsize=(18,9))
widths=[1,1,1]
heights=[1,4]
gs = gridspec.GridSpec(ncols = 3, nrows = 2, width_ratios = widths, height_ratios = heights, wspace=0.00, hspace=0.00,
                       left=0.1, right=0.9, bottom=0.1, top=0.9)

# Simulation names (same as Data directories)
ALL_SIMS = ['TNG',"ORIGINAL",'EAGLE']

ax_big_left = None
labelsOff   = False

redshift = 0
for index, sim in enumerate(ALL_SIMS):
    DATA     = './Data/' + sim + '/'

    ax       = fig.add_subplot(gs[1,index])
    ax_histx = fig.add_subplot(gs[0,index],sharex=ax)

    snapshots, snap2z, BLUE_DIR = switch_sim(sim)

    if (index > 0):
        labelsOff = True
    
    skew(DATA,ax_histx,ax,sim,snapshots[redshift], snap2z, labelsOff)

    if index == 1:
        ax.text( 0.15,0.035,r'${\rm one\!-\!to\!-\!one}$', transform=ax.transAxes,color='gray',fontsize=14 )
    
    ax_histx.axis('off')
    
plt.tight_layout()
plt.savefig(SAVEDIR + 'Figure7.pdf', bbox_inches='tight')