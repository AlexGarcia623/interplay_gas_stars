import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from plot_functions import get_Z_Mstar_SFR, medianZR, ztoSnaps, sSFRcut

sys.path.insert(1,'./Data/')

from additional_data import (G05_masses, G05_metals, C19_masses, C19_metals, K22_masses, K22_metals)

SAVEDIR = './Figures (pdf)/'

SIMS       = ['TNG','ORIGINAL','EAGLE']
SIMS_NAMES = [r'${\rm TNG}$',r'${\rm Illustris}$',r'${\rm EAGLE}$']

redshift = 0

SNAPS     = ztoSnaps[redshift]
CMIN,CMAX = sSFRcut[redshift]
dirs      = ['./Data/%s/snap%s/' %(SIMS[i],SNAPS[i]) for i in range(len(SIMS)) ]


fig, axs = plt.subplots(1, 3, figsize=(10,3.5), sharey=True, sharex=True)

bins = 60

zs=[0,1,2,3,4,5,6,7,8,9,10]
for redshift in zs:
    if (redshift > 8):
        continue
    
    SNAPS = ztoSnaps[redshift]

    dirs = ['./Data/%s/snap%s/' %(SIMS[i],SNAPS[i]) for i in range(len(SIMS)) ]

    for index, ax in enumerate(axs):
        Zstar, Mstar, sSFR = get_Z_Mstar_SFR( dirs[index], which='stars' )

        bin_centers, bin_vals, bin_errs = medianZR( Mstar, Zstar, bins=10, return_err_bars=True )

        if (index==0):
            text = 'TNG'
        elif (index==1):
            text = 'Illustris'
        else:
            text = 'EAGLE'
                
        color = 'C' + str(redshift)
        
        if (redshift == 0):
            ax.plot( G05_masses, G05_metals, color='k', linestyle='--' )
            ax.plot( K22_masses, K22_metals, color='k', linestyle=':' )
            ax.scatter( C19_masses, C19_metals, color='k' )
        
        ax.fill_between( bin_centers, bin_vals - bin_errs, bin_vals + bin_errs, alpha=0.15 )
        ax.plot( bin_centers, bin_vals, lw=3, color=color, label=r'$z=%s$' %redshift )
        
        if (redshift == 0):
            ax.text( 0.075, 0.9, SIMS_NAMES[index], transform=ax.transAxes )
        ax.set_xlabel(r'$\log\left(M_* ~[M_\odot]\right)$')

ymin,ymax = axs[0].get_ylim()
axs[0].set_ylim( ymin*1.1, ymax )
fig.tight_layout()

axs[0].set_ylabel(r'$\log(Z_*~[Z_\odot])$')

leg = axs[2].legend(frameon=False,labelspacing=0.05,
                    handletextpad=0, handlelength=0, markerscale=-1,bbox_to_anchor=(1,1))

for i in range(len(leg.get_texts())): leg.legendHandles[i].set_visible(False)

colors = ['C' + str(i) for i in range(0,9)] + ['k'] + ['teal']
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])

legend_elements = [
    Line2D([0], [0], color='k', linestyle='--', lw=1, 
           label=r'${\rm Gallazzi+(2005)}$'),
    Line2D([0], [0], color='k', linestyle=':', lw=1, 
           label=r'${\rm Kashino+(2022)}$'),
    Line2D([0], [0], marker='o', color='k',
           label=r'${\rm Cullen\,\:\:\;\!+(2019)}$',linestyle='',
           markersize=3)
]
    
legend = axs[0].legend(handles=legend_elements, loc='lower right', fontsize=10,
                       frameon=False,handletextpad=0.5,labelspacing=0.05)

legend_elements = [
    Line2D([0], [0], color='k', linestyle='--', lw=1, 
           label=r'$z=0$'),
    Line2D([0], [0], color='k', linestyle=':', lw=1, 
           label=r'$1.6<z<3.0$'),
    Line2D([0], [0], marker='o', color='k',
           label=r'$2.5<z<5.0$',linestyle='',
           markersize=3),
]

axs[1].legend(handles=legend_elements, loc='lower right', fontsize=10,
              frameon=False,handletextpad=0.5,labelspacing=0.05)

fig.savefig( SAVEDIR + 'Figure2.pdf', bbox_inches="tight" )