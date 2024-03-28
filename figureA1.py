### This is virtually the same as figure3.py

import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm,ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from plot_functions import get_Z_Mstar_SFR, ztoSnaps, sSFRcut


mpl.rcParams['font.size']=18

SAVEDIR = '../Figures (pdf)/'

SIMS       = ['TNG','ORIGINAL','EAGLE']
SIMS_NAMES = [r'${\rm TNG}$',r'${\rm Illustris}$',r'${\rm EAGLE}$']

redshift = 1

SNAPS     = ztoSnaps[redshift]
CMIN,CMAX = sSFRcut[redshift]
dirs      = ['./Data/%s/snap%s/' %(SIMS[i],SNAPS[i]) for i in range(len(SIMS)) ]


fig, axs = plt.subplots(1, 3, figsize=(10,3.5), sharey=True, sharex=True)

bins = 75

spacing = 5
color_bins = np.linspace( CMIN,CMAX,spacing )

newcolors = plt.cm.viridis(np.linspace(0, 1, len(color_bins)))
newcmp = ListedColormap(newcolors)

for index, ax in enumerate(axs):
    Zstar, Mstar, sSFR = get_Z_Mstar_SFR( dirs[index], which="stars" )
    
    Hist1, xedges, yedges = np.histogram2d(Mstar,Zstar,weights=sSFR,bins=(bins,bins))
    Hist2, _     , _      = np.histogram2d(Mstar,Zstar,bins=[xedges,yedges])

    Hist1 = np.transpose(Hist1)
    Hist2 = np.transpose(Hist2)

    hist = Hist1/Hist2
    
    plot = ax.pcolormesh(xedges,yedges,np.log10(hist),cmap=newcmp,vmin=CMIN,vmax=CMAX)
    
    ax.text( 0.075, 0.85, SIMS_NAMES[index], transform=ax.transAxes )
    
    ax.set_xlabel(r'$\log\left(M_* ~[M_\odot]\right)$')
    
fig.tight_layout()

p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten()
p2 = axs[2].get_position().get_points().flatten()
ax_cbar = fig.add_axes([0.2, 1.02, 0.6, 0.05])

cb = plt.colorbar(plot, cax=ax_cbar, ticks=np.linspace(CMIN,CMAX,spacing+1),
                          shrink=0.5,orientation='horizontal')
cb.set_label(r'$\log({\rm sSFR}~[{\rm yr}^{-1}])$')

cb.ax.xaxis.set_label_position('top')

axs[0].set_ylabel(r'$\log(Z_*~[Z_\odot])$')

axs[1].text( 0.75, 0.1, r'$z = %s$' %redshift, transform=axs[1].transAxes )

axs[0].set_xlim(8.001,11.999)

fig.savefig( SAVEDIR + 'FigureA1.pdf', bbox_inches="tight" )