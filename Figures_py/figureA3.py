import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm,ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys, os
sys.path.append(os.path.dirname(os.getcwd()))

from interplay_gas_stars.plot_functions import (
    get_Z_Mstar_SFR, ztoSnaps, sSFRcut
)

SAVEDIR = './Figures (pdf)/'

fig, axs = plt.subplots(1, 3, figsize=(10,3.5), sharey=True, sharex=True)

##### CHANGE ME TO GET DIFFERENT REDSHIFT #####
_100_50_snap_ = 99
_100_50_redshift_ = 0
###############################################

spacing = 5
CMIN,CMAX  = sSFRcut[_100_50_redshift_]
color_bins = np.linspace( CMIN,CMAX,spacing )

newcolors = plt.cm.viridis(np.linspace(0, 1, len(color_bins)))
newcmp    = ListedColormap(newcolors)

_100_50_comp_dirs_ = ['./Data/TNG/snap%s/' %_100_50_snap_,
                      './Data/TNG50-1/snap%s/' %_100_50_snap_,
                      './Data/TNG50-2/snap%s/' %_100_50_snap_]
_100_50_sim_names_ = [r'${\rm TNG}100-1$',r'${\rm TNG}50-1$',r'${\rm TNG}50-2$']

for index, ax in enumerate(axs):
    Zstar, Mstar, sSFR = get_Z_Mstar_SFR( _100_50_comp_dirs_ [index], which='stars')

    Hist1, xedges, yedges = np.histogram2d(Mstar,Zstar,weights=sSFR,bins=(75,75))
    Hist2, _     , _      = np.histogram2d(Mstar,Zstar,bins=[xedges,yedges])

    Hist1 = np.transpose(Hist1)
    Hist2 = np.transpose(Hist2)

    hist = Hist1/Hist2

    _c_min_, _c_max_ = sSFRcut[_100_50_redshift_]
    plot = ax.pcolormesh(xedges,yedges,np.log10(hist),cmap=newcmp,vmin=_c_min_,vmax=_c_max_)

    ax.text( 0.075, 0.9, _100_50_sim_names_[index], transform=ax.transAxes )

    ax.set_xlabel(r'$\log\left(M_* ~[M_\odot]\right)$')

fig.tight_layout()

p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten()
p2 = axs[2].get_position().get_points().flatten()
ax_cbar = fig.add_axes([0.2, 1.02, 0.6, 0.05])

cb = plt.colorbar(plot, cax=ax_cbar, ticks=np.linspace(_c_min_, _c_max_,spacing+1),
                          shrink=0.5,orientation='horizontal')
cb.set_label(r'$\log({\rm sSFR}~[{\rm yr}^{-1}])$')

cb.ax.xaxis.set_label_position('top')

axs[0].set_ylabel(r'$\log(Z_*~[Z_\odot])$')

axs[1].text( 0.75, 0.1, r'$z = %s$' %_100_50_redshift_, transform=axs[1].transAxes )

ymin,ymax = axs[0].get_ylim()
axs[0].set_ylim(ymin,ymax*1.1)

fig.savefig( SAVEDIR + 'FigureA3_z=%s.pdf' %_100_50_redshift_, bbox_inches="tight" )