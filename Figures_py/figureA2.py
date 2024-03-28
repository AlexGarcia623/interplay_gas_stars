## This is very similar to figure4.py
import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm,ListedColormap

import sys, os
sys.path.append(os.path.dirname(os.getcwd()))

from interplay_gas_stars.plot_functions import (
    get_Z_Mstar_SFR, fixed_M_bins, ztoSnaps, sSFRcut
)

SAVEDIR = './Figures (pdf)/'

SIMS       = ['TNG','ORIGINAL','EAGLE']
SIMS_NAMES = [r'${\rm TNG}$',r'${\rm Illustris}$',r'${\rm EAGLE}$']

redshift = 1

SNAPS     = ztoSnaps[redshift]
CMIN,CMAX = sSFRcut[redshift]
dirs      = ['./Data/%s/snap%s/' %(SIMS[i],SNAPS[i]) for i in range(len(SIMS)) ]

fig, axs = plt.subplots(1, 3, figsize=(10,3.5), sharey=True, sharex=True)

mass_spacing = 1.0
mbins = np.arange( 8.0,12.0,mass_spacing )

newcolors = plt.cm.Spectral(np.linspace(0, 1, len(mbins)))
newcmp = ListedColormap(newcolors)

for index, ax in enumerate(axs):
    Zstar, Mstar, sSFR = get_Z_Mstar_SFR( dirs[index], which='stars' )
    
    sSFR = np.log10(sSFR)
    
    Hist1, xedges, yedges = np.histogram2d(sSFR,Zstar,weights=Mstar,bins=(75,75))
    Hist2, _     , _      = np.histogram2d(sSFR,Zstar,bins=[xedges,yedges])

    Hist1 = np.transpose(Hist1)
    Hist2 = np.transpose(Hist2)

    hist = Hist1/Hist2
    
    plot = ax.pcolormesh(xedges,yedges,hist,cmap=newcmp,vmin=8.0,vmax=12.0)
    
    ax.text( 0.075, 0.9, SIMS_NAMES[index], transform=ax.transAxes )
    
    ax.set_xlabel(r'$\log({\rm sSFR~[yr]}^{-1})$')
    
    fixed_mass, fixed_Z, fixed_sSFR = fixed_M_bins(Mstar,Zstar,sSFR,mass_spacing,nbins=10)
        
    for mass in range( len( fixed_mass ) ):
        color = newcmp(mass/len(fixed_mass))
        currentmass = r'$%s$' %( round( (fixed_mass[mass]),2 )+mass_spacing/2 )
        
        mask = (~np.isnan(fixed_sSFR[mass]) & 
                ~np.isnan(fixed_Z[mass]))
        
        x_axis = fixed_sSFR[mass][mask]
        y_axis = fixed_Z   [mass][mask]
            
        if (len(x_axis) != 0):
            polyfit   = np.polyfit( x_axis, y_axis, 1 )

            y_to_plot = np.polyval(polyfit,x_axis)

            ax.plot( x_axis, y_to_plot, label=currentmass, lw=6,
                     color = 'k' )
            ax.plot( x_axis, y_to_plot, label=currentmass, lw=3,
                     color = color )
    
fig.tight_layout()

axs[0].set_ylabel(r'$\log(Z_*~[Z_\odot])$')

ax_cbar = fig.add_axes([0.2, 1., 0.6, 0.05])

cb = plt.colorbar(plot, cax=ax_cbar,
                          shrink=0.5,orientation='horizontal')
cb.set_label(r'$\log (M_*~[M_\odot])$')

cb.ax.xaxis.set_label_position('top')

axs[1].text( 0.75, 0.85, r'$z = %s$' %redshift, transform=axs[1].transAxes )

ymin,ymax = axs[0].get_ylim()
axs[0].set_ylim( ymin, ymax+0.3 )
axs[0].set_xlim( -10.25, -8.5 )

fig.savefig( SAVEDIR + 'FigureA2.pdf', bbox_inches="tight" )