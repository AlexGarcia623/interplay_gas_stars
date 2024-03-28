'''
This file is used to create Figure 1 of "Interplay of Stellar
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
# Import from this library
# from plot_functions import (
#     get_Z_Mstar_SFR, ztoSnaps, sSFRcut
# )
# from Data.additional_data import (
#     G05_masses, G05_metals, P08_masses, P08_metals
# )

SAVEDIR = './Figures (pdf)/' # Where to save files

# Simulation names (same as Data directories)
SIMS = ['TNG','ORIGINAL','EAGLE']
SIMS_NAMES = [r'${\rm TNG}$',r'${\rm Illustris}$',r'${\rm EAGLE}$'] # TeX

# Get constants and set-up directories
redshift = 0
snaps = ztoSnaps[redshift]
CMIN,CMAX = sSFRcut[redshift]
dirs = ['./Data/%s/snap%s/' %(SIMS[i],snaps[i]) for i in range(len(SIMS)) ]

# Create Figure
fig, axs = plt.subplots(1, 3, figsize=(10,3.5), sharey=True, sharex=True)

for index, ax in enumerate(axs):
    Zstar, Mstar, sSFR = get_Z_Mstar_SFR( dirs[index], which='stars' )
        
    Hist1, xedges, yedges = np.histogram2d(Mstar,Zstar,bins=(60,60))
    Hist1 = np.transpose(Hist1)
    
    plot = ax.pcolormesh(xedges,yedges,np.log10(Hist1),cmap=plt.cm.Blues)
    
    ax.text( 0.075, 0.9, SIMS_NAMES[index], transform=ax.transAxes )
    
    ax.set_xlabel(r'$\log\left(M_* ~[M_\odot]\right)$')
    
    ax.set_xlim(7.9,12.1)
    
    ax.plot( G05_masses, G05_metals, color='C2' )
    ax.plot( P08_masses, P08_metals, color='C1' )

fig.tight_layout()

ax_cbar = fig.add_axes([0.2, 1.02, 0.6, 0.05])

cb = plt.colorbar(plot, cax=ax_cbar, shrink=0.5,orientation='horizontal')
cb.set_label(r'$\log N_{\rm galaxies}$')

cb.ax.xaxis.set_label_position('top')

axs[0].set_ylabel(r'$\log(Z_*~[Z_\odot])$')

axs[0].text( 0.25, 0.08, r'${\rm Gallazzi~et~al.~ (2005)}$', transform=axs[0].transAxes,color='C2' )
axs[0].text( 0.3 , 0.18, r'${\rm Panter~et~al.~ (2008)}$'  , transform=axs[0].transAxes,color='C1' )
        
axs[1].text( 0.7, 0.1, r'$z = %s$' %redshift, transform=axs[1].transAxes )

fig.savefig( SAVEDIR + 'Figure1.pdf', bbox_inches="tight" )