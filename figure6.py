'''
This file is used to create Figure 6 of "Interplay of Stellar
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
from getAlpha import get_alpha

SAVEDIR = '../Figures (pdf)/' # Where to save files

# Simulation names (same as Data directories)
sims = ['ORIGINAL','TNG','EAGLE']

m_star_min = 8.0
m_star_max = 12.0
m_gas_min  = 8.5

polyorder = 1 # Polynomial fit of line

EAGLE, EAGLE_lower, EAGLE_upper = get_alpha( 'EAGLE', m_star_min=m_star_min, m_star_max=m_star_max,
                                             polyorder=polyorder )
TNG, TNG_lower, TNG_upper = get_alpha( 'TNG', m_star_min=m_star_min, m_star_max=m_star_max,
                                       polyorder=polyorder )
ORIGINAL, ORIGINAL_lower, ORIGINAL_upper = get_alpha( 'ORIGINAL', m_star_min=m_star_min, m_star_max=m_star_max, 
                                                      polyorder=polyorder )

EAGLE_upper = EAGLE_upper - EAGLE
EAGLE_lower = EAGLE - EAGLE_lower

ORIGINAL_upper = ORIGINAL_upper - ORIGINAL
ORIGINAL_lower = ORIGINAL - ORIGINAL_lower

TNG_upper = TNG_upper - TNG
TNG_lower = TNG - TNG_lower


fig = plt.figure(figsize=(8,3.5))

z = np.arange(0,9)

plt.errorbar( z+0.00, ORIGINAL, label=r'${\rm Illustris}$',
              alpha=0.75, color='C1', yerr = [ORIGINAL_lower, ORIGINAL_upper],
              linestyle='none', marker='^',markersize=7)

plt.errorbar( z+0.00, TNG, label=r'${\rm TNG}$',
              alpha=0.75, color='C2', yerr = [TNG_lower, TNG_upper],
              linestyle='none', marker='*',markersize=7 )

plt.errorbar( z-0.00, EAGLE, label=r'${\rm EAGLE}$',
              alpha=0.75, color='C0', yerr = [EAGLE_lower, EAGLE_upper],
              linestyle='none', marker='o',markersize=7 )

## Make the legend (includes removing error bars) ##
leg  = plt.legend(frameon=True,handletextpad=0, handlelength=0,
                  markerscale=0,loc='lower right',labelspacing=0.05)
lCol = ['C1','C2','C0']
for n, text in enumerate( leg.texts ):
    text.set_color( lCol[n] )   
# get handles
handles, labels = plt.gca().get_legend_handles_labels()
# remove the errorbars
handles = [h[0] for h in handles]
# use them in the legend
leg = plt.legend(frameon=True,handletextpad=0.75, handlelength=0,labelspacing=0.01,
             loc='upper left',fontsize=14)
for n, text in enumerate( leg.texts ):
    text.set_color( lCol[n] )

leg.get_frame().set_alpha(0.0)
leg.get_frame().set_edgecolor('white')
####################################################

plt.xlabel(r'${\rm Redshift}$')
plt.ylabel(r'$\alpha_{\rm min}$')

plt.axhline(0.33 ,color='k',linestyle='--')
plt.text(0.8,0.375,r'${\rm M10 ~(Gas)}$',fontsize=14,alpha=0.5,transform=plt.gca().transAxes)
plt.axhline(0.66 ,color='k',linestyle='-.')
plt.text(0.8,0.7,r'${\rm AM13 ~(Gas)}$',fontsize=14,alpha=0.5,transform=plt.gca().transAxes)
plt.axhline(0.55 ,color='k',linestyle='solid')
plt.text(0.8,0.575,r'${\rm C20~(Gas)}$',fontsize=14,alpha=0.5,transform=plt.gca().transAxes)

ymin, _ = plt.ylim()
plt.ylim(ymin,1.)

plt.tight_layout()
plt.savefig(SAVEDIR+"Figure6.pdf", bbox_inches='tight')