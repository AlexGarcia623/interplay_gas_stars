import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from getAlpha import get_alpha

m_star_min, m_star_max = 8.0, 12.0

SAVEDIR = './Figures (pdf)/'

mpl.rcParams['font.size'] = 20
fig, axs = plt.subplots(3, 1, figsize=(8,8), sharex=True, sharey=True)

thresholds = [-1.00E-01,-5.00E-01,-1.00E+0,-np.inf]
labels = [r'$0.1~{\rm dex}$',r'$0.5~{\rm dex}~({\rm fiducial})$',r'$1.0~{\rm dex}$',r'${\rm SFR} > 0$']
shape = ['^','o','*','v']

redshifts = np.arange(0,9)

for index, THRESHOLD in enumerate(thresholds):
    
    EAGLE, EAGLE_lower, EAGLE_upper = get_alpha( 'EAGLE', THRESHOLD=THRESHOLD,
                                                 m_star_min=m_star_min, m_star_max=m_star_max)
    TNG, TNG_lower, TNG_upper = get_alpha( 'TNG', THRESHOLD=THRESHOLD,
                                           m_star_min=m_star_min, m_star_max=m_star_max)
    ORIGINAL, ORIGINAL_lower, ORIGINAL_upper = get_alpha( 'ORIGINAL', THRESHOLD=THRESHOLD,
                                                          m_star_min=m_star_min, m_star_max=m_star_max)
    
    axs[0].errorbar( redshifts, EAGLE, yerr=[EAGLE - EAGLE_lower,EAGLE_upper - EAGLE], marker=shape[index],
                     markersize=7, capsize=3, label=labels[index] )
    axs[1].errorbar( redshifts, TNG, yerr=[TNG - TNG_lower,TNG_upper - TNG], marker=shape[index],
                     markersize=7, capsize=3, label=labels[index] )
    axs[2].errorbar( redshifts, ORIGINAL, yerr=[ORIGINAL - ORIGINAL_lower,ORIGINAL_upper - ORIGINAL], marker=shape[index],
                     markersize=7, capsize=3, label=labels[index] )

for ax in axs:
    ax.set_ylabel(r'$\alpha_*$')
axs[2].set_xlabel(r'${\rm Redshift}$')
    
axs[0].text(0.5,0.85,r'${\rm EAGLE}$',ha='center', transform=axs[0].transAxes)
axs[1].text(0.5,0.85,r'${\rm TNG}$',ha='center', transform=axs[1].transAxes)
axs[2].text(0.5,0.85,r'${\rm Illustris}$',ha='center', transform=axs[2].transAxes)

axs[0].set_ylim(-0.15,1.15)
    
leg = axs[1].legend(frameon=False, fontsize=15, labelspace=0.0)
colors = ['C' + str(i) for i in range(0,9)] + ['k'] + ['teal']
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])
    
plt.tight_layout()
plt.subplots_adjust(hspace=0.0,wspace=0.0)
plt.savefig(SAVEDIR + 'FigureB1_alpha.pdf')