import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from getSlopes import get_slopes

sys.path.insert(1,'./Data/')

from additional_data import (hlines, vals)

mpl.rcParams['font.size']=15

SAVEDIR = './Figures (pdf)/'


fig = plt.figure(figsize=(7,3.5))

gamma_vals  = vals[np.argsort(vals)]
pred_slopes = hlines[np.argsort(vals)]

plt.plot( np.log10(vals[np.argsort(vals)]), hlines[np.argsort(vals)], color='k', lw=2.5 )

ax = plt.gca()

EAGLE = get_slopes('eagle')
ORIGINAL = get_slopes('original')
TNG = get_slopes('tng')

sims  = [TNG,EAGLE,ORIGINAL]
clrs  = ['C2','C0','C1']
for index,slopes in enumerate(sims):
    
    sim = [
        np.log10(vals[ np.argmin(abs(np.max(slopes) - hlines ) ) ]),
        np.log10(vals[ np.argmin(abs(np.min(slopes) - hlines ) ) ])
    ]
    plt.axvline( np.min(sim), color=clrs[index], lw=2 )
    plt.axvline( np.max(sim), color=clrs[index], lw=2 )
    alpha = 0.3
    if index == 1:
        alpha = alpha+0.2
    plt.axvspan( np.min(sim), np.max(sim), alpha=alpha, color=clrs[index])


plt.text( 0.8, 0.5, r'${\rm EAGLE}$'    , transform=ax.transAxes, color='C0' )
plt.text( 0.8, 0.4, r'${\rm Illustris}$', transform=ax.transAxes, color='C1' )
plt.text( 0.8, 0.3, r'$\bf{\rm TNG}$'   , transform=ax.transAxes, color='C2' )

ymin, ymax = -0.05,1.1
plotlim = (np.min(np.log10(vals)),np.max(np.log10(vals)),
           ymin, ymax)#plt.xlim() + plt.ylim()
plt.imshow([[1,1],[0,0]], cmap=plt.cm.binary, interpolation='bicubic', extent=plotlim, alpha=0.5, aspect='auto')

plt.xlabel(r'$\log\Gamma = \log\left({\tau_{\rm c}}/{\tau_{\rm SF}}\right)$')
plt.ylabel(r'${\rm Slope}$')

plt.text( 0.02, 0.8 , r'${\rm Stronger~ Correlation}$', transform=ax.transAxes)
plt.text( 0.02, 0.21, r'${\rm Weaker~ Correlation}$', transform=ax.transAxes )

plt.ylim(ymin, ymax)

plt.tight_layout()

plt.savefig(SAVEDIR + 'Figure9.pdf', bbox_inches="tight")
plt.show()