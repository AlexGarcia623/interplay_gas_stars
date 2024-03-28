'''
This file is used to create Figure 8 of "Interplay of Stellar
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
from getSlopes import get_slopes

mpl.rcParams['font.size']=16 # Change font size for this file

SAVEDIR = '../Figures (pdf)/' # Where to save files

fig = plt.figure(figsize=(8,4))

EAGLE = get_slopes('eagle')
ORIGINAL = get_slopes('original')
TNG = get_slopes('tng')

zs = np.arange(0,9)

size = 60
plt.scatter(zs, EAGLE   , label=r'${\rm EAGLE}$'    ,s=size,alpha=0.75,marker='o',edgecolor='k',linewidth=1.0)
plt.scatter(zs, ORIGINAL, label=r'${\rm Illustris}$',s=size,alpha=0.75,marker='^',edgecolor='k',linewidth=1.0)
plt.scatter(zs, TNG     , label=r'${\rm TNG}$'      ,s=size,alpha=0.75,marker='*',edgecolor='k',linewidth=1.0)

plt.ylabel(r'${\rm Slope}$')
plt.xlabel(r'${\rm Redshift}$')

## Emperically determined
hlines = [0.11 ,0.67  ,0.38 ,0.75  ,0.56  ,0.9   ,0.96   ,0.02  ]
vals   = [0.1*1,0.1*15,0.1*5,0.1*20,0.1*10,0.1*50,0.1*100,0.1/10]
color  = ['k' for _ in range(len(hlines))]


for index, line in enumerate(hlines):
    plt.axhline(line, color=color[index],linestyle='--',alpha=0.5)
    
    text = r"$%s$" %vals[index]
    fs=13
    if (index == len(hlines)-1):
        text =  r'$\Gamma$ = %s' %text
        plt.text( 1.708+1, line + 0.01 , text, color=color[index], alpha=0.75, fontsize=fs )
    else:
        if (index == 0):
            text += r' ${\rm (T18)}$'
        plt.text( 2.25+1, line + 0.01 , text, color=color[index], alpha=0.75, fontsize=fs )
        
leg = plt.legend(frameon=True,handletextpad=0.75, handlelength=0,labelspacing=0.01,
                 loc='lower right',framealpha=1,edgecolor=(1, 1, 1, 0))

colors = ['C0','C1','C2']
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])

plotlim = (plt.xlim()[0],plt.xlim()[1]+0.25,
           plt.ylim()[0],plt.ylim()[1])#plt.xlim() + plt.ylim()
plt.imshow([[1,1],[0,0]], cmap=plt.cm.binary, 
           interpolation='bicubic', extent=plotlim, alpha=0.5, aspect='auto')

plt.text( 0.02, 0.9 , r'${\rm Stronger~ Correlation}$', transform=plt.gca().transAxes)
plt.text( 0.02, 0.11, r'${\rm Weaker~ Correlation}$'  , transform=plt.gca().transAxes )

plt.xticks(np.arange(0,9))

plt.tight_layout()
plt.savefig( SAVEDIR + "Figure8.pdf", bbox_inches='tight' )