import sys
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm,ListedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from plot_functions import get_Z_Mstar_SFR, ztoSnaps, sSFRcut

sys.path.insert(1,'./Data/')

from additional_data import dr18_masses, dr18_metals, dr18_avg_sSFR

mpl.rcParams['font.size']=18

SAVEDIR = './Figures (pdf)/'

SIMS       = ['TNG','ORIGINAL','EAGLE']
SIMS_NAMES = [r'${\rm TNG}$',r'${\rm Illustris}$',r'${\rm EAGLE}$']

redshift = 0

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

    if ax == axs[2]:
        axins = ax.inset_axes([0.35, 0.125, 0.6, 0.525], transform=None)

        norm = mpl.colors.Normalize(vmin=CMIN, vmax=CMAX)

        axins.plot(dr18_masses[ :4],dr18_metals[ :4], lw=2.25, color='k')
        axins.plot(dr18_masses[4:9],dr18_metals[4:9], lw=2.25, color='k')
        axins.plot(dr18_masses[9: ],dr18_metals[9: ], lw=2.25, color='k')
        axins.plot(dr18_masses[ :4],dr18_metals[ :4], lw=2, color=newcmp(norm(dr18_avg_sSFR[0])))
        axins.plot(dr18_masses[4:9],dr18_metals[4:9], lw=2, color=newcmp(norm(dr18_avg_sSFR[1])))
        axins.plot(dr18_masses[9: ],dr18_metals[9: ], lw=2, color=newcmp(norm(dr18_avg_sSFR[2])))

        xmin, xmax = axins.get_xlim()
        ymin, ymax = axins.get_ylim()

        axins.pcolormesh(xedges, yedges, np.log10(hist),cmap=newcmp,vmin=CMIN,vmax=CMAX, alpha=0.33)

        axins.set_xlim(xmin, xmax)
        axins.set_ylim(ymin, ymax)
        axins.set_xticks([9,10,11])
        axins.set_yticks([])
        t=axins.text(0.60,0.1,r"${\rm De\;Rossi+(2018)}$",fontsize=11,transform=axins.transAxes,ha='center')
        t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))
        ax.indicate_inset_zoom(axins, edgecolor="black", alpha=0.25)
    
fig.tight_layout()

p0 = axs[0].get_position().get_points().flatten()
p1 = axs[1].get_position().get_points().flatten()
p2 = axs[2].get_position().get_points().flatten()
ax_cbar = fig.add_axes([0.2, 1.02, 0.6, 0.05])

cb = plt.colorbar(plot, cax=ax_cbar, ticks=np.linspace(CMIN,CMAX,spacing+1),
                          shrink=0.5,orientation='horizontal')
cb.set_label(r'$\log({\rm sSFR}~[{\rm yr}^{-1}])$')

# cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')

axs[0].set_ylabel(r'$\log(Z_*~[Z_\odot])$')

axs[1].text( 0.75, 0.1, r'$z = %s$' %redshift, transform=axs[1].transAxes )

fig.savefig( SAVEDIR + 'Figure3.pdf', bbox_inches="tight" )



sub_res_ztoSnaps = {
    0 :[99,99,99,135,135,135],
    1 :[50,50,50,86,86,86],
    2 :[33,33,33,68,68,68],
    3 :[25,25,25,60,60,60],
    4 :[21,21,21,54,54,54],
    5 :[17,17,17,49,49,49],
    6 :[13,13,13,45,45,45],
    7 :[11,11,11,41,41,41],
    8 :[8 ,8 ,8 ,38,38,38],
    9 :[6 ,6 ,6 ,35,35,35],
    10:[4 ,4 ,4 ,32,32,32]
}
sub_res_snaps = sub_res_ztoSnaps[redshift]

sub_res = ['TNG','TNG-2','TNG-3','ORIGINAL','ORIGINAL-2','ORIGINAL-3']
sub_res_names = [r'${\rm TNG}$',r'${\rm TNG}-2$',r'${\rm TNG}-3$',
                 r'${\rm Illustris}$',r'${\rm Illustris}-2$',r'${\rm Illustris}-3$']
sub_res_dirs = [ '../blue_FMR/%s/data/snap%s/' %(sub_res[i],sub_res_snaps[i]) for i in range(len(sub_res)) ]