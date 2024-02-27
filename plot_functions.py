import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import h5py

from scipy.interpolate import interp1d
from scipy.stats import norm, skewnorm

from helpers import sfmscut, getMedians

# Custom rcParams

mpl.rcParams['font.size'] = 15
mpl.rcParams['axes.linewidth'] = 1.5
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.minor.visible'] = 'true'
mpl.rcParams['ytick.minor.visible'] = 'true'
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['xtick.minor.width'] = 1.0
mpl.rcParams['ytick.minor.width'] = 1.0
mpl.rcParams['xtick.major.size'] = 7.5
mpl.rcParams['ytick.major.size'] = 7.5
mpl.rcParams['xtick.minor.size'] = 3.5
mpl.rcParams['ytick.minor.size'] = 3.5
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True    
# %set_env MANPATH=/home/paul.torrey/local/texlive/2018/texmf-dist/doc/man:$MANPATH
# %set_env INFOPATH=/home/paul.torrey/local/texlive/2018/texmf-dist/doc/info:$INFOPATH
# %set_env PATH=/home/paul.torrey/local/texlive/2018/bin/x86_64-linux:/home/paul.torrey/local/texlive/2018/texmf-dist:$PATH
mpl.rcParams['text.usetex']        = True
# mpl.rcParams['text.latex.unicode'] = True
mpl.rcParams['font.family']        = 'serif'
mpl.rc('font',**{'family':'sans-serif','serif':['Times New Roman'],'size':15})
mpl.rc('text', usetex=True)

m_star_min = 8.0
m_star_max = 12.0
m_gas_min  = 8.5

h      = 6.774E-01
xh     = 7.600E-01
zo     = 3.500E-01
mh     = 1.6726219E-24
kb     = 1.3806485E-16
mc     = 1.270E-02
Zsun   = 1.27E-02

ztoSnaps = {
    0 :[99,135,28],
    1 :[50,86,19],
    2 :[33,68,15],
    3 :[25,60,12],
    4 :[21,54,10],
    5 :[17,49,8 ],
    6 :[13,45,6 ],
    7 :[11,41,5 ],
    8 :[8 ,38,4 ],
    9 :[6 ,35,3 ],
    10:[4 ,32,2 ]
}
sSFRcut = {
    0 :(-10.5,-9.5),
    1 :(-10,-9),
    2 :(-9.5,-8.5),
    3 :(-9,-8),
    4 :(-10.5,-9.6),
    5 :(-10,-9),
    6 :(-9.5,-8.5),
    7 :(-9,-8),
    8 :(-10.5,-9.6),
    9 :(-10,-9),
    10:(-10,-9),
}

def get_Z_Mstar_SFR(currentDir,which='stars'):
    
    which=which.upper()
    
    h      = 6.774E-01
    xh     = 7.600E-01
    zo     = 3.500E-01
    mh     = 1.6726219E-24
    kb     = 1.3806485E-16
    mc     = 1.270E-02
    Zsun   = 1.27E-02
    
    Zgas      = np.load( currentDir + 'Zgas.npy' )
    Zstar     = np.load( currentDir + 'Zstar.npy' ) 
    star_mass = np.load( currentDir + 'Stellar_Mass.npy'  )
    gas_mass  = np.load( currentDir + 'Gas_Mass.npy' )
    SFR       = np.load( currentDir + 'SFR.npy' )

    sfms_idx = sfmscut(star_mass, SFR)

    desired_mask = ((star_mass > 1.00E+01**(m_star_min)) &
                    (star_mass < 1.00E+01**(m_star_max)) &
                    (gas_mass  > 1.00E+01**(m_gas_min))  &
                    (sfms_idx))

    gas_mass  = gas_mass[desired_mask]
    star_mass = star_mass[desired_mask]
    SFR       = SFR[desired_mask]
    Zstar     = Zstar[desired_mask]
    Zgas      = Zgas[desired_mask]

    Zstar /= Zsun
    OH     = Zgas * (zo/xh) * (1.00/16.00)

    Zgas      = np.log10(OH) + 12

    # Get rid of nans and random values -np.inf
    nonans    = ~(np.isnan(Zgas)) & ~(np.isnan(Zstar)) & (Zstar > 0.0) & (Zgas > 0.0) 

    sSFR      = SFR/star_mass

    gas_mass  = gas_mass[nonans]
    star_mass = star_mass[nonans]
    SFR       = SFR[nonans]
    sSFR      = sSFR[nonans]
    Zstar     = Zstar[nonans]
    Zgas      = Zgas[nonans]

    star_mass     = np.log10(star_mass)
    Zstar         = np.log10(Zstar)
    
    if which == "STARS":
        return Zstar, star_mass, sSFR
    elif which == "GAS":
        return Zgas, star_mass, sSFR
    
    
def fixed_sSFR_bins(mass,Z,sSFR, nBins=4):
    
    start = np.percentile(sSFR,5)
    end   = np.percentile(sSFR,95)
    
    binWidth = (end - start) / nBins
        
    current = start
    
    index = 0
    
    mbin0 = 8.0
    mbin1 = 8.5
    mbin2 = 9.0
    mbin3 = 9.5
    mbin4 = 10.0
    mbin5 = 10.5
    mbin6 = 11.0
    mbins = [mbin0,mbin1,mbin2,mbin3,mbin4,mbin5,mbin6]
    
    fixed_mass = np.ones( (nBins,len(mbins)) )
    fixed_Z    = np.ones( (nBins,len(mbins)) )
    fixed_sSFR = np.ones( nBins )
    while (current < end - binWidth/2):
        within_bin = ( (sSFR > current) & (sSFR < current + binWidth) )
        
        c_mass = mass[within_bin]
        c_Z    =    Z[within_bin]
                
        Z0 = np.nanmedian( c_Z[ (c_mass > mbin0) & (c_mass < mbin1) ] )
        Z1 = np.nanmedian( c_Z[ (c_mass > mbin1) & (c_mass < mbin2) ] )
        Z2 = np.nanmedian( c_Z[ (c_mass > mbin2) & (c_mass < mbin3) ] )
        Z3 = np.nanmedian( c_Z[ (c_mass > mbin3) & (c_mass < mbin4) ] )
        Z4 = np.nanmedian( c_Z[ (c_mass > mbin4) & (c_mass < mbin5) ] )
        Z5 = np.nanmedian( c_Z[ (c_mass > mbin5) & (c_mass < mbin6) ] )
        Z6 = np.nanmedian( c_Z[ (c_mass > mbin6) ] )
        
        fixed_Z[index], fixed_mass[index] = [Z0,Z1,Z2,Z3,Z4,Z5,Z6] , mbins 
        
        fixed_sSFR[index] = current
        current += binWidth
        index += 1
        
    return fixed_mass, fixed_Z, fixed_sSFR

def fixed_M_bins(mass,Z,sSFR,spacing=1.0,nbins=5):
    
    mbins = np.arange( 8.0,12.0,spacing )
    start = 8.0
    end   = 12.0
    
    binWidth = (end - start) / len( mbins )
        
    current = start
    
    n_sSFR_bins = nbins
    sSFR_bins   = np.linspace( np.min(sSFR), np.max(sSFR) , n_sSFR_bins )
        
    fixed_mass = np.ones( len(mbins) )
    fixed_sSFR = np.ones( (len(mbins),len(sSFR_bins)) )
    fixed_Z    = np.ones( (len(mbins),len(sSFR_bins)) )
    for index, current in enumerate(mbins):
        within_bin = ( (mass > current) & (mass < current + binWidth) )
        
        c_sSFR = sSFR[within_bin]
        c_Z    =    Z[within_bin]
        
        Zs = np.ones( n_sSFR_bins )
        
        for i in range(len(sSFR_bins)-1):
            bin0 = sSFR_bins[i]
            bin1 = sSFR_bins[i+1]
            Zs[i] = np.nanmedian( c_Z[ (c_sSFR > bin0) & (c_sSFR < bin1) ] )
        
        Zs[-1] = np.nanmedian( c_Z[ (c_sSFR > sSFR_bins[-1]) ] )
        
        fixed_Z[index], fixed_sSFR[index] = Zs, sSFR_bins
        
        fixed_mass[index] = current
        index += 1
        
    return fixed_mass, fixed_Z, fixed_sSFR

def get_avg_scatter(x,y,bins=5):
    start = np.min(x)
    end   = np.max(x)
    
    bin_centers = np.linspace(start,end,bins)
    binWidth = (end - start) / bins
    
    bin_vals = np.ones(len(bin_centers))
    bin_disp = np.ones(len(bin_centers))
    
    for index, current in enumerate(bin_centers):
        mask = ((x > current - binWidth/2) & (x < current + binWidth/2))
        
        if (sum(mask) < 10):
            bin_vals[index] = np.nan
            bin_disp[index] = np.nan
        else:
            bin_vals[index] = np.nanmedian( y[mask] )
            bin_disp[index] = np.nanstd( y[mask] )
        
    return bin_centers, bin_vals, bin_disp

def medianZR( x, y, bins = 5, return_err_bars = False ):
    
    start = np.min(x)
    end   = np.max(x)
    
    bin_centers = np.linspace(start,end,bins)
    binWidth = (end - start) / bins
    
    bin_vals = np.ones(len(bin_centers))
    bin_errs = np.ones(len(bin_centers))
    
    for index, current in enumerate(bin_centers):
        mask = ((x > current - binWidth/2) & (x < current + binWidth/2))
        
        if (sum(mask) < 10):
            bin_vals[index] = np.nan
            bin_errs[index] = np.nan
        else:
            bin_vals[index] = np.nanmedian( y[mask] )
            bin_errs[index] = np.nanstd( y[mask] )
    
    if return_err_bars:
        return bin_centers, bin_vals, bin_errs
    else:
        return bin_centers, bin_vals
    
def getScatter(x, y, nbins = 10, func=np.median, percLow = 1, percHigh=99):
    
    binWidth = ( np.percentile(x,percHigh) - np.percentile(x,percLow) ) / nbins
    
    start = np.percentile(x,percLow) + binWidth/2
    end   = np.percentile(x,percHigh)
    
    current = start
    
    devs       = np.ones(nbins) * np.nan
    weightDevs = np.ones(nbins) * np.nan
    binLocs    = np.ones(nbins) * np.nan
    medians    = np.ones(nbins) * np.nan
    
    index = 0
    while (current < end):
        
        mask = ( (x < current + binWidth/2) &
                 (x > current - binWidth/2) )
                
        medians[index]    = np.median( y[mask] )
        devs[index]       = np.std( y[mask] )
        weightDevs[index] = np.std( y[mask] * sum(mask) ) 
        binLocs[index]    = current
        
        index += 1
        current += binWidth
    
    return binLocs, medians,devs#func(devs), func(weightDevs)/len(y)

def skew(DATA, ax_histx, ax, WHICH_SIM, snap, snap2z, labelsOff=False):
    
    currentDir = DATA + 'snap%s/' %snap

    Zgas      = np.load( currentDir + 'Zgas.npy' )
    Zstar     = np.load( currentDir + 'Zstar.npy' ) 
    star_mass = np.load( currentDir + 'Stellar_Mass.npy'  )
    gas_mass  = np.load( currentDir + 'Gas_Mass.npy' )
    SFR       = np.load( currentDir + 'SFR.npy' )

    sfms_idx = sfmscut(star_mass, SFR)

    desired_mask = ((star_mass > 1.00E+01**(m_star_min)) &
                    (star_mass < 1.00E+01**(m_star_max)) &
                    (gas_mass  > 1.00E+01**(m_gas_min))  &
                    (sfms_idx))

    gas_mass  = gas_mass[desired_mask]
    star_mass = star_mass[desired_mask]
    SFR       = SFR[desired_mask]
    Zstar     = Zstar[desired_mask]
    Zgas      = Zgas[desired_mask]
    
    Zstar /= Zsun
    OH     = Zgas * (zo/xh) * (1.00/16.00)
    
    Zgas      = np.log10(OH) + 12
    
    # Get rid of nans and random values -np.inf
    nonans    = ~(np.isnan(Zgas)) & ~(np.isnan(Zstar)) & (Zstar > 0.0) & (Zgas > 0.0) 
    
    sSFR      = SFR/star_mass

    sSFR[~(sSFR > 0.0)] = 1e-16

    star_mass = star_mass[nonans]
    sSFR      = sSFR[nonans]
    Zstar     = Zstar[nonans]
    Zgas      = Zgas[nonans]
    
    star_mass     = np.log10(star_mass)
    Zstar         = np.log10(Zstar)
    
    masses, starMedians = getMedians(star_mass,Zstar)
    
    filter_fit_nans_star = ~(np.isnan(starMedians))
    
    masses      = masses[filter_fit_nans_star]
    starMedians = starMedians[filter_fit_nans_star]
    
    MZsR = interp1d(masses,starMedians,fill_value='extrapolate')
    
    starOffset = Zstar - MZsR(star_mass)

    masses, gasMedians = getMedians(star_mass,Zgas)
    
    filter_fit_nans_gas = ~(np.isnan(gasMedians))
    
    masses     = masses[filter_fit_nans_gas]
    gasMedians = gasMedians[filter_fit_nans_gas]
    
    MZgR = interp1d(masses, gasMedians, fill_value='extrapolate')
    
    gasOffset = Zgas - MZgR(star_mass)
    
    nbins = 100

    fig = plt.gcf()

    gaussX = np.linspace(np.min(gasOffset), np.max(gasOffset),100)
    
    gasMu    = np.mean(gasOffset)
    gasSigma = np.std(gasOffset)
    
    mean_guess, std_guess = norm.fit(gasOffset)
    
    gaussY = norm.pdf( gaussX, mean_guess, std_guess )
    
    a_guess, loc_guess, scale_guess = skewnorm.fit(gasOffset)
    
    gaussSkewY = skewnorm.pdf( gaussX, a_guess, loc = loc_guess, scale = scale_guess )
    
    bins, edges, _ = ax_histx.hist( gasOffset, bins=nbins)
    ax_histx.plot( gaussX, gaussY    *np.max(bins)/np.max(gaussY)    , lw=2, color='red' )
    ax_histx.plot( gaussX, gaussSkewY*np.max(bins)/np.max(gaussSkewY), lw=2, color='b')
#         ax_histx.axvline(0,color='k',linestyle='--')
#         ax_histx.axvline(gasMu,color='red',linestyle='--')
    
    ax_histx.text(0.05,0.7,r'${\rm Skewness:}~ %s$' %round(a_guess,2), transform=ax_histx.transAxes, color='b' )
    ax_histx.text(0.05,0.5,r'${\rm Mean:}~ %s$' %round(mean_guess,2) , transform=ax_histx.transAxes, color='r' )
    ax_histx.text(0.05,0.3,r'${\rm Std:}~ %s$' %round(std_guess,2)   , transform=ax_histx.transAxes, color='r' )
    
    plotXmin = -1
    plotXmax = 0.49
    plotYmin = -0.45
    plotYmax = 0.29
    
    ax.set_xlim(plotXmin,plotXmax)
    ax.set_ylim(plotYmin,plotYmax)
    
    Hist1, xedges, yedges = np.histogram2d(gasOffset,starOffset,bins=(nbins,nbins),
                                           range=[[plotXmin, plotXmax],[plotYmin,plotYmax]])
    
    Hist1 = np.transpose(Hist1)
    ax.axhline(0,color='gray',linestyle='--',alpha=0.75)
    ax.axvline(0,color='gray',linestyle='--',alpha=0.75)
    
    ax.pcolormesh(xedges,yedges,np.log10(Hist1),cmap=plt.cm.viridis)

    ax.set_xlabel(r'$\Delta Z_{\rm gas}$')


    x = gasOffset[:,np.newaxis]
    y = starOffset

    a, _, _, _ = np.linalg.lstsq(x,y)

    ax.plot(x,a*x,color='k',lw=4)
        
    if (WHICH_SIM=="ORIGINAL"):
        sim = r'${\rm Illustris}$'
    if (WHICH_SIM=='TNG'):
        sim = r'${\rm TNG}$'
    if (WHICH_SIM == "EAGLE"):
        sim = r'${\rm EAGLE}$'
    
    ax.text(0.05,0.9 , r'%s' %(sim) + ' ' + r'$%s$' %(snap2z[snap]), transform=ax.transAxes)
    ax.text(0.05,0.85, r'${\rm Slope:}~ %s$' %(round(float(a),2)), transform=ax.transAxes)
    
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    one_to_one = np.arange(-2,2,0.01)
    
    ax.plot( one_to_one, one_to_one, color='gray', linestyle='solid', lw=3 )

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    if (labelsOff):
        ax.set_yticklabels([])
    else:
        ax.set_ylabel(r'$\Delta Z_*$')