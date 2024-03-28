'''
This module sets matplotlib rcParams and defines several 
functions useful for calculating alpha values
#
Functions:
    - switch_sim(WHICH_SIM)
        Get constants in order to switch which simulation to analyze
        
    - get_alpha(sim, m_star_min, m_star_max, m_gas_min=8.5, STARS_OR_GAS='stars',
                polyorder=1,THRESHOLD=-5.00E-01)
        Get the projection of minimum scatter
        
    - get_alpha_scatter(sim, m_star_min, m_star_max, m_gas_min=8.5,
                        STARS_OR_GAS='gas',polyorder=1)
        Same as get_alpha(...) but returns scatter values
        
    - plot(sim, name, redshift, m_star_min=8.0, m_star_max=12.0,
           m_gas_min=8.5, STARS_OR_GAS='stars',polyorder=1,
           THRESHOLD=-5.00E-01)
        Used to create Figure 5
        
Code written by: Alex Garcia, 2023-24
'''
# Standard Imports
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from scipy.stats import ks_2samp, iqr
# Import from this library
import sys, os
sys.path.append(os.path.dirname(os.getcwd()))

from interplay_gas_stars.helpers import sfmscut

###### My Custom rcParams ######
fontsize = 20
mpl.rcParams['font.size'] = fontsize
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
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rc('font',**{'family':'sans-serif','serif':['Times New Roman'],'size':15})
mpl.rc('text', usetex=True)
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.size'] = 15
###############################


BLUE = './Data/'
run, base, out_dir, snapshots = None, None, None, []
snap2z = {}

whichSim2Tex = {
    'TNG'     :r'${\rm TNG}$',
    'ORIGINAL':r'${\rm Illustris}$',
    'EAGLE'   :r'${\rm EAGLE}$'
}

h      = 6.774E-01
xh     = 7.600E-01
zo     = 3.500E-01
mh     = 1.6726219E-24
kb     = 1.3806485E-16
mc     = 1.270E-02
Zsun   = 1.27E-02

def switch_sim(WHICH_SIM):
    '''Get constants in order to switch which simulation to analyze
    
    Inputs:
    - WHICH_SIM: any of (eagle, original, TNG, TNG50-1, TNG50-2])
    
    Returns
    - (list): snapshots
    - (dict): snap2z
    - (String): location in Data
    '''
    BLUE_DIR = BLUE + WHICH_SIM + "/"
    if (WHICH_SIM.upper() == "TNG"):
        # TNG
        run       = 'L75n1820TNG'
        base      = '/orange/paul.torrey/IllustrisTNG/Runs/' + run + '/' 
        out_dir   = base 
        snapshots = [99,50,33,25,21,17,13,11,8] # 6,4
        snap2z = {
            99:'z=0',
            50:'z=1',
            33:'z=2',
            25:'z=3',
            21:'z=4',
            17:'z=5',
            13:'z=6',
            11:'z=7',
            8 :'z=8',
            6 :'z=9',
            4 :'z=10',
        }
    elif (WHICH_SIM.upper() == "TNG50-1"):
        # TNG
        run       = 'L35n2160TNG'
        base      = '/orange/paul.torrey/IllustrisTNG/Runs/' + run + '/' 
        out_dir   = base 
        snapshots = [99,50,33,25,21,17,13,11,8] # 6,4
        snap2z = {
            99:'z=0',
            50:'z=1',
            33:'z=2',
            25:'z=3',
            21:'z=4',
            17:'z=5',
            13:'z=6',
            11:'z=7',
            8 :'z=8',
            6 :'z=9',
            4 :'z=10',
        }
    elif (WHICH_SIM.upper() == "TNG50-2"):
        # TNG
        run       = 'L35n1080TNG'
        base      = '/orange/paul.torrey/IllustrisTNG/Runs/' + run + '/' 
        out_dir   = base 
        snapshots = [99,50,33,25,21,17,13,11,8] # 6,4
        snap2z = {
            99:'z=0',
            50:'z=1',
            33:'z=2',
            25:'z=3',
            21:'z=4',
            17:'z=5',
            13:'z=6',
            11:'z=7',
            8 :'z=8',
            6 :'z=9',
            4 :'z=10',
        }
    elif (WHICH_SIM.upper() == "ORIGINAL"):
        # Illustris
        run       = 'L75n1820FP'
        base      = '/orange/paul.torrey/Illustris/Runs/' + run + '/'
        out_dir   = base
        snapshots = [135,86,68,60,54,49,45,41,38] # 35,32
        snap2z = {
            135:'z=0',
            86 :'z=1',
            68 :'z=2',
            60 :'z=3',
            54 :'z=4',
            49 :'z=5',
            45 :'z=6',
            41 :'z=7',
            38 :'z=8',
            35 :'z=9',
            32 :'z=10',
        }
    elif (WHICH_SIM.upper() == "EAGLE"):
        EAGLE_DATA = BLUE_DIR + 'data/'
        snapshots = [28,19,15,12,10,8,6,5,4] # 3,2
        snap2z = {
            28:'z=0',
            19:'z=1',
            15:'z=2',
            12:'z=3',
            10:'z=4',
             8:'z=5',
             6:'z=6',
             5:'z=7',
             4:'z=8',
             3:'z=9',
             2:'z=10'
        }
    return snapshots, snap2z, BLUE_DIR

def get_alpha(sim, m_star_min, m_star_max, m_gas_min=8.5, STARS_OR_GAS='stars',
              polyorder=1,THRESHOLD=-5.00E-01):
    '''Get the projection of minimum scatter
    
    Inputs:
    - sim (String): name of simulation (eagle, original, tng)
    - m_star_min (float): minimum stellar mass 
    - m_star_max (float): maximum stellar mass
    - m_gas_min (float): minimum gas mass
    - STARS_OR_GAS (String): Get stellar or gas-phase metallicity
    - polyorder (int): Order of fitting polynomial
    - THRESHOLD (float): threshold for sSFMS (see appendix of paper)
    
    Returns:
    - (ndarray): all alpha values for the simulation
    - (ndarray): all lower estimates for alpha
    - (ndarray): all upper estimates for alpha
    '''
    STARS_OR_GAS = STARS_OR_GAS.upper()
    
    snapshots, snap2z, BLUE_DIR = switch_sim(sim)
    
    min_alphas = np.zeros(len(snapshots))
    low_errbar = np.zeros(len(snapshots))
    hi_errbar  = np.zeros(len(snapshots))
    
    for gbl_index, snap in enumerate(snapshots):

        currentDir = BLUE_DIR  + 'snap%s/' %snap

        Zgas      = np.load( currentDir + 'Zgas.npy' )
        Zstar     = np.load( currentDir + 'Zstar.npy' ) 
        star_mass = np.load( currentDir + 'Stellar_Mass.npy'  )
        gas_mass  = np.load( currentDir + 'Gas_Mass.npy' )
        SFR       = np.load( currentDir + 'SFR.npy' )

        # Nominal threshold = -5.000E-01
        sfms_idx = sfmscut(star_mass, SFR, THRESHOLD=THRESHOLD,
                           m_star_min=m_star_min, m_star_max=m_star_max)

        desired_mask = ((star_mass > 1.00E+01**(m_star_min)) &
                        (star_mass < 1.00E+01**(m_star_max)) &
                        (gas_mass  > 1.00E+01**(m_gas_min))  &
                        (sfms_idx))

        gas_mass  = gas_mass [desired_mask]
        star_mass = star_mass[desired_mask]
        SFR       = SFR      [desired_mask]
        Zstar     = Zstar    [desired_mask]
        Zgas      = Zgas     [desired_mask]
        
        Zstar /= Zsun
        OH     = Zgas * (zo/xh) * (1.00/16.00)
        Zgas   = np.log10(OH) + 12
        
        nonans    = ~(np.isnan(Zgas)) & ~(np.isnan(Zstar)) & (Zstar > 0.0) & (Zgas > 0.0) 
        
        gas_mass  = gas_mass [nonans]
        star_mass = star_mass[nonans]
        SFR       = SFR      [nonans]
        Zstar     = Zstar    [nonans]
        Zgas      = Zgas     [nonans]
        
        star_mass = np.log10(star_mass)
        Zstar     = np.log10(Zstar)
                
        if (STARS_OR_GAS == "GAS"):
            Z_use = Zgas
        elif (STARS_OR_GAS == "STARS"):
            Z_use = Zstar
        else:
            break
            
        alphas = np.linspace(0,1,100)        
        disps = np.zeros(len(alphas))
        
        for index, alpha in enumerate(alphas):
            mu_fit  = star_mass - alpha*np.log10( SFR )

            popt = np.polyfit(mu_fit, Z_use, polyorder)
            interp = np.polyval( popt, mu_fit )

            disps[index] = np.std( np.abs(Z_use) - np.abs(interp) ) 

        argmin = np.argmin(disps)
        min_alpha = alphas[argmin]
        min_disp  = disps[argmin]

        width = min_disp * 1.01

        within_uncertainty = alphas[ (disps < width) ]

        min_uncertain = within_uncertainty[0]
        max_uncertain = within_uncertainty[-1] 
        
        min_alphas[gbl_index] = min_alpha
        low_errbar[gbl_index] = min_uncertain
        hi_errbar [gbl_index] = max_uncertain
            
    return min_alphas, low_errbar, hi_errbar
    
def get_alpha_scatter(sim, m_star_min, m_star_max, m_gas_min=8.5, STARS_OR_GAS='gas',
                      polyorder=1):
    '''Same as get_alpha(...) but returns scatter values
    
    Inputs:
    - sim (String): name of simulation (eagle, original, tng)
    - m_star_min (float): minimum stellar mass 
    - m_star_max (float): maximum stellar mass
    - m_gas_min (float): minimum gas mass
    - STARS_OR_GAS (String): Get stellar or gas-phase metallicity
    - polyorder (int): Order of fitting polynomial
    
    Returns:
    - (ndarray): all alpha values for the simulation
    - (ndarray): all lower estimates for alpha
    - (ndarray): all upper estimates for alpha
    '''
    STARS_OR_GAS = STARS_OR_GAS.upper()
    
    snapshots, snap2z, BLUE_DIR = switch_sim(sim)
    
    scatter_global = np.zeros( len(snapshots) )
    scatter_local  = np.zeros( len(snapshots) )
    scatter_MZR    = np.zeros( len(snapshots) )
    
    z0_alpha = None
    
    for snap_index, snap in enumerate(snapshots):        
        currentDir = BLUE_DIR + 'data/' + 'snap%s/' %snap

        Zgas      = np.load( currentDir + 'Zgas.npy' )
        Zstar     = np.load( currentDir + 'Zstar.npy' ) 
        star_mass = np.load( currentDir + 'Stellar_Mass.npy'  )
        gas_mass  = np.load( currentDir + 'Gas_Mass.npy' )
        SFR       = np.load( currentDir + 'SFR.npy' )
        R_gas     = np.load( currentDir + 'R_gas.npy' )
        R_star    = np.load( currentDir + 'R_star.npy' )

        THRESHOLD = -5.00E-01 # nominal = -5.00E-01
        sfms_idx = sfmscut(star_mass, SFR, THRESHOLD=THRESHOLD,
                       m_star_min=m_star_min, m_star_max=m_star_max)

        desired_mask = ((star_mass > 1.00E+01**(m_star_min)) &
                        (star_mass < 1.00E+01**(m_star_max)) &
                        (gas_mass  > 1.00E+01**(m_gas_min))  &
                        (sfms_idx))

        gas_mass  = gas_mass [desired_mask]
        star_mass = star_mass[desired_mask]
        SFR       = SFR      [desired_mask]
        Zstar     = Zstar    [desired_mask]
        Zgas      = Zgas     [desired_mask]
        R_gas     = R_gas    [desired_mask]
        R_star    = R_star   [desired_mask]

        Zstar /= Zsun
        OH     = Zgas * (zo/xh) * (1.00/16.00)

        Zgas      = np.log10(OH) + 12

        # Get rid of nans and random values -np.inf
        nonans    = ~(np.isnan(Zgas)) & ~(np.isnan(Zstar)) & (Zstar > 0.0) & (Zgas > 0.0) 

        sSFR      = SFR/star_mass

        sSFR[~(sSFR > 0.0)] = 1e-16

        star_mass = star_mass[nonans]
        sSFR      = sSFR     [nonans]
        Zstar     = Zstar    [nonans]
        Zgas      = Zgas     [nonans]
        R_gas     = R_gas    [nonans]
        R_star    = R_star   [nonans]

        gas_mass      = np.log10(gas_mass)
        star_mass     = np.log10(star_mass)
        Zstar         = np.log10(Zstar)

        if (STARS_OR_GAS == "GAS"):
            Z_use = Zgas
        elif (STARS_OR_GAS == "STARS"):
            Z_use = Zstar

        alphas = np.linspace(0,1,100)
        disp   = np.zeros( len(alphas) )

        for index, alpha in enumerate(alphas):

            muCurrent = star_mass - alpha*np.log10(SFR) 

            popt   = np.polyfit( muCurrent, Z_use, polyorder )            
            interp = np.polyval( popt, muCurrent )

            disp[index] = np.std( np.abs(Z_use) - np.abs(interp) )
        
        if (snap_index == 0):
            # If z=0, save alpha for comparisons
            z0_alpha = alphas[ np.argmin(disp) ]
        
        scatter_local[snap_index]  = np.min(disp)
        ## [Get the location of the z0 alpha]
        scatter_global[snap_index] = disp[np.argmin(np.abs(alphas - z0_alpha ))]
        ## MZR Scatter = alpha=0.0
        scatter_MZR[snap_index]    = disp[ 0 ]
    
    return scatter_local, scatter_global, scatter_MZR   
    
def plot(sim, name, redshift, m_star_min=8.0, m_star_max=12.0,
         m_gas_min=8.5, STARS_OR_GAS='stars',polyorder=1,
         THRESHOLD=-5.00E-01):
    '''Used to create Figure 5
    
    Inputs:
    - sim (String): name of simulation (eagle, original, tng)
    - name (String): name and directory of output file
    
    Optional Inputs:
    - m_star_min (float): minimum stellar mass = 8.0
    - m_star_max (float): maximum stellar mass = 12.0
    - m_gas_min (float): minimum gas mass = 8.5
    - STARS_OR_GAS (String): Get stellar or gas-phase metallicity = "stars"
    - polyorder (int): Order of fitting polynomial = 1
    - THRESHOLD (float): threshold for sSFMS (see appendix of paper)
    
    Returns:
    - (none): saves figure
    '''
    STARS_OR_GAS = STARS_OR_GAS.upper()
    
    snapshots, snap2z, BLUE_DIR = switch_sim(sim)
    
    if redshift < 0 or redshift > 8:
        return
    
    snap = snapshots[redshift]

    currentDir = BLUE_DIR + 'snap%s/' %snap

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

    gas_mass  = gas_mass [desired_mask]
    star_mass = star_mass[desired_mask]
    SFR       = SFR      [desired_mask]
    Zstar     = Zstar    [desired_mask]
    Zgas      = Zgas     [desired_mask]

    Zstar /= Zsun
    OH     = Zgas * (zo/xh) * (1.00/16.00)

    Zgas      = np.log10(OH) + 12

    nonans    = ~(np.isnan(Zgas)) & ~(np.isnan(Zstar)) & (Zstar > 0.0) & (Zgas > 0.0) 

    sSFR      = SFR/star_mass

    gas_mass  = gas_mass [nonans]
    star_mass = star_mass[nonans]
    SFR       = SFR      [nonans]
    sSFR      = sSFR     [nonans]
    Zstar     = Zstar    [nonans]
    Zgas      = Zgas     [nonans]

    star_mass = np.log10(star_mass)
    Zstar     = np.log10(Zstar)

    alphas = np.linspace(0,1,100)        
    disps = np.ones(len(alphas)) * np.nan
    slopes = np.ones(len(alphas)) * np.nan
    intercepts = np.ones(len(alphas)) * np.nan

    if (STARS_OR_GAS == "GAS"):
        Z_use = Zgas
    elif (STARS_OR_GAS == "STARS"):
        Z_use = Zstar

    for index, alpha in enumerate(alphas):

        muCurrent  = star_mass - alpha*np.log10( SFR )

        mu_fit = muCurrent
        Z_fit  = Z_use

        popt = np.polyfit(mu_fit, Z_fit, polyorder)
        slopes[index], intercepts[index] = popt
        interp = np.polyval( popt, mu_fit )

        disps[index] = np.std( np.abs(Z_fit) - np.abs(interp) ) 

    argmin = np.argmin(disps)
    min_alpha = round(alphas[argmin],2)
    min_disp  = disps[argmin]

    print('%s: %s, best alpha_%s = %s' %(sim, snap2z[snap], STARS_OR_GAS.lower(), min_alpha))
    
    fig, axs = plt.subplots(2,1,figsize=(5,8))
    
    axs[0].plot( alphas, disps, lw=2 )
    axs[0].axhline( min_disp, color='r', linestyle='--' )
    axs[0].axvline( min_alpha, color='r',linestyle='--' )
    axs[0].set_xlim(0,1)
    axs[0].set_ylabel(r'${\rm Dispersion}$')
    axs[0].set_xlabel(r'$\alpha$')


    WHICH_SIM_TEX = {
        "TNG":r"${\rm TNG}$",
        "ORIGINAL":r"${\rm Illustris}$",
        "EAGLE":r"${\rm EAGLE}$"
    }

    axs[0].text( 0.4, 0.6 , r'%s' %WHICH_SIM_TEX[sim], transform=axs[0].transAxes )
    axs[0].text( 0.4, 0.5 , r'$%s$' %snap2z[snap], transform=axs[0].transAxes )

    best_mu = star_mass - min_alpha * np.log10(SFR)
    best_line = np.polyval( [slopes[argmin],intercepts[argmin]], best_mu )
    
    axs[1].hist2d( best_mu, Z_use, bins=(100,100), norm=LogNorm() )
    axs[1].plot( best_mu, best_line, color='k', lw=2 )
    axs[1].set_ylabel(r'$\log(Z_* [Z_\odot])$')
    axs[1].set_xlabel(r'$\mu_{%s} = \log M_* - %s\log{\rm SFR}$' %(min_alpha,min_alpha))

    axs[1].text( 0.375 , 0.2  , r'${\rm Minimum~dispersion~fit}$', 
                transform=axs[1].transAxes, fontsize=14 )
    axs[1].text( 0.55  , 0.125, r'$\alpha=%s$' %min_alpha        , 
                transform=axs[1].transAxes, fontsize=14 )

    plt.tight_layout()
    plt.savefig(name + '.pdf',bbox_inches='tight')
    
if __name__ == "__main__":
    print('Hello World')