'''
This module defines a function for calculating slopes
#
Functions:
    - get_slopes(sim,m_star_min=8.0, m_star_max=12.0,
                 m_gas_min=8.5, THRESHOLD=-5.00E-01)
        Get slope of offset in MZ*R vs MZgR 
        
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

from interplay_gas_stars.helpers import (
    sfmscut, getMedians
)
from interplay_gas_stars.getAlpha import (
    switch_sim, whichSim2Tex
)

h      = 6.774E-01
xh     = 7.600E-01
zo     = 3.500E-01
mh     = 1.6726219E-24
kb     = 1.3806485E-16
mc     = 1.270E-02
Zsun   = 1.27E-02

def get_slopes(sim,m_star_min=8.0, m_star_max=12.0,
               m_gas_min=8.5, THRESHOLD=-5.00E-01):
    '''Get slope of offset in MZ*R vs MZgR 
    
    Inputs
    - sim (String): simulation name
    - m_star_min (float): minimum stellar mass
    - m_star_max (float): maximum stellar mass
    - m_gas_min (float): minimum gas mass
    - THRESHOLD (float): value for sSFMS cut
    
    Return:
    - (ndarray): all slopes for current simulation
    '''
    sim = sim.upper()
    
    snapshots, snap2z, DATA = switch_sim(sim)
    
    slopes = np.zeros(len(snapshots))
    
    for gbl_index, snap in enumerate(snapshots):
        currentDir = DATA + 'snap%s/' %snap

        Zgas      = np.load( currentDir + 'Zgas.npy' )
        Zstar     = np.load( currentDir + 'Zstar.npy' ) 
        star_mass = np.load( currentDir + 'Stellar_Mass.npy'  )
        gas_mass  = np.load( currentDir + 'Gas_Mass.npy' )
        SFR       = np.load( currentDir + 'SFR.npy' )

        sfms_idx = sfmscut(star_mass, SFR, THRESHOLD=THRESHOLD)

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
        Zgas   = np.log10(OH) + 12

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

        x = gasOffset[:,np.newaxis]
        y = starOffset

        slopes[gbl_index], _, _, _ = np.linalg.lstsq(gasOffset[:,np.newaxis],starOffset)
        
    return slopes
        