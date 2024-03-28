'''
This module defines several useful functions for the analysis
of this work
#
Functions:
    - getMedians(mass,metals,width=0.1,step=0.05)
        Get the medians metallicity within fixed mass bins
        
    - estimateSkew(data)
        Use scipy's skewnorm.fit
        
    - line(data, p1, p2)
        Defines a line
        
    - sfmscut(m0, sfr0, THRESHOLD=-5.00E-01,
              m_star_min=8.0, m_star_max=12.0)
        Compute specific star formation main sequence
        
Code written by: Alex Garcia, 2023-24
'''
# Standard Imports
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm, skewnorm

m_star_min = 8.0
m_star_max = 12.0
m_gas_min  = 8.5

def getMedians(mass,metals,width=0.1,step=0.05):
    '''Get the medians metallicity within fixed mass bins
    
    Inputs:
    - mass (ndarray): masses
    - metals (ndarray): metallicities
    - width (float): mass bin width
    - step (float): mass bin step size
    
    Returns:
    - (ndarray): median mass bins
    - (ndarray): corresponding metallicity bins
    '''
    start = np.min(mass)
    end   = np.max(mass)
    
    current = start
    
    medians = []
    masses  = []
    
    while (current < end + 2*step):
        mask = ((mass > (current)) & (mass < (current + width)))
        if (len(metals[mask]) > 10):
            medians.append( np.median( metals[mask] ) )
        else:
            medians.append( np.nan )
        masses.append( current )
        current += step
    return np.array(masses),np.array(medians)
    
def estimateSkew(data):
    '''Use scipy's skewnorm.fit
    
    Inputs:
    - data (ndarray): data to test
    
    Returns:
    - (float, 3): skew parameters
    '''
    a_estimate, loc_estimate, scale_estimate = skewnorm.fit(data)
    
    return a_estimate, loc_estimate, scale_estimate
    
def line(data, p1, p2):
    '''Defines a line
    
    Inputs:
    - data (ndarray): x values
    - p1 (float): slope
    - p2 (float): intercept
    
    Returns:
    - (ndarray) y values
    '''
    return p1*data + p2  

def sfmscut(m0, sfr0, THRESHOLD=-5.00E-01,
            m_star_min=8.0, m_star_max=12.0):
    '''Compute specific star formation main sequence
    
    Adapted from Z.S.Hemler+(2021)
    
    Inputs:
    - m0 (ndarray): mass array
    - sfr0 (ndarray): SFR array
    - THRESHOLD (float): value below which galaxies omitted
    - m_star_min (float): minimum stellar mass
    - m_star_max (float): maximum stellar mass
    
    Returns:
    - (ndarray): boolean array of systems that meet criteria
    '''
    nsubs = len(m0)
    idx0  = np.arange(0, nsubs)
    non0  = ((m0   > 0.000E+00) & 
             (sfr0 > 0.000E+00) )
    m     =    m0[non0]
    sfr   =  sfr0[non0]
    idx0  =  idx0[non0]
    ssfr  = np.log10(sfr/m)
    sfr   = np.log10(sfr)
    m     = np.log10(  m)

    idxbs   = np.ones(len(m), dtype = int) * -1
    cnt     = 0
    mbrk    = 1.0200E+01
    mstp    = 2.0000E-01
    mmin    = m_star_min
    mbins   = np.arange(mmin, mbrk + mstp, mstp)
    rdgs    = []
    rdgstds = []


    for i in range(0, len(mbins) - 1):
        idx   = (m > mbins[i]) & (m < mbins[i+1])
        idx0b = idx0[idx]
        mb    =    m[idx]
        ssfrb = ssfr[idx]
        sfrb  =  sfr[idx]
        rdg   = np.median(ssfrb)
        idxb  = (ssfrb - rdg) > THRESHOLD
        lenb  = np.sum(idxb)
        idxbs[cnt:(cnt+lenb)] = idx0b[idxb]
        cnt += lenb
        rdgs.append(rdg)
        rdgstds.append(np.std(ssfrb))

    rdgs       = np.array(rdgs)
    rdgstds    = np.array(rdgstds)
    mcs        = mbins[:-1] + mstp / 2.000E+00
    
    nonans = (~(np.isnan(mcs)) &
              ~(np.isnan(rdgs)) &
              ~(np.isnan(rdgs)))
    parms, cov = curve_fit(line, mcs[nonans], rdgs[nonans], sigma = rdgstds[nonans])
    mmin    = mbrk
    mmax    = m_star_max
    mbins   = np.arange(mmin, mmax + mstp, mstp)
    mcs     = mbins[:-1] + mstp / 2.000E+00
    ssfrlin = line(mcs, parms[0], parms[1])
        
    for i in range(0, len(mbins) - 1):
        idx   = (m > mbins[i]) & (m < mbins[i+1])
        idx0b = idx0[idx]
        mb    =    m[idx]
        ssfrb = ssfr[idx]
        sfrb  =  sfr[idx]
        idxb  = (ssfrb - ssfrlin[i]) > THRESHOLD
        lenb  = np.sum(idxb)
        idxbs[cnt:(cnt+lenb)] = idx0b[idxb]
        cnt += lenb
    idxbs    = idxbs[idxbs > 0]
    sfmsbool = np.zeros(len(m0), dtype = int)
    sfmsbool[idxbs] = 1
    sfmsbool = (sfmsbool == 1)
    return sfmsbool

if __name__ == "__main__":
    print("Hello World!")