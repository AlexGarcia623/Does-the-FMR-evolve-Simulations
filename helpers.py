'''
This module defines several useful functions for the analysis
of this work
#
Functions:
    - getMedians(x,y,width=0.1,step=0.05,
                 return_masks=False,percentile=50,
                 min_samp=10)
        Get the medians metallicity within fixed mass bins
        
    - line(data, a, b)
        Defines a line

    - fourth-order(data, a, b, c, d, e)
        Defines a fourth-order polynomial

    - ttest(hypothesized_value,measurements,errors):
        Perform 1 sample t-test

    - estimate_symmetric_error(lower, upper):
        Estimate symmetic error bars given upper and lower
        
    - sfmscut(m0, sfr0, THRESHOLD=-5.00E-01,
              m_star_min=8.0, m_star_max=12.0)
        Compute specific star formation main sequence
#     
Code written by: Alex Garcia, 2023-24
'''
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats


def getMedians(x,y,width=0.1,step=0.05,return_masks=False,percentile=50,min_samp=10):
    '''Get the medians metallicity within fixed mass bins
    
    Inputs:
    - x (ndarray): masses
    - y (ndarray): metallicities
    - width (float): mass bin width
    - step (float): mass bin step size
    - return_masks (bool): Flag to return mass bin masks
    - percentile (int): percentile to return (default is 50, median)
    - min_sample (int): minimum number of galaxies in a bin
    
    Returns:
    - (ndarray): median mass bins
    - (ndarray): corresponding metallicity bins
    - IF return_masks... (ndarray): corresponding mass bin masks
    '''
    start = np.min(x)
    end   = np.max(x)
    
    current = start
    
    medians = []
    xs      = []
    if (return_masks):
        masks = []
    
    while (current < end + 2*step):
        
        mask = ((x > (current)) & (x < (current + width)))
        if (return_masks):
            masks.append( mask )
        
        if (len(y[mask]) > min_samp):
            medians.append( np.percentile( y[mask], percentile ) )
        else:
            medians.append( np.nan )
            
        xs.append( current )
    
        current += step
    
    medians = np.array(medians)
    xs      = np.array(xs)
    
    nonans = ~(np.isnan(medians))
    
    xs      = xs[nonans] 
    medians = medians[nonans]

    if (return_masks):
        masks = np.array(masks)
        masks = masks[nonans]
        masks = list(masks)
        return xs, medians, masks
    else:
        return xs, medians
    
def line(data, a, b):
    '''Creates linear regression
    
    Inputs:
    - data (ndarray): x-axis data
    - a (float): slope
    - b (float): intercept

    Returns:
    - (ndarray): a*data + b
    '''
    return a*data + b

def fourth_order(data, a, b, c, d, e):
    '''Creates fourth-order regression

    Inputs:
    - data (ndarray): x-axis data
    - a, b, c, d, e (float): regression parameters

    Returns:
    - (ndarray): a + b*data + c*data**2 + d*data**3 + e*data**4
    '''
    return a + b*data + c*data**2 + d*data**3 + e*data**4

def ttest(hypothesized_value,measurements,errors):
    '''Perform 1 sample t-test

    Input:
    - hypothesized_value (float): z=0 alpha value to compare against
    - measuements (ndarry): array of alpha values
    - errors (ndarray): uncertainty on alpha values
    '''
    l = 20
    # Calculate weighted mean and standard error
    weighted_mean = np.sum(measurements / errors**2) / np.sum(1 / errors**2)
    weighted_std_error = np.sqrt(1 / np.sum(1 / errors**2))

    # Calculate t-statistic
    t_stat = (weighted_mean - hypothesized_value) / weighted_std_error

    # Degrees of freedom
    degrees_freedom = len(measurements) - 1

    # Calculate p-value (two-tailed)
    p_val = 2 * stats.t.sf(np.abs(t_stat), degrees_freedom)

    print(f"\t{'Weighted Mean':<{l}}: {weighted_mean:0.3f}")
    print(f"\t{'ref val':<{l}}: {hypothesized_value:0.3f}")
    print("\t\tStatistical Test")
    print(f"\t{'T-statistic':<{l}}: {t_stat:0.3f}")
    print(f"\t{'P-value':<{l}}: {p_val:0.3E}")
    print(f"\t{'Reject (0.05 level)':<{l}}: {p_val < 0.05}")
    
def estimate_symmetric_error(lower, upper):
    '''Errors are non symmetric, but not by much.
    I am just estimating them here
    
    Inputs:
    - lower (ndarray): all lower bound uncertainties
    - upper (ndarray): all upper bound uncertainties

    Returns:
    - (ndarray): average of uncertainties
    '''
    return (lower + upper) / 2

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

if __name__ == '__main__':
    print('Hello World')