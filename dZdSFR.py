'''
This module sets defines a useful function for
the creation of Figure 4
#
Functions:
    - thin_mass_bin(sim, ax, m_star_min=8.0, m_star_max=12.0, m_gas_min=8.5,
                  STARS_OR_GAS='gas', polyorder=1, THRESHOLD=-5.00E-01,
                  thin_low=8.0,thin_high=12.0)
        Get the (offset from MZR) metallicity and SFRs in a thin mass bin
#
Code written by: Alex Garcia, 2023-24
'''
### Standard Imports
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
### Imports From this library
import sys, os
sys.path.append(os.path.dirname(os.getcwd()))
from does_the_fmr_evolve_simulations.getAlpha import (
    switch_sim, BLUE, whichSim2Tex
)
from does_the_fmr_evolve_simulations.helpers import (
    getMedians, sfmscut
)

run, base, out_dir, snapshots = None, None, None, []
snap2z = {}

h      = 6.774E-01
xh     = 7.600E-01
zo     = 3.500E-01
mh     = 1.6726219E-24
kb     = 1.3806485E-16
mc     = 1.270E-02
Zsun   = 1.27E-02

def thin_mass_bin(sim, ax, m_star_min=8.0, m_star_max=12.0, m_gas_min=8.5,
                  STARS_OR_GAS='gas', polyorder=1, THRESHOLD=-5.00E-01,
                  thin_low=8.0,thin_high=12.0):
    '''Get the (offset from MZR) metallicity and SFRs in a thin mass bin

    Inputs:
    - sim (String): name of simulation
    - ax (plt axis): axis to plot on
    - m_star_min (float): minimum stellar mass 
    - m_star_max (float): maximum stellar mass
    - m_gas_min (float): minimum gas mass
    - STARS_OR_GAS (String): Get stellar or gas-phase metallicity
    - polyorder (int): Order of fitting polynomial
    - THRESHOLD (float): threshold for sSFMS (see appendix of paper)
    - thin_low (float): minimum of mass bin
    - thin_high (float): maximum of mass bin

    Returns:
    - (ndarray): median metallicity bins
    - (ndarray): median SFR bins
    '''
    STARS_OR_GAS = STARS_OR_GAS.upper()
    
    snapshots, snap2z, BLUE_DIR = switch_sim(sim)
    
    worm_x = []
    worm_y = []
    
    for gbl_index, snap in enumerate(snapshots):

        currentDir = BLUE_DIR + 'snap%s/' %snap

        Zgas      = np.load( currentDir + 'Zgas.npy' )
        Zstar     = np.load( currentDir + 'Zstar.npy' ) 
        star_mass = np.load( currentDir + 'Stellar_Mass.npy'  )
        gas_mass  = np.load( currentDir + 'Gas_Mass.npy' )
        SFR       = np.load( currentDir + 'SFR.npy' )
        R_gas     = np.load( currentDir + 'R_gas.npy' )
        R_star    = np.load( currentDir + 'R_star.npy' )

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
        R_gas     = R_gas    [desired_mask]
        R_star    = R_star   [desired_mask]
        
        Zstar /= Zsun
        OH     = Zgas * (zo/xh) * (1.00/16.00)
        
        Zgas      = np.log10(OH) + 12
        
        # Get rid of nans and random values -np.inf
        nonans    = ~(np.isnan(Zgas)) & ~(np.isnan(Zstar)) & (Zstar > 0.0) & (Zgas > 0.0) 
        
        sSFR      = SFR/star_mass

        gas_mass  = gas_mass [nonans]
        star_mass = star_mass[nonans]
        SFR       = SFR      [nonans]
        sSFR      = sSFR     [nonans]
        Zstar     = Zstar    [nonans]
        Zgas      = Zgas     [nonans]
        R_gas     = R_gas    [nonans]
        R_star    = R_star   [nonans]
        
        star_mass = np.log10(star_mass)
        Zstar     = np.log10(Zstar)
        SFR       = np.log10(SFR)
        
        thin_mass_bin = (star_mass > thin_low) & (star_mass < thin_high)
        
        star_mass = star_mass[thin_mass_bin]
        Zgas = Zgas[thin_mass_bin]
        SFR = SFR[thin_mass_bin]
        
        mass_bins, median_bins = getMedians(star_mass,Zgas,width=0.05,step=0.025)
        
        filter_nans = ~(np.isnan(median_bins))
        
        mass_bins   = mass_bins[filter_nans]
        median_bins = median_bins[filter_nans]
        
        MZR = interp1d(mass_bins,median_bins,fill_value='extrapolate')
        
        metal_offset = Zgas - MZR(star_mass)
        
        SFR_bins, metal_bins = getMedians(SFR,metal_offset,width=0.5,step=0.25)
        filter_nans = ~(np.isnan(metal_bins))
        worm_x.append(SFR_bins[filter_nans])
        worm_y.append(metal_bins[filter_nans])
        
    return worm_x, worm_y

if __name__ == "__main__":
    print('Hello World')
