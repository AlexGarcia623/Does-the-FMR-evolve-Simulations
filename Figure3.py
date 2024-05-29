'''
This file is used to create Figure 3 of "Does the fundamental 
metallicity relation evolve with redshift? I: the correlation
between offsets from the mass-metallicity relation and star 
formation rate"
#
Paper: https://academic.oup.com/mnras/article/531/1/1398/7671150
#
Code written by: Alex Garcia, 2023-24
'''
### Standard Imports
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
### Imports From this library
import sys, os
sys.path.append(os.path.dirname(os.getcwd()))
from does_the_fmr_evolve_simulations.getAlpha import (
    get_alpha_scatter, whichSim2Tex
)

mpl.rcParams['font.size'] = 18

### Fiducial Hyperparameters
m_star_min = 8.0
m_star_max = 12.0
m_gas_min  = 8.5
polyorder  = 1

### Make figure
fig, axs = plt.subplots(1,3,figsize=(9,3.5),sharey=True,sharex=True)

sims = ['original','tng','eagle']
col  = ['C1','C2','C0']
mark = ['^','*','o']
linestyles = ['solid','solid','solid']

all_weak, all_strong, all_MZR = [],[],[]

for index, sim in enumerate(sims):
    color = col[index]
    symb  = mark[index]
    sim   = sim.upper()
    
    scatter_weak, scatter_strong, scatter_MZR = get_alpha_scatter(sim, m_star_min, m_star_max, m_gas_min=m_gas_min,
                                                                  STARS_OR_GAS='gas',polyorder=1)
    all_weak.append( scatter_weak )
    all_strong.append( scatter_strong )
    all_MZR.append( scatter_MZR )
    
axs[0].axhline(1, color='gray', linestyle=':', alpha=0.5)

for ax in axs:
    ax.set_xlabel( r'${\rm Redshift}$' )
    ax.axhline(1, color='k', linestyle='solid', alpha=1, lw=3)

redshifts = np.arange(0,9)

weak  = np.array( all_weak )
strong = np.array( all_strong )
MZR  = np.array( all_MZR )

ratios1 = weak   / MZR
ratios2 = strong / MZR
ratios3 = weak   / strong

illustris = ratios1[0,:]
tng       = ratios1[1,:]
eagle     = ratios1[2,:]


ms = 7
lw = 1.5
axs[0].plot( redshifts, illustris, color=col[0],
          marker=mark[0], label=whichSim2Tex['ORIGINAL'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[1], lw=lw
)
axs[0].plot( redshifts, tng      , color=col[1],
          marker=mark[1], label=whichSim2Tex['TNG'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[1], lw=lw
)
axs[0].plot( redshifts, eagle    , color=col[2],
          marker=mark[2], label=whichSim2Tex['EAGLE'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[1], lw=lw
)
    
illustris = ratios2[0,:]
tng       = ratios2[1,:]
eagle     = ratios2[2,:]

axs[1].plot( redshifts, illustris, color=col[0],
          marker=mark[0], label=whichSim2Tex['ORIGINAL'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[2], lw=lw
)
axs[1].plot( redshifts, tng      , color=col[1],
          marker=mark[1], label=whichSim2Tex['TNG'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[2], lw=lw
)
axs[1].plot( redshifts, eagle    , color=col[2],
          marker=mark[2], label=whichSim2Tex['EAGLE'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[2], lw=lw
)

illustris = ratios3[0,:]
tng       = ratios3[1,:]
eagle     = ratios3[2,:]
    
axs[2].plot( redshifts, illustris, color=col[0],
          marker=mark[0], label=whichSim2Tex['ORIGINAL'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[2], lw=lw
)
axs[2].plot( redshifts, tng      , color=col[1],
          marker=mark[1], label=whichSim2Tex['TNG'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[2], lw=lw
)
axs[2].plot( redshifts, eagle    , color=col[2],
          marker=mark[2], label=whichSim2Tex['EAGLE'], 
          alpha=0.75, markersize=ms, linestyle=linestyles[2], lw=lw
)
    
fs = 15
leg = axs[2].legend( frameon=False,handletextpad=0.25, handlelength=0,
                     labelspacing=0.05, loc='upper right', fontsize=fs )
for index, text in enumerate(leg.get_texts()):
    text.set_color(col[index])

xmin, xmax = axs[1].get_xlim()
ymin, ymax = axs[1].get_ylim()

axs[0].set_ylabel( r'${\rm Ratio}$' )

axs[0].text( 0.05, 0.85,r'$\sigma_{\rm weak} / \sigma_{\rm MZR}$',
            transform=axs[0].transAxes, fontsize=fs+2)
axs[1].text( 0.05, 0.85,r'$\sigma_{\rm strong} / \sigma_{\rm MZR}$',
        transform=axs[1].transAxes, fontsize=fs+2)
axs[2].text( 0.05, 0.85, r'$\sigma_{\rm weak} / \sigma_{\rm strong}$',
            transform=axs[2].transAxes, fontsize=fs+2)

xmin, xmax = axs[1].get_xlim()
ymin, ymax = axs[1].get_ylim()

error_bar = 0.05 # Nominal 1% error bars

axs[0].fill_between( np.arange(-1,11), (1 + error_bar), (1 - error_bar), color='gray', alpha=0.5 )
axs[1].fill_between( np.arange(-1,11), (1 + error_bar), (1 - error_bar), color='gray', alpha=0.5 )
axs[2].fill_between( np.arange(-1,11), (1 + error_bar), (1 - error_bar), color='gray', alpha=0.5 )

error_bar = 0.10 # Nominal 1% error bars

axs[0].fill_between( np.arange(-1,11), (1 + error_bar), (1 - error_bar), color='gray', alpha=0.25 )
axs[1].fill_between( np.arange(-1,11), (1 + error_bar), (1 - error_bar), color='gray', alpha=0.25 )
axs[2].fill_between( np.arange(-1,11), (1 + error_bar), (1 - error_bar), color='gray', alpha=0.25 )

axs[0].set_xlim(xmin, xmax)
axs[1].set_xlim(xmin, xmax)
axs[0].set_ylim(ymin, ymax)

axs[0].set_yticks([0.6,0.8,1.0,1.2,1.4])

axs[0].set_xticks([0,2,4,6,8])

plt.tight_layout()
plt.subplots_adjust(wspace=0.0)

plt.savefig('Figures (pdfs)/'+'Figure3' + '.pdf', bbox_inches='tight' )