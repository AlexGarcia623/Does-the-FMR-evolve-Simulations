'''
This file is used to create Figure 4 of "Does the fundamental 
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
import cmasher as cmr
### Imports From this library
import sys, os
sys.path.append(os.path.dirname(os.getcwd()))
from does_the_fmr_evolve_simulations.getAlpha import (
    whichSim2Tex
)
from does_the_fmr_evolve_simulations.dZdSFR import (
    thin_mass_bin
)

sims = ['original','tng','eagle']

mpl.rcParams['font.size']=18

thin_low=8.0
thin_high=8.5

redshifts = np.arange(0,9)

cmap = cmr.get_sub_cmap('cmr.guppy', 0.0, 1.0, N=len(redshifts))
#mpl.colors.LinearSegmentedColormap.from_list("", ["lightseagreen","gold","lightcoral"])

newcolors = np.linspace(0, 1, len(redshifts))
colors = [ cmap(x) for x in newcolors[::-1] ]

fig, axs = plt.subplots(1,3,figsize=(10,3.5),sharey=True)
for sim_index, sim in enumerate(sims):
    ax = axs[sim_index]
    worms_x, worms_y = thin_mass_bin(sim.upper(),plt.gca(),thin_low=thin_low,thin_high=thin_high)

    for index, worm_x in enumerate(worms_x):
        worm_y = worms_y[index]
        
        a,b = np.polyfit(worm_x,worm_y,1)
        xs = np.linspace(np.min(worm_x),np.max(worm_x),100)
        
        color = colors[index]#'C' + str(index)
        
        ax.plot(worm_x,worm_y,lw=1.5,color=color,alpha=0.5)
        ax.plot(xs, a*xs+b,lw=2,linestyle='--',color=color,label=r'$z=%s$' %index)
    
    ax.text(0.03, 0.05, whichSim2Tex[sim.upper()], transform=ax.transAxes)

    ax.set_xlabel(r'$\log({\rm SFR}~[M_\odot\;{\rm yr}^{-1}] )$')
    ax.axhline(0.0,color='gray',linestyle='--')
        
axs[0].set_ylabel(r'$\Delta Z$')

axs[0].text(0.03,0.15,r'$10^{%.1f} < M_* [M_\odot] < 10^{%.1f}$' %(thin_low,thin_high),
            transform=axs[0].transAxes, fontsize=16)

ymin,ymax=axs[0].get_ylim()
# axs[0].set_ylim(-0.25,0.25)

leg = axs[2].legend(frameon=False,labelspacing=0.05,
                    handletextpad=0, handlelength=0, markerscale=-1,bbox_to_anchor=(1,1))
for i in range(len(leg.get_texts())): leg.legendHandles[i].set_visible(False)
for index, text in enumerate(leg.get_texts()):
    text.set_color(colors[index])
plt.tight_layout()
plt.subplots_adjust(wspace=0.0)

plt.savefig('./Figures (pdfs)/' + 'Figure4.pdf',bbox_inches='tight')
