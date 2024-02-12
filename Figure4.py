import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

import illustris_python as il

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
import matplotlib.gridspec as gridspec

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import ks_2samp, iqr

from getAlpha import whichSim2Tex
from dZdSFR import thin_mass_bin

sims = ['original','tng','eagle']

mpl.rcParams['font.size']=18

thin_low=8.0
thin_high=8.5

fig, axs = plt.subplots(1,3,figsize=(10,3.5),sharey=True)
for sim_index, sim in enumerate(sims):
    ax = axs[sim_index]
    worms_x, worms_y = thin_mass_bin(sim.upper(),plt.gca(),
                                     thin_low=thin_low,thin_high=thin_high)

    for index, worm_x in enumerate(worms_x):
        worm_y = worms_y[index]
        
        a,b = np.polyfit(worm_x,worm_y,1)
        xs = np.linspace(np.min(worm_x),np.max(worm_x),100)
        
        color = 'C' + str(index)
        
        ax.plot(worm_x,worm_y,lw=1.5,label=r'$z=%s$' %index,color=color,alpha=0.5)
        ax.plot(xs, a*xs+b,lw=2,linestyle='--',color=color)
    
    ax.text(0.03, 0.05, whichSim2Tex[sim.upper()], transform=ax.transAxes)

    ax.set_xlabel(r'$\log({\rm SFR}~[M_\odot\;{\rm yr}^{-1}] )$')
    ax.axhline(0.0,color='gray',linestyle='--')
        
axs[0].set_ylabel(r'$\Delta Z$')

axs[0].text(0.03,0.15,r'$10^{8.0} < M_* [M_\odot] < 10^{8.5}$', transform=axs[0].transAxes, fontsize=16)

ymin,ymax=axs[0].get_ylim()
# axs[0].set_ylim(-0.25,0.25)

leg = axs[2].legend(frameon=False,bbox_to_anchor=(1,1),
                    handletextpad=0.25, handlelength=0,
                    labelspacing=0.05,fontsize=18)
col = ['C' + str(i) for i in range(0,9)]
for index, text in enumerate(leg.get_texts()):
    text.set_color(col[index])
plt.tight_layout()
plt.subplots_adjust(wspace=0.0)

plt.savefig('./Figures (pdfs)/' + 'dZdSFR.pdf',bbox_inches='tight')