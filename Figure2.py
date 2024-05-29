'''
This file is used to create Figure 2 of "Does the fundamental 
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
    get_alpha
)
from does_the_fmr_evolve_simulations.helpers import (
    ttest, estimate_symmetric_error
)

sims = ['ORIGINAL','TNG','EAGLE']

### Fiducial Hyperparameters
m_star_min = 8.0
m_star_max = 12.0
m_gas_min  = 8.5
polyorder  = 1

### Get alpha values
EAGLE, EAGLE_lower, EAGLE_upper = get_alpha( 'EAGLE',
                                             m_star_min=m_star_min,
                                             m_star_max=m_star_max,
                                             polyorder=polyorder )
TNG, TNG_lower, TNG_upper = get_alpha( 'TNG',
                                       m_star_min=m_star_min,
                                       m_star_max=m_star_max,
                                       polyorder=polyorder )
ORIGINAL, ORIGINAL_lower, ORIGINAL_upper = get_alpha( 'ORIGINAL',
                                                      m_star_min=m_star_min,
                                                      m_star_max=m_star_max, 
                                                      polyorder=polyorder )


EAGLE_upper = EAGLE_upper - EAGLE
EAGLE_lower = EAGLE - EAGLE_lower

ORIGINAL_upper = ORIGINAL_upper - ORIGINAL
ORIGINAL_lower = ORIGINAL - ORIGINAL_lower

TNG_upper = TNG_upper - TNG
TNG_lower = TNG - TNG_lower

### Make figure
fig = plt.figure(figsize=(9,4))

z = np.arange(0,9)

ms = 7
plt.errorbar( z+0.00, ORIGINAL, label=r'${\rm Illustris}$',
              alpha=0.75, color='C1', yerr = [ORIGINAL_lower, ORIGINAL_upper],
              linestyle='none', marker='^',markersize=ms)

plt.errorbar( z+0.05, TNG, label=r'${\rm TNG}$',
              alpha=0.75, color='C2', yerr = [TNG_lower, TNG_upper],
              linestyle='none', marker='*',markersize=ms )

plt.errorbar( z-0.05, EAGLE, label=r'${\rm EAGLE}$',
              alpha=0.75, color='C0', yerr = [EAGLE_lower, EAGLE_upper],
              linestyle='none', marker='o',markersize=ms )

## Legend
leg  = plt.legend(frameon=True,handletextpad=0, handlelength=0,
                  markerscale=0,loc='lower right',labelspacing=0.05)
lCol = ['C1','C2','C0']
for n, text in enumerate( leg.texts ):
    text.set_color( lCol[n] )
handles, labels = plt.gca().get_legend_handles_labels()
handles = [h[0] for h in handles]
leg = plt.legend(frameon=True,handletextpad=0.75, handlelength=0,labelspacing=0.01,
             loc='lower center')
for n, text in enumerate( leg.texts ):
    text.set_color( lCol[n] )

leg.get_frame().set_alpha(1.0)
leg.get_frame().set_edgecolor('white')

plt.xlabel(r'${\rm Redshift}$')
plt.ylabel(r'$\alpha_{\rm min}$')

plt.scatter(0.25,0.33 ,color='k',alpha=0.5,marker='s')
plt.text(0.35,0.3,r'${\rm M10}$',fontsize=14,alpha=0.5)

plt.scatter(0.25,0.66 ,color='k',alpha=0.5,marker='s')
plt.text(0.35,0.63,r'${\rm AM13}$',fontsize=14,alpha=0.5)

plt.scatter(0.25,0.55 ,color='k',alpha=0.5,marker='s')
plt.text(0.35,0.52,r'${\rm C20}$',fontsize=14,alpha=0.5)
ymin, _ = plt.ylim()

plt.ylim(ymin,1.)

plt.tight_layout()

plt.savefig('Figures (pdfs)/'+"Figure2.pdf", bbox_inches='tight')
plt.show()

### Do reference value t-test
all_alpha = [EAGLE, TNG, ORIGINAL]
all_lower = [EAGLE_lower, TNG_lower, ORIGINAL_lower]
all_upper = [EAGLE_upper, TNG_upper, ORIGINAL_upper]
    
for index,alphas in enumerate(all_alpha):
    lower = all_lower[index]
    upper = all_upper[index]
    
    which_redshift_compare = 0
    hypothesized_value = alphas[which_redshift_compare]
    
    est_error = upper#estimate_symmetric_error( lower, upper )

    print('\n')
    print(f'{sims[index]}, compared to z={which_redshift_compare} alpha value')
    ttest(hypothesized_value, alphas, est_error)
