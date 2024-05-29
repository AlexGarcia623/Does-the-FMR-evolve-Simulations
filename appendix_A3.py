###########
# Not used in the final paper. Minimal documentation provided
###########
assert(1==0)

import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from getAlpha import get_alpha_handle_mass_bins as get_alpha

mpl.rcParams['font.size'] = 22

sims = ['ORIGINAL','TNG','EAGLE']

m_star_min = 8.0
m_gas_min = 8.5

polyorder=1

fig, axs = plt.subplots(4,1,figsize=(8,16),sharey=True,sharex=True,
                       gridspec_kw={'height_ratios': [1, 1, 1, 0.66]})

z = np.arange(0,9)
ms = 10

offset = [-0.15,-0.1,-0.05,0.05,0.1,0.15]
m_star_mins = np.array([8.0,9.0,10.0,11.0])
m_star_maxs = m_star_mins + 1.0#np.array([8.5,9.0,10.0,12.0])
labels = []
for index, Max in enumerate(m_star_maxs):
    Min = m_star_mins[index]
    labels.append(r'$10^{%.1f} < M_* (M_\odot) < 10^{%.1f}$' %(Min, Max) )
colors = ['C%s' %i for i in range(len(m_star_mins))]
markers = ['o','^','<','>','v','*','s']

illustris_full_pop = [0.23,0.33,0.39,0.45,0.49,0.52,0.53,0.56,0.59]
tng_full_pop = [0.31,0.61,0.60,0.65,0.68,0.70,0.70,0.70,0.70]
eagle_full_pop = [0.74,0.73,0.65,0.59,0.53,0.46,0.44,0.40,0.31]

zs = np.arange(0,9)

axs[0].scatter(zs, illustris_full_pop, color='k', marker='x',
               s=ms*10, label=r'${\rm Full~Population}$')
axs[1].scatter(zs, tng_full_pop, color='k', marker='x',
               s=ms*10, label=r'${\rm Full~Population}$')
axs[2].scatter(zs, eagle_full_pop, color='k', marker='x',
               s=ms*10, label=r'${\rm Full~Population}$')

for index, m_star_max in enumerate(m_star_maxs):
    m_star_min = m_star_mins[index]
    print(f'min: {m_star_min:.1f}, max: {m_star_max:.1f}')
    EAGLE, EAGLE_lower, EAGLE_upper = get_alpha( 'EAGLE', m_star_min=m_star_min, m_star_max=m_star_max,
                                                 m_gas_min=m_gas_min )
    
    TNG, TNG_lower, TNG_upper = get_alpha( 'TNG', m_star_min=m_star_min, m_star_max=m_star_max,
                                           m_gas_min=m_gas_min )
    ORIGINAL, ORIGINAL_lower, ORIGINAL_upper = get_alpha( 'ORIGINAL', m_star_min=m_star_min, m_star_max=m_star_max, 
                                                          m_gas_min=m_gas_min )


    EAGLE_upper = EAGLE_upper - EAGLE
    EAGLE_lower = EAGLE - EAGLE_lower

    ORIGINAL_upper = ORIGINAL_upper - ORIGINAL
    ORIGINAL_lower = ORIGINAL - ORIGINAL_lower

    TNG_upper = TNG_upper - TNG
    TNG_lower = TNG - TNG_lower
    
    axs[0].errorbar( z+offset[index], ORIGINAL, label=labels[index],
                  alpha=0.75, color=colors[index], yerr = [ORIGINAL_lower, ORIGINAL_upper],
                  linestyle='none', marker=markers[index],markersize=ms,capsize=5 )

    axs[1].errorbar( z+offset[index], TNG, label=labels[index],
                  alpha=0.75, color=colors[index], yerr = [TNG_lower, TNG_upper],
                  linestyle='none', marker=markers[index],markersize=ms,capsize=5 )

    axs[2].errorbar( z+offset[index], EAGLE, label=labels[index],
                  alpha=0.75, color=colors[index], yerr = [EAGLE_lower, EAGLE_upper],
                  linestyle='none', marker=markers[index],markersize=ms,capsize=5 )
    

colors = ['k']
colors += ['C%s' %i for i in range(len(m_star_mins))]
handles, labels = axs[1].get_legend_handles_labels()
    
leg  = axs[3].legend(handles, labels, frameon=True,
                     handlelength=0,loc='lower right',
                     labelspacing=0.05)
for n, text in enumerate( leg.texts ):
    text.set_color( colors[n] )

axs[3].axis('off')
axs[3].text(0.5,0.8,r'${\rm Redshift}$',ha='center',transform=axs[3].transAxes)
for i in range(0,9):
    axs[3].text(i, 0.95, i, ha='center')

leg.get_frame().set_alpha(1.0)
leg.get_frame().set_edgecolor('white')

axs[0].text( 0.5,0.9, r'${\rm Illustris}$', transform=axs[0].transAxes, ha='center' )
axs[1].text( 0.5,0.9, r'${\rm TNG}$', transform=axs[1].transAxes, ha='center' )
axs[2].text( 0.5,0.9, r'${\rm EAGLE}$', transform=axs[2].transAxes, ha='center' )

plt.xlabel(r'${\rm Redshift}$')
for ax in axs:
    ax.set_ylabel(r'$\alpha_{\rm min}$')
    
axs[0].set_ylim(-0.05,1.05)

plt.tight_layout()
plt.subplots_adjust(hspace=0.0)

plt.savefig('Figures (pdfs)/'+"FigureA3.pdf", bbox_inches='tight')
plt.show()
