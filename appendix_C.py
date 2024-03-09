import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from getAlpha import get_alpha

mpl.rcParams['font.size'] = 22

sims = ['ORIGINAL','TNG','EAGLE']

m_star_min = 8.0
m_star_max = 12.0

polyorder=1

fig, axs = plt.subplots(3,1,figsize=(8,12),sharey=True,sharex=True)

z = np.arange(0,9)
ms = 10

offset = [-0.075,-0.1,-0.05,0.0,0.05,0.1]
m_gas_mins = [8.0,8.1,8.25,8.3,8.4,8.5]
labels = []
for m in m_gas_mins:
    labels.append(r'$M_{\rm min,\;gas} = 10^{%.2f}$' %m)
colors = ['C%s' %i for i in range(len(m_gas_mins))]
markers = ['o','^','<','>','v','*']
    
for index, m_gas_min in enumerate(m_gas_mins):
    EAGLE, EAGLE_lower, EAGLE_upper = get_alpha( 'EAGLE', m_star_min=m_star_min, m_star_max=m_star_max,
                                                 m_gas_min=m_gas_min )
    print('')
    continue
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

leg  = axs[1].legend(frameon=True, handlelength=0,loc='lower right',labelspacing=0.05)
for n, text in enumerate( leg.texts ):
    text.set_color( colors[n] )

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

plt.savefig('Figures (pdfs)/'+"FigureC1.pdf", bbox_inches='tight')
plt.show()