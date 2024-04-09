import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from getAlpha import get_alpha

mpl.rcParams['font.size'] = 22

sims = ['ORIGINAL','TNG','EAGLE']

m_star_min = 8.0
m_star_max = 12.0
m_gas_min  = 8.5

polyorder=1

nominal_sSFMS = -5.00E-01
low_sSFMS = -1.00E-01
hi_sSFMS = -1.00E00
no_sSFMS = -np.inf

thres=[nominal_sSFMS, low_sSFMS, hi_sSFMS, no_sSFMS]

fig, axs = plt.subplots(3,1,figsize=(8,12),sharey=True,sharex=True)

z = np.arange(0,9)
ms = 10

colors = ['C0','C1','C2','C3']
markers = ['^','o','s','*']
labels = [r'$0.1~{\rm dex}$',r'$0.5~{\rm dex~(fiducial)}$',r'$1.0~{\rm dex}$',r'${\rm SFR}>0$']
offset = [-0.05,-0.025,0.025,0.05]

for index, thres_val in enumerate(thres):
    EAGLE, EAGLE_lower, EAGLE_upper = get_alpha( 'EAGLE', m_star_min=m_star_min, m_star_max=m_star_max,
                                                 polyorder=polyorder, THRESHOLD=thres_val )
    TNG, TNG_lower, TNG_upper = get_alpha( 'TNG', m_star_min=m_star_min, m_star_max=m_star_max,
                                           polyorder=polyorder, THRESHOLD=thres_val )
    ORIGINAL, ORIGINAL_lower, ORIGINAL_upper = get_alpha( 'ORIGINAL', m_star_min=m_star_min, m_star_max=m_star_max, 
                                                          polyorder=polyorder, THRESHOLD=thres_val )

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

# leg.get_frame().set_alpha(1.0)
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

plt.savefig('Figures (pdfs)/' + "FigureA1.pdf", bbox_inches='tight')
plt.show()