import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

from getAlpha import get_alpha

sims = ['ORIGINAL','TNG','EAGLE']


m_star_min = 8.0
m_star_max = 12.0
m_gas_min  = 8.5

polyorder=1

EAGLE, EAGLE_lower, EAGLE_upper = get_alpha( 'EAGLE', m_star_min=m_star_min, m_star_max=m_star_max,
                                             polyorder=polyorder )
TNG, TNG_lower, TNG_upper = get_alpha( 'TNG', m_star_min=m_star_min, m_star_max=m_star_max,
                                       polyorder=polyorder )
ORIGINAL, ORIGINAL_lower, ORIGINAL_upper = get_alpha( 'ORIGINAL', m_star_min=m_star_min, m_star_max=m_star_max, 
                                                      polyorder=polyorder )

EAGLE_upper = EAGLE_upper - EAGLE
EAGLE_lower = EAGLE - EAGLE_lower

ORIGINAL_upper = ORIGINAL_upper - ORIGINAL
ORIGINAL_lower = ORIGINAL - ORIGINAL_lower

TNG_upper = TNG_upper - TNG
TNG_lower = TNG - TNG_lower


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

leg  = plt.legend(frameon=False,handletextpad=0, handlelength=0,
                  markerscale=0,loc='upper left',labelspacing=0.05)

leg  = plt.legend(frameon=True,handletextpad=0, handlelength=0,
                  markerscale=0,loc='lower right',labelspacing=0.05)
lCol = ['C1','C2','C0']
for n, text in enumerate( leg.texts ):
    text.set_color( lCol[n] )
    
# get handles
handles, labels = plt.gca().get_legend_handles_labels()
# remove the errorbars
handles = [h[0] for h in handles]
# use them in the legend
leg = plt.legend(frameon=True,handletextpad=0.75, handlelength=0,labelspacing=0.01,
             loc='lower center')
for n, text in enumerate( leg.texts ):
    text.set_color( lCol[n] )

leg.get_frame().set_alpha(1.0)
leg.get_frame().set_edgecolor('white')

plt.xlabel(r'${\rm Redshift}$')
plt.ylabel(r'$\alpha_{\rm min}$')

plt.scatter(0.25,0.33 ,color='k',alpha=0.5,marker='s')
plt.text(0.35,0.3,r'${\rm M10}$',fontsize=14,alpha=0.5)#,transform=plt.gca().transAxes)

plt.scatter(0.25,0.66 ,color='k',alpha=0.5,marker='s')
plt.text(0.35,0.63,r'${\rm AM13}$',fontsize=14,alpha=0.5)#,transform=plt.gca().transAxes)

plt.scatter(0.25,0.55 ,color='k',alpha=0.5,marker='s')
plt.text(0.35,0.52,r'${\rm C20}$',fontsize=14,alpha=0.5)#,transform=plt.gca().transAxes)
ymin, _ = plt.ylim()

plt.ylim(ymin,1.)

plt.tight_layout()

plt.savefig('Figures (pdfs)/'+"Figure2.pdf", bbox_inches='tight')
plt.show()
