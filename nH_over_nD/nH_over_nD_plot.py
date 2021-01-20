# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 17:38:47 2020

@author: JH218595
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
from pppat.control_room.signals import *

#%% loading data
data = np.load('data.npy')
means, mins, maxs = data
pulses = range(56287, 56489)

errors = np.array([means - mins, maxs - means])

#%%
with plt.style.context('ggplot'):
    fig, ax = plt.subplots()
    ax.axhspan(5,10, color='green', alpha=0.2, label='Optimal')
    # ax.fill_between(pulses, mins, maxs, alpha=0.4)
    # ax.plot(pulses, means, 'o', markeredgecolor='darkred', 
    # label='Mean values')
    markers, caps, bars = ax.errorbar(pulses, means, yerr=errors, 
                marker='o', markeredgecolor='darkred', ls='',
                label='Average (min/max) - LODIVOU15')
    ax.set_title('$n_H/n_D$ ($I_p>300$kA)')
    ax.set_xlabel('WEST Pulse #')
    ax.set_ylabel('$n_H/n_D$ [%]')
    ax.set_ylim(bottom=0, top=80)
    ax.legend()
    ax.grid(True)
    plt.xticks(rotation=40)
    fig.tight_layout()
    # loop through bars and caps and set the alpha value
    [bar.set_alpha(0.4) for bar in bars]
    ax.patch.set_alpha(0.5)
    # [cap.set_alpha(0.5) for cap in caps]
    fig.savefig('2020-12-20_nH_over_nD_plot.png', dpi=150)