# -*- coding: utf-8 -*-
"""
WEST C5 ICRH and LH Commissioning figures

@author: JH218595
"""
#%%
import matplotlib.pyplot as plt
from scipy.io import loadmat
from pppat import *  # import all WEST signals stuff

#%%  2 x 1 MW Q2 and Q4
pulses = [56301, 56302]
fig, axes = scope(pulses, 
                  [signals['Ip'], signals['nl'], signals['IC_P_tot']],
                  style_label='seaborn-talk')
axes[-1].set_xlim(0,12)
fig.suptitle(f'WEST pulses {pulses}')
fig.tight_layout()
fig.subplots_adjust(hspace=0)
axes[0].legend()
axes[-1].set_ylim(0,2)
[ax.set_facecolor('.95') for ax in axes]

#%% 56302 with details on Q2 and Q4 powers
pulses = [56302]
fig, axes = scope(pulses, 
                  [signals['Ip'], signals['nl'], 
                   [signals['IC_P_Q2'], signals['IC_P_Q4']]],
                  style_label='seaborn-talk', cycling_mode='color')
axes[-1].set_xlim(0,12)
fig.suptitle(f'WEST pulses {pulses}')
fig.tight_layout()
fig.subplots_adjust(hspace=0)
axes[0].legend()
axes[-1].set_ylim(0,2)
[ax.set_facecolor('.95') for ax in axes]

P_tot, t_P_tot = get_sig(pulses[0], signals['IC_P_tot'])
axes[-1].fill_between(t_P_tot, np.squeeze(P_tot), alpha=0.3, label='Total', color='C3')
axes[-1].legend(('Q2', 'Q4', 'Total'))

#%% 56340 : 2x0.5 MW Q2 and Q4, with LH before
pulses = [56340]

fig, axes = scope(pulses, 
                  [signals['Ip'], 
                   signals['nl'], 
                  [signals['LH_P_LH1'], signals['LH_P_LH2'], signals['IC_P_Q2'], signals['IC_P_Q4']], 
                  signals['Cu']],
                  style_label='seaborn-talk', cycling_mode='color')
axes[-1].set_xlim(0,12)
fig.suptitle(f'WEST pulses {pulses}')
fig.tight_layout()
fig.subplots_adjust(hspace=0)
axes[0].legend()
axes[-2].set_ylim(0,2)

[ax.set_facecolor('.95') for ax in axes]

P_tot, t_P_tot = get_sig(pulses[0], signals['RF_P_tot'])
axes[-2].fill_between(t_P_tot, np.squeeze(P_tot), alpha=0.3, label='Total', color='C6')
axes[-2].legend(('LH1', 'LH2', 'Q2', 'Q4', 'Total'))

#%% 56344: 1 MW Q2, 56345: 1 MW Q4, with LH before
pulses = [56344]

fig, axes = scope(pulses, 
                  [signals['Ip'], 
                   signals['nl'], 
                  [signals['LH_P_LH1'], signals['LH_P_LH2'], signals['IC_P_tot']], 
                  signals['Cu']],
                  style_label='seaborn-talk', cycling_mode='color')
axes[-1].set_xlim(0,12)
fig.suptitle(f'WEST pulses {pulses}')
fig.tight_layout()
fig.subplots_adjust(hspace=0)
axes[0].legend()
axes[-2].set_ylim(0,2)
axes[-2].legend(('LH1', 'LH2', 'Q2'))

[ax.set_facecolor('.95') for ax in axes]

#%% 56345: 1 MW Q4, with LH before
pulses = [56345]

fig, axes = scope(pulses, 
                  [signals['Ip'], 
                   signals['nl'], 
                  [signals['LH_P_LH1'], signals['LH_P_LH2'], signals['IC_P_tot']], 
                  signals['Cu']],
                  style_label='seaborn-talk', cycling_mode='color')
axes[-1].set_xlim(0,12)
fig.suptitle(f'WEST pulses {pulses}')
fig.tight_layout()
fig.subplots_adjust(hspace=0)
axes[0].legend()
axes[-2].set_ylim(0,2)
axes[-2].legend(('LH1', 'LH2', 'Q4'))

[ax.set_facecolor('.95') for ax in axes]

#%% 56349 : LH only
pulses = [56349]

fig, axes = scope(pulses, 
                  [signals['Ip'], 
                   signals['nl'], 
                  [signals['LH_P_LH1'], signals['LH_P_LH2']]],
                  style_label='seaborn-talk', cycling_mode='color')
axes[-1].set_xlim(0,12)
fig.suptitle(f'WEST pulses {pulses}')
fig.tight_layout()
fig.subplots_adjust(hspace=0)
axes[0].legend()
axes[-1].set_ylim(0,3.5)
[ax.set_facecolor('.95') for ax in axes]

P_tot, t_P_tot = get_sig(pulses[0], signals['LH_P_tot'])
axes[-1].fill_between(t_P_tot, np.squeeze(P_tot), alpha=0.3, label='Total', color='C4')
axes[-1].legend(('LH1', 'LH2', 'Total'))
