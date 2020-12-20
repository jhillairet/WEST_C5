# -*- coding: utf-8 -*-
"""
WEST EUROfusion session 18/12/200

ICRH Antenna Toggling

In these three pulses, the Pr/Pi interlock was removed at generator side.

"""
#%%
from pppat.control_room.signals import *

#%%
sigs = [[signals['IC_P_Q1'], signals['IC_P_Q2'], signals['IC_P_Q4']],
        [signals['IC_Phase_Q1'], signals['IC_Phase_Q2'], signals['IC_Phase_Q4']]]


#%% Q1, Q2 and Q4 
pulses = [56487]
fig, axes = scope(pulses, signames=sigs, cycling_mode='color', alpha=.7)
axes[0].set_xlim(3.8, 5)
fig.suptitle(f'WEST #{pulses}')

#%% Q1, Q2 and Q4
pulses = [56489]
fig, axes = scope(pulses, signames=sigs, cycling_mode='color', alpha=.7)
axes[0].set_xlim(3.8, 5)
fig.suptitle(f'WEST #{pulses}')

#%% Q1 and Q2 only (to check if Q4 perturbs)
pulses = [56488]
sigs = [[signals['IC_P_Q1'], signals['IC_P_Q2']],
        [signals['IC_Phase_Q1'], signals['IC_Phase_Q2']]]

fig, axes = scope(pulses, signames=sigs, cycling_mode='color', alpha=0.7)
axes[0].set_xlim(3.8, 5)
fig.suptitle(f'WEST #{pulses}')
#%%


