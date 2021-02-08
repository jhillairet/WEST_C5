# -*- coding: utf-8 -*-
"""
WEST session 21/01/2021

ICRH Antenna Toggling

In these three pulses, the Pr/Pi interlock was removed at generator side.

"Successful" pulses: 56819, 56816, 56811, 56808, 56807, 56805

"""
#%%
from pppat.control_room.signals import *
from pppat.libpulse.pulse_settings import PulseSettings
from pppat.libpulse.waveform import get_waveform

#%%
# pulses = range(56462, 56491+1)  # 18/12/2020 W Sources
# pulses = range(56771, 56802+1)  # 20/01/2021 HC
# pulses = range(56861, 56891+1)  # 26/01/2021 HC
# pulses = range(56892, 56908+1)  # 27/01/2021 HC
pulses = range(56909, 56927+1)  # 27/01/2021 Boron Powder 2
# pulses = range(56803, 56819+1)  # 21/01/2021 Weight of divertor vs main chamber W sources

pulses = [56464]


for pulse in pulses:
    #%% making snapshot images
    # PCS waveforms expected
    ps = PulseSettings(pulse)
    # retrieve the default trajectory waveform
    wf_Q1=get_waveform('rts:WEST_PCS/Actuators/Heating/ICRH/power/1/waveform.ref', ps.waveforms)
    wf_Q2=get_waveform('rts:WEST_PCS/Actuators/Heating/ICRH/power/2/waveform.ref', ps.waveforms)
    wf_Q4=get_waveform('rts:WEST_PCS/Actuators/Heating/ICRH/power/3/waveform.ref', ps.waveforms)
    wf_LH1=get_waveform('rts:WEST_PCS/Actuators/Heating/LHCD/power/1/waveform.ref', ps.waveforms)
    wf_LH2=get_waveform('rts:WEST_PCS/Actuators/Heating/LHCD/power/2/waveform.ref', ps.waveforms)
    
    # retrieve actual achieved and zeros if does not exist
    Ip, t_Ip = get_sig(pulse, signals['Ip'])
    if np.all(np.isnan(Ip)):
        Ip = np.zeros(100)
        t_Ip = np.linspace(0, 10, 100)

    P_Q1, t_Q1 = get_sig(pulse, signals['IC_P_Q1'])
    P_Q2, t_Q2 = get_sig(pulse, signals['IC_P_Q2'])
    P_Q4, t_Q4 = get_sig(pulse, signals['IC_P_Q4'])
    if np.all(np.isnan(t_Q1)):
        P_Q1 = np.zeros_like(Ip)
        t_Q1 = t_Ip
    if np.all(np.isnan(t_Q2)):
        P_Q2 = np.zeros_like(Ip)
        t_Q2 = t_Ip
    if np.all(np.isnan(t_Q4)):
        P_Q4 = np.zeros_like(Ip)
        t_Q4 = t_Ip

    IC_P_tot = np.interp(t_Ip, t_Q1, P_Q1 + P_Q2 + P_Q4)
       
    P_LH1, t_LH1 = get_sig(pulse, signals['LH_P_LH1'])
    P_LH2, t_LH2 = get_sig(pulse, signals['LH_P_LH2'])
    if np.all(np.isnan(t_LH1)):
        P_LH1 = np.zeros_like(Ip)
        t_LH1 = t_Ip
    if np.all(np.isnan(t_LH2)):
        P_LH2 = np.zeros_like(Ip)
        t_LH2 = t_Ip

    P_LH1 = np.interp(t_Ip, t_LH1, P_LH1)
    P_LH2 = np.interp(t_Ip, t_LH2, P_LH2)
    LH_P_tot = P_LH1 + P_LH2
    
    # total IC Power request interpolation
    wf_Q1 = np.interp(t_Ip, wf_Q1.times-40+0.02, wf_Q1.values/1e6) 
    wf_Q2 = np.interp(t_Ip, wf_Q2.times-40+0.02, wf_Q2.values/1e6)
    wf_Q4 = np.interp(t_Ip, wf_Q4.times-40+0.02, wf_Q4.values/1e6)
    wf_IC = wf_Q1 + wf_Q2 + wf_Q4
    wf_LH1 = np.interp(t_Ip, wf_LH1.times-40+0.02, wf_LH1.values/1e6)    
    wf_LH2 = np.interp(t_Ip, wf_LH2.times-40+0.02, wf_LH2.values/1e6)
    wf_LH = wf_LH1 + wf_LH2
    
    # date for the filename
    date, time = get_sig(pulse, signals['Datetime'])
    
    #%%
    with plt.style.context('bmh'):
        fig, ax = plt.subplots()
        ax.plot(t_Ip, wf_LH1, label='LH1 Request', color='darkblue', ls='--')
        ax.plot(t_Ip, wf_LH2, label='LH2 Request', color='lightblue', ls='--')
        ax.plot(t_Ip, P_LH1, color='darkblue', label='LH1')
        ax.plot(t_Ip, P_LH2, color='lightblue', label='LH2')
        
        ax.plot(t_Ip, wf_Q1-0.01, label='Q1 Request', color='C1', ls='--')
        ax.plot(t_Q1, P_Q1, color='C1', label='Q1')
        ax.plot(t_Ip, wf_Q2, label='Q2 Request', color='C2', ls='--')
        ax.plot(t_Q2, P_Q2, color='C2', label='Q2')
        ax.plot(t_Ip, wf_Q4+0.01, label='Q4 Request', color='C3', ls='--')
        ax.plot(t_Q4, P_Q4, color='C3', label='Q4')
        
        ax.legend(ncol=5, loc='upper left', fontsize=7)
        ax.set_xlabel('time [s]')
        ax.set_ylabel('Coupled Power [MW]')
        ax.set_title(f'WEST pulse #{pulse}')
        # determine the figure xlim best fit
        if (np.argwhere(wf_IC > 0.2).size > 0) & (np.argwhere(wf_LH > 0.2).size > 0):  # both LH & IC
            idx_start = np.amin([np.argwhere(wf_LH > 0.2)[0], np.argwhere(wf_IC > 0.2)[0]])
            idx_stop = np.amax([np.argwhere(wf_LH > 0.2)[-1], np.argwhere(wf_IC > 0.2)[-1]])
        elif np.argwhere(wf_IC > 0.2).size > 0:  # IC only
            idx_start = np.argwhere(wf_IC > 0.2)[0]
            idx_stop = np.argwhere(wf_IC > 0.2)[-1]
        elif np.argwhere(wf_LH > 0.2).size > 0:  # LH only
            idx_start = np.argwhere(wf_LH > 0.2)[0]
            idx_stop = np.argwhere(wf_LH > 0.2)[-1]
        else:  # no heating, use Ip
            try:
                idx_start = np.argwhere(Ip > 0.2)[0]
                idx_stop = np.argwhere(Ip > 0.2)[-1]
            except IndexError as e:
                idx_start = 1
                idx_stop = 99

        ax.set_xlim(left=t_Ip[idx_start] - 1, right=t_Ip[idx_stop] + 1)
        ax.set_ylim(bottom=0, top=np.amax([wf_LH.max(), wf_IC.max()])+0.8)    
        fig.tight_layout()
        
        fig.savefig(f'snapshots/WEST_{pulse}_{date[0]}-{date[1]:02}-{date[2]:02}_snapshot_RF.png')
