import pywed as pw
import numpy as np
import matplotlib.pyplot as plt
from signals import *
from tqdm import tqdm

import imas_west

#%%
def get_isotopic_ratio(pulse, channel_name='LODIVOU15'):
    '''
    Get isotopic ratio nH/nD from a given channel name

    Parameters
    ----------
    pulse: int
        WEST pulse number
    channel_name: str
        Spectroscopic channel name. Default is 'LODIVOU15'

    Returns
    -------
    time: np.array
        time array
    isotopic_ratio: np.array
        isotopic ratio nH/nD 
    '''
    density_ratio, time = None, None
    try:
        specv = imas_west.get(pulse, 'spectrometer_visible')
        # find the channel index corresponding to the channel to use for isotopic ratio measurement
        channel_idx = [channel_idx for (channel_idx, channel) in enumerate(specv.channel) if channel_name in channel.name]

        if channel_idx:
            # the specific channel exists, return the data
            channel_idx = channel_idx[0]
            time = specv.channel[channel_idx].isotope_ratios.time - 32
            density_ratio = 100*specv.channel[channel_idx].isotope_ratios.isotope[0].density_ratio
        else:
           time = None
           density_ratio = None
    except Exception as e:
        pass
    return density_ratio, time

#%%
def truncate_using_ip(pulse, t, y, ip_thres = 0.3):
    '''
    truncate y(t) for values of plasma current larger than a threshold
    '''
    ip, ip_t = get_sig(pulse, signals['Ip'])
    y_interpolated = np.interp(ip_t, t, y)
    # filter data
    y_keep = y_interpolated[ip > ip_thres]
    t_keep = ip_t[ip > ip_thres]
    return y_keep, t_keep

def isotopic_ratio_mean_min_max(pulse, ip_thres=0.4):
    ''' 
    return the mean, min and max of isotopic ratio for a given pulse. Measurements are truncated for a given Ip threshold
    '''
    mean, min, max = None, None, None

    density_ratio, time = get_isotopic_ratio(pulse)
    try:
        if np.any(density_ratio):
            nH_over_nD, t_nH_over_nD = truncate_using_ip(pulse, time, density_ratio, ip_thres=ip_thres)
            mean = np.mean(nH_over_nD)
            min = np.amin(nH_over_nD)
            max = np.amax(nH_over_nD)
    except Exception as e:
        pass
    
    return mean, min, max    
     
#%%
pulses = range(56287, 56489)

#%%
means, mins, maxs = [], [], []
for pulse in tqdm(pulses):
   mean, min, max = isotopic_ratio_mean_min_max(pulse)
   means.append(mean)
   mins.append(min)
   maxs.append(max)

_mins = np.array(mins)
_mins[_mins == None] = np.nan
_mins = np.array(_mins, dtype='float')

_maxs = np.array(maxs)
_maxs[_maxs == None] = np.nan
_maxs = np.array(_maxs, dtype='float')

_means = np.array(means)
_means[_means == None] = np.nan
_means = np.array(_means, dtype='float')

np.save('data', [_means, _mins, _maxs])





