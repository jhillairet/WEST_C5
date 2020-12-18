import pywed as pw
import numpy as np
import matplotlib.pyplot as plt
from signals import *

import imas_west

CHANNEL_NAME = 'LODIVOU15'

pulse = 56485

specv = imas_west.get(pulse, 'spectrometer_visible')
# find the channel index corresponding to the channel to use for isotopic ratio measurement
channel_idx = [channel_idx for (channel_idx, channel) in enumerate(specv.channel) if CHANNEL_NAME in channel.name]

if channel_idx:
    # the specific channel exists, return the data
    channel_idx = channel_idx[0]
    time = specv.channel[channel_idx].isotope_ratios.time - 32
    density_ratio = specv.channel[channel_idx].isotope_ratios.isotope[0].density_ratio

fix, ax = plt.subplots()
ax.plot(time, density_ratio*100)
ax.set_xlabel('time [s]')
ax.set_ylabel('nH/nD [%]')




