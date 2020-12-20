import h5py
import numpy as np

# Assume that WEST libraries pywed, IRFMtb and PPPAT are available in the Python Path
import pywed as pw
from IRFMtb import tsdernier_choc
from pppat.control_room.signals import get_sig, signals, is_sig
from tqdm import tqdm

class PulseDB():
    def __init__(self, hdf5_filename):
        self.hdf5_filename = hdf5_filename   
        
    @property
    def pulse_list(self) -> list:
        """
        Return the WEST pulse numbers as a list
        
        Returns
        -------
        pulse_list: list
            List of the WEST pulse numbers
        """
        with h5py.File(self.hdf5_filename, 'a') as f:
            return list(f.attrs['pulse_list'])

    @pulse_list.setter
    def pulse_list(self, pl: list):
        """
        Set the list of WEST pulse numbers
        
        Argument
        --------
        pl: list
            List of the WEST pulse numbers
        """
        with h5py.File(self.hdf5_filename, 'a') as f:
            f.attrs.create('pulse_list', pl)
            
    def update_pulse_list(self, ref_sig=signals['RF_P_tot'], thres=0.2, from_pulse=None):
        """
        Update the database by adding new pulses if any.
        
        Arguments
        ---------
        ref_sig: signal dict (any from pppat.control_room.signals.signals)
            Reference signal to use to discriminate if pulse must be keept or not
        thres: float
            Threshold to check the reference signal with.
        from_pulse: int
            Start from the indicated WEST number. Default is None (starts from latest entry)
        """
        # if pulse_list already exist
        try:
            pulse_list = self.pulse_list
        except KeyError as e:
            pulse_list = []
        
        if pulse_list:  # if not empty
            pulse_start = pulse_list[-1]
        else:
            pulse_start = 56287  # ~ first plasma with IC for C5

        if from_pulse:  # force starting from a given pulse number if indicated
            pulse_start = from_pulse
            
        pulse_end = tsdernier_choc()
        pulses = range(pulse_start, pulse_end+1)
        print(f'Going from {pulse_start} to {pulse_end}')

        # then append from the latest pulse
        for pulse in tqdm(pulses):
            # keep only pulses where threshold is fit
            res = is_sig(pulse, ref_sig, thres)
            if res:
                print(f'Keeping {pulse}')
                pulse_list.append(pulse)
            else:
                print(f'Discarding {pulse}')
        # update the db
        self.pulse_list = pulse_list
            

    def delete_signals(self, signal_names, pulses=None):
        """
        Delete signals (ie a hdf5 group containing both (y,t) dataset)
        for the given pulses. If no pulse list given, all pulses considered

        Arguments
        ---------
        signal_names: list of string
            List of signal to remove from the database
        pulses: list of int
            List of pulses to consider. Default value is None (all pulses)
        """
        if type(signal_names) is not list:
            raise ValueError(f'signal_names must be a list')
        
        with h5py.File(self.hdf5_filename, 'a') as f:
            if not pulses:
                pulses = self.pulse_list

            for pulse in pulses:
                for signal_name in signal_names:
                    sig_path = f'{pulse}/{signal_name}'
                    if sig_path in f:
                        del f[sig_path]
                    else:
                        print(f'nothing to supress on {sig_path}')

    def add_signals(self, pulse, signal_names, force_rewrite=False, do_smooth=False):
        """
        Add signal data into the database
        f : hdf5 file
        pulse : int
        signal_names : list of string
        """
        with h5py.File(self.hdf5_filename, 'a') as f:
            for sig_name in signal_names:
                # do not rewrite if allready exist
                sig_path = f'{pulse}/{sig_name}'
                if (sig_path not in f) or force_rewrite:
                    print(f'Getting {sig_name} for #{pulse}')
                    y, t = get_sig(pulse, signals[sig_name], do_smooth=do_smooth)
                    if pulse == 54719:
                        print(y)
                    # store y,t in the database
                    if sig_path+"/y" in f:
                        f[sig_path+"/y"][...] = y
                    else:
                        f[sig_path+"/y"] = y

                    if sig_path+"/t" in f:
                        f[sig_path+"/t"][...] = t
                    else:
                        f[sig_path+"/t"] = t

                else:
                    print(f'{sig_name} already exist in database for #{pulse}: passing...')

    def add_attr(self, pulse, signal_name, attr_name, attr_value, force_rewrite=False):
        """
        Write an HDF5 attribute to a signal group
        """
        with h5py.File(self.hdf5_filename, 'a') as f:
            sig_path = f'{pulse}/{signal_name}'
            if sig_path in f:
                f[sig_path].attrs.create(attr_name, attr_value)
            else:
                raise KeyError(f'{sig_path} does not exist')

    def get_attr(self, pulse, signal_name, attr_name):
        """
        Return an HDF5 attribute
        """
        with h5py.File(self.hdf5_filename, 'r') as f:
            sig_path = f'{pulse}/{signal_name}'
            if sig_path in f:
                return f[sig_path].attrs[attr_name]
            else:
                raise KeyError(f'{sig_path} does not exist')

    def get_signal(self, pulse, signal_name):
        """
        Return the (y,t) of signal for a given pulse
        
        Arguments
        ---------
        pulse: int
            WEST pulse number
        signal_name: str
            Name of the signal
        
        Returns
        -------
        y: np.array
            value array
        t: np.array
            time array
        
        """
        with h5py.File(self.hdf5_filename, 'r') as f:
            sig_path = f'{pulse}/{signal_name}'
            if sig_path in f:
                return f[sig_path+'/y'][:], f[sig_path+'/t'][:]
            else:
                raise KeyError(f'signal {sig_path} does not exist')
                
    def list_signal(self, pulse):
        """
        List all the signals available for a pulse.
        
        Argument
        --------
        pulse: int
            WEST pulse number
        
        """
        with h5py.File(self.hdf5_filename, 'r') as f:
            if str(pulse) in f:
                return list(f[f'{pulse}'].keys())
            else:
                raise KeyError(f'Pulse {pulse} does not exist')