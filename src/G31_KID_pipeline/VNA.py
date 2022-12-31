import numpy as np
from . import datapaths
from . import functions as fc
import sys
from glob import glob

'''
            VNA class
'''
class VNA():
    def __init__(self, filename, temperature=None, build_dataset=False, remove_baseline=True, label=None):
        '''
        Class for reading the ROACH VNA sweeps.

        Parameters
        ----------
        filename : string
            The filename of the VNA sweep.
        temperature : float, optional
            The temperature of the VNA sweep in mK. The default is None.
        build_dataset : boolean, optional
            If True the raw data from the ROACH is converted in .npy files. 
            The default is False.
        remove_baseline : boolean, optional
            If True the baselines from amplitude and phase data will be removed. 
            The default is True.
        label : string, optional
            Label string of the VNA sweep.
            The default is None.

        Returns
        -------
        None.

        '''
        print('++++++++++++++++++++++++++++++++++++++++++++++++')
        print('Reding VNA sweep {:s}...'.format(filename))
        self.filename = filename
        self.temperature = temperature
        self.entry = []
        
        try:
            self.bb_freqs = np.load(datapaths.vna / self.filename / 'bb_freqs.npy')
        except:
            print("Cannot find the target directory '{:s}'.".format(self.filename))
            sys.exit()
        
        if not (datapaths.vna_S21 / self.filename).exists() or build_dataset:
            fc.buildS21Dataset(self)
        
        # remove baseline
        if remove_baseline:
            self.removeBaseline()
            
        # check for extracted target and load data
        if (datapaths.vna_S21 / self.filename / 'extracted_target').exists():
            print("Extracted target found. Reading S21 data...")
            entries = list((datapaths.vna_S21 / self.filename / 'extracted_target').glob('[0123456789]*'))
            print("Found", len(entries), "entries.")
            self.readS21Data(len(entries))
        
        if label == None:
            self.label = filename
        else:
            self.label = label
            
            
    def removeBaseline(self, mag_filt='airp_lss', phase_filt='lowpass_cosine'):
        '''
        Remove the baseline of the VNA sweep in order to find resonances

        Parameters
        ----------
        mag_filt : string, optional
            Filter algorithm for the amplitude data. It can be:
                - 'lowpass_cosine' for a low pass filter algorithm;
                - 'a_lss' for an asymmetric leas square smoothing algorithm;
                - 'airp_lss' for an adaptive iteratively reweighted penalized
                    least square smoothing algorithm.
            default is 'airp_lss'.
        phase_filt : string, optional
            Same description as the mag_filt parameter. 'lowpass_cosine' works
            better. The default is 'lowpass_cosine'.

        Returns
        -------
        mag : numpy array
            Filtered amplitude.
        phase : numpy arra
            Filtered phase.

        '''
        if mag_filt==None:
            mag_filt = 'None'
        if phase_filt==None:
            phase_filt = 'None'
        print("Removing baselines...")
        print('Magnitude baseline correction algorithm: '+mag_filt)
        print('Phase baseline correction algorithm: '+phase_filt)
    
        # how many channels?
        n_chan = self.bb_freqs.size
        
        freqs = []
        mag = []
        phase = []
        I = []
        Q = []
        
        mag_channel_last = 0.0
        phase_channel_last = 0.0
    
        for chan in range(n_chan):
            
            # read one sweep at a time
            freqs_channel = np.load(datapaths.vna_S21 / self.filename / "{:03d}".format(chan) / "freqs.npy")
            mag_channel = np.load(datapaths.vna_S21 / self.filename / "{:03d}".format(chan) / "mag.npy")
            phase_channel = np.load(datapaths.vna_S21 / self.filename / "{:03d}".format(chan) / "phase.npy")
            I_channel = np.load(datapaths.vna_S21 / self.filename / "{:03d}".format(chan) / "I.npy")
            Q_channel = np.load(datapaths.vna_S21 / self.filename / "{:03d}".format(chan) / "Q.npy")
            
            # remove offsets
            if chan != 0:
                delta = mag_channel_last-mag_channel[0]
                mag_channel = [y+delta for y in mag_channel]
                
                delta = phase_channel_last-phase_channel[0]
                phase_channel = [y+delta for y in phase_channel]
            
            freqs.append(freqs_channel)
            mag.append(mag_channel)
            phase.append(phase_channel)
            I.append(I_channel)
            Q.append(Q_channel)
            
            mag_channel_last = mag_channel[-1]
            phase_channel_last = phase_channel[-1]
            
        freqs = np.hstack(freqs)
        mag = np.hstack(mag)
        mag -= mag[0]
        phase = np.hstack(phase)
        phase = np.unwrap(phase)
        phase -= phase[0]
        
        I = np.hstack(I)
        Q = np.hstack(Q)
        
        self.freqs = freqs
        self.mag = mag
        self.phase = phase
        self.I = I
        self.Q = Q
        
        # phase detrend
        from scipy.signal import detrend
        detrend(phase, axis=-1, type='linear', overwrite_data=True)
        phase -= phase[0]
        
        # mag filter
        if mag_filt == 'lowpass_cosine':
            sweep_step = 1.25 # kHz
            smoothing_scale = 2500.0 # kHz
            filtered = fc.lowpass_cosine(mag, sweep_step, 1.0/smoothing_scale, 0.1 * (1.0/smoothing_scale))
            mag -= filtered    
            
        if mag_filt == 'a_lss':
            baseline = fc.asymmetric_least_squares_smoothing(data=mag, lam=1e6, p=8e-4, N_iter=5)
            mag -= baseline
            
        if mag_filt == 'airp_lss':
            baseline = fc.adaptive_iteratively_reweighted_penalized_least_squares_smoothing(data=mag, lam=1e6, N_iter=5)
            mag -= baseline
        
        
        # phase filter
        if phase_filt == 'lowpass_cosine':
            sweep_step = 1.25 # kHz
            smoothing_scale = 2500.0 # kHz
            filtered = fc.lowpass_cosine(phase, sweep_step, 1./smoothing_scale, 0.1 * (1.0/smoothing_scale))
            phase -= filtered
    
        return mag, phase
    
    def plotSweep(self, xlim=None, mag_filt='airp_lss', phase_filt='lowpass_cosine'):
        '''
        This function plots both the amplitude and phase data of the VNA sweep.

        Parameters
        ----------
        xlim : tuple, optional
            x axis limits. The default is [0, 700].
        mag_filt : string, optional
            Filter algorithm for the amplitude data. It can be:
                - 'lowpass_cosine' for a low pass filter algorithm;
                - 'a_lss' for an asymmetric leas square smoothing algorithm;
                - 'airp_lss' for an adaptive iteratively reweighted penalized
                    least square smoothing algorithm.
            default is 'airp_lss'.
        phase_filt : string, optional
            Same description as the mag_filt parameter. 'lowpass_cosine' works
            better. The default is 'lowpass_cosine'.

        Returns
        -------
        None.

        '''
        print("Plot VNA sweep...")
        mag, phase = self.removeBaseline(mag_filt, phase_filt)
        
        from matplotlib import pyplot as plt
        fig = plt.figure()
        fig.set_size_inches(7.5, 8)
        plt.subplots_adjust(bottom=0.15, right=0.98, top=0.95, left=0.1, wspace=0.35)
        ax0 = plt.subplot(211)
        ax1 = plt.subplot(212, sharex=ax0)
        
        ax0.plot(self.freqs, mag, color='black', linewidth=1)
        ax1.plot(self.freqs, phase, color='black', linewidth=1)
        
        ax0.grid(linestyle='-', alpha=0.5)
        ax0.set_xlim([self.freqs[0], self.freqs[-1]])
        ax0.set_ylabel('Mag [dB]')
        ax0.set_xlabel('Frequency [MHz]')
        
        ax1.grid(linestyle='-', alpha=0.5)
        ax1.set_xlim([self.freqs[0], self.freqs[-1]])
        ax1.set_ylabel('Phase [rad]')
        ax1.set_xlabel('Frequency [MHz]')
        
        plt.show()
        
    
    # read S21 data for all the tones
    def readS21Data(self, N_channels):
        from uncertainties import ufloat
        known_resonaces_not_fitted = 0
        for i in range(N_channels):
            try:
                file_path = datapaths.vna_S21 / self.filename / 'extracted_target' / '{:03d}'.format(i)
                [Rea_n, Rea_s, Ima_n, Ima_s, Q_tot_n, Q_tot_s, Q_c_n, Q_c_s, 
                 Q_i_n, Q_i_s, nu_r_n, nu_r_s, phi_0_n, phi_0_s, tau] = np.load(file_path / "complex_parameters.npy", allow_pickle=False)
                
                freqs = np.load(file_path/'freqs.npy')
                mag = np.load(file_path/'mag.npy')
                
                self.entry.append({'target_freq': (freqs[0]+freqs[-1])*0.5,
                                  'depth': min(mag),
                                  'reduced_chi2': float(np.load(file_path / "reduced_chi2.npy", allow_pickle=False)),
                                  'Re[a]': ufloat(Rea_n, Rea_s),
                                  'Im[a]': ufloat(Ima_n, Ima_s),
                                  'Q_tot': ufloat(Q_tot_n, Q_tot_s),
                                  'Q_c': ufloat(Q_c_n, Q_c_s),
                                  'Q_i': ufloat(Q_i_n, Q_i_s),
                                  'nu_r': ufloat(nu_r_n, nu_r_s),
                                  'phi_0': ufloat(phi_0_n, phi_0_s)})
            except: 
                self.entry.append({'target_freq': (freqs[0]+freqs[-1])*0.5,
                                  'depth': min(mag),
                                  'reduced_chi2': None,
                                  'Re[a]': None,
                                  'Im[a]': None,
                                  'Q_tot': None,
                                  'Q_c': None,
                                  'Q_i': None,
                                  'nu_r': None,
                                  'phi_0': None})
                known_resonaces_not_fitted += 1
    
        if known_resonaces_not_fitted>0:
            print("Complex fit parameters not available for {:d}/{:d} resonances.".format(known_resonaces_not_fitted, N_channels))
    
        
        
    def findPeaks(self, xlim=[0, 700], mag_filt='airp_lss', phase_filt='lowpass_cosine', peak_width=(1.0, 150.0), peak_height=1.0, peak_prominence=(1.0, 30.0)):
        '''
        A function that searh for resonances in the VNA sweep

        Parameters
        ----------
        xlim : tuple, optional
            x axis limits. The default is [0, 700].
        mag_filt : string, optional
            Filter algorithm for the amplitude data. It can be:
                - 'lowpass_cosine' for a low pass filter algorithm;
                - 'a_lss' for an asymmetric leas square smoothing algorithm;
                - 'airp_lss' for an adaptive iteratively reweighted penalized
                    least square smoothing algorithm.
            default is 'airp_lss'.
        phase_filt : string, optional
            Same description as the mag_filt parameter. 'lowpass_cosine' works
            better. The default is 'lowpass_cosine'.
        peak_width : float of tuple, optional
            Minimum width of peaks or range of values. The default is 
            (1.0, 150.0).
        peak_height : float or tuple, optional
            Minimum height of peaks or range of values. The default is 1.0.
        peak_prominence : float or tuple, optional
            Minimum prominence of peaks or range of values. The default is (1.0, 30.0).
            
        For other information see scipy.signal.find_peaks() help.

        Returns
        -------
        peaks : nunpy array
            index of peaks.
        peaks_info : dictionary
            info on peaks.

        '''
        print(self.filename+': peak finding...')
        
        mag, phase = self.removeBaseline(mag_filt, phase_filt)
        
        # peak finding
        from scipy.signal import find_peaks
        peaks, peaks_info = find_peaks(-mag, 
                                       width=peak_width, 
                                       height=peak_height, 
                                       prominence=peak_prominence)
        print("Found ", len(peaks), " resonances")
        
        from matplotlib import pyplot as plt
        fig = plt.figure()
        fig.set_size_inches(7.5, 8)
        plt.subplots_adjust(bottom=0.15, right=0.98, top=0.95, left=0.1, wspace=0.35)
        ax0 = plt.subplot(211)
        ax1 = plt.subplot(212, sharex=ax0)
        
        ax0.plot(self.freqs[peaks], mag[peaks], "x", label="Resonance")
        ax0.vlines(x=self.freqs[peaks], ymin=mag[peaks], ymax=mag[peaks]+peaks_info["prominences"], color='red', ls='--', lw=1)
        
        for left, right, height in zip(peaks_info["left_ips"],peaks_info["right_ips"], peaks_info["width_heights"]):
            xmin = self.freqs[int(left)]
            xmax = self.freqs[int(right)]
            ax0.hlines(y=-height, xmin=xmin, xmax=xmax, color='blue', ls='--', lw=1)
        
        ax0.legend(loc='best')

        ax0.plot(self.freqs, mag, color='black', linewidth=1)
        ax1.plot(self.freqs, phase, color='black', linewidth=1)
        
        ax0.grid(linestyle='-', alpha=0.5)
        ax0.set_xlim([self.freqs[0], self.freqs[-1]])
        ax0.set_ylabel('Mag [dB]')
        ax0.set_xlabel('Frequency [MHz]')
        
        ax1.grid(linestyle='-', alpha=0.5)
        ax1.set_xlim([self.freqs[0], self.freqs[-1]])
        ax1.set_ylabel('Phase [rad]')
        ax1.set_xlabel('Frequency [MHz]')
        
        plt.show()
        
        return peaks, peaks_info
        
        
    def extractTarget(self, peaks, peaks_info, extr_width=2.5):
        print("Extracting target information...")
        from matplotlib import pyplot as plt
        from matplotlib import cm
        
        cmap = cm.get_cmap('tab20b', lut=None)
        
        fig = plt.figure()
        fig.set_size_inches(7.5, 8)
        plt.subplots_adjust(bottom=0.15, right=0.98, top=0.95, left=0.1, wspace=0.35)
        ax0 = plt.subplot(211)
        ax1 = plt.subplot(212, sharex=ax0)
        
        
        I = []
        Q = []
        mag = []
        phase = []
        freqs = []
        self.entry = []
        
        
        for i,(peak, left_ips, right_ips) in enumerate(zip(peaks, peaks_info["left_ips"], peaks_info["right_ips"])):
            xmin = self.freqs[int(left_ips)]
            xmax = self.freqs[int(right_ips)]
            
            width = xmax - xmin
            
            keep = np.logical_and(self.freqs>=xmin-extr_width*width, self.freqs<=xmax+extr_width*width)
            
            # save I and Q data
            output_path = datapaths.vna_S21 / self.filename / 'extracted_target' / '{:03d}'.format(i)
            if not output_path.exists():
                output_path.mkdir(parents=True)
                
            np.save(output_path / "I.npy", self.I[keep])
            np.save(output_path / "Q.npy", self.Q[keep])
            np.save(output_path / "freqs.npy", self.freqs[keep])
            np.save(output_path / "mag.npy", self.mag[keep])
            
            dictionary = {'channel': i,
                          'target_freq': self.freqs[peak],
                          'depth': self.mag[peak],
                          'Re[a]': None,
                          'Im[a]': None,
                          'Q_tot': None,
                          'Q_c': None,
                          'Q_i': None,
                          'nu_r': None,
                          'phi_0': None,
                          'reduced_chi2': None}
            
            self.entry.append(dictionary)
            
            I.append(np.asarray(self.I[keep]))
            Q.append(np.asarray(self.Q[keep]))
            mag.append(np.asarray(self.mag[keep]))
            phase.append(np.asarray(self.phase[keep]))
            freqs.append(np.asarray(self.freqs[keep]))
            
            color = cmap(np.random.rand())
            
            ax0.plot(self.freqs[keep], self.mag[keep], color=color, linewidth=1)
            ax1.plot(self.freqs[keep], self.phase[keep], color=color, linewidth=1)
            
        
        ax0.grid(linestyle='-', alpha=0.5)
        ax0.set_xlim([self.freqs[0], self.freqs[-1]])
        ax0.set_ylabel('Mag [dB]')
        ax0.set_xlabel('Frequency [MHz]')
        
        ax1.grid(linestyle='-', alpha=0.5)
        ax1.set_xlim([self.freqs[0], self.freqs[-1]])
        ax1.set_ylabel('Phase [rad]')
        ax1.set_xlabel('Frequency [MHz]')
        
        plt.show()
    
    
    def fitS21(self, channel, DATAPOINTS=50):
        print("")
        from tqdm import tqdm
        
        if channel == 'all':
            pbar = tqdm(self.entry, position=0, leave=True)
            for e in pbar:
                pbar.set_description(self.filename+" complex fit... ")
                out_path = datapaths.vna_S21 / self.filename / 'extracted_target' /  "{:03d}".format(e['channel'])
                try:
                    I = np.load(out_path / 'I.npy')
                    Q = np.load(out_path / 'Q.npy')
                    freqs = np.load(out_path / 'freqs.npy')
                    params, chi2 = fc.complexS21Fit(I=I, Q=Q, freqs=freqs, res_freq=e['target_freq'], 
                                           output_path=out_path, DATAPOINTS=DATAPOINTS)
                    
                    e['Re[a]'] = params['Re[a]']
                    e['Im[a]'] = params['Im[a]']
                    e['Q_tot'] = params['Q_tot']
                    e['Q_c'] = params['Q_c']
                    e['Q_i'] = params['Q_i']
                    e['nu_r'] = params['nu_r']
                    e['phi_0'] = params['phi_0']
                    e['reduced_chi2'] = float(chi2)
                except:
                    pass
        else:
            out_path = datapaths.vna_S21 / self.filename / 'extracted_target' /  "{:03d}".format(channel)
            I = np.load(out_path / 'I.npy')
            Q = np.load(out_path / 'Q.npy')
            freqs = np.load(out_path / 'freqs.npy')
            params, chi2 = fc.complexS21Fit(I=I, Q=Q, freqs=freqs, res_freq=self.entry[channel]['target_freq'], 
                          output_path=out_path, DATAPOINTS=DATAPOINTS, verbose=True)
                
            self.entry[channel]['Re[a]'] = params['Re[a]']
            self.entry[channel]['Im[a]'] = params['Im[a]']
            self.entry[channel]['Q_tot'] = params['Q_tot']
            self.entry[channel]['Q_c'] = params['Q_c']
            self.entry[channel]['Q_i'] = params['Q_i']
            self.entry[channel]['nu_r'] = params['nu_r']
            self.entry[channel]['phi_0'] = params['phi_0']
            self.entry[channel]['reduced_chi2'] = float(chi2)
            
        return
            
        
    def plotS21(self, channel):
        extracted_target_path = datapaths.vna_S21 / self.filename / 'extracted_target' / '{:03d}'.format(channel)
        
        fc.complexS21Plot(extracted_target_path)
        
        
        