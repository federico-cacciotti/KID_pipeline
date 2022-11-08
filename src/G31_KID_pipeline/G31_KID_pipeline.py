import numpy as np
from matplotlib.lines import Line2D
import sys
from pathlib import Path


# plot configuration
#plt.rc('text', usetex=False)                        #
#plt.rc('font', family='serif')
#plt.rcParams['text.usetex'] = False                 #
#plt.rcParams['text.latex.preamble'] = r"\usepackage{amssymb} \usepackage{siunitx}"
#plt.rcParams.update({'figure.max_open_warning': 0})
#plt.rcParams["axes.formatter.use_mathtext"] = True
#plt.minorticks_on()
#plt.tick_params(axis='both', which='both', direction='in', labelleft='on', labelright='off', bottom=True, top=True, left=True, right=True)
#plt.xkcd()
#plt.close('all')


#NEP_photon = 5.e-15 #W/Hz^1/2
#NEP_photon_noatm = 8.e-17 #W/Hz^1/2


'''
            VNA class
'''
class VNA():
    def __init__(self, filename, temperature=None, build_dataset=False):
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

        Returns
        -------
        None.

        '''
        print('Reding VNA sweep {:s}...'.format(filename))
        self.filename = filename
        self.temperature = temperature
        
        try:
            self.bb_freqs = np.load(paths.vna / self.filename / 'bb_freqs.npy')
        except:
            print("Cannot find the target directory '{:s}'.".format(self.filename))
            sys.exit()
        
        if not (paths.vna_S21 / self.filename).exists() or build_dataset:
            buildS21Dataset(self)
            
            
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
            freqs_channel = np.load(paths.vna_S21 / self.filename / "{:03d}".format(chan) / "freqs.npy")
            mag_channel = np.load(paths.vna_S21 / self.filename / "{:03d}".format(chan) / "mag.npy")
            phase_channel = np.load(paths.vna_S21 / self.filename / "{:03d}".format(chan) / "phase.npy")
            I_channel = np.load(paths.vna_S21 / self.filename / "{:03d}".format(chan) / "I.npy")
            Q_channel = np.load(paths.vna_S21 / self.filename / "{:03d}".format(chan) / "Q.npy")
            
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
            filtered = lowpass_cosine(mag, sweep_step, 1.0/smoothing_scale, 0.1 * (1.0/smoothing_scale))
            mag -= filtered    
            
        if mag_filt == 'a_lss':
            baseline = asymmetric_least_squares_smoothing(data=mag, lam=1e6, p=8e-4, N_iter=5)
            mag -= baseline
            
        if mag_filt == 'airp_lss':
            baseline = adaptive_iteratively_reweighted_penalized_least_squares_smoothing(data=mag, lam=1e6, N_iter=5)
            mag -= baseline
        
        
        # phase filter
        if phase_filt == 'lowpass_cosine':
            sweep_step = 1.25 # kHz
            smoothing_scale = 2500.0 # kHz
            filtered = lowpass_cosine(phase, sweep_step, 1./smoothing_scale, 0.1 * (1.0/smoothing_scale))
            phase -= filtered
    
        return mag, phase
    
    def plotVNA(self, xlim=[0, 700], mag_filt='airp_lss', phase_filt='lowpass_cosine'):
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
        print('Peak finding...')
        
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
        
        xmin = (peaks_info["left_ips"]/(1+self.freqs.size)) * (self.freqs[-1]-self.freqs[0]) + self.freqs[0]
        xmax = (peaks_info["right_ips"]/(1+self.freqs.size)) * (self.freqs[-1]-self.freqs[0]) + self.freqs[0]
        ax0.hlines(y=-peaks_info["width_heights"], xmin=xmin, xmax=xmax, color='blue', ls='--', lw=1)
        
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
        from pathlib import Path
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
        res_freq = []
        
        
        for i,(peak, left_ips, right_ips) in enumerate(zip(peaks, peaks_info["left_ips"], peaks_info["right_ips"])):
            xmin = (left_ips/(1+self.freqs.size)) * (self.freqs[-1]-self.freqs[0]) + self.freqs[0]
            xmax = (right_ips/(1+self.freqs.size)) * (self.freqs[-1]-self.freqs[0]) + self.freqs[0]
            
            width = xmax - xmin
            
            keep = np.logical_and(self.freqs>=xmin-extr_width*width, self.freqs<=xmax+extr_width*width)
            
            # save I and Q data
            output_path = paths.vna_S21 / self.filename / 'extracted_target' / '{:03d}'.format(i)
            if not output_path.exists():
                output_path.mkdir(parents=True)
                
            np.save(output_path / "I.npy", self.I[keep])
            np.save(output_path / "Q.npy", self.Q[keep])
            np.save(output_path / "freqs.npy", self.freqs[keep])
            
            I.append(np.asarray(self.I[keep]))
            Q.append(np.asarray(self.Q[keep]))
            mag.append(np.asarray(self.mag[keep]))
            phase.append(np.asarray(self.phase[keep]))
            freqs.append(np.asarray(self.freqs[keep]))
            res_freq.append(self.freqs[peak])
            
            color = cmap(np.random.rand())
            
            ax0.plot(self.freqs[keep], self.mag[keep], color=color, linewidth=1)
            ax1.plot(self.freqs[keep], self.phase[keep], color=color, linewidth=1)
            
        
        output = {"I": I, "Q": Q, "mag": mag, "phase":phase, "freqs": freqs, "res_freq": res_freq}
        
        ax0.grid(linestyle='-', alpha=0.5)
        ax0.set_xlim([self.freqs[0], self.freqs[-1]])
        ax0.set_ylabel('Mag [dB]')
        ax0.set_xlabel('Frequency [MHz]')
        
        ax1.grid(linestyle='-', alpha=0.5)
        ax1.set_xlim([self.freqs[0], self.freqs[-1]])
        ax1.set_ylabel('Phase [rad]')
        ax1.set_xlabel('Frequency [MHz]')
        
        plt.show()
        
        return output
    
    
    def fitS21(self, extracted_target, channel):
        print("")
        from tqdm import tqdm
        
        I = extracted_target['I']
        Q = extracted_target['Q']
        freqs = extracted_target['freqs']
        res_freq = extracted_target['res_freq']
        
        if channel == 'all':
            pbar = tqdm(res_freq, position=0, leave=True)
            for i,e in enumerate(pbar):
                pbar.set_description("Complex fit... ")
                out_path = paths.vna_S21 / self.filename / 'extracted_target' /  "{:03d}".format(i)
                if not out_path.exists():
                    out_path.mkdir(parents=True)
                    
                try:
                    params, chi2 = complexS21Fit(I=I[i], Q=Q[i], freqs=freqs[i], res_freq=res_freq[i], 
                                           output_path=out_path, DATAPOINTS=20)
                except:
                    pass
        else:
            out_path = paths.vna_S21 / self.filename / 'extracted_target' /  "{:03d}".format(channel)
            if not out_path.exists():
                out_path.mkdir(parents=True)
            
            try:
                params, chi2 = complexS21Fit(I=I[channel], Q=Q[channel], freqs=freqs[channel], res_freq=res_freq[channel], 
                              output_path=out_path, DATAPOINTS=20, verbose=True)
            except:
                print("Not able to perform a complex fit.")
                return
            
        
    def plotS21(self, channel):
        extracted_target_path = paths.vna_S21 / self.filename / 'extracted_target' / '{:03d}'.format(channel)
        
        complexS21Plot(extracted_target_path)
        


'''
            Target class
'''
class Target():
    def __init__(self, filename, temperature=None, build_dataset=False, label=None, out_of_resonance_parameter=2.0):
        print("+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +")
        print('Reding target sweep {:s}...'.format(filename))
        self.filename = filename
        self.temperature = temperature
        
        if label != None:
            self.label = label
        else:
            self.label = filename
        
        self.in_res_target_idx_aligned = []
        
        try:
            self.target_freqs_new = np.load(paths.target / self.filename / 'target_freqs_new.npy')
        except:
            print("Cannot find the target directory '{:s}'.".format(self.filename))
            sys.exit()
            
        self.entries = self.target_freqs_new.size
        print('{:d} entries found.'.format(self.entries))
        
        # building S21 dataset
        if not (paths.target_S21 / self.filename).exists() or build_dataset:
            buildS21Dataset(self)
        
        # building data dictonaries
        self.entry = []
        for channel,f in enumerate(self.target_freqs_new):
            mag = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "mag.npy")
            
            self.entry.append({'target_freq': f,
                               'channel': channel,
                               'depth': max(mag)-min(mag),
                               'is_out_of_res': False,
                               'number_of_peaks': 0,
                               'Re[a]': None, 
                               'Im[a]': None, 
                               'Q_tot': None, 
                               'Q_i': None, 
                               'Q_c': None,
                               'nu_r': None, 
                               'phi_0': None,
                               'reduced_chi2': None})
        
        self.out_of_res = self.filterOutOfResTones(std_mult=out_of_resonance_parameter)
        
        self.readS21Data()
        
        self.findDouble()
        print("+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +")
    

    # filter out of resonance tones refusing them, sigma_mult allows to refuse out of resonance tones
    def filterOutOfResTones(self, std_mult=3.0):
        print("\nChecking for out of resonance tones...")
        
        depths = [self.entry[i]['depth'] for i in range(self.entries)]
        average_depth = np.average(depths)
        print("\taverage depth: {:f}dB".format(average_depth))
        std_depth = np.std(depths)
        print("\tstd depth: {:f}dB".format(std_depth))
        
        threshold = average_depth-std_mult*std_depth
        print("\tdiscarding tones with depth less than {:f}dB".format(threshold))
        
        out_of_res_number = 0
        for e in self.entry:
            if e['depth']<=threshold:
                e['is_out_of_res'] = True
                out_of_res_number += 1
        
        print("\tFound {:d} out of resonance tones: ".format(out_of_res_number))
        if out_of_res_number>0:
            for i,e in enumerate(self.entry):
                if e['is_out_of_res']:
                    print("\t\tChannel {:d} - {:f} MHz".format(i, e['target_freq']))
                    
        return out_of_res_number

    def fitS21(self, channel, DATAPOINTS=70):
        print("")
        from tqdm import tqdm
        
        path = paths.target_S21 / self.filename
        
        if channel == 'all':
            pbar = tqdm(self.entry, position=0, leave=True)
            for e in pbar:
                pbar.set_description("Complex fit... ")
                if not e['is_out_of_res']:
                    c = e['channel']
                    I = np.load(path / "{:03d}".format(c) / "I.npy", allow_pickle=True)
                    Q = np.load(path / "{:03d}".format(c) / "Q.npy", allow_pickle=True)
                    freqs = np.load(path / "{:03d}".format(c) / "freqs.npy", allow_pickle=True)
                    in_path = paths.target / self.filename / "{:03d}".format(c)
                    out_path = paths.target_S21 / self.filename / "{:03d}".format(c)
                    
                    try:
                        params, chi2 = complexS21Fit(I=I, Q=Q, freqs=freqs, res_freq=e['target_freq'], 
                                               output_path=out_path, DATAPOINTS=DATAPOINTS)
                    
                        if params != None:
                            e['Re[a]'] = params['Re[a]']
                            e['Im[a]'] = params['Im[a]']
                            e['Q_tot'] = params['Q_tot']
                            e['Q_c'] = params['Q_c']
                            e['Q_i'] = params['Q_i']
                            e['nu_r'] = params['nu_r']
                            e['phi_0'] = params['phi_0']
                            e['reduced_chi2'] = chi2
                    except:
                        pass
        else:
            I = np.load(path / "{:03d}".format(channel) / "I.npy", allow_pickle=True)
            Q = np.load(path / "{:03d}".format(channel) / "Q.npy", allow_pickle=True)
            freqs = np.load(path / "{:03d}".format(channel) / "freqs.npy", allow_pickle=True)
            in_path = paths.target / self.filename / "{:03d}".format(channel)
            out_path = paths.target_S21 / self.filename / "{:03d}".format(channel)
            
            try:
                params, chi2 = complexS21Fit(I=I, Q=Q, freqs=freqs, res_freq=self.entry[channel]['target_freq'], 
                              output_path=out_path, DATAPOINTS=DATAPOINTS, verbose=True)
                
                self.entry[channel]['Re[a]'] = params['Re[a]']
                self.entry[channel]['Im[a]'] = params['Im[a]']
                self.entry[channel]['Q_tot'] = params['Q_tot']
                self.entry[channel]['Q_c'] = params['Q_c']
                self.entry[channel]['Q_i'] = params['Q_i']
                self.entry[channel]['nu_r'] = params['nu_r']
                self.entry[channel]['phi_0'] = params['phi_0']
                self.entry[channel]['reduced_chi2'] = chi2
            
            except:
                print("Not able to perform a complex fit.")
                return
        
    # read S21 data for all the tones
    def readS21Data(self):
        from uncertainties import ufloat
        known_resonaces_not_fitted = 0
        for i,e in enumerate(self.entry):
            if not e['is_out_of_res']:
                try:
                    file_path = paths.target_S21 / self.filename / '{:03d}'.format(i)
                    [Rea_n, Rea_s, Ima_n, Ima_s, Q_tot_n, Q_tot_s, Q_c_n, Q_c_s, 
                     Q_i_n, Q_i_s, nu_r_n, nu_r_s, phi_0_n, phi_0_s, tau] = np.load(file_path / "complex_parameters.npy", allow_pickle=False)
                    
                    e['reduced_chi2'] = np.load(file_path / "reduced_chi2.npy", allow_pickle=False)
                    
                    e['Re[a]'] = ufloat(Rea_n, Rea_s)
                    e['Im[a]'] = ufloat(Ima_n, Ima_s)
                    e['Q_tot'] = ufloat(Q_tot_n, Q_tot_s)
                    e['Q_c'] = ufloat(Q_c_n, Q_c_s)
                    e['Q_i'] = ufloat(Q_i_n, Q_i_s)
                    e['nu_r'] = ufloat(nu_r_n, nu_r_s)
                    e['phi_0'] = ufloat(phi_0_n, phi_0_s)
                except: 
                    known_resonaces_not_fitted += 1
    
        if known_resonaces_not_fitted>0:
            print("\nComplex fit parameters not available for {:d}/{:d} resonances.".format(known_resonaces_not_fitted, self.entries-self.out_of_res))
    
    
    def findDouble(self):
        print("\nSearching for double resonances...")
        double_found = False
        
        from scipy.signal import find_peaks
        
        for e in self.entry:
            if not e['is_out_of_res']:
                # read one sweep at a time
                channel = e['channel']
                #x_data_chan = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "freqs.npy")
                y_data_chan = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "mag.npy")
            
                peaks, info = find_peaks(-y_data_chan, width=1, height=2.0, prominence=2.0)
                
                n_peaks = len(peaks)
                
                if n_peaks>1:
                    print("\tFound ", n_peaks, " resonances in channel ", channel)
                
                e['number_of_peaks'] = n_peaks
                if n_peaks>1:
                    double_found = True
        
        if not double_found:
            print("\tNo double resonance found!")
    
    
    def entry_parameters_to_table(self):
        print('temperature target_freq channel depth is_out_of_res number_of_peaks Re[a] Re[a]_err Im[a] Im[a]_err Q_tot Q_tot_err Q_i Q_i_err Q_c Q_c_err nu_r nu_r_err phi_0 phi_0_err reduced_chi2')
        for e in self.entry:
            try:
                print(self.temperature, e['target_freq'], e['channel'], e['depth'], e['is_out_of_res'], e['number_of_peaks'], e['Re[a]'].n, e['Re[a]'].s, e['Im[a]'].n, e['Im[a]'].s, e['Q_tot'].n, e['Q_tot'].s, e['Q_i'].n, e['Q_i'].s, e['Q_c'].n, e['Q_c'].s, e['nu_r'].n, e['nu_r'].s, e['phi_0'].n, e['phi_0'].s, e['reduced_chi2'])
            except:
                print(self.temperature, e['target_freq'], e['channel'], e['depth'], e['is_out_of_res'], e['number_of_peaks'], None, None, None, None, None, None, None, None, None, None, None, None, None, None, None)
    
    def plotTarget(self, flat_at_0db=False):
        from matplotlib import pyplot as plt
        fig,ax = plt.subplots()
        
        artists = []
        
        for e in self.entry:
            # read one sweep at a time
            channel = e['channel']
            x_data_chan = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "freqs.npy")
            y_data_chan = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "mag.npy")
            
            if flat_at_0db:
                y_data_chan -= max(y_data_chan)
            
            if e['is_out_of_res']:
                line, = plt.plot(x_data_chan, y_data_chan, linewidth=1, alpha=0.2, color='black', gid=e['channel'])
            else:
                line, = plt.plot(x_data_chan, y_data_chan, linewidth=1, alpha=1.0, gid=e['channel'])
        
            artists.append(line)
        
        annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        
        def on_plot_hover(event):
            if event.inaxes != ax: #exit if mouse is not on figure
                return
            is_vis = annot.get_visible() #check if an annotation is visible
            # x,y = event.xdata,event.ydata #coordinates of mouse in graph
            for ii, artist in enumerate(artists):
                is_contained, dct = artist.contains(event)
        
                if(is_contained):
                    if('get_data' in dir(artist)): #for plot
                        data = list(zip(*artist.get_data()))
                    elif('get_offsets' in dir(artist)): #for scatter
                        data = artist.get_offsets().data
        
                    #___ Set Annotation settings
                    pos = data[0]
                    annot.xy = pos
                    annot.set_text("{:d}".format(artist.get_gid()))
                    annot.get_bbox_patch().set_edgecolor('black')
                    annot.get_bbox_patch().set_alpha(0.7)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                     if is_vis:
                         annot.set_visible(False) #disable when not hovering
                         fig.canvas.draw_idle()
        
        fig.canvas.mpl_connect("motion_notify_event", on_plot_hover)
        
        plt.grid(linestyle='-', alpha=0.5)
        plt.ylabel('Mag [dB]')
        plt.xlabel('Frequency [MHz]')
        
        plt.show()
        
        
    def plotChannel(self, channel):
        from matplotlib import pyplot as plt
        fig = plt.figure()
        fig.set_size_inches(12, 4.5)
        
        plt.ticklabel_format(axis='both', style='sci', useMathText=True)
        ax0 = plt.subplot(131)
        ax0.set_aspect('equal')
        ax1 = plt.subplot(132)
        ax2 = plt.subplot(133)
        plt.subplots_adjust(bottom=0.15, right=0.98, top=0.95, left=0.09, wspace=0.35)
        
        # read one sweep at a time
        freqs = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "freqs.npy")
        
        I = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "I.npy")
        Q = np.load(paths.target_S21 / self.filename / "{:03d}".format(channel) / "Q.npy")
        
        accumulation_length = 2**20 #for fs=244.14Hz #2**21 for fs=122.07Hz
        fft_len = 1024
        I /= (2**31-1)    # mistral client
        I /= (accumulation_length-1)/(0.5*fft_len)    # mistral client
        Q /= (2**31-1)    # mistral client
        Q /= (accumulation_length-1)/(0.5*fft_len)    # mistral client
        
        # va capito come mettere I e Q in dB
        #I = 20*np.log10(I)
        #Q = 20*np.log10(Q)
        
        amp = np.sqrt(I*I+Q*Q)
        amp = 20*np.log10(amp)
        
        ph = np.arctan2(Q, I)
        ph = np.unwrap(ph, period=np.pi)
        
        ax0.plot(I, Q, linewidth=1, color='k')
        ax1.plot(freqs, amp, linewidth=1, color='k')
        ax2.plot(freqs, ph, linewidth=1, color='k')
    
        ax0.grid(linestyle='-', alpha=0.5)
        ax1.grid(linestyle='-', alpha=0.5)
        ax2.grid(linestyle='-', alpha=0.5)
        ax0.set_ylabel('Q [arbitrary units]')
        ax0.set_xlabel('I [arbitrary units]')
        ax1.set_ylabel('Mag [dB]')
        ax1.set_xlabel('Frequency [MHz]')
        ax2.set_ylabel('Phase [rad]')
        ax2.set_xlabel('Frequency [MHz]')
        
        plt.show()
        
        
    def plotS21(self, channel):
        target_path = paths.target_S21 / self.filename / '{:03d}'.format(channel)
        
        complexS21Plot(target_path)
        

'''
            Othe functions
'''

'''
            readMS2034B function
'''   
class readMS2034B():
    def __init__(self, filename, label=None, mode='real_imag'):
        self.filename = filename
        self.mode = mode
        
        if label != None:
            self.label = label
        else:
            self.label = filename
        
        if self.mode == 'real_imag':
            [self.freqs, self.ReS11, self.ImS11, self.ReS21, self.ImS21, 
             self.ReS12, self.ImS12, self.ReS22, self.ImS22] = np.loadtxt(fname=paths.anritsuMS2034B / filename, dtype=float, skiprows=23, unpack=True)
        
        if self.mode == 'log_mag_phase':
            [self.freqs, self.S11DB, self.S11A, self.S21DB, self.S21A, 
             self.S12DB, self.S12A, self.S22DB, self.S22A] = np.loadtxt(fname=paths.anritsuMS2034B / filename, dtype=float, skiprows=23, unpack=True)
        
        self.freqs *= 1e3 # GHz to MHz
    
    
    def fitS21(self):
        if self.mode != 'real_imag':
            print('S21 fit only available for "real_imag" .s2p file format.')
            return
        
        I = self.ReS21
        Q = self.ImS21
        
        amp = np.sqrt(I*I + Q*Q)
        res_freq = self.freqs[np.argmin(amp)]
        out_path = paths.anritsuMS2034B
        
        np.save(out_path/"I.npy", arr=I)
        np.save(out_path/"Q.npy", arr=Q)
        np.save(out_path/"freqs.npy", arr=self.freqs)
        np.save(out_path/"mag.npy", arr=amp)
        
        try:
            params, chi2 = complexS21Fit(I=I, Q=Q, freqs=self.freqs, res_freq=res_freq, 
                                   output_path=out_path, verbose=True, DATAPOINTS=3500)
        
            if params != None:
                self.Rea = params['Re[a]']
                self.Ima = params['Im[a]']
                self.Qtot = params['Q_tot']
                self.Qc = params['Q_c']
                self.Qi = params['Q_i']
                self.nur = params['nu_r']
                self.phi0 = params['phi_0']
                self.chi2 = chi2
        except:
            pass
        
    
    def plotS21(self):
        if self.mode != 'real_imag':
            print('S21 plot only available for "real_imag" .s2p file format.')
            return
        
        target_path = paths.anritsuMS2034B
        complexS21Plot(target_path)
    
    def plotVNA(self):
        from matplotlib import pyplot as plt
        fig = plt.figure()
        fig.set_size_inches(12, 6)
        ax0 = plt.subplot(221)
        ax1 = plt.subplot(223)
        ax2 = plt.subplot(222)
        
        if self.mode == 'log_mag_phase':
            amp = self.S21DB
            ph = self.S21A
            ph = np.unwrap(ph)
        if self.mode == 'real_imag':
            amp = 20*np.log10(np.sqrt(self.ReS21**2 + self.ImS21**2))
            ph = np.arctan2(self.ImS21, self.ReS21)
            ph = np.unwrap(ph)
        
        ax0.plot(self.freqs, amp, color='black', linewidth=1)
        ax1.plot(self.freqs, ph, color='black', linewidth=1)
        ax2.plot(self.ReS21*1e3, self.ImS21*1e3, color='black', linewidth=1)
        
        ax0.yaxis.set_ticks_position('both')
        ax0.xaxis.set_ticks_position('both')
        ax0.minorticks_on()
        ax0.yaxis.set_tick_params(direction='in', which='both')
        ax0.xaxis.set_tick_params(direction='in', which='both')
        ax0.grid(linestyle='-', alpha=0.5)
        ax0.set_ylabel('Mag [dB]')
        ax0.set_xlabel('Frequency [MHz]')
        
        ax1.yaxis.set_ticks_position('both')
        ax1.xaxis.set_ticks_position('both')
        ax1.minorticks_on()
        ax1.yaxis.set_tick_params(direction='in', which='both')
        ax1.xaxis.set_tick_params(direction='in', which='both')
        ax1.grid(linestyle='-', alpha=0.5)
        ax1.set_ylabel('Phase [rad]')
        ax1.set_xlabel('Frequency [MHz]')
        
        ax2.set_aspect('equal')
        ax2.yaxis.set_ticks_position('both')
        ax2.xaxis.set_ticks_position('both')
        ax2.minorticks_on()
        ax2.yaxis.set_tick_params(direction='in', which='both')
        ax2.xaxis.set_tick_params(direction='in', which='both')
        ax2.grid(linestyle='-', alpha=0.5)
        ax2.set_ylabel(r'Im[S$_{{21}}$] $\times 10^{{-3}}$')
        ax2.set_xlabel(r'Re[S$_{{21}}$] $\times 10^{{-3}}$')
        
        plt.show()


'''
        Function for overplotting target data
'''
def overplotTarget(targets=None, ms2034b_data_list=None, complex_fit_above=False, flat_at_0db=True, colormap='coolwarm'):
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import cm
    cmap = cm.get_cmap(colormap, lut=None)
    
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax0 = plt.subplot(111)
    
    ax0.yaxis.set_ticks_position('both')
    ax0.xaxis.set_ticks_position('both')
    ax0.minorticks_on()
    ax0.yaxis.set_tick_params(direction='in', which='both')
    ax0.xaxis.set_tick_params(direction='in', which='both')
    ax0.grid(linestyle='-', alpha=0.5)
    ax0.set_ylabel('Mag [dB]')
    ax0.set_xlabel('Frequency [MHz]')
    
    # plot roach target sweeps
    handles = []
    if targets != None:
        for i,target in enumerate(targets):
            
            color = cmap(i/(len(targets)-1))
            
            for e in target.entry:
                # read one sweep at a time
                if e['is_out_of_res']==False:
                    channel = e['channel']
                    x_data_chan = np.load(paths.target_S21 / target.filename / "{:03d}".format(channel) / "freqs.npy")
                    y_data_chan = np.load(paths.target_S21 / target.filename / "{:03d}".format(channel) / "mag.npy")
                    
                    if flat_at_0db:
                        y_offset = y_data_chan[0]
                        y_data_chan -= y_offset
                
                    ax0.plot(x_data_chan, y_data_chan, color=color, linestyle='', marker='o', markersize=5)
                    
                    if complex_fit_above:
                        nu_linear = np.linspace(x_data_chan[0], x_data_chan[-1], num=200) # linear sample
                        nu_peack = np.random.normal(e['nu_r'].n, 0.001, 1000) # peak sample
                        nu = np.concatenate([nu_linear, nu_peack])
                        nu = np.sort(nu)
                        Z = S_21(nu, e['Re[a]'].n, e['Im[a]'].n, e['Q_tot'].n, e['Q_c'].n, e['nu_r'].n, e['phi_0'].n, tau=0.04)
                        mag = 20*np.log10(np.abs(Z))
                        if flat_at_0db:
                            mag -= y_offset
                        ax0.plot(nu, mag, linestyle='solid', color=color, alpha=0.6, linewidth=4)
                
            handles.append(Line2D([0], [0], label=target.label, color=color))
            
    # plot ms2034b VNA sweeps
    if ms2034b_data_list != None:
        for i,sweep in enumerate(ms2034b_data_list):
            if sweep.mode == 'log_mag_phase':
                amp = sweep.S21DB
                ph = sweep.S21A
                ph = np.unwrap(ph)
            if sweep.mode == 'real_imag':
                amp = 20*np.log10(np.sqrt(sweep.ReS21**2 + sweep.ImS21**2))
                ph = np.arctan2(sweep.ImS21, sweep.ReS21)
                ph = np.unwrap(ph)
            
            if flat_at_0db:
                amp -= amp[0]
            
            color = cmap(i/(len(ms2034b_data_list)-1))
            ax0.plot(sweep.freqs, amp, color=color, linestyle='dotted', linewidth=1)
            handles.append(Line2D([0], [0], linestyle='dotted', label=sweep.label, color=color))
    
    ax0.legend(loc='best', handles=handles, fontsize=8)
    
    plt.show()



def overplotTargetCircles(targets=None, ms2034b_data_list=None, complex_fit_above=False, colormap='coolwarm'):
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import cm
    cmap = cm.get_cmap(colormap, lut=None)
    
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax0 = plt.subplot(111)
    
    ax0.yaxis.set_ticks_position('both')
    ax0.xaxis.set_ticks_position('both')
    ax0.minorticks_on()
    ax0.yaxis.set_tick_params(direction='in', which='both')
    ax0.xaxis.set_tick_params(direction='in', which='both')
    ax0.grid(linestyle='-', alpha=0.5)
    ax0.set_ylabel('Q')
    ax0.set_xlabel('I')
    ax0.set_aspect('equal')
    
    # plot roach target sweeps
    handles = []
    if targets != None:
        for i,target in enumerate(targets):
            
            color = cmap(i/(len(targets)-1))
            
            for e in target.entry:
                # read one sweep at a time
                if e['is_out_of_res']==False:
                    channel = e['channel']
                    x_data_chan = np.load(paths.target_S21 / target.filename / "{:03d}".format(channel) / "freqs.npy")
                    
                    I = np.load(paths.target_S21 / target.filename / "{:03d}".format(channel) / "I_prime.npy")
                    Q = np.load(paths.target_S21 / target.filename / "{:03d}".format(channel) / "Q_prime.npy")
                
                    ax0.plot(I, Q, color=color, linestyle='', marker='o', markersize=5)
                
                    if complex_fit_above:
                        nu_linear = np.linspace(x_data_chan[0], x_data_chan[-1], num=200) # linear sample
                        nu_peack = np.random.normal(e['nu_r'].n, 0.001, 1000) # peak sample
                        nu = np.concatenate([nu_linear, nu_peack])
                        nu = np.sort(nu)
                        Z = S_21(nu, e['Re[a]'].n, e['Im[a]'].n, e['Q_tot'].n, e['Q_c'].n, e['nu_r'].n, e['phi_0'].n, tau=0.04)
                        
                        ax0.plot(np.real(Z), np.imag(Z), linestyle='solid', color=color, alpha=0.6, linewidth=4)
                
            handles.append(Line2D([0], [0], label=target.label, color=color))
            
    # plot ms2034b VNA sweeps
    if ms2034b_data_list != None:
        for i,sweep in enumerate(ms2034b_data_list):
            amp = 20*np.log10(np.sqrt(sweep.ReS21**2.0 + sweep.ImS21**2.0))
            
            color = cmap(i/(len(ms2034b_data_list)-1))
            ax0.plot(sweep.freqs, amp, color=color, linestyle='dotted', linewidth=1)
            handles.append(Line2D([0], [0], linestyle='dotted', label=sweep.label, color=color))
    
    ax0.legend(loc='best', handles=handles, fontsize=8)
    
    plt.show()


        
'''
            lowpass_cosine function
'''
def lowpass_cosine(y, tau, f_3db, width, padd_data=True):
    
    import numpy as nm
        # padd_data = True means we are going to symmetric copies of the data to the start and stop
    # to reduce/eliminate the discontinuities at the start and stop of a dataset due to filtering
    #
    # False means we're going to have transients at the start and stop of the data

    # kill the last data point if y has an odd length
    if nm.mod(len(y),2):
        y = y[0:-1]

    # add the weird padd
    # so, make a backwards copy of the data, then the data, then another backwards copy of the data
    if padd_data:
        y = nm.append( nm.append(nm.flipud(y),y) , nm.flipud(y) )

    # take the FFT
        import scipy
        import scipy.fftpack
    ffty=scipy.fftpack.fft(y)
    ffty=scipy.fftpack.fftshift(ffty)

    # make the companion frequency array
    delta = 1.0/(len(y)*tau)
    nyquist = 1.0/(2.0*tau)
    freq = nm.arange(-nyquist,nyquist,delta)
    # turn this into a positive frequency array
    pos_freq = freq[(len(ffty)//2):]

    # make the transfer function for the first half of the data
    i_f_3db = min( nm.where(pos_freq >= f_3db)[0] )
    f_min = f_3db - (width/2.0)
    i_f_min = min( nm.where(pos_freq >= f_min)[0] )
    f_max = f_3db + (width/2);
    i_f_max = min( nm.where(pos_freq >= f_max)[0] )

    transfer_function = nm.zeros(int(len(y)//2))
    transfer_function[0:i_f_min] = 1
    transfer_function[i_f_min:i_f_max] = (1 + nm.sin(-nm.pi * ((freq[i_f_min:i_f_max] - freq[i_f_3db])/width)))/2.0
    transfer_function[i_f_max:(len(freq)//2)] = 0

    # symmetrize this to be [0 0 0 ... .8 .9 1 1 1 1 1 1 1 1 .9 .8 ... 0 0 0] to match the FFT
    transfer_function = nm.append(nm.flipud(transfer_function),transfer_function)

    # apply the filter, undo the fft shift, and invert the fft
    filtered=nm.real(scipy.fftpack.ifft(scipy.fftpack.ifftshift(ffty*transfer_function)))

    # remove the padd, if we applied it
    if padd_data:
        filtered = filtered[(len(y)//3):(2*(len(y)//3))]

    # return the filtered data
    return filtered

'''
            asymmetric_least_squares_smoothing function
'''
def asymmetric_least_squares_smoothing(data, lam, p, N_iter):
    '''
    lam: adjusting parameter
    p: asymmetry parameter
    N_iter: number of iteration
    '''
    from scipy.sparse import spdiags, linalg, diags
    L = len(data)
    D = diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np.ones(L)
    for i in range(N_iter):
        W = spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = linalg.spsolve(Z, w*data)
        w = p*(data < z) + (1.0-p)*(data >= z)
    return z

'''
            adaptive_iteratively_reweighted_penalized_least_squares_smoothing function
'''
def adaptive_iteratively_reweighted_penalized_least_squares_smoothing(data, lam, N_iter):
    '''
    lam: adjusting parameter
    N_iter: number of iteration
    '''
    from scipy.sparse import spdiags, linalg, diags
    from scipy.linalg import norm
    L = len(data)
    D = diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np.ones(L)
    for i in range(N_iter):
        W = spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = linalg.spsolve(Z, w*data)
        d_mod = norm((z-data)[z>data])
        if d_mod < 0.001 * norm(data):
            return z
        p = np.exp(i*(data-z)/d_mod)
        w = 0.0*(data < z) + p*(data >= z)
    return z
        
        
'''
            buildS21Dataset function
'''      
def buildS21Dataset(sweep, ROACH='MISTRAL'):
    '''
    This function converts the raw data sweep into different .npy files
    each for a single channel

    Parameters
    ----------
    sweep : Target object or VNA object
        A Target object or VNA object to be converted.
    ROACH : string, optional
        If 'MISTRAL' the tone frequencies are computed as freqs = LO + bb,
        if 'OLIMPO' the tone frequencies are computed as freqs = LO/2 + bb. 
        The default is 'MISTRAL'.

    Returns
    -------
    None.

    '''
    from  tqdm import tqdm
    
    filename = sweep.filename
    
    if  (paths.target / filename).exists():
        sweep_path = paths.target / filename
        output_path = paths.target_S21 / filename
    elif (paths.vna / filename).exists():
        sweep_path = paths.vna / filename
        output_path = paths.vna_S21 / filename
    else:
        print("Sweep file not found.")
        sys.exit()
    
    # make the S21 directory
    if not output_path.exists():
        output_path.mkdir()
        
    print("\nBuilding S21 dataset from raw data sweep...")
    
    if ROACH == "OLIMPO":
        print("\tOLIMPO ROACH selected (freqs = LO/2 + bb)")
    elif ROACH == "MISTRAL":
        print("\tMISTRAL ROACH selected (freqs = LO + bb)")
    
    LO_freqs = np.load(sweep_path / "sweep_freqs.npy")
    bb_freqs = np.load(sweep_path / "bb_freqs.npy")
    
    n_res = bb_freqs.size
    n_files = LO_freqs.size
    
    I = np.zeros([n_res, n_files])
    Q = np.zeros([n_res, n_files])
    
    pbar = tqdm(LO_freqs, position=0, leave=True)
    for i_file, LO_freq in enumerate(pbar):
        pbar.set_description("\tReading LO data... ")
        I_all = np.load(sweep_path / ("I"+str(LO_freq)+".npy"))
        Q_all = np.load(sweep_path / ("Q"+str(LO_freq)+".npy"))
    
        for i_res in range(n_res):
            I[i_res][i_file] = I_all[i_res]
            Q[i_res][i_file] = Q_all[i_res]
    
    # data from mistral client to compute amplitude in dBm
    accumulation_length = 2**21
    fft_length = 1024
    pbar = tqdm(range(n_res), position=0, leave=True)
    for i in pbar:
        pbar.set_description("\tComputing frequencies... ")
        
        Q[i] *= (0.5*fft_length) / ((2**31-1)*(accumulation_length-1))
        I[i] *= (0.5*fft_length) / ((2**31-1)*(accumulation_length-1))
        
        mag = np.sqrt(Q[i]*Q[i] + I[i]*I[i])
        #mag /= (2**31-1)    # mistral client
        #mag /= (accumulation_length-1)/(0.5*fft_length)    # mistral client
        mag = 20*np.log10(mag)
        phase = np.arctan2(Q[i], I[i])
    
        if ROACH == "OLIMPO":
            freqs = LO_freqs*0.5+bb_freqs[i]
        elif ROACH == "MISTRAL":
            freqs = LO_freqs+bb_freqs[i]
    
        if not (output_path / "{:03d}".format(i)).exists():
            (output_path / "{:03d}".format(i)).mkdir()
        np.save(output_path / "{:03d}".format(i) / "I.npy", I[i])
        np.save(output_path / "{:03d}".format(i) / "Q.npy", Q[i])
        np.save(output_path / "{:03d}".format(i) / "mag.npy", mag)
        np.save(output_path / "{:03d}".format(i) / "phase.npy", phase)
        np.save(output_path / "{:03d}".format(i) / "freqs.npy", freqs/1e6)




'''
            S_21 function
'''
def S_21(nu, Rea, Ima, Q_tot, Q_c, nu_r, phi_0, tau):
    '''
    This function returns the S21 scattering parameter.

    Parameters
    ----------
    nu : number, usually is a numpy array or listo of number
        Frequency in [Hz].
    Rea : number
        Real part of the amplitude parameter.
    Ima : number
        Imaginary part of the amplitude parameter.
    Q_tot : number
        Total quality factor.
    Q_c : number
        Coupling quality factor.
    nu_r : number
        Resonant frequency in [Hz].
    phi_0 : number, between [0; 2pi]
        Phase parameter decribing the rotation of the resonant circle.
    tau : number
        Time delay due to the transmission line length in [s].

    Returns
    -------
    comlex number or comlex numpy array
        The computed S21 scattering parameter.

    '''
    a = Rea + Ima*1j
    return np.exp(1j*2.0*np.pi*tau*nu)*a*(1.0 - (Q_tot/Q_c)*np.exp(1j*phi_0) / (1+2*1j*Q_tot*(nu-nu_r)/nu_r) )



'''
            jointTargetSweeps function
'''
def jointTargetSweeps(targets, exclude_channels=[[]], flat_at_0db=False):
    '''
    This function generates a plot with jointed target sweeps
    
    Parameters
    ----------
    targets : list of target objecs
        A list of target objects to be jointed.
    exclude_channels : list of list with channel indexes, optional
        Channel lists to be excluded in the plot. The default is [[]].
    flat_at_0db : boolean, optional
        If True, each resonance will be vertically shifted such that 
        its maximum is 0dB. The default is False.

    Returns
    -------
    None.

    '''
    from matplotlib import pyplot as plt
    fig = plt.figure()
    fig.set_size_inches(6, 4)
    
    for t,exclude_channels in zip(targets,exclude_channels):
        # how many channels?
        n_chan = list(range(t.target_freqs_new.size))
        for c in exclude_channels:
            n_chan.remove(c)
        
        for chan in n_chan:
            # read one sweep at a time
            x_data_chan = np.load(paths.target_S21 / t.filename / "{:03d}".format(chan) / "freqs.npy")
            y_data_chan = np.load(paths.target_S21 / t.filename / "{:03d}".format(chan) / "mag.npy")
            
            if flat_at_0db == True:
                y_data_chan -= max(y_data_chan)
            
            plt.plot(x_data_chan, y_data_chan, linewidth=1)
        
    plt.grid(linestyle='-', alpha=0.5)
    plt.ylabel('Mag [dB]')
    plt.xlabel('Frequency [MHz]')
    
    plt.show()
    #fig.savefig('vna.png', dpi=300)
        


'''
                complexS21Fit function
'''
def complexS21Fit(I, Q, freqs, res_freq, output_path, DATAPOINTS=100, verbose=False):
    from lmfit import Parameters, minimize, fit_report
    from uncertainties import ufloat
    
    # DATAPOINTS: number of left and right datapoints with respect to the resonance frequency 
    # to which the complex fit is computed
    
    from pathlib import Path
    output_path = Path(output_path)
    
    if not output_path.exists():
        print("\n{:s} does not exists.".format(output_path.as_posix()))
        return 0
    
    if verbose:
        print("\nFit S21 complex function for {:s}...".format(output_path.as_posix()))
    
    A = np.sqrt(I**2 + Q**2)
    phase = np.arctan2(Q, I)
    
    ARG_RESFREQ = np.argmin(A)
    ARG_MIN = ARG_RESFREQ-DATAPOINTS
    ARG_MAX = ARG_RESFREQ+DATAPOINTS
    # check if there are enough datapoints at the left of the resonance frequency
    if ARG_MIN < 0:
        ARG_MIN = 0
    # check if there are enough datapoints at the right of the resonance frequency
    if ARG_MAX >= A.size:
        ARG_MAX = A.size-1
        
    
    # removing the cable delay
    tau = 0.04 # microsec
    phase += 2.0*np.pi*tau*freqs
    phase -= int(phase[0]/np.pi)*np.pi
    phase = np.unwrap(phase)
    I = A*np.cos(phase)
    Q = A*np.sin(phase)
    
    # errors on I and Q data is 1%
    IErr = I*0.01
    QErr = Q*0.01
    
    # compute the center coordinates by averaging max and min data values
    I_m, I_M = I.min(), I.max()
    Q_m, Q_M = Q.min(), Q.max()
    Xc = (I_m+I_M)*0.5
    Yc = (Q_m+Q_M)*0.5
    
    # compute the radius by averaging the 'x radius' and the 'y radius'
    _rx = 0.5*np.abs(I_M-I_m)
    _ry = 0.5*np.abs(Q_M-Q_m)
    radius = 0.5*(_rx+_ry)
    
    if verbose:
        print("\n\tcenterd circle parameters")
        print("\tXc = ", Xc)
        print("\tYc = ", Yc)
        print("\tradius = ", radius)
    
    # translation
    I_rot = I-Xc
    Q_rot = Q-Yc
    
    # rotation
    I_mid = 0.5*(I_rot[-1]+I_rot[0])
    Q_mid = 0.5*(Q_rot[-1]+Q_rot[0])
    rotAngle = np.pi-np.arctan2(Q_mid, I_mid)
    
    if verbose:
        print("\n\trotation angle")
        print("\tangle = ", rotAngle)
    
    file = open(output_path / "log.txt","w")
    file.write("centerd circle parameters\n")
    file.write("Xc = {:.3e}\n".format(Xc))
    file.write("Yc = {:.3e}\n".format(Yc))
    file.write("radius = {:.3e}\n".format(radius))
    file.write("\n")
    file.write("rotation angle\n")
    file.write("angle = {:.3e}\n".format(rotAngle))
    file.write("\n")
    file.close()
    
    I_second = np.cos(rotAngle)*I_rot-np.sin(rotAngle)*Q_rot
    Q_second = np.sin(rotAngle)*I_rot+np.cos(rotAngle)*Q_rot
    
    
    # determine the phase and amplitude (prime)
    ph_second = np.arctan2(Q_second, I_second)
    mag_second = np.sqrt(I_second**2 + Q_second**2)
    
    
    # phase fit
    def phaseFunction(nu, Q_tot, nu_r):
        return -2*np.arctan2((2.0*Q_tot*(nu/nu_r-1.0)), 1.0)
    
    def phaseFunctionResiduals(params, nu, data, uncertainty):
        Q_tot = params['Q_tot']
        nu_r = params['nu_r']
        return (data-phaseFunction(nu, Q_tot, nu_r))/uncertainty
        
    params = Parameters()
    params.add('Q_tot', value=5000, min=0)
    params.add('nu_r', value=res_freq, min=res_freq*0.99, max=res_freq*1.01)
    
    x_data = freqs[ARG_MIN:ARG_MAX]
    y_data = ph_second[ARG_MIN:ARG_MAX]
    uncertainty = 0.01*ph_second[ARG_MIN:ARG_MAX]
    
    out = minimize(phaseFunctionResiduals, params, args=(x_data, y_data, uncertainty), method='leastsq')
    
    if verbose:
        print("\n\tphase fit")
        print("\t"+fit_report(out))
    file = open(output_path / "log.txt","a")
    file.write(fit_report(out))
    file.close()
    
    Q_totPhFit = out.params['Q_tot'].value
    nu_rPhFit = out.params['nu_r'].value
    
    # compute Qc and Phi0 from fit parameters
    # guess the value of a
    Ima = 0.5*(Q[-1]+Q[0])
    Rea = 0.5*(I[-1]+I[0])
    a_norm = np.sqrt(Ima**2.0 + Rea**2.0)
    Q_c = 0.5*a_norm*Q_totPhFit/radius
    # let's compute phi0
    alpha_angle = np.angle(Rea+1j*Ima)
    beta_angle = np.angle(Xc-Rea+1j*(Yc-Ima))
    phi_0 = np.pi+beta_angle-alpha_angle
    if phi_0<0.0:
        phi_0 += 2.0*np.pi
    if phi_0>2.0*np.pi:
        phi_0 %= (2.0*np.pi)
    
    # 6-parameters complex fit (fine tuning)
    # complex residuals over errorbars
    def complexResiduals(params, freq, z_val, z_err):
        Rea = params['Rea']
        Ima = params['Ima']
        Q_tot = params['Q_tot']
        Q_c = params['Q_c']
        nu_r = params['nu_r']
        phi_0 = params['phi_0']
        tau = params['tau']
        return np.abs( (z_val - S_21(freq, Rea, Ima, Q_tot, Q_c, nu_r, phi_0, tau))/z_err )
    
    # chi2 minimization
    params = Parameters()
    params.add('Rea', value=Rea)
    params.add('Ima', value=Ima)
    params.add('Q_tot', value=Q_totPhFit, min=0)
    params.add('Q_c', value=Q_c, min=0)
    params.add('nu_r', value=nu_rPhFit, min=nu_rPhFit*0.99, max=nu_rPhFit*1.01)
    params.add('phi_0', value=phi_0, min=0.0, max=2.0*np.pi)
    params.add('tau', value=tau, vary=False)
    
    z_data = I[ARG_MIN:ARG_MAX]+1j*Q[ARG_MIN:ARG_MAX]
    z_err = IErr[ARG_MIN:ARG_MAX]+1j*QErr[ARG_MIN:ARG_MAX]
    freqs = freqs[ARG_MIN:ARG_MAX]
    
    
    out = minimize(complexResiduals, params, args=(freqs, z_data, z_err), method='leastsq')
    
    if verbose:
        print("\n\tcomplex fit results")
        print("\t"+fit_report(out))
    file = open(output_path / "log.txt","a")
    file.write(fit_report(out))
    file.close()
    
    if out.message != 'Fit succeeded.':
        return None
    
    Rea = ufloat(out.params['Rea'].value, out.params['Rea'].stderr)
    Ima = ufloat(out.params['Ima'].value, out.params['Ima'].stderr)
    Q_tot = ufloat(out.params['Q_tot'].value, out.params['Q_tot'].stderr)
    Q_c = ufloat(out.params['Q_c'].value, out.params['Q_c'].stderr)
    nu_r = ufloat(out.params['nu_r'].value, out.params['nu_r'].stderr)
    phi_0 = ufloat(out.params['phi_0'].value, out.params['phi_0'].stderr)

    # compute the Q_i value
    Q_i = Q_c*Q_tot/(Q_c-Q_tot)
    
    if verbose:
        print("\nQ_i: {:uf}".format(Q_i))
    file = open(output_path / "log.txt","a")
    file.write("\nQ_i: {:uf}".format(Q_i))
    file.close()

    # save data into .npy files
    np.save(output_path / "I_prime.npy", I)
    np.save(output_path / "Q_prime.npy", Q)
    np.save(output_path / "phase_prime.npy", phase)
    np.save(output_path / "I_second.npy", I_second)
    np.save(output_path / "Q_second.npy", Q_second)
    np.save(output_path / "mag_second.npy", mag_second)
    np.save(output_path / "phase_second.npy", ph_second)
    np.save(output_path / "transformation_parameters.npy", arr=[Xc, Yc, rotAngle])
    np.save(output_path / "phase_fit_parameters.npy", arr=[radius, Q_totPhFit, nu_rPhFit])
    np.save(output_path / "reduced_chi2.npy", arr=out.redchi)
    np.save(output_path / "complex_parameters.npy", arr=[Rea.n, Rea.s, Ima.n, Ima.s, Q_tot.n, Q_tot.s, Q_c.n, Q_c.s, 
                                                         Q_i.n, Q_i.s, nu_r.n, nu_r.s, phi_0.n, phi_0.s, tau])
    
    return {'Re[a]': Rea, 'Im[a]': Ima, 'Q_tot': Q_tot, 'Q_c': Q_c, 'Q_i': Q_i, 'nu_r': nu_r, 'phi_0': phi_0, 'tau': tau}, out.redchi


'''
            complexS21Plot function
'''
def complexS21Plot(complex_fit_data_path):
    from pathlib import Path
    complex_fit_data_path = Path(complex_fit_data_path)
    
    I = np.load(complex_fit_data_path / "I.npy")
    Q = np.load(complex_fit_data_path / "Q.npy")
    mag = np.load(complex_fit_data_path / "mag.npy")
    #ph = np.load(output_path / "phi.npy")
    freqs = np.load(complex_fit_data_path / "freqs.npy")
        
    A = np.sqrt(I**2 + Q**2)
    
    #IErr = I*0.01
    #QErr = Q*0.01
    
    I_prime = np.load(complex_fit_data_path / "I_prime.npy")
    Q_prime = np.load(complex_fit_data_path / "Q_prime.npy")
    phase_prime = np.load(complex_fit_data_path / "phase_prime.npy")
    phase_prime = np.unwrap(phase_prime)
    
    I_second = np.load(complex_fit_data_path / "I_second.npy")
    Q_second = np.load(complex_fit_data_path / "Q_second.npy")
    #mag_second = np.load(output_path / "mag_second.npy")
    phase_second = np.load(complex_fit_data_path / "phase_second.npy")
    if phase_second[0] <= 0.0:
        phase_second += 2.0*np.pi
    phase_second = np.unwrap(phase_second)
    
    [Xc, Yc, rotAngle] = np.load(complex_fit_data_path / "transformation_parameters.npy")
    [radius, Q_phfit, nu_r_ph_fit] = np.load(complex_fit_data_path / "phase_fit_parameters.npy")
    [Rea, ReaErr, Ima, ImaErr, Qt, QtErr, Qc, QcErr, Qi, QiErr, nu_r, nu_rErr, phi0, phi0Err, tau] = np.load(complex_fit_data_path / "complex_parameters.npy")
    
    def phaseFunction(nu, Q, nu_r):
        return -2*np.arctan2((2.0*Q*(nu/nu_r-1.0)), 1.0)
    
    # PLOT
    from matplotlib import pyplot as plt
    # create a figure
    # text settings
    fig = plt.figure()
    fig.set_size_inches(15, 7)
    plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.07)
    plt.subplots_adjust(hspace=0.0)
    circlePlot = plt.subplot(121, aspect='equal')
    phasePlot1 = plt.subplot(322)
    phasePlot2 = plt.subplot(324)
    amplitudePlot = plt.subplot(326)
    
    color_raw = (17/255, 157/255, 164/255, 1.0)
    color_raw_alpha = (17/255, 157/255, 164/255, 0.3)
    color_notau = (0, 0, 0, 1.0)
    color_notau_alpha = (0, 0, 0, 0.3)
    color_centered = (0, 210/255, 249/255, 1.0)
    color_centered_alpha = (0, 210/255, 249/255, 0.3)
    
    # raw IQ data plot
    circlePlot.plot(I, Q, color=color_raw,linestyle='-', marker='o', markersize=4, markerfacecolor=color_raw_alpha)#, label='Raw data')
    circlePlot.plot(I_prime, Q_prime, color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)#, label=r'$\tau$ removed')
    #circlePlot.errorbar(I, Q, xerr=IErr, yerr=QErr, marker='.', linestyle='--', linewidth=1.0, markersize=1.0, color='green', alpha=1.0, label='Raw data')
    #circlePlot.errorbar(I_prime, Q_prime, xerr=IErr, yerr=QErr, marker='.', linestyle='--', linewidth=1.0, markersize=1.0, color='black', alpha=1.0, label='Tau removed')
    
    # Complex plot 
    Z_freqs = np.linspace(freqs[0], freqs[-1], num=2000)
    Z = S_21(Z_freqs, Rea, Ima, Qt, Qc, nu_r, phi0, tau)
    circlePlot.plot(np.real(Z), np.imag(Z), linestyle='-', color='red', alpha=0.5, linewidth=3.0)#, label='S$_{21}$ Fit')
    
    # resonance after translation and rotation
    #circlePlot.plot(I_second, Q_second, marker='.', markersize=1.0, color='black', alpha=0.5, label='Centered data')
    circlePlot.plot(I_second, Q_second, color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha)#, label='Centered data')
    # circle fit plot
    #phiTrain = np.linspace(0.0, 2.0*np.pi, 100)
    #circlePlot.plot(radius*np.cos(phiTrain), radius*np.sin(phiTrain), linestyle='-', color='blue', alpha=0.3, linewidth=3.0)#, label='Circle Fit')
    circlePlot.axvline(x=0.0, ymin=-10.0, ymax=10.0, linewidth='1.0', linestyle='--', color='black', alpha=0.3)
    circlePlot.axhline(y=0.0, xmin=-10.0, xmax=10.0, linewidth='1.0', linestyle='--', color='black', alpha=0.3)
    
    # phase plot
    phasePlot1.plot(freqs, phase_second, color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha)#, label='Centered data')
    #phasePlot1.plot(freqs, ph_second, marker='.', markersize=1.0, color=color_centered, alpha=0.5, label='Phase')
    phasePlot1.plot(freqs, phaseFunction(freqs, Q_phfit, nu_r_ph_fit), linestyle='-', color='blue', alpha=0.3, linewidth=3.0)#, label='Phase fit')
    
    #phasePlot2.plot(freqs, ph_prime, marker='.', markersize=1.0, color=color_notau, alpha=0.5, label='Phase')
    phasePlot2.plot(freqs, phase_prime, color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)#, label='Phase')
    phasePlot2.plot(Z_freqs, np.unwrap(np.angle(Z)), linestyle='-', color='red', alpha=0.5, linewidth=3.0)#, label='S$_{21}$ Fit')
    
    # amplitude plot
    #amplitudePlot.plot(freqs, A, marker='.', markersize=1.0, color=color_notau, alpha=0.5, label='Amplitude')
    amplitudePlot.plot(freqs, mag-mag[0], color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)#, label='Amplitude')
    amplitudePlot.plot(Z_freqs, 20*np.log10(np.abs(Z))-mag[0], linestyle='-', color='red', alpha=0.5, linewidth=3.0)#, label='S$_{21}$ Fit')
    
    handle = [Line2D([0], [0], color=color_raw,linestyle='-', marker='o', markersize=4, markerfacecolor=color_raw_alpha, label='Raw data'),
              Line2D([0], [0], color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha, label=r'$\tau$ removed'),
              Line2D([0], [0], color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha, label='Centered data'),
              #Line2D([0], [0], linestyle='-', color='blue', alpha=0.3, linewidth=3.0, label='Circle Fit'),
              Line2D([0], [0], linestyle='-', color='blue', alpha=0.3, linewidth=3.0, label='Phase fit'),
              Line2D([0], [0], linestyle='-', color='red', alpha=0.5, linewidth=3.0, label='S$_{21}$ Fit')]
    
    # graphics
    circlePlot.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
    circlePlot.set_xlabel('$I$')
    circlePlot.set_ylabel('$Q$')
    #circlePlot.legend(loc='best')
    circlePlot.grid(alpha=0.5)
    circlePlot.yaxis.set_ticks_position('both')
    circlePlot.xaxis.set_ticks_position('both')
    circlePlot.minorticks_on()
    circlePlot.yaxis.set_tick_params(direction='in', which='both')
    circlePlot.xaxis.set_tick_params(direction='in', which='both')
    
    phasePlot1.legend(loc='upper right', handles=handle)
    phasePlot1.grid(axis="x")
    phasePlot1.set_ylabel('Phase [rad]')
    phasePlot1.set_xlabel('Frequency [MHz]')
    phasePlot1.grid(alpha=0.5)
    phasePlot1.yaxis.set_ticks_position('both')
    phasePlot1.xaxis.set_ticks_position('both')
    phasePlot1.minorticks_on()
    phasePlot1.yaxis.set_tick_params(direction='in', which='both')
    phasePlot1.xaxis.set_tick_params(direction='in', which='both')
    #plt.setp(phasePlot1.get_xticklabels(), visible=False)
    
    #phasePlot2.legend(loc='best')
    phasePlot2.grid(axis="x")
    phasePlot2.set_ylabel('Phase [rad]')
    phasePlot2.set_xlabel('Frequency [MHz]')
    phasePlot2.grid(alpha=0.5)
    phasePlot2.yaxis.set_ticks_position('both')
    phasePlot2.xaxis.set_ticks_position('both')
    phasePlot2.minorticks_on()
    phasePlot2.yaxis.set_tick_params(direction='in', which='both')
    phasePlot2.xaxis.set_tick_params(direction='in', which='both')
    #plt.setp(phasePlot2.get_xticklabels(), visible=False)
    
    #amplitudePlot.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #amplitudePlot.legend(loc='best')
    amplitudePlot.set_xlabel('Frequency [MHz]')
    amplitudePlot.set_ylabel(r'Amplitude [dBm]')
    amplitudePlot.grid(alpha=0.5)
    amplitudePlot.yaxis.set_ticks_position('both')
    amplitudePlot.xaxis.set_ticks_position('both')
    amplitudePlot.minorticks_on()
    amplitudePlot.yaxis.set_tick_params(direction='in', which='both')
    amplitudePlot.xaxis.set_tick_params(direction='in', which='both')
    
    # print the figures
    plt.show()



'''
            Plot picolog data
'''
def plotPicolog(filename, OFFSET_TIME=False):
    '''
    This function returns a plot of a given picolog file with three channels

    Parameters
    ----------
    filename : string
        The picolog filename.
    OFFSET_TIME : boolean, optional
        If True, ch1 (usually time) will start with a zero. The default is 
        False.

    Returns
    -------
    None.

    '''
    import matplotlib.pyplot as plt
    from G31_thermometry import Thermometer
    thermometer = Thermometer(model='DT670', serial_no='D6068043')
    
    # read data
    ch1,ch2,ch3 = np.loadtxt(paths.picolog/filename, dtype=float, skiprows=1, unpack=True, delimiter=',')
    
    if OFFSET_TIME:
        ch1 -= ch1[0]
    
    # convert thermometer voltage to temperature
    temperature = thermometer.temperature(ch2/1000.)
    
    fig = plt.figure(figsize=(9,9))
    ax0 = fig.add_subplot(111)
    ax0.plot(ch1, temperature, c='black')
    ax0.set_ylabel('T$^{Cu}$ [K]', color='black')
    ax0.tick_params(axis='y', colors='black')
    ax0.set_xlim([ch1[0], ch1[-1]])
    ax0.xaxis.grid()
    
    ax0_ = fig.add_subplot(111, sharex = ax0, frameon = False)
    ax0_.fill_between(ch1, ch3, 0.0, alpha=0.20, color='red')
    ax0_.yaxis.set_label_position("right")
    ax0_.tick_params(axis='y', colors='red')
    ax0_.yaxis.tick_right()
    delta_ch3 = ch3.max() - ch3.min()
    ax0_.set_ylim([ch3.min()-(delta_ch3*0.2), ch3.max()+(delta_ch3*0.2)])
    ax0_.set_ylabel('Heater [V]', color='red')
    
    plt.show()

def Delta(T_c):
    from scipy.constants import k as kb
    return 1.764*kb*T_c

def n_qp(T, T_c, N_0):
    from scipy.constants import k as kb
    return 2.0*N_0*np.sqrt(2.0*np.pi*kb*T*Delta(T_c)) * np.exp(-Delta(T_c)/(kb*T))

def electrical_phase_responsivity_linear_fit(nu_r, base_nu_r, T, T_c, N_0, V_abs, label='pixel', color='black', axis=None):
    from lmfit import Minimizer, Parameters
    from uncertainties import ufloat
    
    # converts lists to np arrays
    delta_x = [(n-base_nu_r)/base_nu_r for n in nu_r]
    error_max = max([d.s for d in delta_x])
    delta_x = [(n-base_nu_r)/base_nu_r + 0.5*ufloat(0.0, error_max) for n in nu_r]
    
    T = np.asarray(T)
    
    # function used for the fit procedure
    def linear_function(x, m, q):
        return x*m + q
    
    def fcn2min(params, x, data, errs=None):
        m = params['m']
        q = params['q']
        
        model = linear_function(x, m, q)
        
        if errs is None:
            return data-model
        
        return (data-model)/errs
    
    # quasiparticle number density
    N_qp = V_abs * n_qp(T, T_c, N_0)
    
    # linear fit
    params = Parameters()
    params.add('m', value=(min(delta_x).n-max(delta_x).n)/(max(N_qp)-min(N_qp)), min=0, max=-np.inf)
    params.add('q', value=0.0, min=-1, max=1)
    
    try:
        result = Minimizer(fcn2min, params, fcn_args=(N_qp, [d.n for d in delta_x], [d.s for d in delta_x])).minimize(method='least_squares')
        if axis is not None:
            axis.ticklabel_format(axis='both', style='sci', useMathText=True, scilimits=(0,0))
            
            axis.set_xlabel("Quasiparticle number")
            axis.set_ylabel("Relative resonant frequency shift")
            
            axis.errorbar(N_qp, [d.n for d in delta_x], yerr=[d.s for d in delta_x], color=color, linestyle='', fmt='o', capsize=2, markersize=3)
            axis.plot(N_qp, linear_function(N_qp, result.params['m'].value, result.params['q'].value), color=color, label=label)
            
            axis.grid(color='gray', alpha=0.4)
            
            return result.params
    except:
        print("Cannot find fit parameters.")
        return
        
    return result.params

'''
noise_type = 'chp' or 'cha'
'''
def plot_extracted_stream(filename, channel, axis, noise_type='chp', label=None, color=None, linestyle='solid', linewidth=1):
    path = paths.output_noise / filename

    if label is None:
        label = noise_type+'_{:03d}'.format(channel)

    # open time streams
    time_stream = np.load(path / 'time_stream.npy')
    ch = np.load(path / (noise_type+'_{:03d}.npy'.format(channel)) )

    axis.plot(time_stream, ch, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    axis.grid(color="gray", alpha=0.5)
    axis.set_xlabel("Time [s]")
    axis.set_ylabel(noise_type)


'''
noise_type = 'chp' or 'cha'
'''
def plot_extracted_noise_spectral_density(filename, channel, sampling_frequency, axis, noise_type='chp', label=None, color=None, linestyle='solid', linewidth=1):
    from scipy.signal import periodogram
    
    path = paths.output_noise / filename
    
    if label is None:
        label = noise_type+'_{:03d}'.format(channel)

    # open time streams
    time_stream = np.load(path / 'time_stream.npy')
    ch = np.load(path / (noise_type+'_{:03d}.npy'.format(channel)) )
    
    freqs, noise_spectral_density = periodogram(ch, fs=sampling_frequency, window="han")

    axis.plot(freqs, noise_spectral_density, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    axis.grid(color="gray", alpha=0.5)
    axis.set_xlabel("Frequency [Hz]")
    axis.set_xlim([freqs[0], freqs[-1]])
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylabel("Noise spectral density [1/Hz]")
    
    return freqs, noise_spectral_density
    
    
    
'''

'''
def plot_S21_MS2034B(ms2034b_object, axis, label=None, linestyle='solid', linewidth=1, color='black'):
    if label is None:
        label = ms2034b_object.label

    if ms2034b_object.mode == 'real_imag':
        S21 = np.sqrt(ms2034b_object.ReS21**2.0 + ms2034b_object.ImS21**2.0)
    if ms2034b_object.mode == 'log_mag_phase':
        S21 = ms2034b_object.S21DB

    axis.plot(ms2034b_object.freqs, S21, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    axis.grid(color="gray", alpha=0.5)
    axis.set_xlabel("Frequency [MHz]")
    axis.set_ylabel("$|S_{21}|$ [dBm]")



def electrical_phase_noise_equivalent_power(responsivity, noise_spectral_density, freqs, axis=None, label='ch', color='black', linestyle='solid', linewidth=1):

    NEP = np.sqrt(noise_spectral_density)/responsivity  

    axis.plot(freqs, NEP, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    axis.grid(color="gray", alpha=0.5)
    axis.set_xlabel("Frequency [Hz]")
    axis.set_xlim([freqs[0], freqs[-1]])
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylabel("NEP$_{el}^{ph}$ [W/$\sqrt{Hz}$]")
    
    return freqs, NEP


'''
            lsTarget function
'''
def lsTarget():
    '''
    This function lists all the target sweep directories

    Returns
    -------
    target_list : list of strings
        A list of strings with the name of target sweep directories.
    target_S21_list : list of strings
        A list of strings with the name of converted target sweep directories.

    '''
    import os
    print('List of target sweeps:')
    target_list = sorted([t for t in os.listdir(paths.target) if t[0] != '.'])
    for t in target_list:
        print(t)
    print('')
    print('List of target sweeps converted into different channels:')
    target_S21_list = sorted([t for t in os.listdir(paths.target_S21) if t[0] != '.'])
    for t in target_S21_list:
        print(t)
    print('')
    return target_list, target_S21_list
    

'''
            lsVNA function
'''
def lsVNA():
    '''
    This function lists all the VNA sweep directories
    
    Returns
    -------
    vna_list : list of strings
        A list of strings with the name of VNA sweep directories.
    vna_S21_list : list of strings
        A list of strings with the name of converted VNA sweep directories.

    '''
    import os
    print('List of VNA sweeps:')
    vna_list = sorted([vna for vna in os.listdir(paths.vna) if vna[0] != '.'])
    for vna in vna_list:
        print(vna)
    print('')
    print('List of VNA sweeps converted into different channels:')
    vna_S21_list = sorted([vna for vna in os.listdir(paths.vna_S21) if vna[0] != '.'])
    for vna in vna_S21_list:
        print(vna)
    print('') 
    return vna_list, vna_S21_list
    