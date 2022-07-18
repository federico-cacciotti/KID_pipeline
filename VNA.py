import numpy as np
import paths
import pipeline as pl

class VNA():
    def __init__(self, filename, temperature=None, build_dataset=False):
        self.filename = filename
        self.temperature = temperature
        
        self.bb_freqs = np.load(paths.vna / self.filename / 'bb_freqs.npy')
        
        if not (paths.vna_S21 / self.filename).exists() or build_dataset:
            pl.buildS21Dataset(self)
        
        
    def findPeaks(self, xlim=[0, 700], mag_filt='airp_lss', phase_filt='lowpass_cosine', phase_detrend=True, peak_width=(1.0, 150.0), peak_height=1.0, peak_prominence=(1.0, 30.0)):
        
        print('Peak finding...')
        print('Magnitude baseline correction algorithm: '+mag_filt)
        if phase_detrend:
            print('Phase detrend: True')
        else:
            print('Phase detrend: False')
        print('Phase baseline correction algorithm: '+phase_filt)
        
        from matplotlib import pyplot as plt
        # how many channels?
        n_chan = self.bb_freqs.size
        
        fig = plt.figure()
        fig.set_size_inches(7.5, 8)
        plt.subplots_adjust(bottom=0.15, right=0.98, top=0.95, left=0.1, wspace=0.35)
        ax0 = plt.subplot(211)
        ax1 = plt.subplot(212, sharex=ax0)
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
        if phase_detrend:
            from scipy.signal import detrend
            detrend(phase, axis=-1, type='linear', overwrite_data=True)
            phase -= phase[0]
        
        # mag filter
        if mag_filt == 'lowpass_cosine':
            sweep_step = 1.25 # kHz
            smoothing_scale = 2500.0 # kHz
            filtered = pl.lowpass_cosine(mag, sweep_step, 1.0/smoothing_scale, 0.1 * (1.0/smoothing_scale))
            mag -= filtered    
            
        if mag_filt == 'a_lss':
            baseline = pl.asymmetric_least_squares_smoothing(data=mag, lam=1e6, p=8e-4, N_iter=5)
            mag -= baseline
            
        if mag_filt == 'airp_lss':
            baseline = pl.adaptive_iteratively_reweighted_penalized_least_squares_smoothing(data=mag, lam=1e6, N_iter=5)
            mag -= baseline
        
        
        # phase filter
        if phase_filt == 'lowpass_cosine':
            sweep_step = 1.25 # kHz
            smoothing_scale = 2500.0 # kHz
            filtered = pl.lowpass_cosine(phase, sweep_step, 1./smoothing_scale, 0.1 * (1.0/smoothing_scale))
            phase -= filtered
            
            
        # peak finding
        from scipy.signal import find_peaks
        peaks, peaks_info = find_peaks(-mag, 
                                       width=peak_width, 
                                       height=peak_height, 
                                       prominence=peak_prominence)
        print("Found ", len(peaks), " resonances")
        
        ax0.plot(freqs[peaks], mag[peaks], "x", label="Resonance")
        ax0.legend(loc='best')
            

        ax0.plot(freqs, mag, color='black', linewidth=1)
        ax1.plot(freqs, phase, color='black', linewidth=1)
        
        ax0.grid(linestyle='-', alpha=0.5)
        ax0.set_xlim([freqs[0], freqs[-1]])
        ax0.set_ylabel('Mag [dB]')
        ax0.set_xlabel('Frequency [MHz]')
        
        ax1.grid(linestyle='-', alpha=0.5)
        ax1.set_xlim([freqs[0], freqs[-1]])
        ax1.set_ylabel('Phase [rad]')
        ax1.set_xlabel('Frequency [MHz]')
        
        plt.show()
        
        return peaks, peaks_info
        
        