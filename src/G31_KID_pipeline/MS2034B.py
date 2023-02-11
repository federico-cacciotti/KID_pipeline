import numpy as np
from . import datapaths
from . import functions as fc

'''
            readMS2034B function
'''   
class MS2034B():
    def __init__(self, filename, label=None, mode='real_imag'):
        self.filename = filename
        self.mode = mode
        
        if label != None:
            self.label = label
        else:
            self.label = filename
        
        if self.mode == 'real_imag':
            [self.freqs, self.ReS11, self.ImS11, self.ReS21, self.ImS21, 
             self.ReS12, self.ImS12, self.ReS22, self.ImS22] = np.loadtxt(fname=datapaths.anritsuMS2034B / filename, dtype=float, skiprows=23, unpack=True)
        
        if self.mode == 'log_mag_phase':
            [self.freqs, self.S11DB, self.S11A, self.S21DB, self.S21A, 
             self.S12DB, self.S12A, self.S22DB, self.S22A] = np.loadtxt(fname=datapaths.anritsuMS2034B / filename, dtype=float, skiprows=23, unpack=True)
        
        self.freqs *= 1e3 # GHz to MHz
    
    
    def S21mag_dB(self):
        if self.mode == 'log_mag_phase':
            return self.S12DB
        if self.mode == 'real_imag':
            return 20.0*np.log10(np.sqrt(self.ReS21**2.0+self.ImS21**2.0))
        
    def S21mag(self):
        if self.mode == 'log_mag_phase':
            return 10.0**(self.S12DB/20.0)
        if self.mode == 'real_imag':
            return np.log10(np.sqrt(self.ReS21**2.0+self.ImS21**2.0))
                            
    def S21_I(self):
        if self.mode == 'log_mag_phase':
            return self.S21mag()*np.cos(self.S21A)
        if self.mode == 'real_imag':
            return self.ReS21
    
    def S21_Q(self):
        if self.mode == 'log_mag_phase':
            return self.S21mag()*np.sin(self.S21A)
        if self.mode == 'real_imag':
            return self.ImS21
        
    def S21phase(self):
        if self.mode == 'log_mag_phase':
            return self.S21A
        if self.mode == 'real_imag':
            return np.arctan2(self.ImS21, self.ReS21)
    
    def fitS21(self, DATAPOINTS=3500):
        I = self.S21_I()
        Q = self.S21_Q()
        
        amp = self.S21mag_dB()
        res_freq = self.freqs[np.argmin(amp)]
        out_path = datapaths.anritsuMS2034B
        
        np.save(out_path/"I.npy", arr=I)
        np.save(out_path/"Q.npy", arr=Q)
        np.save(out_path/"freqs.npy", arr=self.freqs)
        np.save(out_path/"mag.npy", arr=amp)
        
        try:
            params, chi2 = fc.complexS21Fit(I=I, Q=Q, freqs=self.freqs, res_freq=res_freq, 
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
        
        target_path = datapaths.anritsuMS2034B
        fc.complexS21Plot(target_path)
    
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