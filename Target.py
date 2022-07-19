import numpy as np
import paths
import sys
import pipeline as pl

class Target():
    def __init__(self, filename, temperature=None, build_dataset=False):
        print('Reding target sweep {:s}...'.format(filename))
        self.filename = filename
        self.temperature = temperature
        
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
            pl.buildS21Dataset(self)
        
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
        
        self.out_of_res = self.filterOutOfResTones(std_mult=3.0)
        
        self.readS21Data()
        
        self.findDouble()
    

    # filter out of resonance tones refusing them, sigma_mult allows to refuse out of resonance tones
    def filterOutOfResTones(self, std_mult=3.0):
        print("\nChecking if there are out of resonance tones...")
        
        depths = [self.entry[i]['depth'] for i in range(self.entries)]
        average_depth = np.average(depths)
        print("average depth: {:f}".format(average_depth))
        std_depth = np.std(depths)
        print("std depth: {:f}".format(std_depth))
            
        threshold = average_depth-std_mult*std_depth
        
        out_of_res_number = 0
        for e in self.entry:
            if e['depth']<=threshold:
                e['is_out_of_res'] = True
                out_of_res_number += 1
        
        print("Found {:d} out of resonance tones: ".format(out_of_res_number))
        if out_of_res_number>0:
            for i,e in enumerate(self.entry):
                if e['is_out_of_res']:
                    print("\tChannel {:d} - {:f} MHz".format(i, e['target_freq']))
                    
        return out_of_res_number

    def fitS21(self, channel):
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
                        params, chi2 = pl.complexS21Fit(I=I, Q=Q, freqs=freqs, res_freq=e['target_freq'], 
                                               output_path=out_path, DATAPOINTS=70)
                    
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
                params = pl.complexS21Fit(I=I, Q=Q, freqs=freqs, res_freq=self.entry[channel]['target_freq'], 
                              output_path=out_path, DATAPOINTS=70, verbose=True)
                
                self.entry[channel]['Re[a]'] = params['Re[a]']
                self.entry[channel]['Im[a]'] = params['Im[a]']
                self.entry[channel]['Q_tot'] = params['Q_tot']
                self.entry[channel]['Q_c'] = params['Q_c']
                self.entry[channel]['Q_i'] = params['Q_i']
                self.entry[channel]['nu_r'] = params['nu_r']
                self.entry[channel]['phi_0'] = params['phi_0']
                self.entry[channel]['reduced_chi2'] = params['reduced_chi2']
            
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
        
        accumulation_length = 2**21
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
