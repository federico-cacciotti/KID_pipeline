from . import datapaths

class Dirfile():
    def __init__(self, filename, fs, first_frame, num_frames):
        self.filename = filename
        self.fs = fs
        self.first_frame = first_frame
        self.num_frames = num_frames
  
    
    '''
    noise_type = 'p', 'a', 'I' or 'Q'
    '''
    def plot_stream(self, channel, data_type, first_frame, num_frames, axis, xdata='frame', label=None, color=None, linestyle='solid', linewidth=1):
        import pygetdata

        if label is None:
            label = "ch{:s}_{:03d}".format(data_type, channel)

        # open time streams
        ctime = pygetdata.dirfile((datapaths.dirfile/self.filename).as_posix()).getdata("time", first_frame=first_frame, num_frames=num_frames)
        ctime -= ctime[0]
        ch = pygetdata.dirfile((datapaths.dirfile/self.filename).as_posix()).getdata("ch{:s}_{:03d}".format(data_type, channel), first_frame=first_frame, num_frames=num_frames)
        
        if xdata == 'time':
            axis.plot(ctime, ch, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
            axis.set_xlabel("Time [s]")
        if xdata == 'frame':
            axis.plot(range(first_frame, first_frame+ch.size), ch, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
            axis.set_xlabel("Frames")
        
        axis.grid(color="gray", alpha=0.5)
        
        axis.set_ylabel(data_type)



    '''
    noise_type = 'p', 'a', 'I' or 'Q'
    '''
    def plot_psd(self, channel, data_type, axis, label=None, color=None, linestyle='solid', linewidth=1):
        import pygetdata
        from scipy.signal import periodogram, welch
        import numpy as np
        
        if label is None:
            label = "ch{:s}_{:03d}".format(data_type, channel)

        # open time streams
        ctime = pygetdata.dirfile((datapaths.dirfile/self.filename).as_posix()).getdata("time", first_frame=self.first_frame, num_frames=self.num_frames)
        ctime -= ctime[0]
        ch = pygetdata.dirfile((datapaths.dirfile/self.filename).as_posix()).getdata("ch{:s}_{:03d}".format(data_type, channel), first_frame=self.first_frame, num_frames=self.num_frames)
        
        freqs, spectrum = periodogram(ch, fs=self.fs, scaling="density", window="hann")
        #freqs, spectrum = welch(ch, fs=fs, window='hann', scaling='density')
        
        psd = np.sqrt(spectrum)
        
        axis.plot(freqs, psd, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.grid(color="gray", alpha=0.5)
        axis.set_xlabel("Frequency [s]")
        axis.set_ylabel("PSD [1/$\sqrt{Hz}$]")