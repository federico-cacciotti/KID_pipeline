from . import datapaths

class Dirfile():
    def __init__(self, filename, fs, first_frame, num_frames=10000, label=None):
        self.filename = filename
        self.fs = fs
        self.first_frame = first_frame
        self.num_frames = num_frames
        
        if label == None:
            self.label = filename
        else:
            self.label = label
  
    
    '''
    noise_type = 'p', 'a', 'I' or 'Q'
    '''
    def get_stream(self, channel, data_type, xdata_type='time'):
        import pygetdata
        
        # open time streams
        xdata = pygetdata.dirfile((datapaths.dirfile/self.filename).as_posix()).getdata("time", first_frame=self.first_frame, num_frames=self.num_frames)
        xdata -= xdata[0]
        ydata = pygetdata.dirfile((datapaths.dirfile/self.filename).as_posix()).getdata("ch{:s}_{:03d}".format(data_type, channel), first_frame=self.first_frame, num_frames=self.num_frames)
        
        if xdata_type == 'time':
            return xdata, ydata
        
        if xdata_type == 'index':
            xdata = range(self.first_frame, self.first_frame+ydata.size)
            return xdata, ydata
    
    def plot_stream(self, xdata, ydata, axis=None, label=None, color=None, linestyle='solid', linewidth=1, alpha=1.0):
        
        if axis == None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(7,7))
            axis = fig.gca()
        axis.grid(color="gray", alpha=0.5)
        
        axis.plot(xdata, ydata, color=color, label=label, linestyle=linestyle, linewidth=linewidth, alpha=alpha)