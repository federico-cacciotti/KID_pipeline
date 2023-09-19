from . import datapaths
import pygetdata
from pathlib import Path

class Dirfile():
    def __init__(self, filename, fs, label=None):
        self.filename = Path(filename)
        self.fs = fs
        
        if label == None:
            self.label = filename
        else:
            self.label = label
  
        self.dirfile = pygetdata.dirfile((datapaths.dirfile/self.filename).as_posix())
        self.nframes = self.dirfile.nframes

        print("Selected dirfile: {:s}".format(self.filename.as_posix()))
        print("Found {:d} frames".format(self.nframes))


    def get_time(self, first_frame=0, num_frames=None, remove_offset=False):
        if num_frames == None:
            num_frames = self.nframes

        time = self.dirfile.getdata("time", first_frame=first_frame, num_frames=num_frames)
        
        if remove_offset:
            time -= time[0]

        return time


    '''
    stream_type = 'p', 'a', 'I' or 'Q'
    '''

    def get_stream(self, channel, stream_type='p', first_frame=0, num_frames=None):
        if num_frames == None:
            num_frames = self.nframes
        
        stream = self.dirfile.getdata("ch{:s}_{:03d}".format(stream_type, channel), first_frame=first_frame, num_frames=num_frames)
        return stream
    

    def plot_stream(self, xdata, ydata, axis=None, label=None, color=None, linestyle='solid', linewidth=1, alpha=1.0):
        
        if axis == None:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(7,7))
            axis = fig.gca()
        axis.grid(color="gray", alpha=0.5)
        
        axis.plot(xdata, ydata, color=color, label=label, linestyle=linestyle, linewidth=linewidth, alpha=alpha)
