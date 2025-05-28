import pygetdata

dirfile_name = "dirfile_20250410_111700"

# open a dirfile
dirfile = pygetdata.dirfile(dirfile_name)


# load stream
first_frame = 100 # skip first N frames
num_frames = 10000 # number of frames to load

# data_string must be in the form "chx_nnn"
# where
# x = p -> normalized loop phase [rad]
# x = a -> normalized loop amplitude [ADC units]
# x = I -> raw in-phase [ADC units]
# x = Q -> raw quadrature [ADC units]
# and (example)
# nnn = 004 -> is the readout channel 4 (always with three digits)

channel = 0
stream_type = 'p'
data_string = "ch{:s}_{:03d}".format(stream_type, channel)

stream = dirfile.getdata(data_string, first_frame=first_frame, num_frames=num_frames)
time = dirfile.getdata("time", first_frame=first_frame, num_frames=num_frames) # unix timestamp


# plots, fft, ...

# for the fft there exist many routines 
# see scipy.fft.fft, scipy.signal.periodogram, scipy.signal.welch, numpy.fft.fft, ...


# always close the dirfile to avoid memory issues
dirfile.close()