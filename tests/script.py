import G31_KID_pipeline as pl
import os

# specify the working directory where data files are stored
pl.paths(working_directory=os.getcwd())

'''
v140 = pl.VNA(filename='20220531_01', temperature=140)
peaks, peaks_info = v140.findPeaks(mag_filt='airp_lss', 
                    peak_width=(1.0, 100.0), 
                    peak_height=1.0, 
                    peak_prominence=(2.0, 80.0),
                    phase_filt='lowpass_cosine', phase_detrend=True)
'''


'''
v270 = pl.VNA(filename='20220531_04', temperature=270)
v270.plotVNA(mag_filt=None, phase_filt=None)

peaks, peaks_info = v270.findPeaks(peak_width=(1.0, 150.0), peak_height=1.0, peak_prominence=(1.0, 30.0))

v270_target = v270.extractTarget(peaks, peaks_info, extr_width=3.0)
v270.fitS21(v270_target, channel=1)
v270.plotS21(1)
'''

t = pl.Target(filename='20220530_145142')

t.fitS21(channel=1)

t.plotS21(1)
