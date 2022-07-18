import pipeline as pl

'''
v140 = pl.VNA.VNA(filename='20220531_01', temperature=140)
peaks, peaks_info = v140.findPeaks(mag_filt='airp_lss', 
                    peak_width=(1.0, 100.0), 
                    peak_height=1.0, 
                    peak_prominence=(2.0, 80.0),
                    phase_filt='lowpass_cosine', phase_detrend=True)
'''


v270 = pl.VNA.VNA(filename='20220531_04', temperature=270)
peaks, peaks_info = v270.findPeaks(mag_filt='airp_lss', 
                    peak_width=(1.0, 150.0), 
                    peak_height=1.0, 
                    peak_prominence=(1.0, 30.0),
                    phase_filt='lowpass_cosine', phase_detrend=True)
