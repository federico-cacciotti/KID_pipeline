# Example of reading a VNA sweep file, finding resonances, extracting a single target resonance,
# and complex fitting the chosen resonance

import G31_KID_pipeline as pl
import os

# specify the working directory where data files are stored
pl.paths(working_directory=os.getcwd())


# open a VNA sweep file
vna = pl.VNA(filename='20220531_04', temperature=270)

# Figure 1: a plot of the raw VNA sweep
vna.plotVNA(xlim=[0, 700], mag_filt=None, phase_filt=None)

# Figure 2: a plot of the VNA sweep without the baselines
vna.plotVNA(xlim=[0, 700], mag_filt='airp_lss', phase_filt='lowpass_cosine')

# Figure 3: find resonances
peaks, peaks_info = vna.findPeaks(mag_filt='airp_lss', 
                    peak_width=(1.0, 100.0), 
                    peak_height=1.0, 
                    peak_prominence=(2.0, 80.0),
                    phase_filt='lowpass_cosine')

# Figure 4: extracted target data
target_extr = vna.extractTarget(peaks, peaks_info, extr_width=3.0)

# complex fit on the resonance 1
vna.fitS21(target_extr, channel=1)

# Figure 5: plot of the complex fit
vna.plotS21(channel=1)