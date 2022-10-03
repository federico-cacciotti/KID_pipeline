# Example of reading a Target sweep file

import G31_KID_pipeline as pl
import os

# specify the working directory where data files are stored
pl.paths(working_directory=os.getcwd())

# open a Target sweep file
target = pl.Target(filename='20220530_145142', temperature=140)

# Figure 1: interactive plot of the Target sweep
target.plotTarget()

# Figure 2: plot of channel 1
target.plotChannel(channel=1)

# complex fit on the resonance 1
target.fitS21(channel=1)

# Figure 3: plot of the complex fit
target.plotS21(1)

