# KID_pipeline package
A pipeline for Kinetic Inductance Detectors data analysis and characterization that I developed during my master thesis and PhD.

# Overview
A short overview of the package is given here. To import the package just include it in your Python code as follows:

```Python
import pipeline as pl
```
Once the package is imported in your code, it automatically checks if all the data directories exist. If it doesn't find one or more data directories, it creates them in the same path of your Python script. At the end of this operation you should have the following directories in your script path:

```sh
anritsu_MS2034B
cosmic_rays
data_logger
noise
picolog
target
target_S21
vna
vna_S21
your_python_script.py
```

- `anritsu_MS2034B` should contain the output files from the Vector Network Analyzer Anritsu MS2034B;
- `cosmic_rays` should contain the output files of cosmic ray events taken using a Tektronix DPO3054 digital oscilloscope;
- `data_logger` should contain the dirfiles coming from the ROACH readout;
- `noise` should contain the extracted noise time streams in `.npy` file format;
- `picolog` should contain the PicoLog output files;
- `target` should contain the target sweep files coming from the ROACH readout;
- `target_S21` should contain the target sweep files converted into `.npy` files for each channel;
- `vna` should contain the VNA sweep files coming from the ROACH readout;
- `vna_S21` should contain the VNA sweep files converted into `.npy` files for each channel.



# `VNA` object

You can open a `VNA` file by defining a `VAN` object as

```Python
vna = pl.VNA.VNA(filename='20220531_04', temperature=270, build_dataset=True)
```
where the optional `temperature` parameter (default is `None`) defines the temperature of the VNA sweep expressed in milli kelvin and the optional `build_dataset` parameter (default is `False`) forces the conversion the raw VNA sweep files into `.npy` files by creating a subfolder under the `VNA_S21` directory with the same name of the `filename` parameter.

## Operations with the `VNA` object
#### `vna.plotVNA(xlim, mag_filt, phase_filt)` <br>
This is a function that plots the amplitude (in dB) and phase (in radians) of a VNA sweep in the frequency range defined by the `xlim` parameter. The optional parameters `mag_filt` and `phase_filt` allow to chose an algorithm to remove the baseline from the VNA sweep magnitude and phase data respectively. They can be equal to one of the following values:
- `None`: no baseline will be removed;
- `'lowpass_cosine'`: the function `pl.lowpass_cosine()` will be called;
- `'a_lss'`: the function `pl.asymmetric_least_squares_smoothing()` will be called (this is not available for `phase_filt`);
- `'airp_lss'`: the function `pl.adaptive_iteratively_reweighted_penalized_least_squares_smoothing()` will be called (this is not available for `phase_filt`).

This is an example showing the original (without filtering out the baselines) and the flattened magnitude and phase data using an `airp_lss` algorithm for the `mag_filt` and a `lowpass_cosine` for the `phase_filt`. <br>
<image src="images/vna_plot.png" width="100%">
<image src="images/vna_plot_filt.png" width="100%">

#### `vna.findPeaks(xlim, mag_filt, phase_filt, peak_width, peak_height, peak_prominence)` <br>
This function find the resonances dips in the VNA sweep and marks the found resonances with a cross. `xlim` is the frequency range of the plot, `mag_filt` and `phase_filt` are described in the previous function, `peak_width`, `peak_height` and `peak_prominence` are used to recognise local minima of the sweep data (see `scipy.signal.find_peaks()`). Currently only the magnitude data is used to find resonances. Following an example of application of this function where it found 353 resonances.
<image src="images/vna_plot_findpeaks.png" width="100%">


#### `vna.extractTarget(peaks, peaks_info, extr_width)` <br>
With this function it is possible to extract a target sweep from a VNA sweep in order to fit the resonance circles and measure the electrical parameters of a detector. This function takes in input the `peaks` and `peaks_info` parameters returned by the `vna.findPeaks()` function and returns a dictionary with I, Q, magnitude, phase and frequency for all the resonances found. The optional `extr_width` parameter (default is 2.5) is the frequency range in unit of FWHM of each extracted resonance.
<image src="images/vna_plot_targetext.png" width="100%">

# `target` object

You can open a `target` file by defining a `target` object as

```Python
target = pl.Target.Target(filename='20220530_145142', temperature=150, build_dataset=True)
```
where the optional `temperature` parameter (default is `None`) defines the temperature of the target sweep expressed in milli kelvin and the optional `build_dataset` parameter (default is `False`) forces the conversion the raw target sweep files into `.npy` files by creating a subfolder under the `target_S21` directory with the same name of the `filename` parameter.

## Operations with the `target` object



# Other operations with the `pipeline` class
#### `pl.buildS21Dataset(sweep, ROACH='MISTRAL')` <br>
Sweep files are generated by the ROACH readout electronics as a list of `.npy` files each one storing the in-quadrature or the in-phase (Q or I respectively) data relative to a single base band frequency. This data is converted into a more readable way in which a single `.npy` file stores a single frequency sweep data as I, Q, magnitude in dB, phase in radians and frequency. The `sweep` parameter can be both a `Target` or a `VNA` object and the optional `ROACH` parameter (default is `'MISTRAL'`) determines how the sweep frequencies are computed starting from the Local Oscillator frequency and the base bands frequencies:
- `ROACH='MISTRAL'` --> `freqs = LO + bb_feqs`;
- `ROACH='OLIMPO'` --> `freqs = 0.5*LO + bb_feqs`.

When defining a new `VNA` or `Target` class, the function ```pl.buildS21Dataset(sweep, ROACH='MISTRAL')``` will be called automatically only if the S21 dataset is not already present. The conversion of the raw sweep dataset can be forced by setting the `build_dataset` parameter to `True`.
