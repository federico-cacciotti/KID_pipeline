# G31_KID_pipeline package
A pipeline for Kinetic Inductance Detectors data analysis and characterization that I developed during my master thesis and PhD.

# How to install
The use of this package is recommended under Linux or Mac OS. This software is developed under Mac OS.

### Using the pip package
Make sure you are using a `>=3.7` Python version and the latest `pip` version available by executing, for Mac Os / Linux users, the following command
```shell
python3 -m pip install --upgrade pip
```
or, for Windows users,
```shell
py -m pip install --upgrade pip
```
Then proceed with the package installation by typing
```shell
pip install G31_KID_pipeline
```
and you are done.

# Stable versions and changelog
- 1.0.3 (October 4rd, 2022) - added `lsTarget()` and `lsVNA()` functions and function descriptions
- 1.0.2 (October 3rd, 2022) - first release

# Overview
A short overview of the package is given here. To import the package just include it in your Python code as follows:

```Python
import G31_KID_pipeline as pl
```
Once the package is imported in your code, you must specify the path to the working directory where all the data files are stored. If the data is in the same folder of your Python script than the easiest way is to use the `os.getcwd()` function:
```python
import os
# specify the working directory where data files are stored
pl.paths(working_directory=os.getcwd())
```
In this way it automatically checks if all the data directories exist. If it doesn't find one or more data directories, it creates them in the same path of your Python script. At the end of this operation you should have the following directories in your script path:

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

You can open a `VNA` file by defining a `VNA` object as

```Python
vna = pl.VNA(filename='20220531_04', temperature=270, build_dataset=True)
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
With this function it is possible to extract a target sweep from a VNA sweep in order to fit the resonance circles and measure the electrical parameters of a detector. This function takes in input the `peaks` and `peaks_info` parameters returned by the `vna.findPeaks()` function and returns a dictionary with I, Q, magnitude, phase and frequency for all the resonances found. The optional `extr_width` parameter (default is 2.5) is the frequency range in unit of FWHM of each extracted resonance. The following plot shows the extracted target sweep from a VNA sweep.
<image src="images/vna_plot_targetext.png" width="100%">


#### `vna.fitS21(extracted_target, channel)` <br>
This function performs a complex fit of the S21 scattering parameter on an extracted target resonance from a target sweep. The fitting routine is described under the `pl.complexS21Fit()` function section. The `channel` parameter can be an integer value between 0 and the number of extracted resonances -1 or `'all'` in order to perform the complex on all the resonances at once.

This is an example of this function, the typical output will be:

```shell
Fit S21 complex function for vna_S21/20220531_04/extracted_target/001...

	centerd circle parameters
	Xc =  -4106.814783255634
	Yc =  13468.54613369031
	radius =  4041.616970887738

	rotation angle
	angle =  1.4642370255123929

	phase fit
	[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 21
    # data points      = 40
    # variables        = 2
    chi-square         = 3300.41786
    reduced chi-square = 86.8531016
    Akaike info crit   = 180.516996
    Bayesian info crit = 183.894755
[[Variables]]
    Q_tot:  14783.0227 +/- 415.342964 (2.81%) (init = 5000)
    nu_r:   201.026147 +/- 5.3325e-05 (0.00%) (init = 201.0256)
[[Correlations]] (unreported correlations are < 0.100)
    C(Q_tot, nu_r) = -0.105

	complex fit results
	[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 221
    # data points      = 40
    # variables        = 6
    chi-square         = 64.3328892
    reduced chi-square = 1.89214380
    Akaike info crit   = 31.0076617
    Bayesian info crit = 41.1409384
[[Variables]]
    Rea:   -42.5742370 +/- 69.2660770 (162.69%) (init = -4509.955)
    Ima:    18188.3458 +/- 143.319915 (0.79%) (init = 17237.46)
    Q_tot:  13858.0282 +/- 491.054533 (3.54%) (init = 14783.02)
    Q_c:    30879.1224 +/- 910.735366 (2.95%) (init = 32585.86)
    nu_r:   201.025926 +/- 1.6956e-04 (0.00%) (init = 201.0261)
    phi_0:  6.15323124 +/- 0.02365153 (0.38%) (init = 6.133844)
    tau:    0.04 (fixed)
[[Correlations]] (unreported correlations are < 0.100)
    C(nu_r, phi_0) = -0.929
    C(Ima, Q_c)    = -0.901
    C(Q_tot, Q_c)  = 0.892
    C(Ima, Q_tot)  = -0.782
    C(Rea, phi_0)  = -0.510
    C(Rea, nu_r)   = 0.454
    C(Rea, Ima)    = 0.308
    C(Rea, Q_c)    = -0.277
    C(Q_tot, nu_r) = 0.163
    C(Rea, Q_tot)  = -0.147
    C(Ima, nu_r)   = -0.130
    C(Q_c, nu_r)   = 0.108

Q_i: 25141+/-1725
```

This function will save the fitted parameters in a `.npy` file inside the VNA sweep folder inside the `vna_S21` directory.


#### `vna.plotS21(channel)` <br>
This function is useful to visualize the result of the complex fit on the extracted resonances of a VNA sweep. The `channel` parameter must be an integer number between 0 and the number of extracted resonances (-1). This function calls the `pl.complexS21Plot()` routine.

This is an example of a complex fit on an extracted resonance from a VNA sweep.
<image src="images/complex_vna_fit.png" width="100%">


# `target` object

You can open a `target` file by defining a `target` object as

```Python
target = pl.Target(filename='20220530_145142', temperature=150, build_dataset=True)
```
where the optional `temperature` parameter (default is `None`) defines the temperature of the target sweep expressed in milli kelvin and the optional `build_dataset` parameter (default is `False`) forces the conversion of the raw target sweep files into `.npy` files by creating a subfolder under the `target_S21` directory with the same name of the `filename` parameter.

When a `target` object is defined a series of operations are performed in background. A dictionary `entry` will be created as an attribute of the `target` object. The dictionary has the following items:
- `'target_freq'` a list with the target frequencies of the sweep;
- `'channel'` a list indexes that label the resonances;
- `'depth'` a list of the resonance depths in dB;
- `'is_out_of_res'` a list of boolean values for each tone (`True` is the tone is out of resonance and `False` if it is a resonance, `False` by default);
- `'number_of_peaks'` a list of integer number representing the number of peaks for each tone (useful to determine if a single tone sees more than one pixel);

the next few items are used to store the output values of the complex fit routine for each resonance:
- `'Re[a]'`;
- `'Im[a]'`;
- `'Q_tot'`;
- `'Q_i'`;
- `'Q_c'`;
- `'nu_r'`;
- `'phi_0'`;
- `'reduced_chi2'`.

The next operation concerns the search for out of resonance tones with the function `target.filterOutOfResTones()`, the search in memory for the already computed S21 complex fit parameters with the function `target.readS21Data()` (if present) and the search for double resonances in the target sweep with the `target.findDouble()` function. All these functions are explained in the next section.


## Operations with the `target` object

#### `target.filterOutOfResTones(std_mult)` <br>
This function is called when a `target` object is defined. It selects and labels all the out of resonance tones (if present). The selection is done by comparing each tone depth with the target sweep average depth. If the depth of a tone is lower than `std_mult` times the standard deviation over the target sweep depths then it is labeled as an out of resonance tone and its corresponding value in the list `'is_out_of_res'` will be turned to `True`.


#### `target.readS21Data()` <br>
This function is called when a `target` object is defined. It will try to search for complex fit parameters of each resonance already present in the memory.


#### `target.findDouble()` <br>
This function is called when a `target` object is defined. It will try to search for double (or more) resonances in each tone of the target sweep. If it find a double resonance the corresponding value in the list `'number_of_peaks'` will be turned to the number of resonances found.

#### `target.plotTarget(flat_at_0db)` <br>
This is a function that plots the amplitude (in dB) of a target sweep. The optional boolean parameter `flat_at_0db` (`False` by default) if `True`, allows to plot all the tones with their highest point at the 0dB level. The output plot will be an interactive plot that will show the channel number of a tone when the mouse is over it, like in the example image below.
<image src="images/target_plot.png" width="100%">


#### `target.plotChannel(channel)` <br>
This function is similar to the previous one with the exception that it plots the IQ data, the amplitude and the phase of a tone indexed by the value of the `channel` parameter. The following image shows an example.
<image src="images/target_plot_single_channel.png" width="100%">


#### `target.plotS21(channel)` <br>
This function is useful to visualize the result of the complex fit on a channel of a target sweep. The `channel` parameter must be an integer number between 0 and the number of channels (-1). This function calls the `pl.complexS21Plot()` routine.

This is an example of a complex fit on a target sweep channel.
```shell
Fit S21 complex function for target_S21/20220530_145142/001...

	centerd circle parameters
	Xc =  -350852.7043057035
	Yc =  436318.853189466
	radius =  337206.6320412996

	rotation angle
	angle =  0.9912039193143776

	phase fit
	[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 19
    # data points      = 31
    # variables        = 2
    chi-square         = 19.1696365
    reduced chi-square = 0.66102195
    Akaike info crit   = -10.9004479
    Bayesian info crit = -8.03247352
[[Variables]]
    Q_tot:  25316.3945 +/- 78.2077610 (0.31%) (init = 5000)
    nu_r:   201.232202 +/- 3.0830e-06 (0.00%) (init = 201.2318)
[[Correlations]] (unreported correlations are < 0.100)
    C(Q_tot, nu_r) = -0.170

	complex fit results
	[[Fit Statistics]]
    # fitting method   = leastsq
    # function evals   = 484
    # data points      = 31
    # variables        = 6
    chi-square         = 8.37844104
    reduced chi-square = 0.33513764
    Akaike info crit   = -28.5580856
    Bayesian info crit = -19.9541623
[[Variables]]
    Rea:   -285363.955 +/- 1722.87122 (0.60%) (init = -520071.1)
    Ima:    847947.008 +/- 2399.48949 (0.28%) (init = 694830.9)
    Q_tot:  25690.3985 +/- 205.204873 (0.80%) (init = 25316.39)
    Q_c:    33902.2276 +/- 160.239752 (0.47%) (init = 32579.87)
    nu_r:   201.232238 +/- 3.6001e-05 (0.00%) (init = 201.2322)
    phi_0:  6.21302085 +/- 0.00710776 (0.11%) (init = 6.220246)
    tau:    0.04 (fixed)
[[Correlations]] (unreported correlations are < 0.100)
    C(nu_r, phi_0)  = -0.892
    C(Q_tot, Q_c)   = 0.754
    C(Ima, Q_c)     = -0.524
    C(Ima, Q_tot)   = -0.502
    C(Ima, nu_r)    = 0.493
    C(Q_tot, nu_r)  = -0.487
    C(Ima, phi_0)   = -0.408
    C(Rea, phi_0)   = -0.396
    C(Rea, nu_r)    = 0.342
    C(Rea, Ima)     = -0.330
    C(Rea, Q_c)     = 0.215
    C(Q_tot, phi_0) = 0.145
    C(Q_c, nu_r)    = -0.136

Q_i: 106062+/-3833
```
This function will save the fitted parameters in a `.npy` file inside the target sweep folder inside the `target_S21` directory.


#### `target.plotS21(channel)` <br>
This function is useful to visualize the result of the complex fit on the extracted resonances of a VNA sweep. The `channel` parameter must be an integer number between 0 and the number of extracted resonances (-1). This function calls the `pl.complexS21Plot()` routine.

This is an example of a complex fit on an extracted resonance from a VNA sweep.
<image src="images/complex_target_fit.png" width="100%">


# Other operations with the `pipeline` class
#### `pl.buildS21Dataset(sweep, ROACH='MISTRAL')` <br>
Sweep files are generated by the ROACH readout electronics as a list of `.npy` files each one storing the in-quadrature or the in-phase (Q or I respectively) data relative to a single base band frequency. This data is converted into a more readable way in which a single `.npy` file stores a single frequency sweep data as I, Q, magnitude in dB, phase in radians and frequency. The `sweep` parameter can be both a `Target` or a `VNA` object and the optional `ROACH` parameter (default is `'MISTRAL'`) determines how the sweep frequencies are computed starting from the Local Oscillator frequency and the base bands frequencies:
- `ROACH='MISTRAL'` --> `freqs = LO + bb_feqs`;
- `ROACH='OLIMPO'` --> `freqs = 0.5*LO + bb_feqs`.

When defining a new `VNA` or `Target` class, the function ```pl.buildS21Dataset(sweep, ROACH='MISTRAL')``` will be called automatically only if the S21 dataset is not already present. The conversion of the raw sweep dataset can be forced by setting the `build_dataset` parameter to `True`.
