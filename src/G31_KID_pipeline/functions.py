import numpy as np
from matplotlib.lines import Line2D
import sys
from . import datapaths


#NEP_photon = 5.e-15 #W/Hz^1/2
#NEP_photon_noatm = 8.e-17 #W/Hz^1/2


'''
        Function for overplotting targets and MS2034B data
'''
def overplotTargetSweeps(targets=None, ms2034b_data_list=None, channel_index=False, add_out_of_res_plot=False, complex_fit_above=False, flat_at_0db=True, colormap='coolwarm', markers=True, only_idxs=None):
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import colormaps as cm
    cmap = cm.get_cmap(colormap)
    
    fig = plt.figure()
    fig.set_size_inches(7, 7)
    ax0 = plt.subplot(111)
    
    ax0.yaxis.set_ticks_position('both')
    ax0.xaxis.set_ticks_position('both')
    ax0.minorticks_on()
    ax0.yaxis.set_tick_params(direction='in', which='both')
    ax0.xaxis.set_tick_params(direction='in', which='both')
    ax0.grid(linestyle='-', alpha=0.5)
    ax0.set_ylabel('Mag [dB]')
    ax0.set_xlabel('Frequency [MHz]')
    
    # plot roach target sweeps
    handles = []
    if targets != None:
        for i,target in enumerate(targets):
            try:
                color = cmap(i/(len(targets)-1))
            except ZeroDivisionError:
                color = cmap(0)
            if only_idxs != None:
                target_entries = [target.entry[idx] for idx in only_idxs]
            else:
                target_entries = target.entry
            for e in target_entries:
                # read one sweep at a time
                if not e['is_out_of_res'] or add_out_of_res_plot:
                    channel = e['channel']
                    x_data_chan = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "freqs.npy")
                    y_data_chan = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "mag.npy")
                    
                    if flat_at_0db:
                        hist, bin_edges = np.histogram(y_data_chan, bins=70)
                        hist_arg_max = hist.argmax()
                        y_offset = 0.5*(bin_edges[hist_arg_max]+bin_edges[hist_arg_max+1])
                        y_data_chan -= y_offset
                
                    if markers:
                        ax0.plot(x_data_chan, y_data_chan, color=color, linestyle='-', linewidth=1, marker='o', markersize=3.5)
                    else:
                        ax0.plot(x_data_chan, y_data_chan, color=color, linestyle='solid')
                    
                    # plot the channel index above the correspondend sweep
                    if channel_index:
                        ax0.text(x_data_chan[0], y_data_chan[0], str(e['channel']), color=color)
                    
                    if complex_fit_above:
                        try:
                            nu_linear = np.linspace(x_data_chan[0], x_data_chan[-1], num=200) # linear sample
                            nu_peack = np.random.normal(e['nu_r'].n, 0.001, 1000) # peak sampling
                            nu = np.concatenate([nu_linear, nu_peack])
                            nu = np.sort(nu)
                            Z = S_21(nu, e['Re[a]'].n, e['Im[a]'].n, e['Q_tot'].n, e['Q_c'].n, e['nu_r'].n, e['phi_0'].n)
                            mag = 20*np.log10(np.abs(Z))
                            if flat_at_0db:
                                mag -= y_offset
                            ax0.plot(nu, mag, linestyle='solid', color=color, alpha=0.6, linewidth=4)
                        except:
                            pass
                
            handles.append(Line2D([0], [0], label=target.label, color=color))
            
    # plot ms2034b VNA sweeps
    if ms2034b_data_list != None:
        for i,sweep in enumerate(ms2034b_data_list):
            if sweep.mode == 'log_mag_phase':
                amp = sweep.S21DB
                ph = sweep.S21A
                ph = np.unwrap(ph)
            if sweep.mode == 'real_imag':
                amp = 20*np.log10(np.sqrt(sweep.ReS21**2 + sweep.ImS21**2))
                ph = np.arctan2(sweep.ImS21, sweep.ReS21)
                ph = np.unwrap(ph)
            
            if flat_at_0db:
                amp -= amp[0]
            
            color = cmap(i/(len(ms2034b_data_list)-1))
            ax0.plot(sweep.freqs, amp, color=color, linestyle='dotted', linewidth=1)
            handles.append(Line2D([0], [0], linestyle='dotted', label=sweep.label, color=color))
    
    ax0.legend(loc='best', handles=handles, fontsize=8)
    
    plt.show()



def overplotTargetCircles(targets=None, ms2034b_data_list=None, complex_fit_above=False, colormap='coolwarm', markers=True, only_idxs=None, force_raw_data=False):
    from matplotlib import pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib import colormaps as cm
    cmap = cm.get_cmap(colormap)
    
    fig = plt.figure(figsize=(7,7))
    ax0 = fig.gca()
    
    ax0.yaxis.set_ticks_position('both')
    ax0.xaxis.set_ticks_position('both')
    ax0.minorticks_on()
    ax0.yaxis.set_tick_params(direction='in', which='both')
    ax0.xaxis.set_tick_params(direction='in', which='both')
    ax0.grid(linestyle='-', alpha=0.5)
    ax0.set_ylabel('Q')
    ax0.set_xlabel('I')
    ax0.set_aspect('equal')
    
    # plot roach target sweeps
    handles = []
    if targets != None:
        for i,target in enumerate(targets):
            
            color = cmap(i/(len(targets)-1))
            
            if only_idxs != None:
                target_entries = [target.entry[idx] for idx in only_idxs]
            else:
                target_entries = target.entry
            for e in target_entries:
                # read one sweep at a time
                if e['is_out_of_res']==False:
                    channel = e['channel']
                    x_data_chan = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "freqs.npy")
                    
                    try:
                        if force_raw_data:
                            raise FileNotFoundError
                        I = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "I_prime.npy")
                        Q = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "Q_prime.npy")
                    except:
                        I = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "I.npy")
                        Q = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "Q.npy")
                    
                    if markers:
                        ax0.plot(I, Q, color=color, linestyle='', marker='o', markersize=5)
                    else:
                        ax0.plot(I, Q, color=color, linestyle='-')
                    
                    ax0.axvline(0.0, linestyle='solid', color='gray')
                    ax0.axhline(0.0, linestyle='solid', color='gray')
                
                    if complex_fit_above:
                        try:
                            nu_linear = np.linspace(x_data_chan[0], x_data_chan[-1], num=200) # linear sample
                            nu_peack = np.random.normal(e['nu_r'].n, 0.001, 1000) # peak sample
                            nu = np.concatenate([nu_linear, nu_peack])
                            nu = np.sort(nu)
                            Z = S_21(nu, e['Re[a]'].n, e['Im[a]'].n, e['Q_tot'].n, e['Q_c'].n, e['nu_r'].n, e['phi_0'].n)
                            
                            ax0.plot(np.real(Z), np.imag(Z), linestyle='solid', color=color, alpha=0.6, linewidth=4)
                        except:
                            pass
                
            handles.append(Line2D([0], [0], label=target.label, color=color))
            
    # plot ms2034b VNA sweeps
    if ms2034b_data_list != None:
        for i,sweep in enumerate(ms2034b_data_list):
            amp = 20*np.log10(np.sqrt(sweep.ReS21**2.0 + sweep.ImS21**2.0))
            
            color = cmap(i/(len(ms2034b_data_list)-1))
            ax0.plot(sweep.freqs, amp, color=color, linestyle='dotted', linewidth=1)
            handles.append(Line2D([0], [0], linestyle='dotted', label=sweep.label, color=color))
    
    ax0.legend(loc='best', handles=handles, fontsize=8)
    
    plt.show()


        
'''
            lowpass_cosine function
'''
def lowpass_cosine(y, tau, f_3db, width, padd_data=True):
    
    import numpy as nm
        # padd_data = True means we are going to symmetric copies of the data to the start and stop
    # to reduce/eliminate the discontinuities at the start and stop of a dataset due to filtering
    #
    # False means we're going to have transients at the start and stop of the data

    # kill the last data point if y has an odd length
    if nm.mod(len(y),2):
        y = y[0:-1]

    # add the weird padd
    # so, make a backwards copy of the data, then the data, then another backwards copy of the data
    if padd_data:
        y = nm.append( nm.append(nm.flipud(y),y) , nm.flipud(y) )

    # take the FFT
        import scipy
        import scipy.fftpack
    ffty=scipy.fftpack.fft(y)
    ffty=scipy.fftpack.fftshift(ffty)

    # make the companion frequency array
    delta = 1.0/(len(y)*tau)
    nyquist = 1.0/(2.0*tau)
    freq = nm.arange(-nyquist,nyquist,delta)
    # turn this into a positive frequency array
    pos_freq = freq[(len(ffty)//2):]

    # make the transfer function for the first half of the data
    i_f_3db = min( nm.where(pos_freq >= f_3db)[0] )
    f_min = f_3db - (width/2.0)
    i_f_min = min( nm.where(pos_freq >= f_min)[0] )
    f_max = f_3db + (width/2);
    i_f_max = min( nm.where(pos_freq >= f_max)[0] )

    transfer_function = nm.zeros(int(len(y)//2))
    transfer_function[0:i_f_min] = 1
    transfer_function[i_f_min:i_f_max] = (1 + nm.sin(-nm.pi * ((freq[i_f_min:i_f_max] - freq[i_f_3db])/width)))/2.0
    transfer_function[i_f_max:(len(freq)//2)] = 0

    # symmetrize this to be [0 0 0 ... .8 .9 1 1 1 1 1 1 1 1 .9 .8 ... 0 0 0] to match the FFT
    transfer_function = nm.append(nm.flipud(transfer_function),transfer_function)

    # apply the filter, undo the fft shift, and invert the fft
    filtered=nm.real(scipy.fftpack.ifft(scipy.fftpack.ifftshift(ffty*transfer_function)))

    # remove the padd, if we applied it
    if padd_data:
        filtered = filtered[(len(y)//3):(2*(len(y)//3))]

    # return the filtered data
    return filtered

'''
            asymmetric_least_squares_smoothing function
'''
def asymmetric_least_squares_smoothing(data, lam, p, N_iter):
    '''
    lam: adjusting parameter
    p: asymmetry parameter
    N_iter: number of iteration
    '''
    from scipy.sparse import spdiags, linalg, diags
    L = len(data)
    D = diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np.ones(L)
    for i in range(N_iter):
        W = spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = linalg.spsolve(Z, w*data)
        w = p*(data < z) + (1.0-p)*(data >= z)
    return z

'''
            adaptive_iteratively_reweighted_penalized_least_squares_smoothing function
'''
def adaptive_iteratively_reweighted_penalized_least_squares_smoothing(data, lam, N_iter):
    '''
    lam: adjusting parameter
    N_iter: number of iteration
    '''
    from scipy.sparse import spdiags, linalg, diags
    from scipy.linalg import norm
    L = len(data)
    D = diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np.ones(L)
    for i in range(N_iter):
        W = spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = linalg.spsolve(Z, w*data)
        d_mod = norm((z-data)[z>data])
        if d_mod < 0.001 * norm(data):
            return z
        p = np.exp(i*(data-z)/d_mod)
        w = 0.0*(data < z) + p*(data >= z)
    return z
        
        
'''
            buildS21Dataset function
'''      
def buildS21Dataset(sweep, ROACH='MISTRAL'):
    '''
    This function converts the raw data sweep into different .npy files
    each for a single channel

    Parameters
    ----------
    sweep : Target object or VNA object
        A Target object or VNA object to be converted.
    ROACH : string, optional
        If 'MISTRAL' the tone frequencies are computed as freqs = LO + bb,
        if 'OLIMPO' the tone frequencies are computed as freqs = LO/2 + bb. 
        The default is 'MISTRAL'.

    Returns
    -------
    None.

    '''
    from  tqdm import tqdm
    
    filename = sweep.filename
    
    if  (datapaths.target / filename).exists():
        sweep_path = datapaths.target / filename
        output_path = datapaths.target_S21 / filename
    elif (datapaths.vna / filename).exists():
        sweep_path = datapaths.vna / filename
        output_path = datapaths.vna_S21 / filename
    else:
        print("Sweep file not found.")
        sys.exit()
    
    # make the S21 directory
    if not output_path.exists():
        output_path.mkdir()
        
    print("\nBuilding S21 dataset from raw data sweep...")
    
    if ROACH == "OLIMPO":
        print("\tOLIMPO ROACH selected (freqs = LO/2 + bb)")
    elif ROACH == "MISTRAL":
        print("\tMISTRAL ROACH selected (freqs = LO + bb)")
    
    try:
        LO_freqs = np.loadtxt(sweep_path / "sweep_freqs.dat")
    except FileNotFoundError:
        LO_freqs = np.load(sweep_path / "sweep_freqs.npy")
        pass
    try:
        bb_freqs = np.loadtxt(sweep_path / "bb_freqs.dat")
    except FileNotFoundError:
        bb_freqs = np.load(sweep_path / "bb_freqs.npy")
    
    n_res = bb_freqs.size
    n_files = LO_freqs.size
    
    I = np.zeros([n_res, n_files])
    Q = np.zeros([n_res, n_files])
    
    pbar = tqdm(LO_freqs, position=0, leave=True)
    for i_file, LO_freq in enumerate(pbar):
        pbar.set_description("\tReading LO data... ")
        I_all = np.load(sweep_path / ("I"+str(LO_freq)+".npy"))
        Q_all = np.load(sweep_path / ("Q"+str(LO_freq)+".npy"))
    
        for i_res in range(n_res):
            I[i_res][i_file] = I_all[i_res]
            Q[i_res][i_file] = Q_all[i_res]
    
    # data from mistral client to compute amplitude in dBm
    accumulation_length = 2**21
    fft_length = 1024
    pbar = tqdm(range(n_res), position=0, leave=True)
    for i in pbar:
        pbar.set_description("\tComputing frequencies... ")
        
        Q[i] *= (0.5*fft_length) / ((2**31-1)*(accumulation_length-1))
        I[i] *= (0.5*fft_length) / ((2**31-1)*(accumulation_length-1))
        
        mag = np.sqrt(Q[i]*Q[i] + I[i]*I[i])
        #mag /= (2**31-1)    # mistral client
        #mag /= (accumulation_length-1)/(0.5*fft_length)    # mistral client
        mag = 20*np.log10(mag)
        phase = np.arctan2(Q[i], I[i])
    
        if ROACH == "OLIMPO":
            freqs = LO_freqs*0.5+bb_freqs[i]
        elif ROACH == "MISTRAL":
            freqs = LO_freqs+bb_freqs[i]
    
        if not (output_path / "{:03d}".format(i)).exists():
            (output_path / "{:03d}".format(i)).mkdir()
        np.save(output_path / "{:03d}".format(i) / "I.npy", I[i])
        np.save(output_path / "{:03d}".format(i) / "Q.npy", Q[i])
        np.save(output_path / "{:03d}".format(i) / "mag.npy", mag)
        np.save(output_path / "{:03d}".format(i) / "phase.npy", phase)
        np.save(output_path / "{:03d}".format(i) / "freqs.npy", freqs/1e6)




'''
            S_21 function
'''
def S_21(nu, Rea, Ima, Q_tot, Q_c, nu_r, phi_0, tau=0.0):
    '''
    This function returns the S21 scattering parameter.

    Parameters
    ----------
    nu : float, usually is a numpy array or list of number
        Frequency in [Hz].
    Rea : float
        Real part of the amplitude parameter.
    Ima : float
        Imaginary part of the amplitude parameter.
    Q_tot : float
        Total quality factor.
    Q_c : float
        Coupling quality factor.
    nu_r : float
        Resonant frequency in [Hz].
    phi_0 : number, between [0; 2pi]
        Phase parameter decribing the rotation of the resonant circle.
    tau : float, optional
        Time delay due to the transmission line length in [s].
        The default is 0.0

    Returns
    -------
    comlex number or comlex numpy array
        The computed S21 scattering parameter.

    '''
    a = Rea + Ima*1j
    return np.exp(-1j*2.0*np.pi*tau*nu)*a*(1.0 - (Q_tot/Q_c)*np.exp(1j*phi_0) / (1+2*1j*Q_tot*(nu-nu_r)/nu_r) )



'''
            jointTargetSweeps function
'''
def jointTargetSweeps(targets, exclude_channels=[[]], flat_at_0db=False):
    '''
    This function generates a plot with jointed target sweeps
    
    Parameters
    ----------
    targets : list of target objecs
        A list of target objects to be jointed.
    exclude_channels : list of list with channel indexes, optional
        Channel lists to be excluded in the plot. The default is [[]].
    flat_at_0db : boolean, optional
        If True, each resonance will be vertically shifted such that 
        its maximum is 0dB. The default is False.

    Returns
    -------
    None.

    '''
    from matplotlib import pyplot as plt
    fig = plt.figure()
    fig.set_size_inches(6, 4)
    
    for t,exclude_channels in zip(targets,exclude_channels):
        # how many channels?
        n_chan = list(range(t.target_freqs_new.size))
        for c in exclude_channels:
            n_chan.remove(c)
        
        for chan in n_chan:
            # read one sweep at a time
            x_data_chan = np.load(datapaths.target_S21 / t.filename / "{:03d}".format(chan) / "freqs.npy")
            y_data_chan = np.load(datapaths.target_S21 / t.filename / "{:03d}".format(chan) / "mag.npy")
            
            if flat_at_0db == True:
                y_data_chan -= max(y_data_chan)
            
            plt.plot(x_data_chan, y_data_chan, linewidth=1)
        
    plt.grid(linestyle='-', alpha=0.5)
    plt.ylabel('Mag [dB]')
    plt.xlabel('Frequency [MHz]')
    
    plt.show()
    #fig.savefig('vna.png', dpi=300)
        


'''
                complexS21Fit function
'''
def complexS21Fit(I, Q, freqs, output_path, RESFREQ=None, DATAPOINTS=100, verbose=False, force_emcee=False):
    '''
    Returns the complex fit of the S21 transfer function

    Parameters
    ----------
    I : TYPE
        DESCRIPTION.
    Q : TYPE
        DESCRIPTION.
    freqs : TYPE
        DESCRIPTION.
    output_path : TYPE
        DESCRIPTION.
    RESFREQ : TYPE, optional
        DESCRIPTION. The default is None.
    DATAPOINTS : TYPE, optional
        DESCRIPTION. The default is 100.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    force_emcee : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    from lmfit import Parameters, minimize, fit_report
    from uncertainties import ufloat
    
    # DATAPOINTS: number of left and right datapoints with respect to the resonance frequency 
    # to which the complex fit is computed
    
    # normalization of I and Q to max()
    IQ_MAX = np.max([np.abs(I).max(), np.abs(Q).max()])
    I /= IQ_MAX
    Q /= IQ_MAX
    
    from pathlib import Path
    output_path = Path(output_path)
    
    if not output_path.exists():
        print("\n{:s} does not exists.".format(output_path.as_posix()))
        return 0
    
    if verbose:
        print("\nFit S21 complex function for {:s}...".format(output_path.as_posix()))
    
    A = np.sqrt(I**2.0 + Q**2.0)
    phase = np.arctan2(Q, I)
    
    # a stupid guess of the resonant frequency if not given
    if RESFREQ != None:
        try:
            ARG_RESFREQ = np.argwhere(freqs >= RESFREQ)[0][0]
            RESFREQ = freqs[ARG_RESFREQ]
        except IndexError: 
            print('Given value of RESFREQ is out of bounds. Setting res_feq to min(A).')
    else:
        #ARG_RESFREQ = np.argmin(A[int(A.size*0.33):int(A.size*0.66)]) + int(A.size*0.33)
        ARG_RESFREQ = np.argmin(A)
        RESFREQ = freqs[ARG_RESFREQ]
    
    ARG_MIN = ARG_RESFREQ-int(0.5*DATAPOINTS)
    ARG_MAX = ARG_RESFREQ+int(0.5*DATAPOINTS)
    # check if there are enough datapoints at the left of the resonant frequency
    if ARG_MIN < 0:
        ARG_MIN = 0
    # check if there are enough datapoints at the right of the resonant frequency
    if ARG_MAX >= A.size:
        ARG_MAX = A.size-1
        
    
    # removing the cable delay
    tau = 0.05 # microsec
    phase += 2.0*np.pi*tau*freqs
    phase -= int(phase[0]/np.pi)*np.pi
    phase = np.unwrap(phase)
    I = A*np.cos(phase)
    Q = A*np.sin(phase)
    
    # errors on I and Q data is 1%
    IErr = I*0.01
    QErr = Q*0.01
    
    # compute the center coordinates by averaging max and min data values
    I_m, I_M = I[ARG_MIN:ARG_MAX].min(), I[ARG_MIN:ARG_MAX].max()
    Q_m, Q_M = Q[ARG_MIN:ARG_MAX].min(), Q[ARG_MIN:ARG_MAX].max()
    Xc = (I_m+I_M)*0.5
    Yc = (Q_m+Q_M)*0.5
    
    # compute the radius by averaging the 'x radius' and the 'y radius'
    _rx = 0.5*np.abs(I_M-I_m)
    _ry = 0.5*np.abs(Q_M-Q_m)
    radius = 0.5*(_rx+_ry)
    
    if verbose:
        print("\n\tcenterd circle parameters")
        print("\tXc = ", Xc)
        print("\tYc = ", Yc)
        print("\tradius = ", radius)
    
    # translation
    I_rot = I-Xc
    Q_rot = Q-Yc
    
    # rotation
    I_mid = 0.5*(I_rot[ARG_MAX]+I_rot[ARG_MIN])
    Q_mid = 0.5*(Q_rot[ARG_MAX]+Q_rot[ARG_MIN])
    rotAngle = np.pi-np.arctan2(Q_mid, I_mid)
    
    if verbose:
        print("\n\trotation angle")
        print("\tangle = ", rotAngle)
    
    file = open(output_path / "log.txt","w")
    file.write("centerd circle parameters\n")
    file.write("Xc = {:.3e}\n".format(Xc))
    file.write("Yc = {:.3e}\n".format(Yc))
    file.write("radius = {:.3e}\n".format(radius))
    file.write("\n")
    file.write("rotation angle\n")
    file.write("angle = {:.3e}\n".format(rotAngle))
    file.write("\n")
    file.close()
    
    I_second = np.cos(rotAngle)*I_rot-np.sin(rotAngle)*Q_rot
    Q_second = np.sin(rotAngle)*I_rot+np.cos(rotAngle)*Q_rot
    
    
    # determine the phase and amplitude (prime)
    ph_second = np.arctan2(Q_second, I_second)
    mag_second = np.sqrt(I_second**2 + Q_second**2)
    
    
    # phase fit
    def phaseFunction(nu, Q_tot, nu_r, y0):
        return -2*np.arctan2((2.0*Q_tot*(nu/nu_r-1.0)), 1.0) - y0
    
    def phaseFunctionResiduals(params, nu, data, uncertainty):
        Q_tot = params['Q_tot']
        nu_r = params['nu_r']
        y0 = params['y0']
        return (data-phaseFunction(nu, Q_tot, nu_r, y0))/uncertainty
        
    params = Parameters()
    params.add('Q_tot', value=5000, min=0)
    params.add('nu_r', value=RESFREQ, min=RESFREQ*0.8, max=RESFREQ*1.2)
    params.add('y0', value=-ph_second[ARG_RESFREQ], min=-5.0*np.pi, max=5.0*np.pi)
    
    x_data = freqs[ARG_MIN:ARG_MAX]
    y_data = ph_second[ARG_MIN:ARG_MAX] # offset here!
    uncertainty = 0.1*ph_second[ARG_MIN:ARG_MAX]
    
    out = minimize(phaseFunctionResiduals, params, args=(x_data, y_data, uncertainty), method='leastsq', max_nfev=1000)
    
    if verbose:
        print("\n\tphase fit")
        print("\t"+fit_report(out))
    file = open(output_path / "log.txt","a")
    file.write(fit_report(out))
    file.close()
    
    Q_totPhFit = out.params['Q_tot'].value
    nu_rPhFit = out.params['nu_r'].value
    y0_ph_fit = out.params['y0'].value
    
    # compute Qc and Phi0 from fit parameters
    # guess the value of a
    Ima = 0.5*(Q[ARG_MAX]+Q[ARG_MIN])
    Rea = 0.5*(I[ARG_MAX]+I[ARG_MIN])
    a_norm = np.sqrt(Ima**2.0 + Rea**2.0)
    Q_c = 0.5*a_norm*Q_totPhFit/radius
    try:
        Q_i = Q_totPhFit*Q_c/(Q_c-Q_totPhFit)
        if Q_i<0.0:
            Q_i = 1.0e7 # this should be a sufficiently high value
    except ZeroDivisionError:
        Q_i = 1.0e7 # this should be a sufficiently high value
        pass
        
    # let's compute phi0
    alpha_angle = np.angle(Rea+1j*Ima)
    beta_angle = np.angle(Xc-Rea+1j*(Yc-Ima))
    phi_0 = np.pi+beta_angle-alpha_angle
    if phi_0<0.0:
        phi_0 += 2.0*np.pi
    if phi_0>2.0*np.pi:
        phi_0 %= (2.0*np.pi)
    
    # 6-parameters complex fit (fine tuning)
    # complex residuals over errorbars
    def complexResiduals(params, freq, z_val, z_err=None):
        Rea = params['Rea']
        Ima = params['Ima']
        Q_tot = params['Q_tot']
        Q_c = params['Q_c']
        nu_r = params['nu_r']
        phi_0 = params['phi_0']
        tau = params['tau']
        if z_err is None:
            return np.abs(z_val - S_21(freq, Rea, Ima, Q_tot, Q_c, nu_r, phi_0, tau)) 
        else:
            return np.abs( (z_val - S_21(freq, Rea, Ima, Q_tot, Q_c, nu_r, phi_0, tau))/z_err) 
    
    # chi2 minimization
    params = Parameters()
    params.add('Rea', value=Rea, min=-2.0, max=2.0)
    params.add('Ima', value=Ima, min=-2.0, max=2.0)
    params.add('Q_c', value=Q_c, min=0.0, max=1.0e12)
    params.add('Q_i', value=Q_i, min=0.0, max=1.0e12)
    params.add('Q_tot', expr='(Q_c*Q_i)/(Q_c+Q_i)')
    params.add('nu_r', value=nu_rPhFit, min=nu_rPhFit*0.8, max=nu_rPhFit*1.2)
    params.add('phi_0', value=phi_0, min=0.0, max=2.0*np.pi)
    params.add('tau', value=0.0, vary=False)
    
    z_data = I[ARG_MIN:ARG_MAX]+1j*Q[ARG_MIN:ARG_MAX]
    z_err = IErr[ARG_MIN:ARG_MAX]+1j*QErr[ARG_MIN:ARG_MAX]
    freqs = freqs[ARG_MIN:ARG_MAX]
    
    # try a leastsquare fit otherwise an MCMC
    try:
        out = minimize(complexResiduals, params, args=(freqs, z_data, z_err), method='leastsq', max_nfev=1000)
        if verbose:
            print("\n\tcomplex fit results")
            print("\t"+fit_report(out))
        file = open(output_path / "log.txt","a")
        file.write(fit_report(out))
        file.close()
        Rea = ufloat(out.params['Rea'].value, out.params['Rea'].stderr)*IQ_MAX
        Ima = ufloat(out.params['Ima'].value, out.params['Ima'].stderr)*IQ_MAX
        Q_c = ufloat(out.params['Q_c'].value, out.params['Q_c'].stderr)
        Q_i = ufloat(out.params['Q_i'].value, out.params['Q_i'].stderr)
        Q_tot = ufloat(out.params['Q_tot'].value, out.params['Q_tot'].stderr)
        nu_r = ufloat(out.params['nu_r'].value, out.params['nu_r'].stderr)
        phi_0 = ufloat(out.params['phi_0'].value, out.params['phi_0'].stderr)
    except: 
        if force_emcee:
            emcee_kws = dict(steps=5500, burn=500, nwalkers=100, progress=True)
            out = minimize(complexResiduals, params, args=(freqs, z_data, z_err), method='emcee', **emcee_kws)
            if verbose:
                print("\n\tcomplex fit results")
                print("\t"+fit_report(out))
            file = open(output_path / "log.txt","a")
            file.write(fit_report(out))
            file.close()
            Rea = ufloat(out.params['Rea'].value, out.params['Rea'].stderr)
            Ima = ufloat(out.params['Ima'].value, out.params['Ima'].stderr)
            Q_tot = ufloat(out.params['Q_tot'].value, out.params['Q_tot'].stderr)
            Q_c = ufloat(out.params['Q_c'].value, out.params['Q_c'].stderr)
            nu_r = ufloat(out.params['nu_r'].value, out.params['nu_r'].stderr)
            phi_0 = ufloat(out.params['phi_0'].value, out.params['phi_0'].stderr)
        pass
    
    
    '''
    # compute the Q_i value
    Q_i = Q_c*Q_tot/(Q_c-Q_tot)
    
    if verbose:
        print("\nQ_i: {:uf}".format(Q_i))
    file = open(output_path / "log.txt","a")
    file.write("\nQ_i: {:uf}".format(Q_i))
    file.close()
    '''

    # save data into .npy files
    np.save(output_path / "I_prime.npy", I*IQ_MAX)
    np.save(output_path / "Q_prime.npy", Q*IQ_MAX)
    np.save(output_path / "phase_prime.npy", phase)
    np.save(output_path / "I_second.npy", I_second*IQ_MAX)
    np.save(output_path / "Q_second.npy", Q_second*IQ_MAX)
    np.save(output_path / "mag_second.npy", mag_second)
    np.save(output_path / "phase_second.npy", ph_second)
    np.save(output_path / "transformation_parameters.npy", arr=[Xc, Yc, rotAngle])
    np.save(output_path / "phase_fit_parameters.npy", arr=[radius, Q_totPhFit, nu_rPhFit, y0_ph_fit])
    np.save(output_path / "reduced_chi2.npy", arr=out.redchi)
    np.save(output_path / "complex_parameters.npy", arr=[Rea.n, Rea.s, Ima.n, Ima.s, Q_tot.n, Q_tot.s, Q_c.n, Q_c.s, 
                                                         Q_i.n, Q_i.s, nu_r.n, nu_r.s, phi_0.n, phi_0.s, tau])
    np.save(output_path / "fit_interval_extrema.npy", arr=[ARG_MIN, ARG_MAX])
    
    return {'Re[a]': Rea, 'Im[a]': Ima, 'Q_tot': Q_tot, 'Q_c': Q_c, 'Q_i': Q_i, 'nu_r': nu_r, 'phi_0': phi_0, 'tau': tau}, out.redchi


'''
            complexS21Plot function
'''
def complexS21Plot(complex_fit_data_path):
    from pathlib import Path
    complex_fit_data_path = Path(complex_fit_data_path)
    
    I = np.load(complex_fit_data_path / "I.npy")
    Q = np.load(complex_fit_data_path / "Q.npy")
    mag = np.sqrt(Q*Q + I*I)
    mag = 20*np.log10(mag)
    
    #mag = np.load(complex_fit_data_path / "mag.npy")
    #ph = np.load(output_path / "phi.npy")
    freqs = np.load(complex_fit_data_path / "freqs.npy")
        
    A = np.sqrt(I**2 + Q**2)
    
    #IErr = I*0.01
    #QErr = Q*0.01
    
    I_prime = np.load(complex_fit_data_path / "I_prime.npy")
    Q_prime = np.load(complex_fit_data_path / "Q_prime.npy")
    phase_prime = np.load(complex_fit_data_path / "phase_prime.npy")
    phase_prime = np.unwrap(phase_prime)
    
    I_second = np.load(complex_fit_data_path / "I_second.npy")
    Q_second = np.load(complex_fit_data_path / "Q_second.npy")
    #mag_second = np.load(output_path / "mag_second.npy")
    phase_second = np.load(complex_fit_data_path / "phase_second.npy")
    #if phase_second[0] <= 0.0:
    #    phase_second += 2.0*np.pi
    #phase_second = np.unwrap(phase_second)
    
    try:
        [Xc, Yc, rotAngle] = np.load(complex_fit_data_path / "transformation_parameters.npy")
    except FileNotFoundError:
        print("Phase fit parameters not found!")
        pass
    try:
        [radius, Q_phfit, nu_r_ph_fit, y0_ph_fit] = np.load(complex_fit_data_path / "phase_fit_parameters.npy")
        [Rea, ReaErr, Ima, ImaErr, Qt, QtErr, Qc, QcErr, Qi, QiErr, nu_r, nu_rErr, phi0, phi0Err, tau] = np.load(complex_fit_data_path / "complex_parameters.npy")
    except ValueError:
        [radius, Q_phfit, nu_r_ph_fit] = np.load(complex_fit_data_path / "phase_fit_parameters.npy")
        y0_ph_fit = 0.0
        pass
    except FileNotFoundError:
        print("S21 fit parameters not found!")
        pass
    try:
        [ARG_MIN, ARG_MAX] = np.load(complex_fit_data_path / "fit_interval_extrema.npy")
        ARG_MIN = int(ARG_MIN)
        ARG_MAX = int(ARG_MAX)
    except FileNotFoundError:
        pass
    
    def phaseFunction(nu, Q, nu_r, y0):
        return -2*np.arctan2((2.0*Q*(nu/nu_r-1.0)), 1.0) - y0
    
    # PLOT
    from matplotlib import pyplot as plt
    # create a figure
    # text settings
    fig = plt.figure()
    fig.set_size_inches(15, 7)
    plt.subplots_adjust(left=0.07, right=0.93, top=0.93, bottom=0.07)
    plt.subplots_adjust(hspace=0.0)
    circlePlot = plt.subplot(121, aspect='equal')
    phasePlot1 = plt.subplot(322)
    phasePlot2 = plt.subplot(324)
    amplitudePlot = plt.subplot(326)
    
    color_raw = (17/255, 157/255, 164/255, 1.0)
    color_raw_alpha = (17/255, 157/255, 164/255, 0.3)
    color_notau = (0, 0, 0, 1.0)
    color_notau_alpha = (0, 0, 0, 0.3)
    color_centered = (0, 210/255, 249/255, 1.0)
    color_centered_alpha = (0, 210/255, 249/255, 0.3)
    color_fit_interval = (250/255, 253/255, 15/255)
    
    # axis
    circlePlot.axvline(x=0.0, ymin=-10.0, ymax=10.0, linewidth='1.0', linestyle='--', color='black', alpha=0.3)
    circlePlot.axhline(y=0.0, xmin=-10.0, xmax=10.0, linewidth='1.0', linestyle='--', color='black', alpha=0.3)
    
    # raw IQ data
    circlePlot.plot(I, Q, color=color_raw,linestyle='-', marker='o', markersize=4, markerfacecolor=color_raw_alpha)
    # centerd data
    circlePlot.plot(I_second, Q_second, color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha)
    phasePlot1.plot(freqs, phase_second, color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha)
    
    # IQ data to be fitted
    circlePlot.plot(I_prime, Q_prime, color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)
    phasePlot2.plot(freqs, phase_prime, color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)
    amplitudePlot.plot(freqs, mag-mag[0], color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)
    
    # plot of the fit interval points
    try:
        circlePlot.plot(I_prime[ARG_MIN:ARG_MAX], Q_prime[ARG_MIN:ARG_MAX], color=color_fit_interval, linestyle='', marker='o', markersize=2, markerfacecolor=color_fit_interval)
        phasePlot2.plot(freqs[ARG_MIN:ARG_MAX], phase_prime[ARG_MIN:ARG_MAX], color=color_fit_interval, linestyle='', marker='o', markersize=2, markerfacecolor=color_fit_interval)
        amplitudePlot.plot(freqs[ARG_MIN:ARG_MAX], mag[ARG_MIN:ARG_MAX]-mag[0], color=color_fit_interval, linestyle='', marker='o', markersize=2, markerfacecolor=color_fit_interval)
    except:
        pass
    
    # phase fit
    try:
        phasePlot1.plot(freqs, phaseFunction(freqs, Q_phfit, nu_r_ph_fit, y0_ph_fit), linestyle='-', color='blue', alpha=0.3, linewidth=3.0)#, label='Phase fit')
    except:
        pass
    # Complex fit
    try:
        Z_freqs = np.linspace(freqs[0], freqs[-1], num=2000)
        Z = S_21(Z_freqs, Rea, Ima, Qt, Qc, nu_r, phi0)
        circlePlot.plot(np.real(Z), np.imag(Z), linestyle='-', color='red', alpha=0.5, linewidth=3.0)
        phasePlot2.plot(Z_freqs, np.unwrap(np.angle(Z)), linestyle='-', color='red', alpha=0.5, linewidth=3.0)
        amplitudePlot.plot(Z_freqs, 20*np.log10(np.abs(Z))-mag[0], linestyle='-', color='red', alpha=0.5, linewidth=3.0)
    except:
        pass
    
    
    handle = [Line2D([0], [0], color=color_raw,linestyle='-', marker='o', markersize=4, markerfacecolor=color_raw_alpha, label='Raw data'),
              Line2D([0], [0], color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha, label=r'$\tau$ removed'),
              Line2D([0], [0], color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha, label='Centered data'),
              #Line2D([0], [0], linestyle='-', color='blue', alpha=0.3, linewidth=3.0, label='Circle Fit'),
              Line2D([0], [0], linestyle='-', color='blue', alpha=0.3, linewidth=3.0, label='Phase fit'),
              Line2D([0], [0], linestyle='-', color='red', alpha=0.5, linewidth=3.0, label='S$_{21}$ Fit')]
    
    # graphics
    circlePlot.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
    circlePlot.set_xlabel('$I$')
    circlePlot.set_ylabel('$Q$')
    #circlePlot.legend(loc='best')
    circlePlot.grid(alpha=0.5)
    circlePlot.yaxis.set_ticks_position('both')
    circlePlot.xaxis.set_ticks_position('both')
    circlePlot.minorticks_on()
    circlePlot.yaxis.set_tick_params(direction='in', which='both')
    circlePlot.xaxis.set_tick_params(direction='in', which='both')
    
    phasePlot1.legend(loc='upper right', handles=handle)
    phasePlot1.grid(axis="x")
    phasePlot1.set_ylabel('Phase [rad]')
    phasePlot1.set_xlabel('Frequency [MHz]')
    phasePlot1.grid(alpha=0.5)
    phasePlot1.yaxis.set_ticks_position('both')
    phasePlot1.xaxis.set_ticks_position('both')
    phasePlot1.minorticks_on()
    phasePlot1.yaxis.set_tick_params(direction='in', which='both')
    phasePlot1.xaxis.set_tick_params(direction='in', which='both')
    #plt.setp(phasePlot1.get_xticklabels(), visible=False)
    
    #phasePlot2.legend(loc='best')
    phasePlot2.grid(axis="x")
    phasePlot2.set_ylabel('Phase [rad]')
    phasePlot2.set_xlabel('Frequency [MHz]')
    phasePlot2.grid(alpha=0.5)
    phasePlot2.yaxis.set_ticks_position('both')
    phasePlot2.xaxis.set_ticks_position('both')
    phasePlot2.minorticks_on()
    phasePlot2.yaxis.set_tick_params(direction='in', which='both')
    phasePlot2.xaxis.set_tick_params(direction='in', which='both')
    #plt.setp(phasePlot2.get_xticklabels(), visible=False)
    
    #amplitudePlot.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #amplitudePlot.legend(loc='best')
    amplitudePlot.set_xlabel('Frequency [MHz]')
    amplitudePlot.set_ylabel(r'Amplitude [dB]')
    amplitudePlot.grid(alpha=0.5)
    amplitudePlot.yaxis.set_ticks_position('both')
    amplitudePlot.xaxis.set_ticks_position('both')
    amplitudePlot.minorticks_on()
    amplitudePlot.yaxis.set_tick_params(direction='in', which='both')
    amplitudePlot.xaxis.set_tick_params(direction='in', which='both')
    
    # print the figures
    plt.show()



'''
            Plot picolog data
'''
def plotPicolog(filename, OFFSET_TIME=False):
    '''
    This function returns a plot of a given picolog file with three channels

    Parameters
    ----------
    filename : string
        The picolog filename.
    OFFSET_TIME : boolean, optional
        If True, ch1 (usually time) will start with a zero. The default is 
        False.

    Returns
    -------
    None.

    '''
    import matplotlib.pyplot as plt
    from G31_thermometry import Thermometer
    thermometer = Thermometer(model='DT670', serial_no='D6068043')
    
    # read data
    ch1,ch2,ch3 = np.loadtxt(datapaths.picolog/filename, dtype=float, skiprows=1, unpack=True, delimiter=',')
    
    if OFFSET_TIME:
        ch1 -= ch1[0]
    
    # convert thermometer voltage to temperature
    temperature = thermometer.temperature(ch2/1000.)
    
    fig = plt.figure(figsize=(9,9))
    ax0 = fig.add_subplot(111)
    ax0.plot(ch1, temperature, c='black')
    ax0.set_ylabel('T$^{Cu}$ [K]', color='black')
    ax0.tick_params(axis='y', colors='black')
    ax0.set_xlim([ch1[0], ch1[-1]])
    ax0.xaxis.grid()
    
    ax0_ = fig.add_subplot(111, sharex = ax0, frameon = False)
    ax0_.fill_between(ch1, ch3, 0.0, alpha=0.20, color='red')
    ax0_.yaxis.set_label_position("right")
    ax0_.tick_params(axis='y', colors='red')
    ax0_.yaxis.tick_right()
    delta_ch3 = ch3.max() - ch3.min()
    ax0_.set_ylim([ch3.min()-(delta_ch3*0.2), ch3.max()+(delta_ch3*0.2)])
    ax0_.set_ylabel('Heater [V]', color='red')
    
    plt.show()
    
    


'''
        
'''
def plot_target(target, axis, color=None, linestyle='solid', linewidth=1, flat_at_0db=False):
    for e in target.entry:
        # read one sweep at a time
        channel = e['channel']
        x_data_chan = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "freqs.npy")
        y_data_chan = np.load(datapaths.target_S21 / target.filename / "{:03d}".format(channel) / "mag.npy")
        
        if flat_at_0db:
            y_data_chan -= max(y_data_chan)
        
        if e['is_out_of_res']:
            axis.plot(x_data_chan, y_data_chan, linewidth=linewidth, alpha=0.2, color=color)
        else:
            axis.plot(x_data_chan, y_data_chan, linewidth=linewidth, linestyle=linestyle, color=color)
    


def Delta(T_c):
    from scipy.constants import k as kb
    return 1.764*kb*T_c

def n_qp(T, T_c, N_0):
    from scipy.constants import k as kb
    from uncertainties import unumpy as unp
    return 2.0*N_0*unp.sqrt(2.0*np.pi*kb*T*Delta(T_c)) * unp.exp(-Delta(T_c)/(kb*T))

def electrical_phase_responsivity_linear_fit(nu_r, base_nu_r, T, T_c, N_0, V_abs, label=None, color='black', axis=None):
    '''
    Returns the slopes of the dx vs dN_qp trend

    Parameters
    ----------
    nu_r : numpy array
        resonant frequencies in MHz.
    base_nu_r : float
        base resonant frequency in MHz.
    T : numpy array
        sweep temperatures in K.
    T_c : float
        critical temperature in K.
    N_0 : float
        density of states at the Fermi surface in 1/um^3 1/J.
    V_abs : float
        absorber volume in um^3.
    label : list, optional
        list of labels. The default is 'pixel'.
    color : list, optional
        list of colors. The default is 'black'.
    axis : matplotlib.pyplot axis, optional
        axis for plot. The default is None.

    Returns
    -------
    lmfit.params
        dictionary of fit paramters.

    '''
    from lmfit import Minimizer, Parameters
    from uncertainties import ufloat, unumpy as unp
    
    # converts lists to np arrays
    nu0 = ufloat(base_nu_r.n, base_nu_r.s)
    delta_x = [(n-nu0)/nu0 for n in nu_r]
    
    T = np.asarray(T)
    
    # function used for the fit procedure
    def linear_function(x, slope, intercept):
        return x*slope + intercept
    
    def fcn2min(params, x, data, errs=None):
        slope = params['slope']
        intercept = params['intercept']
        
        model = linear_function(x, slope, intercept)
        
        if errs is None:
            return data-model
        
        return (data-model)/errs
    
    # quasiparticle number density
    N_qp = V_abs * n_qp(T, T_c, N_0)
    
    # linear fit
    params = Parameters()
    params.add('slope', value=(min(delta_x).n-max(delta_x).n)/(max(unp.nominal_values(N_qp))-min(unp.nominal_values(N_qp))), min=0, max=-np.inf)
    params.add('intercept', value=0.0, min=-1, max=1)
    
    try:
        result = Minimizer(fcn2min, params, fcn_args=(unp.nominal_values(N_qp), [d.n for d in delta_x], [d.s for d in delta_x])).minimize(method='least_squares')
        if axis is not None:
            axis.ticklabel_format(axis='both', style='sci', useMathText=True, scilimits=(0,0))
            axis.tick_params(axis='both', which='both', direction='in', bottom=True, left=True, top=True, right=True)
            
            axis.set_xlabel("$N_{qp}$")
            axis.set_ylabel("$\delta x$")
            
            axis.errorbar(unp.nominal_values(N_qp), [d.n for d in delta_x], yerr=[d.s for d in delta_x], xerr=unp.std_devs(N_qp), color=color, linestyle='', fmt='o', capsize=2, markersize=3)
            axis.plot(unp.nominal_values(N_qp), linear_function(unp.nominal_values(N_qp), result.params['slope'].value, result.params['intercept'].value), color=color, label=label)
            
            axis.grid(color='gray', alpha=0.4)
            
            return result.params
    except ValueError:
        print("Cannot find fit parameters.")
        return
    
    
    

def electrical_amplitude_responsivity_linear_fit(Q_i, T, T_c, N_0, V_abs, label=None, color='black', axis=None):
    '''
    Returns the slopes of the dx vs dN_qp trend

    Parameters
    ----------
    Q_i : numpy array
        internal quality factors.
    T : numpy array
        sweep temperatures in K.
    T_c : float
        critical temperature in K.
    N_0 : float
        density of states at the Fermi surface in 1/um^3 1/J.
    V_abs : float
        absorber volume in um^3.
    label : list, optional
        list of labels. The default is 'pixel'.
    color : list, optional
        list of colors. The default is 'black'.
    axis : matplotlib.pyplot axis, optional
        axis for plot. The default is None.

    Returns
    -------
    lmfit.params
        dictionary of fit paramters.

    '''
    from lmfit import Minimizer, Parameters
    from uncertainties import unumpy as unp
    
    # converts lists to np arrays
    uQ_i = unp.uarray([x.n for x in Q_i], [x.s for x in Q_i])
    one_over_Q_i = 1.0/uQ_i
    
    T = np.asarray(T)
    
    # function used for the fit procedure
    def linear_function(x, slope, intercept):
        return x*slope + intercept
    
    def fcn2min(params, x, data, errs=None):
        slope = params['slope']
        intercept = params['intercept']
        
        model = linear_function(x, slope, intercept)
        
        if errs is None:
            return data-model
        
        return (data-model)/errs
    
    # quasiparticle number density
    N_qp = V_abs * n_qp(T, T_c, N_0)
    
    # linear fit
    params = Parameters()
    params.add('slope', value=(max(unp.nominal_values(one_over_Q_i))-min(unp.nominal_values(one_over_Q_i)))/(max(unp.nominal_values(N_qp))-min(unp.nominal_values(N_qp))), min=0, max=np.inf)
    params.add('intercept', value=0.0, min=-np.inf, max=np.inf)
    
    try:
        result = Minimizer(fcn2min, params, fcn_args=(unp.nominal_values(N_qp), unp.nominal_values(one_over_Q_i),  unp.std_devs(one_over_Q_i))).minimize(method='least_squares')
        if axis is not None:
            axis.ticklabel_format(axis='both', style='sci', useMathText=True, scilimits=(0,0))
            axis.tick_params(axis='both', which='both', direction='in', bottom=True, left=True, top=True, right=True)
            
            axis.set_xlabel("$N_{qp}$")
            axis.set_ylabel("$1/Q_i$")
            
            axis.errorbar(unp.nominal_values(N_qp), [Q.n for Q in one_over_Q_i], yerr=[Q.s for Q in one_over_Q_i], xerr=unp.std_devs(N_qp), color=color, linestyle='', fmt='o', capsize=2, markersize=3)
            axis.plot(unp.nominal_values(N_qp), linear_function(unp.nominal_values(N_qp), result.params['slope'].value, result.params['intercept'].value), color=color, label=label)
            
            axis.grid(color='gray', alpha=0.4)
            
            return result.params
    except ValueError:
        print("Cannot find fit parameters.")
        return    




def electrical_phase_responsivity(slope, T_c, Q_tot, tau_qp, eta_qp):
    """
    This functions returns the electrical phase responsivity

    Parameters
    ----------
    slope : (u)float
        The slope of the delta nu_r / nu_r0 vs N_qp linear trend.
    T_c : (u)float
        Critical temperature in Kelvin.
    Q_tot : (u)float
        Total quality factor.
    tau_qp : (u)float
        Quasiparticle recombination time in seconds.
    eta_qp : float
        Pair breaking efficiency.

    Returns
    -------
    (u)float
        The electrical phase responsivity in rad/Watt.

    """
    return -slope*4.0*Q_tot*eta_qp*tau_qp/Delta(T_c) # rad/W


def electrical_amplitude_responsivity(slope, T_c, Q_tot, tau_qp, eta_qp):
    """
    This functions returns the electrical amplitude responsivity

    Parameters
    ----------
    slope : (u)float
        The slope of the 1/Q_i vs N_qp linear trend.
    T_c : (u)float
        Critical temperature in Kelvin.
    Q_tot : (u)float
        Total quality factor.
    tau_qp : (u)float
        Quasiparticle recombination time in seconds.
    eta_qp : float
        Pair breaking efficiency.

    Returns
    -------
    (u)float
        The electrical amplitude responsivity in rad/Watt.

    """
    return slope*Q_tot*eta_qp*tau_qp/Delta(T_c) # rad/W



'''
noise_type = 'chp' or 'cha'
'''
def plot_extracted_stream(filename, channel, axis, noise_type='chp', label=None, color=None, linestyle='solid', linewidth=1):
    path = datapaths.output_noise / filename

    if label is None:
        label = noise_type+'_{:03d}'.format(channel)

    # open time streams
    time_stream = np.load(path / 'time_stream.npy')
    ch = np.load(path / (noise_type+'_{:03d}.npy'.format(channel)) )

    axis.plot(time_stream, ch, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    axis.grid(color="gray", alpha=0.5)
    axis.set_xlabel("Time [s]")
    axis.set_ylabel(noise_type)


'''
noise_type = 'chp' or 'cha'
'''
def plot_extracted_noise_spectral_density(filename, channel, sampling_frequency, axis=None, noise_type='chp', label=None, color=None, linestyle='solid', linewidth=1, ylim=None):
    from scipy.signal import periodogram
    
    path = datapaths.output_noise / filename
    
    if label is None:
        label = noise_type+'_{:03d}'.format(channel)

    # open time streams
    #time_stream = np.load(path / 'time_stream.npy')
    ch = np.load(path / (noise_type+'_{:03d}.npy'.format(channel)) )
    
    # when scaling = "density" the function returns the power spectral density in 1/Hz
    freqs, noise_spectral_density = periodogram(ch, fs=sampling_frequency, window="han", scaling="density")

    if axis is not None:
        axis.plot(freqs, np.sqrt(noise_spectral_density), color=color, label=label, linestyle=linestyle, linewidth=linewidth, alpha=0.6)
        axis.tick_params(axis='both', which='both', direction='in', bottom=True, left=True, top=True, right=True)
        axis.grid(color="gray", alpha=0.5)
        axis.set_xlabel("Frequency [Hz]")
        axis.set_xlim([freqs[1], freqs[-1]])
        axis.set_xscale('log')
        axis.set_yscale('log')
        axis.set_ylabel("Noise spectral density [1/$\sqrt{Hz}$]")
        
        if ylim is not None:
            axis.set_ylim(ylim)
    
    return freqs, noise_spectral_density
    
    
    
'''

'''
def plot_S21_MS2034B(ms2034b_object, axis, label=None, linestyle='solid', linewidth=1, color='black'):
    if label is None:
        label = ms2034b_object.label

    if ms2034b_object.mode == 'real_imag':
        S21 = np.sqrt(ms2034b_object.ReS21**2.0 + ms2034b_object.ImS21**2.0)
    if ms2034b_object.mode == 'log_mag_phase':
        S21 = ms2034b_object.S21DB

    axis.plot(ms2034b_object.freqs, S21, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    axis.grid(color="gray", alpha=0.5)
    axis.set_xlabel("Frequency [MHz]")
    axis.set_ylabel("$|S_{21}|$ [dBm]")



def electrical_phase_noise_equivalent_power(responsivity, noise_spectral_density, freqs, axis=None, label='ch', color='black', linestyle='solid', linewidth=1):

    NEP = np.sqrt(noise_spectral_density)/responsivity  

    axis.plot(freqs, NEP, color=color, label=label, linestyle=linestyle, linewidth=linewidth)
    axis.grid(color="gray", alpha=0.5)
    axis.set_xlabel("Frequency [Hz]")
    axis.set_xlim([freqs[0], freqs[-1]])
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_ylabel("NEP$_{el}^{ph}$ [W/$\sqrt{Hz}$]")
    
    return freqs, NEP



def psd(ydata, fs, window='hann', nfft=None):
    from scipy.signal import periodogram, welch
    import numpy as np
    freqs, spectrum = periodogram(ydata, fs=fs, scaling="density", window=window, nfft=nfft)
    #freqs, spectrum = welch(ch, fs=fs, window='hann', scaling='density')
    psd = np.sqrt(spectrum)
    return freqs, psd



def filt(ydata, Wn, fs, order=2, filter_type='highpass'):
    from scipy.signal import butter, sosfilt
    sos = butter(N=order, Wn=Wn, btype=filter_type, fs=fs, output='sos')
    filtered = sosfilt(sos, ydata)
    return filtered




'''
            lsTarget function
'''
def lsTarget():
    '''
    This function lists all the target sweep directories

    Returns
    -------
    target_list : list of strings
        A list of strings with the name of target sweep directories.
    target_S21_list : list of strings
        A list of strings with the name of converted target sweep directories.

    '''
    import os
    print('List of target sweeps:')
    target_list = sorted([t for t in os.listdir(datapaths.target) if t[0] != '.'])
    for t in target_list:
        print(t)
    print('')
    print('List of target sweeps converted into different channels:')
    target_S21_list = sorted([t for t in os.listdir(datapaths.target_S21) if t[0] != '.'])
    for t in target_S21_list:
        print(t)
    print('')
    return target_list, target_S21_list
    

'''
            lsVNA function
'''
def lsVNA():
    '''
    This function lists all the VNA sweep directories
    
    Returns
    -------
    vna_list : list of strings
        A list of strings with the name of VNA sweep directories.
    vna_S21_list : list of strings
        A list of strings with the name of converted VNA sweep directories.

    '''
    import os
    print('List of VNA sweeps:')
    vna_list = sorted([vna for vna in os.listdir(datapaths.vna) if vna[0] != '.'])
    for vna in vna_list:
        print(vna)
    print('')
    print('List of VNA sweeps converted into different channels:')
    vna_S21_list = sorted([vna for vna in os.listdir(datapaths.vna_S21) if vna[0] != '.'])
    for vna in vna_S21_list:
        print(vna)
    print('') 
    return vna_list, vna_S21_list
    