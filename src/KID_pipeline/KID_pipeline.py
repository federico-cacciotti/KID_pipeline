import numpy as np
from matplotlib.lines import Line2D
import sys

# check if data directories exist
import paths
paths.check_if_dirs_exist(mkdir=True)

import Target
import VNA

# plot configuration
#plt.rc('text', usetex=False)                        #
#plt.rc('font', family='serif')
#plt.rcParams['text.usetex'] = False                 #
#plt.rcParams['text.latex.preamble'] = r"\usepackage{amssymb} \usepackage{siunitx}"
#plt.rcParams.update({'figure.max_open_warning': 0})
#plt.rcParams["axes.formatter.use_mathtext"] = True
#plt.minorticks_on()
#plt.tick_params(axis='both', which='both', direction='in', labelleft='on', labelright='off', bottom=True, top=True, left=True, right=True)
#plt.xkcd()
#plt.close('all')


#NEP_photon = 5.e-15 #W/Hz^1/2
#NEP_photon_noatm = 8.e-17 #W/Hz^1/2
        
        
class readMS2034B():
    def __init__(self, filename):
        self.filename = filename
        [self.freqs, self.ReS11, self.ImS11, self.ReS21, self.ImS21, 
         self.ReS12, self.ImS12, self.ReS22, self.ImS22] = np.loadtxt(fname=paths.anritsuMS2034B / filename, dtype=float, skiprows=23, unpack=True)
        self.freqs *= 1e3 # GHz to MHz
    
    
    def plotS21(self):
        from matplotlib import pyplot as plt
        fig = plt.figure()
        fig.set_size_inches(6, 12)
        ax0 = plt.subplot(311)
        ax1 = plt.subplot(312)
        ax2 = plt.subplot(313)
        
        amp = np.sqrt(self.ReS21**2 + self.ImS21**2)
        ph = np.arctan2(self.ImS21, self.ReS21)
        ph = np.unwrap(ph)
        
        ax0.plot(self.freqs, amp, color='black', linewidth=1)
        ax1.plot(self.freqs, ph, color='black', linewidth=1)
        ax2.plot(self.ReS21, self.ImS21, color='black', linewidth=1)
        
        ax0.yaxis.set_ticks_position('both')
        ax0.xaxis.set_ticks_position('both')
        ax0.minorticks_on()
        ax0.yaxis.set_tick_params(direction='in', which='both')
        ax0.xaxis.set_tick_params(direction='in', which='both')
        ax0.grid(linestyle='-', alpha=0.5)
        ax0.set_ylabel('Mag [dB]')
        ax0.set_xlabel('Frequency [MHz]')
        
        ax1.yaxis.set_ticks_position('both')
        ax1.xaxis.set_ticks_position('both')
        ax1.minorticks_on()
        ax1.yaxis.set_tick_params(direction='in', which='both')
        ax1.xaxis.set_tick_params(direction='in', which='both')
        ax1.grid(linestyle='-', alpha=0.5)
        ax1.set_ylabel('Phase [rad]')
        ax1.set_xlabel('Frequency [MHz]')
        
        ax2.set_aspect('equal')
        ax2.yaxis.set_ticks_position('both')
        ax2.xaxis.set_ticks_position('both')
        ax2.minorticks_on()
        ax2.yaxis.set_tick_params(direction='in', which='both')
        ax2.xaxis.set_tick_params(direction='in', which='both')
        ax2.grid(linestyle='-', alpha=0.5)
        ax2.set_ylabel(r'Im[S_{{21}}]')
        ax2.set_xlabel(r'Re[S_{{21}}]')
        
        plt.show()
        
        
        
        
        
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
        
        
        
def buildS21Dataset(sweep, ROACH='MISTRAL'):
    from  tqdm import tqdm
    
    filename = sweep.filename
    
    if  (paths.target / filename).exists():
        sweep_path = paths.target / filename
        output_path = paths.target_S21 / filename
    elif (paths.vna / filename).exists():
        sweep_path = paths.vna / filename
        output_path = paths.vna_S21 / filename
    else:
        print("Sweep file not found.")
        sys.exit()
    
    # make the S21 directory
    if not output_path.exists():
        output_path.mkdir()
        
    print("\nBuilding S21 dataset from raw data sweep...")
    
    if ROACH == "OLIMPO":
        print("OLIMPO ROACH selected (freqs = LO/2 + bb)")
    elif ROACH == "MISTRAL":
        print("MISTRAL ROACH selected (freqs = LO + bb)")
    
    LO_freqs = np.load(sweep_path / "sweep_freqs.npy")
    bb_freqs = np.load(sweep_path / "bb_freqs.npy")
    
    n_res = bb_freqs.size
    n_files = LO_freqs.size
    
    I = np.zeros([n_res, n_files])
    Q = np.zeros([n_res, n_files])
    
    pbar = tqdm(LO_freqs, position=0, leave=True)
    for i_file, LO_freq in enumerate(pbar):
        pbar.set_description("Reading LO data... ")
        I_all = np.load(sweep_path / ("I"+str(LO_freq)+".npy"))
        Q_all = np.load(sweep_path / ("Q"+str(LO_freq)+".npy"))
    
        for i_res in range(n_res):
            I[i_res][i_file] = I_all[i_res]
            Q[i_res][i_file] = Q_all[i_res]
    
    # data from mistral client to compute amplitude in dBm
    accumulation_length = 2**21
    fft_len = 1024
    pbar = tqdm(range(n_res), position=0, leave=True)
    for i in pbar:
        pbar.set_description("Computing frequencies... ")
        mag = np.sqrt(Q[i]*Q[i] + I[i]*I[i])
        mag /= (2**31-1)    # mistral client
        mag /= (accumulation_length-1)/(0.5*fft_len)    # mistral client
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





def S_21(nu, Rea, Ima, Q_tot, Q_c, nu_r, phi_0, tau):
        a = Rea + Ima*1j
        return np.exp(1j*2.0*np.pi*tau*nu)*a*(1.0 - (Q_tot/Q_c)*np.exp(1j*phi_0) / (1+2*1j*Q_tot*(nu-nu_r)/nu_r) )






def jointTargetSweeps(targets, exclude_channels=[[]], flat_at_0db=False):
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
            x_data_chan = np.load(paths.target_S21 / t.filename / "{:03d}".format(chan) / "freqs.npy")
            y_data_chan = np.load(paths.target_S21 / t.filename / "{:03d}".format(chan) / "mag.npy")
            
            if flat_at_0db == True:
                y_data_chan -= max(y_data_chan)
            
            plt.plot(x_data_chan, y_data_chan, linewidth=1)
        
    plt.grid(linestyle='-', alpha=0.5)
    plt.ylabel('Mag [dB]')
    plt.xlabel('Frequency [MHz]')
    
    plt.show()
    #fig.savefig('vna.png', dpi=300)
        





def complexS21Fit(I, Q, freqs, res_freq, output_path, DATAPOINTS=100, verbose=False):
    from lmfit import Parameters, minimize, fit_report
    from uncertainties import ufloat
    
    # DATAPOINTS: number of left and right datapoints with respect to the resonance frequency 
    # to which the complex fit is computed
    
    from pathlib import Path
    output_path = Path(output_path)
    
    if not output_path.exists():
        print("\n{:s} does not exists.".format(output_path.as_posix()))
        return 0
    
    if verbose:
        print("\nFit S21 complex function for {:s}...".format(output_path.as_posix()))
    
    A = np.sqrt(I**2 + Q**2)
    phase = np.arctan2(Q, I)
    
    ARG_RESFREQ = np.argmin(A)
    ARG_MIN = ARG_RESFREQ-DATAPOINTS
    ARG_MAX = ARG_RESFREQ+DATAPOINTS
    # check if there are enough datapoints at the left of the resonance frequency
    if ARG_MIN < 0:
        ARG_MIN = 0
    # check if there are enough datapoints at the right of the resonance frequency
    if ARG_MAX >= A.size:
        ARG_MAX = A.size-1
        
    
    # removing the cable delay
    tau = 0.04 # microsec
    phase += 2.0*np.pi*tau*freqs
    phase -= int(phase[0]/np.pi)*np.pi
    phase = np.unwrap(phase)
    I = A*np.cos(phase)
    Q = A*np.sin(phase)
    
    # errors on I and Q data is 1%
    IErr = I*0.01
    QErr = Q*0.01
    
    # compute the center coordinates by averaging max and min data values
    I_m, I_M = I.min(), I.max()
    Q_m, Q_M = Q.min(), Q.max()
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
    I_mid = 0.5*(I_rot[-1]+I_rot[0])
    Q_mid = 0.5*(Q_rot[-1]+Q_rot[0])
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
    def phaseFunction(nu, Q_tot, nu_r):
        return -2*np.arctan2((2.0*Q_tot*(nu/nu_r-1.0)), 1.0)
    
    def phaseFunctionResiduals(params, nu, data, uncertainty):
        Q_tot = params['Q_tot']
        nu_r = params['nu_r']
        return (data-phaseFunction(nu, Q_tot, nu_r))/uncertainty
        
    params = Parameters()
    params.add('Q_tot', value=5000, min=0)
    params.add('nu_r', value=res_freq, min=res_freq*0.99, max=res_freq*1.01)
    
    x_data = freqs[ARG_MIN:ARG_MAX]
    y_data = ph_second[ARG_MIN:ARG_MAX]
    uncertainty = 0.01*ph_second[ARG_MIN:ARG_MAX]
    
    out = minimize(phaseFunctionResiduals, params, args=(x_data, y_data, uncertainty), method='leastsq')
    
    if verbose:
        print("\n\tphase fit")
        print("\t"+fit_report(out))
    file = open(output_path / "log.txt","a")
    file.write(fit_report(out))
    file.close()
    
    Q_totPhFit = out.params['Q_tot'].value
    nu_rPhFit = out.params['nu_r'].value
    
    # compute Qc and Phi0 from fit parameters
    # guess the value of a
    Ima = 0.5*(Q[-1]+Q[0])
    Rea = 0.5*(I[-1]+I[0])
    a_norm = np.sqrt(Ima**2.0 + Rea**2.0)
    Q_c = 0.5*a_norm*Q_totPhFit/radius
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
    def complexResiduals(params, freq, z_val, z_err):
        Rea = params['Rea']
        Ima = params['Ima']
        Q_tot = params['Q_tot']
        Q_c = params['Q_c']
        nu_r = params['nu_r']
        phi_0 = params['phi_0']
        tau = params['tau']
        return np.abs( (z_val - S_21(freq, Rea, Ima, Q_tot, Q_c, nu_r, phi_0, tau))/z_err )
    
    # chi2 minimization
    params = Parameters()
    params.add('Rea', value=Rea)
    params.add('Ima', value=Ima)
    params.add('Q_tot', value=Q_totPhFit, min=0)
    params.add('Q_c', value=Q_c, min=0)
    params.add('nu_r', value=nu_rPhFit, min=nu_rPhFit*0.99, max=nu_rPhFit*1.01)
    params.add('phi_0', value=phi_0, min=0.0, max=2.0*np.pi)
    params.add('tau', value=tau, vary=False)
    
    z_data = I[ARG_MIN:ARG_MAX]+1j*Q[ARG_MIN:ARG_MAX]
    z_err = IErr[ARG_MIN:ARG_MAX]+1j*QErr[ARG_MIN:ARG_MAX]
    freqs = freqs[ARG_MIN:ARG_MAX]
    
    
    out = minimize(complexResiduals, params, args=(freqs, z_data, z_err), method='leastsq')
    
    if verbose:
        print("\n\tcomplex fit results")
        print("\t"+fit_report(out))
    file = open(output_path / "log.txt","a")
    file.write(fit_report(out))
    file.close()
    
    if out.message != 'Fit succeeded.':
        return None
    
    Rea = ufloat(out.params['Rea'].value, out.params['Rea'].stderr)
    Ima = ufloat(out.params['Ima'].value, out.params['Ima'].stderr)
    Q_tot = ufloat(out.params['Q_tot'].value, out.params['Q_tot'].stderr)
    Q_c = ufloat(out.params['Q_c'].value, out.params['Q_c'].stderr)
    nu_r = ufloat(out.params['nu_r'].value, out.params['nu_r'].stderr)
    phi_0 = ufloat(out.params['phi_0'].value, out.params['phi_0'].stderr)

    # compute the Q_i value
    Q_i = Q_c*Q_tot/(Q_c-Q_tot)
    
    if verbose:
        print("\nQ_i: {:uf}".format(Q_i))
    file = open(output_path / "log.txt","a")
    file.write("\nQ_i: {:uf}".format(Q_i))
    file.close()

    # save data into .npy files
    np.save(output_path / "I_prime.npy", I)
    np.save(output_path / "Q_prime.npy", Q)
    np.save(output_path / "phase_prime.npy", phase)
    np.save(output_path / "I_second.npy", I_second)
    np.save(output_path / "Q_second.npy", Q_second)
    np.save(output_path / "mag_second.npy", mag_second)
    np.save(output_path / "phase_second.npy", ph_second)
    np.save(output_path / "transformation_parameters.npy", arr=[Xc, Yc, rotAngle])
    np.save(output_path / "phase_fit_parameters.npy", arr=[radius, Q_totPhFit, nu_rPhFit])
    np.save(output_path / "reduced_chi2.npy", arr=out.redchi)
    np.save(output_path / "complex_parameters.npy", arr=[Rea.n, Rea.s, Ima.n, Ima.s, Q_tot.n, Q_tot.s, Q_c.n, Q_c.s, 
                                                         Q_i.n, Q_i.s, nu_r.n, nu_r.s, phi_0.n, phi_0.s, tau])
    
    return {'Re[a]': Rea, 'Im[a]': Ima, 'Q_tot': Q_tot, 'Q_c': Q_c, 'Q_i': Q_i, 'nu_r': nu_r, 'phi_0': phi_0, 'tau': tau}, out.redchi



def complexS21Plot(complex_fit_data_path):
    from pathlib import Path
    complex_fit_data_path = Path(complex_fit_data_path)
    
    I = np.load(complex_fit_data_path / "I.npy")
    Q = np.load(complex_fit_data_path / "Q.npy")
    #mag = np.load(output_path / "mag.npy")
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
    if phase_second[0] <= 0.0:
        phase_second += 2.0*np.pi
    phase_second = np.unwrap(phase_second)
    
    [Xc, Yc, rotAngle] = np.load(complex_fit_data_path / "transformation_parameters.npy")
    [radius, Q_phfit, nu_r_ph_fit] = np.load(complex_fit_data_path / "phase_fit_parameters.npy")
    [Rea, ReaErr, Ima, ImaErr, Qt, QtErr, Qc, QcErr, Qi, QiErr, nu_r, nu_rErr, phi0, phi0Err, tau] = np.load(complex_fit_data_path / "complex_parameters.npy")
    
    def phaseFunction(nu, Q, nu_r):
        return -2*np.arctan2((2.0*Q*(nu/nu_r-1.0)), 1.0)
    
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
    
    # raw IQ data plot
    circlePlot.plot(I, Q, color=color_raw,linestyle='-', marker='o', markersize=4, markerfacecolor=color_raw_alpha)#, label='Raw data')
    circlePlot.plot(I_prime, Q_prime, color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)#, label=r'$\tau$ removed')
    #circlePlot.errorbar(I, Q, xerr=IErr, yerr=QErr, marker='.', linestyle='--', linewidth=1.0, markersize=1.0, color='green', alpha=1.0, label='Raw data')
    #circlePlot.errorbar(I_prime, Q_prime, xerr=IErr, yerr=QErr, marker='.', linestyle='--', linewidth=1.0, markersize=1.0, color='black', alpha=1.0, label='Tau removed')
    
    # Complex plot plot
    Z = S_21(freqs, Rea, Ima, Qt, Qc, nu_r, phi0, tau)
    circlePlot.plot(np.real(Z), np.imag(Z), linestyle='-', color='red', alpha=0.5, linewidth=3.0)#, label='S$_{21}$ Fit')
    
    # resonance after translation and rotation
    #circlePlot.plot(I_second, Q_second, marker='.', markersize=1.0, color='black', alpha=0.5, label='Centered data')
    circlePlot.plot(I_second, Q_second, color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha)#, label='Centered data')
    # circle fit plot
    #phiTrain = np.linspace(0.0, 2.0*np.pi, 100)
    #circlePlot.plot(radius*np.cos(phiTrain), radius*np.sin(phiTrain), linestyle='-', color='blue', alpha=0.3, linewidth=3.0)#, label='Circle Fit')
    circlePlot.axvline(x=0.0, ymin=-10.0, ymax=10.0, linewidth='1.0', linestyle='--', color='black', alpha=0.3)
    circlePlot.axhline(y=0.0, xmin=-10.0, xmax=10.0, linewidth='1.0', linestyle='--', color='black', alpha=0.3)
    
    # phase plot
    phasePlot1.plot(freqs, phase_second, color=color_centered,linestyle='-', marker='o', markersize=4, markerfacecolor=color_centered_alpha)#, label='Centered data')
    #phasePlot1.plot(freqs, ph_second, marker='.', markersize=1.0, color=color_centered, alpha=0.5, label='Phase')
    phasePlot1.plot(freqs, phaseFunction(freqs, Q_phfit, nu_r_ph_fit), linestyle='-', color='blue', alpha=0.3, linewidth=3.0)#, label='Phase fit')
    
    #phasePlot2.plot(freqs, ph_prime, marker='.', markersize=1.0, color=color_notau, alpha=0.5, label='Phase')
    phasePlot2.plot(freqs, phase_prime, color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)#, label='Phase')
    phasePlot2.plot(freqs, np.unwrap(np.angle(Z)), linestyle='-', color='red', alpha=0.5, linewidth=3.0)#, label='S$_{21}$ Fit')
    
    # amplitude plot
    #amplitudePlot.plot(freqs, A, marker='.', markersize=1.0, color=color_notau, alpha=0.5, label='Amplitude')
    amplitudePlot.plot(freqs, A/1e6, color=color_notau,linestyle='-', marker='o', markersize=4, markerfacecolor=color_notau_alpha)#, label='Amplitude')
    amplitudePlot.plot(freqs, np.abs(Z)/1e6, linestyle='-', color='red', alpha=0.5, linewidth=3.0)#, label='S$_{21}$ Fit')
    
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
    
    amplitudePlot.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #amplitudePlot.legend(loc='best')
    amplitudePlot.set_xlabel('Frequency [MHz]')
    amplitudePlot.set_ylabel(r'Amplitude $\times 10^{{6}}$[au]')
    amplitudePlot.grid(alpha=0.5)
    amplitudePlot.yaxis.set_ticks_position('both')
    amplitudePlot.xaxis.set_ticks_position('both')
    amplitudePlot.minorticks_on()
    amplitudePlot.yaxis.set_tick_params(direction='in', which='both')
    amplitudePlot.xaxis.set_tick_params(direction='in', which='both')
    
    # print the figures
    plt.show()

