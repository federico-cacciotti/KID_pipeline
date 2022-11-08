import numpy as np
from pathlib import Path

paths_cosmic_rays = Path("/Users/federicocacciotti/Documents/PhD/dati/cosmic_rays")


'''
        Impulse analysis for quasiparticle recombination time measurements
'''
    

def impulse_response(t, tau_rise, tau_fall, A, t_offset, v_offset):
    y = v_offset + A*(np.exp(-(t-t_offset)/tau_fall) - np.exp(-(t-t_offset)/tau_rise))
    y = [v_offset if y_<=v_offset else y_ for y_ in y]
    return y

def double_exponential(t, tau_rise, tau_fall, tau_therm, A1, A2, t_offset):
    y = A1*(np.exp(-(t-t_offset)/tau_fall) - np.exp(-(t-t_offset)/tau_rise)) + A2*np.exp(-(t-t_offset)/tau_therm)
    y = [0.0 if y_<0 else y_ for y_ in y]
    return y

def exponential_decay(t, A, tau, t_offset):
    return A*np.exp(-(t-t_offset)/tau)


class Event():
    
    def __init__(self, filename, CH1=None, CH2=None, CH3=None, CH4=None, label='Event'):
        
        self.label = label
        self.time = {'label': 'TIME', 'data': []}
        self.channels = {}
        
        for i,CHi in enumerate([CH1, CH2, CH3, CH4]):
            if CHi != None:
                channel_data = {'label': CHi, 'data': []}
                self.channels['CH{:d}'.format(i+1)] = channel_data
                if CHi == 'Q':
                    self.Q_channel = 'CH{:d}'.format(i+1)
                if CHi == 'I':
                    self.I_channel = 'CH{:d}'.format(i+1)
        
        # read data from .csv file
        n_channels = len(self.channels)
        data = np.loadtxt(paths_cosmic_rays/filename, dtype=float, delimiter=',', skiprows=21, unpack=True, usecols=range(n_channels+1))
        self.vertical_units = str(np.loadtxt(paths_cosmic_rays/filename, dtype=str, delimiter=',', skiprows=12, usecols=1, max_rows=1))
        self.horizontal_units = str(np.loadtxt(paths_cosmic_rays/filename, dtype=str, delimiter=',', skiprows=5, usecols=1, max_rows=1))
        
        self.time['data'] = data[0]
        self.time['data'] -= self.time['data'][0]   # remove time offset
        
        j = 1
        for CH in ['CH1', 'CH2', 'CH3', 'CH4']:
            try:
                self.channels[CH]['data'] = data[j]
                j += 1
            except:
                pass
            
        
    def plot(self):
        import matplotlib.pyplot as plt
        
        colors = ['gray', 'red', 'blue', 'green']
        
        fig = plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        
        for CHi, color in zip(self.channels, colors):
            ax0.plot(self.time['data'], self.channels[str(CHi)]['data'], linestyle='solid', label=self.channels[str(CHi)]['label'], color=color)
    
        # compute amplitude and phase of the signal
        amplitude = np.sqrt(self.channels[self.I_channel]['data']**2. + self.channels[self.Q_channel]['data']**2.)
        phase = np.arctan2(self.channels[self.Q_channel]['data']+5, self.channels[self.I_channel]['data']+5)
        phase = np.unwrap(phase)
        
        # plot the amplitude on the same y axis
        ax0.plot(self.time['data'], amplitude, color='orange', label='Amplitude', linestyle='solid')
    
        ax0.legend(loc='best')
        ax0.grid(color='gray', alpha=0.4)
        ax0.set_xlabel('Time ['+self.horizontal_units+']')
        ax0.set_ylabel('Voltage ['+self.vertical_units+']')
        plt.title(self.label)
        
        # plot the phase on a different y axis
        ax1 = ax0.twinx()
        ax1.plot(self.time['data'], phase, color='violet', linestyle='solid', label='Phase')
        ax1.set_ylabel('Phase [rad]')
        
        plt.show()
        
    
    def fit(self):
        '''
        This function computes the magnitude of the modulated IQ signal and 
        performs the best fit with a double exponential function

        Returns
        -------
        list
            DESCRIPTION.

        '''
        from scipy.signal import correlate, correlation_lags, medfilt
        from lmfit import Minimizer, Parameters, report_fit
        from uncertainties import ufloat
        
        
        self.A = np.sqrt(self.channels[self.I_channel]['data']**2. + self.channels[self.Q_channel]['data']**2.)
        
        A_filt = medfilt(self.A, kernel_size=51)
            
        # normalization
        self.A /= max(A_filt)
        
        # fit
        # errorbars are the std dev of the initial plateau
        keep = np.array(self.A) >= 0.15
        std_dev = np.std(self.A[keep])
        v_offset = np.average(self.A[keep])
        
        keep = np.array(self.A) >= 0.0 # all the stream
        #keep = np.array(self.A) >= 0.15*max(self.A)
        #keep = np.array([True if i%10==0 else False for i,k in enumerate(keep)])
        
        errors = np.ones(len(keep))*std_dev
        
        def fcn2min(params, x, data, errs=None):
            model = impulse_response(x, params['tau_rise'].value, params['tau_fall'].value, params['A'].value, params['t_offset'].value, params['v_offset'].value)
            
            if errs is None:
                return model-data
            
            return (model-data)/errs
        
        # parameters initial guess
        mask = self.A >= 0.5#0.36
        t_offset = self.time['data'][mask][0]
        t_fall_offset = self.time['data'][mask][-1]
        
        mask = self.A >= 0.95
        time_at_max = 0.5*(self.time['data'][mask][0] + self.time['data'][mask][-1])
        
        params = Parameters()
        params.add('tau_rise', value=time_at_max-t_offset, min=0)
        params.add('tau_fall', value=t_fall_offset-time_at_max, min=0)
        params.add('A', value=2., min=0)
        params.add('t_offset', value=t_offset)
        params.add('v_offset', value=v_offset)
        
        try:
            result = Minimizer(fcn2min, params, fcn_args=(self.time['data'][keep], self.A[keep], errors)).minimize()
            self.par = result.params
            
            
            # print data
            report_fit(result)
        except:
            print("Optimal parameters not found")
    
    
    def plot_fit(self):
        import matplotlib.pyplot as plt
        
        fig = plt.figure(figsize=(5, 5))
        ax0 = plt.subplot(111)
        
        ax0.plot(self.time['data']*1e6, self.A, linestyle='solid', label='Data', color='gray')
        
        FIR = impulse_response(self.time['data']*1e6, self.par['tau_rise'].value*1e6, self.par['tau_fall'].value*1e6, self.par['A'].value, self.par['t_offset'].value*1e6, self.par['v_offset'].value)
        ax0.plot(self.time['data']*1e6, FIR, linestyle='solid', linewidth=4, label='Fit', color='red', alpha=0.5)
        
        ax0.legend(loc='best')
        ax0.grid(color='gray', alpha=0.4)
        ax0.set_xlabel('Time [$\mu$'+self.horizontal_units+']')
        ax0.set_ylabel('Amplitude')
        plt.title(self.label)
        
        plt.show()
        


def fitAverage(IQ_timestreams, bounds=[[0, 0, 0, -50], [400, 400, 100, 10]], p0=[100, 150, 10, 0], time_unit='us', voffset=False):
    from scipy.signal import correlate, correlation_lags, medfilt
    from scipy.optimize import curve_fit
    from uncertainties import ufloat
    
    amplitudes = []
    amplitudes_filt = []
    times = []
    amps_max = []
    
    # read data
    for i,file in enumerate(IQ_timestreams):
        print(file)
        t, chI, chQ = np.loadtxt(paths_cosmic_rays / file, float, delimiter=',', skiprows=21, unpack=True)
        
        if time_unit == 'ms':
            t *= 1e3 # in millisec
        else:
            t *= 1e6 # in microsec
        t -= t[0]
        
        if np.abs(chI.max()) < np.abs(chI.min()):
            chI = -chI
        if np.abs(chQ.max()) < np.abs(chQ.min()):
            chQ = -chQ
        A = np.sqrt((chI+1)**2. + (chQ+1)**2.) - np.sqrt(2)
        
        A_filt = medfilt(A, kernel_size=151)
        
        amp_max_filt = max(A_filt)
        
        A /= amp_max_filt
        A_filt /= amp_max_filt
        
        amplitudes_filt.append(A_filt)
    
        amplitudes.append(A)
        times.append(t)
        amps_max.append(amp_max_filt)
        
    # horizontal offset
    for i,(t, A, amp_max) in enumerate(zip(times, amplitudes_filt, amps_max)):
        
        # the two datased must have the same horizontal scale!
        correlation = correlate(amplitudes[1], A)
        lag = correlation_lags(amplitudes[1].size, A.size)
        
        max_lag = lag[np.argmax(correlation)]
        
        time_offset = max_lag*(t[1]-t[0])
        
        times[i] += time_offset
    
    # merging of the data
    times = np.hstack(times)
    amplitudes = np.hstack(amplitudes)
    
    times_idx_sorted = np.flip(np.argsort(times))
    times = times[times_idx_sorted[::-1]]
    amplitudes = amplitudes[times_idx_sorted[::-1]]
    
    amplitudes_filt = medfilt(amplitudes, kernel_size=151)
    
    # vertical offset
    if voffset == True:
        keep = times > times[-1]-0.25*(times[-1]-times[0])
        average = np.average(amplitudes_filt[keep])
        amplitudes_filt -= average
        amplitudes -= average
        # second normalization
        max_amp = max(amplitudes_filt)
        amplitudes_filt /= max_amp
        amplitudes /= max_amp
    
    # guess time shift
    time_shift = times[np.argwhere(amplitudes_filt >= 0.2)[0]]
    print('Time shift = ', time_shift)
    times -= time_shift
    
    
    # fit
    # MISTRAL
    '''
    par, cov = curve_fit(impulseResponse, times[12000:100000], amplitudes[12000:100000], p0=[50, 50, 10, -10], sigma=np.ones(times[12000:100000].size)*0.1,
                             absolute_sigma=True, bounds=[[0, 0, 0, -200], [100, 100, 200, 200]])
    '''
    
    # COSMO
    par, cov = curve_fit(impulse_response, times[0:37000], amplitudes[0:37000], p0=p0, sigma=np.ones(times[0:37000].size)*0.2,
                             absolute_sigma=True, bounds=bounds)

    parErr = np.sqrt(np.diag(cov))
    
    print("Tau_rise = {:.3u} ".format(ufloat(par[0], parErr[0]))+time_unit)
    print("Tau_fall = {:.3u} ".format(ufloat(par[1], parErr[1]))+time_unit)
    print("A = {:.3u}".format(ufloat(par[2], parErr[2])))
    print("t_offset = {:.3u} ".format(ufloat(par[3], parErr[3]))+time_unit)
    
    return [times, amplitudes, amplitudes_filt, par, parErr]
    

def plotAverage(times, amplitudes, amplitudes_filt, par, parErr, title=r'$\tau_{{qp}}$', time_unit='um'):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    
    # plot
    fig = plt.figure()
    fig.set_size_inches(4, 4)
    ax0 = plt.subplot(111)
    fig.subplots_adjust(top=0.9, right=0.99, bottom=0.1, left=0.18)
    
    max_time = max(times)
    
    ax0.plot(times, amplitudes, marker='.', linestyle='', markersize=1, color='gray', alpha=0.05)
    ax0.plot(times, amplitudes_filt, marker='.', linestyle='', markersize=1, color='k', alpha=0.5)
    time_fit = np.linspace(par[3], max_time, num=500)
    ax0.plot(time_fit, impulse_response(time_fit, par[0], par[1], par[2], par[3]), 
             color='blue', linewidth=3, alpha=1)
    
    #ax0.fill_between(times, amplitudes_filt-0.1, amplitudes_filt+0.1, alpha=0.4, color='blue')
    
    labels = ['Data', 'Average', 'Fit']
    #labels = ['Data', 'Fit']
    handles = [Line2D([0], [0],color='gray', ls='', markersize=3, marker='.'),
               Line2D([0], [0],color='black', ls='', markersize=3, marker='.'),
               Line2D([0], [0],color='blue', lw=3)]
    ax0.legend(labels=labels, handles=handles, loc='upper right', prop={'size': 10})
    ax0.set_ylabel("Scaled amplitude")
    if time_unit == 'um':
        ax0.set_xlabel(r"Time [$\mu s$]")
    else:
        ax0.set_xlabel(r"Time [$ms$]")
    ax0.set_xlim([times[0], times[-1]])
    ax0.grid(alpha=0.5)
    title = title
    ax0.set_title(title)
    ax0.yaxis.set_ticks_position('both')
    ax0.xaxis.set_ticks_position('both')
    ax0.minorticks_on()
    ax0.yaxis.set_tick_params(direction='in', which='both')
    ax0.xaxis.set_tick_params(direction='in', which='both')
    
    #ax0.plot(times, amplitudes_filt-impulseResponse(times, par[0], par[1], par[2], par[3]))
    
    #np.save(paths.cosmic_rays / "138mK_377MHz.npy", [times, amplitudes, par, parErr])
    
    fig.savefig(title+'.png', dpi=330)
    
    plt.show()
    

def fitAverageDoubleExp(IQ_timestreams, idx_lim=35000, bounds=[[0, 0, 0, 0, 0, -5], [400, 400, 1200, 100, 100, 5]], p0=[100, 150, 750, 1, 0.5, 0], time_unit='us', voffset=False):
    from scipy.signal import correlate, correlation_lags, medfilt
    from scipy.optimize import curve_fit
    from uncertainties import ufloat
    
    amplitudes = []
    amplitudes_filt = []
    times = []
    amps_max = []
    
    # read data
    for i,file in enumerate(IQ_timestreams):
        print(file)
        t, chI, chQ = np.loadtxt(paths_cosmic_rays / file, float, delimiter=',', skiprows=21, unpack=True)
        
        if time_unit == 'ms':
            t *= 1e3 # in millisec
        else:
            t *= 1e6 # in microsec
        t -= t[0]
        
        if np.abs(chI.max()) < np.abs(chI.min()):
            chI = -chI
        if np.abs(chQ.max()) < np.abs(chQ.min()):
            chQ = -chQ
        A = np.sqrt((chI+1)**2. + (chQ+1)**2.) - np.sqrt(2)
        
        A_filt = medfilt(A, kernel_size=151)
        
        amp_max_filt = max(A_filt)
        
        A /= amp_max_filt
        A_filt /= amp_max_filt
        
        amplitudes_filt.append(A_filt)
    
        amplitudes.append(A)
        times.append(t)
        amps_max.append(amp_max_filt)
        
    # horizontal offset
    for i,(t, A, amp_max) in enumerate(zip(times, amplitudes_filt, amps_max)):
        
        # the two datased must have the same horizontal scale!
        correlation = correlate(amplitudes[1], A)
        lag = correlation_lags(amplitudes[1].size, A.size)
        
        max_lag = lag[np.argmax(correlation)]
        
        time_offset = max_lag*(t[1]-t[0])
        
        times[i] += time_offset
    
    # merging of the data
    times = np.hstack(times)
    amplitudes = np.hstack(amplitudes)
    
    times_idx_sorted = np.flip(np.argsort(times))
    times = times[times_idx_sorted[::-1]]
    amplitudes = amplitudes[times_idx_sorted[::-1]]
    
    amplitudes_filt = medfilt(amplitudes, kernel_size=151)
    
    # vertical offset
    if voffset == True:
        keep = times > times[-1]-0.25*(times[-1]-times[0])
        average = np.average(amplitudes_filt[keep])
        amplitudes_filt -= average
        amplitudes -= average
        # second normalization
        max_amp = max(amplitudes_filt)
        amplitudes_filt /= max_amp
        amplitudes /= max_amp
    
    # guess time shift
    time_shift = times[np.argwhere(amplitudes_filt >= 0.2)[0]]
    print('Time shift = ', time_shift)
    times -= time_shift
    
    # thermal tail fil
    parTail, covTail = curve_fit(exponential_decay, times[idx_lim:-1], amplitudes_filt[idx_lim:-1], p0=[p0[2], p0[4], p0[5]])
    print('Thermal tail fit: ', parTail)
    
    parIR, covIR = curve_fit(impulse_response, times[0:idx_lim], amplitudes[0:idx_lim], p0=[p0[0], p0[1], p0[3], p0[5]], sigma=np.ones(times[0:idx_lim].size)*0.2,
                             absolute_sigma=True)
    print('Impulse response 1: ', parIR)


    par, cov = curve_fit(double_exponential, times, amplitudes, p0=[parIR[0], parIR[1], parTail[0], parIR[2], parTail[1], parIR[3]], sigma=np.ones(times.size)*0.2,
                             absolute_sigma=True, bounds=bounds)
    parErr = np.sqrt(np.diag(cov))
    
    print("Tau_rise = {:.3u} ".format(ufloat(par[0], parErr[0]))+time_unit)
    print("Tau_fall = {:.3u} ".format(ufloat(par[1], parErr[1]))+time_unit)
    print("Tau_therm = {:.3u} ".format(ufloat(par[2], parErr[2]))+time_unit)
    print("A1 = {:.3u}".format(ufloat(par[3], parErr[3])))
    print("A2 = {:.3u}".format(ufloat(par[4], parErr[4])))
    print("t_offset = {:.3u} ".format(ufloat(par[5], parErr[5]))+time_unit)
    
    return [times, amplitudes, amplitudes_filt, par, parErr]
    

def plotAverageDoubleExp(times, amplitudes, amplitudes_filt, par, parErr, title=r'$\tau_{{qp}}$', time_unit='um'):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    
    # plot
    fig = plt.figure()
    fig.set_size_inches(4, 4)
    ax0 = plt.subplot(111)
    fig.subplots_adjust(top=0.9, right=0.99, bottom=0.1, left=0.18)
    
    max_time = max(times)
    
    ax0.plot(times, amplitudes, marker='.', linestyle='', markersize=1, color='gray', alpha=0.05)
    ax0.plot(times, amplitudes_filt, marker='.', linestyle='', markersize=1, color='k', alpha=0.5)
    time_fit = np.linspace(par[3], max_time, num=500)
    ax0.plot(time_fit, double_exponential(time_fit, *par), 
             color='blue', linewidth=3, alpha=1)
    
    #ax0.fill_between(times, amplitudes_filt-0.1, amplitudes_filt+0.1, alpha=0.4, color='blue')
    
    labels = ['Data', 'Average', 'Fit']
    #labels = ['Data', 'Fit']
    handles = [Line2D([0], [0],color='gray', ls='', markersize=3, marker='.'),
               Line2D([0], [0],color='black', ls='', markersize=3, marker='.'),
               Line2D([0], [0],color='blue', lw=3)]
    ax0.legend(labels=labels, handles=handles, loc='upper right', prop={'size': 10})
    ax0.set_ylabel("Scaled amplitude")
    if time_unit == 'um':
        ax0.set_xlabel(r"Time [$\mu s$]")
    else:
        ax0.set_xlabel(r"Time [$ms$]")
    ax0.set_xlim([times[0], times[-1]])
    ax0.grid(alpha=0.5)
    title = title
    ax0.set_title(title)
    ax0.yaxis.set_ticks_position('both')
    ax0.xaxis.set_ticks_position('both')
    ax0.minorticks_on()
    ax0.yaxis.set_tick_params(direction='in', which='both')
    ax0.xaxis.set_tick_params(direction='in', which='both')
    
    #ax0.plot(times, amplitudes_filt-impulseResponse(times, par[0], par[1], par[2], par[3]))
    
    #np.save(paths.cosmic_rays / "138mK_377MHz.npy", [times, amplitudes, par, parErr])
    
    fig.savefig(title+'.png', dpi=330)
    
    plt.show()

    
def fitEvents(IQ_timestreams):
    from scipy.signal import correlate, correlation_lags, medfilt
    from scipy.optimize import curve_fit
    
    amplitudes = []
    amplitudes_filt = []
    times = []
    amps_max = []
    pars = []
    parErrs = []
    fitFile = []
    
    # read data
    for i,file in enumerate(IQ_timestreams):
        print(file)
        t, chI, chQ = np.loadtxt(paths_cosmic_rays / file, float, delimiter=',', skiprows=21, unpack=True)
        t *= 1e6
        t -= t[0]
        A = np.sqrt((chI+1)**2. + (chQ+1)**2.) - np.sqrt(2)
        
        A_filt = medfilt(A, kernel_size=151)
        amp_max_filt = max(A_filt)
        
        # normalization
        A /= amp_max_filt
    
        amplitudes.append(A)
        times.append(t)
        amps_max.append(amp_max_filt)
        
    # horizontal offset
    for i,(t, A, amp_max) in enumerate(zip(times, amplitudes, amps_max)):
        correlation = correlate(amplitudes[2], A)
        lag = correlation_lags(amplitudes[2].size, A.size)
        max_lag = lag[np.argmax(correlation)]
        time_offset = max_lag*(t[1]-t[0])
        times[i] += time_offset
    
    for i,(t, A, amp_max) in enumerate(zip(times, amplitudes, amps_max)):
        # guess and subtract the time shift
        A_filt = medfilt(A, kernel_size=101)
        time_shift = t[np.argwhere(A_filt >= 0.2)[0]]
        print('Time shift = ', time_shift)
        t -= time_shift
    
    # fit
    for t, A, file in zip(times, amplitudes, IQ_timestreams):
        keep = np.array(A) >= 0.15*max(A)
        keep = np.array([True if i%10==0 else False for i,k in enumerate(keep)])
        try:
            par, cov = curve_fit(impulse_response, t[keep], A[keep], p0=[20, 20, 100, 0], sigma=np.ones(len(t[keep]))*0.1,
                                 absolute_sigma=True, bounds=[[1, 1, 0, -30], [300, 300, 200, 0]])
            parErr = np.sqrt(np.diag(cov))
            # print data
            print(file)
            print("Tau_rise = {:.3e} +/- {:.3e}".format(par[0], parErr[0]))
            print("Tau_fall = {:.3e} +/- {:.3e}".format(par[1], parErr[1]))
            print("A = {:.3e} +/- {:.3e}".format(par[2], parErr[2]))
            print("t_offset = {:.3e} +/- {:.3e}".format(par[3], parErr[3]))
            print("")
            pars.append(par)
            parErrs.append(parErr)
            fitFile.append(file)
        except:
            print(file)
            print("Optimal parameters not found")
            print("")
    
    return [times, amplitudes, pars, parErrs, fitFile]
    

def plotEvents(fitFile, times, amplitudes, pars, parErrs):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import get_cmap
    
    # plot
    fig = plt.figure(figsize=(6, 4))
    ax0 = plt.subplot(111)
    cmap = get_cmap('jet', lut=None)
    
    offset = 0.1
    
    for i,(file, t, A, par, parErr) in enumerate(zip(fitFile, times, amplitudes, pars, parErrs)):
        offset_i = i*offset
        
        color = cmap(i/len(fitFile))
        ax0.plot(t, A+offset_i, marker='.', linestyle='', markersize=1, color=color, alpha=0.03)
        ax0.plot(t, np.array(impulse_response(t, par[0], par[1], par[2], par[3]))+offset_i, 
                 color=color, linewidth=2, label=file[6:12]+r" $\tau_{{qp}}$ = {:.2f}".format(par[1]))
        
        
    handles, labels = ax0.get_legend_handles_labels()
    handles = np.flip(handles)
    labels = np.flip(labels)
    
    #ax0.set_xlim([-10, 250])
    #ax0.set_ylim([-0.2, 2.2])
    ax0.grid()
    ax0.set_ylabel('Relative amplitude')
    ax0.set_xlabel('Time [microsec]')
    ax0.legend(handles, labels, loc='best')
    plt.show()
    
def plotAmpTauCorrelation(IQ_timestreams, times, amplitudes, pars, parErrs):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import get_cmap
    from scipy.signal import medfilt
    
    max = []
    
    # plot
    fig = plt.figure(figsize=(9, 9))
    ax0 = plt.subplot(111)
    cmap = get_cmap('jet', lut=None)
    
    for i,(file, t, A, par, parErr) in enumerate(zip(IQ_timestreams, times, amplitudes, pars, parErrs)):
        #current_max = fmin(lambda x,a,b,c,d: -np.array(impulseResponse(x,a,b,c,d)), x0=0, args=(par[0], par[1], par[2], par[3]))
        current_max = np.max(medfilt(A, kernel_size=151))
        max.append(current_max)
        
        color = cmap(i/len(IQ_timestreams))
        
        ax0.plot(par[1], current_max, marker='o', linestyle='', markersize=5, color=color, alpha=1., label=file)
        
    handles, labels = ax0.get_legend_handles_labels()
    handles = np.flip(handles)
    labels = np.flip(labels)
    
    #ax0.set_xlim([-10, 250])
    #ax0.set_ylim([-0.2, 2.2])
    ax0.grid()
    ax0.set_xlabel(r'$\tau_{{qp}}$ [$\mu$s]')
    ax0.set_ylabel(r'Maximum [V]')
    #ax0.legend(handles, labels, loc='best')
    plt.show()