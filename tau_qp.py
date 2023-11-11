import numpy as np
from pathlib import Path

#paths_cosmic_rays = Path("/Users/federicocacciotti/Documents/PhD/dati/cosmic_rays")


'''
        Impulse analysis for quasiparticle recombination time measurements
'''
    

def impulse_response(t, params):
    tau_rise = params['tau_rise']
    tau_fall = params['tau_fall']
    A = params['A']
    t_offset = params['t_offset']
    v_offset = params['v_offset']
    
    y = v_offset + A*(np.exp(-(t-t_offset)/tau_fall) - np.exp(-(t-t_offset)/tau_rise))
    y = [v_offset if y_<=v_offset else y_ for y_ in y]
    return np.array(y)

def double_exponential(t, tau_rise, tau_fall, tau_therm, A1, A2, t_offset):
    y = A1*(np.exp(-(t-t_offset)/tau_fall) - np.exp(-(t-t_offset)/tau_rise)) + A2*np.exp(-(t-t_offset)/tau_therm)
    y = [0.0 if y_<0 else y_ for y_ in y]
    return y

def exponential_decay(t, A, tau, t_offset):
    return A*np.exp(-(t-t_offset)/tau)


class Event():
    
    def __init__(self, filename, CH1=None, CH2=None, CH3=None, CH4=None, I_DC=None, Q_DC=None, label='Event', n_sigma=2.0):
        
        self.label = label
        self.time = {'label': 'TIME', 'data': []}
        self.channels = {}
        self.filename = Path(filename)
        self.n_sigma = n_sigma
        self.I_DC = I_DC
        self.Q_DC = Q_DC
        print("Reading "+self.filename.stem+"...")
        
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
        data = np.loadtxt(self.filename, dtype=float, delimiter=',', skiprows=21, unpack=True, usecols=range(n_channels+1))
        self.vertical_units = str(np.loadtxt(self.filename, dtype=str, delimiter=',', skiprows=12, usecols=1, max_rows=1))
        self.horizontal_units = str(np.loadtxt(self.filename, dtype=str, delimiter=',', skiprows=5, usecols=1, max_rows=1))
        
        
        self.time['data'] = data[0]
        self.time['data'] -= self.time['data'][0]   # remove time offset
        
        j = 1
        for CH in ['CH1', 'CH2', 'CH3', 'CH4']:
            try:
                self.channels[CH]['data'] = data[j]
                j += 1
            except:
                pass
        
        if self.I_DC == None or self.Q_DC==None:
            # add a constant >> Q and I
            C = 1000.0
            if np.abs(self.channels[self.I_channel]['data'].min()) > np.abs(self.channels[self.I_channel]['data'].max()):
                self.I = self.channels[self.I_channel]['data'] - C
            else:
                self.I = self.channels[self.I_channel]['data'] + C
            if np.abs(self.channels[self.Q_channel]['data'].min()) > np.abs(self.channels[self.Q_channel]['data'].max()):
                self.Q = self.channels[self.Q_channel]['data'] - C
            else:
                self.Q = self.channels[self.Q_channel]['data'] + C
            self.A = np.sqrt(self.I*self.I + self.Q*self.Q) - np.sqrt(2.0*C**2.0)
        else:
            self.I = self.channels[self.I_channel]['data'] + self.I_DC
            self.Q = self.channels[self.Q_channel]['data'] + self.Q_DC
            self.A = np.sqrt(self.I*self.I + self.Q*self.Q)
        self.phase = np.arctan2(self.Q, self.I)
        self.phase = np.unwrap(self.phase)
            
        # check if fit parameters already exists
        self.fit_result = None
        self.par = None
        try:
            self.fit_result = np.load(self.filename.parent / (self.filename.stem+'.npy'), allow_pickle=True).tolist()
            self.par = self.fit_result.params
            print("Found fit parameters!")
        except FileNotFoundError:
            print("Fit parameters do not exist.")
            
        
    def plot(self):
        import matplotlib.pyplot as plt
        
        colors = ['gray', 'red', 'blue', 'green']
        
        fig, axis = plt.subplots(1, 3, sharex=True, figsize=(18,7))
        
        plt.subplots_adjust(hspace=0.2)
        
        for CHi, color in zip(self.channels, colors):
            axis[0].plot(self.time['data'], self.channels[str(CHi)]['data'], linestyle='solid', label=self.channels[str(CHi)]['label'], color=color)
    
        fig.suptitle(self.label)
        
        # plot amplitude and phase
        axis[1].plot(self.time['data'], self.A, color='black', label='Amplitude', linestyle='solid')
        axis[2].plot(self.time['data'], self.phase, color='black', label='Phase', linestyle='solid')
        
        for ax in axis:
            ax.legend(loc='best')
            ax.grid(color='gray', alpha=0.4)
            ax.set_xlabel('Time ['+self.horizontal_units+']')
        axis[0].set_ylabel('Voltage ['+self.vertical_units+']')
        axis[1].set_ylabel('Voltage ['+self.vertical_units+']')
        axis[2].set_ylabel('Phase [rad]')
        
        plt.show()
        
    def compute_errorbars(self):
        bin_heights, bin_position = np.histogram(self.A, 200)
        bin_max_height = np.max(bin_heights)
        bin_max_argument = np.argmax(bin_heights)
        bin_max_position = bin_position[bin_max_argument]
        
        from lmfit import Minimizer, Parameters
        def gaussian(x, params):
            return params['A']*np.exp(- 0.5 * (x-params['mu'])**2.0/params['sigma']**2.0)
        def fcn2min(params, x, data):
            return data-gaussian(x, params)
        
        params = Parameters()
        params.add('A', value=bin_max_height, min=bin_max_height*0.7, max=bin_max_height*1.3)
        params.add('mu', value=bin_max_position, min=bin_max_position*0.9, max=bin_max_position*1.1)
        params.add('sigma', value=1.0e-3, min=0.0, max=1.0)
        
        try:
            keep = np.logical_and(bin_position[0:-1] <= self.n_sigma*(bin_max_position-bin_position[0]), bin_heights != 0.0)
            result = Minimizer(fcn2min, params, fcn_args=(bin_position[0:-1][keep], bin_heights[keep])).minimize()
        except TypeError:
            keep = np.logical_and(bin_position[0:-1] <= 4.0*(bin_max_position-bin_position[0]), bin_heights != 0.0)
            result = Minimizer(fcn2min, params, fcn_args=(bin_position[0:-1][keep], bin_heights[keep])).minimize()

        return bin_heights, bin_position, result
        
    
    def fit(self, force_fit=False, method='least_squares', time_mask=None):
        # check if fit parameters already exist
        if not np.any(self.par == None) and force_fit == False:
            print("Fit parameters already exist. Try to force the fitting routine by passing 'force_fit=True' to the fit function.")
            return 0
        
        from lmfit import Minimizer, Parameters, report_fit
        
        # fit
        # errorbars are the std dev of the lower part of the stream
        bin_heights, bin_position, result = self.compute_errorbars()
        
        errors = np.ones(len(self.time['data']))*result.params['sigma'].value
        
        def fcn2min(params, x, data, errs=None):
            model = impulse_response(x, params)
            
            if errs is None:
                return model-data
            
            return (model-data)/errs
        
        # parameters initial guess
        mask = self.A-self.A.min() >= 0.6*(self.A.max()-self.A.min())
        t_offset = self.time['data'][mask][0]
        t_fall_offset = self.time['data'][mask][-1]
        
        mask = self.A-self.A.min() >= 0.90*(self.A.max()-self.A.min())
        time_at_max = 0.5*(self.time['data'][mask][0] + self.time['data'][mask][-1])
        
        params = Parameters()
        params.add('tau_rise', value=time_at_max-t_offset, min=0.0, max=1.0e-2)
        params.add('tau_fall', value=t_fall_offset-time_at_max, min=0.0, max=1.0e-2)
        params.add('A', value=self.A.max(), min=0.0, max=np.inf)
        params.add('t_offset', value=t_offset, min=-1.0, max=1.0)
        params.add('v_offset', value=result.params['mu'].value, min=-1.0, max=1.0)
        
        if np.any(time_mask) == None:
            time_mask = self.time['data'] < np.inf
        
        try:
            result = Minimizer(fcn2min, params, fcn_args=(self.time['data'][time_mask], self.A[time_mask], errors[time_mask])).minimize(method=method)
            self.fit_result = result
            self.par = result.params
            np.save(self.filename.parent / (self.filename.stem+'.npy'), result)
            
            
            # print data
            report_fit(result, show_correl=False)
        except ValueError:
            print("Optimal parameters not found")
    
    
    def plotFit(self, save_figure=True, show_figure=False, title=None):
        import matplotlib.pyplot as plt
        
        fig = plt.figure(figsize=(7, 7))
        ax0 = plt.subplot(111)
        
        ax0.plot(self.time['data'], self.A, linestyle='solid', label='Data', color='gray')
        
        FIR = impulse_response(self.time['data'], self.par)
        ax0.plot(self.time['data']*1e6, FIR, linestyle='solid', linewidth=4, label='Fit', color='red', alpha=0.5)
        
        ax0.legend(loc='best')
        ax0.grid(color='gray', alpha=0.4)
        ax0.set_xlabel('Time [s]')
        ax0.set_ylabel('Amplitude [V]')
        
        if title != None:
            plt.title(title)
        else:
            plt.title(self.label)
        
        if save_figure:
            fig.savefig(self.filename.parents[1]/(self.label+'.png'), dpi=250)
        if show_figure:
            plt.show()
        
    def fancyPlotFit(self, save_figure=True, show_figure=False, title=None):
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        
        fig = plt.figure(figsize=(7, 7))
        plt.clf()
        fig.clear()
        plt.subplots_adjust(hspace = 0.0, wspace=0.0)
        gs0 = gridspec.GridSpec(nrows=5, ncols=5, figure=fig)
        
        if title != None:
            fig.suptitle(title)
        else:
            fig.suptitle(self.label)
        
        # FIR PLOT
        ax0 = fig.add_subplot(gs0[0:4, 0:4])
        ax0.tick_params(axis='both', which='both', direction='in', bottom=True, left=True, top=True, right=True)
        ax0.set_xticklabels([])
        ax_res = fig.add_subplot(gs0[4, 0:4])
        ax_res.tick_params(axis='both', which='both', direction='in', bottom=True, left=True, top=True, right=True)
        ax_hist = fig.add_subplot(gs0[0:4, 4])
        ax_hist.tick_params(axis='both', which='both', direction='in', bottom=True, left=True, top=True, right=True)
        ax_hist.set_yticklabels([])
        
        FIR = impulse_response(self.time['data'], self.par)
        ax0.plot(self.time['data'], self.A, linestyle='solid', label='Data', color='gray', linewidth=1)
        ax0.plot(self.time['data'], FIR, linestyle='solid', linewidth=4, label='Fit', color='red', alpha=0.5)
        ax0.legend(loc='best')
        ax0.grid(color='gray', alpha=0.4)
        ax0.set_xlabel('Time [s]')
        ax0.set_ylabel('Amplitude [V]')
        
        # RESIDUALS PLOT
        res = self.A-FIR
        ax_res.plot(self.time['data'], res, linestyle='solid', linewidth=1, label='Residuals', color='black')
        ax_res.legend(loc='best')
        ax_res.grid(color='gray', alpha=0.4)
        ax_res.set_xlabel('Time [s]')
        
        # HISTOGRAM PLOT
        bin_heights, bin_position, result = self.compute_errorbars()
        plt.stairs(bin_heights, bin_position, orientation='horizontal', fill=False, color='black')
        
        def gaussian(x, params):
            return params['A']*np.exp(- 0.5 * (x-params['mu'])**2.0/params['sigma']**2.0)
            
        x_data = np.linspace(result.params['mu']-5.0*result.params['sigma'], result.params['mu']+5.0*result.params['sigma'], num=500)
        ax_hist.plot(gaussian(x_data, result.params), x_data, color='blue')
        ax_hist.grid(color='gray', alpha=0.4)
        ax_hist.set_xlabel('Counts')
        ax_hist.set_ylim(ax0.get_ylim())
        ax_hist.xaxis.tick_top()
        ax_hist.xaxis.set_label_position('top') 
        
        if save_figure:
            fig.savefig(self.filename.parents[1]/(self.label+'.png'), dpi=250)
        if show_figure:
            plt.show()
    
    def printParameters(self):
        try:
            string = "{:s},{:.4e},{:.4e},{:.4e},{:.4e},{:.4e},{:.4e},{:.4e},{:.4e},{:.4e},{:.4e}".format(self.label, self.par['tau_rise'].value*1e6, self.par['tau_rise'].stderr*1e6, self.par['tau_fall'].value*1e6, self.par['tau_fall'].stderr*1e6, self.par['A'].value, self.par['A'].stderr, self.par['t_offset'].value*1e6, self.par['t_offset'].stderr*1e6, self.par['v_offset'].value, self.par['v_offset'].stderr)

            print(string)
        except:
            print("--")
            pass

        
        
        
def add_batch(tek_numbers, path_to_tek_files, CH1=None, CH2=None, CH3=None, CH4=None, label_prefix='Event_', tek_prefix='tek', tek_extension='.csv', digits=4):
    batch = []
    for i,tek_number in enumerate(tek_numbers):
        batch.append(Event(path_to_tek_files+'/'+tek_prefix+'{:04d}'.format(tek_number)+tek_extension, CH1=CH1, CH2=CH2, CH3=CH3, CH4=CH4, label=label_prefix+'{:d}'.format(i)))
    return batch
        

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
        t, chI, chQ = np.loadtxt(file, float, delimiter=',', skiprows=21, unpack=True)
        
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
        t, chI, chQ = np.loadtxt(file, float, delimiter=',', skiprows=21, unpack=True)
        
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
        t, chI, chQ = np.loadtxt(file, float, delimiter=',', skiprows=21, unpack=True)
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
    

def plotEvents(events, y_offset=0.3, colormap='coolwarm', plot_fit=False):
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import get_cmap
    
    # plot
    fig = plt.figure(figsize=(7, 7))
    ax = fig.gca()
    cmap = get_cmap(colormap, lut=None)
    
    for i,e in enumerate(events):
        offset_i = i*y_offset
        
        color = cmap(i/len(events))
        ax.plot(e.time['data']*1e6, e.A+offset_i, linestyle='solid', color=color, label=e.label, alpha=1.0)
        if plot_fit:
            FIR = impulse_response(e.time['data'], e.par['tau_rise'].value, e.par['tau_fall'].value, e.par['A'].value, e.par['t_offset'].value, e.par['v_offset'].value)
            ax.plot(e.time['data']*1e6, FIR+offset_i, color=color, linewidth=2)
        
    ax.grid()
    ax.set_ylabel('Amplitude + offset')
    ax.set_xlabel('Time [microsec]')
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
