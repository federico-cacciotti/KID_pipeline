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
    