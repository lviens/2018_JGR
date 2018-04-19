def stretching_with_output(ref, cur, t, dvmin, dvmax, nbtrial, window, fmin, fmax, tmin, tmax):
"""
    Stretching function: 
    This function compares the Reference waveform to stretched/compressed current waveforms to get the relative seismic velocity variation (and associated error).
    It also computes the correlation coefficient between the Reference waveform and the current waveform.

    INPUTS:
        - ref = Reference waveform (size N)
        - cur = Current waveform (size N)
        - t = time vector, common to both ref and cur (size N)
        - dvmin = minimum bound for the velocity variation (example: dvmin=-0.03 for -3% of relative velocity change)
        - dvmax = maximum bound for the velocity variation (example: dvmin=0.03 for 3% of relative velocity change)
        - nbtrial = number of stretching coefficient between dvmin and dvmax (no need to be > 100)
        - window = vector of the indices of the cur and ref windows on wich you want to do the measurements
        For error computation:
            - fmin = minimum frequency of the data
            - fmax = maximum frequency of the data
            - tmin = minimum time window where the dv/v is computed 
            - tmax = maximum time window where the dv/v is computed 

    OUTPUTS:
        - dv = Relative velocity change dv/v (in %)
        - cc = correlation coefficient between the reference waveform and the best stretched/compressed current waveform
        - cdp = correlation coefficient between the reference waveform and the initial current waveform
        - Eps = Vector of Epsilon values (Epsilon =-dt/t = dv/v)
        - error = Errors in the dv/v measurements based on Weaver, R., C. Hadziioannou, E. Larose, and M. Camnpillo (2011), On the precision of noise-correlation interferometry, Geophys. J. Int., 185(3), 1384?1392
        - C = vector of the correlation coefficient between the reference waveform and every stretched/compressed current waveforms

    The code first finds the best correlation coefficient between the Reference waveform and the stretched/compressed current waveform among the "nbtrial" values. 
    A refined analysis is then performed around this value to obtain a more precise dv/v measurement .

    """     
    Eps = np.asmatrix(np.linspace(dvmin, dvmax, nbtrial))
    
    L = 1 + Eps
    tt = np.matrix.transpose(np.asmatrix(t))

    tau = tt.dot(L)  # stretched/compressed time axis
    print(tau.shape)
    # print(tau)
    C = np.zeros((1, np.shape(Eps)[1]))
    print(np.shape(Eps)[1])
    # blif_ref = []
    blif_ref = np.zeros((np.shape(Eps)[1], 200))
    blif_curr_stretc = np.zeros((np.shape(Eps)[1],200))
    # Set of stretched/compressed ref GF
    for j in np.arange(np.shape(Eps)[1]):
        s = np.interp(x=np.ravel(tt), xp=np.ravel(tau[:, j]), fp=cur) # np.ravel converts matrix back to 1D array
        # print(s)
        # print(ref.shape)
        bli1 = ref[window]
        bli2 = s[window]
        print(bli1.shape)
        blif_ref[j,:] = bli1
        blif_curr_stretc[j,:] = bli2
        C[0, j] = np.corrcoef(bli1, bli2)[0, 1]

    cdp = np.corrcoef(cur[window], ref[window])[0, 1]

    imax = np.nanargmax(C)
    # print(imax)
    if imax == np.shape(Eps)[1]:
        imax = imax - 2
    if imax <= 2:
        imax = imax + 2

    dtfiner = np.linspace(Eps[0, imax-2], Eps[0, imax+1], 500)
    # np.interpolation around the max to find the real max
    # 29/11/17 Chris altered this line so that spline np.interpolation takes two points either side of imax, instead of one
    func = scipy.interpolate.interp1d(np.ravel(Eps[0, np.arange(imax-3, imax+2)]), np.ravel(C[0,np.arange(imax-3, imax+2)]), kind='cubic')
    CCfiner = func(dtfiner)
    # CCfiner = scipy.np.interpolate.spline(xk=np.ravel(dt[0, np.arange(imax-3, imax+2)]),
    #                                       yk=np.ravel(C[0,np.arange(imax-3, imax+2)]), xnew=dtfiner, order=3)
    cc = np.max(CCfiner)
    # Multiply by 100 to convert to percentage (Epsilon = -dt/t - dv/v)
    dv = -100. * dtfiner[np.argmax(CCfiner)]

    # Error computation
    T = 1 / (fmax - fmin)
    X = cc
    wc = np.pi * (fmin + fmax)
    t1 = np.min([tmin, tmax])
    t2 = np.max([tmin, tmax])
    # Based on Weaver, R., C. Hadziioannou, E. Larose, and M. Camnpillo (2011), On the precision of noise-correlation interferometry, Geophys. J. Int., 185(3), 1384?1392
    error = 100*(np.sqrt(1-X**2)/(2*X)*np.sqrt((6* np.sqrt(np.pi/2)*T)/(wc**2*(t2**3-t1**3))))

    return dv, cc, cdp, dt, error, C, blif_curr_stretc , blif_ref
