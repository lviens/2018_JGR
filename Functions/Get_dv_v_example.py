"""
L. Viens 04/26/2018
This code computes relative homogeneous velocity changes (dv/v) using the stretching technique from single-station cross-correlation functions
All the data for Viens et al. (2018, submitted to JGR) can be download at: https://doi.org/10.7910/DVN/TPDPHA
To run this code:
   - Download SC function data for at least one station
   - Set the "dir_ini" variable to the path where the .mat files are.
   - Set the name of the station(s) in the "virt" variable.
   - Set the filter bandwidth (cut and cut2 in Hz, with cut<cut2)
   - Set the number of hours for the stack (hr_to_stack)
   - Set the parameter "Epsilon" (Epsilon = -dt/t = dv/v), which determines the range of values to look for dv/v changes
   - Set "t_ini" and "t_length" variables:
      - t_ini = 0  and t_length = 1  -> the dv/v is computed over the first second of the causal part of the SC function is selected
      - t_ini = 1  and t_length = 2  -> the dv/v is computed from the 1 s to 3 s after the zero lag time of the causal part of the SC function
The different steps are:
    - Reads the mat file (E.ABHM.mat) that contains the single station cross-correlation functions stacked over 6 h from February 20 to May 11, 2011 at the E.ABHM station.
    - Computes dv/v and correlation coefficient measurements for the vertical/North-south and vertical/East-west components and average them.
    - Plot the results and save the figure in a pdf file in the "dir_out" folder

"""
from __future__ import division
import scipy.io 
import numpy as np
import scipy
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
import os
from scipy.signal import butter, filtfilt
from Stretching_current import *


def main():
    dir_ini = './' # Directory folder 
    dir_out = './' # Output directory to save pdf file
    virt = ['E.ABHM'] # Station name
    year = 2011  # Year
    cut  = 1     # Pass band low corner frequency of the filter in Hz
    cut2 = 20    # Pass band high corner frequency of the filter in Hz
    hr_to_stack = 6 # number of hours to stack
    av_wind = 60/15*hr_to_stack # 60/15 as the SC functions are computed every 15 min.

    # Stretching parameters
    Epsilon = .10  # Stretching between -Epsilon to +Epsilon (multiply by 100 to get the dv in %) STRETCHING LIMITS
    t_ini = 0      # Time to start computing the dv/v (in second). Zero is at the zero lag time
    t_length = 1   # Length of the signal over which the dv/v is computed (in second)

    for sta in virt:    # Loop over the different stations
        print(sta + ": Data filter: " + str(cut) + "-" + str(cut2) + " Hz, dv/v computed between " + str(t_ini) + " s and " + str(t_length) + " s of the causal part, with Epsilon = " + str(Epsilon*100)+" %")

        fid = dir_ini + os.sep + sta + '.mat'
        [UN_ref , UN , UE_ref , UE, delta, day, stack_hr ]= read_mat(fid) # Reads the file using the function read_hdf5 (see below)
        # Get parameters for the stretching
        t_vec = np.arange(-(UE.shape[0]-1)/delta/2, UE.shape[0]/delta/2, 1/delta) # Time axis
        zero_lag_ind = round(len(UE[:,1])/2) # To compute the dv/v on the causal part of the SC functions.
        t_ini_d = t_ini*delta          # start dv/v computation at t_ini_d/delta seconds from the signal begining
        t_length_d = int(t_length*delta)  # dv/v computation over  t_length_d/delta seconds after t_ini_d/delta
        dvE, ccE, errorE = np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)) # Set variables
        dvN, ccN, errorN = np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)) # Set variables

        datUE = []
        datUN = []
        a=1
        cp =0
        for i in range(len(UE[1,:])):
            UE[:,i] = butter_bandpass_filter(data=UE[:,i] , lowcut=cut, highcut=cut2, fs=delta, order=4)
            # Stacking over av_wind hours of data
            if a<av_wind:
                if a==1:
                    datUE = UE[:,i] 
                    datUN = UN[:,i]
                else:
                    datUE = datUE + UE[:,i]  
                    datUN = datUN + UN[:,i]
                a+=1;
            else:
                a=1;

                # Perform stretching
                [dvE[cp], ccE[cp], cdpE, dtE, errorE[cp], CE ] = Stretching_current(ref = UE_ref, cur = datUE, t = t_vec, dvmin = -Epsilon, dvmax = Epsilon, nbtrial = 50,
                                                                  window = np.arange(int(t_ini_d+zero_lag_ind),int(t_ini_d+t_length_d+zero_lag_ind)), 
                                                                  fmin = cut, fmax = cut2, tmin = t_ini_d, tmax = t_length_d)

                [dvN[cp], ccN[cp], cdpN, dtN, errorN[cp], CN] = Stretching_current(ref = UN_ref, cur = datUN, t = t_vec, dvmin = -Epsilon, dvmax = Epsilon, nbtrial = 50,
                                                              window = np.arange(int(t_ini_d+zero_lag_ind),int(t_ini_d+t_length_d+zero_lag_ind)), 
                                                              fmin = cut, fmax = cut2, tmin = t_ini_d, tmax = t_length_d)
                cp+=1
                datUE = []
                datUN = []

        ponddv = (ccE**2 * dvE + ccN**2 * dvN)/(ccE**2 + ccN**2) # Combine NZ and EZ dv/v measurements, based on Hobiger et al. (2014, GJI)
        pondcc = (ccE**3 + ccN**3)/(ccE**2 + ccN**2)             # Combine NZ and EZ correlation coefficients, based on Hobiger et al. (2014, GJI)

        # Plot with the function plot_dv_v
        [day_unit, dat_plot] = plot_dv_v(pondcc, ponddv, day, hr_to_stack , sta, cut, cut2, t_ini_d, delta, t_length_d, dir_out)




def plot_dv_v(pondcc, ponddv, day, stack_hr, virt, cut, cut2, tbis, delta, length_t, directory):
    """
    Function to plot the correlation coefficients and dv/v measurements through time and save the result in the directory folder as: Fig_dv_E.OMNM.pdf'
    """
    wind_per_day = 24/stack_hr
    day_unit = np.arange(day[0],day[0]+len(ponddv)/wind_per_day, 1/wind_per_day)
    x_ind = np.arange(int(day[0]) , int(day[-1]),10)
    mm, dd = ddd2mmdd(2011, x_ind)
    dat_plot = []
    for e in range(len(dd)):
        dat_plot.append(str(mm[e]) + '/' + str(dd[e]))

    # Figure
    fig = plt.figure(figsize=(6, 8))
    fig.subplots_adjust(hspace=0.6)

    # PLOT CC
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel("Correlation coefficient")
    ax1.set_title(virt + '  (Filter: '+("%1.0f" % cut)+'-'+("%1.0f" % cut2)+' Hz)\n'+
                  'Stretching computed between '+("%2.1f" % (tbis/delta))+' and '+
                  ("%2.1f" % (tbis/delta+length_t/delta))+' s\n')
    ax1.set_xlim(day[0], day[-1])
    ax1.set_xticks(x_ind)
    ax1.axes.xaxis.set_ticklabels(dat_plot)
    ax1.set_ylim(0, 1)
    ax1.grid(which='major', linewidth=1, alpha=0.5)
    ax1.plot(day_unit, pondcc)

    ax1.set_xlabel("Month / Day 2011")

    # Plot dv/v
    ax2 = fig.add_subplot(212)
    ax2.set_ylabel("dv/v (%)")
    ax2.set_xlim(day[0], day[-1])
    ax2.set_ylim(1.1*np.min((ponddv)), 1.1*np.max((ponddv)))
    ax2.grid(which='major', linewidth=1, alpha=0.5)
    ax2.scatter(x=day_unit, y=ponddv, s=15, c=np.abs(pondcc), cmap='hot_r')

    ax2.set_xticks(x_ind)
    ax2.axes.xaxis.set_ticklabels(dat_plot)
    ax2.set_xlabel("Month / Day 2011")

    plt.show()
    # Save figure
    outfname = directory + '/' + 'Fig_dv_' + virt + '.pdf'
    fig.savefig(outfname, format='pdf', dpi=400)
    plt.close(fig)

    return day_unit, dat_plot



def read_mat(fid):
    """
    Function to read the .mat file data.
    INPUT: 
        - fid = file folder and name 
    OUTPUT:
        - UN_ref = Reference waveform for the vertical and North-South component of the single-station cross correlation
        - UN = Vertical and North-South matrix of single-station cross correlation
        - UE_ref = Reference waveform for the vertical and East-West component of the single-station cross correlation
        - UE = Vertical and East-West matrix of single-station cross correlation
        - delta = Sampling frequency in Hz
        - day = Day of the year 2011 for which the single-station cross-correlations were computed
        - stack_hr = The initial 15 min SC functions are stacked over "stack_hr" hours (in this example, stack_hr is 6 hours -> 4 waveforms per day)
    """
    mat = scipy.io.loadmat(fid)
    data = mat['dat']
    UN = data[0,0]['UN']
    UE = data[0,0]['UE']
    UN_ref = data[0,0]['ref_UN']
    UE_ref = data[0,0]['ref_UE'] 
    delta = np.squeeze(data[0,0]['delta'] )
    day =np.squeeze(data[0,0]['day'] )
    stack_hr = 60/15 #data[0,0]['t_step'] 
    UN = UN.transpose()
    UE = UE.transpose()
    UE_ref= np.squeeze(UE_ref)
    UN_ref= np.squeeze(UN_ref)

    return UN_ref , UN , UE_ref , UE, delta, day, stack_hr


def butter_bandpass(lowcut, highcut, fs, order=4):
    """
    Function to define the filter
    INPUT: 
        - lowcut = lower frequency of the filter (in Hz)
        - highcut = lower frequency of the filter (in Hz)
        - fs = Sampling rate in Hz
        - order = order of the filter
    OUTPUT:
        - [a,b] =  Numerator (b) and denominator (a) polynomials of the IIR filter.
    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    """
    Function to apply a Butterworth bandpass-filter, two-pass to the data
    INPUT: 
        - data = data time series 
        - lowcut = lower frequency of the filter (in Hz)
        - highcut = lower frequency of the filter (in Hz)
        - fs = Sampling rate in Hz
        - order = order of the filter
    OUTPUT:
        - y = filtered time series
    """
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def ddd2mmdd(year, ddd):
    """
    Function to get the day and month of the year from the day of the year.
    INPUT: 
        - year = year 
        - ddd = day of the year ( can be a vector)
    OUTPUT:
        - mm = month of the year for the input days (ddd)
        - dd = day of the month for the input days
    """
    v = [datetime(year, 1, 1) + timedelta(days=int(dy-1)) for dy in ddd]
    mm = [v[ind].month for ind in range(len(v))]
    dd = [v[ind].day for ind in range(len(v))]
    return mm, dd


if __name__=="__main__":
   main()
