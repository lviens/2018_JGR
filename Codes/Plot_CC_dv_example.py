'''
Written by L. Viens (initial release: 04/26/2018)
Cleaned by L. Viens (22/07/2020)

This code computes relative homogeneous velocity changes (dv/v) using the stretching technique from single-station cross-correlation functions
All the data for Viens et al. (2018, JGR) can be download at: https://doi.org/10.7910/DVN/TPDPHA
To run this code:
   - Download SC function data for at least one station 
   - Set the 'dir_ini' variable to the path where the .mat files are ((No need to change it if the data are in the Data folder)
   - Set the name of the station(s) in the 'virt' variable.
   - Set the filter bandwidth (cut and cut2 in Hz, with cut<cut2)
   - Set the number of hours for the stack (hr_to_stack)
   - Set the parameter 'Epsilon' (Epsilon = -dt/t = dv/v), which determines the range of values to look for dv/v changes
   - Set 't_ini' and 't_length' variables:
      - t_ini = 0  and t_length = 1  -> the dv/v is computed over the first second of the causal part of the SC function is selected
      - t_ini = 1  and t_length = 2  -> the dv/v is computed from the 1 s to 3 s after the zero lag time of the causal part of the SC function

The different steps are:
    - The code reads the .mat file (E.ABHM.mat) that contains the single station cross-correlation functions stacked over 6 h from February 20 to May 11, 2011 at the station.
    - Computes dv/v and correlation coefficient measurements for the vertical/North-south and vertical/East-west components and average them.
    - Plot the results and save the figure in a png file in the 'dir_out' folder

'''
from __future__ import division
import os
import numpy as np
import matplotlib.pyplot as plt
from Stretching_current import Stretching_current
from Additional_Functions import read_mat, butter_bandpass_filter, ddd2mmdd

#%% Main parameters
dir_ini = '../Data' # Directory folder 
dir_out = '../Figures' # Output directory to save the .png figure
virt = ['E.NS7M'] # Station name
cut  = 1     # Filter lower corner frequency (in Hz)
cut2 = 20    # Filter higher corner frequency (in Hz)
hr_to_stack = 6 # number of hours to stack
av_wind = 60 / 15 * hr_to_stack # 60/15 as the SC functions are computed every 15 min.

#%% Stretching parameters
Epsilon = .10  # Stretching between -Epsilon to +Epsilon (multiply by 100 to get the dv in %) STRETCHING LIMITS
t_ini = 0      # Time to start computing the dv/v (in second). Zero is at the zero lag time
t_length = 1   # Length of the signal over which the dv/v is computed (in second)

#%%
for sta in virt:    # Loop over the different stations
    print(sta + ': Data filter: ' + str(cut) + '-' + str(cut2) + ' Hz, dv/v computed between ' + str(t_ini) + ' s and ' + str(t_length) + ' s of the causal part, with Epsilon = ' + str(Epsilon*100)+' %')

    fid = dir_ini + os.sep + sta + '.mat'
    [UN_ref , UN , UE_ref , UE, delta, day ]= read_mat(fid) # Load the data
    
    # Get parameters for the stretching
    t_vec = np.arange(-(UE.shape[0]-1)/delta/2, UE.shape[0]/delta/2, 1/delta) # SCF time vector
    zero_lag_ind = round(len(UE[:,1])/2) # To compute the dv/v on the causal part of the SC functions.
    t_ini_d = t_ini*delta          # start dv/v computation at t_ini_d/delta seconds from the signal begining
    t_length_d = int(t_length*delta)  # dv/v computation over t_length_d/delta seconds after t_ini_d/delta
    dvE, ccE, errorE = np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)) # Set variables
    dvN, ccN, errorN = np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)), np.zeros((int(len(UE[1,:])/av_wind), 1)) # Set variables


    stck_ind = 1
    cp = 0
    datUE = np.zeros(len(UE[:,0]) )
    datUN = np.zeros(len(UN[:,0]) )
    for i in range(len(UE[1,:])):
        UE[:,i] = butter_bandpass_filter(data = UE[:,i] , lowcut = cut, highcut = cut2, fs = delta, order = 4)
        UN[:,i] = butter_bandpass_filter(data = UN[:,i] , lowcut = cut, highcut = cut2, fs = delta, order = 4)

        
        if stck_ind<av_wind: # Stack SC functions over av_wind hours
            datUE += UE[:,i]  
            datUN += UN[:,i]
            stck_ind+=1
        else: # At the last SC function to consider, perform stretching
            datUE += UE[:,i]  
            datUN += UN[:,i]
            datUE/=av_wind
            datUN/=av_wind
            stck_ind=1

            # Perform stretching
            [dvE[cp], ccE[cp], cdpE, dtE, errorE[cp], CE ] = Stretching_current(ref = UE_ref, cur = datUE, t = t_vec, dvmin = -Epsilon, dvmax = Epsilon, nbtrial = 50, window = np.arange(int(t_ini_d+zero_lag_ind), int(t_ini_d+t_length_d+zero_lag_ind)), fmin = cut, fmax = cut2, tmin = t_ini_d, tmax = t_length_d)

            [dvN[cp], ccN[cp], cdpN, dtN, errorN[cp], CN] = Stretching_current(ref = UN_ref, cur = datUN, t = t_vec, dvmin = -Epsilon, dvmax = Epsilon, nbtrial = 50,  window = np.arange(int(t_ini_d+zero_lag_ind),int(t_ini_d+t_length_d+zero_lag_ind)), fmin = cut, fmax = cut2, tmin = t_ini_d, tmax = t_length_d)
            
            cp += 1
            datUE = np.zeros(len(UE[:,0]) )
            datUN = np.zeros(len(UN[:,0]) )

    # Average the Z-E and Z-N dv/v and CC values using the Hobiger et al. (2014, GJI) method
    ponddv = (ccE**2 * dvE + ccN**2 * dvN)/(ccE**2 + ccN**2) # Combine Z-N and Z-E dv/v measurements
    pondcc = (ccE**3 + ccN**3)/(ccE**2 + ccN**2)             # Combine Z-N and Z-E correlation coefficients

        
    #%% Figure
    # Create day/month axis
    wind_per_day = 24/hr_to_stack
    day_unit = np.arange(day[0],day[0]+len(ponddv)/wind_per_day, 1/wind_per_day)
    x_ind = np.arange(int(day[0]) , int(day[-1]),10)
    mm, dd = ddd2mmdd(2011, x_ind)
    dat_plot = []
    for e in range(len(dd)):
        dat_plot.append(str(mm[e]) + '/' + str(dd[e]))
       
    # Plot
    fig = plt.figure(figsize=(8, 8))
    # Plot Correlation coefficient
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('Correlation coefficient')
    ax1.set_title(sta[2:] + '  (Filter: '+ str(cut) +'-'+str(cut2) + ' Hz)\n' +
                   'Stretching between '+ str(t_ini_d/delta) + ' and '+ str( t_length_d/delta) + ' s', fontsize = 15)              
    ax1.set_xlim(day[0], day[-1])
    ax1.plot(day_unit, pondcc, linewidth = 3)
    ax1.set_xticks(x_ind)
    ax1.axes.xaxis.set_ticklabels(dat_plot)
    ax1.set_ylim(0, 1)
    ax1.grid(which='major', linewidth=1, alpha=1)
    ax1.set_xlabel('Month / Day 2011')

    # Plot dv/v
    ax2 = fig.add_subplot(212)
    ax2.scatter(x = day_unit, y = ponddv, s = 45, c = np.abs(pondcc), cmap = 'hot_r', edgecolors = 'k')
    ax2.set_ylabel('dv/v (%)')
    ax2.set_ylim(1.1*np.nanmin((ponddv)), 1.5*np.nanmax((ponddv)))
    ax2.grid(which='major', linewidth = 1, alpha=0.5)
    ax2.set_xticks(x_ind)
    ax2.axes.xaxis.set_ticklabels(dat_plot)
    ax2.set_xlim(day[0], day[-1])
    ax2.set_xlabel('Month / Day 2011')
    
    plt.tight_layout()
    plt.show()
    
    # Save figure
    outfname = dir_out + '/' + 'Fig_dv_' + sta + '.png'
    fig.savefig(outfname, format='png', dpi=100)
    plt.close(fig)
