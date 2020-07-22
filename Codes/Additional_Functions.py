#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 12:42:20 2020

@author: loic
"""
from __future__ import division
import scipy.io 
import numpy as np
import scipy
from datetime import datetime,timedelta
from scipy.signal import butter, filtfilt


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
    """
    mat = scipy.io.loadmat(fid)
    data = mat['dat']
    UN = data[0,0]['UN']
    UE = data[0,0]['UE']
    UN_ref = data[0,0]['ref_UN']
    UE_ref = data[0,0]['ref_UE'] 
    delta = np.squeeze(data[0,0]['delta'] )
    day = np.squeeze(data[0,0]['day'] )
    UN = UN.transpose()
    UE = UE.transpose()
    UE_ref = np.squeeze(UE_ref)
    UN_ref = np.squeeze(UN_ref)

    return UN_ref , UN , UE_ref , UE, delta, day


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
    b, a = butter_bandpass(lowcut, highcut, fs, order = order)
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
    v = [datetime(year, 1, 1) + timedelta(days = int(dy-1)) for dy in ddd]
    mm = [v[ind].month for ind in range(len(v))]
    dd = [v[ind].day for ind in range(len(v))]
    return mm, dd

