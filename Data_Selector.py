import os
import numpy as np
import matplotlib.pyplot as plt
import wfdb # read Physionet/picsdb file format
from picsdb import *
from wfdb import processing
from numba import jit

def getdata(infant,x_ecg,x_resp,fs_ecg,fs_resp):
    time_ecg = np.arange(x_ecg.shape[0])*dt_ecg
    time_resp = np.arange(x_resp.shape[0])*dt_resp
    freq_lo_ecg = 20
    if infant==1 or infant==5:
        freq_hi_ecg = 120
    else:
        freq_hi_ecg = 125
    freq_lo_resp = 0.5
    freq_hi_resp = 5 # 5
    x_ecg_filt  = bp_filter(x_ecg, fs_ecg, freq_lo_ecg, freq_hi_ecg)
    x_resp_filt = bp_filter(x_resp, fs_resp, freq_lo_resp, freq_hi_resp)
    if infant==1 or infant==5:
        x_ecg=processing.resample_sig(x_ecg_filt,fs=250,fs_target=500)[0]
    else:
        x_resp_filt=processing.resample_sig(x_resp_filt,fs=250,fs_target=500)[0]
        x_ecg=x_ecg_filt
    x_resp_flow=[]
    floewer=[]

    x_ecg_flow=[]
    ecg_flow=[]
    if infant==1 or infant==5:
        time_resp_temp=time_resp
        time_ecg_temp=time_resp
    else:
        time_resp_temp=time_ecg
        time_ecg_temp=time_ecg

    x_temp=[]
    temp=[]
    for i in range(0,len(x_resp_filt)-1,1):
        if i==0:
            floewer.append(x_resp_filt[i])
            ecg_flow.append(x_ecg[i])
            temp.append(time_resp_temp[i])
        else:
            floewer.append(x_resp_filt[i])
            temp.append(time_resp_temp[i])
            ecg_flow.append(x_ecg[i])
            if not (x_resp_filt[i]<x_resp_filt[i+1] and x_resp_filt[i]<=x_resp_filt[i-1]):
                pass
            else:
                x_resp_flow.append(floewer)
                x_ecg_flow.append(ecg_flow)
                x_temp.append(temp)
                floewer=[]
                ecg_flow=[]
                temp=[]
    for i in x_resp_flow:
        if(len(i)<200):
            x_resp_flow.remove(i)

    for j in x_ecg_flow:
        if(len(j)<200):
            x_ecg_flow.remove(j)

    for z in x_temp:
        if(len(z)<200):
            x_temp.remove(z)
    lomx=[]
    for z in x_ecg_flow:
        c=locmax(z)
        ecg_peaks = np.array([np.mean(z[i-1:i+1]) for i in c])
        x0, p_ecg_peaks = p_hat(ecg_peaks, n_bins=50, method='histogram')  #histogram
        c = np.array([j for i, j in enumerate(c) if ecg_peaks[i] > 0.12]) #ecg_peacks > 0.12 / 0.2
        lomx.append(c)
    timespan=2 #ms
    rate=[]
    for i in lomx:
        if not len(i)==0:
            span=2*(max(i)-min(i))/(len(i)-1)
            rate.append(60000/span)
        else:
            rate.append(0)
    ## do more start here