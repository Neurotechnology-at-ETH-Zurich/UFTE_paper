#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:13:38 2024

@author: baran
"""

import numpy as np
from matplotlib import pyplot as plt
from pandas import read_csv

siteMap = [27,28,12,29,11,20,7,21,10,1,15,8,26,18,31,23,9,2,14,6,25,17,30,22,0,19,13,24,4,16,3,5,43,32,56,45,60,35,52,40,41,33,46,38,57,49,62,54,42,34,47,39,58,50,63,55,59,44,61,36,53,37,51,48]

#%%
LL2_impedances = np.transpose(np.array(read_csv('mouseLL2_impedances_magnitudes.csv'))[:,1:])
LL2_phases = np.transpose(np.array(read_csv('mouseLL2_impedances_phases.csv'))[:,1:])

days_LL2 = [21,24,28,29,30,38,44,50,134,142,147,253,254,255,260,261,268,274,280,282,283,286,288,293,295,297,301,302,304,307]
LL2_nonbroken_impedances = {}
LL2_nonbroken_phases = {}

for day in range(len(days_LL2)):
    broken_channels = np.where(LL2_impedances[day] > 5e6)[0]
    non_broken_channels = np.delete(LL2_impedances[day], broken_channels)
    non_broken_channels_phases = np.delete(LL2_phases[day], broken_channels)
    non_broken_channel_idx = np.where(np.isin(LL2_impedances[day], non_broken_channels))[0]
    LL2_nonbroken_impedances[days_LL2[day]] = LL2_impedances[day,non_broken_channel_idx]
    LL2_nonbroken_phases[days_LL2[day]] = LL2_phases[day,non_broken_channel_idx]
    
LL2_impedance_stats = np.zeros((6,len(days_LL2)))
LL2_phase_stats = np.zeros((6,len(days_LL2)))

for i in range(len(days_LL2)):
    LL2_impedance_stats[0,i] = np.median(LL2_nonbroken_impedances[days_LL2[i]])
    LL2_impedance_stats[1,i] = np.quantile(LL2_nonbroken_impedances[days_LL2[i]],0.25)
    LL2_impedance_stats[2,i] = np.quantile(LL2_nonbroken_impedances[days_LL2[i]],0.75)
    LL2_impedance_stats[3,i] = np.mean(LL2_nonbroken_impedances[days_LL2[i]])
    LL2_impedance_stats[4,i] = np.std(LL2_nonbroken_impedances[days_LL2[i]])
    LL2_impedance_stats[5,i] = len(LL2_nonbroken_impedances[days_LL2[i]])
    
    LL2_phase_stats[0,i] = np.median(LL2_nonbroken_phases[days_LL2[i]])
    LL2_phase_stats[1,i] = np.quantile(LL2_nonbroken_phases[days_LL2[i]],0.25)
    LL2_phase_stats[2,i] = np.quantile(LL2_nonbroken_phases[days_LL2[i]],0.75)
    LL2_phase_stats[3,i] = np.mean(LL2_nonbroken_phases[days_LL2[i]])
    LL2_phase_stats[4,i] = np.std(LL2_nonbroken_phases[days_LL2[i]])
    LL2_phase_stats[5,i] = len(LL2_nonbroken_phases[days_LL2[i]])    
    

#%%
LL2_spike_sorting_results = ['mouseLL2_200818_134119.mat',
           'mouseLL2_200821_100626.mat',
           'mouseLL2_200825_132217.mat',
           'mouseLL2_200826_163558.mat',
           'mouseLL2_200827_154304.mat',
           'mouseLL2_200904_142528.mat',
           'mouseLL2_200910_170840.mat',
           'mouseLL2_200916_155143.mat',
           'mouseLL2_201209_144624.mat',
           'mouseLL2_201210_152147.mat',
           'mouseLL2_201217_163659.mat',
           'mouseLL2_201222_152605.mat',
           'mouseLL2_210407_135410.mat',
           'mouseLL2_210408_145158.mat',
           'mouseLL2_210409_130407.mat',
           'mouseLL2_210422_161822.mat',
           'mouseLL2_210517_101440.mat',
           'mouseLL2_210519_113603.mat',
           'mouseLL2_210520_123525.mat',
           'mouseLL2_210521_121322.mat',
           'mouseLL2_210525_110414.mat',
           'mouseLL2_210526_112217.mat',
           'mouseLL2_210528_100147.mat',
           'mouseLL2_210531_094913.mat']

LL2_mean_SNRs = np.zeros(len(LL2_unit_folders))
LL2_SNR_errs = np.zeros(len(LL2_unit_folders))
single_unit_counts = np.zeros(len(LL2_unit_folders))
LL2_SNRs = {}

"""
Run this part on the spike sorting result (res.mat) files to extract unit SNRs. 
Otherwise the SNRs of individual single units are already provided in the mouseLL2_single_unit_SNRs.csv

for i,session in enumerate(LL2_spike_sorting_results):
    file = 'spike_sorting_results/' + session
    data_dict = mat73.loadmat(file)
    single_unit = np.zeros(len(data_dict['clusterNotes']))
    for j in range(len(data_dict['clusterNotes'])):
        if data_dict['clusterNotes'][j] == ['single']:
            single_unit[j] = 1
    single_unit_counts[i] = len(np.where(single_unit)[0]) 
    non_zero_SNRs = data_dict['unitSNR'] * single_unit
    non_zero_SNRs = non_zero_SNRs[non_zero_SNRs != 0]
    LL2_mean_SNRs[i] = np.mean(non_zero_SNRs)
    LL2_SNR_errs[i] = np.std(non_zero_SNRs) / np.sqrt(len(non_zero_SNRs))
    LL2_SNRs[session] = non_zero_SNRs
"""

SNRs = read_csv('mouseLL2_single_unit_SNRs.csv')
for i,session in enumerate(LL2_spike_sorting_results):
    session_SNRs = np.nan_to_num(np.array(SNRs[str(days_LL2[i])]))
    non_zero_SNRs = session_SNRs[session_SNRs != 0]
    
    LL2_mean_SNRs[i] = np.mean(non_zero_SNRs)
    LL2_SNR_errs[i] = np.std(non_zero_SNRs) / np.sqrt(len(non_zero_SNRs))
    LL2_SNRs[session] = non_zero_SNRs

#%%
mouseLL2_lifetime_array = read_csv('mouseLL2_single_unit_lifetimes.csv')
mouseLL2_lifetime_array = np.array(mouseLL2_lifetime_array)[:,1:]

mouseLL2_lifetimes = np.zeros(mouseLL2_lifetime_array.shape[0])
for i in range(len(mouseLL2_lifetimes)):
    last_day = days_LL2[np.where(mouseLL2_lifetime_array[i])[0][-1]]
    first_day = days_LL2[np.where(mouseLL2_lifetime_array[i])[0][0]]
    mouseLL2_lifetimes[i] = last_day - first_day
    
#%%
y,x = np.histogram(mouseLL2_lifetimes,bins=10,range=(0,300))
y = y/np.sum(y)
plt.figure('Fig. S9d')
plt.bar(np.arange(1,10),np.cumsum(np.flip(y[1:])),width=0.9)
plt.show()
   
plt.figure('Fig. S9c')
plt.plot(days_LL2, LL2_mean_SNRs)
plt.fill_between(days_LL2, LL2_mean_SNRs-LL2_SNR_errs, LL2_mean_SNRs+LL2_SNR_errs,alpha=0.3)
plt.ylim([0,20])
plt.show()

plt.figure('Fig. S9a')
plt.plot(days_LL2[1:],LL2_impedance_stats[0,1:],color='b')
plt.fill_between(days_LL2[1:],LL2_impedance_stats[1,1:],LL2_impedance_stats[2,1:],alpha=0.1,color='b')
plt.yscale('log')
plt.ylim([1e5,1e7])
plt.ylabel('Impedances (Ohm)')
plt.xlabel('Days')
plt.show()    
    
plt.figure('Fig. S9b')
plt.plot(days_LL2[1:],LL2_phase_stats[0,1:],color='b')
plt.fill_between(days_LL2[1:],LL2_phase_stats[1,1:],LL2_phase_stats[2,1:],alpha=0.1,color='b')
plt.ylim([-90,0])
plt.ylabel('Phase')
plt.xlabel('Days')
plt.show()  
