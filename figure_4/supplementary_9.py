#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 14:13:38 2024

@author: baran
"""

from load_intan_rhd_format import read_data
import numpy as np
from matplotlib import pyplot as plt
import mat73

siteMap = [27,28,12,29,11,20,7,21,10,1,15,8,26,18,31,23,9,2,14,6,25,17,30,22,0,19,13,24,4,16,3,5,43,32,56,45,60,35,52,40,41,33,46,38,57,49,62,54,42,34,47,39,58,50,63,55,59,44,61,36,53,37,51,48]

def read_impedances_from_metadata(folder, siteMap=siteMap):
    """This function reads the impedance measurements stored in the Intan meta data and exports it into a .csv file, with channels ordered according to their spatial organization. This should be used in case the impedance measurements are not stored in a .csv file.

    Inputs:
        -folder (string): Directory to the folder containing the impedance measurements
        -siteMap (1xN list where N is the number of channels in the electrode array):
        the list that contains the information on how the Intan amplifier channels (e.g. 0-63 for port A-000 to A-063; 64-127 for port B-000 to B-063 etc.) are organized spatially (i.e. 0-63 for sites on first shank with 64 channels dorsal to ventral, 64-128 for sites on the second shank with 64 channels dorsal to ventral etc.)
    """

    meta_data = read_data(folder+'/info.rhd')
    amplifier_channels = meta_data['amplifier_channels']

    channel_ids = []
    impedances = []
    phases = []
    native_orders = []

    for i in range(len(amplifier_channels)):
        channel_ids.append(amplifier_channels[i]['native_channel_name'])
        impedances.append(amplifier_channels[i]['electrode_impedance_magnitude'])
        phases.append(amplifier_channels[i]['electrode_impedance_phase'])
        native_orders.append(amplifier_channels[i]['native_order'])

    impedances_full = []
    phases_full = []
    count = 0
    for site in siteMap:
        if site in native_orders:
            impedances_full.append(impedances[count])
            phases_full.append(phases[count])
            count = count + 1
        else:
            impedances_full.append(10000000.0)
            phases_full.append(0)
    return impedances_full, phases_full

#%%
LL2_impedance_folders = ['/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200818_134119',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200821_100626',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200825_132217',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200826_163558',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200827_154304',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200904_142528',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200910_170840',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200916_155143',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201209_144624',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201210_152147',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201217_163659',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201222_152605',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210407_135410',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210408_145158',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210409_130407',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210414_131258',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210415_131337',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210422_161822',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210428_152350',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210504_112529',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210506_113444',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210507_110813',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210510_152309',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210512_122400',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210517_101440',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210519_113603',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210520_123525',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210521_121322',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210525_110414',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210526_112217',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210528_100147',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210531_094913']
        
   
LL2_impedances = np.zeros((len(LL2_impedance_folders),64))
LL2_phases = np.zeros((len(LL2_impedance_folders),64))
for i,folder in enumerate(LL2_impedance_folders):
    print(folder)
    impedances, phases = read_impedances_from_metadata(folder)
    LL2_impedances[i,:] = impedances
    LL2_phases[i,:] = phases

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
LL2_unit_folders = ['/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200818_134119',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200821_100626',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200825_132217',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200826_163558',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200827_154304',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200904_142528',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200910_170840',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_200916_155143',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201209_144624',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201210_152147',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201217_163659',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_201222_152605',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210407_135410',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210408_145158',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210409_130407',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210422_161822',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210517_101440',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210519_113603',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210520_123525',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210521_121322',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210525_110414',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210526_112217',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210528_100147',
           '/media/baran/Elements/Data/Chris_data_sorted/mouseLL2/mouseLL2_210531_094913'
    ]

days_LL2 = [21,24,28,29,30,38,44,50,134,135,142,147,253,254,255,268,293,295,296,297,301,302,304,307]

LL2_mean_SNRs = np.zeros(len(LL2_unit_folders))
LL2_SNR_errs = np.zeros(len(LL2_unit_folders))
single_unit_counts = np.zeros(len(LL2_unit_folders))

for i,folder in enumerate(LL2_unit_folders):
    file = folder + '/amplifier_data_res.mat'
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
    

mouseLL2_lifetimes = [54,277,11,39,257,173,263,160,8,5,34,33,48,33,150,42,218,2,108,13,97,13,13,118,92,8,6,12,12,14,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

LL2_impedances = np.zeros((len(LL2_unit_folders),64))
LL2_phases = np.zeros((len(LL2_unit_folders),64))
for i,folder in enumerate(LL2_unit_folders):
    print(folder)
    impedances, phases = read_impedances_from_metadata(folder)
    LL2_impedances[i,:] = impedances
    LL2_phases[i,:] = phases

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
