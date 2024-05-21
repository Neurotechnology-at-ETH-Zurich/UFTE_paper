#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 09:04:04 2022

@author: baran
"""
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

impedances = pd.read_csv('rat1_impedances.csv')
phases = pd.read_csv('rat1_impedance_phases.csv')
impedances = impedances.to_numpy()
phases = phases.to_numpy()
impedances = impedances[1:,:]
phases = phases[1:,:]
days = impedances[:,0]
days_35 = days
impedances = impedances[:,1:]
phases = phases[:,1:]

impedances = np.delete(impedances, np.arange(128,192), axis=1)
phases = np.delete(phases, np.arange(128,192), axis=1)

broken_channels = np.where(impedances[0] > 5e6)[0]
non_broken_channels = np.delete(impedances[0], broken_channels)
non_broken_channels_phases = np.delete(phases[0], broken_channels)
bad_channels = np.where(np.abs(stats.zscore(non_broken_channels))>3)[0][0]
good_channels = np.delete(non_broken_channels, bad_channels)
good_channels_phases = np.delete(non_broken_channels_phases, bad_channels)
good_channel_idx = np.where(np.isin(impedances[0], good_channels))[0]

post_surgery_impedances = impedances[1:,good_channel_idx]
post_surgery_phases = phases[1:,good_channel_idx]
post_surgery_nonbroken_impedances = {}
post_surgery_nonbroken_phases = {}

for day in range(len(days[1:])):
    broken_channels = np.where(post_surgery_impedances[day] > 5e6)[0]
    non_broken_channels = np.delete(post_surgery_impedances[day], broken_channels)
    non_broken_channels_phases = np.delete(post_surgery_phases[day], broken_channels)
    non_broken_channel_idx = np.where(np.isin(post_surgery_impedances[day], non_broken_channels))[0]
    post_surgery_nonbroken_impedances[days[day+1]] = post_surgery_impedances[day,non_broken_channel_idx]
    post_surgery_nonbroken_phases[days[day+1]] = post_surgery_phases[day,non_broken_channel_idx]
    
impedance_stats_35 = np.zeros((6,len(days)))
phase_stats_35 = np.zeros((6,len(days)))

for i in range(len(days)):
    if i == 0:
        impedance_stats_35[0,i] = np.median(good_channels)
        impedance_stats_35[1,i] = np.quantile(good_channels,0.25)
        impedance_stats_35[2,i] = np.quantile(good_channels,0.75)
        impedance_stats_35[3,i] = np.mean(good_channels)
        impedance_stats_35[4,i] = np.std(good_channels)
        impedance_stats_35[5,i] = len(good_channels)
        
        phase_stats_35[0,i] = np.median(good_channels_phases)
        phase_stats_35[1,i] = np.quantile(good_channels_phases,0.25)
        phase_stats_35[2,i] = np.quantile(good_channels_phases,0.75)
        phase_stats_35[3,i] = np.mean(good_channels_phases)
        phase_stats_35[4,i] = np.std(good_channels_phases)
        phase_stats_35[5,i] = len(good_channels_phases)
        
    else:
        impedance_stats_35[0,i] = np.median(post_surgery_nonbroken_impedances[days[i]])
        impedance_stats_35[1,i] = np.quantile(post_surgery_nonbroken_impedances[days[i]],0.25)
        impedance_stats_35[2,i] = np.quantile(post_surgery_nonbroken_impedances[days[i]],0.75)
        impedance_stats_35[3,i] = np.mean(post_surgery_nonbroken_impedances[days[i]])
        impedance_stats_35[4,i] = np.std(post_surgery_nonbroken_impedances[days[i]])
        impedance_stats_35[5,i] = len(post_surgery_nonbroken_impedances[days[i]])
    
        phase_stats_35[0,i] = np.median(post_surgery_nonbroken_phases[days[i]])
        phase_stats_35[1,i] = np.quantile(post_surgery_nonbroken_phases[days[i]],0.25)
        phase_stats_35[2,i] = np.quantile(post_surgery_nonbroken_phases[days[i]],0.75)
        phase_stats_35[3,i] = np.mean(post_surgery_nonbroken_phases[days[i]])
        phase_stats_35[4,i] = np.std(post_surgery_nonbroken_phases[days[i]])
        phase_stats_35[5,i] = len(post_surgery_nonbroken_phases[days[i]])
    
#%%   

impedances = pd.read_csv('rat2_impedances.csv')
phases = pd.read_csv('rat2_impedance_phases.csv')
impedances = impedances.to_numpy()
phases = phases.to_numpy()
impedances = impedances[1:,:]
phases = phases[1:,:]
days = impedances[:,0]
days_37 = days
impedances = impedances[:,1:]
phases = phases[:,1:]

broken_channels = np.where(impedances[0] > 5e6)[0]
non_broken_channels = np.delete(impedances[0], broken_channels)
non_broken_channels_phases = np.delete(phases[0], broken_channels)
bad_channels = np.where(np.abs(stats.zscore(non_broken_channels))>3)[0][0]
good_channels = np.delete(non_broken_channels, bad_channels)
good_channels_phases = np.delete(non_broken_channels_phases, bad_channels)
good_channel_idx = np.where(np.isin(impedances[0], good_channels))[0]

post_surgery_impedances = impedances[1:,good_channel_idx]
post_surgery_phases = phases[1:,good_channel_idx]
post_surgery_nonbroken_impedances = {}
post_surgery_nonbroken_phases = {}

for day in range(len(days[1:])):
    broken_channels = np.where(post_surgery_impedances[day] > 5e6)[0]
    non_broken_channels = np.delete(post_surgery_impedances[day], broken_channels)
    non_broken_channels_phases = np.delete(post_surgery_phases[day], broken_channels)
    non_broken_channel_idx = np.where(np.isin(post_surgery_impedances[day], non_broken_channels))[0]
    post_surgery_nonbroken_impedances[days[day+1]] = post_surgery_impedances[day,non_broken_channel_idx]
    post_surgery_nonbroken_phases[days[day+1]] = post_surgery_phases[day,non_broken_channel_idx]
    
impedance_stats_37 = np.zeros((6,len(days)))
phase_stats_37 = np.zeros((6,len(days)))

for i in range(len(days)):
    if i == 0:
        impedance_stats_37[0,i] = np.median(good_channels)
        impedance_stats_37[1,i] = np.quantile(good_channels,0.25)
        impedance_stats_37[2,i] = np.quantile(good_channels,0.75)
        impedance_stats_37[3,i] = np.mean(good_channels)
        impedance_stats_37[4,i] = np.std(good_channels)
        impedance_stats_37[5,i] = len(good_channels)
        
        phase_stats_37[0,i] = np.median(good_channels_phases)
        phase_stats_37[1,i] = np.quantile(good_channels_phases,0.25)
        phase_stats_37[2,i] = np.quantile(good_channels_phases,0.75)
        phase_stats_37[3,i] = np.mean(good_channels_phases)
        phase_stats_37[4,i] = np.std(good_channels_phases)
        phase_stats_37[5,i] = len(good_channels_phases)
        
    else:
        impedance_stats_37[0,i] = np.median(post_surgery_nonbroken_impedances[days[i]])
        impedance_stats_37[1,i] = np.quantile(post_surgery_nonbroken_impedances[days[i]],0.25)
        impedance_stats_37[2,i] = np.quantile(post_surgery_nonbroken_impedances[days[i]],0.75)
        impedance_stats_37[3,i] = np.mean(post_surgery_nonbroken_impedances[days[i]])
        impedance_stats_37[4,i] = np.std(post_surgery_nonbroken_impedances[days[i]])
        impedance_stats_37[5,i] = len(post_surgery_nonbroken_impedances[days[i]])
    
        phase_stats_37[0,i] = np.median(post_surgery_nonbroken_phases[days[i]])
        phase_stats_37[1,i] = np.quantile(post_surgery_nonbroken_phases[days[i]],0.25)
        phase_stats_37[2,i] = np.quantile(post_surgery_nonbroken_phases[days[i]],0.75)
        phase_stats_37[3,i] = np.mean(post_surgery_nonbroken_phases[days[i]])
        phase_stats_37[4,i] = np.std(post_surgery_nonbroken_phases[days[i]])
        phase_stats_37[5,i] = len(post_surgery_nonbroken_phases[days[i]])

#%%    
    
color1 = 'r'
color2 = 'b'    
    
plt.figure()
plt.plot(days_35[1:],impedance_stats_35[0,1:],color=color1)
plt.fill_between(days_35[1:],impedance_stats_35[1,1:],impedance_stats_35[2,1:],alpha=0.1,color=color1)
plt.plot(days_37[1:],impedance_stats_37[0,1:],color=color2)
plt.fill_between(days_37[1:],impedance_stats_37[1,1:],impedance_stats_37[2,1:],alpha=0.1,color=color2)
plt.yscale('log')
plt.ylim([1e5,1e7])
plt.ylabel('Impedances (Ohm)')
plt.xlabel('Days')
plt.show()    

color1 = 'r'
color2 = 'b'    
    
plt.figure()
plt.plot(days_35[1:],phase_stats_35[0,1:],color=color1)
plt.fill_between(days_35[1:],phase_stats_35[1,1:],phase_stats_35[2,1:],alpha=0.1,color=color1)
plt.plot(days_37[1:],phase_stats_37[0,1:],color=color2)
plt.fill_between(days_37[1:],phase_stats_37[1,1:],phase_stats_37[2,1:],alpha=0.1,color=color2)
plt.ylim([-90,0])
plt.ylabel('Phase')
plt.xlabel('Days')
plt.show()    

