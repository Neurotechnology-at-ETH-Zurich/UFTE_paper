#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:44:13 2024

@author: baran
"""

import mat73
import numpy as np 
from matplotlib import pyplot as plt
from bisect import bisect_left
import matplotlib

def return_unit_lifetimes(spikeTimes, spikesByCluster, num_units, days, time_boundaries, to_delete): 
    """
    This function calculates the tracking durations of each cluster based on the spike sorting results. 

    Parameters
    ----------
    spikeTimes : Array of int32
        Array from JRClust containing times of each spike in samples (loaded from res.mat file and converted to Python indexing by subtracting 1).
    spikesByCluster : list
        List of arrays from JRClust containing the ids of spikes belonging to each cluster.
    num_units : int
        Total number of units.
    days : Array of int64
        Array containing the post-implantation days of the recording sessions.
    time_boundaries : Array of float64
        Array containing the boundaries between consecutive sessions in the concatenated data.
    to_delete : Array of int64
        Array containing the ids of isi-violating single units to ignore.

    Returns
    -------
    unit_lifetimes : Array of int
        (1xN) array containing the tracked durations of each of the N units in days .
    unit_lifetime_matrix : Array of int
        (NxM) array showing in which of the M sessions each of the N units were detected .

    """    
    unit_lifetime_matrix = np.zeros((num_units,len(days)),dtype=int)
    unit_lifetimes = np.zeros(num_units, dtype=int)
    
    #Iterating over each unit and counting their spikes in each of the sessions 
    for unit in range(num_units):
        unit_spike_idx = spikesByCluster[unit][0].astype(int) - 1 #converting to Python by subtracting minus 1
        unit_spike_times = spikeTimes[unit_spike_idx] / (sr)
        unit_spike_sessions = np.zeros(len(unit_spike_times))
        
        #Classifying the unit spikes into sessions 
        for i, spike_time in enumerate(unit_spike_times):
            unit_spike_sessions[i] = bisect_left(time_boundaries[1:],spike_time)
        unit_spike_sessions = unit_spike_sessions.astype('int')
        counts = np.bincount(unit_spike_sessions)
        counts[counts<10] = 0 #A unit is considered to be appearing in a session only if it has more than 10 spikes
        
        #Constructing the unit_lifetime_matrix
        unit_lifetime_matrix[unit][np.where(counts)[0]] = 1
        if to_delete != []:
            unit_lifetime_matrix[to_delete[:,0],to_delete[:,1]]=0
        unit_days = days[np.where(unit_lifetime_matrix[unit]==1)[0]]
        try:
            unit_lifetimes[unit] = unit_days[-1] - unit_days[0] + 1
        except IndexError:
            unit_lifetimes[unit] = 0 
            
    return unit_lifetimes, unit_lifetime_matrix

#Defining parameters 
sr = 20000 #Sample rate (Hz)
mins = 20 #Minutes taken from each session 

#Post-implantation days of the recording sessions 
#35 and 37 refer to the rat IDs throughout the code
days_35 = np.array([11,13,14,18,20,25,27,33,35,41,43,47,49,54,57,61,64,68,71,75,77,82,84,89,92,95])
days_37 = np.array([3,7,9,14,17,21,24,28,31,35,37,42,44,49,52,56,59,63,66,70,76,80,84,87,90])

#Boundaries between consecutive sessions in the concatenated files. 
time_boundaries_35 = np.arange(0, 60*mins*(len(days_35)+1), 60*mins)
time_boundaries_37 = np.arange(984.9152, 984.9152+60*mins*len(days_37), 60*mins)
time_boundaries_37 = np.insert(time_boundaries_37, 0, 0)

#Single units violating ISI 
to_delete_35 = np.array([[28,0],[28,2],[28,3],[28,6],[39,4],[39,5],[40,7],[41,11],[117,13],[139,7],[143,13],[241,17],[267,18],[337,19]])
to_delete_37 = np.array([[48,12],[48,13],[79,12],[82,7],[98,0],[98,13],[126,8],[130,21],[151,21],[156,18],[252,0],[264,5],[265,18],[273,18],[283,0],[283,3],[283,4],[283,9],[305,1],[309,5],[317,21],[320,10],[354,5],[354,6],[354,11]])

#loading spike sorting results
data_dict_35 = mat73.loadmat('/home/baran/Multiarea_bundle_data/rTBY35/combined/combined_res.mat')
data_dict_37 = mat73.loadmat('/home/baran/Multiarea_bundle_data/rTBY37/combined/except_RSC/combined_res.mat')
data_dict_37_rsc = mat73.loadmat('/home/baran/Multiarea_bundle_data/rTBY37/combined/RSC/combined_res.mat')

#Loading JRClust parameter files
prm_file_35 = '/home/baran/Multiarea_bundle_data/rTBY35/combined/combined.prm'
prm_file_37 = '/home/baran/Multiarea_bundle_data/rTBY37/combined/except_RSC/combined.prm'
prm_file_37_rsc = '/home/baran/Multiarea_bundle_data/rTBY37/combined/RSC/combined.prm'

#Loading the parameters from the JRClust result files
spikeTimes_35 = data_dict_35['spikeTimes'] - 1 #Times of each spike in samples (python indexing)
spikesByCluster_35 = data_dict_35['spikesByCluster'] #spike IDs belonging to each unit 
num_units_35 = len(spikesByCluster_35) #Number of units detected in this animal

#Loading the same parameters for Rat 37
spikeTimes_37 = data_dict_37['spikeTimes'] - 1
spikesByCluster_37 = data_dict_37['spikesByCluster']
num_units_37 = len(spikesByCluster_37)

#Loading the parameters for Rat 37 RSC since it was done separately
spikeTimes_37_rsc = data_dict_37_rsc['spikeTimes'] - 1
spikesByCluster_37_rsc = data_dict_37_rsc['spikesByCluster']
num_units_37_rsc = len(spikesByCluster_37_rsc)
clusterNotes_37_rsc = np.array(data_dict_37_rsc['clusterNotes'])

unit_lifetimes_35, unit_lifetime_matrix_35 = return_unit_lifetimes(spikeTimes_35, spikesByCluster_35, num_units_35, days_35, time_boundaries_35, to_delete_35)
unit_lifetimes_37, unit_lifetime_matrix_37 = return_unit_lifetimes(spikeTimes_37, spikesByCluster_37, num_units_37, days_37, time_boundaries_37, to_delete_37)
unit_lifetimes_37_rsc, unit_lifetime_matrix_37_rsc = return_unit_lifetimes(spikeTimes_37_rsc, spikesByCluster_37_rsc, num_units_37_rsc, days_37,time_boundaries_37, [])
unit_lifetime_matrix_37 = np.vstack((unit_lifetime_matrix_37, unit_lifetime_matrix_37_rsc))

#%%
#Loading the files with the ensemble activation and further ensemble information
ens_info_35 = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Matlab_35.mat')
ens_info_37 = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Matlab_37.mat')
ens_activations_35_info_1 = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Matlab_35_Assembly_activations.mat')
ens_activations_37_info_1 = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Matlab_37_Assembly_activations.mat')
ens_activations_35_info = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Matlab_35_Assembly_activations_bigger_bin_smoothed.mat')
ens_activations_37_info = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Matlab_37_Assembly_activations_bigger_bin_smoothed.mat')
ens_activations_ripple_35 = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Assemby_Activation_During_Ripples_35.mat')['Assembly_Activation_During_Ripple']
ens_activations_ripple_37 = mat73.loadmat('/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 5/Assemby_Activation_During_Ripples_37.mat')['Assembly_Activation_During_Ripple']

#Activations of ensembles over time
ens_act_35 = ens_activations_35_info['Activity_strength_smoothed']
ens_act_37 = ens_activations_37_info['Activity_strength_smoothed']

#IDs of units that are part of ensembles
ens_cellID_35 = (ens_info_35['Assembly_cellID'] - 1).astype('int')
ens_cellID_37 = (ens_info_37['Assembly_cellID'] - 1).astype('int')

#Participations of units in ensembles
ens_cell_part_35 = ens_info_35['Significant_Neurons_Assemblies']
ens_cell_part_37 = ens_info_37['Significant_Neurons_Assemblies']

#Putative brain areas of units in ensembles
ens_cell_structures_35 = np.array(ens_info_35['PutitativeStructure'])
ens_cell_structures_37 = np.array(ens_info_37['PutitativeStructure'])

#Lifetimes of units that are parts of ensembles
ens_cell_lifetimes_35 = unit_lifetime_matrix_35[ens_cellID_35]
ens_cell_lifetimes_37 = unit_lifetime_matrix_37[ens_cellID_37]

#What percentage of units are detected in each ensembe in each session
ens_cell_part_perc_35 = np.zeros((ens_cell_part_35.shape[1], len(days_35)))
ens_cell_part_perc_37 = np.zeros((ens_cell_part_37.shape[1], len(days_37)))

#Distribution of units in ensembles based on their brain regions
ens_cell_structure_perc_35 = np.zeros((ens_cell_part_35.shape[1],2))
ens_cell_structure_perc_37 = np.zeros((ens_cell_part_37.shape[1],3))

#Activation of ensembles in each session
ens_act_ses_35 = np.zeros((len(days_35), ens_cell_part_35.shape[1]))
ens_act_ses_37 = np.zeros((len(days_37), ens_cell_part_37.shape[1]))

#Figuring out which bin in the ensemble activation traces correspond to which sessions
rel_times_edges_35 = ens_activations_35_info_1['center_of_relative_times_edges']
acu_time_35 = ens_activations_35_info['bin_activities']
bin_sessions_35 = np.arange(len(acu_time_35))
for i, sec in enumerate(acu_time_35):
    bin_sessions_35[i] = bisect_left(rel_times_edges_35[1:],sec)

rel_times_edges_37 = ens_activations_37_info_1['center_of_relative_times_edges']
acu_time_37 = ens_activations_37_info['bin_activities']
bin_sessions_37 = np.arange(len(acu_time_37))
for i, sec in enumerate(acu_time_37):
    bin_sessions_37[i] = bisect_left(rel_times_edges_37[1:],sec)

#Iterating over the ensembles to populate the arrays for ensemble cell participation, ensemble activation by session and ensemble brain structure compositions
for ensemble in range(ens_cell_part_35.shape[1]):
    ens_cell_part_perc_35[ensemble] = np.mean(ens_cell_lifetimes_35[ens_cell_part_35[:,ensemble]], axis=0)
    ens_act_ses_35[np.unique(bin_sessions_35[np.where(ens_act_35[ensemble,:] > 2 * np.std(ens_act_35[ensemble,:]))[0]]), ensemble] = 1
    x = ens_cell_structures_35[ens_cell_part_35[:,ensemble]]
    ens_cell_structure_perc_35[ensemble,0] = len(np.where(np.logical_or((x=='dHP'),(x=='iHP')))[0]) / len(x)
    ens_cell_structure_perc_35[ensemble,1] = len(np.where(np.logical_or((x=='PrL'),(x=='Cg1')))[0]) / len(x)
    
for ensemble in range(ens_cell_part_37.shape[1]):
    ens_cell_part_perc_37[ensemble] = np.mean(ens_cell_lifetimes_37[ens_cell_part_37[:,ensemble]], axis=0)
    ens_act_ses_37[np.unique(bin_sessions_37[np.where(ens_act_37[ensemble,:] > 2 * np.std(ens_act_37[ensemble,:]))[0]]), ensemble] = 1
    x = ens_cell_structures_37[ens_cell_part_37[:,ensemble]]
    ens_cell_structure_perc_37[ensemble,0] = len(np.where(np.logical_or((x=='dHP'),(x=='iHP')))[0]) / len(x)
    ens_cell_structure_perc_37[ensemble,1] = len(np.where(np.logical_or(np.logical_or((x=='PrL'),(x=='Cg1')),(x=='IL')))[0]) / len(x)
    ens_cell_structure_perc_37[ensemble,2] = len(np.where(x=='RSC')[0]) / len(x)

ens_thresh = 0.66 #the threshold of the percentage of tracked units at which the ensemble is considered to be "alive"

#getting the info of in how many sessions each ensemble has more members detected than the ens_thresh
ens_cell_part_perc_35_bin = np.zeros((ens_cell_part_perc_35.shape))
ens_cell_part_perc_35_bin[np.where(ens_cell_part_perc_35 > ens_thresh)] = 1
ens_cell_part_perc_35_bin = np.transpose(ens_cell_part_perc_35_bin)
ens_cell_part_perc_37_bin = np.zeros((ens_cell_part_perc_37.shape))
ens_cell_part_perc_37_bin[np.where(ens_cell_part_perc_37 > ens_thresh)] = 1
ens_cell_part_perc_37_bin = np.transpose(ens_cell_part_perc_37_bin)

#Calculating the lifetimes of ensembles based on activation and the trackability of neurons 
ens_cell_lifetimes_35 = np.zeros(ens_cell_part_35.shape[1])
ens_cell_lifetimes_37 = np.zeros(ens_cell_part_37.shape[1])
ens_lifetimes_35 = np.zeros(ens_cell_part_35.shape[1])
ens_lifetimes_37 = np.zeros(ens_cell_part_37.shape[1])

for ensemble in range(ens_cell_part_35.shape[1]):
    if len(np.where(ens_cell_part_perc_35_bin[:,ensemble])[0]) > 1:
        ens_cell_lifetimes_35[ensemble] = days_35[np.where(ens_cell_part_perc_35_bin[:,ensemble])[0][-1]] - days_35[np.where(ens_cell_part_perc_35_bin[:,ensemble])[0][0]]
    if len(np.where(ens_act_ses_35[:,ensemble])[0]) > 1:
        ens_lifetimes_35[ensemble] = days_35[np.where(ens_act_ses_35[:,ensemble])[0][-1]] - days_35[np.where(ens_act_ses_35[:,ensemble])[0][0]]
    
for ensemble in range(ens_cell_part_37.shape[1]):
    if len(np.where(ens_cell_part_perc_37_bin[:,ensemble])[0]) > 1:
        ens_cell_lifetimes_37[ensemble] = days_37[np.where(ens_cell_part_perc_37_bin[:,ensemble])[0][-1]] - days_37[np.where(ens_cell_part_perc_37_bin[:,ensemble])[0][0]]
    if len(np.where(ens_act_ses_37[:,ensemble])[0]) > 1:
        ens_lifetimes_37[ensemble] = days_37[np.where(ens_act_ses_37[:,ensemble])[0][-1]] - days_37[np.where(ens_act_ses_37[:,ensemble])[0][0]]    
        
ens_cell_lifetimes = np.concatenate((ens_cell_lifetimes_35, ens_cell_lifetimes_37))
ens_lifetimes = np.concatenate((ens_lifetimes_35, ens_lifetimes_37))

#Calculating the percentage of hippocampal neurons in the ensembles
ens_hpc_part = np.zeros(len(ens_cell_structure_perc_35)+len(ens_cell_structure_perc_37))
for i in range(len(ens_cell_structure_perc_35)):
    hpc_part = ens_cell_structure_perc_35[i][0]
    ens_hpc_part[i] = hpc_part
        
for i in range(len(ens_cell_structure_perc_37)):
    hpc_part = ens_cell_structure_perc_37[i][0]
    ens_hpc_part[i+len(ens_cell_structure_perc_35)] = hpc_part

#Calculating the activations of the ensembles during ripples where the ripples activate some ensembles
a = np.where(ens_activations_ripple_35 > np.mean(ens_activations_ripple_35)+2*np.std(ens_activations_ripple_35))
b = np.zeros(ens_activations_ripple_35.shape)
b[a[0],a[1]] = 1
non_zero_ripples_35 = np.where(np.sum(b, axis= 0))[0]
ens_activations_ripple_35_st = ens_activations_ripple_35[:,non_zero_ripples_35]

a = np.where(ens_activations_ripple_37 > np.mean(ens_activations_ripple_37) + 2*np.std(ens_activations_ripple_37))
b = np.zeros(ens_activations_ripple_37.shape)
b[a[0],a[1]] = 1
non_zero_ripples_37 = np.where(np.sum(b, axis= 0))[0]
ens_activations_ripple_37_st = ens_activations_ripple_37[:,non_zero_ripples_37]

ens_activations_ripple_35_mean = np.mean(ens_activations_ripple_35_st, axis=1)
ens_activations_ripple_37_mean = np.mean(ens_activations_ripple_37_st, axis=1)
ens_activations_ripple_mean = np.concatenate((ens_activations_ripple_35_mean, ens_activations_ripple_37_mean))

#%%

#Generating the figures
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("",["#4C4C4C","#76AC42"])

plt.figure("Fig. 5D")
plt.violinplot(ens_lifetimes, positions=[1], showmeans=True)
plt.violinplot(ens_cell_lifetimes, positions=[2], showmeans=True)
plt.show()

plt.figure("Fig. 5F-scatter")
plt.scatter(ens_lifetimes, ens_activations_ripple_mean, c=ens_hpc_part, cmap=cmap, s=50)
plt.axvline(x=np.mean(ens_lifetimes) + np.std(ens_lifetimes))
plt.axhline(y=np.mean(ens_activations_ripple_mean) + np.std(ens_activations_ripple_mean))
plt.colorbar()
plt.show()

sl = ens_lifetimes < (np.mean(ens_lifetimes) + np.std(ens_lifetimes)) #short life ensembles 
ll = ens_lifetimes > (np.mean(ens_lifetimes) + np.std(ens_lifetimes)) #long life ensembles
ha = ens_activations_ripple_mean > (np.mean(ens_activations_ripple_mean) + np.std(ens_activations_ripple_mean)) #high ripple activated ensembles
la = ens_activations_ripple_mean < (np.mean(ens_activations_ripple_mean) + np.std(ens_activations_ripple_mean)) #low ripple activated ensembles

plt.figure("Fig. 5F-boxplot left")
plt.boxplot(ens_activations_ripple_mean[ll], positions = [1])
plt.boxplot(ens_activations_ripple_mean[np.logical_and(sl,ha)], positions = [2])
plt.show()

plt.figure("Fig. 5F-boxplot right")
plt.boxplot(ens_lifetimes[np.logical_and(la,ll)], positions = [1])
plt.boxplot(ens_lifetimes[ha], positions = [2])
plt.show()
