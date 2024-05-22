#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:42:12 2024

@author: baran
"""

#Importing packages
import numpy as np
import mat73
import pandas as pd
from tqdm import tqdm
from sklearn.decomposition import PCA
from scipy import stats
from scipy import signal
from matplotlib import pyplot as plt
from bisect import bisect_left

#Defining functions
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

def calculate_global_PCA(spikeSites, spikesFilt_rs, global_neigh, total_num_channels):
    """
    This function reduces the dimensionality of all spike waveforms detected in the concatenated session.

    Parameters
    ----------
    spikeSites : Array of int32
        Array from JRClust indicating the primary sites of each spike (subtracted 1 to convert to Python indexing).
    spikesFilt_rs : Array of float64
        KxMxN array from JRClust file (_filt.jrc) containing the filtered waveforms of each spike in K waveform samples, M recording contacts, N spikes
    global_neigh : Array of int32
        MxN array containing the neighboring channels of each recording contact, starting from the recording contact. M = number of contacts, N = number of neighbors
    total_num_channels : int
        Total number of recording contacts in the electrode arrays (in all bundles).

    Returns
    -------
    ch_spikePCA_all : Array of float64
        MxN array containing the dimensionality-reduced spike waveforms where M corresponds to number of PCsxnumber of recording contacts in bundle, N corresponds to numberr of spikes.

    """
    num_PCs = 3  
    ch_spikePCA_all = np.zeros((num_PCs*total_num_channels,len(spikeSites)))
    chs = np.arange(1,total_num_channels+1)
    for ch in tqdm(chs):
        #For each spike, figuring out which relative channel index in JRC nomenclature corresponds to the channel of 
        #interest. If the spike is outside of the neighborhood of "ch", then the index is -1
        rel_spikeSites = np.ones(len(spikeSites), dtype='int32') * -1
        rel_spikeSites[np.where(global_neigh[spikeSites] == ch)[0]] = np.where(global_neigh[spikeSites] == ch)[1]
        
        #determining which spikes fall under the neighborhood of ch, and extracting their waveforms at site "ch"
        ch_spikes = np.where(rel_spikeSites != -1)[0]
        ch_spikeWaveforms = spikesFilt_rs[:,rel_spikeSites[ch_spikes],ch_spikes]
        
        #performing PCA on the spike waveforms detected in site "ch" and calculating the projections of the spike 
        #waveforms on these PC axes
        if ch_spikes.size > (num_PCs * total_num_channels):
            pca = PCA(n_components=num_PCs)
            ch_spikePCA = pca.fit_transform(ch_spikeWaveforms.T)
            ch_spikePCA_all[((ch-1)*num_PCs):(ch*num_PCs),ch_spikes] = ch_spikePCA.T
            
    return ch_spikePCA_all

def get_global_neigh(total_num_channels, rec_sites, prm_file, num_shanks, channel_per_shank, ignoreSites):
    """
    This function gets the "neighboring channels" of each recording contact as defined by the JRClust nomenclature, skipping broken channels. All calculation is done based on MATLAB indexing (starting from 1) in this function.

    Parameters
    ----------
    total_num_channels : int
        Total number of recording contacts in the electrode arrays (in all bundles).
    rec_sites : int
        Number of recording contacts in the neighboring of primary contact.
    prm_file : str
        Location of the .prm file.
    num_shanks : int
        Number of bundles.
    channel_per_shank : int
        Number of recording contacts in bundle.
    ignoreSites : Array of int64
        Array containing the IDs of broken channels.

    Returns
    -------
    global_neigh : Array of int32
        MxN array containing the neighboring channels of each recording contact, starting from the recording contact. M = number of contacts, N = number of neighbors

    """
    global_neigh = np.zeros((total_num_channels, rec_sites))
    for shank in range(num_shanks): 
        shank_global_neigh = np.zeros((channels_per_shank,rec_sites))
        first_ch = shank * channels_per_shank + 1
        last_ch = (shank+1) * channels_per_shank 
        chs = np.arange(first_ch, last_ch+1)
        dead_sites = ignoreSites[np.logical_and(np.greater_equal(ignoreSites,first_ch), np.less_equal(ignoreSites,last_ch))] #read the dead channel info for the shank

        for i, ch in enumerate(chs):
            #populate a matrix with the channel neighborhood matrix, starting from first channel and moving to the immediate neighbors
            #e.g. for channel 25, it goes as follows: [25,24,26,23,27,22,28...]
            neigh_array = np.zeros(rec_sites*4+1)
            neigh_array[0] = ch
            neigh_array[1::2] = np.arange(ch-1,ch-(rec_sites*2+1),-1)
            neigh_array[2::2] = np.arange(ch+1,ch+(rec_sites*2+1),+1)
            
            #clean up the channels that are outside of the shank range or were ignored by JRClust
            neigh_array = np.delete(neigh_array, np.where(neigh_array < first_ch)[0])
            neigh_array = np.delete(neigh_array, np.where(neigh_array > last_ch)[0])
            neigh_array = np.delete(neigh_array, np.isin(neigh_array, dead_sites))
            if neigh_array.size > 0:
                shank_global_neigh[i] = neigh_array[0:rec_sites]
        
        global_neigh[first_ch-1:last_ch] = shank_global_neigh
            
    global_neigh = global_neigh.astype('int32')
    return global_neigh

def calculate_meanWf_meanPC(days, num_units, spikesFilt_rs, ch_spikePCA_all, unit_lifetime_matrix, spikesByCluster, spike_sessions, clusterSites, global_neigh):
    """
    This function gets the mean waveforms and the mean and variance of PCs for each unit in each session

    Parameters
    ----------
    days : Array of int64
        Array containing the post-implantation days of the recording sessions.
    num_units : int
       Total number of units.
     spikesFilt_rs : Array of float64
         KxMxN array from JRClust file (_filt.jrc) containing the filtered waveforms of each spike in K waveform samples, M recording contacts, N spikes
    ch_spikePCA_all : Array of float64
        MxN array containing the dimensionality-reduced spike waveforms where M corresponds to number of PCsxnumber of recording contacts in bundle, N corresponds to numberr of spikes.
    unit_lifetime_matrix : Array of int
        (NxM) array showing in which of the M sessions each of the N units were detected.
    spikesByCluster : list
        List of arrays from JRClust containing the ids of spikes belonging to each cluster.
    spike_sessions : Array of int
        1xN array showing to which session each spike belongs to (N=number of spikes).
    clusterSites : Array of float64
        1xN array showing the primary site of each unit, N=number of units.
    global_neigh : Array of int32
        MxN array containing the neighboring channels of each recording contact, starting from the recording contact. M = number of contacts, N = number of neighbors

    Returns
    -------
    meanWfs : Array of float64
        KxLxMxN array containing the mean waveforms of each unit in each recording contact in each session. K: number of units, L: number of sessions, M: number of recording contacts, N: waveform samples
    meanPCs :  Array of float64
        KxLxM array containing the mean PCs of each unit in each session. K: number of units, L: number of sessions, M: number of total PCs=number of PC per recording contact x number of recording contacts
    var_PCs :  Array of float64
        KxLxM array containing the variances of PCs of each unit in each session. K: number of units, L: number of sessions, M: number of total PCs=number of PC per recording contact x number of recording contacts.

    """
    meanWfs = np.zeros((unit_lifetime_matrix.shape[0], unit_lifetime_matrix.shape[1], total_num_channels, 41))
    meanPCs = np.zeros((unit_lifetime_matrix.shape[0], unit_lifetime_matrix.shape[1], total_num_channels*3))
    var_PCs = np.zeros((unit_lifetime_matrix.shape[0], unit_lifetime_matrix.shape[1], total_num_channels*3))
    for session in tqdm(range(len(days))):
        first_spike_session = np.where(spike_sessions == session)[0][0]
        last_spike_session = np.where(spike_sessions == session)[0][-1]
        spikesFilt_rs_session = spikesFilt_rs[:,:,first_spike_session:(last_spike_session+1)]
        ch_spikePCA_all_session = ch_spikePCA_all[:,first_spike_session:(last_spike_session+1)]
        for unit in range(num_units):
            if unit_lifetime_matrix[unit][session]==1:
                unit_spike_idx = spikesByCluster[unit][0].astype(int) - 1
                session_unit_spikes = np.where(np.isin(np.where(spike_sessions == session)[0], unit_spike_idx))[0]
                unit_session_PC = ch_spikePCA_all_session[:,session_unit_spikes]
                meanPCs[unit][session] = np.mean(unit_session_PC, axis=1)
                var_PCs[unit][session] = np.var(unit_session_PC, axis=1)
                unitWfs_session_mean = np.mean(spikesFilt_rs_session[:,:,session_unit_spikes], axis=2)
                for i,ch in enumerate(global_neigh[int(clusterSites[unit])]):
                    meanWfs[unit,session,(ch-1),:] = unitWfs_session_mean[:,i]   
    return meanWfs, meanPCs, var_PCs
#%%

#Defining parameters 
sr = 20000 #Sample rate (Hz)
mins = 20 #Minutes taken from each session 
total_num_channels = 256 #total number of recording contacts
num_shanks = 4 #Number of bundles in each animal 
channels_per_shank = 64 #Number of recording contacts per bundle 
rec_sites = 7 #Number of recording contacts considered in the neighborhood of each channel 
num_PCs = 3 #Number of principal components to be taken while reducing the dimensionality of the spike waveforms 

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
data_dict_35 = mat73.loadmat('rTBY35_combined_res.mat')
data_dict_37 = mat73.loadmat('rTBY37_except_RSC_combined_res.mat')
data_dict_37_rsc = mat73.loadmat('rTBY37_RSC_combined_res.mat')

#Loading JRClust parameter files
prm_file_35 = 'rTBY35_combined.prm'
prm_file_37 = 'rTBY37_except_RSC_combined.prm'
prm_file_37_rsc = 'rTBY37_RSC_combined.prm'

#Loading the parameters from the JRClust result files
spikeTimes_35 = data_dict_35['spikeTimes'] - 1 #Times of each spike in samples (python indexing)
spikesByCluster_35 = data_dict_35['spikesByCluster'] #spike IDs belonging to each unit 
clusterSites_35 = data_dict_35['clusterSites'] - 1 #Primary sites of each unit (python indexing)
spikeSites_35 = data_dict_35['spikeSites'] - 1 #Primary sites of each spike (python indexing)
num_units_35 = len(spikesByCluster_35) #Number of units detected in this animal
spike_sessions_35 = np.zeros(len(spikeTimes_35), dtype=int) #Identifying the sessions of each spike
for i, spike_time in enumerate(spikeTimes_35):
    spike_sessions_35[i] = bisect_left(time_boundaries_35[1:],spike_time/sr)
clusterNotes_35 = np.array(data_dict_35['clusterNotes']) #Array containing info on whether each unit is single or multi

#Loading the same parameters for Rat 37
spikeTimes_37 = data_dict_37['spikeTimes'] - 1
spikesByCluster_37 = data_dict_37['spikesByCluster']
clusterSites_37 = data_dict_37['clusterSites'] - 1
num_units_37 = len(spikesByCluster_37)
spikeSites_37 = data_dict_37['spikeSites'] - 1
spike_sessions_37 = np.zeros(len(spikeTimes_37),dtype=int)
for i, spike_time in enumerate(spikeTimes_37):
    spike_sessions_37[i] = bisect_left(time_boundaries_37[1:],spike_time/sr)
clusterNotes_37 = np.array(data_dict_37['clusterNotes'])

#Loading the parameters for Rat 37 RSC since it was done separately
spikeTimes_37_rsc = data_dict_37_rsc['spikeTimes'] - 1
clusterSites_37_rsc = data_dict_37_rsc['clusterSites'] - 1
spikesByCluster_37_rsc = data_dict_37_rsc['spikesByCluster']
num_units_37_rsc = len(spikesByCluster_37_rsc)
spikeSites_37_rsc = data_dict_37_rsc['spikeSites'] - 1
spike_sessions_37_rsc = np.zeros(len(spikeTimes_37_rsc),dtype=int)
for i, spike_time in enumerate(spikeTimes_37_rsc):
    spike_sessions_37_rsc[i] = bisect_left(time_boundaries_37[1:],spike_time/sr)
clusterNotes_37_rsc = np.array(data_dict_37_rsc['clusterNotes'])

total_num_channels = num_shanks * channels_per_shank

#Defining sites that were excluded from the spike sorting
ignoreSites_35 = np.r_[[12,52,64],np.arange(129,193),[246,251,256]]
ignoreSites_37 = np.r_[[1, 2, 9, 53], np.arange(129,193), [196, 210]]
ignoreSites_37_rsc = np.r_[np.arange(1,129),[131,145,147,157,182,183,184,185,186],np.arange(193,257)]

#Generating the global_neigh matrices
global_neigh_35 = get_global_neigh(total_num_channels, rec_sites, prm_file_35, num_shanks, channels_per_shank, ignoreSites_35)
global_neigh_37 = get_global_neigh(total_num_channels, rec_sites, prm_file_37, num_shanks, channels_per_shank, ignoreSites_37)
global_neigh_37_rsc = get_global_neigh(total_num_channels, rec_sites, prm_file_37_rsc, num_shanks, channels_per_shank, ignoreSites_37_rsc)

#Loading the unit quality metrics from the files generated by unit_quality_analysis.py
r35_quality = pd.read_csv('rTBY35_combined.csv')
r37_quality = pd.read_csv('rTBY37_combined.csv')
r37_quality_rsc = pd.read_csv('rTBY37_combined_RSC.csv')

#Gathering isi violation percentages
isi35 = np.array(r35_quality['ISI violation percentage'])
isi37 = np.array(r37_quality['ISI violation percentage'])
isi37_rsc = np.array(r37_quality_rsc['ISI violation percentage'])
isi37 = np.concatenate((isi37, isi37_rsc))
isi_all = np.concatenate((isi35,isi37))

#Gathering nearest neighbor distances
nn_distance_35 = np.array(r35_quality['Unit dists from nearest unit'])
nn_distance_37 = np.array(r37_quality['Unit dists from nearest unit'])
nn_distance_37_rsc = np.array(r37_quality_rsc['Unit dists from nearest unit'])
nn_distance_37 = np.concatenate((nn_distance_37, nn_distance_37_rsc))
nn_distances_all = np.concatenate((nn_distance_35, nn_distance_37))

#Identifying single and multi units based on manual labeling
clusterNotes_37 = np.concatenate((clusterNotes_37, clusterNotes_37_rsc))
single_units35 = np.squeeze(np.array(clusterNotes_35 == 'single'))
single_units35[np.where(single_units35)[0][np.where(isi35[single_units35] > 2)[0]]] = False
single_units35[28] = False #eliminated by to_delete
single_units37 = np.squeeze(np.array(clusterNotes_37 == 'single'))
single_units37[np.where(single_units37)[0][np.where(isi37[single_units37] > 2)[0]]] = False
single_units37[[48,98]] = False #eliminated by to_delete
single_units_all = np.concatenate((single_units35, single_units37))
multi_units_all = np.logical_not(single_units_all)
multi_units_all[[28,len(single_units35)+48,len(single_units35)+98]] = False #eliminated by to_delete

#Generating the unit tracking durations and the unit trackability matrices
unit_lifetimes_35, unit_lifetime_matrix_35 = return_unit_lifetimes(spikeTimes_35, spikesByCluster_35, num_units_35, days_35, time_boundaries_35, to_delete_35)
unit_lifetimes_37, unit_lifetime_matrix_37 = return_unit_lifetimes(spikeTimes_37, spikesByCluster_37, num_units_37, days_37, time_boundaries_37, to_delete_37)
unit_lifetimes_37_rsc, unit_lifetime_matrix_37_rsc = return_unit_lifetimes(spikeTimes_37_rsc, spikesByCluster_37_rsc, num_units_37_rsc, days_37,time_boundaries_37, [])
unit_lifetimes_37 = np.concatenate((unit_lifetimes_37, unit_lifetimes_37_rsc))
unit_lifetimes_all = np.concatenate((unit_lifetimes_35,unit_lifetimes_37))

single_unit_lifetimes35 = unit_lifetimes_35[single_units35]
single_unit_lifetimes37 = unit_lifetimes_37[single_units37]

#Generating the mean of waveforms and mean/variance of principal components
spikesFilt_dir = 'rTBY35_combined_filt.jrc'
spikesFilt = np.fromfile(spikesFilt_dir, 'int16') * 0.195
spikesFilt_rs = np.reshape(spikesFilt, [41,7,len(spikeTimes_35)], 'F')  
ch_spikePCA_all = calculate_global_PCA(spikeSites_35, spikesFilt_rs, global_neigh_35, total_num_channels)
meanWfs_35, meanPCs_35, varPCs_35 = calculate_meanWf_meanPC(days_35, num_units_35, spikesFilt_rs, ch_spikePCA_all, unit_lifetime_matrix_35, spikesByCluster_35, spike_sessions_35, clusterSites_35, global_neigh_35)

spikesFilt_dir = 'rTBY37_except_RSC_combined_filt.jrc'
spikesFilt = np.fromfile(spikesFilt_dir, 'int16') * 0.195
spikesFilt_rs = np.reshape(spikesFilt, [41,7,len(spikeTimes_37)], 'F')  
ch_spikePCA_all = calculate_global_PCA(spikeSites_37, spikesFilt_rs, global_neigh_37, total_num_channels)
meanWfs_37, meanPCs_37, varPCs_37 = calculate_meanWf_meanPC(days_37, num_units_37, spikesFilt_rs, ch_spikePCA_all, unit_lifetime_matrix_37, spikesByCluster_37, spike_sessions_37, clusterSites_37, global_neigh_37)

spikesFilt_dir = 'rTBY37_RSC_combined_filt.jrc'
spikesFilt = np.fromfile(spikesFilt_dir, 'int16') * 0.195
spikesFilt_rs = np.reshape(spikesFilt, [41,11,len(spikeTimes_37_rsc)], 'F')  
ch_spikePCA_all = calculate_global_PCA(spikeSites_37_rsc, spikesFilt_rs, global_neigh_37_rsc, total_num_channels)
meanWfs_37_rsc, meanPCs_37_rsc, varPCs_37_rsc = calculate_meanWf_meanPC(days_37, num_units_37_rsc, spikesFilt_rs, ch_spikePCA_all, unit_lifetime_matrix_37_rsc, spikesByCluster_37_rsc, spike_sessions_37_rsc, clusterSites_37_rsc, global_neigh_37_rsc)

#bundle identity of each unit
clusterShanks_35 = np.floor(clusterSites_35 / 64) 
clusterSites_37 = np.concatenate((clusterSites_37, clusterSites_37_rsc))
clusterShanks_37 = np.floor(clusterSites_37 / 64)

#%%
#Generating the correlation and SMD stats for each tracked unit
corr_stats = {}
smd_stats = {}
corrs_unit_all = np.zeros(1)
corrs_neigh_nonunit_all = np.zeros(1)
corrs_nonunit_all = np.zeros(1)
smds_unit_all = np.zeros(1)
smds_neigh_nonunit_all = np.zeros(1)
smds_nonunit_all = np.zeros(1)

#First going over the bundles in Rat 35
for shank in [0,1,3]:
    shank_spikes = np.where(clusterShanks_35 == shank)[0] #IDs of spikes in the bundle
    
    #Gathering the mean waveforms and mean/variance of PCs belonging to this bundle
    meanWf_shank = meanWfs_35[:,:,(shank*64):((shank+1)*64),:]
    meanWf_shank = meanWf_shank[shank_spikes[0]:(shank_spikes[-1]+1),:,:,:]
    
    meanPC_shank = meanPCs_35[:,:,(shank*64*num_PCs):((shank+1)*64*num_PCs)]
    meanPC_shank = meanPC_shank[shank_spikes[0]:(shank_spikes[-1]+1),:,:]
    
    varPC_shank = varPCs_35[:,:,(shank*64*num_PCs):((shank+1)*64*num_PCs)] 
    varPC_shank = varPC_shank[shank_spikes[0]:(shank_spikes[-1]+1),:,:]
    
    #Identifying the single units, their sites and lifetimes in the bundle
    shank_single_units = single_units35[clusterShanks_35 == shank]
    shank_single_unit_idx = np.where(shank_single_units)[0]
    num_shank_single_units = len(shank_single_unit_idx)
    clusterSites_35_shank_single_units = clusterSites_35[shank_spikes][shank_single_units]
    single_unit_lifetime_matrix_shank = unit_lifetime_matrix_35[shank_spikes[shank_single_units]]
    
    #Taking the mean waveforms of only the single units and reshaping the matrix to perform correlation
    meanWfs_by_shank_rs = meanWf_shank[shank_single_units,:,:,:]
    meanWfs_by_shank_rs = np.reshape(meanWfs_by_shank_rs,[num_shank_single_units,len(days_35),64*41])
    meanWfs_by_shank_rs = np.reshape(meanWfs_by_shank_rs,[num_shank_single_units*len(days_35), 64*41])
    
    #Calculating the correlation coefficients of mean waveforms of single units in the bundle 
    corrs = np.corrcoef(meanWfs_by_shank_rs)
    corrs = np.nan_to_num(corrs)
    
    #Taking the mean/variance PCs of only the single units and reshaping the matrix to perform correlation
    meanPCs_by_shank_rs = meanPC_shank[shank_single_units,:,:]
    varPCs_by_shank_rs = varPC_shank[shank_single_units,:,:]
    meanPCs_by_shank_rs = np.reshape(meanPCs_by_shank_rs,[num_shank_single_units*len(days_35),64*num_PCs])
    varPCs_by_shank_rs = np.reshape(varPCs_by_shank_rs,[num_shank_single_units*len(days_35), 64*num_PCs])
    
    #Calculating the standard mean differences of PCs of single units in the bundle
    smds = np.zeros((num_shank_single_units*len(days_35),num_shank_single_units*len(days_35)))
    for i in tqdm(range(num_shank_single_units*len(days_35))):
        for j in range(num_shank_single_units*len(days_35)):
            m1 = meanPCs_by_shank_rs[i]
            m2 = meanPCs_by_shank_rs[j]
            s1 = varPCs_by_shank_rs[i]
            s2 = varPCs_by_shank_rs[j]
            smds[i][j] = np.sqrt(np.sum(np.square(np.nan_to_num((m1-m2)/np.sqrt((s1+s2)/2)))))
    
    #Calculating the statistics of correlation and SMD for each unit
    corrs_stats_shank = np.zeros((num_shank_single_units,9))
    smd_stats_shank = np.zeros((num_shank_single_units,9))

    for unit in range(num_shank_single_units):
        unit_id_in_shank = shank_single_unit_idx[unit]
        #Correlations and SMDs of clusters belonging to the same neuron across sessions  
        corrs_unit = corrs[unit*len(days_35):((unit+1)*len(days_35)),unit*len(days_35):((unit+1)*len(days_35))]
        corrs_unit = np.triu(corrs_unit,1)
                
        smds_unit = smds[unit*len(days_35):((unit+1)*len(days_35)),unit*len(days_35):((unit+1)*len(days_35))]
        smds_unit = np.triu(smds_unit,1)
        try: 
            adjacent_pairs = np.lib.stride_tricks.sliding_window_view(np.nonzero(single_unit_lifetime_matrix_shank[unit])[0], 2)
            corrs_unit = corrs_unit[adjacent_pairs[:,0],adjacent_pairs[:,1]]
            corrs_unit_all = np.concatenate((corrs_unit_all, corrs_unit))
            smds_unit = smds_unit[adjacent_pairs[:,0],adjacent_pairs[:,1]]
            smds_unit_all = np.concatenate((smds_unit_all, smds_unit))
        except ValueError:
            corrs_unit = []
            smds_unit = []
        
        #Correlations and SMDs of clusters belonging to different neuron than the neuron of interest in the neighboring channels
        #triu refers to the units whose ids are larger than the unit of interest, so that in the grand population we don't repeat the same stats twice
        unit_site = int(clusterSites_35_shank_single_units[unit])
        unit_neigh = global_neigh_35[unit_site]-1
        neigh_nonunit = []
        for i in range(len(unit_neigh)):
            if neigh_nonunit == []:
                neigh_nonunit = np.where(clusterSites_35_shank_single_units == unit_neigh[i])[0]
            else:
                neigh_nonunit = np.append(neigh_nonunit, np.where(clusterSites_35_shank_single_units == unit_neigh[i])[0])
        neigh_nonunit = np.sort(np.delete(neigh_nonunit,np.where(neigh_nonunit==unit)))
        
        neigh_nonunit_idx = []
        neigh_nonunit_idx_triu = []
        for i in neigh_nonunit:
            if neigh_nonunit_idx == []:
                neigh_nonunit_idx = np.arange(i*len(days_35),(i+1)*len(days_35))
            else:
                neigh_nonunit_idx = np.r_[neigh_nonunit_idx, i*len(days_35):(i+1)*len(days_35)]
            if i > unit:
                if neigh_nonunit_idx_triu == []:
                    neigh_nonunit_idx_triu = np.arange(i*len(days_35),(i+1)*len(days_35))
                else: 
                    neigh_nonunit_idx_triu = np.r_[neigh_nonunit_idx_triu, i*len(days_35):(i+1)*len(days_35)]
        
        corrs_neigh_nonunit = corrs[unit*len(days_35):((unit+1)*len(days_35)), neigh_nonunit_idx]        
        corrs_neigh_nonunit_triu = corrs[unit*len(days_35):((unit+1)*len(days_35)), neigh_nonunit_idx_triu]
        
        smds_neigh_nonunit = smds[unit*len(days_35):((unit+1)*len(days_35)), neigh_nonunit_idx]        
        smds_neigh_nonunit_triu = smds[unit*len(days_35):((unit+1)*len(days_35)), neigh_nonunit_idx_triu]
        
        neigh_nonunit_nonzero = np.nonzero(corrs_neigh_nonunit)
        neigh_nonunit_triu_nonzero = np.nonzero(corrs_neigh_nonunit_triu)
        
        corrs_neigh_nonunit = corrs_neigh_nonunit[neigh_nonunit_nonzero]
        corrs_neigh_nonunit_triu = corrs_neigh_nonunit_triu[neigh_nonunit_triu_nonzero]
        corrs_neigh_nonunit_all = np.concatenate((corrs_neigh_nonunit_all, corrs_neigh_nonunit_triu))
        
        smds_neigh_nonunit = smds_neigh_nonunit[neigh_nonunit_nonzero]
        smds_neigh_nonunit_triu = smds_neigh_nonunit_triu[neigh_nonunit_triu_nonzero]
        smds_neigh_nonunit_all = np.concatenate((smds_neigh_nonunit_all, smds_neigh_nonunit_triu))
        
        #Correlations and SMDs of clusters belonging to different neurons than the neuron of interest
        #triu refers to the units whose ids are larger than the unit of interest, so that in the grand population we don't repeat the same stats twice
        corrs_nonunit = corrs[unit*len(days_35):((unit+1)*len(days_35)),np.r_[0:unit*len(days_35),(unit+1)*len(days_35):num_shank_single_units*len(days_35)]]
        corrs_nonunit_triu = corrs[unit*len(days_35):((unit+1)*len(days_35)),(unit+1)*len(days_35):num_shank_single_units*len(days_35)]        
        
        smds_nonunit = smds[unit*len(days_35):((unit+1)*len(days_35)),np.r_[0:unit*len(days_35),(unit+1)*len(days_35):num_shank_single_units*len(days_35)]]
        smds_nonunit_triu = smds[unit*len(days_35):((unit+1)*len(days_35)),(unit+1)*len(days_35):num_shank_single_units*len(days_35)]             
        
        nonunit_nonzero = np.nonzero(corrs_nonunit)
        nonunit_triu_nonzero = np.nonzero(corrs_nonunit_triu)
        
        corrs_nonunit = corrs_nonunit[nonunit_nonzero]
        corrs_nonunit_triu = corrs_nonunit_triu[nonunit_triu_nonzero]
        corrs_nonunit_all = np.concatenate((corrs_nonunit_all, corrs_nonunit_triu))
        
        smds_nonunit = smds_nonunit[nonunit_nonzero]
        smds_nonunit_triu = smds_nonunit_triu[nonunit_triu_nonzero]
        smds_nonunit_all = np.concatenate((smds_nonunit_all, smds_nonunit_triu))
        
        #Calculating all the corr and SMD stats
        corrs_stats_shank[unit] = [np.mean(corrs_unit), np.std(corrs_unit), len(corrs_unit), np.mean(corrs_nonunit), np.std(corrs_nonunit), len(corrs_nonunit), np.mean(corrs_neigh_nonunit), np.std(corrs_neigh_nonunit), len(corrs_neigh_nonunit)]
        smd_stats_shank[unit] = [np.mean(smds_unit), np.std(smds_unit), len(smds_unit), np.mean(smds_nonunit), np.std(smds_nonunit), len(smds_nonunit), np.mean(smds_neigh_nonunit), np.std(smds_neigh_nonunit), len(smds_neigh_nonunit)]
    corr_stats[shank] = corrs_stats_shank
    smd_stats[shank] = smd_stats_shank

#Now doing the same for Rat 37 
unit_lifetime_matrix_37 = np.vstack((unit_lifetime_matrix_37, unit_lifetime_matrix_37_rsc))

for shank in range(4):
    shank_spikes = np.where(clusterShanks_37 == shank)[0]
    if shank == 2:
        meanWf_shank = meanWfs_37_rsc[:,:,(shank*64):((shank+1)*64),:]
        meanPC_shank = meanPCs_37_rsc[:,:,(shank*64*num_PCs):((shank+1)*64*num_PCs)]
        varPC_shank = varPCs_37_rsc[:,:,(shank*64*num_PCs):((shank+1)*64*num_PCs)]     
        global_neigh = global_neigh_37_rsc
        
    else:
        meanWf_shank = meanWfs_37[:,:,(shank*64):((shank+1)*64),:]
        meanPC_shank = meanPCs_37[:,:,(shank*64*num_PCs):((shank+1)*64*num_PCs)]
        varPC_shank = varPCs_37[:,:,(shank*64*num_PCs):((shank+1)*64*num_PCs)]
        
        meanWf_shank = meanWf_shank[shank_spikes[0]:(shank_spikes[-1]+1),:,:,:]
        meanPC_shank = meanPC_shank[shank_spikes[0]:(shank_spikes[-1]+1),:,:]
        varPC_shank = varPC_shank[shank_spikes[0]:(shank_spikes[-1]+1),:,:]
        global_neigh = global_neigh_37

    shank_single_units = single_units37[clusterShanks_37 == shank]
    shank_single_unit_idx = np.where(shank_single_units)[0]
    num_shank_single_units = len(shank_single_unit_idx)
    
    clusterSites_37_shank_single_units = clusterSites_37[shank_spikes][shank_single_units]
    single_unit_lifetime_matrix_shank = unit_lifetime_matrix_37[shank_spikes[shank_single_units]]
    
    meanWfs_by_shank_rs = meanWf_shank[shank_single_units,:,:,:]
    meanWfs_by_shank_rs = np.reshape(meanWfs_by_shank_rs,[num_shank_single_units,len(days_37),64*41])
    meanWfs_by_shank_rs = np.reshape(meanWfs_by_shank_rs,[num_shank_single_units*len(days_37), 64*41])
    corrs = np.corrcoef(meanWfs_by_shank_rs)
    corrs = np.nan_to_num(corrs)
    corrs_stats_shank = np.zeros((num_shank_single_units,9))
    
    meanPCs_by_shank_rs = meanPC_shank[shank_single_units,:,:]
    varPCs_by_shank_rs = varPC_shank[shank_single_units,:,:]
    meanPCs_by_shank_rs = np.reshape(meanPCs_by_shank_rs,[num_shank_single_units*len(days_37),64*num_PCs])
    varPCs_by_shank_rs = np.reshape(varPCs_by_shank_rs,[num_shank_single_units*len(days_37), 64*num_PCs])
    smds = np.zeros((num_shank_single_units*len(days_37),num_shank_single_units*len(days_37)))
    smd_stats_shank = np.zeros((num_shank_single_units,9))
    
    for i in tqdm(range(num_shank_single_units*len(days_37))):
        for j in range(num_shank_single_units*len(days_37)):
            m1 = meanPCs_by_shank_rs[i]
            m2 = meanPCs_by_shank_rs[j]
            s1 = varPCs_by_shank_rs[i]
            s2 = varPCs_by_shank_rs[j]
            smds[i][j] = np.sqrt(np.sum(np.square(np.nan_to_num((m1-m2)/np.sqrt((s1+s2)/2)))))
              
    for unit in range(num_shank_single_units):
        unit_id_in_shank = shank_single_unit_idx[unit]
        corrs_unit = corrs[unit*len(days_37):((unit+1)*len(days_37)),unit*len(days_37):((unit+1)*len(days_37))]
        corrs_unit = np.triu(corrs_unit,1)
        
        smds_unit = smds[unit*len(days_37):((unit+1)*len(days_37)),unit*len(days_37):((unit+1)*len(days_37))]
        smds_unit = np.triu(smds_unit,1)
        try: 
            adjacent_pairs = np.lib.stride_tricks.sliding_window_view(np.nonzero(single_unit_lifetime_matrix_shank[unit])[0], 2)
            corrs_unit = corrs_unit[adjacent_pairs[:,0],adjacent_pairs[:,1]]
            corrs_unit_all = np.concatenate((corrs_unit_all, corrs_unit))
            smds_unit = smds_unit[adjacent_pairs[:,0],adjacent_pairs[:,1]]
            smds_unit_all = np.concatenate((smds_unit_all, smds_unit))
        except ValueError:
            corrs_unit = []
            smds_unit = []
        
        unit_site = int(clusterSites_37_shank_single_units[unit])
        unit_neigh = global_neigh[unit_site]-1
        neigh_nonunit = []
        for i in range(len(unit_neigh)):
            if neigh_nonunit == []:
                neigh_nonunit = np.where(clusterSites_37_shank_single_units == unit_neigh[i])[0]
            else:
                neigh_nonunit = np.append(neigh_nonunit, np.where(clusterSites_37_shank_single_units == unit_neigh[i])[0])
        neigh_nonunit = np.sort(np.delete(neigh_nonunit,np.where(neigh_nonunit==unit)))
        
        neigh_nonunit_idx = []
        neigh_nonunit_idx_triu = []
        for i in neigh_nonunit:
            if neigh_nonunit_idx == []:
                neigh_nonunit_idx = np.arange(i*len(days_37),(i+1)*len(days_37))
            else:
                neigh_nonunit_idx = np.r_[neigh_nonunit_idx, i*len(days_37):(i+1)*len(days_37)]
            if i > unit:
                if neigh_nonunit_idx_triu == []:
                    neigh_nonunit_idx_triu = np.arange(i*len(days_37),(i+1)*len(days_37))
                else: 
                    neigh_nonunit_idx_triu = np.r_[neigh_nonunit_idx_triu, i*len(days_37):(i+1)*len(days_37)]
        
        corrs_neigh_nonunit = corrs[unit*len(days_37):((unit+1)*len(days_37)), neigh_nonunit_idx]        
        corrs_neigh_nonunit_triu = corrs[unit*len(days_37):((unit+1)*len(days_37)), neigh_nonunit_idx_triu]
        
        smds_neigh_nonunit = smds[unit*len(days_37):((unit+1)*len(days_37)), neigh_nonunit_idx]        
        smds_neigh_nonunit_triu = smds[unit*len(days_37):((unit+1)*len(days_37)), neigh_nonunit_idx_triu]
        
        neigh_nonunit_nonzero = np.nonzero(corrs_neigh_nonunit)
        neigh_nonunit_triu_nonzero = np.nonzero(corrs_neigh_nonunit_triu)
        
        corrs_neigh_nonunit = corrs_neigh_nonunit[neigh_nonunit_nonzero]
        corrs_neigh_nonunit_triu = corrs_neigh_nonunit_triu[neigh_nonunit_triu_nonzero]
        corrs_neigh_nonunit_all = np.concatenate((corrs_neigh_nonunit_all, corrs_neigh_nonunit_triu))
        
        smds_neigh_nonunit = smds_neigh_nonunit[neigh_nonunit_nonzero]
        smds_neigh_nonunit_triu = smds_neigh_nonunit_triu[neigh_nonunit_triu_nonzero]
        smds_neigh_nonunit_all = np.concatenate((smds_neigh_nonunit_all, smds_neigh_nonunit_triu))
        
        corrs_nonunit = corrs[unit*len(days_37):((unit+1)*len(days_37)),np.r_[0:unit*len(days_37),(unit+1)*len(days_37):num_shank_single_units*len(days_37)]]
        corrs_nonunit_triu = corrs[unit*len(days_37):((unit+1)*len(days_37)),(unit+1)*len(days_37):num_shank_single_units*len(days_37)]
        
        smds_nonunit = smds[unit*len(days_37):((unit+1)*len(days_37)),np.r_[0:unit*len(days_37),(unit+1)*len(days_37):num_shank_single_units*len(days_37)]]
        smds_nonunit_triu = smds[unit*len(days_37):((unit+1)*len(days_37)),(unit+1)*len(days_37):num_shank_single_units*len(days_37)]
        
        nonunit_nonzero = np.nonzero(corrs_nonunit)
        nonunit_triu_nonzero = np.nonzero(corrs_nonunit_triu)
        
        corrs_nonunit = corrs_nonunit[nonunit_nonzero]
        corrs_nonunit_triu = corrs_nonunit_triu[nonunit_triu_nonzero]
        corrs_nonunit_all = np.concatenate((corrs_nonunit_all, corrs_nonunit_triu))
        
        smds_nonunit = smds_nonunit[nonunit_nonzero]
        smds_nonunit_triu = smds_nonunit_triu[nonunit_triu_nonzero]
        smds_nonunit_all = np.concatenate((smds_nonunit_all, smds_nonunit_triu))
        
        corrs_stats_shank[unit] = [np.mean(corrs_unit), np.std(corrs_unit), len(corrs_unit), np.mean(corrs_nonunit), np.std(corrs_nonunit), len(corrs_nonunit), np.mean(corrs_neigh_nonunit), np.std(corrs_neigh_nonunit), len(corrs_neigh_nonunit)]
        smd_stats_shank[unit] = [np.mean(smds_unit), np.std(smds_unit), len(smds_unit), np.mean(smds_nonunit), np.std(smds_nonunit), len(smds_nonunit), np.mean(smds_neigh_nonunit), np.std(smds_neigh_nonunit), len(smds_neigh_nonunit)]
    corr_stats[shank+4] = corrs_stats_shank
    smd_stats[shank+4] = smd_stats_shank 

#Deleting the zeros from the beginning of the arrays below   
corrs_unit_all = np.delete(corrs_unit_all,0)
corrs_neigh_nonunit_all = np.delete(corrs_neigh_nonunit_all,0)
corrs_nonunit_all = np.delete(corrs_nonunit_all, 0)
smds_unit_all = np.delete(smds_unit_all,0)
smds_neigh_nonunit_all = np.delete(smds_neigh_nonunit_all,0)
smds_nonunit_all = np.delete(smds_nonunit_all,0)

#comparing the correlation and SMD values of unit vs. neigh. non-unit for each single unit
corr_stats_all_units = np.concatenate((corr_stats[0],corr_stats[1],corr_stats[3],corr_stats[4],corr_stats[5],corr_stats[7],corr_stats[6]))
smd_stats_all_units = np.concatenate((smd_stats[0],smd_stats[1],smd_stats[3],smd_stats[4],smd_stats[5],smd_stats[7],smd_stats[6]))

corr_stats_p_values = np.zeros(len(corr_stats_all_units))
smd_stats_p_values = np.zeros(len(smd_stats_all_units))
for i in range(len(corr_stats_all_units)):
    corr_stats_p_values[i] = stats.ttest_ind_from_stats(corr_stats_all_units[i][0],corr_stats_all_units[i][1], corr_stats_all_units[i][2], corr_stats_all_units[i][6], corr_stats_all_units[i][7], corr_stats_all_units[i][8])[1]
    smd_stats_p_values[i] = stats.ttest_ind_from_stats(smd_stats_all_units[i][0], smd_stats_all_units[i][1], smd_stats_all_units[i][2], smd_stats_all_units[i][6], smd_stats_all_units[i][7], smd_stats_all_units[i][8])[1]
    
corr_stats_all_units = np.hstack((corr_stats_all_units, np.expand_dims(corr_stats_p_values, axis=1)))
smd_stats_all_units = np.hstack((smd_stats_all_units, np.expand_dims(smd_stats_p_values, axis=1)))
corr_stats_all_units = np.nan_to_num(corr_stats_all_units)
smd_stats_all_units = np.nan_to_num(smd_stats_all_units)

#Generating Supplementary Figure 5
plt.figure("Fig. S5A", figsize=[5,15])
plt.title('Correlations')
unit_violin = plt.violinplot(corrs_unit_all,positions=[1], showmeans=True, showextrema=False)
unit_violin['bodies'][0].set_facecolor('green')
neigh_nonunit_violin = plt.violinplot(corrs_neigh_nonunit_all,positions=[2], showmeans=True, showextrema=False)
neigh_nonunit_violin['bodies'][0].set_facecolor('orange')
nonunit_violin = plt.violinplot(corrs_nonunit_all,positions=[3], showmeans=True, showextrema=False)
nonunit_violin['bodies'][0].set_facecolor('red')
plt.ylim([-0.2,1.0])
plt.show()

plt.figure("Fig. S5B", figsize=[5,15])
plt.title('SMDs')
unit_violin = plt.violinplot(smds_unit_all,positions=[1], showmeans=True, showextrema=False)
unit_violin['bodies'][0].set_facecolor('green')
neigh_nonunit_violin = plt.violinplot(smds_neigh_nonunit_all,positions=[2], showmeans=True, showextrema=False)
neigh_nonunit_violin['bodies'][0].set_facecolor('orange')
nonunit_violin = plt.violinplot(smds_nonunit_all,positions=[3], showmeans=True, showextrema=False)
nonunit_violin['bodies'][0].set_facecolor('red')
plt.ylim([0,35])
plt.show()

#%%
#Calculating the SNRs across sessions 
#Location of the raw data 
r35_filedir = '/home/baran/data/rTBY35_combined/combined.dat'
r37_filedir = '/home/baran/data/rTBY37_combined/combined.dat'

#Arrays containing the root mean square of baseline activity for each session and each recording contact
rms_35 = np.zeros((len(days_35), total_num_channels))
rms_37 = np.zeros((len(days_37), total_num_channels))

#Designing the bandpass filter
b, a = signal.butter(4, [500,5000], 'bandpass', fs=sr)

#Iterating over each recording session, taking a 10 minute snippet from each session, filtering it and calculating the RMS of each channel of the bandpassed data
for i in tqdm(range(len(days_35))):
    data = np.fromfile(r35_filedir, 'int16', count=sr*60*total_num_channels, offset=int(time_boundaries_35[i]*sr*total_num_channels+10*60*sr*total_num_channels)) * 0.195
    data_rs = np.reshape(data, [total_num_channels, sr*60], 'F')
    data_bp = signal.filtfilt(b,a,data_rs,axis=1)
    rms_35[i] = np.sqrt(np.mean(np.square(data_bp), axis=1))

for i in tqdm(range(len(days_37))):
    data = np.fromfile(r37_filedir, 'int16', count=sr*60*total_num_channels, offset=int(time_boundaries_37[i]*sr*total_num_channels+10*60*sr*total_num_channels)) * 0.195
    data_rs = np.reshape(data, [total_num_channels, sr*60], 'F')
    data_bp = signal.filtfilt(b,a,data_rs,axis=1)
    rms_37[i] = np.sqrt(np.mean(np.square(data_bp), axis=1))  

#Calculating the amplitudes of each mean waveform in each session 
meanWfs_35_max = np.max(np.abs(meanWfs_35),axis=3)
meanWfs_37_max = np.max(np.abs(np.concatenate((meanWfs_37, meanWfs_37_rsc), axis=0)),axis=3)

#Calculating the unit SNRs in each session by dividing their max amplitude in their primary recording contact by the RMS in that contact 
unitSNRs_35 = np.zeros((meanWfs_35_max.shape[0], meanWfs_35_max.shape[1]))
unitSNRs_37 = np.zeros((meanWfs_37_max.shape[0], meanWfs_37_max.shape[1]))
for i in range(meanWfs_35_max.shape[0]):
    unitSNRs_35[i] = np.max(np.divide(meanWfs_35_max[i,:,:], rms_35), axis=1)
for i in range(meanWfs_37_max.shape[0]):
    unitSNRs_37[i] = np.max(np.divide(meanWfs_37_max[i,:,:], rms_37), axis=1)
single_unit_SNRs_35 = unitSNRs_35[single_units35]
single_unit_SNRs_37 = unitSNRs_37[single_units37]   

#calculating the mean SNR for each unit across sessions where that unit appeared
SNRs_all = np.zeros((unitSNRs_35.shape[0] + unitSNRs_37.shape[0]))
for unit in range(unitSNRs_35.shape[0]):
    unitSNRs = unitSNRs_35[unit]
    unitSNRs = unitSNRs[np.nonzero(unitSNRs)]
    SNRs_all[unit] = np.mean(unitSNRs)
for unit in range(unitSNRs_37.shape[0]):
    unitSNRs = unitSNRs_37[unit]
    unitSNRs = unitSNRs[np.nonzero(unitSNRs)]
    SNRs_all[unit+unitSNRs_35.shape[0]] = np.mean(unitSNRs)

#Generating supplementary Figure 4
plt.figure("Fig. S4A", figsize=[5,15])
plt.title('SNR')
single_violin = plt.violinplot(SNRs_all[single_units_all],positions=[1], showmeans=True, showextrema=False)
single_violin['bodies'][0].set_facecolor('green')
multi_violin = plt.violinplot(SNRs_all[multi_units_all],positions=[2], showmeans=True, showextrema=False)
multi_violin['bodies'][0].set_facecolor('red')
plt.ylim([0,75])
plt.axhline(y=np.min(SNRs_all[single_units_all]), color='k', linestyle='--')
plt.axhline(y=np.min(SNRs_all[multi_units_all]), color='k', linestyle='--')
plt.show()

plt.figure("Fig. S4B", figsize=[5,15])
plt.title('ISI_vio')
single_violin = plt.violinplot(isi_all[single_units_all],positions=[1], showmeans=True, showextrema=False)
single_violin['bodies'][0].set_facecolor('green')
multi_violin = plt.violinplot(isi_all[multi_units_all],positions=[2], showmeans=True, showextrema=False)
multi_violin['bodies'][0].set_facecolor('red')
plt.ylim([0,10])
plt.axhline(y=2, color='k', linestyle='--')
plt.show()

plt.figure("Fig. S4C", figsize=[5,15])
plt.title('dist_from_nn')
single_violin = plt.violinplot(nn_distances_all[single_units_all],positions=[1], showmeans=True, showextrema=False)
single_violin['bodies'][0].set_facecolor('green')
multi_violin = plt.violinplot(nn_distances_all[multi_units_all],positions=[2], showmeans=True, showextrema=False)
multi_violin['bodies'][0].set_facecolor('red')
plt.ylim([0,20])
plt.show()

#%%
#Tracking the single unit SNRs and unit yields across sessions 
#Weeks of each recording session 
weeks_35 = np.floor(days_35/7).astype('int16')
weeks_37 = np.floor(days_37/7).astype('int16')

single_unit_SNR_means = np.zeros(14)
single_unit_SNR_errs = np.zeros(14)
single_unit_SNR_weeks = {}

unit_counts = np.zeros(14)
unit_counts_hpc = np.zeros(14)
unit_counts_ctx = np.zeros(14)

#Identifying which units are HPC and which are cortex
hpc_end_35_all_units = np.where(clusterShanks_35 == 1)[0][-1]
unit_lifetime_matrix_35_hpc = unit_lifetime_matrix_35[0:hpc_end_35_all_units,:]
unit_lifetime_matrix_35_ctx = unit_lifetime_matrix_35[hpc_end_35_all_units+1:,:]

hpc_end_37_all_units = np.where(clusterShanks_37 == 1)[0][-1]
unit_lifetime_matrix_37_hpc = unit_lifetime_matrix_35[0:hpc_end_37_all_units,:]
unit_lifetime_matrix_37_ctx = unit_lifetime_matrix_35[hpc_end_37_all_units+1:,:]

#Iterating over the weeks to calculate the unit counts and mean single unit SNRs
for wk in range(14):
    #Taking all the single unit mean SNRs during the week 
    single_unit_SNRs_week = []
    for i in np.where(weeks_35 == wk)[0]:
        single_unit_SNRs_session = single_unit_SNRs_35[:,i]
        single_unit_SNRs_session = single_unit_SNRs_session[np.nonzero(single_unit_SNRs_session)]
        if single_unit_SNRs_week == []:
            single_unit_SNRs_week = single_unit_SNRs_session 
        else:
            single_unit_SNRs_week = np.concatenate((single_unit_SNRs_week, single_unit_SNRs_session))
    #Taking the number of units occurred during the week         
    if len(np.where(weeks_35==wk)[0]) == 1:
        unit_count_hpc = len(np.where(unit_lifetime_matrix_35_hpc[:,np.where(weeks_35==wk)[0]])[0])
        unit_count_ctx = len(np.where(unit_lifetime_matrix_35_ctx[:,np.where(weeks_35==wk)[0]])[0])
    elif np.where(weeks_35 == wk)[0].size == 0:
        unit_count_hpc = 0
        unit_count_ctx = 0
    else:
        unit_count_hpc = len(np.where(np.sum(unit_lifetime_matrix_35_hpc[:,np.where(weeks_35 == wk)[0]],axis=1))[0])
        unit_count_ctx = len(np.where(np.sum(unit_lifetime_matrix_35_ctx[:,np.where(weeks_35 == wk)[0]],axis=1))[0])
    
    #Repeating the same for rat 37
    for i in np.where(weeks_37 == wk)[0]:
        single_unit_SNRs_session = single_unit_SNRs_37[:,i]
        single_unit_SNRs_session = single_unit_SNRs_session[np.nonzero(single_unit_SNRs_session)]
        if single_unit_SNRs_week == []:
            single_unit_SNRs_week = single_unit_SNRs_session 
        else:
            single_unit_SNRs_week = np.concatenate((single_unit_SNRs_week, single_unit_SNRs_session))
    
    if len(np.where(weeks_37==wk)[0]) == 1:
        unit_count_hpc = unit_count_hpc + len(np.where(unit_lifetime_matrix_37_hpc[:,np.where(weeks_37 == wk)[0]])[0])
        unit_count_ctx = unit_count_ctx + len(np.where(unit_lifetime_matrix_37_ctx[:,np.where(weeks_37 == wk)[0]])[0])
    elif np.where(weeks_37 == wk)[0].size == 0:
        pass
    else:
        unit_count_hpc = unit_count_hpc + len(np.where(np.sum(unit_lifetime_matrix_37_hpc[:,np.where(weeks_37 == wk)[0]],axis=1))[0])
        unit_count_ctx = unit_count_ctx + len(np.where(np.sum(unit_lifetime_matrix_37_ctx[:,np.where(weeks_37 == wk)[0]],axis=1))[0])
    
    #Calculating the mean single unit SNRs and unit counts for the week 
    single_unit_SNR_means[wk] = np.mean(single_unit_SNRs_week)
    single_unit_SNR_errs[wk] = np.std(single_unit_SNRs_week) / np.sqrt(len(single_unit_SNRs_week))
    single_unit_SNR_weeks[wk] = single_unit_SNRs_week
    unit_counts_hpc[wk] = unit_count_hpc
    unit_counts_ctx[wk] = unit_count_ctx

#Calculating the unit yields by dividing the unit counts by the number of functional recording contacts    
"""nr of functional channels: 
rat 35 cortex = 61
rat 35 HPC = 125
rat 37 cortex = 117
rat 37 HPC = 124"""

unit_yields_hpc = np.zeros(14)
unit_yields_hpc[0] = unit_counts_hpc[0]/124
unit_yields_hpc[1:13] = unit_counts_hpc[1:13]/(124+125)
unit_yields_hpc[13] = unit_counts_hpc[13]/(125)

unit_yields_ctx = np.zeros(14)
unit_yields_ctx[0] = unit_counts_ctx[0]/117
unit_yields_ctx[1:13] = unit_counts_ctx[1:13]/(117+61)
unit_yields_ctx[13] = unit_counts_ctx[13]/(61)

#Generating Figure 4D/E
plt.figure("Fig. 4D-SNRs")
plt.plot(np.arange(1,13),single_unit_SNR_means[1:13])
plt.fill_between(np.arange(1,13), single_unit_SNR_means[1:13]-single_unit_SNR_errs[1:13], single_unit_SNR_means[1:13]+single_unit_SNR_errs[1:13],alpha=0.3)
plt.ylim([0,20])
plt.show()

plt.figure("Fig. 4D-yields")
plt.plot(np.arange(1,13),unit_yields_ctx[1:13])
plt.plot(np.arange(1,13),unit_yields_hpc[1:13])
plt.ylim([0,2.0])
plt.show()

single_unit_lifetimes_all = np.concatenate((single_unit_lifetimes35, single_unit_lifetimes37))
y,x = np.histogram(single_unit_lifetimes_all,bins=13,range=(0,91))
y = y/np.sum(y)
plt.figure("Fig. 4E")
plt.bar(np.arange(13),np.cumsum(np.flip(y)),width=0.9)
plt.show()
