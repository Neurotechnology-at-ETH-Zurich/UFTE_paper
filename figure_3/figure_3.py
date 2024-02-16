import numpy as np
from scipy.signal import filtfilt, butter,decimate
from matplotlib import pyplot as plt 
import math
import mat73
from scipy.io import loadmat
import scipy

#Defining the probe site map, importing the chunk of data where the SWR from iHP and dHP are extracted

siteMap = [28,27,29,26,30,25,31,24,32,23,34,22,33,21,36,20,35,19,38,18,37,244,40,242,39,240,42,16,41,17,44,14,43,15,63,13,61,
              12,62,11,59,10,60,9,57,8,58,7,55,0,56,2,53,4,54,6,51,5,52,3,49,1,50,239,120,119,84,83,85,82,86,81,87,80,88,79,117,
              122,118,121,89,78,90,77,91,76,92,75,93,74,115,124,116,123,94,73,95,72,96,71,98,70,97,69,113,126,114,125,100,66,102,
              68,104,67,106,65,108,64,110,127,112,47,99,45,101,48,103,46,212,213,211,214,210,215,209,216,208,217,207,218,175,219,
              176,220,173,221,174,222,171,223,172,224,169,226,170,225,167,228,168,227,165,230,166,229,163,232,164,231,161,234,162,
              233,160,236,159,235,158,237,157,249,156,248,154,246,152,111,150,109,148,107,146,105,184,199,183,200,136,245,135,137,
              247,138,198,181,197,182,186,201,185,202,134,243,133,139,250,140,196,179,195,180,188,203,187,204,132,241,131,141,251,
              142,252,177,193,178,190,205,189,206,130,238,129,143,194,145,253,147,254,149,255,151,192,153,191,155,128,144]


sample_rate = 20000
num_channels = 256
chunk_size = 60
chunk = 39

folder = '/home/baran/Desktop/rTBY37_session8/'
amplifier = folder + 'amplifier.dat'
start_idx = chunk * chunk_size * sample_rate * num_channels
start_offset = int(start_idx*2)
end_idx  = (chunk+1) * chunk_size * sample_rate * num_channels
data = np.fromfile(amplifier, 'int16', end_idx-start_idx, offset=start_offset)
data = np.reshape(data,(num_channels,int(len(data)/num_channels)),'F') * 0.195
data = data[siteMap]

time = np.arange(0,60,1/20000)
b,a = butter(4,500/10000,'lowpass')
data_filtered = filtfilt(b,a,data,axis=1)

#Plotting a sample SWR from multiple channels of dHPC and saving the figure

start_time = 28.2
end_time = 28.3
start_idx = int(start_time * sample_rate)
end_idx = int(end_time * sample_rate) 

plt.figure(figsize=[2,10])
for i in np.arange(65,80):
    plt.plot(time[start_idx:end_idx],data_filtered[i][start_idx:end_idx]+i*1000, c='b') 
plt.savefig('rTBY37_session8_dHP_ripple_snapshot.svg',dpi=500,format='svg')

#Plotting a sample SWR from multiple channels of iHPC and saving the figure

start_time = 28.1
end_time = 28.2
start_idx = int(start_time * sample_rate)
end_idx = int(end_time * sample_rate) 

plt.figure(figsize=[2,15])
for i in np.arange(42,57):
    if i not in [50]: 
        plt.plot(time[start_idx:end_idx],data_filtered[i][start_idx:end_idx]+i*1500, c='r') 
        
plt.savefig('rTBY37_session8_iHP_ripple_snapshot.svg',dpi=500,format='svg')

#Loading the results of spike sorting

data_dict = mat73.loadmat(folder+'/spike_sorting_results_except_RSC/amplifier_res.mat')
spikesFilt_dir = folder+'/spike_sorting_results_except_RSC/amplifier_filt.jrc'
num_clusters = len(data_dict['clusterNotes'])

single_unit = np.zeros(num_clusters)
for i in range(num_clusters):
    if data_dict['clusterNotes'][i] == ['single']:
        single_unit[i] = 1
      
spikeTimes = data_dict['spikeTimes']
spikesByCluster = data_dict['spikesByCluster']
clusterSites = data_dict['clusterSites'].astype('int')
spikeSites = data_dict['spikeSites']
spikesFilt = np.fromfile(spikesFilt_dir, 'int16') * 0.195
spikesFilt_rs = np.reshape(spikesFilt, [41,7,len(spikeTimes)], 'F') 
spikeClusters = data_dict['spikeClusters']

#Plotting the spike waveforms of sample units from dHPC 
cluster = 14
cluster_spikes = spikesByCluster[cluster][0].astype(int) - 1
cluster_spike_wavs = spikesFilt_rs[:,0,cluster_spikes]
plt.figure(figsize=[5,5])
plt.plot(cluster_spike_wavs[:,0:500],'g',alpha=0.05)
plt.show()

cluster = 18
cluster_spikes = spikesByCluster[cluster][0].astype(int) - 1
cluster_spike_wavs = spikesFilt_rs[:,0,cluster_spikes]
plt.figure(figsize=[5,5])
plt.plot(cluster_spike_wavs[:,0:500],'g',alpha=0.05)
plt.show()

cluster = 19
cluster_spikes = spikesByCluster[cluster][0].astype(int) - 1
cluster_spike_wavs = spikesFilt_rs[:,0,cluster_spikes]
plt.figure(figsize=[5,5])
plt.plot(cluster_spike_wavs[:,0:500],'g',alpha=0.05)
plt.show()

#Plotting the spike waveforms of sample units from iHPC
cluster = 6
cluster_spikes = spikesByCluster[cluster][0].astype(int) - 1
cluster_spike_wavs = spikesFilt_rs[:,0,cluster_spikes]
plt.figure(figsize=[5,5])
plt.plot(cluster_spike_wavs[:,0:500],'b',alpha=0.05)
plt.show()

cluster = 3
cluster_spikes = spikesByCluster[cluster][0].astype(int) - 1
cluster_spike_wavs = spikesFilt_rs[:,0,cluster_spikes]
plt.figure(figsize=[5,5])
plt.plot(cluster_spike_wavs[:,0:500],'b',alpha=0.05)
plt.show()

cluster = 
cluster_spikes = spikesByCluster[cluster][0].astype(int) - 1
cluster_spike_wavs = spikesFilt_rs[:,0,cluster_spikes]
plt.figure(figsize=[5,5])
plt.plot(cluster_spike_wavs[:,0:500],'b',alpha=0.05)
plt.show()

#Generating the mPFC units subfigure 

#Reading the data from the chunk of rTBY35/Sesson 18 which contains the DOWN -> UP state transition#plt.plot(time,data[213,:]-2000)
chunk = 30
folder = '/home/baran/Desktop/rTBY35_session18/'
amplifier = folder + 'amplifier.dat'
start_idx = chunk * chunk_size * sample_rate * num_channels
start_offset = int(start_idx*2)
end_idx  = (chunk+1) * chunk_size * sample_rate * num_channels
data = np.fromfile(amplifier, 'int16', end_idx-start_idx, offset=start_offset)
data = np.reshape(data,(num_channels,int(len(data)/num_channels)),'F') * 0.195
data = data[siteMap]
time = np.arange(0,chunk_size,1/sample_rate)

#Generating the plot which shows multiple traces during this transition
plt.figure()
plt.plot(time,data[130,:]) #patch trace
plt.plot(time,data[128,:]-1000) #another extracellular trace from RSC
plt.plot(time,data[94,:]-2000) #extracellular trace from the cortex above dHP
plt.plot(time,data[212,:]-3000) #extracellular trace from mPFC
plt.xlim([17.5,18.2])
plt.show()
plt.savefig('rTBY35_session18_up-down-state.svg')

#Bandpass filtering the patch trace and doing a simple thresholding at 250 uV to detect spikes of the neuron. 
#Generating the plot with the spike waveforms of this neuron afterwards.
b,a = butter(4,[300/10000,5000/10000], 'bandpass')
data_filtered = filtfilt(b,a,data,axis=1)
peak_idx = scipy.signal.find_peaks(data_filtered[130,:],prominence=250)

num_spike_samples = 80
pre_peak_idx = 30
spike_wavs = np.zeros((len(peak_idx[0]),num_spike_samples))

for i in range(num_spike_samples):
    spike_wavs[:,i] = data_filtered[130,(peak_idx[0]-(pre_peak_idx-i))]

plt.figure()
plt.plot(np.transpose(spike_wavs),'b',alpha=0.3)
plt.show()
#plt.savefig('rTBY35_session18_site130_unit_wav.svg')

#code for generating panel B
folder = '/home/baran/Desktop/rTBY34_session5/'
amplifier = folder + 'amplifier.dat'
chunk = 58
start_idx = chunk * chunk_size * sample_rate * num_channels
start_offset = int(start_idx*2)
end_idx  = (chunk+1) * chunk_size * sample_rate * num_channels
data = np.fromfile(amplifier, 'int16', end_idx-start_idx, offset=start_offset)
data = np.reshape(data,(num_channels,int(len(data)/num_channels)),'F') * 0.195
data = data[siteMap]
time = np.arange(0,chunk_size,1/sample_rate)
b,a = butter(4, [5000/10000], 'lowpass')
data_filtered = filtfilt(b,a,data,axis=1)
data_dict = mat73.loadmat(folder+'amplifier_res.mat')
spikeTimes = data_dict['spikeTimes']-1
spikesByCluster = data_dict['spikesByCluster']
clusterSites = data_dict['clusterSites'].astype('int')-1
spikeSites = data_dict['spikeSites']-1
spikeClusters = data_dict['spikeClusters']
clusterNotes = data_dict['clusterNotes']

iHPC_chs = [36,37,38,39]
dHPC_chs = [70,72,74,76]
RSC_chs = [137,139,140,141]
mPFC_chs = [197,223,224,225]

time_begin = 58*60+38.4
time_end = 58*60+38.6

def plot_spikes(chs,color):
    units = np.where(np.isin(clusterSites, chs))[0]
    unit_sites = clusterSites[units]

    for i, unit in enumerate(units):
        unit_site = np.where(unit_sites[i] == chs)[0][0]
        cluster_spikeTimes = (spikeTimes[spikesByCluster[unit][0].astype('int32')])/20000
        unit_spikes = np.where(np.logical_and((cluster_spikeTimes>time_begin),(cluster_spikeTimes<time_end)))[0]
        for spike in unit_spikes:
            spike_time = np.arange((cluster_spikeTimes[spike]-58*60)-5e-4,(cluster_spikeTimes[spike]-58*60)+1.5e-3,1/20000)
            spike_idx = (spike_time*20000).astype('int')
            plt.plot(spike_time,data_filtered[chs[unit_site],spike_idx]+unit_site*100,color,linewidth=1)

#Plotting the traces
plt.figure(figsize=[15,5])
for i in range(4):
    plt.plot(time, data_filtered[iHPC_chs[i]]+i*100, 'r', linewidth=0.2)
plt.xlim(38.4,38.6)
plt.ylim([-900,750])
plt.show()

plt.figure(figsize=[15,5])
for i in range(4):
    plt.plot(time, data_filtered[dHPC_chs[i]]+i*100, 'b', linewidth=0.2)
plot_spikes(dHPC_chs,'b')        
plt.xlim(38.4,38.6)
plt.ylim([-1600,1000])
plt.show()

plt.figure(figsize=[15,5])
for i in range(4):
    plt.plot(time, data_filtered[RSC_chs[i]]+i*100, 'green', linewidth=0.2)
plot_spikes(RSC_chs,'g')        
plt.xlim(38.4,38.6)
plt.ylim([-250,750])
plt.show()

plt.figure(figsize=[15,5])
for i in range(4):   
    plt.plot(time, data_filtered[mPFC_chs[i]]+i*100, 'orange', linewidth=0.2)
plot_spikes(mPFC_chs,'orange')  
plt.xlim(38.4,38.6)
plt.ylim([-200,500])
plt.show()
