import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats

impedances = pd.read_csv('256ch_device2_impedance_spectroscopy.csv')
phases = pd.read_csv('256ch_device_phase_spectroscopy.csv') 

non_broken_phases = phases.drop(labels = impedances.index[impedances['1000'] > 5e6].tolist())
non_broken_impedances = impedances.drop(labels = impedances.index[impedances['1000'] > 5e6].tolist())

good_phases = non_broken_phases.drop(labels = non_broken_impedances.index[stats.zscore(non_broken_impedances['1000'])>3].tolist())
good_impedances = non_broken_impedances.drop(labels = non_broken_impedances.index[stats.zscore(non_broken_impedances['1000'])>3].tolist())

freqs = np.array(impedances.keys()).astype('int')
good_impedances_means = np.mean(good_impedances)
good_impedances_std = np.std(good_impedances)
good_phases_means = np.mean(good_phases)
good_phases_std = np.std(good_phases)

plt.figure(figsize=[3.7,2.2])
ax = plt.gca()
ax.errorbar(freqs, good_impedances_means, good_impedances_std, marker='o', markersize=4, capsize=4, ls='none')
ax.set_yscale('log')
ax.set_xscale('log')
ax.grid(True, 'major')
ax.grid(True, 'minor', alpha=0.1)
ax.set_xlabel('Frequency (Hz)',fontdict={'fontsize':12})
ax.set_ylabel(r'Impedance ($\Omega$)',fontdict={'fontsize':12})
ax.set_xticks([1,10,100,1000,10000])
ax.set_yticks([1e4,1e5,1e6,1e7])
ax.set_yticklabels(['10k','100k','1M','10M'],fontdict={'fontsize':10})
ax.set_xticklabels(['1','10','100','1000','10000'],fontdict={'fontsize':10})
ax.set_xlim([1,1e4])
ax.set_ylim([1e4,5e7])
plt.show()

plt.figure(figsize=[3.7,2.2])
ax = plt.gca()
ax.errorbar(freqs, good_phases_means, good_phases_std, marker='o', markersize=4, capsize=4, ls='none')
ax.set_xscale('log')
ax.grid(True, 'major')
ax.grid(True, 'minor', alpha=0.1)
ax.set_xlabel('Frequency (Hz)',fontdict={'fontsize':12})
ax.set_ylabel('Phase (Degrees)',fontdict={'fontsize':12})
ax.set_xticks([1,10,100,1000,10000])
ax.set_xlim([1,1e4])
ax.set_ylim([90,-90])
plt.savefig('Phase_spectroscopy.svg', dpi=300)
plt.show()
