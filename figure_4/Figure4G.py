#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 15:19:43 2023

@author: baran
"""
from matplotlib import pyplot as plt
import numpy as np
from PIL import Image
from tqdm import tqdm
from scipy.signal import convolve 
from scipy import stats

Image.MAX_IMAGE_PIXELS = 500000000

def neigh_ct(imarray):
    neigh = np.zeros(imarray.shape)
    kernel = np.ones((3,3),dtype=np.int8)
    kernel[[0,0,1,2,2],[0,2,1,0,2]] = 0
    neigh = convolve(imarray, kernel, 'same')
    return neigh

def calculate_fluo_intensities(imarray, radius_images):
    fluo = {} 
    for r in range(1,len(radii)):
        non_zero_points = np.nonzero(imarray*radius_images[r])
        fluo[r] = imarray[non_zero_points[1],non_zero_points[0]]
    return fluo 

neuron_image_dir = '/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 4/material/AVG_TBY37_s1_n1_1776_1_12_neuron_4x4bins.tif'
gfap_image_dir = '/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 4/material/AVG_TBY37_s1_n1_1776_1_12_GFAP_4x4bins.tif'
iba_image_dir = '/home/baran/Dropbox (Yanik Lab)/Multiarea rat recordings/Figures/Figure 4/material/AVG_TBY37_s1_n1_1776_1_12_IBA_4x4bins.tif'

dr = 25
line_begin = [14189/4, 2922/4]
line_end = [16215/4, 7625/4]
line_points = 1000
radii = np.arange(0,500+dr,dr)
um_per_pix = 0.325 * 4

neuron_image = Image.open(neuron_image_dir)
neuron_imarray = np.asarray(neuron_image)

gfap_image = Image.open(gfap_image_dir)
gfap_imarray = np.asarray(gfap_image)

iba_image = Image.open(iba_image_dir)
iba_imarray = np.asarray(iba_image)

line_x = np.linspace(line_begin[0],line_end[0],line_points,dtype=np.single)
line_y = np.linspace(line_begin[1],line_end[1],line_points,dtype=np.single)

plt.figure()
plt.imshow(neuron_image)
plt.plot(line_x,line_y,color='white')
plt.show()

x = np.arange(neuron_image.size[0])
y = np.arange(neuron_image.size[1])
xx,yy = np.meshgrid(x,y)
xx = xx.astype('int16')
yy = yy.astype('int16')

a,b = np.polyfit(line_x, line_y, 1)

points_in_line = {}
points_in_radius = {}
for i in tqdm(range(len(line_x))):
    xi = line_x[i]
    y_perp = a*xi + b + (xi-x) / a
    y_perp = np.around(y_perp)
    points_in_line[i] = np.vstack((x,y_perp)) 
    dist_to_line = np.sqrt((x - line_x[i])**2 + (y_perp-line_y[i])**2) * um_per_pix
    for r in range(1,len(radii)):
        close_points = np.where(np.logical_and((dist_to_line > radii[r-1]),(dist_to_line < radii[r])))[0]
        close_points = points_in_line[i][:,close_points]
        if i == 0:
            points_in_radius[r] = close_points
        else: 
            points_in_radius[r] = np.hstack((points_in_radius[r], close_points))
            
radius_images = np.zeros((len(radii),neuron_imarray.shape[0],neuron_imarray.shape[1]),dtype='int8')
for r in tqdm(range(1,len(radii))):
    points_in_radius[r] = points_in_radius[r].astype('int16')
    radius_images[r][points_in_radius[r][1],points_in_radius[r][0]] = 1
    neigh = neigh_ct(radius_images[r])
    radius_images[r][np.logical_and((neigh > 2),(radius_images[r] == 0))] = 1
    
neuron_fluo = calculate_fluo_intensities(neuron_imarray, radius_images)
gfap_fluo = calculate_fluo_intensities(gfap_imarray, radius_images)
iba_fluo = calculate_fluo_intensities(iba_imarray, radius_images)

n_rand_sam = 40000

neuron_fluo_ds = np.zeros((len(radii)-1,n_rand_sam))
iba_fluo_ds = np.zeros((len(radii)-1,n_rand_sam))
gfap_fluo_ds = np.zeros((len(radii)-1,n_rand_sam))

for i in range(len(radii)-1):
    neuron_fluo_ds[i] = np.random.choice(neuron_fluo[i+1],n_rand_sam)
    iba_fluo_ds[i] = np.random.choice(iba_fluo[i+1],n_rand_sam)
    gfap_fluo_ds[i] = np.random.choice(gfap_fluo[i+1],n_rand_sam)

nbins = len(radii) - 1

neuron_mean_fluo = np.zeros(nbins)
gfap_mean_fluo = np.zeros(nbins)
iba_mean_fluo = np.zeros(nbins)

neuron_err = np.zeros(nbins)
gfap_err = np.zeros(nbins)
iba_err = np.zeros(nbins)

for i in range(0,nbins):
    neuron_mean_fluo[i] = np.mean(neuron_fluo_ds[i])
    gfap_mean_fluo[i] = np.mean(gfap_fluo_ds[i])
    iba_mean_fluo[i] = np.mean(iba_fluo_ds[i])
    
    neuron_err[i] = np.std(neuron_fluo_ds[i]) / np.sqrt(len(neuron_fluo_ds[i]))
    gfap_err[i] = np.std(gfap_fluo_ds[i]) / np.sqrt(len(gfap_fluo_ds[i]))
    iba_err[i] = np.std(iba_fluo_ds[i]) / np.sqrt(len(iba_fluo_ds[i]))

#%%
plt.figure(figsize=[10,15])
plt.errorbar(radii[1:], neuron_mean_fluo/neuron_mean_fluo[-1], color='magenta', yerr=neuron_err/neuron_mean_fluo[-1], capsize=10)
plt.errorbar(-radii[1:], neuron_mean_fluo/neuron_mean_fluo[-1], color='magenta', yerr=neuron_err/neuron_mean_fluo[-1], capsize=10)
plt.ylim([0,1.2])
plt.show()

plt.figure(figsize=[10,15])
plt.errorbar(radii[1:], gfap_mean_fluo/gfap_mean_fluo[-1], color='yellow', yerr=gfap_err/gfap_mean_fluo[-1],capsize=10)
plt.errorbar(-radii[1:], gfap_mean_fluo/gfap_mean_fluo[-1], color='yellow',yerr=gfap_err/gfap_mean_fluo[-1],capsize=10)
plt.ylim([0,1.2])
plt.show()

plt.figure(figsize=[10,15])
plt.errorbar(radii[1:], iba_mean_fluo/iba_mean_fluo[-1], color='cyan',yerr=iba_err/iba_mean_fluo[-1],capsize=10)
plt.errorbar(-radii[1:], iba_mean_fluo/iba_mean_fluo[-1], color='cyan',yerr=iba_err/iba_mean_fluo[-1],capsize=10)
plt.ylim([0,1.2])
plt.show()

