#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 15:23:33 2023

@author: boramert
"""

import numpy as np
import matplotlib
#%matplotlib inline
import matplotlib.pyplot as plt
from itertools import compress
from PIL import Image
import json
import mne
import datetime

datapath = '/Users/boramert/Documents/NIRx/Data/2023-03-02/2023-03-02_003/2023-03-02_003.snirf'
raw_data = mne.io.read_raw_snirf(datapath)

jsonpath = '/Users/boramert/Documents/NIRx/Data/2023-03-02/2023-03-02_003/2023-03-02_003_description.json'
with open(jsonpath, 'r') as f:
  desc = json.load(f)
  
exp_date = raw_data.info['meas_date'].strftime('%m') + "-" + raw_data.info['meas_date'].strftime('%d')
  
name = desc['subject']

raw_data.annotations
raw_data.annotations.rename({
    '2' : 'Std',
    '3' : 'Odd',
    '4' : 'X',
    '5' : 'O (std)',
    '6' : 'O (odd)'
})

# raw_data.plot(n_channels=len(raw_data.ch_names),
#                     duration=500, show_scrollbars=True);

raw_od = mne.preprocessing.nirs.optical_density(raw_data);
# raw_od.plot(n_channels=len(raw_od.ch_names),
#              duration=500, show_scrollbars=False);

sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)
fig, ax = plt.subplots()
ax.hist(sci);
ax.set(xlabel='Scalp Coupling Index', ylabel='Count', xlim=[0, 1]);
plt.show()

raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.5))

#Check PSD for heart beat
raw_od.plot_psd()

#Motion Artifact Removal
corrected_tddr = mne.preprocessing.nirs.temporal_derivative_distribution_repair(raw_od);
#Before Removal
# raw_od.plot(n_channels=40, duration=400, show_scrollbars=False, title = "Before Removal");
# #After Removal
# corrected_tddr.plot(n_channels=40, duration=400, show_scrollbars=False, title = "After Removal");

#Filter Data to Remove Heart Beat
raw_od = raw_od.filter(0.01, 0.1, h_trans_bandwidth=0.08,
                             l_trans_bandwidth=0.001 , method='iir')

# raw_od.plot_psd()

#Apply Beer-Lambert Law to Obtain HbO and HbR
raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)
#raw_haemo.plot(n_channels=len(raw_haemo.ch_names),
#                duration=500, show_scrollbars=True)

events, event_dict = mne.events_from_annotations(raw_haemo)
fig = mne.viz.plot_events(events, event_id=event_dict,
                          sfreq=raw_haemo.info['sfreq'])
fig.subplots_adjust(right=0.7);  # make room for the legend

#Create Epochs
tmin = -5
tmax = 22
picks_fnirs = mne.pick_types(raw_haemo.info, fnirs=True, exclude='bads')
epochs = mne.Epochs(raw_haemo, events, event_id=event_dict, picks = picks_fnirs,
                    tmin=tmin, tmax=tmax,
                    proj=True, baseline=(None, 0), preload=True, event_repeated = 'merge',
                    detrend=1, verbose=True)

epochs.drop_bad()

imgnames = ["" for x in range(len(epochs.ch_names)*2)]
num = 0
for i in epochs.ch_names:
    figs_s = epochs['Std'].plot_image(picks=i,
                             combine=None, vmin=-10,
                             vmax=10,
                             show = False,
                             title = 'Standard ' + i,
                             ts_args=dict(ylim=dict(hbo=[-5, 5],
                                                    hbr=[-5, 5])))

    figs_b = epochs['Odd'].plot_image(picks=i,
                             combine=None, vmin=-10,
                             vmax=10,
                             show = False,
                             title = 'Oddball ' + i,
                             ts_args=dict(ylim=dict(hbo=[-5, 5],
                                                    hbr=[-5, 5])))
    imgnames[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/Tekli/"  + exp_date + "_" + name + "_Std_" + i + ".png"
    figs_s[0].savefig(imgnames[num],dpi=400)
    num = num+1
    imgnames[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/Tekli/"  + exp_date + "_" + name + "_Odd_" + i + ".png"
    figs_b[0].savefig(imgnames[num],dpi=400)
    num = num+1

# Load images
images = [Image.open(f) for f in imgnames]

# Create figure and subplots
coll = 2
rowl = int(len(images)/coll)
fig, axs = plt.subplots(rowl, coll, figsize=(10,40))

num = 0
for row in range(rowl):
    for col in range(coll):
        axs[row,col].imshow(images[num])
        axs[row,col].axis('off')
        num = num+1
        
fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/" + exp_date + "_" + name + "_Oddball_Epoch_Avg_Per_Channel.png",dpi=400)

plt.tight_layout()
plt.show()


color_dict = dict(HbO='#AA3377', HbR='b')
styles_dict = dict(Odd=dict(linestyle='dashed'))

imgnames2 = ["" for x in range(int(len(epochs.ch_names)/2))]
num = 0
for i in epochs.ch_names:
    evoked_dict = {'Std/HbO': epochs['Std'].average(picks=epochs.ch_names[num]),
                   'Std/HbR': epochs['Std'].average(picks=epochs.ch_names[num+int(len(epochs.ch_names)/2)]),
                   'Odd/HbO': epochs['Odd'].average(picks=epochs.ch_names[num]),
                   'Odd/HbR': epochs['Odd'].average(picks=epochs.ch_names[num+int(len(epochs.ch_names)/2)])}

    # Rename channels until the encoding of frequency in ch_name is fixed
    for condition in evoked_dict:
        evoked_dict[condition].rename_channels(lambda x: x[:-4])
        
    figs = mne.viz.plot_compare_evokeds(evoked_dict, combine="mean", ci=0.95,show=False,
                                        title=i, colors=color_dict, styles=styles_dict)
    
    imgnames2[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/Tekli/"  + exp_date + "_" + name + i + ".png"
    figs[0].savefig(imgnames2[num],dpi=400)
    num = num+1
    if num==int(len(epochs.ch_names)/2):
        break

images = [Image.open(f) for f in imgnames2]
coll = 2
rowl = int(len(images)/coll)

fig, axs = plt.subplots(rowl, coll, figsize=(20,40))
num = 0
for row in range(rowl):
    for col in range(coll):
        axs[row,col].imshow(images[num])
        axs[row,col].axis('off')
        num = num+1

fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/" + exp_date + "_" + name + "_Oddball_HbO_HbR_Avg.png",dpi=400)

plt.tight_layout()
plt.show()

times = np.arange(-3, 22, 3.0)
topomap_args = dict(extrapolate='local')
fig = epochs['Std'].average(picks='hbo').plot_joint(
    times=times, topomap_args=topomap_args, title='Standard HbO');

fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/"  + exp_date + "_" + name + "_Joint_Std" + ".png", dpi = 400)


times = np.arange(-3, 22, 3.0)
topomap_args = dict(extrapolate='local')
fig = epochs['Odd'].average(picks='hbo').plot_joint(
    times=times, topomap_args=topomap_args, title='Oddball Hb');

fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/"  + exp_date + "_" + name + "_Joint_Odd" + ".png", dpi = 400)


