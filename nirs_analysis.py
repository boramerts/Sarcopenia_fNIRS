#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 18:40:29 2023

@author: boramert
"""

import nirs
import mne

# nirs.checkSubject("/Users/boramert/Documents/NIRx/Data/2023-04-14/2023-04-14_005")

rootdir = '/Users/boramert/Documents/NIRx/Data/'
subjects = nirs.LoadNirX(rootdir)
print("Data import complete...")
subjects = nirs.set_annotations(subjects,0)
print("Annotations complete...")

subjects = nirs.optical_density(subjects, False)
print("OD complete...")

subjects = nirs.beer_lambert(subjects, True, l_freq = 0.01, 
                            h_freq = 0.1, 
                            l_trans_bandwidth = 0.08, 
                            h_trans_bandwidth = 0.001)
print("Beer Lambert complete...")

# subjects = nirs.get_events(subjects)
# print("Events complete...")

figs = []
for num,haemo in enumerate(subjects[2]):
    experiment = haemo.experiment
    intensity = subjects[0][num]
    nirs.GLM(haemo,intensity,experiment)
    # figs.append(nirs.GLM(haemo,intensity,experiment))
    print("Processed subject ", num)


# for num, i in enumerate(subjects[3]):
#     filename = '/Users/boramert/Desktop/Yüksek Lisans/Matlab_Kod/' + i.info['subject_info']['his_id'] + '_' + i.experiment + '_epochs.fif' 
#     i.save(filename, overwrite = True)
#     mne.write_events('/Users/boramert/Desktop/Yüksek Lisans/Matlab_Kod/' + i.info['subject_info']['his_id'] + '_' + i.experiment + '_events.fif'
#                      , i.events
#                      , overwrite = True)

# for num, i in enumerate(subjects[3]):
#     if i.experiment == 'Oddball' and i.info['subject_info']['his_id'] == '009':
#         a = i

#nirs.plot_mean_resp(subjects)

#nirs.show_channel_locations(subjects[0][1])



# epochs = subjects[3]

# grip_epochs = []
# nback_epochs = []
# oddball_epochs = []

# for i in range(len(epochs)):
#     if epochs[i].experiment == 'Grip':
#         grip_epochs.append(epochs[i])
#     elif epochs[i].experiment == 'Nback':
#         nback_epochs.append(epochs[i])
#     elif epochs[i].experiment == 'Oddball':
#         oddball_epochs.append(epochs[i])
        
# channels = grip_epochs[0].ch_names
# grip_evokeds = defaultdict(list)

# for ch in range(40):
#     for i in range(len(grip_epochs)):
#         if grip_epochs[i].ch_names[ch] in channels: 
            