#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 14:48:29 2023

@author: boramert
"""

import mne
import json
import os
import nirs_analysis
from collections import defaultdict
import matplotlib.pyplot as plt

rootdir = '/Users/boramert/Documents/NIRx/Data/'

got_snirf = 0 
got_json = 0

grip_evokeds = defaultdict(list)
nback_evokeds = defaultdict(list)
oddball_evokeds = defaultdict(list)

channels = list(nirs_analysis.get_channels(rootdir))
channels.sort()

for subdir, dirs, files in os.walk(rootdir):
    got_snirf = 0
    got_json = 0
    for file in files:
        datapath = subdir + "/" + file
        if datapath.find(".snirf") != -1:
            snirf_path = datapath
            got_snirf = 1
        if datapath.find("_description.json") != -1:
            json_path = datapath
            got_json = 1
    if got_json==1 and got_snirf==1:
        with open(json_path, 'r') as f:
          desc = json.load(f)
          if desc['experiment'] == "Grip":
              raw_haemo, epochs = nirs_analysis.set_grips(snirf_path, json_path)
              for cidx, condition in enumerate(epochs.event_id):
                  grip_evokeds[condition].append(epochs[condition].average(picks=channels))
          elif desc['experiment'] == "Nback":
              raw_haemo, epochs = nirs_analysis.set_nbacks(snirf_path, json_path)
              for cidx, condition in enumerate(epochs.event_id):
                  if condition == 'Nback' or condition == 'Oback':
                      nback_evokeds[condition].append(epochs[condition].average(picks=channels))
          elif desc['experiment'] == "Oddball":
              raw_haemo, epochs = nirs_analysis.set_oddballs(snirf_path, json_path)
              for cidx, condition in enumerate(epochs.event_id):
                  if condition == 'Odd' or condition == 'Std':
                      oddball_evokeds[condition].append(epochs[condition].average(picks=channels))
    elif (got_json==1 and got_snirf==0) or (got_json==0 and got_snirf==1):
        text = ".json or .snirf file is missing in " + datapath
        print(text)
        break

num = 0
for n,i in enumerate(channels):
    if num == len(channels):
        break
    fig = nirs_analysis.plot_averages(grip_evokeds,channels[num], channels[num+1])
    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/Grip_Average_" + str(num) + ".png",dpi=400)

    fig = nirs_analysis.plot_averages(nback_evokeds,channels[num], channels[num+1])
    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/Nback_Average_" + str(num) + ".png",dpi=400)

    fig = nirs_analysis.plot_averages(oddball_evokeds,channels[num], channels[num+1])
    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Oddball/Oddball_Average_" + str(num) + ".png",dpi=400)
    num = num+2
    
#Kanalları düzelt
