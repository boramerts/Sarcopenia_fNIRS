#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:24:24 2023

@author: boramert
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from itertools import compress
from PIL import Image
import json
import mne
import datetime
import os
import nirs_analysis

rootdir = '/Users/boramert/Documents/NIRx/Data/'

got_snirf = 0
got_json = 0

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
              nirs_analysis.gripAnalysis(snirf_path, json_path)
          elif desc['experiment'] == "Nback":
              nirs_analysis.nbackAnalysis(snirf_path, json_path)
          elif desc['experiment'] == "Oddball":
              nirs_analysis.oddballAnalysis(snirf_path, json_path)
    else:
        text = ".json or .snirf file is missing in " + datapath
        print(text)
        break