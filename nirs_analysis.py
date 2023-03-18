#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:47:05 2023

@author: boramert
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import compress
from PIL import Image
import json
import mne
import datetime
import os
from collections import defaultdict

def get_channels(rootdir):
    channels = set(['S1_D1 hbo', 'S1_D2 hbo', 'S2_D1 hbo', 'S2_D2 hbo', 
                'S2_D4 hbo', 'S3_D2 hbo', 'S3_D3 hbo', 'S4_D2 hbo', 
                'S4_D3 hbo', 'S4_D4 hbo', 'S5_D5 hbo', 'S5_D6 hbo', 
                'S6_D5 hbo', 'S6_D6 hbo', 'S6_D8 hbo', 'S7_D6 hbo', 
                'S7_D7 hbo', 'S8_D6 hbo', 'S8_D7 hbo', 'S8_D8 hbo', 
                'S1_D1 hbr', 'S1_D2 hbr', 'S2_D1 hbr', 'S2_D2 hbr', 
                'S2_D4 hbr', 'S3_D2 hbr', 'S3_D3 hbr', 'S4_D2 hbr', 
                'S4_D3 hbr', 'S4_D4 hbr', 'S5_D5 hbr', 'S5_D6 hbr', 
                'S6_D5 hbr', 'S6_D6 hbr', 'S6_D8 hbr', 'S7_D6 hbr', 
                'S7_D7 hbr', 'S8_D6 hbr', 'S8_D7 hbr', 'S8_D8 hbr'])
    ch_names = []
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
                  raw_haemo, epochs = set_grips(snirf_path, json_path)
                  ch_names.append(epochs.ch_names)
              elif desc['experiment'] == "Nback":
                  raw_haemo, epochs = set_nbacks(snirf_path, json_path)
                  ch_names.append(epochs.ch_names)
              elif desc['experiment'] == "Oddball":
                  raw_haemo, epochs = set_oddballs(snirf_path, json_path)
                  ch_names.append(epochs.ch_names)
              
        elif (got_json==1 and got_snirf==0) or (got_json==0 and got_snirf==1):
            text = ".json or .snirf file is missing in " + datapath
            print(text)
            break
    for i in range(len(ch_names)):
        channels = channels.intersection(ch_names[i])
    
    return channels

def gripAnalysis(snirf_path, json_path):
    raw_data = mne.io.read_raw_snirf(snirf_path)

    with open(json_path, 'r') as f:
      desc = json.load(f)
      
    exp_date = raw_data.info['meas_date'].strftime('%m') + "-" + raw_data.info['meas_date'].strftime('%d')
      
    name = desc['subject']

    raw_data.annotations

    raw_data.annotations.rename({
        '2' : 'Block',
        '3' : 'Rest'
    })

    raw_od = mne.preprocessing.nirs.optical_density(raw_data);

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

    #Filter Data to Remove Heart Beat
    raw_od = corrected_tddr.filter(0.01, 0.1, h_trans_bandwidth=0.08,
                                 l_trans_bandwidth=0.001 , method='iir')

    # raw_od.plot_psd()

    #Apply Beer-Lambert Law to Obtain HbO and HbR
    raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)

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
        figs_s = epochs['Block'].plot_image(picks=i,
                                 combine=None, vmin=-10,
                                 vmax=10,
                                 show = False,
                                 title = 'Block ' + i,
                                 ts_args=dict(ylim=dict(hbo=[-5, 5],
                                                        hbr=[-5, 5])))

        figs_b = epochs['Rest'].plot_image(picks=i,
                                 combine=None, vmin=-10,
                                 vmax=10,
                                 show = False,
                                 title = 'Rest ' + i,
                                 ts_args=dict(ylim=dict(hbo=[-5, 5],
                                                        hbr=[-5, 5])))
        imgnames[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/Tekli/"  + exp_date + "_" + name + "_S_" + i + ".png"
        figs_s[0].savefig(imgnames[num],dpi=400)
        num = num+1
        imgnames[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/Tekli/"  + exp_date + "_" + name + "_B_" + i + ".png"
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
            
    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/" + exp_date + "_" + name + "_Grip_Epoch_Avg_Per_Channel.png",dpi=400)

    plt.tight_layout()
    plt.show()


    color_dict = dict(HbO='#AA3377', HbR='b')
    styles_dict = dict(Rest=dict(linestyle='dashed'))

    imgnames2 = ["" for x in range(int(len(epochs.ch_names)/2))]
    num = 0
    for i in epochs.ch_names:
        evoked_dict = {'Block/HbO': epochs['Block'].average(picks= epochs.ch_names[num]),
                       'Block/HbR': epochs['Block'].average(picks= epochs.ch_names[num+int(len(epochs.ch_names)/2)]),
                       'Rest/HbO': epochs['Rest'].average(picks=epochs.ch_names[num]),
                       'Rest/HbR': epochs['Rest'].average(picks=epochs.ch_names[num+int(len(epochs.ch_names)/2)])}

        # Rename channels until the encoding of frequency in ch_name is fixed
        for condition in evoked_dict:
            evoked_dict[condition].rename_channels(lambda x: x[:-4])
            
        figs = mne.viz.plot_compare_evokeds(evoked_dict, combine="mean", ci=0.95,show=False,
                                            title=i, colors=color_dict, styles=styles_dict)
        
        imgnames2[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/Tekli/"  + exp_date + "_" + name + i + ".png"
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

    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/" + exp_date + "_" + name + "_Grip_HbO_HbR_Avg.png",dpi=400)

    plt.tight_layout()
    plt.show()

    times = np.arange(-3, 22, 3.0)
    topomap_args = dict(extrapolate='local')
    fig = epochs['Block'].average(picks='hbo').plot_joint(
        times=times, topomap_args=topomap_args, title='Block HbO');

    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/"  + exp_date + "_" + name + "_Joint_Block" + ".png", dpi = 400)


    times = np.arange(-3, 22, 3.0)
    topomap_args = dict(extrapolate='local')
    fig = epochs['Rest'].average(picks='hbo').plot_joint(
        times=times, topomap_args=topomap_args, title='Rest HbO');

    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Grip/"  + exp_date + "_" + name + "_Joint_Rest" + ".png", dpi = 400)

def nbackAnalysis(snirf_path, json_path):
    raw_data = mne.io.read_raw_snirf(snirf_path)

    with open(json_path, 'r') as f:
      desc = json.load(f)
      
    exp_date = raw_data.info['meas_date'].strftime('%m') + "-" + raw_data.info['meas_date'].strftime('%d')
      
    name = desc['subject']

    raw_data.annotations

    raw_data.annotations.rename({
        '2' : 'Oback',
        '3' : 'Nback',
        '4' : 'Button'
    })

    raw_od = mne.preprocessing.nirs.optical_density(raw_data);

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

    #Filter Data to Remove Heart Beat
    raw_od = corrected_tddr.filter(0.01, 0.1, h_trans_bandwidth=0.08,
                                 l_trans_bandwidth=0.001 , method='iir')

    # raw_od.plot_psd()

    #Apply Beer-Lambert Law to Obtain HbO and HbR
    raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)

    events, event_dict = mne.events_from_annotations(raw_haemo)
    fig = mne.viz.plot_events(events, event_id=event_dict,
                              sfreq=raw_haemo.info['sfreq'])
    fig.subplots_adjust(right=0.7);  # make room for the legend

    #Create Epochs
    tmin = -5
    tmax = 30
    picks_fnirs = mne.pick_types(raw_haemo.info, fnirs=True, exclude='bads')
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict, picks = picks_fnirs,
                        tmin=tmin, tmax=tmax,
                        proj=True, baseline=(None, 0), preload=True, event_repeated = 'merge',
                        detrend=1, verbose=True)

    epochs.drop_bad()

    imgnames = ["" for x in range(len(epochs.ch_names)*2)]
    num = 0
    for i in epochs.ch_names:
        figs_s = epochs['Oback'].plot_image(picks=i,
                                 combine=None, vmin=-10,
                                 vmax=10,
                                 show = False,
                                 title = 'Oback ' + i,
                                 ts_args=dict(ylim=dict(hbo=[-5, 5],
                                                        hbr=[-5, 5])))

        figs_b = epochs['Nback'].plot_image(picks=i,
                                 combine=None, vmin=-10,
                                 vmax=10,
                                 show = False,
                                 title = 'Nback ' + i,
                                 ts_args=dict(ylim=dict(hbo=[-5, 5],
                                                        hbr=[-5, 5])))
        imgnames[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/Tekli/"  + exp_date + "_" + name + "_Oback_" + i + ".png"
        figs_s[0].savefig(imgnames[num],dpi=400)
        num = num+1
        imgnames[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/Tekli/"  + exp_date + "_" + name + "_Nback_" + i + ".png"
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
            
    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/" + exp_date + "_" + name + "_Nback_Epoch_Avg_Per_Channel.png",dpi=400)

    plt.tight_layout()
    plt.show()


    color_dict = dict(HbO='#AA3377', HbR='b')
    styles_dict = dict(Nback=dict(linestyle='dashed'))

    imgnames2 = ["" for x in range(int(len(epochs.ch_names)/2))]
    num = 0
    for i in epochs.ch_names:
        evoked_dict = {'Oback/HbO': epochs['Oback'].average(picks=epochs.ch_names[num]),
                       'Oback/HbR': epochs['Oback'].average(picks=epochs.ch_names[num+int(len(epochs.ch_names)/2)]),
                       'Nback/HbO': epochs['Nback'].average(picks=epochs.ch_names[num]),
                       'Nback/HbR': epochs['Nback'].average(picks=epochs.ch_names[num+int(len(epochs.ch_names)/2)])}

        # Rename channels until the encoding of frequency in ch_name is fixed
        for condition in evoked_dict:
            evoked_dict[condition].rename_channels(lambda x: x[:-4])
            
        figs = mne.viz.plot_compare_evokeds(evoked_dict, combine="mean", ci=0.95,show=False,
                                            title=i, colors=color_dict, styles=styles_dict)
        
        imgnames2[num] = "/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/Tekli/"  + exp_date + "_" + name + i + ".png"
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

    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/" + exp_date + "_" + name + "_Nback_HbO_HbR_Avg.png",dpi=400)

    plt.tight_layout()
    plt.show()

    times = np.arange(-3, 22, 3.0)
    topomap_args = dict(extrapolate='local')
    fig = epochs['Oback'].average(picks='hbo').plot_joint(
        times=times, topomap_args=topomap_args, title='Nback HbO');

    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/"  + exp_date + "_" + name + "_Joint_Oback" + ".png", dpi = 400)


    times = np.arange(-3, 22, 3.0)
    topomap_args = dict(extrapolate='local')
    fig = epochs['Nback'].average(picks='hbo').plot_joint(
        times=times, topomap_args=topomap_args, title='Nback Hb');

    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/Nback/"  + exp_date + "_" + name + "_Joint_Nback" + ".png", dpi = 400)

def oddballAnalysis(snirf_path, json_path):
    raw_data = mne.io.read_raw_snirf(snirf_path)

    with open(json_path, 'r') as f:
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

    raw_od = mne.preprocessing.nirs.optical_density(raw_data);

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

    #Filter Data to Remove Heart Beat
    raw_od = corrected_tddr.filter(0.01, 0.1, h_trans_bandwidth=0.08,
                                 l_trans_bandwidth=0.001 , method='iir')

    # raw_od.plot_psd()

    #Apply Beer-Lambert Law to Obtain HbO and HbR
    raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)

    events, event_dict = mne.events_from_annotations(raw_haemo)
    fig = mne.viz.plot_events(events, event_id=event_dict,
                              sfreq=raw_haemo.info['sfreq'])
    fig.subplots_adjust(right=0.7);  # make room for the legend

    #Create Epochs
    tmin = -5
    tmax = 40
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
    
def set_grips(snirf_path, json_path):
    raw_data = mne.io.read_raw_snirf(snirf_path)

    with open(json_path, 'r') as f:
      desc = json.load(f)
      
    exp_date = raw_data.info['meas_date'].strftime('%m') + "-" + raw_data.info['meas_date'].strftime('%d')
      
    name = desc['subject']

    raw_data.annotations

    raw_data.annotations.rename({
        '2' : 'Block',
        '3' : 'Rest'
    })

    raw_od = mne.preprocessing.nirs.optical_density(raw_data);

    sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)

    raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.5))

    #Check PSD for heart beat
    raw_od.plot_psd()

    #Motion Artifact Removal
    corrected_tddr = mne.preprocessing.nirs.temporal_derivative_distribution_repair(raw_od);

    #Filter Data to Remove Heart Beat
    raw_od = corrected_tddr.filter(0.01, 0.1, h_trans_bandwidth=0.08,
                                 l_trans_bandwidth=0.001 , method='iir')

    # raw_od.plot_psd()

    #Apply Beer-Lambert Law to Obtain HbO and HbR
    raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)

    events, event_dict = mne.events_from_annotations(raw_haemo)

    #Create Epochs
    tmin = -5
    tmax = 22
    picks_fnirs = mne.pick_types(raw_haemo.info, fnirs=True, exclude='bads')
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict, picks = picks_fnirs,
                        tmin=tmin, tmax=tmax,
                        proj=True, baseline=(None, 0), preload=True, event_repeated = 'merge',
                        detrend=1, verbose=True)

    epochs.drop_bad()
    return raw_haemo, epochs

def set_nbacks(snirf_path, json_path):
    raw_data = mne.io.read_raw_snirf(snirf_path)

    with open(json_path, 'r') as f:
      desc = json.load(f)
      
    exp_date = raw_data.info['meas_date'].strftime('%m') + "-" + raw_data.info['meas_date'].strftime('%d')
      
    name = desc['subject']

    raw_data.annotations

    raw_data.annotations.rename({
        '2' : 'Oback',
        '3' : 'Nback',
        '4' : 'Button'
    })

    raw_od = mne.preprocessing.nirs.optical_density(raw_data);

    sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)

    raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.5))

    #Check PSD for heart beat
    raw_od.plot_psd()

    #Motion Artifact Removal
    corrected_tddr = mne.preprocessing.nirs.temporal_derivative_distribution_repair(raw_od);

    #Filter Data to Remove Heart Beat
    raw_od = corrected_tddr.filter(0.01, 0.1, h_trans_bandwidth=0.08,
                                 l_trans_bandwidth=0.001 , method='iir')

    # raw_od.plot_psd()

    #Apply Beer-Lambert Law to Obtain HbO and HbR
    raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)

    events, event_dict = mne.events_from_annotations(raw_haemo)

    #Create Epochs
    tmin = -5
    tmax = 30
    picks_fnirs = mne.pick_types(raw_haemo.info, fnirs=True, exclude='bads')
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict, picks = picks_fnirs,
                        tmin=tmin, tmax=tmax,
                        proj=True, baseline=(None, 0), preload=True, event_repeated = 'merge',
                        detrend=1, verbose=True)

    epochs.drop_bad()
    return raw_haemo, epochs
    
def set_oddballs(snirf_path, json_path):
    raw_data = mne.io.read_raw_snirf(snirf_path)

    with open(json_path, 'r') as f:
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

    raw_od = mne.preprocessing.nirs.optical_density(raw_data);

    sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)

    raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.5))

    #Check PSD for heart beat
    raw_od.plot_psd()

    #Motion Artifact Removal
    corrected_tddr = mne.preprocessing.nirs.temporal_derivative_distribution_repair(raw_od);

    #Filter Data to Remove Heart Beat
    raw_od = corrected_tddr.filter(0.01, 0.1, h_trans_bandwidth=0.08,
                                 l_trans_bandwidth=0.001 , method='iir')

    # raw_od.plot_psd()

    #Apply Beer-Lambert Law to Obtain HbO and HbR
    raw_haemo = mne.preprocessing.nirs.beer_lambert_law(raw_od, ppf=0.1)

    events, event_dict = mne.events_from_annotations(raw_haemo)

    #Create Epochs
    tmin = -5
    tmax = 40
    picks_fnirs = mne.pick_types(raw_haemo.info, fnirs=True, exclude='bads')
    epochs = mne.Epochs(raw_haemo, events, event_id=event_dict, picks = picks_fnirs,
                        tmin=tmin, tmax=tmax,
                        proj=True, baseline=(None, 0), preload=True, event_repeated = 'merge',
                        detrend=1, verbose=True)

    epochs.drop_bad()
    return raw_haemo, epochs
    
def plot_averages(evokeds, hbo, hbr):
    fig, axes = plt.subplots(nrows=1, ncols=len(evokeds), figsize=(17, 5))
    lims = dict(hbo=[-5, 6], hbr=[-5, 6])
    picks = [hbo,hbr]
    print(picks)
    for (pick, color) in zip([hbo, hbr], ['r', 'b']):
        for idx, evoked in enumerate(evokeds):
            print(idx)
            print(pick)
            mne.viz.plot_compare_evokeds({evoked: evokeds[evoked]}, combine='mean',
                                 picks=pick, axes=axes[idx], show=False,
                                 colors=[color], legend=False, ylim=lims, ci=0.95,
                                 show_sensors=idx == 2)
            axes[idx].set_title('{} {}'.format(evoked, hbo))
    axes[0].legend(["HbO (Mean)", "HbO (Std)","", "Hb (Mean)", "Hb (Std)"])
    return fig
    
    
    
    
    