#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:34:40 2023

@author: boramert
"""
import numpy as np
import matplotlib.pyplot as plt
from itertools import compress
import json
import mne
import os
import mne_nirs
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm
from nilearn.plotting import plot_design_matrix

reject = ["004", "005", "008", "010", "011"]

def LoadNirX(rootdir):
    subjects = [[],[],[],[]]
    
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            datapath = subdir+"/"+file
            if datapath.find("_description.json") != -1:
                with open(datapath, 'r') as f:
                    desc = json.load(f)
        try:
            raw_data = mne.io.read_raw_nirx(subdir)
            raw_data.experiment = desc['experiment'] # type: ignore
            subjects[0].append(raw_data)
        except:
            continue
    return subjects
        
def show_channel_locations(raw_data):
    subjects_dir = mne.datasets.sample.data_path() / 'subjects'

    brain = mne.viz.Brain(
        'fsaverage', subjects_dir=subjects_dir, background='w', cortex='0.5')

    brain.add_sensors(
        raw_data.info, trans='fsaverage',
        fnirs = ['channels', 'pairs', 'sources', 'detectors']) # type: ignore
    brain.show_view(azimuth=20, elevation=60, distance=400)
        
def set_annotations(subjects, a):
    for num, i in enumerate(subjects[0]):
        if i.experiment == 'Grip':
            i.annotations.set_durations({
                '2.0': 35,
                '3.0': 20})
            i.annotations.rename({
                '2.0' : 'Block',
                '3.0' : 'Rest'
            })
            unwanted = np.nonzero(i.annotations.description == '1.0')
            i.annotations.delete(unwanted)
        elif i.experiment == 'Nback':
            i.annotations.rename({
                '2.0' : 'Oback',
                '3.0' : 'Nback'
            })
            i = remove_extras(i,a)
            unwanted = np.nonzero(i.annotations.description == '1.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '4.0')
            i.annotations.delete(unwanted)
            i.annotations.set_durations(27.5)
        elif i.experiment == 'Oddball':
            i.annotations.rename({
                '2.0' : 'Std',
                '3.0' : 'Odd'
            })
            i = remove_extras(i,a)
            unwanted = np.nonzero(i.annotations.description == '1.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '4.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '5.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '6.0')
            i.annotations.delete(unwanted)
            i.annotations.set_durations(37.4)
    return subjects

def remove_extras(subject,x):
    if subject.experiment == 'Nback':
        unwanted = np.nonzero(subject.annotations.description == 'Oback')
        unwanted =np.delete(unwanted, [a for a in range(x,len(unwanted[0]),2)])
        subject.annotations.delete(unwanted)
        unwanted = [a for a in range(len(subject.annotations.description)-3,len(subject.annotations.description))]
        subject.annotations.delete(unwanted)
        unwanted = np.nonzero(subject.annotations.description == 'Nback')
        unwanted = np.delete(unwanted, [a for a in range(x,len(unwanted[0]),2)])
        subject.annotations.delete(unwanted)
    elif subject.experiment == 'Oddball':
        unwanted = np.nonzero(subject.annotations.description == 'Std')
        unwanted = np.delete(unwanted, [a for a in range(x,len(unwanted[0]),2)])
        subject.annotations.delete(unwanted)
        unwanted = [a for a in range(len(subject.annotations.description)-3,len(subject.annotations.description))]
        subject.annotations.delete(unwanted)
        unwanted = np.nonzero(subject.annotations.description == 'Odd')
        unwanted = np.delete(unwanted, [a for a in range(x,len(unwanted[0]),2)])
        subject.annotations.delete(unwanted)
    
    return subject
    
def optical_density(subjects, p):
    print("_____________________________________")
    print("Converting to Optical Density")
    print("_____________________________________")
    for num, i in enumerate(subjects[0]):
        raw_od = mne.preprocessing.nirs.optical_density(i)
        
        sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)
        
        if p:
            fig, ax = plt.subplots()
            ax.hist(sci)
            ax.set(xlabel='Scalp Coupling Index', ylabel='Count', xlim=[0, 1])
            ax.set_title(i.info['subject_info']['his_id'])
            fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/SCI/" + i.info['subject_info']['his_id'] + ".png",dpi=400)
        
        raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.5))
        
        subjects[1].append(raw_od)
    
    return subjects
        
def beer_lambert(subjects, repair, 
                  l_freq, h_freq, 
                  h_trans_bandwidth,
                  l_trans_bandwidth):
    
    print("_____________________________________")
    print("Converting to Hbo-Hb")
    print("_____________________________________")
    
    for num, i in enumerate(subjects[1]):
        corrected_tddr = mne.preprocessing.nirs.temporal_derivative_distribution_repair(i)
        
        raw_haemo = mne.preprocessing.nirs.beer_lambert_law(corrected_tddr, ppf=0.1)
        
        raw_haemo = raw_haemo.filter(l_freq, h_freq,
                                     h_trans_bandwidth=h_trans_bandwidth,
                                     l_trans_bandwidth=l_trans_bandwidth , method='iir')
        subjects[2].append(raw_haemo)
    return subjects

def get_events(subjects):
    
    print("_____________________________________")
    print("Getting Epochs")
    print("_____________________________________")
    for num, i in enumerate(subjects[2]):
        events, event_dict = mne.events_from_annotations(i)
        
    
        picks_fnirs = mne.pick_types(i.info, fnirs=True, exclude='bads')
        if i.experiment == 'Grip': 
            tmin, tmax = -2, 20
            experiment = 'Grip'
        elif i.experiment == 'Nback': 
            tmin,tmax = -2, 27
            experiment = 'Nback'
        elif i.experiment == 'Oddball': 
            tmin,tmax = -2, 37
            experiment = 'Oddball'
        
        reject_criteria = dict(hbo=40e-6)
        epochs = mne.Epochs(i, events, event_id=event_dict, picks = picks_fnirs, reject=reject_criteria,
                            tmin=tmin, tmax=tmax,
                            proj=True, baseline=(None,0), preload=True,
                            detrend=1, verbose=True)
        epochs.plot_drop_log(subject = epochs.info['subject_info']['his_id'])
        
        baseline = (None,0)
        
        epochs = epochs.apply_baseline(baseline)
        
        if i.experiment == 'Grip' and (len(epochs['Rest']) == 0 or len(epochs['Block']) == 0):
            epochs.status = 'reject' # type: ignore
        elif i.experiment == 'Nback' and (len(epochs['Nback']) == 0 or len(epochs['Oback']) == 0):
            epochs.status = 'reject' # type: ignore
        elif i.experiment == 'Oddball' and (len(epochs['Std']) == 0 or len(epochs['Odd']) == 0):
            epochs.status = 'reject' # type: ignore
        else:
            epochs.status = 'normal'
        
        epochs.experiment = experiment
        
        subjects[3].append(epochs)
    return subjects
        
def GLM(haemo,intensity,experiment):
    if experiment == 'Grip':
            stim_dur = 20
            d_matrix = 'Block'
    elif experiment == 'Nback':
            stim_dur = 27.5
            d_matrix = 'Nback'
    elif experiment == 'Oddball':
            stim_dur = 37.4
            d_matrix = 'Odd'
    design_matrix = make_first_level_design_matrix(haemo,
                                                    drift_model='cosine',
                                                    high_pass=0.005,  # Must be specified per experiment
                                                    hrf_model='spm',
                                                    stim_dur=stim_dur)
        
    # fig, ax1 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
    # fig = plot_design_matrix(design_matrix, ax=ax1)

    s = mne_nirs.experimental_design.create_boxcar(intensity, stim_dur=stim_dur)
    plt.plot(intensity.times, s[:, 1])
    plt.plot(design_matrix[d_matrix])
    plt.xlim(0, 200)
    plt.legend(["Stimulus", "Expected Response"])
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.show()

    # return fig

def checkSubject(folderpath):
    raw_data = mne.io.read_raw_nirx(folderpath)
    print(raw_data.annotations)
    raw_data.plot(n_channels=len(raw_data.ch_names), show_scrollbars=False)
    return
