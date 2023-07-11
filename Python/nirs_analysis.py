import mne
import mne_nirs
from mne_nirs.experimental_design import make_first_level_design_matrix
from mne_nirs.statistics import run_glm
from mne_nirs.channels import picks_pair_to_idx
from nilearn.plotting import plot_design_matrix

import pandas as pd
import numpy as np
import json
import os
import matplotlib.pyplot as plt
from itertools import compress

rootdir = '/Users/boramert/Documents/NIRx/Data/'
poststim = 5
prestim = -2

#Load Data
def LoadNirX(rootdir):
    # 0: raw, 1 : OD, 2: raw_hb, 3: filt_hb, 4: Epochs, 5: ROI GLM
    subjects = [[], [], [], [], [], []]
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            datapath = subdir+"/"+file
            if datapath.find("_description.json") != -1:
                with open(datapath, 'r') as f:
                    desc = json.load(f)
        try:
            raw_data = mne.io.read_raw_nirx(subdir)
            raw_data.experiment = desc['experiment']  # type: ignore
            subjects[0].append(raw_data)
        except:
            continue
    subjects = set_annotations(subjects)
    return subjects

#Configure Annotations
def set_annotations(subjects):
    for num, i in enumerate(subjects[0]):
        if i.experiment == 'Grip':
            i.annotations.set_durations({
                '2.0': 30,
                '3.0': 20})
            i.annotations.rename({
                '2.0' : 'Block',
                '3.0' : 'Rest'
            })
            unwanted = np.nonzero(i.annotations.description == '1.0')
            i.annotations.delete(unwanted)
        elif i.experiment == 'Nback':
            i.annotations.set_durations({
                 '2.0': 27.5,
                 '3.0': 27.5
            })
            i.annotations.rename({
                '2.0' : 'Oback',
                '3.0' : 'Nback'
            })
            i = remove_extras(i,0)
            unwanted = np.nonzero(i.annotations.description == '1.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '4.0')
            i.annotations.delete(unwanted)
            i.annotations.set_durations(27)
        elif i.experiment == 'Oddball':
            i.annotations.set_durations({
                 '2.0': 37.4,
                 '3.0': 37.4
            })
            i.annotations.rename({
                '2.0' : 'Std',
                '3.0' : 'Odd'
            })
            i = remove_extras(i,0)
            unwanted = np.nonzero(i.annotations.description == '1.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '4.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '5.0')
            i.annotations.delete(unwanted)
            unwanted = np.nonzero(i.annotations.description == '6.0')
            i.annotations.delete(unwanted)
            i.annotations.set_durations(37)
    return subjects

#Remove Extra Triggers
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

subjects = LoadNirX(rootdir)

for num,i in enumerate(subjects[0]):
    #Get Optical Density
    raw_od = mne.preprocessing.nirs.optical_density(i)
    #Get SCI
    sci = mne.preprocessing.nirs.scalp_coupling_index(raw_od)
    #Plot SCI
    fig, ax = plt.subplots()
    ax.hist(sci)
    ax.set(xlabel='Scalp Coupling Index', ylabel='Count', xlim=[0, 1])
    ax.set_title(i.info['subject_info']['his_id'])
    fig.savefig("/Users/boramert/Desktop/Yüksek Lisans/Python_Kod/figs/SCI/" + i.info['subject_info']['his_id'] + ".png",dpi=400)
    #Mark Bads
    raw_od.info['bads'] = list(compress(raw_od.ch_names, sci < 0.5))
    #TDDR
    corrected_tddr = mne.preprocessing.nirs.temporal_derivative_distribution_repair(raw_od)
    #Beer Lambert
    raw_haemo = mne.preprocessing.nirs.beer_lambert_law(corrected_tddr)
    #Filter Data
    filt_haemo = raw_haemo.filter(0.01, 0.1,
                                  h_trans_bandwidth = 0.08,
                                  l_trans_bandwidth = 0.001,
                                  method = 'iir', verbose = False)
    #Get Events
    events, event_dict = mne.events_from_annotations(filt_haemo)
    picks_fnirs = mne.pick_types(filt_haemo.info, fnirs=True, exclude='bads')
    reject_criteria = dict(hbo=6e-7)
    #Get Epochs
    if filt_haemo.experiment == 'Grip':
            stim = 20
            conditions = ['Rest', 'Block']
    elif filt_haemo.experiment == 'Nback':
            stim = 27.5
            conditions = ['Oback', 'Nback']
    elif filt_haemo.experiment == 'Oddball':
            stim = 37.4
            conditions = ['Std', 'Odd']
    tmin = prestim
    tmax = stim+poststim
    epochs = mne.Epochs(filt_haemo, events, event_id = event_dict, picks = picks_fnirs,
                        reject = reject_criteria, tmin = tmin, tmax = tmax,
                        proj = True, baseline = (prestim, 0) , preload = True,
                        detrend = 1, verbose = False)
    #Plot Removed Epochs
    epochs.plot_drop_log(subject = epochs.info['subject_info']['his_id'])
    #______________________________GLM______________________________
    #Create Design Matrix
    design_matrix = make_first_level_design_matrix(filt_haemo,
                                                   drift_model = 'cosine',
                                                   high_pass = 0.005,
                                                   hrf_model = 'spm',
                                                   stim_dur = stim)
    #Plot Design Matrix
    fig, ax1 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
    fig = plot_design_matrix(design_matrix, ax=ax1)
    #Plot Expected Response
    fig, ax1 = plt.subplots(figsize=(10, 6), nrows=1, ncols=1)
    s = mne_nirs.experimental_design.create_boxcar(filt_haemo, stim_dur=stim) # type: ignore
    plt.plot(filt_haemo.times, s[:, 0])
    plt.plot(filt_haemo.times, s[:, 1])
    plt.plot(design_matrix[conditions[0]]) # type: ignore
    plt.plot(design_matrix[conditions[1]]) # type: ignore
    plt.xlim(0, 200)
    ttl = filt_haemo.info['subject_info']['his_id'] + ' ' + filt_haemo.experiment
    plt.title(ttl)
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.show()
    #Fit GLM
    glm_est = run_glm(filt_haemo.pick(picks='all', exclude='bads'), design_matrix)
    #Get 3D view
    # glm_est.copy().surface_projection(condition="Block", view="dorsal",
    #                                   chroma="hbo",
    #                                   subjects_dir=mne.datasets.sample.data_path() / "subjects")
    #Set ROI
    
    # Primary Motor Cortex (BA4)
    Primary_Motor_Cortex_L = [[1,1],[2,1]]
    Primary_Motor_Cortex_R = [[5,5],[6,5]]

    # Premotor and Supplementary Motor Cortex (BA6)
    Pre_Supplementary_Motor_Cortex_L = [[2,2],[2,4],[4,2],[4,4],[3,2]]
    Pre_Supplementary_Motor_Cortex_R = [[5,5],[6,5],[5,6],[6,6],[7,6]]

    # Dorselateral Prefrontal Cortex (BA46 & BA9)
    Dorsolateral_Prefrontal_Cortex_L = [[4,4],[4,3],[4,2],[2,4]]
    Dorsolateral_Prefrontal_Cortex_R = [[8,8],[8,7],[8,6],[6,8]]

    groups = dict(Primary_Motor_Cortex_L =picks_pair_to_idx(filt_haemo, Primary_Motor_Cortex_L, on_missing='ignore'),
                  Primary_Motor_Cortex_R =picks_pair_to_idx(filt_haemo, Primary_Motor_Cortex_R, on_missing='ignore'),
                  Pre_Supplementary_Motor_Cortex_L =picks_pair_to_idx(filt_haemo, Pre_Supplementary_Motor_Cortex_L, on_missing='ignore'),
                  Pre_Supplementary_Motor_Cortex_R =picks_pair_to_idx(filt_haemo, Pre_Supplementary_Motor_Cortex_R, on_missing='ignore'),
                  Dorsolateral_Prefrontal_Cortex_L =picks_pair_to_idx(filt_haemo, Dorsolateral_Prefrontal_Cortex_L, on_missing='ignore'),
                  Dorsolateral_Prefrontal_Cortex_R =picks_pair_to_idx(filt_haemo, Dorsolateral_Prefrontal_Cortex_R, on_missing='ignore'))
    #Get ROI GLM Data
    df = glm_est.to_dataframe_region_of_interest(groups, conditions) # type: ignore
    df.head()

    #Append to Subjects Array
    subjects[1].append(raw_od)
    subjects[2].append(raw_haemo)
    subjects[3].append(filt_haemo)
    subjects[4].append(epochs)
    subjects[5].append(df)

    #Export Data
    filename = '/Users/boramert/Desktop/Yüksek Lisans/Exports/GLM_Data/' + i.info['subject_info']['his_id'] + '_' + i.experiment + '_GLM.csv' 
    df.to_csv(filename, index = False)

    filename = '/Users/boramert/Desktop/Yüksek Lisans/Exports/fNIRS_Data/' + filt_haemo.info['subject_info']['his_id'] + '_' + filt_haemo.experiment + '_data.fif' 
    filt_haemo.save(filename, overwrite = True)
    mne.write_events('/Users/boramert/Desktop/Yüksek Lisans/Exports/fNIRS_Data/' + filt_haemo.info['subject_info']['his_id'] + '_' + filt_haemo.experiment + '_events.fif'
                     , epochs.events
                     , overwrite = True)
    
    filename = '/Users/boramert/Desktop/Yüksek Lisans/Exports/fNIRS_Data/' + filt_haemo.info['subject_info']['his_id'] + '_' + filt_haemo.experiment + '_epochs.fif'
    epochs.save(filename, overwrite = True)




