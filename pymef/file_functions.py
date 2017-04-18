#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 14:54:39 2017

Lightweight python library for pymef.

Ing.,Mgr. (MSc.) Jan Cimbálník
Biomedical engineering
International Clinical Research Center
St. Anne's University Hospital in Brno
Czech Republic
&
Mayo systems electrophysiology lab
Mayo Clinic
200 1st St SW
Rochester, MN
United States
"""

from .mef_file import pymef3_file
import numpy as np

def uutc_for_sample(sample,channel_info):
    
    # Sampling freq
    fs = channel_info['section_2']['sampling_frequency']
    
    seg_list = list(channel_info['segments'])
    seg_list.sort()
    
    prev_sample = 0
    prev_time = channel_info['segments'][seg_list[0]]['indices'][0]['start_time']
    
    for segment in seg_list: 
        indices = channel_info['segments'][segment]['indices']
        for index in indices:
            index_start_time = index['start_time']
            index_start_sample = index['start_sample']  #TODO: !!!!Start sample might be 0 for each segment!!!! Check this
            
            if index_start_sample > sample:
                return int(prev_time + np.ceil(((sample - prev_sample) / fs)))
            
            prev_sample = index_start_sample
            prev_time = index_start_time
            
    print('Sample number out of file')
    
    return None

def sample_for_uutc(uutc,channel_md,return_discont_distance=False):
    
    # UUTC check
    if not(uutc_check(uutc,channel_md)):
        print('uUTC time out of file')
        return None
    
    # Sampling freq
    fs = channel_md['section_2']['sampling_frequency']
    
    seg_list = list(channel_md['segments'])
    seg_list.sort()
    
    prev_sample = 0
    prev_N_samples = channel_md['segments'][seg_list[0]]['indices'][0]['number_of_samples']
    prev_time = channel_md['segments'][seg_list[0]]['indices'][0]['start_time']
    
    for segment in seg_list: 
        indices = channel_md['segments'][segment]['indices']
        for index in indices:
            index_start_time = index['start_time']
            index_start_sample = index['start_sample']   #TODO: !!!!Start sample might be 0 for each segment!!!! Check this
            index_N_samples = index['number_of_samples']
            
            if index_start_time > uutc:
                # Check for discontinuity
                if ((index_start_time - prev_time) / 1e6)*fs > prev_N_samples:
                    print('uUTC time at discontinuity, returning the closest point in the future')
                    if return_discont_distance:
                        return int(index_start_sample),((index_start_time-uutc) / 1000000) * fs
                    else:
                        return int(index_start_sample)
                
                if return_discont_distance:
                    return int(prev_sample + ((uutc - prev_time) / 1000000) * fs),0
                else:
                    return int(prev_sample + ((uutc - prev_time) / 1000000) * fs)
            
            prev_sample = index_start_sample
            prev_time = index_start_time
            prev_N_samples = index_N_samples
            
    # Ran out of indices
    if return_discont_distance:
        return int(prev_sample + ((uutc - prev_time) / 1000000) * fs),0
    else:
        return None
    
def uutc_check(uUTC,channel_md):
    
    # Get the uUTC start
    chan_uutc_start = channel_md['channel_specific_metadata']['earliest_start_time']
    # Get the uUTC stop
    chan_uutc_stop = channel_md['channel_specific_metadata']['latest_end_time']
    
    if type(channel_md['section_3']['recording_time_offset']) == int:
        chan_uutc_start = -chan_uutc_start + channel_md['section_3']['recording_time_offset']
        chan_uutc_stop = -chan_uutc_stop + channel_md['section_3']['recording_time_offset']
        
    if (uUTC>=chan_uutc_start and uUTC<=chan_uutc_stop):
        return True
    else:
        return False
        
def get_discontinuities(channel_md):
    """
    Processes indices and returns discontinuities accross segments
    
    Parameters:
    -----------
    channel_md - channel metadata dictionart\n
    
    Returns:
    --------
    discont_index_list - list of lists discontinuity start indices for individual segments
    discont_uutc_list - list of lists discontinuity start uutc times for individual segments
    discont_sample_list - list of lists discontinuity start samples for individual segments
    """
    
    discont_index_list = []
    discont_uutc_list = []
    discont_sample_list = []
    
    # Sampling freq
    fs = channel_md['section_2']['sampling_frequency']
    
    # Find out if offset is used
    if type(channel_md['section_3']['recording_time_offset']) == int:
        rt_offset = channel_md['section_3']['recording_time_offset']
        offset_flag = True
    else:
        rt_offset = 0
        offset_flag = False
    
    seg_list = list(channel_md['segments'])
    seg_list.sort()
    
    prev_N_samples = channel_md['segments'][seg_list[0]]['indices'][0]['number_of_samples']
    prev_time = channel_md['segments'][seg_list[0]]['indices'][0]['start_time']
    if offset_flag:
        prev_time = -prev_time + rt_offset
    
    for si,segment in enumerate(seg_list): 
        indices = channel_md['segments'][segment]['indices']
        
        # First block in segment is a discontinuity by definition
        seg_discont_index_list = []
        seg_discont_index_list.append(0)
        seg_discont_uutc_list = []
        seg_discont_uutc_list.append(indices[0]['start_time'])
        seg_discont_sample_list = []
        seg_discont_sample_list.append(indices[0]['start_sample']) #TODO: !!!!Start sample might be 0 for each segment!!!! Check this
        
        for ii,index in enumerate(indices):
            index_start_time = index['start_time']
            if offset_flag:
                index_start_time = -index_start_time + rt_offset
            index_N_samples = index['number_of_samples']
            
            # Previous block time length
            prev_block_span = prev_N_samples * fs * 1e6
            
            # Check for discontinuity
            if index_start_time > prev_time+prev_block_span:
                seg_discont_index_list.append(ii)
                seg_discont_uutc_list.append(index['start_time'])
                seg_discont_sample_list.append(index['start_sample']) #TODO: !!!!Start sample might be 0 for each segment!!!! Check this
                
            prev_time = index_start_time
            prev_N_samples = index_N_samples
            
        discont_index_list.append(seg_discont_index_list)
        discont_uutc_list.append(seg_discont_uutc_list)
        discont_sample_list.append(seg_discont_sample_list)
        
    return discont_index_list,discont_uutc_list,discont_sample_list

def read_ts_data_sample(session_path,
                        password,channel_map,sample_map):
    """
    Reads desired channels in desired sample segment
    
    Parameters:
    -----------
    session_path - path to mef3 session (.mefd)\n
    password - mef3 data password\n
    channel_map - list of channels to be read\n
    sample_map - list of [start,stop] samples to be loaded that correspond\n
        to channel_map. if there is only one entry the same range is applied\n
        to all channels\n
    
    Returns:
    --------
    data - numpy array [channels,samples]\n
    """
    
    data_list = [] # Creating a list since sampling frequency is not garanteed to be the same in all channels
    
    if len(sample_map) == 1:
        sample_map = sample_map*len(channel_map)
        
    if len(sample_map) != len(channel_map):
        print('Length of sample map is not equivalet to the length of channel map')
        return
        
    
    for channel,sample_ss in zip(channel_map,sample_map):
        channel_path = session_path+'/'+channel+'.timd'
        
        data = pymef3_file.read_mef_ts_data(channel_path,password, sample_ss[0], sample_ss[1])
        data_list.append(data)
        
    return data_list

#def read_ts_data_uutc(session_md,session_path,
#                      password,channel_map,uutc_map,
#                      on_out_of_channel_fail = True):
#    """
#    Reads desired channels in desired time segment. Missing data at
#    discontinuities are filled with NaNs.
#    
#    Parameters:
#    -----------
#    session_md - mef3 session metadata dictionary\n
#    session_path - path to mef3 session (.mefd)\n
#    password - mef3 data password\n
#    channel_map - list of channels to be read\n
#    uutc_map - list of [start,stop] uutc times to be loaded that correspond\n
#        to channel_map. if there is only one entry the same range is applied\n
#        to all channels\n
#    on_out_of_channel_fail - fail if uUTC times are out of channel times\n
#        default = True
#    
#    Returns:
#    --------
#    data - numpy array [channels,samples]\n
#    """
#    
#    data_list = []
#    
#    if len(uutc_map) == 1:
#        uutc_map = uutc_map*len(channel_map)
#        
#    if len(uutc_map) != len(uutc_map):
#        print('Length of uutc map is not equivalet to the length of channel map')
#        return
#    
#    # Checkt that uUTCs are in file
#    for channel,uutc_ss in zip(channel_map,uutc_map):
#        channel_path = session_path+'/'+channel+'.timd'
#        channel_md = session_md['time_series_channels'][channel]
#        
#        fs = channel_md['section_2']['sampling_frequency']
#        
#        # Get the size of the output array
#        N_samples_out = int(np.ceil(((uutc_ss[1] - uutc_ss[0]) / 1e6) * fs))
#        out_array = np.zeros(N_samples_out)
#        out_array[:] = np.nan
#        
#        # Get discontinuities
#        disc_uutc_list,disc_sample_list = get_discontinuities(channel_md)[1:]
#        
#        flat_uutcs = [i for sublist in disc_uutc_list for i in sublist]
#        flat_samples = [i for sublist in disc_sample_list for i in sublist]
#        
#        # Find out if we hit any discontinuities, if so - get them
#        disc_in_request = [x for x in zip(flat_uutcs,flat_samples) if (x[0] > uutc_ss[0] and x[0] < uutc_ss[1])]
#        
#        # Get the segment samples in mef and corresponding indices in output array
#        if len(disc_uuts_in_request):
#            out_array_idxs = []
#            mef_idxs = []
#            
#            # Get the first mef sample and check if we hit discontinuity
#            sample, disc_dist = sample_for_uutc(uutc_ss[0],channel_md,True)
#            
#            if disc_dist: #uutc was in discontinuity
#                array_idx = disc_dist
#            
#            disc_sample = disc_in_request[0][1]
#            
#            for disc_uutc in disc_uuts_in_request:
#                
#                
#        else:
#            out_array_idxs = [[0,N_samples_out]]
#            start_samp = sample_for_uutc(uutc_ss[0],channel_md)
#            stop_samp = sample_for_uutc(uutc_ss[1],channel_md)
#            mef_idxs[[start_samp,stop_samp]]
#        
#        
#       
#        
#        
#        
#        
#    return

def read_ts_channel_basic_info(mef3_session_path,password):
    """
    Reads session time series channel names
    
    Parameters:
    -----------
    mef3_session_path - path to mef3 session\n
    password - mef3 data password\n
    
    Returns:
    --------
    channel_list\n
    """
    
    session_ts_metadata_dict = pymef3_file.read_mef_session_metadata(mef3_session_path,password)
    
    channel_list = list(session_ts_metadata_dict['time_series_channels'].keys())
    channel_list.sort()
    
    channel_infos = []
    for channel in channel_list:
        
        fsamp = session_ts_metadata_dict['time_series_channels'][channel]['section_2']['sampling_frequency']
        nsamp = session_ts_metadata_dict['time_series_channels'][channel]['section_2']['number_of_samples']
        ufact = session_ts_metadata_dict['time_series_channels'][channel]['section_2']['units_conversion_factor']
        unit = session_ts_metadata_dict['time_series_channels'][channel]['section_2']['units_description']
        
        channel_infos.append({'name':channel,'fsamp':fsamp,'nsamp':nsamp,
                              'ufact':ufact,'unit':unit})
    
    return channel_infos
