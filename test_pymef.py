#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 16:00:20 2017

Testing script for pymef

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

from pymef import pymef3_file
import os
import numpy as np

# TODO - remember to try things for different password settings
# TODO - tests for memory leaks
# TODO - test the LNTP record where consulted with Matt
# DONE - test CSti record
# DONE - CRC checks - meflib takes care of this when reading!!!


# %% Presets

# Path to test session
#test_session_path = '/home/jan_cimbalnik/Dropbox/Source/C/Mef3_python/test_session/'
#test_session_path = '/Users/jan/Desktop/mef3_test/'
test_session_path = '/home/cimba/Dropbox/Source/C/Mef3_python/test_session/'


secs_to_write = 10
samps_per_mef_block = 5000
sampling_frequency = 5000
secs_to_append = 5
pass_1 = 'chair'
pass_2 = 'table'
# pass_1 = None
# pass_2 = None

# Clean up the session_path
os.system('rm -rf '+test_session_path)

start_time = 946684800000000 # midnight, 1 January 2000 from Dan's code
end_time = int(start_time + 1e6*secs_to_write)
record_time_1 = int(start_time + 1e6)
record_time_2 = int(start_time + 2*1e6)

# %% ---------- Mef write test ----------
print("\n\n---------- Writing mef files ----------\n\n")

# %% Write one data record file (with indices file)
#record_file_path = test_session_path+'c2_2.mefd/O1.timd/O1-000000.segd'
record_file_path = test_session_path+'msel_fnusa.mefd/msel.timd/msel-000000.segd/'

record_list = []

# Create Note
note_dict = {'type_string':'Note',
			 'version_major':1,
			 'version_minor':0,
			 'time':record_time_1,
			 'note':'Note_test'}

# Create SyLg
sylg_dict = {'type_string':'SyLg',
			 'version_major':1,
			 'version_minor':0,
			 'time':record_time_1,
			 'text':'SyLg_test'}

# Create EDFA
edfa_dict = {'type_string':'EDFA',
			 'version_major':1,
			 'version_minor':0,
			 'time':record_time_1,
			 'duration':1000000,
			 'annotation':'EDFA_test'}

# Not sure what is the template, ask Matt
# lntp_dict = {'type_string':'LNTP',
# 			 'version_major':1,
# 			 'version_minor':0,
# 			 'time':int(946684800000000 + 2*1e6),
# 			 'length':1000000,
# 			 'template':20}


# Create Seiz
seiz_chans = []
seiz_chan_dict_1 = {'name':'msel',
                    'onset':record_time_1,
                    'offset':record_time_2}
seiz_chans.append(seiz_chan_dict_1)

seiz_dict = {'type_string':'Seiz',
			 'version_major':1,
			 'version_minor':0,
			 'time':min([x['onset'] for x in seiz_chans]),
			 'earliest_onset':min([x['onset'] for x in seiz_chans]),
			 'latest_offset':max([x['onset'] for x in seiz_chans]),
			 'duration':max([x['onset'] for x in seiz_chans]) - min([x['onset'] for x in seiz_chans]),
			 'number_of_channels':len(seiz_chans),
			 'onset_code':2, # Decide what to do with his when reading - we can translate to string
			 'marker_name_1':'beer_tap',
			 'marker_name_2':'wine_tap',
			 'annotation':'test_seizure',
			 'channels':seiz_chans}

csti_dict = {'type_string':'CSti',
             'version_major':1,
             'version_minor':0,
             'time':record_time_1,
             'task_type':'beerdrink',
             'stimulus_duration':1000000,
             'stimulus_type':'pilsner',
             'patient_response':'hmmm'}

esti_dict = {'type_string':'ESti',
            'version_major':1,
            'version_minor':0,
            'time':record_time_1,
            'amplitude':1.5,
            'frequency':250.5,
            'pulse_width':100,
            'ampunit_code':1,
            'mode_code':2,
            'waveform':'nice',
            'anode':'positive',
            'catode':'negative'}

record_list.append(note_dict)
record_list.append(sylg_dict)
record_list.append(edfa_dict)

# record_list.append(lntp_dict)

record_list.append(seiz_dict)
record_list.append(csti_dict)
record_list.append(esti_dict)




pymef3_file.write_mef_data_records(record_file_path,
                              pass_1,
                              pass_2,
                              start_time,  
                              end_time,
                              record_list)

print("Records written at segment level")

record_file_path_channel = test_session_path+'msel_fnusa.mefd/msel.timd/'
pymef3_file.write_mef_data_records(record_file_path_channel,
                              pass_1,
                              pass_2,
                              start_time,
                              end_time,
                              record_list)

print("Records written at channel level")

record_file_path_channel = test_session_path+'msel_fnusa.mefd/'
pymef3_file.write_mef_data_records(record_file_path_channel,
                              pass_1,
                              pass_2,
                              start_time,
                              end_time,
                              record_list)

print("Records written at segment level")

# %% Write one time series channel metadata file
section3_dict = {'recording_time_offset':start_time,
                 'DST_start_time':0,
				'DST_end_time': 0,
				'GMT_offset':3600,
				'subject_name_1': 'Olaf',
                 	'subject_name_2': 'Mefson',
                 'subject_ID': '2017',
                 'recording_location':'pub'}

section2_ts_dict = {'channel_description':'Test_channel',
                    'session_description':'Test_session',
                    'recording_duration':1,# Test None here
                    'reference_description':'wine',
                    'acquisition_channel_number': 5,
                    'sampling_frequency': sampling_frequency,
                    'notch_filter_frequency_setting':50.0,
                    'low_frequency_filter_setting':1.0,
                    'high_frequency_filter_setting':10.0,
                    'AC_line_frequency':70,
                    'units_conversion_factor': 1.5,
                    'units_description': 'uV', # Not wure if the info below should be specified by user
                    'maximum_native_sample_value':0.0,
                    'minimum_native_sample_value':0.0,
                    'start_sample':0, # WATCH OUT - start sample IS important!!!
                    'number_of_blocks':0,
                    'maximum_block_bytes':0,
                    'maximum_block_samples':0,
                    'maximum_difference_bytes':0,
                    'block_interval':0, # Will be entered when writing ts data
                    'number_of_discontinuities':1,
                    'maximum_contiguous_blocks':0,
                    'maximum_contiguous_block_bytes':0,
                    'maximum_contiguous_samples':0,
                    'number_of_samples':0}

# Note: the rest of the tmd2 fileds is subject to discussion
time_series_metadata_file_path = record_file_path
pymef3_file.write_mef_ts_metadata(time_series_metadata_file_path,
							 pass_1,
							 pass_2,
							 start_time,
							 end_time,
							 section2_ts_dict,
							 section3_dict)

print("Time series metadata written")

# %% Write time series data with indices file


# Generate data
raw_data = np.random.randint(-200,200,sampling_frequency*secs_to_write,dtype='int32')

time_series_data_files_path = record_file_path
pymef3_file.write_mef_ts_data_and_indices(time_series_data_files_path,
                             		 pass_1,
                             		 pass_2,
                             		 samps_per_mef_block,
                             		 raw_data,
                             		 0)

print("Time series data and indices written")

# %% Append time series data and modify files accordinglly, but update metadata first

raw_data_to_append = np.random.randint(-200,200,sampling_frequency*secs_to_append,dtype='int32')

pymef3_file.append_ts_data_and_indices(time_series_data_files_path,
                             		pass_1,
                             		pass_2,
                             		start_time,
                             		int(end_time + (1e6*secs_to_append)),
                             		samps_per_mef_block,
                             		raw_data_to_append,
                             		0)
raw_data = np.concatenate([raw_data,raw_data_to_append])

# %% Write one video channel metadata file
section2_v_dict={'channel_description':'Test_channel',
				'session_description':'Test_session',
				'recording_duration':1,
				'horizontal_resolution':10,
				'vertical_resolution':20,
				'frame_rate':60.0,
				'number_of_clips':1,
				'maximum_clip_bytes':5000,
				'video_format':'mpeg',
				'video_file_CRC':111111} # Some of these fields won't be probably needed

video_metadata_file_path = test_session_path+'c2_2.mefd/movie.vidd/movie-000000.segd'
video_metadata_file_path = test_session_path+'msel_fnusa.mefd/fnusa.vidd/fnusa-000000.segd'
pymef3_file.write_mef_v_metadata(video_metadata_file_path,
							pass_1,
							pass_2,
							start_time,
							end_time,
							section2_v_dict,
							section3_dict)

print("Video metadata written")


# %% Write video indices file

index_entry = {'start_time':start_time,
               'end_time':end_time,
               'start_frame':0,
               'end_frame':300,
               'file_offset':0,
               'clip_bytes':5000}

index_entries=[index_entry]

video_indices_file_path = video_metadata_file_path
pymef3_file.write_mef_v_indices(video_indices_file_path,
                           pass_1,
                           pass_2,
                           start_time, # min from index entries
                           end_time, # max from index entries
                           index_entries)

print("Video indices written")



# %% ---------- Mef read test ----------

print("\n\n---------- Reading mef files ----------\n\n")

# %% Read time series segment metadata (including records)
# segment_directory = '/home/jan_cimbalnik/Desktop/edf2mef3/c2_2.mefd/O1.timd/O1-000000.segd'
segment_directory = test_session_path+'msel_fnusa.mefd/msel.timd/msel-000000.segd'
segment_ts_metadata_dict = pymef3_file.read_mef_segment_metadata(segment_directory,pass_2)

print("Time series segment read")

# %% Read time series channel metadata
channel_directory = test_session_path+'msel_fnusa.mefd/msel.timd'
channel_ts_metadata_dict = pymef3_file.read_mef_channel_metadata(channel_directory,pass_2)

print("Time series channel read")

# %% Read video segment metadata
segment_directory = test_session_path+'msel_fnusa.mefd/fnusa.vidd/fnusa-000000.segd'
segment_v_metadata_dict = pymef3_file.read_mef_segment_metadata(segment_directory,pass_2)

print("Video segment read")

# %% Read video channel metadata
channel_directory = test_session_path+'msel_fnusa.mefd/fnusa.vidd'
channel_v_metadata_dict = pymef3_file.read_mef_channel_metadata(channel_directory,pass_2)

print("Video channel read")

# %% Read session
session_directory = test_session_path+'msel_fnusa.mefd'
session_metadata_dict = pymef3_file.read_mef_session_metadata(session_directory,pass_2)

print("Session read")


# %% Read time series data

data = pymef3_file.read_mef_ts_data(test_session_path+'msel_fnusa.mefd/msel.timd',pass_2, None, None)

print("Time series data read")

# %% ---------- Mef read/write comparison ----------

print("\n\n---------- Comparing input and output ----------\n\n")

# Check the records goind in and out
records = segment_ts_metadata_dict['records_info']['records']

if len(record_list) == len(records):
    print("Number of recrods OK")
else:
    print("MISMATCH! Different number of written and read records")

for rec_id in list(range(len(record_list))):
    write_record = record_list[rec_id]
    read_record = segment_ts_metadata_dict['records_info']['records'][rec_id]
    
    print('Record type: ---'+write_record['type_string']+'---')
    
    # Record header
    if ((write_record['time'] == read_record['time']) &
        (write_record['type_string'] == read_record['type_string'])):
        print('Write = read in record header '+str(rec_id)+' OK')
    else:
        print('MISMATCH! Write != read in record header '+str(rec_id))
        
    # Record body
    if write_record['type_string'] == 'EDFA':
        if ((write_record['duration'] == read_record['record_body']['duration']) &
            (write_record['annotation'] == read_record['record_body']['annotation'])):
            print('Write = read in record body '+str(rec_id)+' OK')
        else:
            print('MISMATCH! Write != read in record body '+str(rec_id))

    if write_record['type_string'] == 'Note':
        if (write_record['note'] == read_record['record_body']['note']):
            print('Write = read in record body '+str(rec_id)+' OK')
        else:
            print('MISMATCH! Write != read in record body '+str(rec_id))
            
    if write_record['type_string'] == 'SyLg':
        if (write_record['text'] == read_record['record_body']['text']):
            print('Write = read in record body '+str(rec_id)+' OK')
        else:
            print('MISMATCH! Write != read in record body '+str(rec_id))

    if write_record['type_string'] == 'Seiz':
        if ((write_record['earliest_onset'] == read_record['record_body']['earliest_onset_uUTC']) &
        		(write_record['latest_offset'] == read_record['record_body']['latest_offset_uUTC']) &
        		(write_record['duration'] == read_record['record_body']['duration']) &
        		(write_record['number_of_channels'] == read_record['record_body']['number_of_channels']) &
        		(write_record['onset_code'] == read_record['record_body']['onset_code']) &
        		(write_record['marker_name_1'] == read_record['record_body']['marker_name_1']) &
        		(write_record['marker_name_2'] == read_record['record_body']['marker_name_2']) &
        		(write_record['annotation'] == read_record['record_body']['annotation'])):
                print('Write = read in record body '+str(rec_id)+' OK')
        else:
            print('MISMATCH! Write != read in record body '+str(rec_id))
            
        # Check the channel entries
        for ci, write_channel in enumerate(write_record['channels']):
            if ((write_channel['name'] == read_record['record_body']['channels'][ci]['name']) &
                (write_channel['onset'] == read_record['record_body']['channels'][ci]['onset_uUTC']) &
                (write_channel['offset'] == read_record['record_body']['channels'][ci]['offset_uUTC'])):
                    print('Write == read for seizure channel '+str(ci)+' OK')
            else:
                print('MISMATCH! Write != read in seizure channel '+str(ci))

    if write_record['type_string'] == 'CSti':
        if ((write_record['task_type'] == read_record['record_body']['task_type']) &
            (write_record['stimulus_duration'] == read_record['record_body']['stimulus_duration']) &
            (write_record['stimulus_type'] == read_record['record_body']['stimulus_type']) &
            (write_record['patient_response'] == read_record['record_body']['patient_response'])):
            print('Write = read in record body '+str(rec_id)+' OK')
        else:
            print('MISMATCH! Write != read in record body '+str(rec_id))

    if write_record['type_string'] == 'ESti':
        if ((write_record['amplitude'] == read_record['record_body']['amplitude']) &
            (write_record['frequency'] == read_record['record_body']['frequency']) &
            (write_record['pulse_width'] == read_record['record_body']['pulse_width']) &
            (write_record['ampunit_code'] == read_record['record_body']['ampunit_code']) &
            (write_record['mode_code'] == read_record['record_body']['mode_code']) &
            (write_record['waveform'] == read_record['record_body']['waveform']) &
            (write_record['anode'] == read_record['record_body']['anode']) &
            (write_record['catode'] == read_record['record_body']['catode'])):
            print('Write = read in record body '+str(rec_id)+' OK')
        else:
            print('MISMATCH! Write != read in record body '+str(rec_id))


# Check the metadata from time series

mismatch_flag = False
for md2_user_key in ['channel_description','session_description','recording_duration',
					'reference_description','acquisition_channel_number','sampling_frequency',
                     'notch_filter_frequency_setting','low_frequency_filter_setting',
					'high_frequency_filter_setting','AC_line_frequency','units_conversion_factor',
                     	'units_description']:
    if section2_ts_dict[md2_user_key] != segment_ts_metadata_dict['section_2'][md2_user_key]:
        print('MISMATCH! Section 2 "'+md2_user_key+'" is different ')
        mismatch_flag = True
        
if mismatch_flag == False: print('User defined time series metadata section 2 OK')

mismatch_flag = False
for md3_user_key in ['recording_time_offset','DST_start_time','DST_end_time',
                     'GMT_offset','subject_name_1','subject_name_2','subject_ID','recording_location']:
    
    if section3_dict[md3_user_key] != segment_ts_metadata_dict['section_3'][md3_user_key]:
        print('MISMATCH! Section 3 "'+md3_user_key+'" is different ')
        mismatch_flag = True

if mismatch_flag == False: print('User defined time series metadata section 3 OK')

mismatch_flag = False

if segment_ts_metadata_dict['section_2']['maximum_native_sample_value'] != np.max(data) * section2_ts_dict['units_conversion_factor']:
    print('MISMATCH! Section 2 "maximum_native_sample_value" is different than it should be')
    mismatch_flag = True

if segment_ts_metadata_dict['section_2']['minimum_native_sample_value'] != np.min(data) * section2_ts_dict['units_conversion_factor']:
    print('MISMATCH! Section 2 "minimum_native_sample_value" is different than it should be')
    mismatch_flag = True
    
# This field should be tested across segments
if segment_ts_metadata_dict['section_2']['start_sample'] != section2_ts_dict['start_sample']:
    print('MISMATCH! Section 2 "start_sample" is different than it should be')
    mismatch_flag = True
    
if segment_ts_metadata_dict['section_2']['number_of_blocks'] != np.ceil((sampling_frequency * (secs_to_write + secs_to_append)) / samps_per_mef_block):
    print('MISMATCH! Section 2 "number_of_blocks" is different than it should be')
    mismatch_flag = True
    
if segment_ts_metadata_dict['section_2']['maximum_block_samples'] != samps_per_mef_block:
    print('MISMATCH! Section 2 "maximum_block_samples" is different than it should be')
    mismatch_flag = True
    
if segment_ts_metadata_dict['section_2']['number_of_samples'] != sampling_frequency * (secs_to_write + secs_to_append):
    print('MISMATCH! Section 2 "number_of_samples" is different than it should be')
    mismatch_flag = True

if mismatch_flag == False: print('Data defined time series metadata section 2 OK')


# Check the metadata from video
mismatch_flag = False
for md2_user_key in section2_v_dict.keys():
    if section2_v_dict[md2_user_key] != segment_v_metadata_dict['section_2'][md2_user_key]:
        print('MISMATCH! Section 2 "'+md2_user_key+'" is different ')
        print(section2_v_dict[md2_user_key])
        print(segment_v_metadata_dict['section_2'][md2_user_key])
        mismatch_flag = True
        
if mismatch_flag == False: print('User defined video metadata section 2 OK')

mismatch_flag = False
for md3_user_key in ['recording_time_offset','DST_start_time','DST_end_time',
                     'GMT_offset','subject_name_1','subject_name_2','subject_ID','recording_location']:
    
    if section3_dict[md3_user_key] != segment_v_metadata_dict['section_3'][md3_user_key]:
        print('MISMATCH! Section 3 "'+md3_user_key+'" is different ')
        print(section3_dict[md3_user_key])
        print(segment_v_metadata_dict['section_3'][md3_user_key])
        mismatch_flag = True

if mismatch_flag == False: print('User defined video metadata section 3 OK')

# Check the time series data
data_diff = np.sum(raw_data - data)
if data_diff:
    print('MISMATCH! Data I/O failed, diff: '+str(data_diff))
else:
    print('Data I/O OK')

# %% ---------- Mef helper function test ----------

print("\n\n---------- Testing helper functions ----------\n\n")

ts_metadata_file = test_session_path+'msel_fnusa.mefd/msel.timd/msel-000000.segd/'+'msel-000000.tmet'

# Test wrong password
result = pymef3_file.check_mef_password(ts_metadata_file,'bu')
if result == 0:
    print('Wrong password check OK!')
else:
    print('Wrong password check failed!')

# Test level 1 password
result = pymef3_file.check_mef_password(ts_metadata_file,pass_1)
if result == 1:
    print('Level 1 password check OK!')
else:
    print('Level 1 password check failed!')

# Test level 2 password
result = pymef3_file.check_mef_password(ts_metadata_file,pass_2)
if result == 2:
    print('Level 2 password check OK!')
else:
    print('Level 2 password check failed!')
