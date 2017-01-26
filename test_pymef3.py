#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 16:00:20 2017

Testing script for pymef3

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

import pymef3
import numpy as np

# TODO - remember to try things for different password settings

# Path to test session
#test_session_path = '/home/jan_cimbalnik/Desktop/edf2mef3/'
test_session_path = '/home/jan_cimbalnik/Dropbox/Source/C/Mef3_python/test_session/'

# %% ---------- Mef write test ----------


# %% Write one data record file (with indeces file)
#record_file_path = test_session_path+'c2_2.mefd/O1.timd/O1-000000.segd'
record_file_path = test_session_path+'msel_fnusa.mefd/msel.timd/msel-000000.segd/'

record_list = []

# Create EDFA
#TODO - test the rest of mef records
edfa_dict = {'type_string':'EDFA',
			 'version_major':1,
			 'version_minor':0,
			 'time':int(946684800000000 + 1e6),
			 'duration':1000000,
			 'annotation':'EDFA_test'}
record_list.append(edfa_dict)

#pymef3.write_mef_data_record(record_file_path,
#                              'chair',
#                              'table',
#                              946684800000000,  # midnight, 1 January 2000 from Dan's code
#                              int(946684800000000 + (1e6*5)),
#                              record_list)



# %% Write one time series channel metadata file
section3_dict = {'recording_time_offset':0,
				 'DST_start_time':946684800000000,
				 'DST_end_time': int(946684800000000 + (1e6*5)),
				 'GMT_offset':0,
				 'subject_name_1': 'Olaf',
             	 'subject_name_2': 'Mefson',
                 'subject_ID': '2017',
                 'recording_location':'pub'}

section2_ts_dict = {'channel_description':'Test_channel',
					'session_description':'Test_session',
					'recording_duration':1,# Test None here
					'reference_description':'wine',
					'acquisition_channel_number': 5,
					'sampling_frequency': 5000,
                    'notch_filter_frequency_setting':50.0,
					'low_frequency_filter_setting':1.0,
					'high_frequency_filter_setting':10.0,
					'AC_line_frequency':70,
                 	'units_conversion_factor': 1.5,
                 	'units_description': 'uV', # Not wure if the info below should be specified by user
                 	'maximum_native_sample_value':0,
                 	'minimum_native_sample_value':0,
                 	'start_sample':0,
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
pymef3.write_mef_ts_metadata(time_series_metadata_file_path,
							 'chair',
							 'table',
							 946684800000000,
							 int(946684800000000 + (1e6*5)),
							 section2_ts_dict,
							 section3_dict)



# %% Write one video channel metadata file
section2_v_dict={'channel_description':'Test_channel',
				 'session_description':'Test_session',
				 'recording_duration':1,
				 'horizontal_resolution':10,
				 'vertical_resolution':20,
				 'frame_rate':2,
				 'number_of_clips':2,
				 'maximum_clip_bytes':50,
				 'video_format':'mpeg',
				 'video_file_CRC':111111} # Some of these fields won't be probably needed

video_metadata_file_path = test_session_path+'c2_2.mefd/movie.vidd/movie-000000.segd'
video_metadata_file_path = test_session_path+'msel_fnusa.mefd/fnusa.vidd/fnusa-000000.segd'
#pymef3.write_mef_v_metadata(video_metadata_file_path,
#							'chair',
#							'table',
#							946684800000000,
#							int(946684800000000 + (1e6*5)),
#							section2_v_dict,
#							section3_dict)


# %% Write time series data with indeces file

# Generate data
raw_data = np.random.randint(-200,200,5000*10,dtype='int32')

time_series_data_files_path = record_file_path
pymef3.write_mef_ts_data_and_indeces(time_series_data_files_path,
                             		 'chair',
                             		 'table',
                             		 5000,
                             		 raw_data,
                             		 0)



# %% ---------- Mef read test ----------


# %% Read time series segment metadata (including records)
# segment_directory = '/home/jan_cimbalnik/Desktop/edf2mef3/c2_2.mefd/O1.timd/O1-000000.segd'
segment_directory = test_session_path+'msel_fnusa.mefd/msel.timd/msel-000000.segd'
segment_ts_metadata_dict = pymef3.read_mef_segment_metadata(segment_directory,'table')


# %% Read time series channel metadata
# channel_directory = test_session_path+'msel_fnusa.mefd/msel.timd'
# channel_ts_metadata_dict = pymef3.read_mef_channel_metadata(channel_directory,'table')

# %% Read video segment metadata
#segment_directory = test_session_path+'msel_fnusa.mefd/fnusa.vidd/fnusa-000000.segd'
#segment_v_metadata_dict = pymef3.read_mef_segment_metadata(segment_directory,'table')


# %% Read video channel metadata



# %% Read session
# session_directory = test_session_path+'msel_fnusa.mefd'
# session_metadata_dict = pymef3.read_mef_session_metadata(session_directory,'table')



# Read time series data
data = pymef3.read_mef_ts_data(test_session_path+'msel_fnusa.mefd/msel.timd','table')

# %% ---------- Mef read/write comparison ----------

# Check the records goind in and out

# Check the metadata from time series

# Check the metadata from video

# Check the time series data
data_diff = np.sum(raw_data - data)
if data_diff:
    print('Data I/O failed, diff: '+str(data_diff))
else:
    print('Data I/O OK')



