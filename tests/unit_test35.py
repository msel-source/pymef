#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 09:23:17 2017

Unit testing for pymef library.

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

# Standard library imports
import unittest
import tempfile
import warnings

# Third party imports
import numpy as np

# Local imports
from pymef import pymef3_file


class TestStringMethods(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        
        self.temp_dir = tempfile.TemporaryDirectory(suffix='.mefd')
        self.mef_session_path = self.temp_dir.name
        
        # Some presets
        self.secs_to_write = 10
        self.samps_per_mef_block = 5000
        self.sampling_frequency = 5000
        self.secs_to_append = 5
        self.discont_length = 2
        self.secs_to_seg2 = 5
        
        self.pwd_1 = 'chair'
        self.pwd_2 = 'table'
        
        self.start_time = 946684800000000
        self.end_time = int(self.start_time + 1e6*self.secs_to_write)
        self.record_time_1 = int(self.start_time + 1e6)
        self.record_time_2 = int(self.start_time + 2*1e6)
        self.rec_offset = int(self.start_time - 1e6)
        
        # Create paths for channels and segments
        self.ts_channel_path = self.mef_session_path+'/ts_channel.timd'
        
        self.ts_seg1_path = self.ts_channel_path+'/ts_channel-000000.segd'
        self.ts_seg2_path = self.ts_channel_path+'/ts_channel-000001.segd'
        
        self.vid_channel_path = self.mef_session_path+'/vid_channel.vidd'
        
        self.vid_seg1_path = self.vid_channel_path+'/vid_channel-000000.segd'
        
        
        # Prepare dummy records
        
        self.record_list = []
        
        # Create Note
        note_dict = {'type_string':'Note',
        			 'version_major':1,
        			 'version_minor':0,
        			 'time':self.record_time_1,
        			 'note':'Note_test'}
        
        # Create SyLg
        sylg_dict = {'type_string':'SyLg',
        			 'version_major':1,
        			 'version_minor':0,
        			 'time':self.record_time_1,
        			 'text':'SyLg_test'}
        
        # Create EDFA
        edfa_dict = {'type_string':'EDFA',
        			 'version_major':1,
        			 'version_minor':0,
        			 'time':self.record_time_1,
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
                            'onset':self.record_time_1,
                            'offset':self.record_time_2}
        seiz_chans.append(seiz_chan_dict_1)
        
        seiz_time = min([x['onset'] for x in seiz_chans])
        earliest_onset = min([x['onset'] for x in seiz_chans])
        latest_offset = max([x['offset'] for x in seiz_chans])
        
        seiz_dict = {'type_string':'Seiz',
        			 'version_major':1,
        			 'version_minor':0,
        			 'time':seiz_time,
        			 'earliest_onset':earliest_onset,
        			 'latest_offset':latest_offset,
        			 'duration':latest_offset - earliest_onset,
        			 'number_of_channels':len(seiz_chans),
        			 'onset_code':2, # Decide what to do with his when reading - we can translate to string
        			 'marker_name_1':'beer_tap',
        			 'marker_name_2':'wine_tap',
        			 'annotation':'test_seizure',
        			 'channels':seiz_chans}
        
        csti_dict = {'type_string':'CSti',
                     'version_major':1,
                     'version_minor':0,
                     'time':self.record_time_1,
                     'task_type':'beerdrink',
                     'stimulus_duration':1000000,
                     'stimulus_type':'pilsner',
                     'patient_response':'hmmm'}
        
        esti_dict = {'type_string':'ESti',
                    'version_major':1,
                    'version_minor':0,
                    'time':self.record_time_1,
                    'amplitude':1.5,
                    'frequency':250.5,
                    'pulse_width':100,
                    'ampunit_code':1,
                    'mode_code':2,
                    'waveform':'nice',
                    'anode':'positive',
                    'catode':'negative'}
        
        self.record_list.append(note_dict)
        self.record_list.append(sylg_dict)
        self.record_list.append(edfa_dict)
        
        # record_list.append(lntp_dict)
        
        self.record_list.append(seiz_dict)
        self.record_list.append(csti_dict)
        self.record_list.append(esti_dict)
        
        # Prepare dummy time series metadata
        self.section3_dict = {'recording_time_offset':self.rec_offset,
                              'DST_start_time':0,
                              'DST_end_time': 0,
                              'GMT_offset':3600,
                              'subject_name_1': 'Olaf',
                              'subject_name_2': 'Mefson',
                              'subject_ID': '2017',
                              'recording_location':'pub'}
        
        self.section2_ts_dict = {'channel_description':'Test_channel',
                                 'session_description':'Test_session',
                                 'recording_duration':1,# TODO: test 0 and None here
                                 'reference_description':'wine',
                                 'acquisition_channel_number': 5,
                                 'sampling_frequency': self.sampling_frequency,
                                 'notch_filter_frequency_setting':50.0,
                                 'low_frequency_filter_setting':1.0,
                                 'high_frequency_filter_setting':10.0,
                                 'AC_line_frequency':70,
                                 'units_conversion_factor': 1.5,
                                 'units_description': 'uV',
                                 'maximum_native_sample_value':0.0,
                                 'minimum_native_sample_value':0.0,
                                 'start_sample':0, # WATCH OUT - start sample IS important!!!
                                 'number_of_blocks':0,
                                 'maximum_block_bytes':0,
                                 'maximum_block_samples':0,
                                 'maximum_difference_bytes':0,
                                 'block_interval':0,
                                 'number_of_discontinuities':1,
                                 'maximum_contiguous_blocks':0,
                                 'maximum_contiguous_block_bytes':0,
                                 'maximum_contiguous_samples':0,
                                 'number_of_samples':0}
        
        # Second segment
        self.section2_ts_dict_seg2 = self.section2_ts_dict.copy()
        seg2_start_sample = self.sampling_frequency*(self.secs_to_write
                                                     + self.secs_to_append)
        self.section2_ts_dict_seg2['start_sample'] = seg2_start_sample
        
        # Prepare dummy time series data
        N = self.sampling_frequency * self.secs_to_write
        self.raw_data = np.random.randint(-200,200,N,dtype='int32')
        
        # Prepare dummy data for appending to test appedn function
        N = self.sampling_frequency * self.secs_to_append
        self.raw_data_to_append = np.random.randint(-200,200,N,dtype='int32')
        
        self.raw_data_seg_1 = np.concatenate([self.raw_data,
                                             self.raw_data_to_append])
        
        # Second segment data
        N = self.sampling_frequency * self.secs_to_seg2
        self.raw_data_seg_2 = np.random.randint(-200,200,N,dtype='int32')
        
        self.raw_data_all = np.concatenate([self.raw_data_seg_1,
                                            self.raw_data_seg_2])
        
        # Preapare dummy video metadata and indices
        
        self.section2_v_dict={'channel_description':'Test_channel',
                              'session_description':'Test_session',
                              'recording_duration':1,
                              'horizontal_resolution':10,
                              'vertical_resolution':20,
                              'frame_rate':60.0,
                              'number_of_clips':1,
                              'maximum_clip_bytes':5000,
                              'video_format':'mpeg',
                              'video_file_CRC':111111}
        
        self.v_index_entry = {'start_time':self.start_time,
                              'end_time':self.end_time,
                              'start_frame':0,
                              'end_frame':300,
                              'file_offset':0,
                              'clip_bytes':5000}
        
        self.v_index_entries=[self.v_index_entry]
        
        
        # Write dummy records
        
        pymef3_file.write_mef_data_records(self.ts_seg1_path,
                                           self.pwd_1,
                                           self.pwd_2,
                                           self.start_time,  
                                           self.end_time,
                                           self.rec_offset,
                                           self.record_list)
        
#        print("Records written at segment level")
        
        pymef3_file.write_mef_data_records(self.ts_channel_path,
                                           self.pwd_1,
                                           self.pwd_2,
                                           self.start_time,  
                                           self.end_time,
                                           self.rec_offset,
                                           self.record_list)
        
#        print("Records written at channel level")        
        pymef3_file.write_mef_data_records(self.mef_session_path,
                                           self.pwd_1,
                                           self.pwd_2,
                                           self.start_time,  
                                           self.end_time,
                                           self.rec_offset,
                                           self.record_list)
        
#        print("Records written at session level")
        
        # Write dummy time series metadata
        # Note: the rest of the tmd2 fileds is subject to discussion
        pymef3_file.write_mef_ts_metadata(self.ts_seg1_path,
                                          self.pwd_1,
                                          self.pwd_2,
                                          self.start_time,  
                                          self.end_time,
                                          self.section2_ts_dict,
                                          self.section3_dict)
        
        # Write second segment
        seg2_start = int(self.end_time
                         + (1e6*self.secs_to_append)
                         + int(1e6*self.discont_length))
        seg2_stop = seg2_start + int(self.secs_to_seg2 * 1e6)
        self.section2_ts_dict_seg2 = self.section2_ts_dict.copy()
        start_samp = self.sampling_frequency * (self.secs_to_write
                                                + self.secs_to_append)
        self.section2_ts_dict_seg2['start_sample'] = start_samp

        pymef3_file.write_mef_ts_metadata(self.ts_seg2_path,
                                          self.pwd_1,
                                          self.pwd_2,
                                          seg2_start,
                                          seg2_stop,
                                          self.section2_ts_dict_seg2,
                                          self.section3_dict)
        
#        print("Time series metadata written")
        
        # Write dummy time series data
        
        # Write first segment
        pymef3_file.write_mef_ts_data_and_indices(self.ts_seg1_path,
                                                  self.pwd_1,
                                                  self.pwd_2,
                                                  self.samps_per_mef_block,
                                                  self.raw_data,
                                                  0)
        
        # Append time series data and modify files accordinglly, but update metadata first
        append_start = self.end_time
        append_stop = int(append_start + (self.secs_to_append * 1e6))
        pymef3_file.append_ts_data_and_indices(self.ts_seg1_path,
                                               self.pwd_1,
                                               self.pwd_2,
                                               append_start,
                                               append_stop,
                                               self.samps_per_mef_block,
                                               self.raw_data_to_append,
                                               0)
        
        # Write second segment
        
        pymef3_file.write_mef_ts_data_and_indices(self.ts_seg2_path,
                                                  self.pwd_1,
                                                  self.pwd_2,
                                                  self.samps_per_mef_block,
                                                  self.raw_data_seg_2,
                                                  0)
                
#        print("Time series data and indices written")
        
        
        # Write dummy video metadata and indices
        pymef3_file.write_mef_v_metadata(self.vid_seg1_path,
                                         self.pwd_1,
                                         self.pwd_2,
                                         self.start_time,
                                         self.end_time,
                                         self.section2_v_dict,
                                         self.section3_dict)
        
#        print("Video metadata written")
        
        pymef3_file.write_mef_v_indices(self.vid_seg1_path,
                                        self.pwd_1,
                                        self.pwd_2,
                                        self.start_time, # min from index entries
                                        self.end_time, # max from index entries
                                        self.v_index_entries)
        
#        print("Video indices written")
        
        # Read back session metadata (avoids reading metadata in each function)
        self.smd = pymef3_file.read_mef_session_metadata(self.mef_session_path,
                                                         self.pwd_2)
        
    # ----- Read metadata tests -----
    
    def test_read_ts_segment_metadata(self):
        pymef3_file.read_mef_segment_metadata(self.ts_seg1_path,
                                              self.pwd_2)
        
    def test_read_ts_channel_metadata(self):
        pymef3_file.read_mef_channel_metadata(self.ts_channel_path,
                                              self.pwd_2)
        
    def test_read_vid_segment_metadata(self):
        pymef3_file.read_mef_segment_metadata(self.vid_seg1_path,
                                              self.pwd_2)
        
    def test_read_vid_channel_metadata(self):
        pymef3_file.read_mef_channel_metadata(self.vid_channel_path,
                                              self.pwd_2)
        
    def test_read_session_metadata(self):
        pymef3_file.read_mef_session_metadata(self.mef_session_path,
                                              self.pwd_2)
        
        
    # ----- Mef write / read comparison -----
    
    def test_record_reading(self):
        
        segments = self.smd['time_series_channels']['ts_channel']['segments']
        seg_md = segments['ts_channel-000000']        
        read_records = seg_md['records_info']['records']
    
        self.assertEqual(len(self.record_list), len(read_records))

        for rec_id in range(len(self.record_list)):
            write_record = self.record_list[rec_id]
            read_record = read_records[rec_id]
            
#            print('Record type: ---'+write_record['type_string']+'---')
            
            # Record header
            self.assertEqual(write_record['time'],
                             read_record['time'])
            self.assertEqual(write_record['type_string'],
                             read_record['type_string'])
                
            # Record body
            if write_record['type_string'] == 'EDFA':
                self.assertEqual(write_record['duration'],
                                 read_record['record_body']['duration'])
                self.assertEqual(write_record['annotation'],
                                 read_record['record_body']['annotation'])

            if write_record['type_string'] == 'Note':
                self.assertEqual(write_record['note'],
                                 read_record['record_body']['note'])
            
            if write_record['type_string'] == 'SyLg':
                self.assertEqual(write_record['text'],
                                 read_record['record_body']['text'])
                
            if write_record['type_string'] == 'Seiz':
                self.assertEqual(write_record['earliest_onset'],
                                 read_record['record_body']['earliest_onset_uUTC'])
                self.assertEqual(write_record['latest_offset'],
                                 read_record['record_body']['latest_offset_uUTC'])
                self.assertEqual(write_record['duration'],
                                 read_record['record_body']['duration'])
                self.assertEqual(write_record['number_of_channels'],
                                 read_record['record_body']['number_of_channels'])
                self.assertEqual(write_record['onset_code'],
                                 read_record['record_body']['onset_code'])
                self.assertEqual(write_record['marker_name_1'],
                                 read_record['record_body']['marker_name_1'])
                self.assertEqual(write_record['marker_name_2'],
                                 read_record['record_body']['marker_name_2'])
                self.assertEqual(write_record['annotation'],
                                 read_record['record_body']['annotation'])
           
                # Check the channel entries
                for ci, write_channel in enumerate(write_record['channels']):
                    read_channel = read_record['record_body']['channels'][ci]
                    
                    self.assertEqual(write_channel['name'],
                                     read_channel['name'])
                    self.assertEqual(write_channel['onset'],
                                     read_channel['onset_uUTC'])
                    self.assertEqual(write_channel['offset'],
                                     read_channel['offset_uUTC'])

        
            if write_record['type_string'] == 'CSti':
                self.assertEqual(write_record['task_type'],
                                 read_record['record_body']['task_type'])
                self.assertEqual(write_record['stimulus_duration'],
                                 read_record['record_body']['stimulus_duration'])
                self.assertEqual(write_record['stimulus_type'],
                                 read_record['record_body']['stimulus_type'])
                self.assertEqual(write_record['patient_response'],
                                 read_record['record_body']['patient_response'])
                
        
            if write_record['type_string'] == 'ESti':
                self.assertEqual(write_record['amplitude'],
                                 read_record['record_body']['amplitude'])
                self.assertEqual(write_record['frequency'],
                                 read_record['record_body']['frequency'])
                self.assertEqual(write_record['pulse_width'],
                                 read_record['record_body']['pulse_width'])
                self.assertEqual(write_record['ampunit_code'],
                                 read_record['record_body']['ampunit_code'])
                self.assertEqual(write_record['mode_code'],
                                 read_record['record_body']['mode_code'])
                self.assertEqual(write_record['waveform'],
                                 read_record['record_body']['waveform'])
                self.assertEqual(write_record['anode'],
                                 read_record['record_body']['anode'])
                self.assertEqual(write_record['catode'],
                                 read_record['record_body']['catode'])
        
    def test_time_series_metadata_section_2_usr(self):
        
        segments = self.smd['time_series_channels']['ts_channel']['segments']
        seg_md = segments['ts_channel-000000']    
        
        # User specified fields section 2
        section_2_usr_field_list = ['channel_description',
                                    'session_description',
                                    'recording_duration',
                                    'reference_description',
                                    'acquisition_channel_number',
                                    'sampling_frequency',
                                    'notch_filter_frequency_setting',
                                    'low_frequency_filter_setting',
                                    'high_frequency_filter_setting',
                                    'AC_line_frequency',
                                    'units_conversion_factor',
                                    'units_description']
        
        for md2_user_key in section_2_usr_field_list:
#            print(md2_user_key,self.section2_ts_dict[md2_user_key])
            self.assertEqual(self.section2_ts_dict[md2_user_key],
                             seg_md['section_2'][md2_user_key])
            
    def test_time_series_metadata_section_2_auto(self):
    
        segments = self.smd['time_series_channels']['ts_channel']['segments']
        seg_md = segments['ts_channel-000000']        
        
        # Fields that are created by C code during data writing
        max_sample_val = (np.max(self.raw_data_seg_1)
                          * self.section2_ts_dict['units_conversion_factor'])
        self.assertEqual(max_sample_val,
                         seg_md['section_2']['maximum_native_sample_value'])
        
        min_sample_val = (np.min(self.raw_data_seg_1)
                          * self.section2_ts_dict['units_conversion_factor'])
        self.assertEqual(min_sample_val,
                         seg_md['section_2']['minimum_native_sample_value'])
        
        # TODO: This field should be tested across segments
        self.assertEqual(self.section2_ts_dict['start_sample'],
                         seg_md['section_2']['start_sample'])
        
        min_sample_val = (np.min(self.raw_data_seg_1)
                          * self.section2_ts_dict['units_conversion_factor'])
        
        N_blocks = np.ceil((self.sampling_frequency 
                            * (self.secs_to_write + self.secs_to_append))
                                / self.samps_per_mef_block)
        self.assertEqual(N_blocks,
                         seg_md['section_2']['number_of_blocks'])
        
        self.assertEqual(self.samps_per_mef_block,
                         seg_md['section_2']['maximum_block_samples'])
        
        N_samples = (self.sampling_frequency 
                     * (self.secs_to_write + self.secs_to_append))
        self.assertEqual(N_samples,
                         seg_md['section_2']['number_of_samples'])
            
    def test_time_series_metadata_section_3(self):
        
        segments = self.smd['time_series_channels']['ts_channel']['segments']
        seg_md = segments['ts_channel-000000']
        
        for md3_key in self.section3_dict.keys():
            self.assertEqual(self.section3_dict[md3_key],
                             seg_md['section_3'][md3_key])
    
    def test_video_metadata_section_2(self):
        
        segments = self.smd['video_channels']['vid_channel']['segments']
        seg_md = segments['vid_channel-000000']

        for md2_key in self.section2_v_dict.keys():
            self.assertEqual(self.section2_v_dict[md2_key],
                             seg_md['section_2'][md2_key])
            
    def test_video_metadata_section_3(self):
        segments = self.smd['video_channels']['vid_channel']['segments']
        seg_md = segments['vid_channel-000000']

        for md3_key in self.section3_dict.keys():
            self.assertEqual(self.section3_dict[md3_key],
                             seg_md['section_3'][md3_key])
          
    def test_time_series_data(self):
                
        read_data = pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                                 self.pwd_2, 
                                                 None,
                                                 None)
        
        # Check the sums
        self.assertEqual(np.sum(self.raw_data_all),
                         np.sum(read_data))


        
    # ----- Data reading tests -----
    
    # Reading by sample
    
    def test_start_sample_out_of_file(self):
        
        ch_md = self.smd['time_series_channels']['ts_channel']
        
        ch_md2 = ch_md['section_2']
        
        start = ch_md2['start_sample'] - self.samps_per_mef_block
        end = ch_md2['start_sample'] + self.samps_per_mef_block
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Start sample smaller than 0. '
                            'Setting start sample to 0')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end)
            
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
        
    def test_end_sample_out_of_file(self):
                
        ch_md = self.smd['time_series_channels']['ts_channel']
        
        ch_md2 = ch_md['section_2']
        
        start = ch_md2['number_of_samples'] - self.samps_per_mef_block
        end = ch_md2['number_of_samples'] + self.samps_per_mef_block
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Stop sample larger than number of samples. '
                            'Setting end sample to number of samples in '
                            'channel')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end)
                
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
    
    def test_start_end_sample_before_file(self):
                
        ch_md = self.smd['time_series_channels']['ts_channel']
        
        ch_md2 = ch_md['section_2']
        
        start = ch_md2['start_sample'] - (self.samps_per_mef_block * 2)
        end = ch_md2['start_sample'] - self.samps_per_mef_block
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Start and stop samples are out of file. '
                            'Returning None')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end)
            
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
        
    def test_start_end_sample_after_file(self):
        
        ch_md = self.smd['time_series_channels']['ts_channel']
        
        ch_md2 = ch_md['section_2']
        
        start = ch_md2['number_of_samples'] + self.samps_per_mef_block
        end = ch_md2['number_of_samples'] + (self.samps_per_mef_block * 2)
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Start and stop samples are out of file. '
                            'Returning None')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end)
            
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
            
        return
    
    def test_start_sample_bigger_than_end_sample(self):
        
        error_text = 'Start sample larger than end sample, exiting...'
        
        try:
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         10,
                                         5)
        except Exception as e:
            self.assertEqual(error_text,str(e))
    
    # Reading by uutc
    
    def test_start_uutc_out_of_file(self):
        
        ch_md = self.smd['time_series_channels']['ts_channel']

        ch_md_spec = ch_md['channel_specific_metadata']
        
        start = int(ch_md_spec['earliest_start_time'] - 1e6)
        end = int(ch_md_spec['earliest_start_time'] + 1e6)
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Start uutc earlier than earliest start time. '
                            'Will insert NaNs')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end,
                                         True)
            
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
        
    
    def test_end_uutc_out_of_file(self):
        
        ch_md = self.smd['time_series_channels']['ts_channel']

        ch_md_spec = ch_md['channel_specific_metadata']
        
        start = int(ch_md_spec['latest_end_time'] - 1e6)
        end = int(ch_md_spec['latest_end_time'] + 1e6)
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Stop uutc later than latest end time. '
                            'Will insert NaNs')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end,
                                         True)
            
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
    
    def test_start_end_uutc_before_file(self):
        
        ch_md = self.smd['time_series_channels']['ts_channel']

        ch_md_spec = ch_md['channel_specific_metadata']
        
        start = int(ch_md_spec['earliest_start_time'] - (1e6 * 2))
        end = int(ch_md_spec['earliest_start_time'] - 1e6)
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Start and stop times are out of file. '
                            'Returning None')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end,
                                         True)
            
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
    
    def test_start_end_uutc_after_file(self):
        ch_md = self.smd['time_series_channels']['ts_channel']

        ch_md_spec = ch_md['channel_specific_metadata']
        
        start = int(ch_md_spec['latest_end_time'] + 1e6)
        end = int(ch_md_spec['latest_end_time'] + (1e6 * 2))
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            warning_text = ('Start and stop times are out of file. '
                            'Returning None')
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         start,
                                         end,
                                         True)
            
            self.assertEqual(len(w),1)
            self.assertEqual(warning_text,str(w[-1].message))
            assert issubclass(w[-1].category, RuntimeWarning)
    
    def test_start_uutc_in_discontinuity(self):
        
        discont_end_time = int(self.end_time + (1e6*self.secs_to_append)
                               + int(1e6*self.discont_length))
                
        
        start = int(discont_end_time - 5e5)
        end = int(discont_end_time + 5e5)
        

        data = pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                            self.pwd_2,
                                            start,
                                            end,
                                            True)
        
        # Check for the number of NaNs
        N_nans = (5e5 / 1e6) * self.sampling_frequency
        read_N_nans = np.sum(np.isnan(data))
        
        self.assertEqual(N_nans,read_N_nans)
        
            
    def test_end_uutc_in_discontinuity(self):
                
        discont_start_time = int(self.end_time + (1e6*self.secs_to_write))
                
        start = int(discont_start_time - self.secs_to_append*1e6)
        end = int(discont_start_time - self.secs_to_append*1e6 + 5e5)

        data = pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                            self.pwd_2,
                                            start,
                                            end,
                                            True)
        
        # Check for the number of NaNs
        N_nans = (5e5 / 1e6) * self.sampling_frequency
        N_nans -= 1 # End is not inclusive 
        read_N_nans = np.sum(np.isnan(data))
        
        self.assertEqual(N_nans,read_N_nans)
    
    def test_start_uutc_bigger_than_end_uutc(self):
        error_text = 'Start time later than end time, exiting...'
        
        try:
            pymef3_file.read_mef_ts_data(self.ts_channel_path,
                                         self.pwd_2,
                                         self.start_time + int(2*1e6),
                                         self.start_time + int(1*1e6),
                                         True)
        except Exception as e:
            self.assertEqual(error_text,str(e))


    # ----- Pymef helpers -----
    
    def test_wrong_password(self):
        ts_metadata_file = self.ts_seg1_path + '/ts_channel-000000.tmet'
        result = pymef3_file.check_mef_password(ts_metadata_file,'bu')
        self.assertEqual(-1,result)
        
    def test_level_1_password(self):
        ts_metadata_file = self.ts_seg1_path + '/ts_channel-000000.tmet'
        result = pymef3_file.check_mef_password(ts_metadata_file,self.pwd_1)
        self.assertEqual(1,result)
        
    def test_level_2_password(self):
        ts_metadata_file = self.ts_seg1_path + '/ts_channel-000000.tmet'
        result = pymef3_file.check_mef_password(ts_metadata_file,self.pwd_2)
        self.assertEqual(2,result)                        

if __name__ == '__main__':
    unittest.main()
