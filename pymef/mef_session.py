#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 5 14:21:39 2019

Mef session object.

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
import os
import struct
import shutil
import warnings
from multiprocessing import Pool

# Third party imports
import numpy as np

# Local imports
from .mef_file.pymef3_file import (read_mef_session_metadata,
                                   read_mef_ts_data,
                                   clean_mef_session_metadata,
                                   write_mef_ts_metadata,
                                   write_mef_v_metadata,
                                   write_mef_ts_data_and_indices,
                                   append_ts_data_and_indices,
                                   write_mef_v_indices,
                                   write_mef_data_records,
                                   create_tmd2_dtype,
                                   create_vmd2_dtype,
                                   create_md3_dtype,
                                   check_mef_password)
from .mef_constants import (
            METADATA_RECORDING_DURATION_NO_ENTRY,
            TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY,
            TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_NO_ENTRY,
            TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY,
            TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY,
            TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY,
            TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY,
            TIME_SERIES_METADATA_AC_LINE_FREQUENCY_NO_ENTRY,
            TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY,
            TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY,
            TIME_SERIES_METADATA_START_SAMPLE_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY,
            TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY,
            TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY,

            VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY,
            VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY,
            VIDEO_METADATA_FRAME_RATE_NO_ENTRY,
            VIDEO_METADATA_VIDEO_FILE_CRC_NO_ENTRY,

            METADATA_RECORDING_TIME_OFFSET_NO_ENTRY,
            METADATA_DST_START_TIME_NO_ENTRY,
            METADATA_DST_END_TIME_NO_ENTRY,
            GMT_OFFSET_NO_ENTRY
        )


class MefSession():
    """
    Basic object for operations with mef sessions.

    Parameters:
    -----------
    session_path - path to mef session
    password - password for mef session
    read_metadata - whether to read metadata (default=True)
    new_session - whether this is a new session for writing (default=False)
    """

    def __init__(self, session_path, password, read_metadata=True,
                 new_session=False):

        if not session_path.endswith('/'):
            session_path += '/'

        if '.mefd' != session_path[-6:-1]:
            raise ValueError("Session path must end with .mefd suffix!")

        if new_session:
            os.makedirs(session_path)
            self.session_md = None
            return

        if not os.path.exists(session_path):
            raise FileNotFoundError(session_path+' does not exist!')

        self.session_path = session_path
        self.password = password

        # Check if path exists
        if not os.path.exists(session_path):
            raise FileNotFoundError(session_path+' does not exist!')
        self._check_password()

        if read_metadata:
            self.session_md = read_mef_session_metadata(session_path,
                                                        password)
        else:
            self.session_md = None

        return

    # ----- Helper functions -----
    def _check_password(self):
        """
        Checks provided password on all files in the session

        Parameters:
        -----------
        mefpath - path to mef3 direcory\n
        password - mef3 data password\n

        Returns:
        --------
        None on success]n
        """
        mef_files = []
        for path, subdirs, files in os.walk(self.session_path):
            for name in files:
                mef_files.append(os.path.join(path, name))

        results = np.zeros(len(mef_files))
        for i, mef_file in enumerate(mef_files):
            results[i] = check_mef_password(mef_file, self.password)

        if np.any(results < 0):
            raise RuntimeError('MEF password is invalid')

        return

    def _get_channel_md(self, channel):
        """
        Paramters:
        ----------
        channel - required channel

        Returns:
        --------
        channel metadata structure
        """
        ts_chs = self.session_md['time_series_channels']
        return ts_chs[channel]['channel_specific_metadata']

    def _arg_merger(self, args):
        return read_mef_ts_data(*args)

    def reload(self):
        self.close()
        self.session_md = read_mef_session_metadata(self.session_path,
                                                    self.password)

    def close(self):
        if self.session_md is not None:
            clean_mef_session_metadata(
                self.session_md['session_specific_metadata'])
            self.session_md = None
        return

    # ----- Data writing functions -----
    def _initialize_tmd2(self):
        """
        Initialize time series metadata section 2 array with default values

        Returns:
        --------
        Initialized numpy array of time series md2 dtype
        """

        tmd2 = np.zeros(1, create_tmd2_dtype())

        tmd2['recording_duration'] = (
            METADATA_RECORDING_DURATION_NO_ENTRY)
        tmd2['acquisition_channel_number'] = (
            TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY)
        tmd2['number_of_samples'] = (
            TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_NO_ENTRY)
        tmd2['sampling_frequency'] = (
            TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY)
        tmd2['low_frequency_filter_setting'] = (
            TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY)
        tmd2['high_frequency_filter_setting'] = (
            TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY)
        tmd2['notch_filter_frequency_setting'] = (
            TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY)
        tmd2['AC_line_frequency'] = (
            TIME_SERIES_METADATA_AC_LINE_FREQUENCY_NO_ENTRY)
        tmd2['units_conversion_factor'] = (
            TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY)
        tmd2['maximum_native_sample_value'] = (
            TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY)
        tmd2['minimum_native_sample_value'] = (
            TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY)
        tmd2['start_sample'] = (
            TIME_SERIES_METADATA_START_SAMPLE_NO_ENTRY)
        tmd2['maximum_difference_bytes'] = (
            TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_NO_ENTRY)
        tmd2['maximum_block_samples'] = (
            TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY)
        tmd2['block_interval'] = (
            TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY)
        tmd2['maximum_contiguous_samples'] = (
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_NO_ENTRY)
        tmd2['maximum_contiguous_block_bytes'] = (
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY)
        tmd2['number_of_discontinuities'] = (
            TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY)

        return tmd2

    def _initialize_vmd2(self):
        """
        Initialize video metadata section 2 array with default values

        Returns:
        --------
        Initialized numpy array of video md2 dtype
        """

        vmd2 = np.zeros(1, create_vmd2_dtype())

        vmd2['horizontal_resolution'] = (
            VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY)
        vmd2['vertical_resolution'] = (
            VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY)
        vmd2['frame_rate'] = (
            VIDEO_METADATA_FRAME_RATE_NO_ENTRY)
        vmd2['video_file_CRC'] = (
            VIDEO_METADATA_VIDEO_FILE_CRC_NO_ENTRY)

        return vmd2

    def _initialize_md3(self):
        """
        Initialize metadta section 3 array with default values

        Returns:
        --------
        Initialized numpy array of md3 dtype
        """

        md3 = np.zeros(1, create_md3_dtype())

        md3['recording_time_offset'] = METADATA_RECORDING_TIME_OFFSET_NO_ENTRY
        md3['DST_start_time'] = METADATA_DST_START_TIME_NO_ENTRY
        md3['DST_end_time'] = METADATA_DST_END_TIME_NO_ENTRY
        md3['GMT_offset'] = GMT_OFFSET_NO_ENTRY
        md3['subject_name_1'] = ""
        md3['subject_name_2'] = ""
        md3['recording_location'] = ""

        return md3

    def write_mef_ts_segment_metadata(self, channel, segment_n,
                                      password_1, password_2,
                                      start_time, end_time,
                                      section_2_dict, section_3_dict):
        """
        Writes new time series metadata in the specified segment

        Parameters:
        -----------
        channel - channel name
        segment_n - segment number
        password_1 - level 1 password
        password_2 - level 2 password
        start_time - start time
        end time - end time
        section_2_dict - dictionary with user specified section_2 fileds
        section_3_dict - dictionary with user specified section_3 fileds
        """

        segment_path = (self.session_path+channel+'.timd/'
                        + channel+'-'+str(segment_n).zfill(6)+'.segd/')

        tmet_path = segment_path+channel+'-'+str(segment_n).zfill(6)+'.tmet'

        if os.path.exists(tmet_path):
            raise RuntimeError('Metadata file '+tmet_path+' already exists!')

        os.makedirs(segment_path, exist_ok=True)

        section_2 = self._initialize_tmd2()
        for key, value in section_2_dict.items():
            section_2[key] = value

        section_3 = self._initialize_md3()
        for key, value in section_3_dict.items():
            section_3[key] = value

        write_mef_ts_metadata(segment_path,
                              password_1,
                              password_2,
                              start_time,
                              end_time,
                              section_2,
                              section_3)

    def write_mef_ts_segment_data(self, channel, segment_n,
                                  password_1, password_2,
                                  samps_per_mef_block,
                                  data):
        """
        Writes new time series metadata in the specified segment

        Parameters:
        -----------
        channel - channel name
        segment_n - segment number
        password_1 - level 1 password
        password_2 - level 2 password
        samps_per_mef_block - number of samples per mef block
        end data - numpy array of type int32
        """

        segment_path = (self.session_path+channel+'.timd/'
                        + channel+'-'+str(segment_n).zfill(6)+'.segd/')

        tdat_path = segment_path+channel+'-'+str(segment_n).zfill(6)+'.tdat'

        if not os.path.exists(segment_path
                              + channel+'-'+str(segment_n).zfill(6)+'.tmet'):
            raise RuntimeError('Please write the metadata file first!')

        if os.path.exists(tdat_path):
            raise RuntimeError('Data file '+tdat_path+' already exists!')

        write_mef_ts_data_and_indices(segment_path,
                                      password_1,
                                      password_2,
                                      samps_per_mef_block,
                                      data,
                                      0)  # lossy compression flag - not used

    def append_mef_ts_segment_data(self, channel, segment_n,
                                   password_1, password_2,
                                   start_time, end_time,
                                   samps_per_mef_block,
                                   data):
        """
        Appends new time series metadata in the specified segment

        Parameters:
        -----------
        channel - channel name
        segment_n - segment number
        password_1 - level 1 password
        password_2 - level 2 password
        start_time - start time
        end_time - end time
        samps_per_mef_block - number of samples per mef block
        end data - numpy array of type int32
        """

        segment_path = (self.session_path+channel+'.timd/'
                        + channel+'-'+str(segment_n).zfill(6)+'.segd/')

        tdat_path = segment_path+channel+'-'+str(segment_n).zfill(6)+'.tdat'

        if not os.path.exists(tdat_path):
            raise RuntimeError('Data file does not exist!')

        append_ts_data_and_indices(segment_path,
                                   password_1,
                                   password_2,
                                   start_time,
                                   end_time,
                                   samps_per_mef_block,
                                   data,
                                   0)  # lossy compression flag - not used

    def write_mef_v_segment_metadata(self, channel, segment_n,
                                     password_1, password_2,
                                     start_time, end_time,
                                     section_2_dict, section_3_dict):
        """
        Writes new video metadata in the specified segment

        Parameters:
        -----------
        channel - channel name
        segment_n - segment number
        password_1 - level 1 password
        password_2 - level 2 password
        start_time - start time
        end time - end time
        section_2_dict - dictionary with user specified section_2 fileds
        section_3_dict - dictionary with user specified section_3 fileds
        """

        segment_path = (self.session_path+channel+'.vidd/'
                        + channel+'-'+str(segment_n).zfill(6)+'.segd/')

        vmet_path = segment_path+channel+'-'+str(segment_n).zfill(6)+'.vmet'

        if os.path.exists(vmet_path):
            raise RuntimeError('Segment '+vmet_path+' already exists!')

        os.makedirs(segment_path, exist_ok=True)

        section_2 = self._initialize_vmd2()
        for key, value in section_2_dict.items():
            section_2[key] = value

        section_3 = self._initialize_md3()
        for key, value in section_3_dict.items():
            section_3[key] = value

        write_mef_v_metadata(segment_path,
                             password_1,
                             password_2,
                             start_time,
                             end_time,
                             section_2,
                             section_3)

    def write_mef_v_segment_indices(self, channel, segment_n,
                                    password_1, password_2,
                                    start_time, end_time,
                                    index_entries):

        """
        Writes new video indices in the specified segment

        Parameters:
        -----------
        channel - channel name
        segment_n - segment number
        password_1 - level 1 password
        password_2 - level 2 password
        start_time - start time
        end time - end time
        index_entries - numpy array with vi_dtype
        """

        segment_path = (self.session_path+channel+'.vidd/'
                        + channel+'-'+str(segment_n).zfill(6)+'.segd/')

        write_mef_v_indices(segment_path,
                            password_1,
                            password_2,
                            start_time,
                            end_time,
                            index_entries)

    def write_mef_records(self, password_1, password_2,
                          start_time, end_time,
                          time_offset, records_list,
                          channel_type='ts',
                          channel=None, segment_n=None):

        """
        TODO: - this funciton should use numpy arrays in the future
        Writes new records on session level. If channel is specified, the
        the records are written at channel level. If segment_n is sepcified,
        the recordes are written at segment level.

        Parameters:
        -----------
        password_1 - level 1 password
        password_2 - level 2 password
        start_time - start time
        end time - end time
        time_offset - time offset for records
        record_list - python list with record dictionaries
        channel_type - 'ts' (time series) or 'v' (video), default='ts'
        channel - channel name
        segment_n - segment number
        """

        if channel is None and segment_n is not None:
            raise ValueError('Channel has to be set if segment is set')

        dir_path = self.session_path
        if channel is not None:
            if channel_type == 'ts':
                dir_path += channel+'.timd/'
            elif channel_type == 'v':
                dir_path += channel+'.vidd/'
            else:
                raise ValueError('Invalid channel_type, allowed options are:'
                                 '"ts" or "v"')
        if segment_n is not None:
            dir_path += channel+'-'+str(segment_n).zfill(6)+'.segd/'

        write_mef_data_records(dir_path,
                               password_1,
                               password_2,
                               start_time,
                               end_time,
                               time_offset,
                               records_list)

    # ----- Data reading functions -----
    def get_channel_toc(self, channel):

        """
        Returns discontinuities accross segments

        Parameters:
        -----------
        channel_md - channel metadata dictionart\n

        Returns:
        --------
        TOC - array with
              - [0,:] = discontinuity flags
              - [1,:] = discont lengths
              - [2,:] = start samples
              - [3,:] = start uutc times
        """

        channel_md = self.session_md['time_series_channels'][channel]
        toc = np.empty([4, 0], dtype='int64')
        fsamp = channel_md['section_2']['sampling_frequency']

        # Sort the segments to eliminate dictionary randomness
        segs = list(channel_md['segments'].keys())
        segs.sort()
        for segment_name in segs:
            seg_toc = channel_md['segments'][segment_name]['TOC']

            # Join into channel TOC
            toc = np.concatenate([toc, seg_toc], axis=1)

        # Once we have all segments get lenghts (differnces between segments)
        toc[1, 1::] = (((np.diff(toc[3, :]) / 1e6)
                        - (np.diff(toc[2, :]) / fsamp)) * 1e6)

        return toc

    def read_ts_channels_sample(self, channel_map, sample_map, process_n=None):
        """
        Reads desired channels in desired sample segment

        Parameters:
        -----------
        channel_map - channel or list of channels to be read
        sample_map - list of [start,stop] samples to be loaded that correspond
            to channel_map. if there is only one entry the same range is
            applied to all channels
        process_n - how many processes use for reading (defualt None)

        Returns:
        --------
        data - numpy array [channels,samples]
        """

        data_list = []

        if not isinstance(channel_map, (list, np.ndarray, str)):
            raise TypeError('Channel map has to be list, array or str')

        if isinstance(channel_map, str):
            is_chan_str = True
            channel_map = [channel_map]
        else:
            is_chan_str = False

        if not isinstance(sample_map[0], (list, np.ndarray)):
            sample_map = [sample_map]

        if len(sample_map) == 1:
            sample_map = sample_map*len(channel_map)

        if len(sample_map) != len(channel_map):
            raise RuntimeError('Length of sample map is not equivalent'
                               'to the length of channel map')

        if process_n is not None and not isinstance(process_n, int):
            raise RuntimeError('Process_n argument must be None or int')

        if process_n is not None:
            mp = Pool(process_n)

            iterator = []
            for channel, sample_ss in zip(channel_map, sample_map):

                iterator.append([self._get_channel_md(channel),
                                 sample_ss[0], sample_ss[1]])

            data_list = mp.map(self._arg_merger, iterator)
            mp.terminate()
            if is_chan_str:
                return data_list[0]
            else:
                return data_list

        for channel, sample_ss in zip(channel_map, sample_map):
            data = read_mef_ts_data(self._get_channel_md(channel),
                                    sample_ss[0], sample_ss[1])
            data_list.append(data)

        if is_chan_str:
            return data_list[0]
        else:
            return data_list

    def read_ts_channels_uutc(self, channel_map, uutc_map, process_n=None):
        """
        Reads desired channels in desired time segment. Missing data at
        discontinuities are filled with NaNs.

        Parameters:
        -----------
        channel_map - channel or list of channels to be read
        uutc_map - list of [start,stop] uutc times to be loaded that correspond
           to channel_map. if there is only one entry the same range is applied
           to all channels
        process_n - how many processes use for reading (defualt None)

        Returns:
        --------
        data - numpy array [channels,samples]
        """

        data_list = []

        if not isinstance(channel_map, (list, np.ndarray, str)):
            raise TypeError('Channel map has to be list, array or str')

        if isinstance(channel_map, str):
            is_chan_str = True
            channel_map = [channel_map]
        else:
            is_chan_str = False

        if not isinstance(uutc_map[0], (list, np.ndarray)):
            uutc_map = [uutc_map]

        if len(uutc_map) == 1:
            uutc_map = uutc_map*len(channel_map)

        if len(uutc_map) != len(uutc_map):
            raise RuntimeError('Length of uutc map is not equivalent'
                               'to the length of channel map')

        if process_n is not None and not isinstance(process_n, int):
            raise RuntimeError('Process_n argument must be Nnoe or int')

        if process_n is not None:
            mp = Pool(process_n)

            iterator = []
            for channel, sample_ss in zip(channel_map, uutc_map):
                iterator.append([self._get_channel_md(channel),
                                 sample_ss[0], sample_ss[1], True])

            data_list = mp.map(self._arg_merger, iterator)
            mp.terminate()
            if is_chan_str:
                return data_list[0]
            else:
                return data_list

        for channel, uutc_ss in zip(channel_map, uutc_map):
            data = read_mef_ts_data(self._get_channel_md(channel),
                                    uutc_ss[0], uutc_ss[1], True)
            data_list.append(data)

        if is_chan_str:
            return data_list[0]
        else:
            return data_list

    def read_ts_channel_basic_info(self):
        """
        Reads session time series channel names

        Parameters:
        -----------
        session_path - path to mef3 session
        password - mef3 data password

        Returns:
        --------
        channel_list
        """

        channel_list = list(self.session_md['time_series_channels'].keys())
        channel_list.sort()

        channel_infos = []
        for channel in channel_list:

            channel_md = self.session_md['time_series_channels'][channel]
            channel_md_spec = channel_md['channel_specific_metadata']
            channel_md_s2 = channel_md['section_2']

            fsamp = channel_md_s2['sampling_frequency']
            nsamp = channel_md_s2['number_of_samples']
            ufact = channel_md_s2['units_conversion_factor']
            unit = channel_md_s2['units_description']
            start_time = channel_md_spec['earliest_start_time']
            end_time = channel_md_spec['latest_end_time']
            ch_desc = channel_md_s2['channel_description']

            channel_infos.append({'name': channel, 'fsamp': fsamp,
                                  'nsamp': nsamp, 'ufact': ufact, 'unit': unit,
                                  'start_time':  start_time,
                                  'end_time': end_time,
                                  'channel_description': ch_desc})

        return channel_infos

    # ----- Session operations -----
    def annonymize_session(self, password_1, password_2,
                           new_name=None, new_id=None):
        """
        Anonymize mef session

        Parameters:
        -----------
        session_path - path to the session.
        password_1 - session password level 1
        password_2 - session password level 2
        new_name - new first name for the subject (default = None)
        new_id - new subject id (default = None)

        Returns:
        --------
        0 - on success
        """

        if self.password == password_1:
            raise ValueError("Password provided for opening the session was \
                              level 1, please provide password level 2")

        # Get metadata files and create a list matching the session md
        md_file_list = []
        for root, _, files in os.walk(self.session_path):
            if root.endswith(".timd"):
                channel = root[root.rindex('/')+1:-5]
                channel_md = self.session_md['time_series_channels'][channel]
            elif root.endswith(".vidd"):
                channel = root[root.rindex('/')+1:-5]
                channel_md = self.session_md['video_channels'][channel]
            elif root.endswith(".segd"):
                segment = root[root.rindex('/')+1:-5]
                segment_md = channel_md['segments'][segment]

            for file in files:
                if file.endswith(".tmet") or file.endswith(".vmet"):
                    md_file_list.append([segment_md, root,
                                         os.path.join(root, file)])

        # Run through the list, modify section_3 and rewrite the files
        for seg_md, seg_path, seg_md_file in md_file_list:
            section_2 = seg_md['section_2']
            section_3 = seg_md['section_3']

            seg_start = seg_md['universal_headers']['metadata']['start_time']
            seg_stop = seg_md['universal_headers']['metadata']['end_time']

            if new_name is None:
                section_3.pop('subject_name_1')
            else:
                section_3['subject_name_1'] = new_name
            section_3.pop('subject_name_2')
            if new_id is None:
                section_3.pop('subject_ID')
            else:
                section_3['subject_ID'] = new_id

            # Remove 'not entered fields' and modify no filter
            items = list(section_3.items())
            for item in items:
                if item[1] is None:
                    section_3.pop(item[0])

            items = list(section_2.items())
            for item in items:
                if item[1] is None:
                    section_2.pop(item[0])
                if item[1] in ('no low frequency filter',
                               'no high frequency filter',
                               'no notch filter'):
                    section_2[item[0]] = 0

            os.remove(seg_md_file)

            if seg_md_file.endswith('.tmet'):
                write_mef_ts_metadata(seg_path,
                                      password_1,
                                      password_2,
                                      seg_start,
                                      seg_stop,
                                      section_2,
                                      section_3)
            elif seg_md_file.endswith('.vmet'):
                write_mef_v_metadata(seg_path,
                                     password_1,
                                     password_2,
                                     seg_start,
                                     seg_stop,
                                     section_2,
                                     section_3)

        # Reload the session metadata
        self.reload()

        return 0

    def detect_corrupt_data(self, repair=False):
        """
        Detects corrupt data

        Parameters:
        -----------
        repair(bool) - whether to try to repair data (default=False)

        Returns:
        --------
        0 - on success
        """

        tsd = self.session_md['time_series_channels']
        channels = list(tsd)
        channels.sort()

        # ----- Check time indices entries -----

        for channel in channels:

            # Get total number of blocks in the channel
            total_analyzed_blocks = 0
            segments = list(tsd[channel]['segments'])
            segments.sort()

            for segment in segments:
                idcs = tsd[channel]['segments'][segment]['indices']
                total_analyzed_blocks += len(idcs)

            current_total_block = 0
            for segment in segments:
                idcs = tsd[channel]['segments'][segment]['indices']

                path_to_data = (self.session_path + '/'
                                + channel + '.timd/'
                                + segment + '.segd/'
                                + segment + '.tdat')

                # Open the data file
                f_dat = open(path_to_data, 'r+b')
                f_dat.seek(1024)  # Skippnig the UH

                if repair:
                    orig_idx_file = path_to_data[:-4] + 'tidx'
                    bup_idx_file = path_to_data[:-4] + 'tidx_bup'
                    orig_dat_file = path_to_data[:-4] + 'tdat'
                    bup_dat_file = path_to_data[:-4] + 'tdat_bup'

                # Check the file offset and block bytres of the previous block
                for i, idx in enumerate(idcs):
                    current_total_block += 1

                    header_bytes = f_dat.read(304)  # Read block header
                    number_of_samples = struct.unpack('I',
                                                      header_bytes[32:36])[0]
                    block_bytes = struct.unpack('I',
                                                header_bytes[36:40])[0]
                    start_time = struct.unpack('q',
                                               header_bytes[40:48])[0]

                    if idx['block_bytes'] != block_bytes:
                        print("Block", i, "/", len(idcs),
                              "(", current_total_block, "/",
                              total_analyzed_blocks,
                              ")", "segment", segment, "channel", channel,
                              "has different block_bytes than block header:",
                              idx['block_bytes'], " X ", block_bytes)

                        if repair:

                            # Index is zeroed out
                            if idx['block_bytes'] == 0:

                                if not os.path.exists(bup_idx_file):
                                    shutil.copyfile(orig_idx_file,
                                                    bup_idx_file)

                                f_idx = open(orig_idx_file, 'r+b')

                                # Skip the UH
                                f_idx.seek(1024)

                                # Seek to the index entry
                                f_idx.seek(i * 56, 1)

                                # Copy the file offset
                                f_idx.write(struct.pack('q', f_dat.tell()-304))

                                # Copy the start time
                                f_idx.write(struct.pack('q', start_time))

                                # Calculate/insert the start sample (in seg)
                                if i == 0:
                                    f_idx.write(struct.pack('q', 0))
                                else:
                                    n_smples = (idcs[i-1]['start_sample']
                                                + idcs[i-1]
                                                    ['number_of_samples'])
                                    f_idx.write(struct.pack('q', n_smples))

                                # Copy number of samples
                                f_idx.write(struct.pack('I',
                                                        number_of_samples))

                                # Copy blok bytes
                                f_idx.write(struct.pack('I',
                                                        block_bytes))

                                # Insert RED_NAN to the sample values
                                f_idx.write(bytearray.fromhex("000000"))
                                f_idx.write(bytearray.fromhex("000000"))

                                # Close the file
                                f_idx.close()

                            # RED_block header is zeroed out - invalidate seg
                            else:
                                path_to_segment = (self.session_path + '/'
                                                   + channel + '.timd/'
                                                   + segment + '.segd')
                                warn_str = ('Data cannot be recovered,'
                                            + ' invalidating segment '
                                            + path_to_segment)
                                warnings.warn(warn_str, RuntimeWarning)
                                os.rename(path_to_segment,
                                          path_to_segment+'_corrupt')
                                break

                    # Skip the block bytes - move to the next block
                    f_dat.seek(block_bytes - 304, 1)

                # Check that we are not out of file (diff from metadata)
                f_size = os.path.getsize(path_to_data)
                f_size_md = f_dat.tell()
                if f_size > f_size_md:

                    if repair:
                        warn_str = ('Data file larger than metadata info,'
                                    + 'cutting file ' + path_to_data)
                        warnings.warn(warn_str, RuntimeWarning)

                        if not os.path.exists(bup_dat_file):
                            shutil.copyfile(orig_dat_file, bup_dat_file)

                        f_dat.seek(0)
                        whole_file = f_dat.read(f_size_md)
                        f_dat.close()

                        f_dat = open(path_to_data, 'wb')
                        f_dat.write(whole_file)

                    else:
                        warn_str = ("Data file larger than metadata info"
                                    + path_to_data)
                        warnings.warn(warn_str, RuntimeWarning)

                f_dat.close()

        # Reload the session metadata
        self.reload()
