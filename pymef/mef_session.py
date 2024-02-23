#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -----------------------------------------------------------------------------
# Copyright (c) Jan Cimbalnik, Matt Stead, Ben Brinkmann, and Dan Crepeau.
# All Rights Reserved.
# Distributed under the (new) BSD License. See LICENSE.txt for more info.
# -----------------------------------------------------------------------------

# Standard library imports
import os
import struct
import shutil
import warnings
from multiprocessing import Pool
from pathlib import Path

# Third party imports
import numpy as np

# Local imports
from pymef.mef_file.pymef3_file import (read_mef_session_metadata,
                                        read_mef_ts_data,
                                        clean_mef_session_metadata,
                                        write_mef_ts_metadata,
                                        write_mef_v_metadata,
                                        write_mef_ts_data_and_indices,
                                        append_ts_data_and_indices,
                                        write_mef_v_indices,
                                        write_mef_data_records,
                                        create_rh_dtype,
                                        create_note_dtype,
                                        create_sylg_dtype,
                                        create_edfa_dtype,
                                        create_lntp_dtype,
                                        create_csti_dtype,
                                        create_esti_dtype,
                                        create_seiz_dtype,
                                        create_seiz_ch_dtype,
                                        create_curs_dtype,
                                        create_epoc_dtype,
                                        create_tmd2_dtype,
                                        create_vmd2_dtype,
                                        create_md3_dtype,
                                        check_mef_password)
from pymef.mef_constants import (
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
            TIME_SERIES_METADATA_MAXIMUM_BLOCK_BYTES_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY,
            TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY,
            TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY,
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCKS_NO_ENTRY,

            VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY,
            VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY,
            VIDEO_METADATA_FRAME_RATE_NO_ENTRY,
            VIDEO_METADATA_NUMBER_OF_CLIPS_NO_ENTRY,
            VIDEO_METADATA_MAXIMUM_CLIP_BYTES_NO_ENTRY,
            VIDEO_METADATA_VIDEO_FILE_CRC_NO_ENTRY,

            METADATA_RECORDING_TIME_OFFSET_NO_ENTRY,
            METADATA_DST_START_TIME_NO_ENTRY,
            METADATA_DST_END_TIME_NO_ENTRY,
            GMT_OFFSET_NO_ENTRY
        )

MEF_FILE_EXTENSIONS = ['mefd', 'segd',
                       'rdat', 'ridx',
                       'vidd', 'vmet', 'vidx',
                       'timd', 'tmet', 'tdat', 'tidx']


class MefSession():
    """
    Basic object for operations with mef sessions.

    Parameters
    ----------
    session_path: str
        path to mef session
    password: str
        password for mef session
    read_metadata: bool
        whether to read metadata (default=True)
    new_session: bool
        whether this is a new session for writing (default=False)
    check_all_passwords: bool
        check all files or just the first one encoutered(default=True)
    """

    def __init__(self, session_path, password, read_metadata=True,
                 new_session=False, check_all_passwords=True):

        if not session_path.endswith('/'):
            session_path += '/'

        if '.mefd' != session_path[-6:-1]:
            raise ValueError("Session path must end with .mefd suffix!")

        self.path = session_path
        self.password = password

        if new_session:
            os.makedirs(session_path)
            self.session_md = None
            return

        if not os.path.exists(session_path):
            raise FileNotFoundError(session_path+' does not exist!')

        # Check if path exists
        if not os.path.exists(session_path):
            raise FileNotFoundError(session_path+' does not exist!')
        self._check_password(check_all_passwords)

        if read_metadata:
            self.session_md = read_mef_session_metadata(session_path,
                                                        password)
        else:
            self.session_md = None

        return

    def __repr__(self):
        return (f"MefSession(session_path={self.path}"
                f", password={self.password})")

    def __str__(self):
        if self.session_md is None:
            return f"Mef session: {self.path}"
        else:
            md = self.session_md['session_specific_metadata']
            ts_ch_n = md['number_of_time_series_channels'][0]
            v_ch_n = md['number_of_video_channels'][0]
            descr_str = (f"Mef session: {self.path}\n"
                         f"Start time: {md['earliest_start_time'][0]}\n"
                         f"End time: {md['latest_end_time'][0]}\n"
                         f"Number of time series channels: {ts_ch_n}\n"
                         f"Number of video channels: {v_ch_n}\n")
        return descr_str

    def __len__(self):
        if self.session_md is None:
            return None
        else:
            md = self.session_md['session_specific_metadata']
            return md['latest_end_time'][0] - md['earliest_start_time'][0]

    # ----- Helper functions -----
    def _check_password(self, check_all=True):
        """
        Checks provided password on all files in the session.

        Parameters
        ----------
        check_all: bool
            whether to check all files or just the first one encoutered,
            checking just one file significantly improves speed.

        Returns
        -------
        result: object
            None on success
        """
        mef_files = []
        for path, subdirs, files in os.walk(self.path):
            for name in files:
                if any([name.endswith(ext) for ext in MEF_FILE_EXTENSIONS]):
                    mef_files.append(os.path.join(path, name))

        results = np.zeros(len(mef_files))
        for i, mef_file in enumerate(mef_files):
            results[i] = check_mef_password(mef_file, self.password)
            if not check_all:
                break

        if np.any(results < 0):
            raise RuntimeError('MEF password is invalid')

        return

    def _get_channel_md(self, channel):
        """
        Paramters
        ---------
        channel: str
            required channel

        Returns
        -------
        channel_md: dict
            Channel metadata structure
        """
        ts_chs = self.session_md['time_series_channels']
        return ts_chs[channel]['channel_specific_metadata']

    def _arg_merger(self, args):
        return read_mef_ts_data(*args)

    def reload(self):
        self.close()
        self.session_md = read_mef_session_metadata(self.path,
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

        Returns
        -------
        init_tmd2: np.array
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
        tmd2['maximum_block_bytes'] = (
            TIME_SERIES_METADATA_MAXIMUM_BLOCK_BYTES_NO_ENTRY)
        tmd2['maximum_block_samples'] = (
            TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY)
        tmd2['maximum_difference_bytes'] = (
            TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_NO_ENTRY)
        tmd2['block_interval'] = (
            TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY)
        tmd2['number_of_discontinuities'] = (
            TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY)
        tmd2['maximum_contiguous_blocks'] = (
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCKS_NO_ENTRY)
        tmd2['maximum_contiguous_block_bytes'] = (
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY)
        tmd2['maximum_contiguous_samples'] = (
            TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_NO_ENTRY)

        return tmd2

    def _initialize_vmd2(self):
        """
        Initialize video metadata section 2 array with default values

        Returns
        -------
        init_vmd2: np.array
            Initialized numpy array of video md2 dtype
        """

        vmd2 = np.zeros(1, create_vmd2_dtype())

        vmd2['recording_duration'] = (
            METADATA_RECORDING_DURATION_NO_ENTRY)
        vmd2['horizontal_resolution'] = (
            VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY)
        vmd2['vertical_resolution'] = (
            VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY)
        vmd2['frame_rate'] = (
            VIDEO_METADATA_FRAME_RATE_NO_ENTRY)
        vmd2['number_of_clips'] = (
            VIDEO_METADATA_NUMBER_OF_CLIPS_NO_ENTRY)
        vmd2['maximum_clip_bytes'] = (
            VIDEO_METADATA_MAXIMUM_CLIP_BYTES_NO_ENTRY)
        vmd2['video_file_CRC'] = (
            VIDEO_METADATA_VIDEO_FILE_CRC_NO_ENTRY)

        return vmd2

    def _initialize_md3(self):
        """
        Initialize metadta section 3 array with default values

        Returns
        -------
        init_md3: np.array
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

    def _create_np_record(self, record):
        """
        Create record dictionary with numpy arrays from python dictionary.

        Parameters
        ----------
        record: dict
            Python dictionary with record entries

        Returns
        -------
        record_dict: dict
            Dictionary with numpy arrays record header, body, subbody
        """

        if not isinstance(record, dict):
            raise ValueError("Record is not a dictionary")

        record_type = record.get('type')

        if record_type is None:
            raise ValueError("Record does not contain 'type' field")

        hdr_dtype = create_rh_dtype()
        hdr_arr = np.zeros(1, hdr_dtype)

        if record_type == 'Note':
            hdr_arr['type_string'] = b'Note'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            if 'text' in record.keys():
                record_str = record['text']+'\0'
                body_dtype = create_note_dtype(len(record_str))
                body_arr = np.zeros(1, body_dtype)
                body_arr['text'] = record_str.encode('utf8') + b'\x00'

        elif record_type == 'SyLg':
            hdr_arr['type_string'] = b'SyLg'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            if 'text' in record.keys():
                record_str = record['text']+'\0'
                body_dtype = create_sylg_dtype(len(record_str))
                body_arr = np.zeros(1, body_dtype)
                body_arr['text'] = record_str.encode('utf8') + b'\x00'

        elif record_type == 'EDFA':
            hdr_arr['type_string'] = b'EDFA'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            if 'text' in record.keys():
                record_str = record['text']+'\0'
                body_dtype = create_edfa_dtype(len(record_str))
                body_arr = np.zeros(1, body_dtype)
                body_arr['text'] = record_str.encode('utf8') + b'\x00'
            if 'duration' in record.keys():
                body_arr['duration'] = record['duration']

        elif record_type == 'LNTP':
            hdr_arr['type_string'] = b'LNTP'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            template = record['template']
            body_dtype = create_lntp_dtype(len(template))
            body_arr = np.zeros(1, body_dtype)
            body_arr['length'] = len(template)
            body_arr['template'] = template

        elif record_type == 'CSti':
            hdr_arr['type_string'] = b'CSti'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            body_dtype = create_csti_dtype()
            d_type_keys = [x[0] for x in body_dtype.descr]
            body_arr = np.zeros(1, body_dtype)
            for key in record.keys():
                if key in ['type', 'time']:
                    continue
                if isinstance(record[key], str):
                    body_arr[key] = record[key].encode('utf8') + b'\x00'
                else:
                    body_arr[key] = record[key]

        elif record_type == 'ESti':
            hdr_arr['type_string'] = b'ESti'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            body_dtype = create_esti_dtype()
            d_type_keys = [x[0] for x in body_dtype.descr]
            body_arr = np.zeros(1, body_dtype)
            for key in record.keys():
                if key in ['type', 'time']:
                    continue
                if isinstance(record[key], str):
                    body_arr[key] = record[key].encode('utf8') + b'\x00'
                else:
                    body_arr[key] = record[key]

        elif record_type == 'Seiz':
            hdr_arr['type_string'] = b'Seiz'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            body_dtype = create_seiz_dtype()
            body_arr = np.zeros(1, body_dtype)
            for key in record.keys():
                if key in ['type', 'time']:
                    continue
                if key == 'channels':
                    channels = record['channels']
                    subbody_dtype = create_seiz_ch_dtype()
                    subbody_arr = np.zeros(len(record['channels']),
                                           subbody_dtype)
                    subbody_arr['name'] = [x['name'] for x in channels]
                    subbody_arr['onset'] = [x['onset'] for x in channels]
                    subbody_arr['offset'] = [x['offset'] for x in channels]
                else:
                    if isinstance(record[key], str):
                        body_arr[key] = record[key].encode('utf8') + b'\x00'
                    else:
                        body_arr[key] = record[key]

        elif record_type == 'Curs':
            hdr_arr['type_string'] = b'Curs'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            body_dtype = create_curs_dtype()
            d_type_keys = [x[0] for x in body_dtype.descr]
            body_arr = np.zeros(1, body_dtype)
            for key in record.keys():
                if key in ['type', 'time']:
                    continue
                if isinstance(record[key], str):
                    body_arr[key] = record[key].encode('utf8') + b'\x00'
                else:
                    body_arr[key] = record[key]

        elif record_type == 'Epoc':
            hdr_arr['type_string'] = b'Epoc'
            if 'time' in record.keys():
                hdr_arr['time'] = record['time']

            body_dtype = create_epoc_dtype()
            d_type_keys = [x[0] for x in body_dtype.descr]
            body_arr = np.zeros(1, body_dtype)
            for key in record.keys():
                if key in ['type', 'time']:
                    continue
                if isinstance(record[key], str):
                    body_arr[key] = record[key].encode('utf8') + b'\x00'
                else:
                    body_arr[key] = record[key]

        else:
            raise ValueError("Unrecognized record type:'%s'" % record_type)

        record_np_dict = {'record_header': hdr_arr,
                          'record_body': body_arr}

        if 'subbody_arr' in locals():
            record_np_dict['record_subbody'] = subbody_arr

        return record_np_dict

    def write_mef_ts_segment_metadata(self, channel, segment_n,
                                      password_1, password_2,
                                      start_time, end_time,
                                      section_2_dict, section_3_dict):
        """
        Writes new time series metadata in the specified segment

        Parameters
        ----------
        channel: str
            Channel name
        segment_n: int
            Segment number
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        start_time: int
            Start time
        end time: int
            End time
        section_2_dict: dict
            Dictionary with user specified section_2 fileds
        section_3_dict: dict
            Dictionary with user specified section_3 fileds
        """

        segment_path = (self.path+channel+'.timd/'
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
        Writes new time series data in the specified segment

        Parameters
        ----------
        channel: str
            Channel name
        segment_n: int
            Segment number
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        samps_per_mef_block: int
            Number of samples per mef block
        data: np.array
            1-D numpy array of type int32
        """

        segment_path = (self.path+channel+'.timd/'
                        + channel+'-'+str(segment_n).zfill(6)+'.segd/')

        tdat_path = segment_path+channel+'-'+str(segment_n).zfill(6)+'.tdat'

        if not os.path.exists(segment_path
                              + channel+'-'+str(segment_n).zfill(6)+'.tmet'):
            raise RuntimeError('Please write the metadata file first!')

        if os.path.exists(tdat_path):
            raise RuntimeError(f"Data file '{tdat_path}' already exists!")

        # lossy compression flag - not used
        write_mef_ts_data_and_indices(segment_path,
                                      password_1,
                                      password_2,
                                      samps_per_mef_block,
                                      data,
                                      0)

    def append_mef_ts_segment_data(self, channel, segment_n,
                                   password_1, password_2,
                                   start_time, end_time,
                                   samps_per_mef_block,
                                   data, discontinuity_flag=False):
        """
        Appends new time series metadata in the specified segment

        Parameters
        ----------
        channel: str
            Channel name
        segment_n: int
            Segment number
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        start_time: int
            Start time of the appended data
        end_time:
            End time of the appended data
        samps_per_mef_block: int
            Number of samples per mef block
        data: np.array
            1-D numpy array of type int32
        discontinuity_flag: bool
            Discontinuity flag for appended data.
        """

        segment_path = (self.path+channel+'.timd/'
                        + channel+'-'+str(segment_n).zfill(6)+'.segd/')

        tdat_path = segment_path+channel+'-'+str(segment_n).zfill(6)+'.tdat'

        if not os.path.exists(tdat_path):
            raise RuntimeError(f"Data file '{tdat_path}' does not exist!")

        # TODO - check for bool in C code
        if discontinuity_flag is True:
            discontinuity_flag = 1
        else:
            discontinuity_flag = 0

        # lossy compression flag - not used
        append_ts_data_and_indices(segment_path,
                                   password_1,
                                   password_2,
                                   start_time,
                                   end_time,
                                   samps_per_mef_block,
                                   data,
                                   discontinuity_flag)

    def write_mef_v_segment_metadata(self, channel, segment_n,
                                     password_1, password_2,
                                     start_time, end_time,
                                     section_2_dict, section_3_dict):
        """
        Writes new video metadata in the specified segment

        Parameters
        ----------
        channel: str
            Channel name
        segment_n: int
            Segment number
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        start_time: int
            Start time of the data
        end time: int
            End time of the data
        section_2_dict: dict
            Dictionary with user specified section_2 fileds
        section_3_dict: dict
            Dictionary with user specified section_3 fileds
        """

        segment_path = (self.path+channel+'.vidd/'
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

        Parameters
        ----------
        channel: str
            Channel name
        segment_n: int
            Segment number
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        start_time: int
            Start time of video indices
        end time: int
            End time of video indices
        index_entries: np.array
            Numpy array with vi_dtype
        """

        segment_path = (self.path+channel+'.vidd/'
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
        Writes new records on session level. If channel is specified, the
        the records are written at channel level. If segment_n is sepcified,
        the recordes are written at segment level.

        Parameters
        ----------
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        start_time: int
            Records start time
        end time: int
            Records end time
        time_offset: int
            Time offset for records
        record_list: list
            List with record dictionaries
        channel_type: str
            Channel type: 'ts' (time series) or 'v' (video), default='ts'
        channel: str
            Channel name
        segment_n: int
            Segment number

        Notes
        -----
        Each entry in record list must be a dictionary and contain the field
        "type". The rest of the entries are optional. All times are in uUTC or
        us.
        The following types are recognized:

        Note - simple note:
            - type - "Note"
            - time - record time (int)
            - text - note text (str)

        SyLg - system log:
            - type - "SyLg"
            - time - record time (int)
            - text - system log text

        EDFA - EDF annotation record:
            - type - "EDFA"
            - time - record time (int)
            - duration - duration of the event (int)
            - text - annotation text (str)

        LNTP - Line noise template:
            - type - "LNTP"
            - time - record time (int)
            - length - template length (int)
            - template - 1D numpy array of template itself (np.arry, dtype=int)

        CSti - Cognitive stimulation:
            - type - "CSti"
            - time - record time (int)
            - task_type - type of the task (str)
            - stimulus_duration - duration of the stimulus (int)
            - stimulus_type - type of the stimulus (str)
            - patient_response - response of the patient (str)

        ESti - Electrical stimulation:
            - type - "Esti"
            - time - record time (int)
            - amplitude - stimulus amplitude (float)
            - frequency - frequency of the stimulus (float)
            - pulse_width - pulse width (int)
            - ampunit_code - code of amplitude unit (int)
                - -1 = no entry
                - 0 = unknown
                - 1 = mA
                - 2 = V
            - mode code - code of the stimulation mode (int)
                - -1 = no entry
                - 0 = unknown
                - 1 = current
                - 2 = voltage
            - anode - stimulation anode (str)
            - catode - stimulation catode (str)

        Seiz - Seizure:
            - type - "Seiz"
            - time - record time (int)
            - earliest_onset - earliest uUTC onset of the seizure (int)
            - latest offset - latest uUTC offset of the seizure (int)
            - duration - duration of the seizure (int)
            - number_of_channels - number of seizure onset channels (int)
            - onset_code - code for the onset
                - -1 = no entry
                - 0 - unknown
                - 1 - docal
                - 2 - generalized
                - 3 -propagated
                - 4 - mixed
            - marker_name_1 - name of the marker 1 (str)
            - marker_name_2 - name of the marker 2 (str)
            - annotation - seizure annotation (str)
            - channels - list of dictionaries with channel entries
                - name - name of the channel (str)
                - onset - seizure onset on this channel (int)
                - offset - seizure offset on this channel (int)

        Curs - Cursor record:
            - type - "Curs"
            - time - record time (int)
            - id_record - id of the record (int)
            - trace_timestamp - timestamp of the trace (int8)
            - latency - latency of the record (int)
            - name - name of the cursor record (str)

        Epoc - Epoch record:
            - type - "Epoc"
            - time - record time (int)
            - id_record - id of the record (int)
            - timestamp - start timestamp of the epoch (int8)
            - end_timestamp - end timestamp of the epoch (int8)
            - duration - duration of the epoch
            - type - epoch type (str)
            - text - descriptive text of the epoch (str)

        """

        if channel is None and segment_n is not None:
            raise ValueError('Channel has to be set if segment is set')

        dir_path = self.path
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

        numpy_records_list = []
        for record in records_list:
            numpy_records_list.append(self._create_np_record(record))

        write_mef_data_records(dir_path,
                               password_1,
                               password_2,
                               start_time,
                               end_time,
                               time_offset,
                               numpy_records_list)

    def create_slice_session(self, slice_session_path, slice_start_stop,
                             password_1, password_2, samps_per_mef_block=None,
                             time_unit='uutc'):
        """
        Function to create slice of the mef session (time series only).

        Parameters
        ----------
        slice_session_path: str
            Path to new sliced session (including .mefd suffix)
        slice_start_stop: list
            List containing 2 int with start and stop uutc times of the slice
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        samps_per_mef_block: int
            Number of samples for mef block. Default=channel sampling
            frequency
        """

        if not slice_session_path.endswith('/'):
            slice_session_path += '/'

        if self.session_md is None:
            raise ValueError("Please read the session metadata first.")

        if not any([isinstance(x, (int, np.int32, np.int64)) for x in slice_start_stop]):
            raise ValueError("Start and stop must be integers.")

        channels_dict = self.session_md['time_series_channels']

        segment_n = 0
        for channel, ch_md in channels_dict.items():
            
            toc = self.get_channel_toc(channel)
            
            toc_start_idx = np.where(toc[0] == 1)[0]
            toc_start_times = toc[-1, toc_start_idx]
            toc_stop_times = toc[-1, toc_start_idx[1:]-1] + ((toc[1, toc_start_idx[1:]-1]/5000)*1e6).astype(int)
            toc_stop_times = np.concatenate([toc_stop_times, self.session_md['session_specific_metadata']['latest_end_time']])
            
            toc_slices = np.stack([toc_start_times, toc_stop_times], axis=1)
            toc_slices = toc_slices[((toc_slices[:, 1] > slice_start_stop[0])
                                     & (toc_slices[:, 0] < slice_start_stop[1]))]

            if toc_slices[0, 0] < slice_start_stop[0]:
                toc_slices[0, 0] = slice_start_stop[0]

            if toc_slices[-1, 1] > slice_start_stop[1]:
                toc_slices[-1, 1] = slice_start_stop[1]

            segment_path = (slice_session_path+channel+'.timd/'
                            + channel+'-'+str(segment_n).zfill(6)+'.segd/')
            
            os.makedirs(segment_path, exist_ok=True)

            tmet_path = (segment_path+channel+'-'+str(segment_n).zfill(6)
                         + '.tmet')

            if os.path.exists(tmet_path):
                raise RuntimeError('Metadata file '+tmet_path
                                   + ' already exists!')

            section_2 = ch_md['section_2'].copy()

            if samps_per_mef_block is None:
                spmb = int(section_2['sampling_frequency'][0])
            else:
                spmb = samps_per_mef_block

            is_first_write = True
            
            for ss in toc_slices:

                if is_first_write is True:
                
                    # Zero out the machine generated fields
                    section_2['recording_duration'] = toc_slices[-1, -1] - toc_slices[0, 0]
                    section_2['maximum_native_sample_value'] = 0.0
                    section_2['minimum_native_sample_value'] = 0.0
                    section_2['number_of_blocks'] = 0
                    section_2['maximum_block_bytes'] = 0
                    section_2['maximum_block_samples'] = 0
                    section_2['maximum_difference_bytes'] = 0
                    section_2['block_interval'] = 0
                    section_2['maximum_contiguous_blocks'] = 0
                    section_2['maximum_contiguous_block_bytes'] = 0
                    section_2['maximum_contiguous_samples'] = 0
                    section_2['number_of_samples'] = 0
                
                    section_3 = ch_md['section_3'].copy()

                    write_mef_ts_metadata(segment_path,
                                          password_1,
                                          password_2,
                                          slice_start_stop[0],
                                          slice_start_stop[1],
                                          section_2,
                                          section_3)

                    data = self.read_ts_channels_uutc(channel, [int(x) for x in ss])
                
                    tdat_path = (segment_path+channel+'-'+str(segment_n).zfill(6)
                                 + '.tdat')
                
                    if os.path.exists(tdat_path):
                        raise RuntimeError('Data file '+tdat_path+' already exists!')

                    if np.sum(np.isnan(data)) == len(data):
                        print(f"Signal of channel {channel} in time {ss} is empty. Skipping...")
                        continue
                
                    # lossy compression flag - not used
                    write_mef_ts_data_and_indices(segment_path,
                                                  password_1,
                                                  password_2,
                                                  spmb,
                                                  data.astype('int32'),
                                                  0)
                    
                    is_first_write = False
                        
                else:
                    
                    data = self.read_ts_channels_uutc(channel, [int(x) for x in ss])
                    if np.sum(np.isnan(data)) == len(data):
                        print(f"Signal of channel {channel} in time {ss} is empty. Skipping...")
                        continue
                    append_ts_data_and_indices(segment_path,
                                               password_1,
                                               password_2,
                                               int(ss[0]),
                                               int(ss[1]),
                                               spmb,
                                               data.astype('int32'),
                                               True)

            # Once finished check if we have written any data at all
            if not os.path.exists(tdat_path):
                print(f"No data written for channel {channel}. Deleting...")
                shutil.rmtree(slice_session_path+channel+'.timd/')


    # ----- Data reading functions -----
    def _create_dict_record(self, np_record):
        """
        Create python dictionary from record dictionary with numpy arrays.

        Parameters
        ----------
        np_record: dict
            Dictionary with numpy arrays record header, body, subbody

        Returns
        -------
        record_dict: dict
            Python dictionary with record entries
        """

        record_header = np_record.get('record_header')
        record_body = np_record.get('record_body')
        record_subbody = np_record.get('record_subbody')

        d_type_keys = [x[0] for x in record_body.dtype.descr]

        rec_dict = {}

        if record_header['type_string'] == b'Note':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        elif record_header['type_string'] == b'SyLg':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        elif record_header['type_string'] == b'EDFA':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        elif record_header['type_string'] == b'LNTP':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        elif record_header['type_string'] == b'CSti':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        elif record_header['type_string'] == b'ESti':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        elif record_header['type_string'] == b'Seiz':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

            if record_subbody is not None:
                ch_list = []
                for ch_record in record_subbody:
                    ch_dict = {'name': ch_record['name'].decode('utf-8'),
                               'onset': ch_record['onset'],
                               'offset': ch_record['offset']}
                    ch_list.append(ch_dict)

                rec_dict['channels'] = ch_list

        elif record_header['type_string'] == b'Curs':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        elif record_header['type_string'] == b'Epoc':
            rec_dict['type'] = record_header['type_string'][0].decode('utf-8')
            rec_dict['time'] = record_header['time'][0]

            for key in d_type_keys:
                value = record_body[key][0]
                if isinstance(value, bytes):
                    rec_dict[key] = value.decode('utf-8')
                else:
                    rec_dict[key] = value

        else:
            warn_string = ('Unrecognized record type: '
                           + record_header['type_string'].decode('utf-8'))
            warnings.warn(warn_string, RuntimeWarning)
            rec_dict['type'] = record_header['type_string'].decode('utf-8')

        return rec_dict

    def read_records(self, channel=None, segment_n=None):
        """
        Returns list of dictionaries with MEF records.

        Parameters
        ----------
        channel: str
            Session channel, if not specified, session records will be read
            (default = None)
        segment_n: int
            Segment number, if not specified, channel records will be read
            (default = None)

        Returns
        -------
        record_list: list
            List of dictionaries with record entries
        """

        if channel is not None:
            if channel in self.session_md['time_series_channels'].keys():
                channel_md = self.session_md['time_series_channels'][channel]
            elif channel in self.session_md['video_channels'].keys():
                channel_md = self.session_md['video_channels'][channel]
            else:
                raise ValueError("No channel %s in this session" % channel)

            if segment_n is not None:
                segment = channel+'-'+str(segment_n).zfill(6)
                if segment in channel_md['segments'].keys():
                    segment_md = channel_md['segments'][segment]
                    records_list = segment_md['records_info']['records']
                else:
                    raise ValueError("No segment %s in this session" % segment)
            else:
                records_list = channel_md['records_info']['records']
        else:
            records_list = self.session_md['records_info']['records']

        python_dict_list = []

        for record in records_list:
            python_dict_list.append(self._create_dict_record(record))

        return python_dict_list

    def get_channel_toc(self, channel):
        """
        Returns discontinuities accross segments.

        Parameters
        ----------
        channel: str
            Channel to calculate TOC on

        Returns
        -------
        TOC: np.array
            Array with
              - [0,:] = discontinuity flags
              - [1,:] = n block samples
              - [2,:] = start samples
              - [3,:] = start uutc times
        """

        channel_md = self.session_md['time_series_channels'][channel]
        toc = np.empty([4, 0], dtype='int64')

        # Sort the segments to eliminate dictionary randomness
        segs = list(channel_md['segments'].keys())
        segs.sort()
        for segment_name in segs:
            seg_toc = channel_md['segments'][segment_name]['TOC']

            # Join into channel TOC
            toc = np.concatenate([toc, seg_toc], axis=1)


        return toc

    def read_ts_channels_sample(self, channel_map, sample_map, process_n=None):
        """
        Reads desired channels in desired sample segment

        Parameters
        ----------
        channel_map: str or list
            Channel or list of channels to be read
        sample_map: list
            List of [start, stop] samples to be loaded that correspond
            to channel_map. if there is only one entry the same range is
            applied to all channels
        process_n: int
            How many processes use for reading (default=None)

        Returns
        -------
        data: np.array(dtype=np.float32)
            Numpy array of numpy array objects [channels,samples] or 1D numpy
            array
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

    def read_ts_channels_uutc(self, channel_map, uutc_map, process_n=None,
                              out_nans=True):
        """
        Reads desired channels in desired time segment. Missing data at
        discontinuities are filled with NaNs.

        Parameters
        ----------
        channel_map: str or list
            Channel or list of channels to be read
        uutc_map: list
            List of [start,stop] uutc times to be loaded that correspond
            to channel_map. if there is only one entry the same range is
            applied to all channels
        process_n: int
            How many processes use for reading (defualt = None)
        out_nans: bool
            Whether to return an array of np.nan if the uutc times for
            channel are completely out of start and end times

        Returns
        -------
        data: np.array(dtype=np.float32)
            Numpy array of numpy array objects [channels,samples] or 1D numpy
            array
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
            if out_nans and data is None:
                channel_md = self.session_md['time_series_channels'][channel]
                size = ((np.diff(uutc_ss) / 1e6)[0] * channel_md['section_2']['sampling_frequency'][0])
                data = np.empty(int(size))
                data[:] = np.nan
            data_list.append(data)

        if is_chan_str:
            return data_list[0]
        else:
            return data_list

    def read_ts_channel_basic_info(self):
        """
        Reads session time series channel names

        Returns
        -------
        channel_list: list
            List of dictionaries with information about channels:
                - Sampling frequency
                - Number of samples
                - Units conversion factor
                - Units description
                - Earliest start time
                - Latest end time
                - Channel description
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

        Parameters
        ----------
        password_1: str
            Level 1 password
        password_2: str
            Level 2 password
        new_name: str
            New first name for the subject (default = None)
        new_id: str
            New subject id (default = None)

        Returns
        -------
        result: object
            None on success
        """

        if self.password == password_1:
            raise ValueError("Password provided for opening the session was \
                              level 1, please provide password level 2")

        # Get metadata files and create a list matching the session md
        md_file_list = []
        for root, _, files in os.walk(self.path):
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

        return None

    def detect_corrupt_data(self, repair=False):
        """
        Detects corrupt data

        Parameters
        ----------
        repair: bool
            Whether to try to repair data (default=False)

        Returns
        -------
        result: object
            None on success
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

                path_to_data = (self.path + '/'
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
                                path_to_segment = (self.path + '/'
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

        return None

    def change_channel_name(self, change_dict):

        """
        Chnges name of one or more channels

        Parameters
        ----------
        change_dict: dict
            Dictionary with channels to be changed {old_name: new_name}

        Returns
        -------
        result: object
            None on success
        """

        # Find channel folder
        for fold in os.listdir(self.path):
            ch_name = os.path.splitext(fold)[0]
            if ch_name in change_dict.keys():
                chan_path = self.path + fold
                
                # Rename files first
                for path in Path(chan_path).rglob('*'):
                    if os.path.splitext(path)[1] == '.segd':
                        continue
                    abs_path = str(path.absolute())
                    path, file_name = os.path.split(abs_path)
                    new_file_name = file_name.replace(ch_name,
                                                      change_dict[ch_name])
                    new_abs_path = '/'.join([path, new_file_name])
                    shutil.move(abs_path, new_abs_path)
                
                # Rename segment folders
                for path in Path(chan_path).rglob('*.segd'):
                    abs_path = str(path.absolute())
                    path, fold_name = os.path.split(abs_path)
                    new_fold_name = fold_name.replace(ch_name,
                                                      change_dict[ch_name])
                    new_abs_path = '/'.join([path, new_fold_name])
                    shutil.move(abs_path, new_abs_path)
                
                # Rename channel folder
                path, fold_name = os.path.split(chan_path)
                new_fold_name = fold_name.replace(ch_name,
                                                  change_dict[ch_name])
                new_abs_path = '/'.join([path, new_fold_name])
                shutil.move(chan_path, new_abs_path)

        self.reload()

        return None
