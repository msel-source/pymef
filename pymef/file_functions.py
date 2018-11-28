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
import os
import struct
import shutil
import warnings
from multiprocessing import Pool

import numpy as np

from .mef_file import pymef3_file


def detect_corrupt_data(session_path, password=None, repair=False):
    """
    Detects corrupt data\n

    Parameters:\n
    -----------\n
    session_path - path to the session.\n
    password - session password (default=None)\n
    repair(bool) - whether to try to repair data (default=False)\n

    Returns:\n
    --------\n
    0 - on success
    """

    session_md = pymef3_file.read_mef_session_metadata(session_path, password)

    tsd = session_md['time_series_channels']
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

            path_to_data = (session_path + '/'
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
                number_of_samples = struct.unpack('I', header_bytes[32:36])[0]
                block_bytes = struct.unpack('I', header_bytes[36:40])[0]
                start_time = struct.unpack('q', header_bytes[40:48])[0]

                if idx['block_bytes'] != block_bytes:
                    print("Block", i, "/", len(idcs),
                          "(", current_total_block, "/", total_analyzed_blocks,
                          ")", "segment", segment, "channel", channel,
                          "has different block_bytes than the block header:",
                          idx['block_bytes'], " X ", block_bytes)

                    if repair:

                        # Index is zeroed out
                        if idx['block_bytes'] == 0:

                            if not os.path.exists(bup_idx_file):
                                shutil.copyfile(orig_idx_file, bup_idx_file)

                            f_idx = open(orig_idx_file, 'r+b')

                            # Skip the UH
                            f_idx.seek(1024)

                            # Seek to the index entry
                            f_idx.seek(i * 56, 1)

                            # Copy the file offset
                            f_idx.write(struct.pack('q', f_dat.tell()-304))

                            # Copy the start time
                            f_idx.write(struct.pack('q', start_time))

                            # Calculate / insert the start sample (in segment)
                            if i == 0:
                                f_idx.write(struct.pack('q', 0))
                            else:
                                n_smples = (idcs[i-1]['start_sample']
                                            + idcs[i-1]['number_of_samples'])
                                f_idx.write(struct.pack('q', n_smples))

                            # Copy number of samples
                            f_idx.write(struct.pack('I', number_of_samples))

                            # Copy blok bytes
                            f_idx.write(struct.pack('I', block_bytes))

                            # Insert RED_NAN to the sample values
                            f_idx.write(bytearray.fromhex("000000"))
                            f_idx.write(bytearray.fromhex("000000"))

                            # Close the file
                            f_idx.close()

                        # RED_block header is zeroed out - invalidate seg
                        else:
                            path_to_segment = (session_path + '/'
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
                    warn_str = ('Data file larger than metadata information,'
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
                    warn_str = ("Data file larger than metadata information"
                                + path_to_data)
                    warnings.warn(warn_str, RuntimeWarning)

            f_dat.close()


def annonimize_session(session_path, password_1, password_2,
                       new_name=None, new_id=None):
    """
    Anonimize mef session\n

    Parameters:
    -----------
    session_path - path to the session.\n
    password_1 - session password level 1\n
    password_2 - session password level 2\n
    new_name - new first name for the subject (default = None)\n
    new_id - new subject id (default = None)\n

    Returns:
    --------
    0 - on success
    """

    # Read the session metadata
    session_md = pymef3_file.read_mef_session_metadata(session_path,
                                                       password_2)

    # Get individual metadata files and create a list matching the session md
    md_file_list = []
    for root, _, files in os.walk(session_path):
        if root.endswith(".timd"):
            channel = root[root.rindex('/')+1:-5]
            channel_md = session_md['time_series_channels'][channel]
        elif root.endswith(".vidd"):
            channel = root[root.rindex('/')+1:-5]
            channel_md = session_md['video_channels'][channel]
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
            pymef3_file.write_mef_ts_metadata(seg_path,
                                              password_1,
                                              password_2,
                                              seg_start,
                                              seg_stop,
                                              section_2,
                                              section_3)
        elif seg_md_file.endswith('.vmet'):
            pymef3_file.write_mef_v_metadata(seg_path,
                                             password_1,
                                             password_2,
                                             seg_start,
                                             seg_stop,
                                             section_2,
                                             section_3)

    return 0


def get_toc(channel_md):

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
          - [2,:] = start uutc times
    """

    toc = np.empty([4, 0], dtype='int64')
    fsamp = channel_md['section_2']['sampling_frequency']

    # Sort the segments to eliminate dictionary randomness
    segs = list(channel_md['segments'].keys())
    segs.sort()
    for segment_name in segs:
        seg_toc = channel_md['segments'][segment_name]['TOC']

        # Join into channel TOC
        toc = np.concatenate([toc, seg_toc], axis=1)

    # Once we have all segments get lenghts (i.e. differnces between segments)
    toc[1, 1::] = (((np.diff(toc[3, :]) / 1e6) - (np.diff(toc[2, :]) / fsamp))
                   * 1e6)

    return toc


def uutc_for_sample(sample, channel_md):
    """
    Convert sample to uUTC time\n

    Parameters:\n
    -------------\n
    sample - sample to be converted\n
    channel_md - channel metadata dictionary\n

    Returns:\n
    ----------\n
    uUTC time
    """

    toc = get_toc(channel_md)
    fsamp = channel_md['section_2']['sampling_frequency']
    nsamp = channel_md['section_2']['number_of_samples']

    if sample > nsamp:
        raise RuntimeError('Sample '+str(sample)+' is too large. '
                           'Number of samples is '+str(nsamp))

    toc_idx = np.where(toc[2] <= sample)[0][-1]
    sample_diff = sample - toc[2][toc_idx]

    return int(toc[3][toc_idx] + (sample_diff/fsamp) * 1e6)


def sample_for_uutc(uutc, channel_md, return_discont_distance=False):
    """
    Convert uUTC to sample time\n

    Parameters:\n
    -------------\n
    uUTC - sample to be converted\n
    channel_md - channel metadata dictionary\n
    return_discont_distance - flag if return distance to neares discontinuity\n

    Returns:\n
    ----------\n
    sample
    """
    # UUTC check
    if not uutc_check(uutc, channel_md):
        raise RuntimeError('uUTC time out of file')

    toc = get_toc(channel_md)
    fsamp = channel_md['section_2']['sampling_frequency']

    toc_idx = np.where(toc[3] <= uutc)[0][-1]
    uutc_diff = uutc - toc[3][toc_idx]
    sample_diff = int(((uutc_diff/1e6)*fsamp) + 0.5)

    # If we are in the last block we do not need to check disconts.
    if toc_idx == len(toc[3])-1:
        return toc[2][toc_idx] + sample_diff

    block_len = toc[2][toc_idx+1] - toc[2][toc_idx]
    if sample_diff >= block_len:
        warn_str = ('uUTC time at discontinuity,'
                    'returning the first sample after'
                    'dicontinuity (not inclusive)')
        warnings.warn(warn_str, RuntimeWarning)
        return toc[2][toc_idx+1]
    else:
        return toc[2][toc_idx] + sample_diff


def uutc_check(uutc, channel_md):
    """
    Check if the uUTC is in the recording.\n

    Parameters:\n
    -----------\n
    uutc - timestamp\n
    channel_md - channel metadata dict\n

    Returns:\n
    --------\n
    True if inside the recording
    """

    channel_spec_md = channel_md['channel_specific_metadata']

    # Get the uUTC start
    chan_uutc_start = channel_spec_md['earliest_start_time']
    # Get the uUTC stop
    chan_uutc_stop = channel_spec_md['latest_end_time']

    return uutc >= chan_uutc_start and uutc <= chan_uutc_stop


def _sample_arg_merger(args):
    return pymef3_file.read_mef_ts_data(*args)


def read_ts_channels_sample(session_path, password, channel_map, sample_map,
                            process_n=None):
    """
    Reads desired channels in desired sample segment

    Parameters:
    -----------
    session_path - path to mef3 session (.mefd)\n
    password - mef3 data password\n
    channel_map - channel or list of channels to be read\n
    sample_map - list of [start,stop] samples to be loaded that correspond\n
        to channel_map. if there is only one entry the same range is applied\n
        to all channels\n
    process_n - how many processes use for reading (defualt None)\n

    Returns:
    --------
    data - numpy array [channels,samples]\n
    """

    # Check if path exists
    if not os.path.exists(session_path):
        raise FileNotFoundError(session_path+' does not exist!')

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
        raise RuntimeError('Process_n argument must be Nnoe or int')

    if process_n is not None:
        mp = Pool(process_n)

        iterator = []
        for channel, sample_ss in zip(channel_map, sample_map):

            channel_path = session_path+'/'+channel+'.timd'
            if not os.path.exists(channel_path):
                raise FileNotFoundError(channel_path+' does not exist!')
            check_password(channel_path, password)
            iterator.append([channel_path, password,
                             sample_ss[0], sample_ss[1]])

        data_list = mp.map(_sample_arg_merger, iterator)
        mp.terminate()
        if is_chan_str:
            return data_list[0]
        else:
            return data_list

    for channel, sample_ss in zip(channel_map, sample_map):
        channel_path = session_path+'/'+channel+'.timd'
        check_password(channel_path, password)
        data = pymef3_file.read_mef_ts_data(channel_path, password,
                                            sample_ss[0], sample_ss[1])
        data_list.append(data)

    if is_chan_str:
        return data_list[0]
    else:
        return data_list


def _uutc_arg_merger(args):
    return pymef3_file.read_mef_ts_data(*args)


def read_ts_channels_uutc(session_path, password, channel_map, uutc_map,
                          process_n=None):
    """
    Reads desired channels in desired time segment. Missing data at
    discontinuities are filled with NaNs.

    Parameters:
    -----------
    session_md - mef3 session metadata dictionary\n
    session_path - path to mef3 session (.mefd)\n
    password - mef3 data password\n
    channel_map - channel or list of channels to be read\n
    uutc_map - list of [start,stop] uutc times to be loaded that correspond\n
       to channel_map. if there is only one entry the same range is applied\n
       to all channels\n
    process_n - how many processes use for reading (defualt None)\n

    Returns:
    --------
    data - numpy array [channels,samples]\n
    """

    # Check if path exists
    if not os.path.exists(session_path):
        raise FileNotFoundError(session_path+' does not exist!')

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

            channel_path = session_path+'/'+channel+'.timd'
            if not os.path.exists(channel_path):
                raise FileNotFoundError(channel_path+' does not exist!')
            check_password(channel_path, password)
            iterator.append([channel_path, password,
                             sample_ss[0], sample_ss[1], True])

        data_list = mp.map(_uutc_arg_merger, iterator)
        mp.terminate()
        if is_chan_str:
            return data_list[0]
        else:
            return data_list

    for channel, uutc_ss in zip(channel_map, uutc_map):
        channel_path = session_path+'/'+channel+'.timd'
        check_password(channel_path, password)
        data = pymef3_file.read_mef_ts_data(channel_path, password,
                                            uutc_ss[0], uutc_ss[1], True)
        data_list.append(data)

    if is_chan_str:
        return data_list[0]
    else:
        return data_list


def change_channel_name(session_path, orig_channel_name, new_channel_name):
    """
    Changes name of channel

    Parameters:
    -----------
    session_path - path to mef session\n
    orig_channel_name - original channel name\n
    new_channel_name - new channel_name\n

    Returns:
    --------
    0 - on success\n
    """

    # Identify the channel folder
    channel_folders = os.listdir(session_path)

    channel_folder = [x for x in channel_folders
                      if x.startswith(orig_channel_name)][0]
    channel_path = session_path + '/' + channel_folder + '/'

    # Run through segments
    segment_folders = [x for x in os.listdir(channel_path) if 'segd' in x]
    channel_files = [x for x in os.listdir(channel_path) if 'segd' not in x]
    for segment_folder in segment_folders:
        segment_path = channel_path + segment_folder+'/'
        # Run through segment_files
        segment_files = os.listdir(segment_path)
        for segment_file in segment_files:
            file_path = segment_path + segment_file
            new_file_name = (new_channel_name
                             + segment_file[segment_file.rindex('-'):])
            new_file_path = segment_path+new_file_name
            os.rename(file_path, new_file_path)

        # Rename the segment folder
        new_segment_folder = (new_channel_name
                              + segment_folder[segment_folder.rindex('-'):])
        new_segment_path = channel_path + new_segment_folder
        os.rename(segment_path, new_segment_path)

    # Run through channel files
    for channel_file in channel_files:
        file_path = channel_path + channel_file
        new_file_name = (new_channel_name
                         + channel_file[channel_file.rindex('.'):])
        new_file_path = channel_path+new_file_name
        os.rename(file_path, new_file_path)

    # Rename the channel folder
    new_channel_folder = (new_channel_name
                          + channel_folder[channel_folder.rindex('.'):])
    new_channel_path = session_path + '/' + new_channel_folder
    os.rename(channel_path, new_channel_path)

    return 0


def read_ts_channel_basic_info(session_path, password):
    """
    Reads session time series channel names

    Parameters:
    -----------
    session_path - path to mef3 session\n
    password - mef3 data password\n

    Returns:
    --------
    channel_list\n
    """

    # Check if path exists
    if not os.path.exists(session_path):
        raise FileNotFoundError(session_path+' does not exist!')

    # Check password
    check_password(session_path, password)

    session_ts_md = pymef3_file.read_mef_session_metadata(session_path,
                                                          password, False)

    channel_list = list(session_ts_md['time_series_channels'].keys())
    channel_list.sort()

    channel_infos = []
    for channel in channel_list:

        channel_md = session_ts_md['time_series_channels'][channel]
        channel_md_spec = channel_md['channel_specific_metadata']
        channel_md_s2 = channel_md['section_2']

        fsamp = channel_md_s2['sampling_frequency']
        nsamp = channel_md_s2['number_of_samples']
        ufact = channel_md_s2['units_conversion_factor']
        unit = channel_md_s2['units_description']
        start_time = channel_md_spec['earliest_start_time']
        ch_desc = channel_md_s2['channel_description']

        channel_infos.append({'name': channel, 'fsamp': fsamp, 'nsamp': nsamp,
                              'ufact': ufact, 'unit': unit,
                              'start_time':  start_time,
                              'channel_description': ch_desc})

    return channel_infos


def check_password(mefpath, password):
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
    for path, subdirs, files in os.walk(mefpath):
        for name in files:
            mef_files.append(os.path.join(path, name))

    results = np.zeros(len(mef_files))
    for i, mef_file in enumerate(mef_files):
        results[i] = pymef3_file.check_mef_password(mef_file, password)

    if np.any(results < 0):
        raise RuntimeError('MEF password is invalid')

    return
