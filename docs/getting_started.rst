Getting started
==================
There are two ways one can interact with PyMef wrapper. Either through High level API which is aimed at users with no experience of MEF format. Or through Low level API which is mainly meant for experienced users.

High level API interaction
----------------------------------
All operations are handled through :class:`pymef.mef_session.MefSession` class. To initialize a MefSession we can use following code:

.. code-block:: python

    import numpy
    from pymef.mef_session import MefSession
	
    session_path = '/path/to/session.mefd'
    password = 'mef_password'
	
    # Either for reading
    ms = MefSession(session_path, password)

    # Or for writing
    ms = MefSession(session_path, password, False, True)
	
Writing time series data
~~~~~~~~~~~~~~~~~~~~~~~~~~~
To write time series data the time series metadata file has to be written first.

.. code-block:: python

    pwd_1 = 'table'
    pwd_2 = 'chair'
	
    start_time = 946684800000000
    end_time = start_time + 1000000

    sampling_frequency = 5000

    channel_name = 'test_channel'

    # Section 3
    section3_dict = {'recording_time_offset': 0,
                     'DST_start_time': 0,
                     'DST_end_time': 0,
                     'GMT_offset': 3600,
                     'subject_name_1': 'Olaf',
                     'subject_name_2': 'Mefson',
                     'subject_ID': '2017',
                     'recording_location': 'pub'}

    # Time series section 2 dictionry with ts section 2 fields
    section2_ts_dict = {'channel_description': 'Test_channel',
                        'session_description': 'Test_session',
                        'recording_duration': 1,
                        'reference_description': 'wine',
                        'acquisition_channel_number': 5,
                        'sampling_frequency': sampling_frequency,
                        'notch_filter_frequency_setting': 50.0,
                        'low_frequency_filter_setting': 1.0,
                        'high_frequency_filter_setting': 10.0,
                        'AC_line_frequency': 70,
                        'units_conversion_factor': 1.5,
                        'units_description': 'uV',
                        'start_sample': 0,  # Different for segments
                        'number_of_discontinuities': 1,
                        # The following entries are filled automatically during data writing
                        'maximum_native_sample_value': 0.0,
                        'minimum_native_sample_value': 0.0,
                        'number_of_blocks': 0,
                        'maximum_block_bytes': 0,
                        'maximum_block_samples': 0,
                        'maximum_difference_bytes': 0,
                        'block_interval': 0,
                        'maximum_contiguous_blocks': 0,
                        'maximum_contiguous_block_bytes': 0,
                        'maximum_contiguous_samples': 0,
                        'number_of_samples': 0}
	
    ms.write_mef_ts_segment_metadata(channel_name,
                                     segment_n,
                                     pwd_1,
                                     pwd_2,
                                     start_time,
                                     end_time,
                                     section2_ts_dict,
                                     section3_dict)
                                          
After the metadata file is in place the data itslef can be written:

.. code-block:: python

    N = 20000  # Number of samples
    samps_per_mef_block = 5000
    data = np.random.randint(-200, 200, N, dtype='int32')

    ms.write_mef_ts_segment_data(channel,
                                 0,
                                 pwd_1,
                                 pwd_2,
                                 samps_per_mef_block,
                                 raw_data)
                                     
Reading time series data
~~~~~~~~~~~~~~~~~~~~~~~~~~~
To obtain information about time series channels the function :meth:`~pymef.mef_session.MefSession.read_ts_channel_basic_info` as in the following code snippet.

.. code-block:: python

	ms.read_ts_channel_basic_info()
	
There are two options for reading time series data:
    - By sample with function :meth:`~pymef.mef_session.MefSession.read_ts_channels_sample`
    - By uUTC with function :meth:`~pymef.mef_session.MefSession.read_ts_channels_uutc`
    
Both functions have similar API. The first argument is either a channel string or a list of channel strings. The secund argument consits of list of lists with start and stop sample or uutc. Both time entries can also be None in which case the recording start / stop is used. If only one entry for start and stop is provided the same span is applied to all channels in the channel list. There is one minor difference in returned data - when reading by uUTC the function returns NaNs if no data is available in the specified time span.

.. code-block:: python

	# Returns 1D numpy array with data from recording start to recording stop
	ms.read_ts_channels_uutc(channel, [[None, None]])
	
	# Returns 1D numpy array with 2 1D numpy arrays with data from recording start to recording stop
	ms.read_ts_channels_uutc([channel, channel], [[None, None]])
	
	# Returns 1D numpy array with data from sample 5 to sample 5000
	ms.read_ts_channels_sample(channel, [[5, 5000]])