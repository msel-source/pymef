

// Python 3 extension for opeartions with Multiscale Electrophysiology Format (MEF) version 3.0
// For licences please see meflib.c

// Written by Jan, Dan, Matt, Ben

#include "pymef3.h"

#include <numpy/arrayobject.h>

#include "meflib.c"
#include "mefrec.c"


/* TODO

DONE - Fix in write_mef_times_series_data_and_indices() - passing pointer to source numpy array causes segfualt when
 freeing RED_procesging_struct

DONE Fix - video metadata writing - encoding problems, file name gets messed up, not sure why because it is pretty
 much he same code as in time series (had to copy file path coming from python - MEF_strncopy(file_path, path_out, MEF_FULL_FILE_NAME_BYTES))

DONE - Fix metadata wtiting - cannot read the written metadata and don't know why. (I had ".tmet" hidden channel
 in the directoty, not woth a day of looking for this if you ask me!!!)

DONE - Fix reading segment time series data - why allocate_file_preocessing_struc() allocating the data struct
 when MEF_FALSE is passed - results in errors if there are no data files. Which makes sense to some extent really but messes up my code :-) (fixed in pymef3.c by directly reading the metadata file instead of read_MEF_segment())

DONE - Fix structure to python mapping - if statements should be XXX != NULL, otherwise 0 is trasnlated as not
 entered which is not the same (coparing with mef "not entered" values)

DONE - Video index write function

DONE - Create decomp mef and substitute read_ts_data(copy-paste from Dan's function)

DONE - Fix occasional segfault when reading ts data (not needed since we will do decomp mef)

DONE - Create appedn ts data function

DONE - Fix uUTC reading (si8) into python. Seems incorrect now or I am missing something. (was likely connected to bad memory management)

DONE - Decomp_mef sort of thing so we can get slices of data - we should modify read_mef_ts_data for this - will
    be handled by start/stop time, if None, the data will be read from beginnig to end, ie whole channel, see Dan's C function for this

DONE - Check if write record function closes the file

DONE - Fix seizures record typ reading compilation problems

DONE - Move Python function declarations and docstrings to the header file.

DONE - Create list for file indices

DONE - Extract segment number when writing or appending data

DONE - Add ability to read universal headers (esp. for CRCs)

DONE - Check times in universal headers - why are they negative? (was connected to GMT offset to being set, corrected by a line after initialize_meflib(), ask MATT about this)

DONE - Check times in time series indices - they are going down instead of up! (solved by the line above)

DONE - Check file closing in all functions (problem was in ts data reading fixed at the end, see git)

DONE - check the times in read_channel/read_session (bug in MEF_LIB, see the notes and git)

DONE - Check maximum number_of_records / maximum_record_bytes in read_seesion - expected behavior, reading records at channel level not segments

DONE - Deal with simple quates in channel names, grrrrrrrr!!!!!!

Address recording_time_offset feature in file writing

Write the help docstrings

Fix the info at the beginning (licence)

Encryption / decryption at all levels - allow for no encryption, level 1 encryption, level 2 encryption - user
 can specify by inserting None into the password field

Address all warnings when compiling

Build in optional lossy compression

Build in opional filtration

*/

/* NOTES

All mkdirs will be in the pure python layer, perhaps apart from the segment number genration?

What is recording time offset and how to use it??? Ask Matt.

Ts and v Metadata write functions could be merged into one!!!

MEF_LIB "bug" - troubles when spaces in file paths

MEF_LIB bug - generate_file_list - had troubles with \' sign. fixed in the bash command (the above could be fix the same way), see git

MEF_LIB bug - read_MEF_session - was reading time_series_channels instead of video_channels, see git for the fix

MEF_LIB bug - read_MEF_channel/read_MEF_session - the earliest_start_time is initialized with zero hence it starts in 1970 and can never be updated, see git fot the fix

*/


/************************************************************************************/
/*********************************  Python stuff  ***********************************/
/************************************************************************************/

/* Module initialisation */
PyObject * PyInit_pymef3(void)
{
    PyObject *m = PyModule_Create(&moduledef);

    if (m == NULL)
        return NULL;

    return m;
}

/************************************************************************************/
/******************************  MEF write functions  *******************************/
/************************************************************************************/

static PyObject *write_mef_data_records(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_file_path;
    si1    *py_level_1_password;
    si1    *py_level_2_password;
    si1    *temp_str_bytes;
    
    si8     recording_start_uutc_time, recording_stop_uutc_time;
    
    PyObject    *py_record_list, *py_record_dict, *py_seiz_chan_list, *temp_o, *temp_UTF_str;
    Py_ssize_t  annot_bytes;

    // Method specific
    FILE_PROCESSING_STRUCT *gen_fps, *rec_data_fps, *rec_idx_fps;
    si8     bytes, rb_bytes, max_rec_bytes, file_offset;
    ui4     type_code, n_records, li, lj;
    ui1     *rd;
    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], record_file_name[MEF_BASE_FILE_NAME_BYTES];
    si1     path_processed;
    
    UNIVERSAL_HEADER        *uh;
    RECORD_HEADER           *rh;
    RECORD_INDEX            *ri;
    PASSWORD_DATA           *pwd;

    MEFREC_EDFA_1_0         *edfa_type;
    MEFREC_LNTP_1_0         *lntp_type;
    MEFREC_Seiz_1_0         *seiz_type;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sssllO!",
                          &py_file_path,
                          &py_level_1_password,
                          &py_level_2_password,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &PyList_Type, &py_record_list)){
        return NULL;
    }

    /// initialize MEF library
    (void) initialize_meflib();   
    MEF_globals->recording_time_offset_mode = RTO_IGNORE;  

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    uh->start_time = recording_start_uutc_time;
    uh->end_time = recording_stop_uutc_time;
    pwd = gen_fps->password_data = process_password_data(NULL, py_level_1_password, py_level_2_password, uh);
    
    // Check for directory type
    // Segment level
    path_processed = 0;
    extract_path_parts(py_file_path, path_out, name, type);
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name);
        if (!path_processed)
            MEF_strncpy(record_file_name, name, MEF_BASE_FILE_NAME_BYTES);
        path_processed = 1;
    }

    // Channel level
    if (path_processed){
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
    }
    if (!strcmp(type,TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING) | !strcmp(type,VIDEO_CHANNEL_DIRECTORY_TYPE_STRING)){
        MEF_strncpy(uh->channel_name, name, MEF_BASE_FILE_NAME_BYTES);
        if (!path_processed)
            MEF_strncpy(record_file_name, name, MEF_BASE_FILE_NAME_BYTES);
        path_processed = 1;
    }

    // Session level
    if (path_processed){
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
    } 
    if (!strcmp(type,SESSION_DIRECTORY_TYPE_STRING)){
        MEF_strncpy(uh->session_name, name, MEF_BASE_FILE_NAME_BYTES);
        if (!path_processed)
            MEF_strncpy(record_file_name, name, MEF_BASE_FILE_NAME_BYTES);
        path_processed = 1;
    }

    // Determine the number and size of records - ASK this is pretty dumb that I am calling this twice, we caould do it piecemeal instead
    n_records = (ui4) PyList_Size(py_record_list);
    bytes = UNIVERSAL_HEADER_BYTES;
    for (li = 0; li<n_records; li++){

        py_record_dict = PyList_GetItem(py_record_list, li);
        temp_o = PyDict_GetItemString(py_record_dict,"type_string");
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        type_code = *((ui4 *) temp_str_bytes);

        // Fork for different record types
        switch (type_code) {
            case MEFREC_EDFA_TYPE_CODE:
                bytes += RECORD_HEADER_BYTES; // + optional annotation
                rb_bytes = MEFREC_EDFA_1_0_BYTES;
                temp_o = PyDict_GetItemString(py_record_dict,"annotation");
                if (temp_o != NULL){
                    temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                    annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    rb_bytes += (si8) annot_bytes + 1; // strncpy copies the null termination as well
                }
                rb_bytes += 16 - (rb_bytes % 16);
                bytes += rb_bytes;
                break;
            case MEFREC_LNTP_TYPE_CODE:
                // TODO - add template but I do not know what that is!!
                bytes += RECORD_HEADER_BYTES + sizeof(MEFREC_LNTP_1_0);
                break;
            case MEFREC_Note_TYPE_CODE:
                bytes += RECORD_HEADER_BYTES;
                rb_bytes = 0;
                temp_o = PyDict_GetItemString(py_record_dict,"note");
                if (temp_o != NULL){
                    temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                    annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    rb_bytes += (si8) annot_bytes + 1;
                }
                rb_bytes += 16 - (rb_bytes % 16);
                bytes += rb_bytes;
                break;
            case MEFREC_Seiz_TYPE_CODE:
                bytes += RECORD_HEADER_BYTES + MEFREC_Seiz_1_0_BYTES;
                py_seiz_chan_list = PyDict_GetItemString(py_record_dict,"channels");
                bytes += (MEFREC_Seiz_1_0_CHANNEL_BYTES * PyList_Size(py_seiz_chan_list));
                break;
            case MEFREC_CSti_TYPE_CODE:
                bytes += RECORD_HEADER_BYTES;
                rb_bytes = MEFREC_CSti_1_0_BYTES;
                rb_bytes += 16 - (rb_bytes % 16);
                bytes += rb_bytes;
                break;
            case MEFREC_ESti_TYPE_CODE:
                bytes += RECORD_HEADER_BYTES;
                rb_bytes = MEFREC_ESti_1_0_BYTES;
                //rb_bytes += 16 - (rb_bytes % 16);// unnecessary but kept for consistency
                bytes += rb_bytes;
                break;
            case MEFREC_SyLg_TYPE_CODE:
                bytes += RECORD_HEADER_BYTES;
                rb_bytes = 0;
                temp_o = PyDict_GetItemString(py_record_dict,"text");
                if (temp_o != NULL){
                    temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                    annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    rb_bytes += (si8) annot_bytes + 1;
                }
                rb_bytes += 16 - (rb_bytes % 16);
                bytes += rb_bytes;
                break;
            default:
                bytes += 0;  // + optional annotation
                break;
                }
    }

    // Create file processing structs for record data and indices files
    rec_data_fps = allocate_file_processing_struct(bytes, RECORD_DATA_FILE_TYPE_CODE, NULL, gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(rec_data_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, record_file_name, RECORD_DATA_FILE_TYPE_STRING);
    generate_UUID(rec_data_fps->universal_header->file_UUID);

    bytes = UNIVERSAL_HEADER_BYTES + (n_records * RECORD_INDEX_BYTES);
    rec_idx_fps = allocate_file_processing_struct(bytes, RECORD_INDICES_FILE_TYPE_CODE, NULL, gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(rec_idx_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, record_file_name, RECORD_INDICES_FILE_TYPE_STRING);
    generate_UUID(rec_idx_fps->universal_header->file_UUID);
    rec_idx_fps->universal_header->maximum_entry_size = RECORD_INDEX_BYTES;
    rec_data_fps->universal_header->number_of_entries = rec_idx_fps->universal_header->number_of_entries = n_records;

    rd = rec_data_fps->records;
    ri = rec_idx_fps->record_indices;
    
    // rh->version_major = ri->version_major = 1;
    // rh->version_minor = ri->version_minor = 0;
    file_offset = ri->file_offset = UNIVERSAL_HEADER_BYTES;

    // Run through the python list, read records and write them
    max_rec_bytes = 0;
    for (li = 0; li<n_records; li++){

        // set up record header
        rh = (RECORD_HEADER *) rd;
        
        rh->encryption = ri->encryption = LEVEL_2_ENCRYPTION_DECRYPTED;  // ASK level 2 because may conatin subject identifying data
        ri->file_offset = file_offset;

        // get info from python dictionary
        py_record_dict = PyList_GetItem(py_record_list, li);
        map_python_rh(py_record_dict, rh);

        ri->time = rh->time;

        if (rh->version_major == NULL)
            rh->version_major = ri->version_major = 1;
        else
            ri->version_major = rh->version_major;
        if (rh->version_minor == NULL)
            rh->version_minor = ri->version_minor = 0;
        else
            ri->version_minor = rh->version_minor;
        
        // Done with record header, do record body
        rd += RECORD_HEADER_BYTES;

        // Fork for different record types
        rh->bytes = 0;
        type_code = *((ui4 *) rh->type_string);

        switch (type_code) {
            case MEFREC_EDFA_TYPE_CODE:

                // ASK should there by types created by this function and passed to subfunctinos?
                map_python_EDFA_type(py_record_dict, (si1 *) rd);

                // Type strings
                rh->bytes = MEFREC_EDFA_1_0_BYTES;
                MEF_strncpy(ri->type_string, MEFREC_EDFA_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_EDFA_TYPE_STRING, TYPE_BYTES);

                // Annotation
                temp_o = PyDict_GetItemString(py_record_dict,"annotation");
                if (temp_o != NULL){
                    temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                    temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char
                    rh->bytes += MEF_strcpy((si1 *) rd + MEFREC_EDFA_1_0_BYTES, temp_str_bytes);
                }

                // Pad to 16
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_LNTP_TYPE_CODE:
                map_python_LNTP_type(py_record_dict, (si1 *) rd);
                rh->bytes = MEFREC_LNTP_1_0_BYTES;
                MEF_strncpy(ri->type_string, MEFREC_LNTP_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_LNTP_TYPE_STRING, TYPE_BYTES);
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_Seiz_TYPE_CODE:
                map_python_Siez_type(py_record_dict, (si1 *) rd);
                rh->bytes = MEFREC_Seiz_1_0_BYTES;
                MEF_strncpy(ri->type_string, MEFREC_Seiz_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Seiz_TYPE_STRING, TYPE_BYTES);
                // Inidividual channels
                py_seiz_chan_list = PyDict_GetItemString(py_record_dict,"channels");
                for (lj = 0; lj < PyList_Size(py_seiz_chan_list); lj++){
                    py_record_dict = PyList_GetItem(py_seiz_chan_list, lj);
                    map_python_Siez_type_channel(py_record_dict, (si1 *) (rd+MEFREC_Seiz_1_0_BYTES+(lj*MEFREC_Seiz_1_0_CHANNEL_BYTES)));
                    rh->bytes += MEFREC_Seiz_1_0_CHANNEL_BYTES;
                }
                // No need to pad, seizure structs are 16 safe
                break;

            case MEFREC_Note_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_Note_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Note_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict,"note");
                if (temp_o != NULL){
                    temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                    temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char
                    rh->bytes += MEF_strcpy((si1 *) rd, temp_str_bytes);
                }
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_CSti_TYPE_CODE:
                map_python_CSti_type(py_record_dict, (si1 *) rd);
                MEF_strncpy(ri->type_string, MEFREC_CSti_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_CSti_TYPE_STRING, TYPE_BYTES);
                rh->bytes = MEFREC_CSti_1_0_BYTES;
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_ESti_TYPE_CODE:
                map_python_ESti_type(py_record_dict, (si1 *) rd);
                MEF_strncpy(ri->type_string, MEFREC_ESti_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_ESti_TYPE_STRING, TYPE_BYTES);
                rh->bytes = MEFREC_ESti_1_0_BYTES;

                //MEF_pad(rd, rh->bytes, 16); // unnecessary but kept for consistency
                break;

            case MEFREC_SyLg_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_SyLg_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_SyLg_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict,"text");
                if (temp_o != NULL){
                    temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                    temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char
                    rh->bytes += MEF_strcpy((si1 *) rd, temp_str_bytes);
                }
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_UnRc_TYPE_CODE:
                rh->bytes = 0;
                break;

            default:
                rh->bytes = 0;
                break;
        }

        if (rh->bytes > max_rec_bytes)
            max_rec_bytes = rh->bytes;

        rd += rh->bytes;
        file_offset += (RECORD_HEADER_BYTES + rh->bytes);
        ++ri;
    }

    rec_data_fps->universal_header->maximum_entry_size = max_rec_bytes + RECORD_HEADER_BYTES;
    rec_data_fps->directives.io_bytes = file_offset;

    write_MEF_file(rec_data_fps);
    free_file_processing_struct(rec_data_fps);

    write_MEF_file(rec_idx_fps);
    free_file_processing_struct(rec_idx_fps);

    //Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_ts_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_file_path;
    si1    *py_level_1_password;
    si1    *py_level_2_password;
    
    si8     recording_start_uutc_time, recording_stop_uutc_time;
    
    PyObject    *py_tmd2_dict, *py_md3_dict, *temp_o;
    Py_ssize_t  li;

    // Method specific
    FILE_PROCESSING_STRUCT *gen_fps, *metadata_fps;
    TIME_SERIES_INDEX               *tsi;
    UNIVERSAL_HEADER        *uh;
    PASSWORD_DATA           *pwd;
    TIME_SERIES_METADATA_SECTION_2  *tmd2;
    METADATA_SECTION_3  *md3;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si8     bytes, max_rec_bytes, file_offset;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sssllO!O!",
                          &py_file_path,
                          &py_level_1_password,
                          &py_level_2_password,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &PyDict_Type, &py_tmd2_dict,
                          &PyDict_Type, &py_md3_dict)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE; 

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    uh->start_time = recording_start_uutc_time;
    uh->end_time = recording_stop_uutc_time;
    pwd = gen_fps->password_data = process_password_data(NULL, py_level_1_password, py_level_2_password, uh);

    // Check for directory type
    extract_path_parts(py_file_path, path_out, name, type);
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name);

        // Copy the segment name for file name construction
        MEF_strncpy(segment_name, name, MEF_BASE_FILE_NAME_BYTES);

        // TODO - extact segment number
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
        if (!strcmp(type,TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING)){

            MEF_strncpy(uh->channel_name, name, MEF_BASE_FILE_NAME_BYTES);
            // Get session name
            MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
            extract_path_parts(path_in, path_out, name, type);
            MEF_strncpy(uh->session_name, name, MEF_BASE_FILE_NAME_BYTES);
        }else{
            //Fire an error that this is not time series directory - hence makes no sense to write metadata
            PyErr_SetString(PyExc_RuntimeError, "Not a time series channel, exiting...");
            PyErr_Occurred();
            return NULL;
        }
    }else{
        //Fire an error that this is not segment directory - hence makes no sense to write metadata
        PyErr_SetString(PyExc_RuntimeError, "Not a segment, exiting...");
        PyErr_Occurred();
        return NULL;
    }


    // generate level UUID into generic universal_header
    generate_UUID(gen_fps->universal_header->level_UUID);

    // set up mef3 time series metadata file
    metadata_fps = allocate_file_processing_struct(METADATA_FILE_BYTES, TIME_SERIES_METADATA_FILE_TYPE_CODE, NULL, gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(metadata_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    uh = metadata_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 1;
    uh->maximum_entry_size = METADATA_FILE_BYTES;
    initialize_metadata(metadata_fps);
    metadata_fps->metadata.section_1->section_2_encryption = LEVEL_1_ENCRYPTION_DECRYPTED;
    metadata_fps->metadata.section_1->section_3_encryption = LEVEL_2_ENCRYPTION_DECRYPTED;

    // Get time series metadata section 2 from python dict
    map_python_tmd2(py_tmd2_dict, metadata_fps->metadata.time_series_section_2);

    // Get time series metadata section 2 from python dict
    map_python_md3(py_md3_dict, metadata_fps->metadata.section_3);

    write_MEF_file(metadata_fps);
    free_file_processing_struct(metadata_fps); // TODO Commented for now because it is likly to close Python objects - solve later

    //Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_v_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_file_path;
    si1    *py_level_1_password;
    si1    *py_level_2_password;
    
    si8     recording_start_uutc_time, recording_stop_uutc_time;
    
    PyObject    *py_vmd2_dict, *py_md3_dict, *temp_o;
    Py_ssize_t  li;

    // Method specific
    FILE_PROCESSING_STRUCT *gen_fps, *metadata_fps;
    UNIVERSAL_HEADER        *uh;
    PASSWORD_DATA           *pwd;
    VIDEO_METADATA_SECTION_2  *vmd2;
    METADATA_SECTION_3  *md3;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si8     bytes, max_rec_bytes, file_offset;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sssllO!O!",
                          &py_file_path,
                          &py_level_1_password,
                          &py_level_2_password,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &PyDict_Type, &py_vmd2_dict,
                          &PyDict_Type, &py_md3_dict)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE; 

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    uh->start_time = recording_start_uutc_time;
    uh->end_time = recording_stop_uutc_time;
    pwd = gen_fps->password_data = process_password_data(NULL, py_level_1_password, py_level_2_password, uh);

    // Check for directory type
    extract_path_parts(py_file_path, path_out, name, type);
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name);

        // Copy the segment name for later use
        MEF_strncpy(segment_name, name, MEF_BASE_FILE_NAME_BYTES);

        // TODO - extact segment number
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
        if (!strcmp(type,VIDEO_CHANNEL_DIRECTORY_TYPE_STRING)){
            MEF_strncpy(uh->channel_name, name, MEF_BASE_FILE_NAME_BYTES);
            // Get session name
            MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
            extract_path_parts(path_in, path_out, name, type);
            MEF_strncpy(uh->session_name, name, MEF_BASE_FILE_NAME_BYTES);
        }else{
            //Fire an error that this is not time series directory - hence makes no sense to write metadata
            PyErr_SetString(PyExc_RuntimeError, "Not a video channel, exiting...");
            PyErr_Occurred();
            return NULL;
        }

    }else{
        //Fire an error that this is not segment directory - hence makes no sense to write metadata
        PyErr_SetString(PyExc_RuntimeError, "Not a segment, exiting...");
        PyErr_Occurred();
        return NULL;
    }

    // generate level UUID into generic universal_header
    generate_UUID(gen_fps->universal_header->level_UUID);

    // set up mef3 time series metadata file
    metadata_fps = allocate_file_processing_struct(METADATA_FILE_BYTES, VIDEO_METADATA_FILE_TYPE_CODE, NULL, gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(metadata_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", py_file_path, segment_name, VIDEO_METADATA_FILE_TYPE_STRING);
    uh = metadata_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 1;
    uh->maximum_entry_size = METADATA_FILE_BYTES;
    initialize_metadata(metadata_fps);
    metadata_fps->metadata.section_1->section_2_encryption = LEVEL_1_ENCRYPTION_DECRYPTED;
    metadata_fps->metadata.section_1->section_3_encryption = LEVEL_2_ENCRYPTION_DECRYPTED;

    // Get time series metadata section 2 from python dict
    map_python_vmd2(py_vmd2_dict, metadata_fps->metadata.video_section_2);

    // Get time series metadata section 2 from python dict
    map_python_md3(py_md3_dict, metadata_fps->metadata.section_3);

    write_MEF_file(metadata_fps);

    free_file_processing_struct(metadata_fps);

    //Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_ts_data_and_indices(PyObject *self, PyObject *args)
{
    // Specified by user
    PyObject    *raw_data;
    si1    *py_file_path;
    si1    *py_level_1_password;
    si1    *py_level_2_password;
    ui1    lossy_flag;
    
    si8    samps_per_mef_block, recording_start_uutc_time, recording_stop_uutc_time;

    // Method specific
    PASSWORD_DATA           *pwd;
    UNIVERSAL_HEADER    *uh;
    FILE_PROCESSING_STRUCT  *gen_fps, *metadata_fps, *ts_idx_fps, *ts_data_fps;
    TIME_SERIES_METADATA_SECTION_2  *tmd2;
    TIME_SERIES_INDEX   *tsi;
    RED_PROCESSING_STRUCT   *rps;
    RED_BLOCK_HEADER    *block_header;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     full_file_name[MEF_FULL_FILE_NAME_BYTES], file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si4     max_samp, min_samp, *np_array_ptr;
    si8     start_sample, ts_indices_file_bytes, n_read, samps_remaining, block_samps, file_offset;
    sf8     curr_time, time_inc;


    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"ssslOb",
                          &py_file_path, // full path including segment
                          &py_level_1_password,
                          &py_level_2_password,
                          &samps_per_mef_block,
                          &raw_data,
                          &lossy_flag)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE; 

    // NOTE: gen_fps is unecessart here if the metadata file with the universal header already exists, or is it?

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    pwd = gen_fps->password_data = process_password_data(NULL, py_level_1_password, py_level_2_password, uh);

    // Check for directory type
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    extract_path_parts(file_path, path_out, name, type);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name);

        // Copy the segment name for later use
        MEF_strncpy(segment_name, name, MEF_BASE_FILE_NAME_BYTES);

        // TODO - extact segment number
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
        if (!strcmp(type,TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING)){
            MEF_strncpy(uh->channel_name, name, MEF_BASE_FILE_NAME_BYTES);
            // Get session name
            MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
            extract_path_parts(path_in, path_out, name, type);
            MEF_strncpy(uh->session_name, name, MEF_BASE_FILE_NAME_BYTES);
        }else{
            //Fire an error that this is not time series directory - hence makes no sense to write metadata
            PyErr_SetString(PyExc_RuntimeError, "Not a time series channel, exiting...");
            PyErr_Occurred();
            return NULL;
        }

    }else{
        //Fire an error that this is not segment directory - hence makes no sense to write metadata
        PyErr_SetString(PyExc_RuntimeError, "Not a segment, exiting...");
        PyErr_Occurred();
        return NULL;
    }

    // TODO - take care of different encryptions

    // ASK There can be a fork here - we can create entierly new metadata file if there is none in the segment, python layer can also take care of this (can call write_ts_metadata first)
    // Get the metadata and update some fields - update the rest when processing RED
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    metadata_fps = read_MEF_file(NULL, full_file_name, py_level_1_password, pwd, NULL, USE_GLOBAL_BEHAVIOR);
    // seg = read_MEF_segment(NULL, file_path, TIME_SERIES_CHANNEL_TYPE, py_level_1_password, NULL, MEF_FALSE, MEF_FALSE);
    // metadata_fps = seg->metadata_fps;
    tmd2 = metadata_fps->metadata.time_series_section_2;

    tmd2->number_of_samples = (si8) PyArray_SHAPE(raw_data)[0];
    tmd2->number_of_blocks = (si8) ceil((sf8) tmd2->number_of_samples / (sf8) samps_per_mef_block);
    tmd2->maximum_block_samples = samps_per_mef_block;

    // Get the start time and end time from the metadata file
    uh->start_time = metadata_fps->universal_header->start_time;
    uh->end_time = metadata_fps->universal_header->end_time;
    
    // Set up mef3 time series indices file
    ts_indices_file_bytes = (tmd2->number_of_blocks * TIME_SERIES_INDEX_BYTES) + UNIVERSAL_HEADER_BYTES;
    ts_idx_fps = allocate_file_processing_struct(ts_indices_file_bytes, TIME_SERIES_INDICES_FILE_TYPE_CODE, NULL, metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(ts_idx_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
    uh = ts_idx_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = tmd2->number_of_blocks;
    uh->maximum_entry_size = TIME_SERIES_INDEX_BYTES;

    // Set up mef3 time series data file
    ts_data_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES + RED_MAX_COMPRESSED_BYTES(samps_per_mef_block, 1), TIME_SERIES_DATA_FILE_TYPE_CODE, NULL, metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(ts_data_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);
    uh = ts_data_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = tmd2->number_of_blocks;
    uh->maximum_entry_size = samps_per_mef_block;
    ts_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
    ts_data_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(ts_data_fps);

    // TODO optional filtration
    // use allocation below if lossy
    if (lossy_flag == 1){
        rps = RED_allocate_processing_struct(samps_per_mef_block, 0, samps_per_mef_block, RED_MAX_DIFFERENCE_BYTES(samps_per_mef_block), samps_per_mef_block, samps_per_mef_block, pwd);
        // ASK RED lossy compression user specified???
        rps->compression.mode = RED_MEAN_RESIDUAL_RATIO;
        rps->directives.detrend_data = MEF_TRUE;
        rps->directives.require_normality = MEF_TRUE;
        rps->compression.goal_mean_residual_ratio = 0.10;
        rps->compression.goal_tolerance = 0.01;
    }else{
        rps = RED_allocate_processing_struct(samps_per_mef_block, 0, 0, RED_MAX_DIFFERENCE_BYTES(samps_per_mef_block), 0, 0, pwd);
    }

    rps->block_header = (RED_BLOCK_HEADER *) (rps->compressed_data = ts_data_fps->RED_blocks);

    // create new RED blocks
    curr_time = metadata_fps->universal_header->start_time;
    time_inc = ((sf8) samps_per_mef_block / tmd2->sampling_frequency) * (sf8) 1e6;
    samps_remaining = tmd2->number_of_samples;
    block_header = rps->block_header;
    tsi = ts_idx_fps->time_series_indices;
    start_sample = 0;
    min_samp = RED_POSITIVE_INFINITY;
    max_samp = RED_NEGATIVE_INFINITY;
    block_samps = samps_per_mef_block;
    file_offset = UNIVERSAL_HEADER_BYTES; 

    // Write the data and update the metadata
    while (samps_remaining) {

        // check
        if (samps_remaining < block_samps)
            block_samps = samps_remaining;
        block_header->number_of_samples = block_samps;
        block_header->start_time = (si8) (curr_time + 0.5); // ASK Why 0.5 here?
        curr_time += time_inc;
        
        rps->original_data = rps->original_ptr = (si4 *) PyArray_DATA(raw_data) + (tmd2->number_of_samples - samps_remaining);

        // filter - comment out if don't want
        // filtps->data_length = block_samps;
        // RED_filter(filtps);

        samps_remaining -= block_samps;

        // compress
        (void) RED_encode(rps);
        ts_data_fps->universal_header->body_CRC = CRC_update((ui1 *) block_header, block_header->block_bytes, ts_data_fps->universal_header->body_CRC);
        e_fwrite((void *) block_header, sizeof(ui1), block_header->block_bytes, ts_data_fps->fp, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, EXIT_ON_FAIL);

        // time series indices
        tsi->file_offset = file_offset;
        file_offset += (tsi->block_bytes = block_header->block_bytes);
        tsi->start_time = block_header->start_time;
        tsi->start_sample = start_sample;
        start_sample += (tsi->number_of_samples = block_samps);
        RED_find_extrema(rps->original_ptr, block_samps, tsi);
        if (max_samp < tsi->maximum_sample_value)
            max_samp = tsi->maximum_sample_value;
        if (min_samp > tsi->minimum_sample_value)
            min_samp = tsi->minimum_sample_value;
        tsi->RED_block_flags = block_header->flags;
        ++tsi;

        // update metadata
        if (tmd2->maximum_block_bytes < block_header->block_bytes)
            tmd2->maximum_block_bytes = block_header->block_bytes;
        if (tmd2->maximum_difference_bytes < block_header->difference_bytes)
            tmd2->maximum_difference_bytes = block_header->difference_bytes;
    }

    // update metadata
    tmd2->maximum_contiguous_block_bytes = file_offset - UNIVERSAL_HEADER_BYTES;
    if (tmd2->units_conversion_factor >= 0.0) {
        tmd2->maximum_native_sample_value = (sf8) max_samp * tmd2->units_conversion_factor;
        tmd2->minimum_native_sample_value = (sf8) min_samp * tmd2->units_conversion_factor;
    } else {
        tmd2->maximum_native_sample_value = (sf8) min_samp * tmd2->units_conversion_factor;
        tmd2->minimum_native_sample_value = (sf8) max_samp * tmd2->units_conversion_factor;
    }
    tmd2->maximum_contiguous_blocks = tmd2->number_of_blocks;

    // Write the files
    ts_data_fps->universal_header->header_CRC = CRC_calculate(ts_data_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);
    e_fseek(ts_data_fps->fp, 0, SEEK_SET, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    e_fwrite(uh, sizeof(ui1), UNIVERSAL_HEADER_BYTES, ts_data_fps->fp, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    fclose(ts_data_fps->fp);
    // write out metadata & time series indices files
    write_MEF_file(metadata_fps);
    write_MEF_file(ts_idx_fps);

    // clean up
    free_file_processing_struct(metadata_fps);
    free_file_processing_struct(ts_data_fps);
    free_file_processing_struct(ts_idx_fps);
    free_file_processing_struct(gen_fps);
    rps->block_header = NULL;
    rps->compressed_data = NULL;
    rps->original_data = NULL;
    rps->original_ptr = NULL;
    RED_free_processing_struct(rps);

    return Py_None;
}

static PyObject *write_mef_v_indices(PyObject *self, PyObject *args)
{
    // Specified by user
    PyObject    *vi_list, *vi_dict; 
    si1    *py_file_path;
    si1    *py_level_1_password;
    si1    *py_level_2_password;
    
    si8     samps_per_mef_block, recording_start_uutc_time, recording_stop_uutc_time;

    // Method specific
    SEGMENT    *seg;
    PASSWORD_DATA           *pwd;
    UNIVERSAL_HEADER    *uh;
    FILE_PROCESSING_STRUCT  *gen_fps, *v_idx_fps;
    VIDEO_METADATA_SECTION_2  *vmd2;
    VIDEO_INDEX   *vi;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     full_file_name[MEF_FULL_FILE_NAME_BYTES], file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si4     li;
    si8     v_indices_file_bytes;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sssllO!",
                          &py_file_path, // full path including segment
                          &py_level_1_password,
                          &py_level_2_password,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &PyList_Type, &vi_list)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE; 

    // NOTE: gen_fps is unecessart here if the metadata file with the universal header already exists, or is it?

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    pwd = gen_fps->password_data = process_password_data(NULL, py_level_1_password, py_level_2_password, uh);

    // Check for directory type
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    extract_path_parts(file_path, path_out, name, type);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name);

        // Copy the segment name for later use
        MEF_strncpy(segment_name, name, MEF_BASE_FILE_NAME_BYTES);

        // TODO - extact segment number
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
        if (!strcmp(type,VIDEO_CHANNEL_DIRECTORY_TYPE_STRING)){
            MEF_strncpy(uh->channel_name, name, MEF_BASE_FILE_NAME_BYTES);
            // Get session name
            MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
            extract_path_parts(path_in, path_out, name, type);
            MEF_strncpy(uh->session_name, name, MEF_BASE_FILE_NAME_BYTES);
        }else{
            //Fire an error that this is not time series directory - hence makes no sense to write metadata
            PyErr_SetString(PyExc_RuntimeError, "Not a video channel, exiting...");
            PyErr_Occurred();
            return NULL;
        }

    }else{
        //Fire an error that this is not segment directory - hence makes no sense to write metadata
        PyErr_SetString(PyExc_RuntimeError, "Not a segment, exiting...");
        PyErr_Occurred();
        return NULL;
    }

    // Get the start time and end time from the metadata file
    uh->start_time = recording_start_uutc_time;
    uh->end_time = recording_stop_uutc_time;

    // Set up mef3 video indices file
    v_indices_file_bytes = (PyList_Size(vi_list) * VIDEO_INDEX_BYTES) + UNIVERSAL_HEADER_BYTES;
    v_idx_fps = allocate_file_processing_struct(v_indices_file_bytes, VIDEO_INDICES_FILE_TYPE_CODE, NULL, gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(v_idx_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, VIDEO_INDICES_FILE_TYPE_STRING);
    uh = v_idx_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = li<PyList_Size(vi_list);
    uh->maximum_entry_size = TIME_SERIES_INDEX_BYTES;

    // Run through the python list and create indices
    vi = v_idx_fps->video_indices;
    for (li = 0; li<PyList_Size(vi_list); li++){
        vi_dict = PyList_GetItem(vi_list, li);
        map_python_vi(vi_dict, vi);
        ++vi;
    }

    // write the file
    write_MEF_file(v_idx_fps);

    // clean up
    free_file_processing_struct(v_idx_fps);

    return Py_None;
}

/************************************************************************************/
/*************************  MEF modify/append functions  ****************************/
/************************************************************************************/

static PyObject *append_ts_data_and_indices(PyObject *self, PyObject *args)
{
    // Specified by user
    PyObject    *raw_data;
    si1    *py_file_path;
    si1    *py_level_1_password;
    si1    *py_level_2_password;
    ui1    lossy_flag;
    
    si8    samps_per_mef_block, recording_start_uutc_time, recording_stop_uutc_time;

    // Method specific
    PASSWORD_DATA           *pwd;
    UNIVERSAL_HEADER    *uh;
    FILE_PROCESSING_STRUCT  *gen_fps, *metadata_fps, *ts_idx_fps, *ts_data_fps;
    TIME_SERIES_METADATA_SECTION_2  *tmd2;
    TIME_SERIES_INDEX   *tsi;
    RED_PROCESSING_STRUCT   *rps;
    RED_BLOCK_HEADER    *block_header;
    FILE_PROCESSING_DIRECTIVES     *gen_directives;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     full_file_name[MEF_FULL_FILE_NAME_BYTES], file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si4     i, max_samp, min_samp, *np_array_ptr;
    si8     start_sample, ts_indices_file_bytes, n_read, samps_remaining, block_samps, file_offset, orig_number_of_blocks;
    sf8     curr_time, time_inc;

    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE; 

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"ssslllOb",
                          &py_file_path, // full path including segment
                          &py_level_1_password,
                          &py_level_2_password,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &samps_per_mef_block,
                          &raw_data,
                          &lossy_flag)){
        return NULL;
    }

    // We don't really need this
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    pwd = gen_fps->password_data = process_password_data(NULL, py_level_1_password, py_level_2_password, uh);

    // Check for directory type
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    extract_path_parts(file_path, path_out, name, type);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){

        // Copy the segment name for later use
        MEF_strncpy(segment_name, name, MEF_BASE_FILE_NAME_BYTES);

        // TODO - extact segment number
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
        if (!strcmp(type,TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING)){
            // Get session name
            MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
            extract_path_parts(path_in, path_out, name, type);
        }else{
            //Fire an error that this is not time series directory - hence makes no sense to write metadata
            PyErr_SetString(PyExc_RuntimeError, "Not a time series channel, exiting...");
            PyErr_Occurred();
            return NULL;
        }

    }else{
        //Fire an error that this is not segment directory - hence makes no sense to write metadata
        PyErr_SetString(PyExc_RuntimeError, "Not a segment, exiting...");
        PyErr_Occurred();
        return NULL;
    }

    // Create directives so that we can modify the files, not just read them?????
    gen_directives = initialize_file_processing_directives(NULL);
    gen_directives->open_mode = FPS_R_PLUS_OPEN_MODE;
    gen_directives->close_file = MEF_FALSE;

    // Read in the metadata file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    metadata_fps = read_MEF_file(NULL, full_file_name, py_level_1_password, pwd, NULL, USE_GLOBAL_BEHAVIOR);
    tmd2 = metadata_fps->metadata.time_series_section_2;
    // We are appending so get only the end time
    metadata_fps->universal_header->end_time = recording_stop_uutc_time;

    orig_number_of_blocks = tmd2->number_of_blocks;
    tmd2->number_of_blocks +=  (si8) ceil((sf8) PyArray_SHAPE(raw_data)[0] / (sf8) samps_per_mef_block);
    if (samps_per_mef_block > tmd2->maximum_block_samples)
        tmd2->maximum_block_samples = samps_per_mef_block;

    // Read in the indices file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
    ts_idx_fps = read_MEF_file(NULL, full_file_name, py_level_1_password, pwd, gen_directives, USE_GLOBAL_BEHAVIOR);

    // Read in the time series data file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);
    ts_data_fps = read_MEF_file(NULL, full_file_name, py_level_1_password, pwd, gen_directives, USE_GLOBAL_BEHAVIOR);
    
    // TODO optional filtration
    // use allocation below if lossy
    if (lossy_flag == 1){
        rps = RED_allocate_processing_struct(samps_per_mef_block, 0, samps_per_mef_block, RED_MAX_DIFFERENCE_BYTES(samps_per_mef_block), samps_per_mef_block, samps_per_mef_block, pwd);
        // ASK RED lossy compression user specified???
        rps->compression.mode = RED_MEAN_RESIDUAL_RATIO;
        rps->directives.detrend_data = MEF_TRUE;
        rps->directives.require_normality = MEF_TRUE;
        rps->compression.goal_mean_residual_ratio = 0.10;
        rps->compression.goal_tolerance = 0.01;
    }else{
        rps = RED_allocate_processing_struct(samps_per_mef_block, 0, 0, RED_MAX_DIFFERENCE_BYTES(samps_per_mef_block), 0, 0, pwd);
    }

    //rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;

    rps->block_header = (RED_BLOCK_HEADER *) (rps->compressed_data = ts_data_fps->RED_blocks);

    // TODO - take care of discontinuity flags here!!!
    // create new RED blocks
    curr_time = recording_start_uutc_time;
    time_inc = ((sf8) samps_per_mef_block / tmd2->sampling_frequency) * (sf8) 1e6;
    samps_remaining = (si8) PyArray_SHAPE(raw_data)[0];
    block_header = rps->block_header;
    start_sample = tmd2->number_of_samples + 1;
    tmd2->number_of_samples = tmd2->number_of_samples + (si8) PyArray_SHAPE(raw_data)[0];
    min_samp = RED_POSITIVE_INFINITY;
    max_samp = RED_NEGATIVE_INFINITY;
    block_samps = samps_per_mef_block; 
    np_array_ptr = (si4 *) PyArray_DATA(raw_data);

    //Move file_offset to the end of RED blocks
    file_offset = UNIVERSAL_HEADER_BYTES + ts_data_fps->raw_data_bytes;

    // fseek to the end of data and indices file
    e_fseek(ts_data_fps->fp, 0, SEEK_END, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    e_fseek(ts_idx_fps->fp, 0, SEEK_END, ts_idx_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);

    // allocate time_series_index
    tsi = e_calloc(1, TIME_SERIES_INDEX_BYTES, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);

    // Write the data and update the metadata
    while (samps_remaining) {

        // check
        if (samps_remaining < block_samps)
            block_samps = samps_remaining;
        block_header->number_of_samples = block_samps;
        block_header->start_time = (si8) (curr_time + 0.5); // ASK Why 0.5 here?
        curr_time += time_inc;

        rps->original_data = rps->original_ptr = (si4 *) PyArray_DATA(raw_data) + ((si8) PyArray_SHAPE(raw_data)[0] - samps_remaining);

        // filter - comment out if don't want
        // filtps->data_length = block_samps;
        // RED_filter(filtps);

        samps_remaining -= block_samps;

        // compress and write blocks
        (void) RED_encode(rps);
        ts_data_fps->universal_header->body_CRC = CRC_update((ui1 *) block_header, block_header->block_bytes, ts_data_fps->universal_header->body_CRC);
        e_fwrite((void *) block_header, sizeof(ui1), block_header->block_bytes, ts_data_fps->fp, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, EXIT_ON_FAIL);

        // time series indices
        tsi->file_offset = file_offset;
        file_offset += (tsi->block_bytes = block_header->block_bytes);
        tsi->start_time = block_header->start_time;
        tsi->start_sample = start_sample;
        start_sample += (tsi->number_of_samples = block_samps);
        RED_find_extrema(rps->original_ptr, block_samps, tsi);
        if (max_samp < tsi->maximum_sample_value)
            max_samp = tsi->maximum_sample_value;
        if (min_samp > tsi->minimum_sample_value)
            min_samp = tsi->minimum_sample_value;
        tsi->RED_block_flags = block_header->flags;

        // write the time index entry
        ts_idx_fps->universal_header->body_CRC = CRC_update((ui1 *) tsi, TIME_SERIES_INDEX_BYTES, ts_idx_fps->universal_header->body_CRC);
        e_fwrite((void *) tsi, TIME_SERIES_INDEX_BYTES, 1, ts_idx_fps->fp, ts_idx_fps->full_file_name, __FUNCTION__, __LINE__, EXIT_ON_FAIL);

        // update metadata
        if (tmd2->maximum_block_bytes < block_header->block_bytes)
            tmd2->maximum_block_bytes = block_header->block_bytes;
        if (tmd2->maximum_difference_bytes < block_header->difference_bytes)
            tmd2->maximum_difference_bytes = block_header->difference_bytes;
    }

    // update metadata
    tmd2->maximum_contiguous_block_bytes = file_offset - UNIVERSAL_HEADER_BYTES;
    if (tmd2->units_conversion_factor >= 0.0) {
        tmd2->maximum_native_sample_value = (sf8) max_samp * tmd2->units_conversion_factor;
        tmd2->minimum_native_sample_value = (sf8) min_samp * tmd2->units_conversion_factor;
    } else {
        tmd2->maximum_native_sample_value = (sf8) min_samp * tmd2->units_conversion_factor;
        tmd2->minimum_native_sample_value = (sf8) max_samp * tmd2->units_conversion_factor;
    }

    // Update the header of data file
    uh = ts_data_fps->universal_header;
    uh->number_of_entries = tmd2->number_of_blocks;
    uh->maximum_entry_size = samps_per_mef_block;
    uh->start_time = metadata_fps->universal_header->start_time;
    uh->end_time = metadata_fps->universal_header->end_time;
    ts_data_fps->universal_header->header_CRC = CRC_calculate(ts_data_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);
    e_fseek(ts_data_fps->fp, 0, SEEK_SET, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    e_fwrite(uh, sizeof(ui1), UNIVERSAL_HEADER_BYTES, ts_data_fps->fp, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    
    // Update the header of indices file
    uh = ts_idx_fps->universal_header;
    uh->number_of_entries = tmd2->number_of_blocks;
    uh->maximum_entry_size = TIME_SERIES_INDEX_BYTES;
    uh->start_time = metadata_fps->universal_header->start_time;
    uh->end_time = metadata_fps->universal_header->end_time;
    ts_idx_fps->universal_header->header_CRC = CRC_calculate(ts_idx_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);
    e_fseek(ts_idx_fps->fp, 0, SEEK_SET, ts_idx_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    e_fwrite(uh, sizeof(ui1), UNIVERSAL_HEADER_BYTES, ts_idx_fps->fp, ts_idx_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    
    // Update the metadta file
    write_MEF_file(metadata_fps);

    // Close the file pointers
    fclose(ts_data_fps->fp);
    fclose(ts_idx_fps->fp);

    // clean up
    free_file_processing_struct(metadata_fps);
    free_file_processing_struct(ts_data_fps);
    free_file_processing_struct(ts_idx_fps);
    free_file_processing_struct(gen_fps);
    rps->block_header = NULL;
    rps->compressed_data = NULL;
    rps->original_data = NULL;
    RED_free_processing_struct(rps);

    return Py_None;
}

// ASK No need for modify functions - can be taken care of at python level - just load and rewrite,
// memory load would be minute in thes cases.

// static PyObject *modify_mef_data_record(PyObject *self, PyObject *args)
// {

// }

// static PyObject *modify_mef_ts_metadata(PyObject *self, PyObject *args)
// {

// }

// static PyObject *modify_mef_v_metadata(PyObject *self, PyObject *args)
// {
    
// }

/************************************************************************************/
/******************************  MEF read functions  ********************************/
/************************************************************************************/

static PyObject *read_mef_session_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_session_path;
    si1    *py_level_1_password;
    
    // Dictionaries
    PyObject *ses_metadata_dict;
    
    // Fuction specific
    SESSION *session;
    si1     session_path[MEF_FULL_FILE_NAME_BYTES];
 
    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"ss",
                          &py_session_path,
                          &py_level_1_password)){
        return NULL;
    }
    
    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE;

    MEF_strncpy(session_path, py_session_path, MEF_FULL_FILE_NAME_BYTES);
    session = read_MEF_session(NULL, session_path, py_level_1_password, NULL, MEF_FALSE, MEF_TRUE);    

    // Session info
    ses_metadata_dict = map_mef3_session(session);

    // clean up
    free_session(session, MEF_TRUE); 
    
    return ses_metadata_dict;   
}

static PyObject *read_mef_channel_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_channel_dir;
    si1    *py_level_1_password;
    
    // Dictionaries
    PyObject *ch_metadata_dict;
    
    // Fuction specific
    CHANNEL *channel;
 
    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"ss",
                          &py_channel_dir,
                          &py_level_1_password)){
        return NULL;
    }
    
    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE;
    
    channel = read_MEF_channel(NULL, py_channel_dir, UNKNOWN_CHANNEL_TYPE, py_level_1_password, NULL, MEF_FALSE, MEF_TRUE);    
    
    // map the channel info
    ch_metadata_dict = map_mef3_channel(channel);

    // clean up
    free_channel(channel, MEF_TRUE);
    
    return ch_metadata_dict;   
} 

static PyObject *read_mef_segment_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_segment_dir;
    si1    *py_level_1_password;
    
    // Dictionaries
    PyObject *seg_metadata_dict;
    
    // Fuction specific
    SEGMENT *segment;
  
    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"ss",
                          &py_segment_dir,
                          &py_level_1_password)){
        return NULL;
    }
    
    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE;

    segment = read_MEF_segment(NULL, py_segment_dir, UNKNOWN_CHANNEL_TYPE, py_level_1_password, NULL, MEF_FALSE, MEF_TRUE);    

    // map the segment info
    seg_metadata_dict = map_mef3_segment(segment);

    // clean up
    free_segment(segment, MEF_TRUE);
    
    return seg_metadata_dict; 
}

static PyObject *read_mef_ts_data(PyObject *self, PyObject *args) // This will be subsituted for read_mef_ts_data - do not want to break code now
{
    // Specified by user
    si1    *py_channel_path;
    si1    *py_level_1_password;
    si8    start_samp;
    si8    end_samp;
 
    // Python variables
    PyObject    *py_array_out;
    PyObject    *ostart_samp, *oend_samp;
    
    
    // Method specific variables
    si1             *compressed_data_buffer, *cdp;
    si1             *in_file_name, *out_file_name;
    si1             channel_path[MEF_FULL_FILE_NAME_BYTES];
    ui4             n_segments, max_samps;
    si4             start_segment, end_segment, sample_counter;
    si4             n_blocks_in_segment, offset_to_start_samp;
    si4             order, pass_order, stop_order;
    si4             *raw_data_buffer, *idp;
    si4             *temp_data_buf;
    ui8             n_samples, total_data_bytes, bytes_to_read;
    ui8             start_idx, end_idx, num_blocks, data_len, n_read;
    si8             i, j;
    si8             segment_start_sample, segment_end_sample;
    si8             total_samps, samp_counter_base;
    si8             block_start_time;
    sf8             *filt_data;

    // RED decompression
    ui4             kept_samples, skipped_samples, tot_samples;
    ui4             *dcdp;

    extern MEF_GLOBALS      *MEF_globals;
    FILE                *out_fp, *in_fp, *fp;
    RED_PROCESSING_STRUCT       *rps;
    PASSWORD_DATA           *pwd;
    CHANNEL             *chan;
    SEGMENT             *seg;
    TIME_SERIES_METADATA_SECTION_2  *md2;
    TIME_SERIES_INDEX       *tsi;

    
    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"ssOO",
                          &py_channel_path,
                          &py_level_1_password,
                          &ostart_samp,
                          &oend_samp)){
        return NULL;
    }
        
    // check inputs


    // initialize MEF library
    (void) initialize_meflib();
    MEF_globals->recording_time_offset_mode = RTO_IGNORE;

    MEF_strncpy(channel_path, py_channel_path, MEF_FULL_FILE_NAME_BYTES); // might be unnecesasry
    chan = read_MEF_channel(NULL, channel_path, UNKNOWN_CHANNEL_TYPE, py_level_1_password, NULL, MEF_FALSE, MEF_FALSE);  
    // else pass NULL instead of NULL
    
    if (chan->channel_type != TIME_SERIES_CHANNEL_TYPE) {
        PyErr_SetString(PyExc_RuntimeError, "Not a time series channel, exiting...");
        PyErr_Occurred();
        return NULL;
    }

    // Take care of the input
    if (ostart_samp == Py_None)
        start_samp = 0;
    else
        start_samp = PyLong_AsLong(ostart_samp);

    if (oend_samp == Py_None)
        end_samp = chan->metadata.time_series_section_2->number_of_samples;
    else
        end_samp = PyLong_AsLong(oend_samp);

    // Fire a warning if the stop samp is larger than number of samples
    if (end_samp > chan->metadata.time_series_section_2->number_of_samples){
        PyErr_WarnEx(PyExc_RuntimeWarning, "Stop sample larger than number of samples. Will append zeros", 1);
    }

    // Number of samples to be read
    n_samples = end_samp - start_samp;

    if (n_samples == 0){
        PyErr_SetString(PyExc_RuntimeError, "Number of samples to be read is 0. The end sample is not inclusive!");
        PyErr_Occurred();
        return NULL;
    }

    // Allocate numpy array
    npy_intp dims[1] = {n_samples};
    import_array();
    py_array_out = PyArray_SimpleNew(1, dims, NPY_INT);

    // Iterate through segments, looking for data that matches our criteria 
    n_segments = chan->number_of_segments;
    start_segment = end_segment = -1; // ASK - this shouldn't be unnecessary - taken care of at the beginnig

    // With samples, it's a little easier, since there are no sample number gaps between segments.
    // so just iterate through segments and see which segments the start and end are in.
    for (i = 0; i < n_segments; ++i) {
        
        segment_start_sample = chan->segments[i].metadata_fps->metadata.time_series_section_2->start_sample;
        segment_end_sample   = chan->segments[i].metadata_fps->metadata.time_series_section_2->start_sample +
        chan->segments[i].metadata_fps->metadata.time_series_section_2->number_of_samples;
        
        if ((start_samp >= segment_start_sample) &&
            (start_samp <= segment_end_sample))
            start_segment = i;
        if ((end_samp >= segment_start_sample) &&
            (end_samp <= segment_end_sample))
            end_segment = i;
    }
    
    // find start block in start segment
    samp_counter_base = chan->segments[start_segment].metadata_fps->metadata.time_series_section_2->start_sample;
    for (j = 1; j < chan->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks; j++) {
        
        block_start_time = chan->segments[start_segment].time_series_indices_fps->time_series_indices[j].start_time;
        remove_recording_time_offset( &block_start_time); // ASK why this???
        
        if (chan->segments[start_segment].time_series_indices_fps->time_series_indices[j].start_sample + samp_counter_base > start_samp) {
            start_idx = j - 1;
            break;
        }
        // starting point is in last block in segment
        start_idx = j;
    }
    
    // find stop block in stop segment
    samp_counter_base = chan->segments[end_segment].metadata_fps->metadata.time_series_section_2->start_sample;
    for (j = 1; j < chan->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks; j++) {
        
        block_start_time = chan->segments[end_segment].time_series_indices_fps->time_series_indices[j].start_time;
        remove_recording_time_offset( &block_start_time);
        
        if (chan->segments[start_segment].time_series_indices_fps->time_series_indices[j].start_sample + samp_counter_base > end_samp) {
            end_idx = j - 1;
            break;
        }
        // ending point is in last block in segment
        end_idx = j;
    }

    // find total_samps and total_data_bytes, so we can allocate buffers
    total_samps = 0;
    total_data_bytes = 0;

    // normal case - everything is in one segment
    if (start_segment == end_segment) {
        if (end_idx < (chan->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            total_samps += chan->segments[start_segment].time_series_indices_fps->time_series_indices[end_idx+1].start_sample -
            chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
            //fprintf(stderr, "total_samps = %d\n", total_samps);
            total_data_bytes += chan->segments[start_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        }
        else {
            // case where end_idx is last block in segment
            total_samps += chan->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
            chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
            total_data_bytes += chan->segments[start_segment].time_series_data_fps->file_length -
            chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        }
        num_blocks = end_idx - start_idx + 1;
        // TBD do we car about this?
        // TBD this max/min calculation doesn't take into account filter/re-sampling.  Also, need to make worth with multiple segments.
        //for (i = start_idx;i <= end_idx;i++)
        //{
        //    // fill in max and min values
        //    if (channel->segments[start_segment].time_series_indices_fps->time_series_indices[i].maximum_sample_value > chan->max_sample_value)
        //        chan->max_sample_value = channel->segments[start_segment].time_series_indices_fps->time_series_indices[i].maximum_sample_value;
        //    if (channel->segments[start_segment].time_series_indices_fps->time_series_indices[i].minimum_sample_value < chan->min_sample_value)
        //        chan->min_sample_value = channel->segments[start_segment].time_series_indices_fps->time_series_indices[i].minimum_sample_value;
        //}
    }
    // spans across segments
    else {
        // start with first segment
        n_blocks_in_segment = chan->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        total_samps += chan->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
        chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
        total_data_bytes +=  chan->segments[start_segment].time_series_data_fps->file_length -
        chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        num_blocks = n_blocks_in_segment - start_idx;
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            n_blocks_in_segment = chan->segments[i].metadata_fps->metadata.time_series_section_2->number_of_blocks;
            total_samps += chan->segments[i].metadata_fps->metadata.time_series_section_2->number_of_samples;
            total_data_bytes += chan->segments[i].time_series_data_fps->file_length -
            chan->segments[i].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += n_blocks_in_segment;
        }
        
        // then last segment
        n_blocks_in_segment = chan->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        if (end_idx < (chan->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            total_samps += chan->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].start_sample -
            chan->segments[end_segment].time_series_indices_fps->time_series_indices[0].start_sample;
            total_data_bytes += chan->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            chan->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += end_idx + 1;
        }
        else {
            // case where end_idx is last block in segment
            total_samps += chan->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
            chan->segments[end_segment].time_series_indices_fps->time_series_indices[0].start_sample;
            total_data_bytes += chan->segments[end_segment].time_series_data_fps->file_length -
            chan->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += end_idx + 1;
        }
    }

    // allocate buffers
    data_len = total_samps;
    order = (pass_order > stop_order) ? pass_order : stop_order;
    filt_data = (sf8 *) calloc((size_t) data_len + (6 * order), sizeof(sf8));
    compressed_data_buffer = (si1 *) malloc((size_t) total_data_bytes);
    raw_data_buffer = (si4 *) malloc((size_t) (total_samps * sizeof(si4)));
    cdp = compressed_data_buffer;
    idp = PyArray_GETPTR1(py_array_out, 0);
    
    // read in RED data
    // normal case - everything is in one segment
    if (start_segment == end_segment) {
        fp = chan->segments[start_segment].time_series_data_fps->fp;
        fseek(fp, chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        n_read = fread(cdp, sizeof(si1), (size_t) total_data_bytes, fp);
    }
    // spans across segments
    else {
        // start with first segment
        fp = chan->segments[start_segment].time_series_data_fps->fp;
        fseek(fp, chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        bytes_to_read = chan->segments[start_segment].time_series_data_fps->file_length -
        chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
        cdp += bytes_to_read;
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            fp = chan->segments[i].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = chan->segments[i].time_series_data_fps->file_length;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            cdp += bytes_to_read;
        }
        
        // then last segment
        fp = chan->segments[end_segment].time_series_data_fps->fp;
        fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
        bytes_to_read = chan->segments[end_segment].time_series_data_fps->file_length -
        chan->segments[end_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
        cdp += bytes_to_read;
    }

    // set up RED processing struct
    cdp = compressed_data_buffer;
    max_samps = chan->metadata.time_series_section_2->maximum_block_samples;
    
    // create RED processing struct
    rps = (RED_PROCESSING_STRUCT *) calloc((size_t) 1, sizeof(RED_PROCESSING_STRUCT));
    rps->compression.mode = RED_DECOMPRESSION;
    //rps->directives.return_block_extrema = MEF_TRUE;
    rps->decompressed_ptr = rps->decompressed_data = idp;
    rps->difference_buffer = (si1 *) e_calloc((size_t) RED_MAX_DIFFERENCE_BYTES(max_samps), sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    offset_to_start_samp = start_samp - chan->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
    
    //thread_info->number_of_discontinuities = 0;
    sample_counter = 0;
    
    // decode first block to temp array
    cdp = compressed_data_buffer;  

    // TBD use real max block length
    temp_data_buf = (int *) malloc(33000 * 4);
    rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
    rps->compressed_data = cdp;
    rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
    RED_decode(rps);
    cdp += rps->block_header->block_bytes;

    // copy requested samples from first block to output buffer
    skipped_samples = offset_to_start_samp;
    if (skipped_samples > rps->block_header->number_of_samples) {
        //this is bad- likely means idx data is corrupt 
        printf("[decomp_mef] block indexing error: decoded %d samples, attempting to skip %lu samples\n", rps->block_header->number_of_samples, skipped_samples);
        free (temp_data_buf);
        free (compressed_data_buffer);
        free (rps->difference_buffer);
        free (rps);
        PyErr_SetString(PyExc_RuntimeError, "Oh no, data are likely corrupt...");
        PyErr_Occurred();
        return NULL;
    }
    kept_samples = rps->block_header->number_of_samples - skipped_samples;
    tot_samples = (unsigned int) (end_samp - start_samp + 1);
    if (kept_samples >= tot_samples) { // start and end indices in same block => already done
        memcpy((void *) idp, (void *) (temp_data_buf + skipped_samples), tot_samples * sizeof(int));
        free (temp_data_buf);
        free (compressed_data_buffer);
        free (rps->difference_buffer);
        free (rps);
        return py_array_out;
    }
    memcpy((void *) idp, (void *) (temp_data_buf + skipped_samples), kept_samples * sizeof(int));
    dcdp = idp + kept_samples;
    
    rps->decompressed_ptr = rps->decompressed_data = dcdp;
    
    // decode bytes to samples
    for (i=1;i<num_blocks-1;i++) {
        rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        RED_decode(rps);
        cdp += rps->block_header->block_bytes;
        dcdp += rps->block_header->number_of_samples;
        rps->decompressed_ptr = rps->decompressed_data = dcdp;
        sample_counter += rps->block_header->number_of_samples;
    }
    
        // decode last block to temp array
            rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
        RED_decode(rps);
    
    // copy requested samples from last block to output buffer
    kept_samples = (unsigned int) (end_samp - chan->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx].start_sample + 1);
    memcpy((void *) dcdp, (void *) temp_data_buf, kept_samples * sizeof(int));
    
    // we're done with the compressed data, get rid of it
    free (temp_data_buf);
    free (compressed_data_buffer);
    free (rps->difference_buffer);
    free (rps);
    
    for (i=0;i<chan->number_of_segments;i++){
        chan->segments[i].metadata_fps->directives.free_password_data = MEF_TRUE;
        chan->segments[i].metadata_fps->directives.close_file = MEF_TRUE;
        if (chan->segments[i].record_data_fps != NULL)
            chan->segments[i].record_data_fps->directives.close_file = MEF_TRUE;
        if (chan->segments[i].record_indices_fps != NULL)
            chan->segments[i].record_indices_fps->directives.close_file = MEF_TRUE;
        chan->segments[i].time_series_data_fps->directives.close_file = MEF_TRUE;
        chan->segments[i].time_series_indices_fps->directives.close_file = MEF_TRUE;
    }
    free_channel(chan, MEF_TRUE);

    return py_array_out;
}
 
/************************************************************************************/
/*******************************  Mapper functions  *********************************/
/************************************************************************************/


/* Set of functions that map mef3 structures into python dictionaries and are
hidden from the python api. Structures are in the same order as in meflib.h.*/

/* Notes for self - size_type X Python conversions
mef3    Py_BuildValue   Python2C
si1     "b"             (char) PyLong_AsLong
ui1     "B"             (unsigned char) PyLong_AsLong
si2     "h"             (short) PyLong_AsLong
ui2     "H"             (unsigned short)PyLong_AsLong
si4     "i"             (int) PyLong_AsLong
ui4     "I"             (unsigned int) PyLong_AsLong
si8     "l"             PyLong_AsLong
ui8     "k"             PyLong_AsUnsignedLong
sf4     "f"             (float) PyFloat_AsDouble
sf8     "d"             PyFloat_AsDouble
sf16     ?

---- for strings a bit more complicated

temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
MEF_strcpy(rh->type_string, temp_str_bytes); // assign to char that we want

*/

/*****************************  Python to Mef struct  *******************************/

// UNIVERSAL_HEADER *map_python_uh(PyObject *)
// {
// }

// METADATA_SECTION_1 *map_python_md1(PyObject *md1_dict)
// {
// }

void    map_python_tmd2(PyObject *tmd2_dict, TIME_SERIES_METADATA_SECTION_2 *tmd2)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;
    si1     *temp_str_bytes;

    // Declare tmd2 struct
    //TIME_SERIES_METADATA_SECTION_2  *tmd2;

    // Assign from dict to struct

    // Type independent fields
    temp_o = PyDict_GetItemString(tmd2_dict,"channel_description");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(tmd2->channel_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(tmd2_dict,"session_description");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(tmd2->session_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(tmd2_dict,"recording_duration");
    if (temp_o != NULL)
        tmd2->recording_duration = PyLong_AsLong(temp_o);

    // Time series specific fields
    temp_o = PyDict_GetItemString(tmd2_dict,"reference_description");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(tmd2->reference_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(tmd2_dict,"acquisition_channel_number");
    if (temp_o != NULL)
        tmd2->acquisition_channel_number = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"sampling_frequency");
    if (temp_o != NULL)
        tmd2->sampling_frequency = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"low_frequency_filter_setting");
    if (temp_o != NULL)
        tmd2->low_frequency_filter_setting = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"high_frequency_filter_setting");
    if (temp_o != NULL)
        tmd2->high_frequency_filter_setting = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"notch_filter_frequency_setting");
    if (temp_o != NULL)
        tmd2->notch_filter_frequency_setting = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"AC_line_frequency");
    if (temp_o != NULL)
        tmd2->AC_line_frequency = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"units_conversion_factor");
    if (temp_o != NULL)
        tmd2->units_conversion_factor = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"units_description");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(tmd2->units_description, temp_str_bytes);
    }

    // The fields below will be written when writing data. I am keeping these in case
    // something bad happens and one would want to fix stuff

    temp_o = PyDict_GetItemString(tmd2_dict,"maximum_native_sample_value");
    if (temp_o != NULL)
        tmd2->maximum_native_sample_value = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"minimum_native_sample_value");
    if (temp_o != NULL)
        tmd2->minimum_native_sample_value = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"start_sample");
    if (temp_o != NULL)
        tmd2->start_sample = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"number_of_samples");
    if (temp_o != NULL)
        tmd2->number_of_samples = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"number_of_blocks");
    if (temp_o != NULL)
        tmd2->number_of_blocks = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"maximum_block_bytes");
    if (temp_o != NULL)
        tmd2->maximum_block_bytes = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"maximum_block_samples");
    if (temp_o != NULL)
        tmd2->maximum_block_samples = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"maximum_difference_bytes");
    if (temp_o != NULL)
        tmd2->maximum_difference_bytes = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"block_interval");
    if (temp_o != NULL)
        tmd2->block_interval = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"number_of_discontinuities");
    if (temp_o != NULL)
        tmd2->number_of_discontinuities = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"maximum_contiguous_blocks");
    if (temp_o != NULL)
        tmd2->maximum_contiguous_blocks = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"maximum_contiguous_block_bytes");
    if (temp_o != NULL)
        tmd2->maximum_contiguous_block_bytes = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(tmd2_dict,"maximum_contiguous_samples");
    if (temp_o != NULL)
        tmd2->maximum_contiguous_samples = PyLong_AsLong(temp_o);

    return;
}

void    map_python_vmd2(PyObject *vmd2_dict, VIDEO_METADATA_SECTION_2 *vmd2)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;
    si1     *temp_str_bytes;

    // Declare tmd2 struct
    //VIDEO_METADATA_SECTION_2  *vmd2;

    // Assign from dict to struct

    // Type independent fields
    temp_o = PyDict_GetItemString(vmd2_dict,"channel_description");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(vmd2->channel_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(vmd2_dict,"session_description");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(vmd2->session_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(vmd2_dict,"recording_duration");
    if (temp_o != NULL)
        vmd2->recording_duration = PyLong_AsLong(temp_o);

    // Video specific fields
    temp_o = PyDict_GetItemString(vmd2_dict,"horizontal_resolution");
    if (temp_o != NULL)
        vmd2->horizontal_resolution = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vmd2_dict,"vertical_resolution");
    if (temp_o != NULL)
        vmd2->vertical_resolution = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vmd2_dict,"frame_rate");
    if (temp_o != NULL)
        vmd2->frame_rate = PyFloat_AsDouble(temp_o);

    temp_o = PyDict_GetItemString(vmd2_dict,"number_of_clips");
    if (temp_o != NULL)
        vmd2->number_of_clips = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vmd2_dict,"maximum_clip_bytes");
    if (temp_o != NULL)
        vmd2->maximum_clip_bytes = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vmd2_dict,"video_format");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(vmd2->video_format, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(vmd2_dict,"video_file_CRC");
    if (temp_o != NULL)
        vmd2->video_file_CRC = PyLong_AsLong(temp_o);

    return;
}

void    map_python_vi(PyObject *vi_dict, VIDEO_INDEX *vi)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;
    si1     *temp_str_bytes;

    temp_o = PyDict_GetItemString(vi_dict,"start_time");
    if (temp_o != NULL)
        vi->start_time = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vi_dict,"end_time");
    if (temp_o != NULL)
        vi->end_time = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vi_dict,"start_frame");
    if (temp_o != NULL)
        vi->start_frame = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vi_dict,"end_frame");
    if (temp_o != NULL)
        vi->end_frame = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vi_dict,"file_offset");
    if (temp_o != NULL)
        vi->file_offset = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(vi_dict,"clip_bytes");
    if (temp_o != NULL)
        vi->clip_bytes = PyLong_AsLong(temp_o);

    return;
}

void    map_python_md3(PyObject *md3_dict, METADATA_SECTION_3 *md3)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;
    si1     *temp_str_bytes;

    // Declare tmd2 struct
    //METADATA_SECTION_3  *md3;

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(md3_dict,"recording_time_offset");
    if (temp_o != NULL)
        md3->recording_time_offset = PyLong_AsLong(temp_o);
    else
        md3->recording_time_offset = METADATA_RECORDING_TIME_OFFSET_NO_ENTRY;

    temp_o = PyDict_GetItemString(md3_dict,"DST_start_time");
    if (temp_o != NULL)
        md3->DST_start_time = PyLong_AsLong(temp_o);
    else
        md3->DST_start_time = METADATA_DST_START_TIME_NO_ENTRY;

    temp_o = PyDict_GetItemString(md3_dict,"DST_end_time");
    if (temp_o != NULL)
        md3->DST_end_time = PyLong_AsLong(temp_o);
    else
        md3->DST_end_time = METADATA_DST_END_TIME_NO_ENTRY;

    temp_o = PyDict_GetItemString(md3_dict,"GMT_offset");
    if (temp_o != NULL)
        md3->GMT_offset = PyLong_AsLong(temp_o);
    else
        md3->GMT_offset = GMT_OFFSET_NO_ENTRY;

    temp_o = PyDict_GetItemString(md3_dict,"subject_name_1");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(md3->subject_name_1, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(md3_dict,"subject_name_2");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(md3->subject_name_2, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(md3_dict,"subject_ID");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(md3->subject_ID, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(md3_dict,"recording_location");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(md3->recording_location, temp_str_bytes);
    }

    return;
}

// SEGMENT, CHANNEL, SESSION are constructed from the above when reading

/**************************  Python record struct to Mef  ****************************/

void    map_python_rh(PyObject *rh_dict, RECORD_HEADER  *rh)
{
    // Helpers
    PyObject    *temp_o;
    PyObject    *temp_UTF_str;

    si1     *temp_str_bytes;
    si8     temp_time;

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(rh_dict,"type_string");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(rh->type_string, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(rh_dict,"version_major");
    if (temp_o != NULL)
        rh->version_major = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(rh_dict,"version_minor");
    if (temp_o != NULL)
        rh->version_minor = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(rh_dict,"encryption");
    if (temp_o != NULL)
        rh->encryption = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(rh_dict,"bytes");
    if (temp_o != NULL)
        rh->bytes = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(rh_dict,"time");
    if (temp_o != NULL){
        rh->time = PyLong_AsLong(temp_o);
    }

    return;
}

void    map_python_EDFA_type(PyObject *EDFA_type_dict, MEFREC_EDFA_1_0  *r_type)
{
    // Helpers
    PyObject    *temp_o;

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(EDFA_type_dict,"duration");
    if (temp_o !=NULL)
        r_type->duration = PyLong_AsLong(temp_o);

    return;
}

void    map_python_LNTP_type(PyObject *LNTP_type_dict, MEFREC_LNTP_1_0  *r_type)
{
    // Helpers
    PyObject    *temp_o;

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(LNTP_type_dict,"length");
    if (temp_o != NULL)
        r_type->length = PyLong_AsLong(temp_o);

    return;
}

void    map_python_Siez_type(PyObject *Siez_type_dict, MEFREC_Seiz_1_0  *r_type)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;

    si1     *temp_str_bytes;

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(Siez_type_dict,"earliest_onset");
    if (temp_o != NULL)
        r_type->earliest_onset = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(Siez_type_dict,"latest_offset");
    if (temp_o != NULL)
        r_type->latest_offset = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(Siez_type_dict,"duration");
    if (temp_o != NULL)
        r_type->duration = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(Siez_type_dict,"number_of_channels");
    if (temp_o != NULL)
        r_type->number_of_channels = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(Siez_type_dict,"onset_code");
    if (temp_o != NULL)
        r_type->onset_code = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(Siez_type_dict,"marker_name_1");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->marker_name_1, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(Siez_type_dict,"marker_name_2");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->marker_name_2, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(Siez_type_dict,"annotation");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->annotation, temp_str_bytes);
    }

    return;
}

void    map_python_Siez_type_channel(PyObject *Siez_ch_type_dict, MEFREC_Seiz_1_0_CHANNEL *r_type)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;

    si1     *temp_str_bytes;

    // Assign from dict to struct
    if (temp_o = PyDict_GetItemString(Siez_ch_type_dict,"name"))
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->name, temp_str_bytes);

    if (temp_o = PyDict_GetItemString(Siez_ch_type_dict,"onset"))
        r_type->onset = PyLong_AsLong(temp_o);

    if (temp_o = PyDict_GetItemString(Siez_ch_type_dict,"offset"))
        r_type->offset = PyLong_AsLong(temp_o);

    return;
}   

void    map_python_CSti_type(PyObject *CSti_type_dict, MEFREC_CSti_1_0  *r_type)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;

    si1     *temp_str_bytes;

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(CSti_type_dict,"task_type");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->task_type, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(CSti_type_dict,"stimulus_duration");
    if (temp_o != NULL)
        r_type->stimulus_duration = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(CSti_type_dict,"stimulus_type");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->stimulus_type, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(CSti_type_dict,"patient_response");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->patient_response, temp_str_bytes);
    }

    return;
}

void    map_python_ESti_type(PyObject *ESti_type_dict, MEFREC_ESti_1_0  *r_type)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;

    si1     *temp_str_bytes;

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(ESti_type_dict,"amplitude");
    if (temp_o != NULL)
        r_type->amplitude = PyFloat_AsDouble(temp_o);
    temp_o = PyDict_GetItemString(ESti_type_dict,"frequency");
    if (temp_o != NULL)
        r_type->frequency = PyFloat_AsDouble(temp_o);
    temp_o = PyDict_GetItemString(ESti_type_dict,"pulse_width");
    if (temp_o != NULL)
        r_type->pulse_width = PyLong_AsLong(temp_o);
    temp_o = PyDict_GetItemString(ESti_type_dict,"ampunit_code");
    if (temp_o != NULL)
        r_type->ampunit_code = PyLong_AsLong(temp_o);
    temp_o = PyDict_GetItemString(ESti_type_dict,"mode_code");
    if (temp_o != NULL)
        r_type->mode_code = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(ESti_type_dict,"waveform");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->waveform, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(ESti_type_dict,"anode");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->anode, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(ESti_type_dict,"catode");
    if (temp_o != NULL){
        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        MEF_strcpy(r_type->catode, temp_str_bytes);
    }

    return;
}

/*****************************  Mef struct to Python  *******************************/

PyObject *map_mef3_uh(UNIVERSAL_HEADER *uh)
{
    // Dictionaries
    PyObject *uh_dict;
    
    // Helper variables
    si1   temp_str[256];
    si8   long_file_time;
 
    /* UNIVERSAL HEADER */
    
    // Create output dictionary   
    uh_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");

    // Insert entries into dictionary
    if (uh->header_CRC != UNIVERSAL_HEADER_HEADER_CRC_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "header_CRC",
            Py_BuildValue("I", uh->header_CRC));
    else
        PyDict_SetItemString(uh_dict, "header_CRC",
            Py_BuildValue("s", temp_str));
    
    if (uh->body_CRC != UNIVERSAL_HEADER_BODY_CRC_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "body_CRC",
            Py_BuildValue("I", uh->body_CRC));
    else
        PyDict_SetItemString(uh_dict, "body_CRC",
            Py_BuildValue("s", temp_str));                         
                             
    if (uh->file_type_string[0])
        PyDict_SetItemString(uh_dict, "file_type_string",
            Py_BuildValue("s", uh->file_type_string));
    else
        PyDict_SetItemString(uh_dict, "file_type_string",
            Py_BuildValue("s", temp_str));
                             
    if (uh->mef_version_major != UNIVERSAL_HEADER_MEF_VERSION_MAJOR_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "mef_version_major",
            Py_BuildValue("B", uh->mef_version_major));
    else
        PyDict_SetItemString(uh_dict, "mef_version_major",
            Py_BuildValue("s", temp_str));   

    if (uh->mef_version_minor != UNIVERSAL_HEADER_MEF_VERSION_MINOR_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "mef_version_minor",
            Py_BuildValue("B", uh->mef_version_minor));
    else
        PyDict_SetItemString(uh_dict, "mef_version_minor",
            Py_BuildValue("s", temp_str));

    if (uh->byte_order_code != UNIVERSAL_HEADER_BYTE_ORDER_CODE_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "byte_order_code",
            Py_BuildValue("B", uh->byte_order_code));
    else
        PyDict_SetItemString(uh_dict, "byte_order_code",
            Py_BuildValue("s", temp_str));

    if (uh->start_time != UNIVERSAL_HEADER_START_TIME_NO_ENTRY){
        PyDict_SetItemString(uh_dict, "start_time",
            Py_BuildValue("l", uh->start_time));
    }
    else
        PyDict_SetItemString(uh_dict, "start_time",
            Py_BuildValue("s", temp_str)); 

    if (uh->end_time != UNIVERSAL_HEADER_END_TIME_NO_ENTRY){
        PyDict_SetItemString(uh_dict, "end_time",
            Py_BuildValue("l", uh->end_time));
    }
    else
        PyDict_SetItemString(uh_dict, "end_time",
            Py_BuildValue("s", temp_str));  

    if (uh->number_of_entries != UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "number_of_entries",
            Py_BuildValue("l", uh->number_of_entries));
    else
        PyDict_SetItemString(uh_dict, "number_of_entries",
            Py_BuildValue("s", temp_str));  

    if (uh->maximum_entry_size != UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "maximum_entry_size",
            Py_BuildValue("l", uh->maximum_entry_size));
    else
        PyDict_SetItemString(uh_dict, "maximum_entry_size",
            Py_BuildValue("s", temp_str)); 

    if (uh->segment_number != UNIVERSAL_HEADER_SEGMENT_NUMBER_NO_ENTRY)
        PyDict_SetItemString(uh_dict, "segment_number",
            Py_BuildValue("i", uh->segment_number));
    else
        PyDict_SetItemString(uh_dict, "segment_number",
            Py_BuildValue("s", temp_str)); 

    if (uh->channel_name[0])
        PyDict_SetItemString(uh_dict, "channel_name",
            Py_BuildValue("s", uh->channel_name));
    else
        PyDict_SetItemString(uh_dict, "channel_name",
            Py_BuildValue("s", temp_str));

    if (uh->session_name[0])
        PyDict_SetItemString(uh_dict, "session_name",
            Py_BuildValue("s", uh->session_name));
    else
        PyDict_SetItemString(uh_dict, "session_name",
            Py_BuildValue("s", temp_str));

    if (uh->anonymized_name[0])
        PyDict_SetItemString(uh_dict, "anonymized_name",
            Py_BuildValue("s", uh->anonymized_name));
    else
        PyDict_SetItemString(uh_dict, "anonymized_name",
            Py_BuildValue("s", temp_str));
     
    // The rest will be done if necessary

    return uh_dict;
}

PyObject *map_mef3_md1(METADATA_SECTION_1 *md1)
{
    // Dictionaries
    PyObject *s1_dict;
    
    // Helper variables
    si1   temp_str[256];
 
    /* SECTION 1 */
    
    // Create output dictionary   
    s1_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    // Insert entries into dictionary
    if (md1->section_2_encryption)
        PyDict_SetItemString(s1_dict, "section_2_encryption",
            Py_BuildValue("i", md1->section_2_encryption));
    else
        PyDict_SetItemString(s1_dict, "section_2_encryption",
            Py_BuildValue("s", temp_str));
    
    if (md1->section_3_encryption)
        PyDict_SetItemString(s1_dict, "section_3_encryption",
            Py_BuildValue("i", md1->section_3_encryption));
    else
        PyDict_SetItemString(s1_dict, "section_3_encryption",
            Py_BuildValue("s", temp_str));                         
                             
    if (md1->protected_region[0])
        PyDict_SetItemString(s1_dict, "protected_region",
            Py_BuildValue("s", md1->protected_region));
    else
        PyDict_SetItemString(s1_dict, "protected_region",
            Py_BuildValue("s", temp_str));
                             
    if (md1->discretionary_region[0])
        PyDict_SetItemString(s1_dict, "discretionary_region",
            Py_BuildValue("s", md1->discretionary_region));
    else
        PyDict_SetItemString(s1_dict, "discretionary_region",
            Py_BuildValue("s", temp_str));          

    return s1_dict;
}

PyObject *map_mef3_tmd2(TIME_SERIES_METADATA_SECTION_2 *tmd)
{
    // Dictionaries
    PyObject *s2_dict;
    
    // Helper variables
    si1   *time_str, temp_str[256];
    si8   long_file_time;
 
    /* SECTION 2 */
    
    // Create output dictionary   
    s2_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    // Insert entries into dictionary

    // Channel description
    if (tmd->channel_description[0])
        PyDict_SetItemString(s2_dict, "channel_description",
            Py_BuildValue("s", tmd->channel_description));
    else
        PyDict_SetItemString(s2_dict, "channel_description",
            Py_BuildValue("s", temp_str));

    // Session description                
    if (tmd->session_description[0])
        PyDict_SetItemString(s2_dict, "session_description",
            Py_BuildValue("s", tmd->session_description));
    else
        PyDict_SetItemString(s2_dict, "session_description",
            Py_BuildValue("s", temp_str));
    // Recording duration                    
    // long_file_time = (si8) (tmd->recording_duration + 500000) / 1000000;
    // time_str = ctime((time_t *) &long_file_time); time_str[24] = 0;
    if (tmd->recording_duration)
        PyDict_SetItemString(s2_dict, "recording_duration",
            Py_BuildValue("l", tmd->recording_duration));
    else
        PyDict_SetItemString(s2_dict, "recording_duration",
            Py_BuildValue("s", temp_str));
    
    // Reference description                         
    if (tmd->reference_description[0])
        PyDict_SetItemString(s2_dict, "reference_description",
            Py_BuildValue("s", tmd->reference_description));
    else
        PyDict_SetItemString(s2_dict, "reference_description",
            Py_BuildValue("s", temp_str));
    // Acquisition channel number                        
    if (tmd->acquisition_channel_number != TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "acquisition_channel_number",
            Py_BuildValue("i", tmd->acquisition_channel_number));
    else
        PyDict_SetItemString(s2_dict, "acquisition_channel_number",
            Py_BuildValue("s", temp_str));
    // Sampling frequency                         
    if (tmd->sampling_frequency != TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "sampling_frequency",
            Py_BuildValue("d", tmd->sampling_frequency));
    else
        PyDict_SetItemString(s2_dict, "sampling_frequency",
            Py_BuildValue("s", temp_str));           
    //Low frequency filter setting
    if (FLOAT_EQUAL (tmd->low_frequency_filter_setting, TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY))
       {sprintf(temp_str, "not entered");
        PyDict_SetItemString(s2_dict, "low_frequency_filter_setting",
            Py_BuildValue("s", temp_str));
        }
    else if (tmd->low_frequency_filter_setting < EPSILON)
        {
        sprintf(temp_str, "no low frequency filter");
        PyDict_SetItemString(s2_dict, "low_frequency_filter_setting",
            Py_BuildValue("s", temp_str));
        }
    else
        PyDict_SetItemString(s2_dict, "low_frequency_filter_setting",
            Py_BuildValue("d", tmd->low_frequency_filter_setting));
    // High frequency filter setting
    if (FLOAT_EQUAL (tmd->high_frequency_filter_setting, TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY))
       {sprintf(temp_str, "not entered");
        PyDict_SetItemString(s2_dict, "high_frequency_filter_setting",
            Py_BuildValue("s", temp_str));
        }
    else if (tmd->high_frequency_filter_setting < EPSILON)
        {
        sprintf(temp_str, "no high frequency filter");
        PyDict_SetItemString(s2_dict, "high_frequency_filter_setting",
            Py_BuildValue("s", temp_str));
        }
    else
        PyDict_SetItemString(s2_dict, "high_frequency_filter_setting",
            Py_BuildValue("d", tmd->high_frequency_filter_setting));
    // Notch frequency filter setting
    if (FLOAT_EQUAL (tmd->notch_filter_frequency_setting, TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY)){
        sprintf(temp_str, "not entered");
        PyDict_SetItemString(s2_dict, "notch_filter_frequency_setting",
            Py_BuildValue("s", temp_str));
        }
    else if (tmd->notch_filter_frequency_setting < EPSILON){
        sprintf(temp_str, "no notch filter");
        PyDict_SetItemString(s2_dict, "notch_filter_frequency_setting",
            Py_BuildValue("s", temp_str));
        }
    else
        PyDict_SetItemString(s2_dict, "notch_filter_frequency_setting",
            Py_BuildValue("d", tmd->notch_filter_frequency_setting));
    // AC line frequency
    if (FLOAT_EQUAL (tmd->AC_line_frequency, TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY)){
        sprintf(temp_str, "not entered");
        PyDict_SetItemString(s2_dict, "AC_line_frequency",
            Py_BuildValue("s", temp_str));
        }
    else
        PyDict_SetItemString(s2_dict, "AC_line_frequency",
            Py_BuildValue("d", tmd->AC_line_frequency));                          
    // Units conversion factor
    if (FLOAT_EQUAL (tmd->units_conversion_factor, TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY)){
        sprintf(temp_str, "not entered");
        PyDict_SetItemString(s2_dict, "units_conversion_factor",
            Py_BuildValue("s", temp_str));
        }
    else
        PyDict_SetItemString(s2_dict, "units_conversion_factor",
            Py_BuildValue("d",tmd->units_conversion_factor));
    // Units description
    if (tmd->units_description[0])
        PyDict_SetItemString(s2_dict, "units_description",
            Py_BuildValue("s", tmd->units_description));
    else
        PyDict_SetItemString(s2_dict, "units_description",
            Py_BuildValue("s", temp_str));
                         
    // TODO  - this is kind o a hack!
    // Maximum and minimum data value
    if (tmd->maximum_native_sample_value != RED_NAN)
        PyDict_SetItemString(s2_dict, "maximum_native_sample_value",
            Py_BuildValue("d", tmd->maximum_native_sample_value));
        
    else
        PyDict_SetItemString(s2_dict, "maximum_native_sample_value",
            Py_BuildValue("s", temp_str));
        

    if (tmd->minimum_native_sample_value != RED_NAN)
        PyDict_SetItemString(s2_dict, "minimum_native_sample_value",
            Py_BuildValue("d", tmd->minimum_native_sample_value));
    else
        PyDict_SetItemString(s2_dict, "minimum_native_sample_value",
            Py_BuildValue("s", temp_str));


    

    // Start sample
    if (tmd->start_sample != TIME_SERIES_METADATA_START_SAMPLE_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "start_sample",
            Py_BuildValue("k", tmd->start_sample));
    else
        PyDict_SetItemString(s2_dict, "start_sample",
            Py_BuildValue("s", temp_str));
    // Number of samples
    if (tmd->number_of_samples != TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "number_of_samples",
            Py_BuildValue("k", tmd->number_of_samples));
    else
        PyDict_SetItemString(s2_dict, "number_of_samples",
            Py_BuildValue("s", temp_str));
    // Number of blocks
    if (tmd->number_of_blocks != TIME_SERIES_METADATA_NUMBER_OF_BLOCKS_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "number_of_blocks",
            Py_BuildValue("k", tmd->number_of_blocks));
    else
        PyDict_SetItemString(s2_dict, "number_of_blocks",
            Py_BuildValue("s", temp_str));
    // Maximum block bytes
    if (tmd->maximum_block_bytes != TIME_SERIES_METADATA_MAXIMUM_BLOCK_BYTES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "maximum_block_bytes",
            Py_BuildValue("k", tmd->maximum_block_bytes));
    else
        PyDict_SetItemString(s2_dict, "maximum_block_bytes",
            Py_BuildValue("s", temp_str));
    // Maximum block samples
    if (tmd->maximum_block_samples != TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "maximum_block_samples",
            Py_BuildValue("k", tmd->maximum_block_samples));
    else
        PyDict_SetItemString(s2_dict, "maximum_block_samples",
            Py_BuildValue("s", temp_str));
    // Maximum difference bytes
    if (tmd->maximum_difference_bytes != TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "maximum_difference_bytes",
            Py_BuildValue("k", tmd->maximum_difference_bytes));
    else
        PyDict_SetItemString(s2_dict, "maximum_difference_bytes",
            Py_BuildValue("s", temp_str));
    // Block interval
    if (tmd->block_interval != TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "block_interval",
            Py_BuildValue("i", tmd->block_interval));
    else
        PyDict_SetItemString(s2_dict, "block_interval",
            Py_BuildValue("s", temp_str));
    // Number of discontinuities
    if (tmd->number_of_discontinuities != TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "number_of_discontinuities",
            Py_BuildValue("i", tmd->number_of_discontinuities));
    else
        PyDict_SetItemString(s2_dict, "number_of_discontinuities",
            Py_BuildValue("s", temp_str));
    // Maximum contiguous blocks
    if (tmd->maximum_contiguous_blocks != TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCKS_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "maximum_contiguous_blocks",
            Py_BuildValue("i", tmd->maximum_contiguous_blocks));
    else
        PyDict_SetItemString(s2_dict, "maximum_contiguous_blocks",
            Py_BuildValue("s", temp_str));
    // Maximum contiguous block bytes
    if (tmd->maximum_contiguous_block_bytes != TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "maximum_contiguous_block_bytes",
            Py_BuildValue("i", tmd->maximum_contiguous_block_bytes));
    else
        PyDict_SetItemString(s2_dict, "maximum_contiguous_block_bytes",
            Py_BuildValue("s", temp_str));
    // Maximum contiguous samples
    if (tmd->maximum_contiguous_samples != TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "maximum_contiguous_samples",
            Py_BuildValue("i", tmd->maximum_contiguous_samples));
    else
        PyDict_SetItemString(s2_dict, "maximum_contiguous_samples",
            Py_BuildValue("s", temp_str));
    // Protected region
    if (tmd->protected_region[0])
        PyDict_SetItemString(s2_dict, "protected_region",
            Py_BuildValue("s", tmd->protected_region));
    else
        PyDict_SetItemString(s2_dict, "protected_region",
            Py_BuildValue("s", temp_str));
    // Discretionary region
    if (tmd->discretionary_region[0])
        PyDict_SetItemString(s2_dict, "discretionary_region",
            Py_BuildValue("s", tmd->discretionary_region));
    else
        PyDict_SetItemString(s2_dict, "discretionary_region",
            Py_BuildValue("s", temp_str));

    return s2_dict;
}

PyObject *map_mef3_vmd2(VIDEO_METADATA_SECTION_2 *vmd)
{
    // Dictionaries
    PyObject *s2_dict;
    
    // Helper variables
    si1   *time_str, temp_str[256];
    si8   long_file_time;
 
    /* SECTION 2 */
    
    // Create output dictionary   
    s2_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");

    // Channel description
    if (vmd->channel_description[0])
        PyDict_SetItemString(s2_dict, "channel_description",
            Py_BuildValue("s", vmd->channel_description));
    else
        PyDict_SetItemString(s2_dict, "channel_description",
            Py_BuildValue("s", temp_str));
    // Session description                
    if (vmd->session_description[0])
        PyDict_SetItemString(s2_dict, "session_description",
            Py_BuildValue("s", vmd->session_description));
    else
        PyDict_SetItemString(s2_dict, "session_description",
            Py_BuildValue("s", temp_str));
    // Recording duration                    
    // long_file_time = (si8) (vmd->recording_duration + 500000) / 1000000;
    // time_str = ctime((time_t *) &long_file_time); time_str[24] = 0;
    if (vmd->recording_duration)
        PyDict_SetItemString(s2_dict, "recording_duration",
            Py_BuildValue("k", vmd->recording_duration));
    else
        PyDict_SetItemString(s2_dict, "recording_duration",
            Py_BuildValue("s", temp_str));
    
    // Horizontal resolution
    if (vmd->horizontal_resolution != VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "horizontal_resolution",
            Py_BuildValue("i", vmd->horizontal_resolution));
    else
        PyDict_SetItemString(s2_dict, "horizontal_resolution",
            Py_BuildValue("s", temp_str));
                             
    // Vertical resolution
    if (vmd->vertical_resolution != VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "vertical_resolution",
            Py_BuildValue("i", vmd->vertical_resolution));
    else
        PyDict_SetItemString(s2_dict, "vertical_resolution",
            Py_BuildValue("s", temp_str));
                             
    // Frame rate
    if (vmd->frame_rate != VIDEO_METADATA_FRAME_RATE_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "frame_rate",
            Py_BuildValue("d", vmd->frame_rate));
    else
        PyDict_SetItemString(s2_dict, "frame_rate",
            Py_BuildValue("s", temp_str));
                             
    // Number of clips
    if (vmd->number_of_clips != VIDEO_METADATA_NUMBER_OF_CLIPS_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "number_of_clips",
            Py_BuildValue("i", vmd->number_of_clips));
    else
        PyDict_SetItemString(s2_dict, "number_of_clips",
            Py_BuildValue("s", temp_str));
                             
    // Maximum clip bytes
    if (vmd->maximum_clip_bytes != VIDEO_METADATA_MAXIMUM_CLIP_BYTES_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "maximum_clip_bytes",
            Py_BuildValue("i", vmd->maximum_clip_bytes));
    else
        PyDict_SetItemString(s2_dict, "maximum_clip_bytes",
            Py_BuildValue("s", temp_str));
              
    // Video format
    if (vmd->video_format[0])
        PyDict_SetItemString(s2_dict, "video_format",
            Py_BuildValue("s", vmd->video_format));
    else
        PyDict_SetItemString(s2_dict, "video_format",
            Py_BuildValue("s", temp_str));
                             
    // Video file CRC
    if (vmd->video_file_CRC != VIDEO_METADATA_VIDEO_FILE_CRC_NO_ENTRY)
        PyDict_SetItemString(s2_dict, "video_file_CRC",
            Py_BuildValue("i", vmd->video_file_CRC));
    else
        PyDict_SetItemString(s2_dict, "video_file_CRC",
            Py_BuildValue("s", temp_str));
                             
    // Protected region
    if (vmd->protected_region[0])
        PyDict_SetItemString(s2_dict, "protected_region",
            Py_BuildValue("s", vmd->protected_region));
    else
        PyDict_SetItemString(s2_dict, "protected_region",
            Py_BuildValue("s", temp_str));
    // Discretionary region
    if (vmd->discretionary_region[0])
        PyDict_SetItemString(s2_dict, "discretionary_region",
            Py_BuildValue("s", vmd->discretionary_region));
    else
        PyDict_SetItemString(s2_dict, "discretionary_region",
            Py_BuildValue("s", temp_str));
                                 
    return s2_dict;
}

PyObject *map_mef3_md3(METADATA_SECTION_3 *md3)
{
    // Dictionaries
    PyObject *s3_dict;
    
    // Helper variables
    si1   *time_str, temp_str[256];
    si8   long_file_time;
    
    /* SECTION 3 */

    // Create section 3 dictionary
    s3_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    // Insert entries into dictionary
    
    // Recording time offset
    // long_file_time = (si8) (md3->recording_time_offset + 500000) / 1000000;
    // time_str = ctime((time_t *) &long_file_time); time_str[24] = 0;
    if (md3->recording_time_offset != METADATA_RECORDING_TIME_OFFSET_NO_ENTRY)
        PyDict_SetItemString(s3_dict, "recording_time_offset",
            Py_BuildValue("l", md3->recording_time_offset));
    else
        PyDict_SetItemString(s3_dict, "recording_time_offset",
            Py_BuildValue("s", temp_str));
                             
    // DST start time
    // long_file_time = (si8) (md3->DST_start_time + 500000) / 1000000;
    // time_str = ctime((time_t *) &long_file_time); time_str[24] = 0;
    if (md3->DST_start_time != METADATA_DST_START_TIME_NO_ENTRY){
        PyDict_SetItemString(s3_dict, "DST_start_time",
            Py_BuildValue("l", md3->DST_start_time));
    }
    else
        PyDict_SetItemString(s3_dict, "DST_start_time",
            Py_BuildValue("s", temp_str));
                             
    // DST end time
    // long_file_time = (si8) (md3->DST_end_time + 500000) / 1000000;
    // time_str = ctime((time_t *) &long_file_time); time_str[24] = 0;
    if (md3->DST_end_time != METADATA_DST_END_TIME_NO_ENTRY){
        PyDict_SetItemString(s3_dict, "DST_end_time",
            Py_BuildValue("l", md3->DST_end_time));
    }
    else
        PyDict_SetItemString(s3_dict, "DST_end_time",
            Py_BuildValue("s", temp_str));
                             
    // GMT offset
    // long_file_time = (si8) (md3->GMT_offset + 500000) / 1000000;
    // time_str = ctime((time_t *) &long_file_time); time_str[24] = 0;
    if (md3->GMT_offset != GMT_OFFSET_NO_ENTRY)
        PyDict_SetItemString(s3_dict, "GMT_offset",
            Py_BuildValue("i", md3->GMT_offset));
    else
        PyDict_SetItemString(s3_dict, "GMT_offset",
            Py_BuildValue("s", temp_str));
                             
    // Subject name 1
    if (md3->subject_name_1[0])
        PyDict_SetItemString(s3_dict, "subject_name_1",
            Py_BuildValue("s", md3->subject_name_1));
    else
        PyDict_SetItemString(s3_dict, "subject_name_1",
            Py_BuildValue("s", temp_str));
                             
    // Subject name 2
    if (md3->subject_name_2[0])
        PyDict_SetItemString(s3_dict, "subject_name_2",
            Py_BuildValue("s", md3->subject_name_2));
    else
        PyDict_SetItemString(s3_dict, "subject_name_2",
            Py_BuildValue("s", temp_str));
                             
    // Subject ID
    if (md3->subject_ID[0])
        PyDict_SetItemString(s3_dict, "subject_ID",
            Py_BuildValue("s", md3->subject_ID));
    else
        PyDict_SetItemString(s3_dict, "subject_ID",
            Py_BuildValue("s", temp_str));
                             
    // Recording location
    if (md3->recording_location[0])
        PyDict_SetItemString(s3_dict, "recording_location",
            Py_BuildValue("s", md3->recording_location));
    else
        PyDict_SetItemString(s3_dict, "recording_location",
            Py_BuildValue("s", temp_str));
    // Protected region
    // if (md3->protected_region[0])
    //     PyDict_SetItemString(s3_dict, "protected_region",
    //         Py_BuildValue("s", md3->protected_region));
    // else
    //     PyDict_SetItemString(s3_dict, "protected_region",
    //         Py_BuildValue("s", temp_str));
    // // Discretionary region
    // if (md3->discretionary_region[0])
    //     PyDict_SetItemString(s3_dict, "discretionary_region",
    //         Py_BuildValue("s", md3->discretionary_region));
    // else
    //     PyDict_SetItemString(s3_dict, "discretionary_region",
    //         Py_BuildValue("s", temp_str));

    return s3_dict; 
}

PyObject *map_mef3_ti(TIME_SERIES_INDEX *ti)
{
    // Dictionaries
    PyObject *ti_dict;
    
    // Helper variables
    si1   temp_str[256];
 
    /* TIME SERIES INDEX */
    
    // Create output dictionary   
    ti_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");

    // Insert entries into dictionary                    
    if (ti->file_offset != TIME_SERIES_INDEX_FILE_OFFSET_NO_ENTRY)
        PyDict_SetItemString(ti_dict, "file_offset",
            Py_BuildValue("l", ti->file_offset));
    else
        PyDict_SetItemString(ti_dict, "file_offset",
            Py_BuildValue("s", temp_str));

    if (ti->start_time != TIME_SERIES_INDEX_START_TIME_NO_ENTRY)
        PyDict_SetItemString(ti_dict, "start_time",
            Py_BuildValue("l", ti->start_time));
    else
        PyDict_SetItemString(ti_dict, "start_time",
            Py_BuildValue("s", temp_str));

    if (ti->start_sample != TIME_SERIES_INDEX_START_SAMPLE_NO_ENTRY)
        PyDict_SetItemString(ti_dict, "start_sample",
            Py_BuildValue("l", ti->start_sample));
    else
        PyDict_SetItemString(ti_dict, "start_sample",
            Py_BuildValue("s", temp_str));

    if (ti->number_of_samples != TIME_SERIES_INDEX_NUMBER_OF_SAMPLES_NO_ENTRY)
        PyDict_SetItemString(ti_dict, "number_of_samples",
            Py_BuildValue("I", ti->number_of_samples));
    else
        PyDict_SetItemString(ti_dict, "number_of_samples",
            Py_BuildValue("s", temp_str));    

    if (ti->block_bytes != TIME_SERIES_INDEX_BLOCK_BYTES_NO_ENTRY)
        PyDict_SetItemString(ti_dict, "block_bytes",
            Py_BuildValue("I", ti->block_bytes));
    else
        PyDict_SetItemString(ti_dict, "block_bytes",
            Py_BuildValue("s", temp_str));   

    if (ti->maximum_sample_value != TIME_SERIES_INDEX_MAXIMUM_SAMPLE_VALUE_NO_ENTRY)
        PyDict_SetItemString(ti_dict, "maximum_sample_value",
            Py_BuildValue("i", ti->maximum_sample_value));
    else
        PyDict_SetItemString(ti_dict, "maximum_sample_value",
            Py_BuildValue("s", temp_str)); 

    if (ti->minimum_sample_value != TIME_SERIES_INDEX_MINIMUM_SAMPLE_VALUE_NO_ENTRY)
        PyDict_SetItemString(ti_dict, "minimum_sample_value",
            Py_BuildValue("i", ti->minimum_sample_value));
    else
        PyDict_SetItemString(ti_dict, "minimum_sample_value",
            Py_BuildValue("s", temp_str));  

    // The rest is to be done if neccessary

    return ti_dict;
}

PyObject *map_mef3_vi(VIDEO_INDEX *vi)
{
    // Dictionaries
    PyObject *vi_dict;
    
    // Helper variables
    si1   temp_str[256];
 
    /* TIME SERIES INDEX */
    
    // Create output dictionary   
    vi_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");

    // Insert entries into dictionary                    
    if (vi->start_time != VIDEO_INDEX_START_TIME_NO_ENTRY)
        PyDict_SetItemString(vi_dict, "start_time",
            Py_BuildValue("l", vi->start_time));
    else
        PyDict_SetItemString(vi_dict, "start_time",
            Py_BuildValue("s", temp_str));

    if (vi->end_time != VIDEO_INDEX_END_TIME_NO_ENTRY)
        PyDict_SetItemString(vi_dict, "end_time",
            Py_BuildValue("l", vi->end_time));
    else
        PyDict_SetItemString(vi_dict, "end_time",
            Py_BuildValue("s", temp_str));

    if (vi->start_frame != VIDEO_INDEX_START_FRAME_NO_ENTRY)
        PyDict_SetItemString(vi_dict, "start_frame",
            Py_BuildValue("I", vi->start_frame));
    else
        PyDict_SetItemString(vi_dict, "start_frame",
            Py_BuildValue("s", temp_str));

    if (vi->end_frame != VIDEO_INDEX_END_FRAME_NO_ENTRY)
        PyDict_SetItemString(vi_dict, "end_frame",
            Py_BuildValue("I", vi->end_frame));
    else
        PyDict_SetItemString(vi_dict, "end_frame",
            Py_BuildValue("s", temp_str));

    if (vi->file_offset != VIDEO_INDEX_FILE_OFFSET_NO_ENTRY)
        PyDict_SetItemString(vi_dict, "file_offset",
            Py_BuildValue("l", vi->file_offset));
    else
        PyDict_SetItemString(vi_dict, "file_offset",
            Py_BuildValue("s", temp_str));

    if (vi->clip_bytes != VIDEO_INDEX_CLIP_BYTES_NO_ENTRY)
        PyDict_SetItemString(vi_dict, "clip_bytes",
            Py_BuildValue("l", vi->clip_bytes));
    else
        PyDict_SetItemString(vi_dict, "clip_bytes",
            Py_BuildValue("s", temp_str));

    // The rest if neccessary

    return vi_dict;
}

PyObject *map_mef3_segment(SEGMENT *segment)
{
    // Dictionaries
    PyObject *metadata_dict;
    PyObject *records_dict;
    PyObject *spec_dict;
    PyObject *uh_dict;
    PyObject *uhs_dict;
    PyObject *s1_dict;
    PyObject *s2_dict;
    PyObject *s3_dict;
    PyObject *sdi_dict;
    PyObject *idx_list;
    
    // Helper variables
    si1   *time_str, temp_str[256];
    ui4   i;
    si8   long_file_time, number_of_entries;
    
    METADATA_SECTION_1      *md1;
        TIME_SERIES_METADATA_SECTION_2  *tmd2;
        VIDEO_METADATA_SECTION_2    *vmd2;
    METADATA_SECTION_3      *md3;
    TIME_SERIES_INDEX       *tsi;
    VIDEO_INDEX       *vi;

    // Create python output dictionary
    metadata_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    /* SEGMENT SPECIFIC */
    
    // Create specific dictionary
    PyDict_SetItemString(metadata_dict, "segment_specific_metadata", PyDict_New()); 
    spec_dict = PyDict_GetItemString(metadata_dict, "segment_specific_metadata");
    
    // Insert entries into dictionary
    
    // Segment name             
    if (segment->name[0])
        PyDict_SetItemString(spec_dict, "segment_name",
            Py_BuildValue("s", segment->name));
    else
        PyDict_SetItemString(spec_dict, "segment_name",
            Py_BuildValue("s", temp_str));
                             
    // Read segment records if present and add it to metadata
    if (segment->record_indices_fps != NULL & segment->record_data_fps != NULL){
        records_dict = map_mef3_records(segment->record_indices_fps, segment->record_data_fps);
        PyDict_SetItemString(metadata_dict, "records_info", records_dict);
    }

    // Assign pointers for reading metadata and read universal headers
    md1 = segment->metadata_fps->metadata.section_1;
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            tmd2 = segment->metadata_fps->metadata.time_series_section_2;
            break;
        case VIDEO_CHANNEL_TYPE:
            vmd2 = segment->metadata_fps->metadata.video_section_2; 
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }      
    md3 = segment->metadata_fps->metadata.section_3;

    // Create section 1 dictionary
    s1_dict = map_mef3_md1(md1);
    PyDict_SetItemString(metadata_dict, "section_1", s1_dict);
    
    // Create section 2 dictionary
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            s2_dict = map_mef3_tmd2(tmd2);
            break;
        case VIDEO_CHANNEL_TYPE:
            s2_dict = map_mef3_vmd2(vmd2);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }

    PyDict_SetItemString(metadata_dict, "section_2", s2_dict); 
    
    // Create section 3 dictionary
    s3_dict = map_mef3_md3(md3);
    PyDict_SetItemString(metadata_dict, "section_3", s3_dict);

    // TODO - this should be a list - there is more indices in indices file!!!!

    // Create indices dictionary
    
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            number_of_entries = segment->time_series_indices_fps->universal_header->number_of_entries;
            tsi = segment->time_series_indices_fps->time_series_indices;
            idx_list = PyList_New(number_of_entries);
            for(i = 0; i < number_of_entries; i++){
                sdi_dict = map_mef3_ti(tsi);
                PyList_SET_ITEM(idx_list, i, sdi_dict);
                tsi++;
            }
            break;
        case VIDEO_CHANNEL_TYPE:
            number_of_entries = segment->video_indices_fps->universal_header->number_of_entries;
            vi = segment->video_indices_fps->video_indices;
            idx_list = PyList_New(number_of_entries);
            for(i = 0; i < number_of_entries; i++){
                sdi_dict = map_mef3_vi(vi);
                PyList_SET_ITEM(idx_list, i, sdi_dict);
                vi++;
            }
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }
    PyDict_SetItemString(metadata_dict, "indices", idx_list);

    // Get universal headers
    uhs_dict = PyDict_New();

    // Metadata
    uh_dict = map_mef3_uh(segment->metadata_fps->universal_header);
    PyDict_SetItemString(uhs_dict, "metadata", uh_dict);

    // Data an indices universal headers
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            uh_dict = map_mef3_uh(segment->time_series_data_fps->universal_header);
            PyDict_SetItemString(uhs_dict, "time_series_data", uh_dict);
            uh_dict = map_mef3_uh(segment->time_series_indices_fps->universal_header);
            PyDict_SetItemString(uhs_dict, "time_series_indices", uh_dict);
            break;
        case VIDEO_CHANNEL_TYPE:
            uh_dict = map_mef3_uh(segment->video_indices_fps->universal_header);
            PyDict_SetItemString(uhs_dict, "video_indices", uh_dict);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }

    PyDict_SetItemString(metadata_dict, "universal_headers", uhs_dict);

    return metadata_dict;
}

PyObject *map_mef3_channel(CHANNEL *channel) // This funtion also loops through segments
{
    // Dictionaries
    PyObject *metadata_dict;
    PyObject *records_dict;
    PyObject *spec_dict;
    PyObject *s1_dict;
    PyObject *s2_dict;
    PyObject *s3_dict;
    PyObject *segment_dict;
    PyObject *segments_dict;
    
    // Helper variables
    si1   temp_str[256];
    si4   i;
    
    // Method
    SEGMENT *segment;
    
    METADATA_SECTION_1		*md1;
        TIME_SERIES_METADATA_SECTION_2	*tmd2;
        VIDEO_METADATA_SECTION_2	*vmd2;
    METADATA_SECTION_3		*md3;

    // Create python output dictionary
    metadata_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    /* CHANNEL SPECIFIC */
    
    // Create specific dictionary
    PyDict_SetItemString(metadata_dict, "channel_specific_metadata", PyDict_New()); 
    spec_dict = PyDict_GetItemString(metadata_dict, "channel_specific_metadata");
    
    // Insert entries into dictionary
    
    // Earliest start time
    if (channel->earliest_start_time != UNIVERSAL_HEADER_START_TIME_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "earliest_start_time",
            Py_BuildValue("l", channel->earliest_start_time));
    else
        PyDict_SetItemString(spec_dict, "earliest_start_time",
            Py_BuildValue("s", temp_str));
    // Latest end time                     
    if (channel->latest_end_time != UNIVERSAL_HEADER_END_TIME_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "latest_end_time",
            Py_BuildValue("l", channel->latest_end_time));
    else
        PyDict_SetItemString(spec_dict, "latest_end_time",
            Py_BuildValue("s", temp_str));
    // Anonymized name                     
    if (channel->anonymized_name[0])
        PyDict_SetItemString(spec_dict, "anonymized_name",
            Py_BuildValue("s", channel->anonymized_name));
    else
        PyDict_SetItemString(spec_dict, "anonymized_name",
            Py_BuildValue("s", temp_str));
    // Maximum number of records
    if (channel->maximum_number_of_records != UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "maximum_number_of_records",
            Py_BuildValue("l", channel->maximum_number_of_records));
    else
        PyDict_SetItemString(spec_dict, "maximum_number_of_records",
            Py_BuildValue("s", temp_str));
    // Maximum record bytes                  
    if (channel->maximum_record_bytes != UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "maximum_record_bytes",
            Py_BuildValue("l", channel->maximum_record_bytes));
    else
        PyDict_SetItemString(spec_dict, "maximum_record_bytes",
            Py_BuildValue("s", temp_str));
    // Channel name             
    if (channel->name[0])
        PyDict_SetItemString(spec_dict, "channel_name",
            Py_BuildValue("s", channel->name));
    else
        PyDict_SetItemString(spec_dict, "channel_name",
            Py_BuildValue("s", temp_str));
    
    // Read channel records if present and add it to metadata
    if (channel->record_indices_fps != NULL & channel->record_data_fps != NULL){
        records_dict = map_mef3_records(channel->record_indices_fps, channel->record_data_fps);
        PyDict_SetItemString(metadata_dict, "records", records_dict);
    }

    // Assign pointers for reading metadata
    md1 = channel->metadata.section_1;
    switch (channel->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            tmd2 = channel->metadata.time_series_section_2;
            break;
        case VIDEO_CHANNEL_TYPE:
            vmd2 = channel->metadata.video_section_2;  
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }     
    md3 = channel->metadata.section_3;
    
    // Create section 1 dictionary
    s1_dict = map_mef3_md1(md1);
    PyDict_SetItemString(metadata_dict, "section_1", s1_dict);
    
    // Create section 2 dictionary
    switch (channel->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            s2_dict = map_mef3_tmd2(tmd2);
            break;
        case VIDEO_CHANNEL_TYPE:
            s2_dict = map_mef3_vmd2(vmd2);
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }

    PyDict_SetItemString(metadata_dict, "section_2", s2_dict); 
    
    // Create section 3 dictionary
    s3_dict = map_mef3_md3(md3);
    PyDict_SetItemString(metadata_dict, "section_3", s3_dict);
    
    // Loop over segments
    for (i = 0; i < channel->number_of_segments; ++i){
        if (i == 0){
            PyDict_SetItemString(metadata_dict, "segments", PyDict_New()); 
            segments_dict = PyDict_GetItemString(metadata_dict, "segments");
        }
        // Get the channel pointer
        segment = channel->segments + i;
        // Map the channel
        segment_dict = map_mef3_segment(segment);
        // Put into ditionary
        PyDict_SetItemString(segments_dict, segment->name, segment_dict); 
    }

    return metadata_dict;
}

PyObject *map_mef3_session(SESSION *session) // This funtion also loops through channels
{
    // Dictionaries
    PyObject *metadata_dict;
    PyObject *records_dict;
    PyObject *spec_dict;

    PyObject *ts_metadata;
    PyObject *ts_dict;
    PyObject *s1_ts_dict;
    PyObject *s2_ts_dict;
    PyObject *s3_ts_dict;

    PyObject *v_metadata;
    PyObject *v_dict;
    PyObject *s1_v_dict;
    PyObject *s2_v_dict;
    PyObject *s3_v_dict;

    PyObject *channel_dict;
    
    // Helper variables
    si1   *time_str, temp_str[256];
    si4   i;
    si8   long_file_time;
    
    // Method 
    CHANNEL *channel;

    METADATA_SECTION_1      *md1;
        TIME_SERIES_METADATA_SECTION_2  *tmd2;
        VIDEO_METADATA_SECTION_2    *vmd2;
    METADATA_SECTION_3      *md3;

    // Create python output dictionary
    metadata_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    /* SESSION SPECIFIC */
    
    // Create specific dictionary
    PyDict_SetItemString(metadata_dict, "session_specific_metadata", PyDict_New()); 
    spec_dict = PyDict_GetItemString(metadata_dict, "session_specific_metadata");
    
    // Insert entries into dictionary
    
    // Earliest start time
    if (session->earliest_start_time != UNIVERSAL_HEADER_START_TIME_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "earliest_start_time",
            Py_BuildValue("l", session->earliest_start_time));
    else
        PyDict_SetItemString(spec_dict, "earliest_start_time",
            Py_BuildValue("s", temp_str));
    // Latest end time                     
    if (session->latest_end_time != UNIVERSAL_HEADER_END_TIME_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "latest_end_time",
            Py_BuildValue("l", session->latest_end_time));
    else
        PyDict_SetItemString(spec_dict, "latest_end_time",
            Py_BuildValue("s", temp_str));
    // Anonymized name                     
    if (session->anonymized_name[0])
        PyDict_SetItemString(spec_dict, "anonymized_name",
            Py_BuildValue("s", session->anonymized_name));
    else
        PyDict_SetItemString(spec_dict, "anonymized_name",
            Py_BuildValue("s", temp_str));
    // Maximum number of records
    if (session->maximum_number_of_records != UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "maximum_number_of_records",
            Py_BuildValue("l", session->maximum_number_of_records));
    else
        PyDict_SetItemString(spec_dict, "maximum_number_of_records",
            Py_BuildValue("s", temp_str));
    // Maximum record bytes                  
    if (session->maximum_record_bytes != UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_NO_ENTRY)
        PyDict_SetItemString(spec_dict, "maximum_record_bytes",
            Py_BuildValue("l", session->maximum_record_bytes));
    else
        PyDict_SetItemString(spec_dict, "maximum_record_bytes",
            Py_BuildValue("s", temp_str));
    // Session name             
    if (session->name[0])
        PyDict_SetItemString(spec_dict, "session_name",
            Py_BuildValue("s", session->name));
    else
        PyDict_SetItemString(spec_dict, "session_name",
            Py_BuildValue("s", temp_str));
    // Number of time series channels            
    if (session->number_of_time_series_channels != 0)
        PyDict_SetItemString(spec_dict, "number_of_times_series_channles",
            Py_BuildValue("i", session->number_of_time_series_channels));
    else
        PyDict_SetItemString(spec_dict, "number_of_times_series_channles",
            Py_BuildValue("i", 0));
    // Number of video channels            
    if (session->number_of_video_channels != 0)
        PyDict_SetItemString(spec_dict, "number_of_video_channels",
            Py_BuildValue("i", session->number_of_video_channels));
    else
        PyDict_SetItemString(spec_dict, "number_of_video_channels",
            Py_BuildValue("i", 0));                    

    // Read session records if present and add it to metadata
    if (session->record_indices_fps != NULL & session->record_data_fps != NULL){
        records_dict = map_mef3_records(session->record_indices_fps, session->record_data_fps);
        PyDict_SetItemString(metadata_dict, "records", records_dict);
    }

    // Get time series metadata
    if (session->number_of_time_series_channels > 0)
    {
        // Create dictionary
        PyDict_SetItemString(metadata_dict, "time_series_metadata", PyDict_New()); 
        ts_metadata = PyDict_GetItemString(metadata_dict, "time_series_metadata");

        // Assign pointers for reading time series metadata
        md1 = session->time_series_metadata.section_1;
        tmd2 = session->time_series_metadata.time_series_section_2;
        md3 = session->time_series_metadata.section_3;
        
        // Create section time series 1 dictionary
        s1_ts_dict = map_mef3_md1(md1);
        PyDict_SetItemString(ts_metadata, "section_1", s1_ts_dict);
        
        // Create section time series 2 dictionary
        s2_ts_dict = map_mef3_tmd2(tmd2);
        PyDict_SetItemString(ts_metadata, "section_2", s2_ts_dict); 
        
        // Create section time series 3 dictionary
        s3_ts_dict = map_mef3_md3(md3);
        PyDict_SetItemString(ts_metadata, "section_3", s3_ts_dict);
    }

    // Get video metadata
    if (session->number_of_video_channels > 0)
    {
        // Create dictionary
        PyDict_SetItemString(metadata_dict, "video_metadata", PyDict_New()); 
        v_metadata = PyDict_GetItemString(metadata_dict, "video_metadata");

        // Assign pointers for reading video metadata
        md1 = session->video_metadata.section_1;
        vmd2 = session->video_metadata.video_section_2;
        md3 = session->video_metadata.section_3;
        
        // Create section video 1 dictionary
        s1_v_dict = map_mef3_md1(md1);
        PyDict_SetItemString(v_metadata, "section_1", s1_v_dict);
        
        // Create section video 2 dictionary
        s2_v_dict = map_mef3_vmd2(vmd2);
        PyDict_SetItemString(v_metadata, "section_2", s2_v_dict); 
        
        // Create section tvideo 3 dictionary
        s3_v_dict = map_mef3_md3(md3);
        PyDict_SetItemString(v_metadata, "section_3", s3_v_dict);
    }

    // Loop over time series channels         
    for (i = 0; i < session->number_of_time_series_channels; ++i){
        if (i == 0){
            PyDict_SetItemString(metadata_dict, "ts_channels", PyDict_New()); 
            ts_dict = PyDict_GetItemString(metadata_dict, "ts_channels");
        }
        // Get the channel pointer
        channel = session->time_series_channels + i;
        // Map the channel
        channel_dict = map_mef3_channel(channel);
        // Put into ditionary
        PyDict_SetItemString(ts_dict, channel->name, channel_dict); 

    }
    // Loop over video channels
    for (i = 0; i < session->number_of_video_channels; ++i){
        if (i == 0){
            PyDict_SetItemString(metadata_dict, "v_channels", PyDict_New()); 
            v_dict = PyDict_GetItemString(metadata_dict, "v_channels");
        }
        // Get the channel pointer
        channel = session->video_channels + i;
        // Map the channel
        channel_dict = map_mef3_channel(channel);
        // Put into ditionary
        PyDict_SetItemString(v_dict, channel->name, channel_dict); 
    }

    return metadata_dict;
}

/**************************  Mef record struct to Python  ****************************/

PyObject *map_mef3_records(FILE_PROCESSING_STRUCT *ri_fps, FILE_PROCESSING_STRUCT *rd_fps) // Runs through records and puts them in a dictionary
{

    // Ditionary
    PyObject    *record_info_dict;
    PyObject    *all_record_list;
    PyObject    *record_dict;
    PyObject    *uhs_dict;
    PyObject    *uh_dict;

    void    *rd;
    si4     i;
    si8     number_of_records;

    RECORD_HEADER   *rh;
    RECORD_INDEX    *ri;
 
    // Set up output dictionary
    record_info_dict = PyDict_New();

    // Get universal headers
    uhs_dict = PyDict_New();
    uh_dict = map_mef3_uh(ri_fps->universal_header);
    PyDict_SetItemString(uhs_dict, "record_indices", uh_dict);
    uh_dict = map_mef3_uh(rd_fps->universal_header);
    PyDict_SetItemString(uhs_dict, "record_data", uh_dict);
    PyDict_SetItemString(record_info_dict, "universal_headers", uhs_dict);

    // Create list for records
    number_of_records = ri_fps->universal_header->number_of_entries;
    all_record_list = PyList_New(number_of_records);

    // First entry
    ri = ri_fps->record_indices; // This is unneccesary now but can be used to filter read records. Will be done in python now.
    rd = rd_fps->raw_data + UNIVERSAL_HEADER_BYTES;
    

    for (i=0; i < number_of_records; ++i){

        rh = (RECORD_HEADER *) rd;

        record_dict = map_mef3_rh(rh);
        PyList_SET_ITEM(all_record_list, i, record_dict); // ASK Matt / Dan. Could also be type_string, 

        rd += (RECORD_HEADER_BYTES + rh->bytes);
        
    }

    PyDict_SetItemString(record_info_dict, "records", all_record_list);

    return record_info_dict;
}

PyObject *map_mef3_rh(RECORD_HEADER *rh)
{
    // Dictionaries
    PyObject *rh_dict;
    PyObject *type_dict;

    // Helper variables
    si1   *time_str, temp_str[256];
    si8   long_file_time;
    ui4   type_code;
 
    /* RECORD HEADER */
    
    // Create output dictionary   
    rh_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    // Insert entries into dictionary
    if (rh->record_CRC != RECORD_HEADER_RECORD_CRC_NO_ENTRY)
        PyDict_SetItemString(rh_dict, "record_CRC",
            Py_BuildValue("I", rh->record_CRC));
    else
        PyDict_SetItemString(rh_dict, "record_CRC",
            Py_BuildValue("s", temp_str));                       
                             
    if (rh->type_string[0])
        PyDict_SetItemString(rh_dict, "type_string",
            Py_BuildValue("s", rh->type_string));
    else
        PyDict_SetItemString(rh_dict, "type_string",
            Py_BuildValue("s", temp_str));
                             
    if (rh->version_major != RECORD_HEADER_VERSION_MAJOR_NO_ENTRY)
        PyDict_SetItemString(rh_dict, "version_major",
            Py_BuildValue("B", rh->version_major));
    else
        PyDict_SetItemString(rh_dict, "version_major",
            Py_BuildValue("s", temp_str));   

    if (rh->version_minor != RECORD_HEADER_VERSION_MINOR_NO_ENTRY)
        PyDict_SetItemString(rh_dict, "version_minor",
            Py_BuildValue("B", rh->version_minor));
    else
        PyDict_SetItemString(rh_dict, "version_minor",
            Py_BuildValue("s", temp_str));

    if (rh->encryption != NULL)
        PyDict_SetItemString(rh_dict, "encryption",
            Py_BuildValue("b", rh->encryption));
    else
        PyDict_SetItemString(rh_dict, "encryption",
            Py_BuildValue("s", temp_str));

    if (rh->bytes != RECORD_HEADER_BYTES_NO_ENTRY)
        PyDict_SetItemString(rh_dict, "bytes",
            Py_BuildValue("I", rh->bytes));
    else
        PyDict_SetItemString(rh_dict, "bytes",
            Py_BuildValue("s", temp_str));

    if (rh->time != RECORD_HEADER_TIME_NO_ENTRY){
        PyDict_SetItemString(rh_dict, "time",
            Py_BuildValue("l", rh->time));
    }
    else
        PyDict_SetItemString(rh_dict, "time",
            Py_BuildValue("s", temp_str));

    type_code = *((ui4 *) rh->type_string);

    switch (type_code) {
        case MEFREC_Note_TYPE_CODE:
            type_dict = map_mef3_Note_type(rh);
            break;
        case MEFREC_EDFA_TYPE_CODE:
            type_dict = map_mef3_EDFA_type(rh);
            break;
        case MEFREC_LNTP_TYPE_CODE:
            type_dict = map_mef3_LNTP_type(rh);
            break;
        case MEFREC_Seiz_TYPE_CODE:
            type_dict = map_mef3_Seiz_type(rh);
            break;
        case MEFREC_CSti_TYPE_CODE:
            type_dict = map_mef3_CSti_type(rh);
            break;
        case MEFREC_ESti_TYPE_CODE:
            type_dict = map_mef3_ESti_type(rh);
            break;
        case MEFREC_SyLg_TYPE_CODE:
            type_dict = map_mef3_SyLg_type(rh);
            break;
        case MEFREC_UnRc_TYPE_CODE:
            //PyErr_SetString(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Format(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Occurred();
            type_dict = PyDict_New();
            break;
        default:
            //PyErr_SetString(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Format(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Occurred();
            type_dict = PyDict_New();
            break;
    }

    PyDict_SetItemString(rh_dict, "record_body", type_dict);

    return rh_dict;
}

PyObject *map_mef3_ri(RECORD_INDEX *ri)
{
    // Dictionaries
    PyObject *ri_dict;
    
    // Helper variables
    si1   *time_str, temp_str[256];
    si8   long_file_time;
 
    /* RECORD HEADER */
    
    // Create output dictionary   
    ri_dict = PyDict_New();
    
    // Set not not entered string   
    sprintf(temp_str, "not entered");
    
    // Insert entries into dictionary                    
    if (ri->type_string[0])
        PyDict_SetItemString(ri_dict, "type_string",
            Py_BuildValue("s", ri->type_string));
    else
        PyDict_SetItemString(ri_dict, "type_string",
            Py_BuildValue("s", temp_str));
                             
    if (ri->version_major != RECORD_INDEX_VERSION_MAJOR_NO_ENTRY)
        PyDict_SetItemString(ri_dict, "version_major",
            Py_BuildValue("B", ri->version_major));
    else
        PyDict_SetItemString(ri_dict, "version_major",
            Py_BuildValue("s", temp_str));   

    if (ri->version_minor != RECORD_INDEX_VERSION_MINOR_NO_ENTRY)
        PyDict_SetItemString(ri_dict, "version_minor",
            Py_BuildValue("B", ri->version_minor));
    else
        PyDict_SetItemString(ri_dict, "version_minor",
            Py_BuildValue("s", temp_str));

    if (ri->encryption != NULL)
        PyDict_SetItemString(ri_dict, "encryption",
            Py_BuildValue("b", ri->encryption));
    else
        PyDict_SetItemString(ri_dict, "encryption",
            Py_BuildValue("s", temp_str));

    if (ri->file_offset != RECORD_INDEX_FILE_OFFSET_NO_ENTRY)
        PyDict_SetItemString(ri_dict, "file_offset",
            Py_BuildValue("l", ri->file_offset));
    else
        PyDict_SetItemString(ri_dict, "file_offset",
            Py_BuildValue("s", temp_str));

    if (ri->time != RECORD_INDEX_TIME_NO_ENTRY){
        long_file_time = Py_BuildValue("l", ri->time);
        PyDict_SetItemString(ri_dict, "time",
            ABS(long_file_time));
    }
    else
        PyDict_SetItemString(ri_dict, "time",
            Py_BuildValue("s", temp_str));       

    return ri_dict;
}

PyObject *map_mef3_Note_type(RECORD_HEADER *rh)
{
    // Dictionary
    PyObject *dict_out;

    si1 *Note;
     
    // Dict definition
    dict_out = PyDict_New();

    // Version 1.0
    if (rh->version_major == 1 && rh->version_minor == 0) {
        Note = (si1 *) rh + MEFREC_Note_1_0_TEXT_OFFSET;
        PyDict_SetItemString(dict_out,"note",Py_BuildValue("s", Note));
    }
    // Unrecognized record version
    else {
        PyDict_SetItemString(dict_out,"note_text",Py_BuildValue("s", "Unrecognized Note version."));
    }
    
    return dict_out;
}

PyObject *map_mef3_EDFA_type(RECORD_HEADER *rh)
{
    // Dictionary
    PyObject *dict_out;

    MEFREC_EDFA_1_0 *edfa;
    si1     *annotation;

    // Dict definition
    dict_out = PyDict_New();
       
    // Version 1.0
    if (rh->version_major == 1 && rh->version_minor == 0) {
        edfa = (MEFREC_EDFA_1_0 *) ((ui1 *) rh + MEFREC_EDFA_1_0_OFFSET);
        annotation = (si1 *) rh + MEFREC_EDFA_1_0_ANNOTATION_OFFSET;
        PyDict_SetItemString(dict_out,"annotation",Py_BuildValue("s", annotation));
        PyDict_SetItemString(dict_out,"duration",Py_BuildValue("l", edfa->duration));
    }

    // Unrecognized record version
    else {
        PyDict_SetItemString(dict_out,"note_text",Py_BuildValue("s", "Unrecognized Note version."));
    }
    
    return dict_out;
}

PyObject *map_mef3_LNTP_type(RECORD_HEADER *rh)
{
    // Dictionary
    PyObject *dict_out;
    PyObject *data_list;

    MEFREC_LNTP_1_0 *lntp;
    si4     *template;
    si8     i;
    
    // Dict definition
    dict_out = PyDict_New();
    
    // Version 1.0
    if (rh->version_major == 1 && rh->version_minor == 0) {
        lntp = (MEFREC_LNTP_1_0 *) ((ui1 *) rh + MEFREC_LNTP_1_0_OFFSET);
        template = (si4 *) rh + MEFREC_LNTP_1_0_TEMPLATE_OFFSET;
        printf("Line Noise Template:\n");

        data_list= PyList_New(lntp->length);
        for (i = 0; i < lntp->length; ++i)
            PyList_SET_ITEM(data_list,i,Py_BuildValue("i",template[i]));

        PyDict_SetItemString(dict_out,"line_noise_teamplate",data_list);
    }
    // Unrecognized record version
    else {
        PyDict_SetItemString(dict_out,"note_text",Py_BuildValue("s", "Unrecognized LNTP version."));
    }
    
    return dict_out;
}

PyObject *map_mef3_Seiz_type(RECORD_HEADER *rh)
{
    // Dictionary
    PyObject *dict_out;
    PyObject *dict_channel;
    PyObject *channel_list;

    si4         i, mn1 = MEF_FALSE, mn2 = MEF_FALSE;
    MEFREC_Seiz_1_0     *seizure;
    MEFREC_Seiz_1_0_CHANNEL *channels;
    si1         time_str[32];
        
    // Dict definition
    dict_out = PyDict_New();        


    // Version 1.0
    if (rh->version_major == 1 && rh->version_minor == 0) {
        seizure = (MEFREC_Seiz_1_0 *) ((ui1 *) rh + MEFREC_Seiz_1_0_OFFSET);

        PyDict_SetItemString(dict_out,"earliest_onset_uUTC",Py_BuildValue("l", seizure->earliest_onset));
        local_date_time_string(seizure->earliest_onset, time_str);
        PyDict_SetItemString(dict_out,"earliest_onset_str",Py_BuildValue("s", time_str));

        PyDict_SetItemString(dict_out,"latest_offset_uUTC",Py_BuildValue("l", seizure->latest_offset));
        local_date_time_string(seizure->latest_offset, time_str);
        PyDict_SetItemString(dict_out,"latest_offset_str",Py_BuildValue("s", time_str));

        PyDict_SetItemString(dict_out,"duration",Py_BuildValue("l", seizure->duration));

        PyDict_SetItemString(dict_out,"number_of_channels",Py_BuildValue("i", seizure->number_of_channels));

        PyDict_SetItemString(dict_out,"onset_code",Py_BuildValue("i", seizure->onset_code));

        switch (seizure->onset_code) {
            case MEFREC_Seiz_1_0_ONSET_NO_ENTRY:
                PyDict_SetItemString(dict_out,"seizure_type",Py_BuildValue("s", "no_entry"));
                break;
            case MEFREC_Seiz_1_0_ONSET_UNKNOWN:
                PyDict_SetItemString(dict_out,"seizure_type",Py_BuildValue("s", "unknown"));
                break;
            case MEFREC_Seiz_1_0_ONSET_FOCAL:
                PyDict_SetItemString(dict_out,"seizure_type",Py_BuildValue("s", "focal"));
                break;
            case MEFREC_Seiz_1_0_ONSET_GENERALIZED:
            PyDict_SetItemString(dict_out,"seizure_type",Py_BuildValue("s", "generalized"));
                break;
            case MEFREC_Seiz_1_0_ONSET_PROPAGATED:
            PyDict_SetItemString(dict_out,"seizure_type",Py_BuildValue("s", "propagated"));
                break;
            case MEFREC_Seiz_1_0_ONSET_MIXED:
            PyDict_SetItemString(dict_out,"seizure_type",Py_BuildValue("s", "mixed"));
                break;
            default:
                PyDict_SetItemString(dict_out,"seizure_type",Py_BuildValue("s", "unrecognized"));
                break;
        }
        if (strlen(seizure->marker_name_1))
                mn1 = MEF_TRUE;
        if (strlen(seizure->marker_name_2))
                mn2 = MEF_TRUE;
        if (mn1 == MEF_TRUE && mn2 == MEF_TRUE){
            PyDict_SetItemString(dict_out,"marker_name_1",Py_BuildValue("s", seizure->marker_name_1));
            PyDict_SetItemString(dict_out,"marker_name_2",Py_BuildValue("s", seizure->marker_name_2));

        }else if (mn1 == MEF_TRUE){
            PyDict_SetItemString(dict_out,"marker_name_1",Py_BuildValue("s", seizure->marker_name_1));
            PyDict_SetItemString(dict_out,"marker_name_2",Py_BuildValue("s", "no_entry"));
        }else if (mn2 == MEF_TRUE){
            PyDict_SetItemString(dict_out,"marker_name_1",Py_BuildValue("s", "no_entry"));
            PyDict_SetItemString(dict_out,"marker_name_2",Py_BuildValue("s", seizure->marker_name_2));
        }else{
            PyDict_SetItemString(dict_out,"marker_name_1",Py_BuildValue("s", "no_entry"));
            PyDict_SetItemString(dict_out,"marker_name_2",Py_BuildValue("s", "no_entry"));
        }
        if (strlen(seizure->annotation))
            PyDict_SetItemString(dict_out,"annotation",Py_BuildValue("s", seizure->annotation));
        else
            PyDict_SetItemString(dict_out,"annotation",Py_BuildValue("s", "no_entry"));

        channels = (MEFREC_Seiz_1_0_CHANNEL *) ((ui1 *) rh + MEFREC_Seiz_1_0_CHANNELS_OFFSET);            
        
        channel_list = PyList_New(seizure->number_of_channels);
        for (i = 0; i < seizure->number_of_channels; ++i) {

            dict_channel = PyDict_New();

            if (strlen(channels[i].name))
                PyDict_SetItemString(dict_channel, "name", Py_BuildValue("s", channels[i].name));
            else
                PyDict_SetItemString(dict_channel, "name", Py_BuildValue("s", "no_entry"));

            PyDict_SetItemString(dict_channel,"onset_uUTC",Py_BuildValue("l", channels[i].onset));
            local_date_time_string(channels[i].onset, time_str);
            PyDict_SetItemString(dict_channel,"onset_str",Py_BuildValue("s", time_str));

            PyDict_SetItemString(dict_channel,"offset_uUTC",Py_BuildValue("l", channels[i].offset));
            local_date_time_string(channels[i].offset, time_str);
            PyDict_SetItemString(dict_channel,"offset_str",Py_BuildValue("s", time_str));

            PyList_SET_ITEM(channel_list, i, dict_channel);
        }
        PyDict_SetItemString(dict_out,"channels", channel_list);
    }
    // Unrecognized record version
    else {
        PyDict_SetItemString(dict_out,"note_text",Py_BuildValue("s", "Unrecognized Seiz version."));
    }
     
    return dict_out;
}

PyObject *map_mef3_CSti_type(RECORD_HEADER *rh)
{
    // Dictionary
    PyObject *dict_out;

    si4         i, mn1 = MEF_FALSE, mn2 = MEF_FALSE;
    MEFREC_CSti_1_0     *cog_stim;
    si1         time_str[32];
        
    // Dict definition
    dict_out = PyDict_New();        


    // Version 1.0
    if (rh->version_major == 1 && rh->version_minor == 0) {
        cog_stim = (MEFREC_CSti_1_0 *) ((ui1 *) rh + MEFREC_CSti_1_0_OFFSET);

        if (strlen(cog_stim->task_type))
            PyDict_SetItemString(dict_out,"task_type",Py_BuildValue("s", cog_stim->task_type));
        else
            PyDict_SetItemString(dict_out,"task_type",Py_BuildValue("s", "no_entry"));

        PyDict_SetItemString(dict_out,"stimulus_duration",Py_BuildValue("l", cog_stim->stimulus_duration));

        if (strlen(cog_stim->stimulus_type))
            PyDict_SetItemString(dict_out,"stimulus_type",Py_BuildValue("s", cog_stim->stimulus_type));
        else
            PyDict_SetItemString(dict_out,"stimulus_type",Py_BuildValue("s", "no_entry"));

        if (strlen(cog_stim->patient_response))
            PyDict_SetItemString(dict_out,"patient_response",Py_BuildValue("s", cog_stim->patient_response));
        else
            PyDict_SetItemString(dict_out,"patient_response",Py_BuildValue("s", "no_entry"));
        
    }
    // Unrecognized record version
    else {
        PyDict_SetItemString(dict_out,"note_text",Py_BuildValue("s", "Unrecognized CSti version."));
    }
     
    return dict_out;
}

PyObject *map_mef3_ESti_type(RECORD_HEADER *rh)
{
    // Dictionary
    PyObject *dict_out;

    si4         i, mn1 = MEF_FALSE, mn2 = MEF_FALSE;
    MEFREC_ESti_1_0     *el_stim;
    si1         time_str[32];
        
    // Dict definition
    dict_out = PyDict_New();        


    // Version 1.0
    if (rh->version_major == 1 && rh->version_minor == 0) {
        el_stim = (MEFREC_ESti_1_0 *) ((ui1 *) rh + MEFREC_ESti_1_0_OFFSET);

        PyDict_SetItemString(dict_out,"amplitude",Py_BuildValue("d", el_stim->amplitude));
        PyDict_SetItemString(dict_out,"frequency",Py_BuildValue("d", el_stim->frequency));
        PyDict_SetItemString(dict_out,"pulse_width",Py_BuildValue("l", el_stim->pulse_width));

        PyDict_SetItemString(dict_out,"ampunit_code",Py_BuildValue("i", el_stim->ampunit_code));
        switch (el_stim->ampunit_code) {
            case MEFREC_ESti_1_0_AMPUNIT_NO_ENTRY:
                PyDict_SetItemString(dict_out,"amplitude_unit",Py_BuildValue("s", "no_entry"));
                break;
            case MEFREC_ESti_1_0_AMPUNIT_UNKNOWN:
                PyDict_SetItemString(dict_out,"amplitude_unit",Py_BuildValue("s", "unknown"));
                break;
            case MEFREC_ESti_1_0_AMPUNIT_MA:
                PyDict_SetItemString(dict_out,"amplitude_unit",Py_BuildValue("s", "mA"));
                break;
            case MEFREC_ESti_1_0_AMPUNIT_V:
            PyDict_SetItemString(dict_out,"amplitude_unit",Py_BuildValue("s", "V"));
                break;
            default:
                PyDict_SetItemString(dict_out,"amplitude_unit",Py_BuildValue("s", "unrecognized"));
                break;
        }

        PyDict_SetItemString(dict_out,"mode_code",Py_BuildValue("i", el_stim->mode_code));
        switch (el_stim->mode_code) {
            case MEFREC_ESti_1_0_MODE_NO_ENTRY:
                PyDict_SetItemString(dict_out,"mode",Py_BuildValue("s", "no_entry"));
                break;
            case MEFREC_ESti_1_0_MODE_UNKNOWN:
                PyDict_SetItemString(dict_out,"mode",Py_BuildValue("s", "unknown"));
                break;
            case MEFREC_ESti_1_0_MODE_CURRENT:
                PyDict_SetItemString(dict_out,"mode",Py_BuildValue("s", "current"));
                break;
            case MEFREC_ESti_1_0_MODE_VOLTAGE:
            PyDict_SetItemString(dict_out,"mode",Py_BuildValue("s", "voltage"));
                break;
            default:
                PyDict_SetItemString(dict_out,"mode",Py_BuildValue("s", "unrecognized"));
                break;
        }

        if (strlen(el_stim->waveform))
            PyDict_SetItemString(dict_out,"waveform",Py_BuildValue("s", el_stim->waveform));
        else
            PyDict_SetItemString(dict_out,"waveform",Py_BuildValue("s", "no_entry"));

        if (strlen(el_stim->anode))
            PyDict_SetItemString(dict_out,"anode",Py_BuildValue("s", el_stim->anode));
        else
            PyDict_SetItemString(dict_out,"anode",Py_BuildValue("s", "no_entry"));


        if (strlen(el_stim->catode))
            PyDict_SetItemString(dict_out,"catode",Py_BuildValue("s", el_stim->catode));
        else
            PyDict_SetItemString(dict_out,"catode",Py_BuildValue("s", "no_entry"));
  
    }
    // Unrecognized record version
    else {
        PyDict_SetItemString(dict_out,"note_text",Py_BuildValue("s", "Unrecognized CSti version."));
    }
     
    return dict_out;
}

PyObject *map_mef3_SyLg_type(RECORD_HEADER *rh)
{
    // Dictionary
    PyObject *dict_out;

    si1 *log_entry;
        
    // Dict definition
    dict_out = PyDict_New();
        
    // Version 1.0
    if (rh->version_major == 1 && rh->version_minor == 0) {
        log_entry = (si1 *) rh + MEFREC_SyLg_1_0_TEXT_OFFSET;
        PyDict_SetItemString(dict_out,"text",Py_BuildValue("s", log_entry));
    }
    // Unrecognized record version
    else {
        PyDict_SetItemString(dict_out,"text",Py_BuildValue("s", "Unrecognized SyLg version."));
    }
    
    return dict_out;
}


/**************************  Other helper functions  ****************************/

si4 extract_segment_number(si1 *segment_name)
{
    si1     *c;
    si4     segment_number;
    ui4     i;

    // move pointer to the end of the string
    c = segment_name + strlen(segment_name) - 1;

    // Get to the dash
    while(*--c == '-'){
        if (*c == '/'){
            PyErr_SetString(PyExc_RuntimeError, "Segment name not in valid form XXX-000000");
            PyErr_Occurred();
            return NULL;
        }
    }
    c++;

    segment_number = (si4) strtol(c, NULL, 10);

    return segment_number;
}
 