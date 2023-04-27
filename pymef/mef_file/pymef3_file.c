

/************************************************************************************/
/****************************  MEF 3.0 Library Python Wrapper ***********************/
/************************************************************************************/


// Python wrapper for Multiscale Electrophysiology Format (MEF) version 3.0 library
// Copyright 2021, Mayo Foundation, Rochester MN. All rights reserved.
// Written by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and Dan Crepeau.

// Usage and modification of this source code is governed by the Apache 2.0 license.
// You may not use this file except in compliance with this License.
// A copy of the Apache 2.0 License may be obtained at http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under this License is distributed on an "as is" basis,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Thanks to all who acknowledge the Mayo Systems Electrophysiology Laboratory, Rochester, MN
// in academic publications of their work facilitated by this software.

// For further information about mef_lib.c API see meflib.c or documentation.

#include "pymef3_file.h"

#include <numpy/arrayobject.h>
#include <numpy/npy_math.h>

#include "meflib.c"
#include "mefrec.c"



/************************************************************************************/
/******************************  MEF write functions  *******************************/
/************************************************************************************/

static PyObject *write_mef_data_records(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_file_path;
    si1     level_1_password_arr[PASSWORD_BYTES] = {0};
    si1     level_2_password_arr[PASSWORD_BYTES] = {0};
    si1    *level_1_password;
    si1    *level_2_password;
    si1    *temp_str_bytes;
    
    si8     recording_start_uutc_time, recording_stop_uutc_time, recording_time_offset;
    
    PyObject    *py_record_list, *py_record_dict, *py_pass_1_obj, *py_pass_2_obj;
    PyObject    *temp_o, *temp_UTF_str;

    // Method specific
    FILE_PROCESSING_STRUCT *gen_fps, *rec_data_fps, *rec_idx_fps;
    si8     bytes, rb_bytes, max_rec_bytes, file_offset;
    ui4     type_code, *type_str_int, n_records, li;
    ui1     *rd;
    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], record_file_name[MEF_BASE_FILE_NAME_BYTES];
    si1     path_processed;
    
    UNIVERSAL_HEADER        *uh;
    RECORD_HEADER           *rh;
    RECORD_INDEX            *ri;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOLLLO!",
                          &py_file_path,
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &recording_time_offset,
                          &PyList_Type, &py_record_list)){
        return NULL;
    }

    // Check if the list is empty
    if (PyList_Size(py_record_list) == 0){
        Py_INCREF(Py_None);
        return Py_None;
    }

    // initialize MEF library
    (void) initialize_meflib();  
    // Apply recording offset
    MEF_globals->recording_time_offset = recording_time_offset;

    // tak care of password entries
    if (PyUnicode_Check(py_pass_1_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_1_password = NULL;
        }else{
            level_1_password = strcpy(level_1_password_arr,temp_str_bytes);
        }
    }else{
        level_1_password = NULL;
    }

    if (PyUnicode_Check(py_pass_2_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_2_password = NULL;
        }else{
            level_2_password = strcpy(level_2_password_arr,temp_str_bytes);
        }
    }else{
        level_2_password = NULL;
    }
    

    if ((level_1_password == NULL) && (level_2_password != NULL)){
        PyErr_SetString(PyExc_RuntimeError, "Level 2 password cannot be set without level 1 password.");
        PyErr_Occurred();
        return NULL;
    }

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    uh->start_time = recording_start_uutc_time;
    uh->end_time = recording_stop_uutc_time;

    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    gen_fps->password_data = process_password_data(NULL, level_1_password, level_2_password, uh);
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

    // Check for directory type
    // Segment level
    path_processed = 0;
    extract_path_parts(py_file_path, path_out, name, type);
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name[0]);
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
    rb_bytes = 0;
    for (li = 0; li<n_records; li++){
        py_record_dict = PyList_GetItem(py_record_list, li);
        // Header bytes
        bytes += RECORD_HEADER_BYTES;
        temp_o = PyDict_GetItemString(py_record_dict, "record_body");
        if (temp_o != NULL)
            rb_bytes = PyArray_ITEMSIZE((PyArrayObject *) temp_o);
        temp_o = PyDict_GetItemString(py_record_dict, "record_subbody");
        if (temp_o != NULL)
            rb_bytes += (PyArray_ITEMSIZE((PyArrayObject *) temp_o) * PyArray_SIZE((PyArrayObject *) temp_o));
        // pad if needed
        if (rb_bytes % 16 != 0)
            rb_bytes += 16 - (rb_bytes % 16);
        bytes += rb_bytes;
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
        temp_o = PyDict_GetItemString(py_record_dict, "record_header");
        map_python_rh(temp_o, rh); 

        ri->time = rh->time;
        if (rh->version_major == 0)
            rh->version_major = ri->version_major = 1;
        else
            ri->version_major = rh->version_major;
        if (rh->version_minor == 0)
            rh->version_minor = ri->version_minor = 0;
        else
            ri->version_minor = rh->version_minor;
        
        // Done with record header, do record body
        rd += RECORD_HEADER_BYTES;

        // Fork for different record types
        rh->bytes = 0;
        type_str_int = (ui4 *) rh->type_string;
        type_code = *type_str_int;
        switch (type_code) {
            case MEFREC_EDFA_TYPE_CODE:
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                map_python_EDFA_type(temp_o, (MEFREC_EDFA_1_0 *) rd);

                // Type strings
                rh->bytes += PyArray_ITEMSIZE((PyArrayObject *) temp_o);
                MEF_strncpy(ri->type_string, MEFREC_EDFA_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_EDFA_TYPE_STRING, TYPE_BYTES);
                rh->bytes = (ui4) MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_LNTP_TYPE_CODE:
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                map_python_LNTP_type(temp_o, (MEFREC_LNTP_1_0 *) rd);
                rh->bytes += PyArray_ITEMSIZE((PyArrayObject *) temp_o);
                MEF_strncpy(ri->type_string, MEFREC_LNTP_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_LNTP_TYPE_STRING, TYPE_BYTES);
                rh->bytes = (ui4) MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_Seiz_TYPE_CODE:
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                map_python_Siez_type(temp_o, (MEFREC_Seiz_1_0 *) rd);
                rh->bytes += MEFREC_Seiz_1_0_BYTES;
                MEF_strncpy(ri->type_string, MEFREC_Seiz_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Seiz_TYPE_STRING, TYPE_BYTES);
                // Inidividual channels
                temp_o = PyDict_GetItemString(py_record_dict, "record_subbody");
                if (temp_o != NULL){
                    map_python_Siez_ch_type(temp_o, (si1 *) rd+MEFREC_Seiz_1_0_BYTES);
                    rh->bytes += (PyArray_ITEMSIZE((PyArrayObject *) temp_o) * PyArray_SIZE((PyArrayObject *) temp_o));
                }
                // No need to pad, seizure structs are 16 safe
                break;

            case MEFREC_Note_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_Note_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Note_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict,"record_body");
                rh->bytes += MEF_strcpy((si1 *) rd, (si1 *) PyArray_DATA((PyArrayObject *) temp_o));
                rh->bytes = (ui4) MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_CSti_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_CSti_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_CSti_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                map_python_CSti_type(temp_o, (MEFREC_CSti_1_0 *) rd);
                rh->bytes += MEFREC_CSti_1_0_BYTES;
                rh->bytes = (ui4) MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_ESti_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_ESti_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_ESti_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                map_python_ESti_type(temp_o, (MEFREC_ESti_1_0 *) rd);
                rh->bytes += MEFREC_ESti_1_0_BYTES;

                //(ui4) MEF_pad(rd, rh->bytes, 16); // unnecessary but kept for consistency
                break;

            case MEFREC_SyLg_TYPE_CODE:
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                MEF_strncpy(ri->type_string, MEFREC_SyLg_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_SyLg_TYPE_STRING, TYPE_BYTES);
                rh->bytes += MEF_strcpy((si1 *) rd, (si1 *) PyArray_DATA((PyArrayObject *) temp_o));
                rh->bytes = (ui4) MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_Curs_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_Curs_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Curs_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                map_python_Curs_type(temp_o, (MEFREC_Curs_1_0 *) rd);
                rh->bytes += MEFREC_Curs_1_0_BYTES;
                rh->bytes = (ui4) MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_Epoc_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_Epoc_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Epoc_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict, "record_body");
                map_python_Epoc_type(temp_o, (MEFREC_Epoc_1_0 *) rd);
                rh->bytes += MEFREC_Epoc_1_0_BYTES;
                rh->bytes = (ui4) MEF_pad(rd, rh->bytes, 16);
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
    write_MEF_file(rec_idx_fps);
    free_file_processing_struct(rec_data_fps);
    free_file_processing_struct(rec_idx_fps);
    free_file_processing_struct(gen_fps);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_ts_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_file_path;
    PyObject    *py_pass_1_obj, *py_pass_2_obj;
    
    si8     recording_start_uutc_time, recording_stop_uutc_time;
    
    PyObject    *py_tmd2_dict, *py_md3_dict, *temp_UTF_str;

    // Method specific
    FILE_PROCESSING_STRUCT *gen_fps, *metadata_fps;
    UNIVERSAL_HEADER        *uh;

    si1     level_1_password_arr[PASSWORD_BYTES] = {0};
    si1     level_2_password_arr[PASSWORD_BYTES] = {0};
    si1    *level_1_password;
    si1    *level_2_password;
    si1    *temp_str_bytes;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOLLOO",
                          &py_file_path,
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &py_tmd2_dict,
                          &py_md3_dict)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_pass_1_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_1_password = NULL;
        }else{
            level_1_password = strcpy(level_1_password_arr,temp_str_bytes);
        }
    }else{
        level_1_password = NULL;
    }

    if (PyUnicode_Check(py_pass_2_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_2_password = NULL;
        }else{
            level_2_password = strcpy(level_2_password_arr,temp_str_bytes);
        }
    }else{
        level_2_password = NULL;
    }


    if ((level_1_password == NULL) && (level_2_password != NULL)){
        PyErr_SetString(PyExc_RuntimeError, "Level 2 password cannot be set without level 1 password.");
        PyErr_Occurred();
        return NULL;
    }

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    uh->start_time = recording_start_uutc_time;
    uh->end_time = recording_stop_uutc_time;
    
    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    gen_fps->password_data = process_password_data(NULL, level_1_password, level_2_password, uh);
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

    // Check for directory type
    extract_path_parts(py_file_path, path_out, name, type);
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name[0]);

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

    // Get metadata section 3 from python dict
    map_python_md3(py_md3_dict, metadata_fps->metadata.section_3);

    // Assign recording_time_offset
    MEF_globals->recording_time_offset = metadata_fps->metadata.section_3->recording_time_offset;

    write_MEF_file(metadata_fps);
    free_file_processing_struct(metadata_fps);
    free_file_processing_struct(gen_fps);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_v_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_file_path;
    PyObject    *py_pass_1_obj, *py_pass_2_obj;
    
    si8     recording_start_uutc_time, recording_stop_uutc_time;
    
    PyObject    *py_vmd2_dict, *py_md3_dict, *temp_UTF_str;

    // Method specific
    FILE_PROCESSING_STRUCT *gen_fps, *metadata_fps;
    UNIVERSAL_HEADER        *uh;

    si1     level_1_password_arr[PASSWORD_BYTES] = {0};
    si1     level_2_password_arr[PASSWORD_BYTES] = {0};
    si1    *level_1_password;
    si1    *level_2_password;
    si1    *temp_str_bytes;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOLLOO",
                          &py_file_path,
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &py_vmd2_dict,
                          &py_md3_dict)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_pass_1_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_1_password = NULL;
        }else{
            level_1_password = strcpy(level_1_password_arr,temp_str_bytes);
        }
    }else{
        level_1_password = NULL;
    }

    if (PyUnicode_Check(py_pass_2_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_2_password = NULL;
        }else{
            level_2_password = strcpy(level_2_password_arr,temp_str_bytes);
        }
    }else{
        level_2_password = NULL;
    }
   

    if ((level_1_password == NULL) && (level_2_password != NULL)){
        PyErr_SetString(PyExc_RuntimeError, "Level 2 password cannot be set without level 1 password.");
        PyErr_Occurred();
        return NULL;
    }

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;
    uh->start_time = recording_start_uutc_time;
    uh->end_time = recording_stop_uutc_time;
    
    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    gen_fps->password_data = process_password_data(NULL, level_1_password, level_2_password, uh);
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

    // Check for directory type
    extract_path_parts(py_file_path, path_out, name, type);
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for video
        uh->segment_number = extract_segment_number(&name[0]);

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
            //Fire an error that this is not video directory - hence makes no sense to write metadata
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

    // set up mef3 video metadata file
    metadata_fps = allocate_file_processing_struct(METADATA_FILE_BYTES, VIDEO_METADATA_FILE_TYPE_CODE, NULL, gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(metadata_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", py_file_path, segment_name, VIDEO_METADATA_FILE_TYPE_STRING);
    uh = metadata_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = 1;
    uh->maximum_entry_size = METADATA_FILE_BYTES;
    initialize_metadata(metadata_fps);
    metadata_fps->metadata.section_1->section_2_encryption = LEVEL_1_ENCRYPTION_DECRYPTED;
    metadata_fps->metadata.section_1->section_3_encryption = LEVEL_2_ENCRYPTION_DECRYPTED;

    // Get video metadata section 2 from python dict
    map_python_vmd2(py_vmd2_dict, metadata_fps->metadata.video_section_2);

    // Get metadata section 3 from python dict
    map_python_md3(py_md3_dict, metadata_fps->metadata.section_3);

    // Apply recording offset
    MEF_globals->recording_time_offset = recording_start_uutc_time;

    write_MEF_file(metadata_fps);

    free_file_processing_struct(metadata_fps);
    free_file_processing_struct(gen_fps);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_ts_data_and_indices(PyObject *self, PyObject *args)
{
    // Specified by user
    PyArrayObject    *raw_data;
    si1    *py_file_path;
    PyObject    *py_pass_1_obj, *py_pass_2_obj;
    si4    lossy_flag;
    si4    array_type;
    
    PyObject *temp_UTF_str;
    si8    samps_per_mef_block;

    // Method specific
    PASSWORD_DATA           *pwd;
    UNIVERSAL_HEADER    *ts_data_uh;
    FILE_PROCESSING_STRUCT  *gen_fps, *metadata_fps, *ts_idx_fps, *ts_data_fps;
    TIME_SERIES_METADATA_SECTION_2  *tmd2;
    TIME_SERIES_INDEX   *tsi;
    RED_PROCESSING_STRUCT   *rps;
    RED_BLOCK_HEADER    *block_header;

    si1     level_1_password_arr[PASSWORD_BYTES] = {0};
    si1     level_2_password_arr[PASSWORD_BYTES] = {0};
    si1    *level_1_password;
    si1    *level_2_password;
    si1    *temp_str_bytes;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     full_file_name[MEF_FULL_FILE_NAME_BYTES], file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si4     max_samp, min_samp;
    ui4     block_samps;
    si8     start_sample, ts_indices_file_bytes, samps_remaining, file_offset;
    si8     curr_time, time_inc;

    // Optional arguments
    lossy_flag = 0; // default - no lossy compression

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOLO|i",
                          &py_file_path, // full path including segment
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &samps_per_mef_block,
                          &raw_data,
                          &lossy_flag)){
        return NULL;
    }

    // check raw_data data type, convert if necessary
    array_type = PyArray_TYPE(raw_data);
    if (array_type != NPY_INT32){
        PyErr_SetString(PyExc_RuntimeError, "Incorrect data type. Please convert your NumPy array to Int32 data type!");
        PyErr_Occurred();
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_pass_1_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_1_password = NULL;
        }else{
            level_1_password = strcpy(level_1_password_arr,temp_str_bytes);
        }
    }else{
        level_1_password = NULL;
    }

    if (PyUnicode_Check(py_pass_2_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_2_password = NULL;
        }else{
            level_2_password = strcpy(level_2_password_arr,temp_str_bytes);
        }
    }else{
        level_2_password = NULL;
    }


    if ((level_1_password == NULL) && (level_2_password != NULL)){
        PyErr_SetString(PyExc_RuntimeError, "Level 2 password cannot be set without level 1 password.");
        PyErr_Occurred();
        return NULL;
    }

    
    // set up a generic mef3 fps and process the password data with it
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    pwd = process_password_data(NULL, level_1_password, level_2_password, gen_fps->universal_header);
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

    // extract the segment name and check the firectory-type (if indeed segment)
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    extract_path_parts(file_path, path_out, name, type);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // segment type/directory

        // copy the segment name for file name construction later
        MEF_strncpy(segment_name, name, MEF_BASE_FILE_NAME_BYTES);

        // extract the channel name and check the type (if indeed time-series)
        MEF_strncpy(path_in, path_out, MEF_FULL_FILE_NAME_BYTES);
        extract_path_parts(path_in, path_out, name, type);
        if (!strcmp(type,TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING)){
            // correct/corresponding directory-type
            
            // extract the session name
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


    // 
    // Read the existing time-series metadata file
    //
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    metadata_fps = read_MEF_file(NULL, full_file_name, level_1_password, pwd, NULL, USE_GLOBAL_BEHAVIOR);

    // 
    MEF_globals->recording_time_offset = metadata_fps->metadata.section_3->recording_time_offset;


    //
    // Point to and update the time-series section 2 of the metadata struct (from the .tmet file)
    // 
    // The fields in this section 2 struct will be updated her and later to reflect the
    // data (that we will be writing), in the end the updated metadata will be written (to the .tmet file)
    // 
    //

    // create a pointer to the existing time-series section 2 metadata (from the .tmet file)
    tmd2 = metadata_fps->metadata.time_series_section_2;

    // update fields in the time-series section 2 metadata based on the data (to be written)    
    tmd2->number_of_samples = (si8) PyArray_SHAPE(raw_data)[0];
    tmd2->recording_duration = (si8) (((sf8)tmd2->number_of_samples / (sf8) tmd2->sampling_frequency) * 1e6);
    tmd2->number_of_blocks = (si8) ceil((sf8) tmd2->number_of_samples / (sf8) samps_per_mef_block);
    tmd2->maximum_block_samples = (ui4) samps_per_mef_block;       


    // 
    // Set up a file-processing-struct and universal-header for the time-series indices (file)
    //

    // allocate a fps and univeral header for the ts-indices (file), based on the ts-metadata (copying the directives, password data, and raw data)
    ts_indices_file_bytes = (tmd2->number_of_blocks * TIME_SERIES_INDEX_BYTES) + UNIVERSAL_HEADER_BYTES;
    ts_idx_fps = allocate_file_processing_struct(ts_indices_file_bytes, TIME_SERIES_INDICES_FILE_TYPE_CODE, NULL, metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(ts_idx_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);

    // generate/update the ts-index file uuid and set some of the entries fields
    generate_UUID(ts_idx_fps->universal_header->file_UUID);
    ts_idx_fps->universal_header->number_of_entries = tmd2->number_of_blocks;
    ts_idx_fps->universal_header->maximum_entry_size = TIME_SERIES_INDEX_BYTES;


    // 
    // Set up a file-processing-struct and universal-header for the time-series data and write to a file
    //	

    // allocate a fps and univeral header for the ts-data, based on the ts-metadata (copying the directives, password data, and raw data, including start_)
    ts_data_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES + RED_MAX_COMPRESSED_BYTES(samps_per_mef_block, 1), TIME_SERIES_DATA_FILE_TYPE_CODE, NULL, metadata_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(ts_data_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);

    // point to the universal-header of the time-series data (file)
    ts_data_uh = ts_data_fps->universal_header;

    // generate/update the ts-data file uuid and set some of the entries fields
    generate_UUID(ts_data_uh->file_UUID);
    ts_data_uh->number_of_entries = tmd2->number_of_blocks;
    ts_data_uh->maximum_entry_size = samps_per_mef_block;

    // write the universal header of the ts-data file
    ts_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
    ts_data_fps->directives.close_file = MEF_FALSE;
    write_MEF_file(ts_data_fps);


    //
    //
    //

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
    time_inc = (si8) (((sf8) samps_per_mef_block / tmd2->sampling_frequency) * (sf8) 1e6);
    samps_remaining = tmd2->number_of_samples;
    block_header = rps->block_header;
    tsi = ts_idx_fps->time_series_indices;
    min_samp = RED_POSITIVE_INFINITY;
    max_samp = RED_NEGATIVE_INFINITY;
    block_samps = (ui4) samps_per_mef_block;
    file_offset = UNIVERSAL_HEADER_BYTES;

    start_sample = 0;

    // Write the data and update the metadata
    while (samps_remaining){

        // check
        if (samps_remaining < block_samps)
            block_samps = (ui4) samps_remaining;
		
        block_header->number_of_samples = block_samps;
        block_header->start_time = (si8) (curr_time + 0.5); // ASK Why 0.5 here?
        curr_time += time_inc;
        
        rps->original_data = rps->original_ptr = (si4 *) PyArray_DATA(raw_data) + (tmd2->number_of_samples - samps_remaining);

        // filter - comment out if don't want
        // filtps->data_length = block_samps;
        // RED_filter(filtps);

        samps_remaining -= (si8) block_samps;

        // compress
        (void) RED_encode(rps);
        ts_data_fps->universal_header->body_CRC = CRC_update((ui1 *) block_header, block_header->block_bytes, ts_data_fps->universal_header->body_CRC);
        e_fwrite((void *) block_header, sizeof(ui1), block_header->block_bytes, ts_data_fps->fp, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, EXIT_ON_FAIL);

        // time series indices
        tsi->file_offset = file_offset;
        file_offset += (tsi->block_bytes = block_header->block_bytes);
        tsi->start_time = block_header->start_time;
        tsi->start_sample = start_sample;
        start_sample += (tsi->number_of_samples = (si8) block_samps);
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
    if (tmd2->units_conversion_factor >= 0.0){
        tmd2->maximum_native_sample_value = (sf8) max_samp * tmd2->units_conversion_factor;
        tmd2->minimum_native_sample_value = (sf8) min_samp * tmd2->units_conversion_factor;
    }else{
        tmd2->maximum_native_sample_value = (sf8) min_samp * tmd2->units_conversion_factor;
        tmd2->minimum_native_sample_value = (sf8) max_samp * tmd2->units_conversion_factor;
    }
    tmd2->maximum_contiguous_blocks = tmd2->number_of_blocks;

    // calculate the CRC for the time-series data-file and set in the universal header
    ts_data_fps->universal_header->header_CRC = CRC_calculate(ts_data_fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);

    // re-write the universal header of the ts-data file (which now includes the CRC)
    e_fseek(ts_data_fps->fp, 0, SEEK_SET, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    e_fwrite(ts_data_uh, sizeof(ui1), UNIVERSAL_HEADER_BYTES, ts_data_fps->fp, ts_data_fps->full_file_name, __FUNCTION__, __LINE__, MEF_globals->behavior_on_fail);
    fclose(ts_data_fps->fp);

    // write/update the time-series metadata file
    write_MEF_file(metadata_fps);

    // write time-series indices (file)
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

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_v_indices(PyObject *self, PyObject *args)
{
    // Specified by user
    PyObject    *vi_array; 
    si1    *py_file_path;
    PyObject    *py_pass_1_obj, *py_pass_2_obj;
    
    PyObject *temp_UTF_str;
    si8     recording_start_uutc_time, recording_stop_uutc_time;

    // Method specific
    UNIVERSAL_HEADER    *uh;
    FILE_PROCESSING_STRUCT  *gen_fps, *v_idx_fps;

    si1     level_1_password_arr[PASSWORD_BYTES] = {0};
    si1     level_2_password_arr[PASSWORD_BYTES] = {0};
    si1    *level_1_password;
    si1    *level_2_password;
    si1    *temp_str_bytes;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si8     v_indices_file_bytes;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOLLO",
                          &py_file_path, // full path including segment
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &vi_array)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // Apply recording offset
    MEF_globals->recording_time_offset = recording_start_uutc_time;

    // NOTE: gen_fps is unecessart here if the metadata file with the universal header already exists, or is it?
    // tak care of password entries
    if (PyUnicode_Check(py_pass_1_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_1_password = NULL;
        }else{
            level_1_password = strcpy(level_1_password_arr,temp_str_bytes);
        }
    }else{
        level_1_password = NULL;
    }

    if (PyUnicode_Check(py_pass_2_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_2_password = NULL;
        }else{
            level_2_password = strcpy(level_2_password_arr,temp_str_bytes);
        }
    }else{
        level_2_password = NULL;
    }


    if ((level_1_password == NULL) && (level_2_password != NULL)){
        PyErr_SetString(PyExc_RuntimeError, "Level 2 password cannot be set without level 1 password.");
        PyErr_Occurred();
        return NULL;
    }

    // set up a generic mef3 fps for universal header and password data
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;

    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    gen_fps->password_data = process_password_data(NULL, level_1_password, level_2_password, uh);
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

    // Check for directory type
    MEF_strncpy(file_path, py_file_path, MEF_FULL_FILE_NAME_BYTES);
    extract_path_parts(file_path, path_out, name, type);
    if (!strcmp(type,SEGMENT_DIRECTORY_TYPE_STRING)){
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number((si1 *) &name);

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
    v_indices_file_bytes = (PyArray_SHAPE((PyArrayObject *) vi_array)[0] * VIDEO_INDEX_BYTES) + UNIVERSAL_HEADER_BYTES;
    v_idx_fps = allocate_file_processing_struct(v_indices_file_bytes, VIDEO_INDICES_FILE_TYPE_CODE, NULL, gen_fps, UNIVERSAL_HEADER_BYTES);
    MEF_snprintf(v_idx_fps->full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, VIDEO_INDICES_FILE_TYPE_STRING);
    uh = v_idx_fps->universal_header;
    generate_UUID(uh->file_UUID);
    uh->number_of_entries = PyArray_SHAPE((PyArrayObject *) vi_array)[0];
    uh->maximum_entry_size = TIME_SERIES_INDEX_BYTES;

    // Run through the python list and create indices
    map_python_vi(vi_array, v_idx_fps->video_indices);

    // write the file
    write_MEF_file(v_idx_fps);

    // clean up
    free_file_processing_struct(v_idx_fps);
    free_file_processing_struct(gen_fps);

    Py_INCREF(Py_None);
    return Py_None;
}

/************************************************************************************/
/*************************  MEF modify/append functions  ****************************/
/************************************************************************************/

static PyObject *append_ts_data_and_indices(PyObject *self, PyObject *args)
{
    // Specified by user
    PyArrayObject    *raw_data;
    si1    *py_file_path;
    PyObject    *py_pass_1_obj, *py_pass_2_obj;
    si4    lossy_flag, discontinuity_flag;
    
    PyObject *temp_UTF_str;
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

    si1     level_1_password_arr[PASSWORD_BYTES] = {0};
    si1     level_2_password_arr[PASSWORD_BYTES] = {0};
    si1    *level_1_password;
    si1    *level_2_password;
    si1    *temp_str_bytes;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     full_file_name[MEF_FULL_FILE_NAME_BYTES], file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si4     max_samp, min_samp;
    si8     start_sample, samps_remaining, block_samps, file_offset, ts_indices_file_bytes;
    sf8     curr_time, time_inc;

    // Optional arguments
    discontinuity_flag = 1; // default - appended samples are discontinuity
    lossy_flag = 0; // default - no lossy compression

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOLLLO|ii",
                          &py_file_path, // full path including segment
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &samps_per_mef_block,
                          &raw_data,
                          &discontinuity_flag,
                          &lossy_flag)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_pass_1_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_1_password = NULL;
        }else{
            level_1_password = strcpy(level_1_password_arr,temp_str_bytes);
        }
    }else{
        level_1_password = NULL;
    }

    if (PyUnicode_Check(py_pass_2_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

        if (!*temp_str_bytes){
            level_2_password = NULL;
        }else{
            level_2_password = strcpy(level_2_password_arr,temp_str_bytes);
        }
    }else{
        level_2_password = NULL;
    }


    if ((level_1_password == NULL) && (level_2_password != NULL)){
        PyErr_SetString(PyExc_RuntimeError, "Level 2 password cannot be set without level 1 password.");
        PyErr_Occurred();
        return NULL;
    }

    // We don't really need this
    gen_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES, NO_FILE_TYPE_CODE, NULL, NULL, 0);
    initialize_universal_header(gen_fps, MEF_TRUE, MEF_FALSE, MEF_TRUE);
    uh = gen_fps->universal_header;

    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    pwd = gen_fps->password_data = process_password_data(NULL, level_1_password, level_2_password, uh);
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

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
            free_file_processing_struct(gen_fps);
            return NULL;
        }

    }else{
        //Fire an error that this is not segment directory - hence makes no sense to write metadata
        PyErr_SetString(PyExc_RuntimeError, "Not a segment, exiting...");
        PyErr_Occurred();
        free_file_processing_struct(gen_fps);
        return NULL;
    }

    // Create directives so that we can modify the files, not just read them?????
    gen_directives = initialize_file_processing_directives(NULL);
    gen_directives->open_mode = FPS_R_PLUS_OPEN_MODE;
    gen_directives->close_file = MEF_FALSE;

    // Read in the metadata file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    metadata_fps = read_MEF_file(NULL, full_file_name, level_1_password, pwd, NULL, USE_GLOBAL_BEHAVIOR);
    
    if (metadata_fps == NULL){
        PyErr_SetString(PyExc_FileNotFoundError, "Metadata file does not exist, exiting...");
        PyErr_Occurred();
        free_file_processing_struct(gen_fps);
        return NULL;
    }

    tmd2 = metadata_fps->metadata.time_series_section_2;
    // We are appending so get only the end time
    metadata_fps->universal_header->end_time = recording_stop_uutc_time;

    MEF_globals->recording_time_offset = metadata_fps->metadata.section_3->recording_time_offset;

    tmd2->number_of_blocks +=  (si8) ceil((sf8) PyArray_SHAPE(raw_data)[0] / (sf8) samps_per_mef_block);
    if (samps_per_mef_block > tmd2->maximum_block_samples)
        tmd2->maximum_block_samples = (ui4) samps_per_mef_block;

    // Read only UH of indices and data files
    gen_directives->io_bytes = UNIVERSAL_HEADER_BYTES;

    // Read in the indices file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
    ts_indices_file_bytes = (tmd2->number_of_blocks * TIME_SERIES_INDEX_BYTES) + UNIVERSAL_HEADER_BYTES;
    ts_idx_fps = allocate_file_processing_struct(ts_indices_file_bytes, TIME_SERIES_INDICES_FILE_TYPE_CODE, gen_directives, NULL, 0);
    ts_idx_fps = read_MEF_file(ts_idx_fps, full_file_name, level_1_password, pwd, gen_directives, USE_GLOBAL_BEHAVIOR);

    if (ts_idx_fps == NULL){
        PyErr_SetString(PyExc_FileNotFoundError, "Index file does not exist, exiting...");
        PyErr_Occurred();
        free_file_processing_struct(gen_fps);
        free_file_processing_struct(metadata_fps);
        return NULL;
    }

    // Read in the time series data file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);
    ts_data_fps = allocate_file_processing_struct(UNIVERSAL_HEADER_BYTES + RED_MAX_COMPRESSED_BYTES(samps_per_mef_block, 1), TIME_SERIES_DATA_FILE_TYPE_CODE, gen_directives, NULL, 0);
    ts_data_fps = read_MEF_file(ts_data_fps, full_file_name, level_1_password, pwd, gen_directives, USE_GLOBAL_BEHAVIOR);
    
    if (ts_data_fps == NULL){
        PyErr_SetString(PyExc_FileNotFoundError, "Data file does not exist, exiting...");
        PyErr_Occurred();
        free_file_processing_struct(gen_fps);
        free_file_processing_struct(metadata_fps);
        free_file_processing_struct(ts_idx_fps);
        return NULL;
    }

    // Switch the directives back for wirting
    gen_directives->io_bytes = FPS_FULL_FILE;

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
    curr_time = (sf8) recording_start_uutc_time;
    time_inc = ((sf8) samps_per_mef_block / tmd2->sampling_frequency) * (sf8) 1e6;
    samps_remaining = (si8) PyArray_SHAPE(raw_data)[0];
    block_header = rps->block_header;
    if (discontinuity_flag != 1)
        block_header->flags = 0;
    start_sample = tmd2->number_of_samples;
    tmd2->number_of_samples = tmd2->number_of_samples + (si8) PyArray_SHAPE(raw_data)[0];
    tmd2->recording_duration = (si8) (((sf8)tmd2->number_of_samples / (sf8) tmd2->sampling_frequency) * 1e6);
    min_samp = RED_POSITIVE_INFINITY;
    max_samp = RED_NEGATIVE_INFINITY;
    block_samps = samps_per_mef_block; 
    //Move file_offset to the end of RED blocks
    file_offset = ts_data_fps->file_length;

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
        block_header->number_of_samples = (ui4) block_samps;
        block_header->start_time = (si8) (curr_time + 0.5);
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
        start_sample += (tsi->number_of_samples = (ui4) block_samps);
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
    Py_INCREF(Py_None);
    return Py_None;
}

// ASK No need for modify functions - can be taken care of at python level - just load and rewrite,
// memory load would be minute in these cases.

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
    PyObject    *py_password_obj;
    
    // Dictionaries
    PyObject *ses_metadata_dict;
    
    // Fuction specific
    SESSION *session;
    si1     session_path[MEF_FULL_FILE_NAME_BYTES];
    si1     password_arr[PASSWORD_BYTES] = {0};
    si1     *temp_str_bytes;
    si1     *password;
    PyObject    *temp_UTF_str;

    // Read indices flag
    si1 map_indices_flag = 1;
 
    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sO|b",
                          &py_session_path,
                          &py_password_obj,
                          &map_indices_flag)){
        return NULL;
    }
    
    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_password_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
        if (!*temp_str_bytes){
            password = NULL;
        }else{
            password = strcpy(password_arr,temp_str_bytes);
        }
    }else{
        password = NULL;
    }


    MEF_strncpy(session_path, py_session_path, MEF_FULL_FILE_NAME_BYTES);
    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    session = read_MEF_session(NULL, py_session_path, password, NULL, MEF_FALSE, MEF_TRUE);    
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

    // Session info
    ses_metadata_dict = map_mef3_session(session, map_indices_flag);
    
    return ses_metadata_dict;   
}

static PyObject *read_mef_channel_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_channel_dir;
    PyObject    *py_password_obj;
    
    // Dictionaries
    PyObject *ch_metadata_dict;
    
    // Fuction specific
    CHANNEL *channel;
    si1     password_arr[PASSWORD_BYTES] = {0};
    si1     *temp_str_bytes;
    si1     *password;
    PyObject    *temp_UTF_str;
 
    // Read indices flag
    si1 map_indices_flag = 1;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sO|b",
                          &py_channel_dir,
                          &py_password_obj,
                          &map_indices_flag)){
        return NULL;
    }
    
    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_password_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
        if (!*temp_str_bytes){
            password = NULL;
        }else{
            password = strcpy(password_arr,temp_str_bytes);
        }
    }else{
        password = NULL;
    }
    

    
    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    channel = read_MEF_channel(NULL, py_channel_dir, UNKNOWN_CHANNEL_TYPE, password, NULL, MEF_FALSE, MEF_TRUE);    

    // map the channel info
    ch_metadata_dict = map_mef3_channel(channel, map_indices_flag);
    
    return ch_metadata_dict;
}

static PyObject *read_mef_segment_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_segment_dir;
    PyObject    *py_password_obj;
    
    // Dictionaries
    PyObject *seg_metadata_dict;
    
    // Fuction specific
    SEGMENT *segment;
    si1     password_arr[PASSWORD_BYTES] = {0};
    si1     *temp_str_bytes;
    si1     *password;
    PyObject    *temp_UTF_str;

    // Read indices flag
    si1 map_indices_flag = 1;
  
    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sO|b",
                          &py_segment_dir,
                          &py_password_obj,
                          &map_indices_flag)){
        return NULL;
    }
    
    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_password_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
        if (!*temp_str_bytes){
            password = NULL;
        }else{
            password = strcpy(password_arr,temp_str_bytes);
        }
    }else{
        password = NULL;
    }
    
    segment = read_MEF_segment(NULL, py_segment_dir, UNKNOWN_CHANNEL_TYPE, password, NULL, MEF_FALSE, MEF_TRUE);    
    // map the segment info
    seg_metadata_dict = map_mef3_segment(segment,map_indices_flag);
    
    return seg_metadata_dict; 
}

static PyObject *read_mef_ts_data(PyObject *self, PyObject *args)
{
    // Specified by user
    PyObject    *py_channel_obj;
    PyObject    *ostart, *oend;
    si8     start_time, end_time;
    si8     start_samp, end_samp;
    si4     times_specified;
 
    // Python variables
    PyArrayObject    *py_array_out;

    // Method specific variables
    ui4     i, j;
    // si4     offset_to_start_samp;
    //ui8     data_len;
    ui4 n_segments;
    CHANNEL    *channel;
    ui4 start_segment, end_segment;
    
    si8  total_samps;//, samp_counter_base;
    ui8  total_data_bytes;
    ui8 start_idx, end_idx, num_blocks, num_block_in_segment;
    ui1 *compressed_data_buffer, *cdp;
    si8  segment_start_sample, segment_end_sample;
    si8  segment_start_time, segment_end_time;
    si8  block_start_time;
    FILE *fp;
    ui8 n_read, bytes_to_read;
    RED_PROCESSING_STRUCT   *rps;
    si4 sample_counter;
    ui4 max_samps;
    si4 *temp_data_buf;
    ui4 num_samps;

    si4 *decomp_data;
    sf8 *numpy_arr_data;
    
    si4 offset_into_output_buffer;
    si8 block_start_time_offset;
    
    si1 py_warning_message[256];
    si4 crc_block_failure, blocks_decoded;
    si1 last_block_decoded_flag;

    npy_intp dims[1];
    
    // Optional arguments
    times_specified = 0; // default behavior - read samples

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"OOO|i",
                          &py_channel_obj,
                          &ostart,
                          &oend,
                          &times_specified)){
        return NULL;
    }
        
    // set up mef 3 library
    (void) initialize_meflib();
    MEF_globals->behavior_on_fail = RETURN_ON_FAIL;
    
    // initialize Numpy
    import_array();


    channel = (CHANNEL *) PyArray_DATA((PyArrayObject *) py_channel_obj);
    MEF_globals->recording_time_offset = channel->metadata.section_3->recording_time_offset;

    if (channel->channel_type != TIME_SERIES_CHANNEL_TYPE) {
        PyErr_SetString(PyExc_RuntimeError, "Not a time series channel, exiting...");
        PyErr_Occurred();
        return NULL;
    }

    //TODO - take care of situation when user passes times but does not raise the flag times_specified

    // If None was passed as one of the arguments move to the start or end of the recording
    start_samp = 0;
    start_time = channel->earliest_start_time;
    end_samp = channel->metadata.time_series_section_2->number_of_samples;
    end_time = channel->latest_end_time;
    if (ostart != Py_None)
        start_samp = start_time = PyLong_AsLongLong(ostart);
    
    if (oend != Py_None)
        end_samp = end_time = PyLong_AsLongLong(oend);

    // check if valid data range
    if (times_specified && start_time >= end_time)
    {
        PyErr_SetString(PyExc_RuntimeError, "Start time later than end time, exiting...");
        PyErr_Occurred();
        return NULL;
    }
    if (!times_specified && start_samp >= end_samp)
    {
        PyErr_SetString(PyExc_RuntimeError, "Start sample larger than end sample, exiting...");
        PyErr_Occurred();
        return NULL;
    }    

    // Fire a warnings if start or stop or both are out of file
    if (times_specified){
        if (((start_time < channel->earliest_start_time) & (end_time < channel->earliest_start_time)) |
            ((start_time > channel->latest_end_time) & (end_time > channel->latest_end_time))){
            PyErr_WarnEx(PyExc_RuntimeWarning, "Start and stop times are out of file. Returning None", 1);
            Py_INCREF(Py_None);
            return Py_None;
        }
        if (end_time > channel->latest_end_time)
            PyErr_WarnEx(PyExc_RuntimeWarning, "Stop uutc later than latest end time. Will insert NaNs", 1);
        
        if (start_time < channel->earliest_start_time)
            PyErr_WarnEx(PyExc_RuntimeWarning, "Start uutc earlier than earliest start time. Will insert NaNs", 1);
    }else{
        if (((start_samp < 0) & (end_samp < 0)) |
            ((start_samp > channel->metadata.time_series_section_2->number_of_samples) & (end_samp > channel->metadata.time_series_section_2->number_of_samples))){
            PyErr_WarnEx(PyExc_RuntimeWarning, "Start and stop samples are out of file. Returning None", 1);
            Py_INCREF(Py_None);
            return Py_None;
        }
        if (end_samp > channel->metadata.time_series_section_2->number_of_samples){
            PyErr_WarnEx(PyExc_RuntimeWarning, "Stop sample larger than number of samples. Setting end sample to number of samples in channel", 1);
            end_samp = channel->metadata.time_series_section_2->number_of_samples;
        }
        if (start_samp < 0){
            PyErr_WarnEx(PyExc_RuntimeWarning, "Start sample smaller than 0. Setting start sample to 0", 1);
            start_samp = 0;
        }
    }
    // Determine the number of samples
    num_samps = 0;
    if (times_specified)
        num_samps = (ui4)((((end_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
    else
        num_samps = (ui4) (end_samp - start_samp);
        
    // Allocate numpy array
    dims[0] = num_samps;

    // Integers represent the "real" data but cannot use NaNs. this way can put data directly into numpy array
    // when decompressing - cannot do this with floats - have to be copied
    //py_array_out = PyArray_SimpleNew(1, dims, NPY_INT); // Integers for now but might have to convert to floats

    // Usnig doubles so we can use NaN values for discontinuities
    py_array_out = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if (py_array_out == NULL){
        PyErr_SetString(PyExc_RuntimeError, "Memory allocation error, please try shortening the requested segment.");
        PyErr_Occurred();
        return NULL;
    }
    numpy_arr_data = (sf8 *) PyArray_GETPTR1(py_array_out, 0);
    
    // Iterate through segments, looking for data that matches our criteria
    n_segments = (ui4) channel->number_of_segments;
    start_segment = end_segment = -1;
    
    if (times_specified) {
        start_samp = sample_for_uutc_c(start_time, channel);
        end_samp = sample_for_uutc_c(end_time, channel);
    }else{
        start_time = uutc_for_sample_c(start_samp, channel);
        end_time = uutc_for_sample_c(end_samp, channel);
    }
 
    // Find start and stop segments by uutc time
    // find start segment by finding first segment whose ending is past the start time.
    // then find stop segment by using the previous segment of the (first segment whose start is past the end time)
    for (i = 0; i < n_segments; ++i) {
        
        if (times_specified){
            segment_start_time = channel->segments[i].time_series_data_fps->universal_header->start_time;
            segment_end_time   = channel->segments[i].time_series_data_fps->universal_header->end_time;
            remove_recording_time_offset( &segment_start_time);
            remove_recording_time_offset( &segment_end_time);
                        
            if ((segment_end_time >= start_time) && (start_segment == -1)) {
                start_segment = i;
                end_segment = i;
            }
            if ((end_segment != -1) && (segment_start_time <= end_time))
                end_segment = i;
            
        }else{
            segment_start_sample = channel->segments[i].metadata_fps->metadata.time_series_section_2->start_sample;
            segment_end_sample   = channel->segments[i].metadata_fps->metadata.time_series_section_2->start_sample +
            channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_samples;
            
            if ((start_samp >= segment_start_sample) && (start_samp <= segment_end_sample))
                start_segment = i;
            if ((end_samp >= segment_start_sample) && (end_samp <= segment_end_sample))
                end_segment = i;
            
        }
    }
    
    // DAN_CONSULT - no need to do this? - the array will be filled with NaNs, warning will be fired
    // TBD should this be re-done to include a partial page if partial data is available?
    // if ((start_segment == -1) || (end_segment == -1)) //hit the end of the file, fill out data with zeros
    // {
    //     // TBD do something here?
        
    //     return DECOMP_MEF3_ERROR_NOT_ENOUGH_SAMPLES_IN_CHANNEL;
    // }
    
    // find start block in start segment
    start_idx = end_idx = 0;
    //samp_counter_base = channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->start_sample;
    for (j = 1; j < channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks; j++) {
        
        block_start_time = channel->segments[start_segment].time_series_indices_fps->time_series_indices[j].start_time;
        remove_recording_time_offset( &block_start_time);
        
        if (block_start_time > start_time) {
            start_idx = j - 1;
            break;
        }
        // starting point is in last block in segment
        start_idx = j;
    }
    
    // find stop block in stop segment
    //samp_counter_base = channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->start_sample;
    for (j = 1; j < channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks; j++) {
        
        block_start_time = channel->segments[end_segment].time_series_indices_fps->time_series_indices[j].start_time;
        remove_recording_time_offset( &block_start_time);
        
        if (block_start_time > end_time) {
            end_idx = j - 1;
            break;
        }
        // ending point is in last block in segment
        end_idx = j;
    }
    
    // find total_samps and total_data_bytes, so we can allocate buffers
    total_samps = 0;
    total_data_bytes = 0;
    
    // TBD do we care about this?
    //chan->max_sample_value = RED_MINIMUM_SAMPLE_VALUE;
    //chan->min_sample_value = RED_MAXIMUM_SAMPLE_VALUE;
    
    // normal case - everything is in one segment
    if (start_segment == end_segment) {
        if (end_idx < (ui8) (channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            total_samps += channel->segments[start_segment].time_series_indices_fps->time_series_indices[end_idx+1].start_sample -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
            total_data_bytes += channel->segments[start_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        }
        else {
            // case where end_idx is last block in segment
            total_samps += channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
            total_data_bytes += channel->segments[start_segment].time_series_data_fps->file_length -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
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
        num_block_in_segment = (ui8) channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        total_samps += channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
        total_data_bytes +=  channel->segments[start_segment].time_series_data_fps->file_length -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        num_blocks = num_block_in_segment - start_idx;

        if (channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset < 1024){
            PyErr_SetString(PyExc_RuntimeError, "Invalid index file offset, exiting...");
            PyErr_Occurred();
            return NULL;
        }
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            num_block_in_segment = (ui8) channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_blocks;
            total_samps += channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_samples;
            total_data_bytes += channel->segments[i].time_series_data_fps->file_length -
            channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += num_block_in_segment;

            if (channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset < 1024){
                PyErr_SetString(PyExc_RuntimeError, "Invalid index file offset, exiting...");
                PyErr_Occurred();
                return NULL;
            }
        }
        
        // then last segment
        num_block_in_segment = (ui8) channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        if (end_idx < (ui8) (channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            total_samps += channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].start_sample -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].start_sample;
            total_data_bytes += channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += end_idx + 1;
        }
        else {
            // case where end_idx is last block in segment
            total_samps += channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].start_sample;
            total_data_bytes += channel->segments[end_segment].time_series_data_fps->file_length -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += end_idx + 1;
        }

        if (channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx].file_offset < 1024){
            PyErr_SetString(PyExc_RuntimeError, "Invalid index file offset, exiting...");
            PyErr_Occurred();
            return NULL;
        }
    }
    
    // allocate buffers
    //data_len = total_samps;
    compressed_data_buffer = (ui1 *) malloc((size_t) total_data_bytes);
    cdp = compressed_data_buffer;
    decomp_data = (si4 *) malloc((size_t) (num_samps * sizeof(si4)));
    memset_int(decomp_data, RED_NAN, num_samps);
    // decomp_data = PyArray_GETPTR1(py_array_out, 0);
    // memset(decomp_data,NPY_NAN,sizeof(NPY_FLOAT)*num_samps);
    
    // read in RED data
    // normal case - everything is in one segment
    if (start_segment == end_segment) {
        if (channel->segments[start_segment].time_series_data_fps->fp == NULL){
            channel->segments[start_segment].time_series_data_fps->fp = fopen(channel->segments[start_segment].time_series_data_fps->full_file_name, "rb");
            channel->segments[start_segment].time_series_data_fps->fd = fileno(channel->segments[start_segment].time_series_data_fps->fp);
        }
        fp = channel->segments[start_segment].time_series_data_fps->fp;
        
        #ifdef _WIN32
            _fseeki64(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #else
            fseek(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #endif
        n_read = fread(cdp, sizeof(si1), (size_t) total_data_bytes, fp);
        if (n_read != total_data_bytes){
            sprintf(py_warning_message, "Read in fewer than expected bytes from data file in segment %d.", start_segment);
            PyErr_WarnEx(PyExc_RuntimeWarning, py_warning_message, 1);
        }
        if (channel->segments[start_segment].time_series_data_fps->directives.close_file == MEF_TRUE)
            fps_close(channel->segments[start_segment].time_series_data_fps);
    }
    // spans across segments
    else {
        // start with first segment
        if (channel->segments[start_segment].time_series_data_fps->fp == NULL){
            channel->segments[start_segment].time_series_data_fps->fp = fopen(channel->segments[start_segment].time_series_data_fps->full_file_name, "rb");
            channel->segments[start_segment].time_series_data_fps->fd = fileno(channel->segments[start_segment].time_series_data_fps->fp);
        }
        fp = channel->segments[start_segment].time_series_data_fps->fp;
        #ifdef _WIN32
             _fseeki64(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #else
            fseek(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #endif
        bytes_to_read = channel->segments[start_segment].time_series_data_fps->file_length -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
        if (n_read != bytes_to_read){
            sprintf(py_warning_message, "Read in fewer than expected bytes from data file in segment %d.", start_segment);
            PyErr_WarnEx(PyExc_RuntimeWarning, py_warning_message, 1);
        }
        cdp += n_read;
        if (channel->segments[start_segment].time_series_data_fps->directives.close_file == MEF_TRUE)
            fps_close(channel->segments[start_segment].time_series_data_fps);
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            if (channel->segments[i].time_series_data_fps->fp == NULL){
                channel->segments[i].time_series_data_fps->fp = fopen(channel->segments[i].time_series_data_fps->full_file_name, "rb");
                channel->segments[i].time_series_data_fps->fd = fileno(channel->segments[i].time_series_data_fps->fp);
            }
            fp = channel->segments[i].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[i].time_series_data_fps->file_length - 
            channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read){
                sprintf(py_warning_message, "Read in fewer than expected bytes from data file in segment %d.", i);
                PyErr_WarnEx(PyExc_RuntimeWarning, py_warning_message, 1);
            }
            cdp += n_read;
            if (channel->segments[i].time_series_data_fps->directives.close_file == MEF_TRUE)
                fps_close(channel->segments[i].time_series_data_fps);
        }
        
        // then last segment
        if (channel->segments[end_segment].time_series_data_fps->fp == NULL){
            channel->segments[end_segment].time_series_data_fps->fp = fopen(channel->segments[end_segment].time_series_data_fps->full_file_name, "rb");
            channel->segments[end_segment].time_series_data_fps->fd = fileno(channel->segments[end_segment].time_series_data_fps->fp);
        }
        num_block_in_segment = channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        if (end_idx < (ui8) (channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            fp = channel->segments[end_segment].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read){
                sprintf(py_warning_message, "Read in fewer than expected bytes from data file in segment %d.", end_segment);
                PyErr_WarnEx(PyExc_RuntimeWarning, py_warning_message, 1);
            }
            cdp += n_read;
        }
        else {
            // case where end_idx is last block in segment
            fp = channel->segments[end_segment].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[end_segment].time_series_data_fps->file_length -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read){
                sprintf(py_warning_message, "Read in fewer than expected bytes from data file in segment %d.", end_segment);
                PyErr_WarnEx(PyExc_RuntimeWarning, py_warning_message, 1);
            }
            cdp += n_read;
        }

        if (channel->segments[end_segment].time_series_data_fps->directives.close_file == MEF_TRUE)
            fps_close(channel->segments[end_segment].time_series_data_fps);

        // // then last segment
        // fp = channel->segments[end_segment].time_series_data_fps->fp;
        // fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
        // bytes_to_read = channel->segments[end_segment].time_series_data_fps->file_length -
        // channel->segments[end_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        // n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
        // if (n_read != bytes_to_read)
        //     printf("Error reading file");
        // cdp += bytes_to_read;
    }
        
    // set up RED processing struct
    cdp = compressed_data_buffer;
    max_samps = channel->metadata.time_series_section_2->maximum_block_samples;
    
    // create RED processing struct
    rps = (RED_PROCESSING_STRUCT *) calloc((size_t) 1, sizeof(RED_PROCESSING_STRUCT));
    rps->compression.mode = RED_DECOMPRESSION;
    //rps->directives.return_block_extrema = MEF_TRUE;
    rps->decompressed_ptr = rps->decompressed_data = decomp_data;
    rps->difference_buffer = (si1 *) e_calloc((size_t) RED_MAX_DIFFERENCE_BYTES(max_samps) + 1, sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // offset_to_start_samp = start_samp - channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
    
    sample_counter = 0;
    
    cdp = compressed_data_buffer;
    
    crc_block_failure = 0;
    blocks_decoded = 0;
    
	// decode the first block
    temp_data_buf = (int *) malloc((max_samps * 1.1) * sizeof(si4));
    rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
    rps->compressed_data = cdp;
    rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
    if (!check_block_crc((ui1*)(rps->block_header), max_samps, compressed_data_buffer, total_data_bytes))
    {
        crc_block_failure++;
        last_block_decoded_flag = 0;
        cdp += rps->block_header->block_bytes;
    }
    else
    {
        RED_decode(rps);
        cdp += rps->block_header->block_bytes;
        blocks_decoded++;
        last_block_decoded_flag = 1;
        
        if (times_specified)
        {
            // rps->block_header->start_time is already offset during RED_decode()
            
            if ((rps->block_header->start_time - start_time) >= 0)
                offset_into_output_buffer = (si4) ((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
            else
                offset_into_output_buffer = (si4) ((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) - 0.5);
        }
        else
            offset_into_output_buffer = (si4) (channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->start_sample +
                                               channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample) - start_samp;
        
        // copy requested samples from first block to output buffer
        // TBD this loop could be optimized
        for (i=0;i<rps->block_header->number_of_samples;i++)
        {
            if (offset_into_output_buffer < 0)
            {
                offset_into_output_buffer++;
                continue;
            }
            
            if ((ui4) offset_into_output_buffer >= num_samps)
                break;
            
            *(decomp_data + offset_into_output_buffer) = temp_data_buf[i];
            
            offset_into_output_buffer++;
        }
		sample_counter = offset_into_output_buffer;
    }

    
    // decode blocks in between the first and the last
    for (i=1;i<num_blocks-1;i++)
    {
        rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        // check that block fits fully within output array
        // this should be true, but it's possible a stray block exists out-of-order, or with a bad timestamp
        
        // we need to manually remove offset, since we are using the time value of the block bevore decoding the block
        // (normally the offset is removed during the decoding process)

        if ((rps->block_header->block_bytes == 0) || !check_block_crc((ui1*)(rps->block_header), max_samps, compressed_data_buffer, total_data_bytes))
        {
            crc_block_failure++;
            
            // two-in-a-row bad block CRCs - this is probably an unrecoverable situation, so just stop decoding.
            if (last_block_decoded_flag == 0)
                goto done_decoding;
            
            // set a flag, and keep trying successive blocks.
            last_block_decoded_flag = 0;
        }
        else
        {
            if (times_specified)
            {
                block_start_time_offset = rps->block_header->start_time;
                remove_recording_time_offset( &block_start_time_offset );
                
                // The next two checks see if the block contains out-of-bounds samples.
                // In that case, skip the block and move on
                if (block_start_time_offset < start_time)
                {
                    cdp += rps->block_header->block_bytes;
                    continue;
                }
                if (block_start_time_offset + ((rps->block_header->number_of_samples / channel->metadata.time_series_section_2->sampling_frequency) * 1e6) >= end_time)
                {
                    // Comment this out for now, it creates a strange boundary condition
                    // cdp += rps->block_header->block_bytes;
                    continue;
                }
                
                rps->decompressed_ptr = rps->decompressed_data = decomp_data + (int)((((block_start_time_offset - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
            }
            else
            {
                // prevent buffer overflow
                if ((sample_counter + rps->block_header->number_of_samples) > num_samps)
                    goto done_decoding;
                
                rps->decompressed_ptr = rps->decompressed_data = decomp_data + sample_counter;
            }
            
            RED_decode(rps);
            sample_counter += rps->block_header->number_of_samples;
            blocks_decoded++;
            last_block_decoded_flag = 1;
        }
		
        cdp += rps->block_header->block_bytes;
    }
    
	// decode last block to temp array
    if (num_blocks > 1)
    {
        rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
        if (!check_block_crc((ui1*)(rps->block_header), max_samps, compressed_data_buffer, total_data_bytes))
        {
            crc_block_failure++;
            
            goto done_decoding;
        }
        
        RED_decode(rps);
        blocks_decoded++;
        last_block_decoded_flag = 1;
        
        if (times_specified)
        {
            if ((rps->block_header->start_time  - start_time) >= 0)
                offset_into_output_buffer = (si4) ((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
            else
                offset_into_output_buffer = (si4) ((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) - 0.5);
        }
        else
            offset_into_output_buffer = sample_counter;
        
        // copy requested samples from last block to output buffer
        for (i=0;i<rps->block_header->number_of_samples;i++)
        {
            if (offset_into_output_buffer < 0)
            {
                offset_into_output_buffer++;
                continue;
            }
            
            if ((ui4) offset_into_output_buffer >= num_samps)
                break;
            
            *(decomp_data + offset_into_output_buffer) = temp_data_buf[i];
            
            offset_into_output_buffer++;
        }
    }
    
done_decoding:
    
    if (crc_block_failure > 0)
    {
        if (start_segment != end_segment)
            sprintf(py_warning_message, "CRC data block failure detected, %ld blocks skipped, in segments %d through %d.", num_blocks - blocks_decoded, start_segment, end_segment);
        else
            sprintf(py_warning_message, "CRC data block failure detected, %ld blocks skipped, in segment %d.", num_blocks - blocks_decoded, start_segment);
        
        PyErr_WarnEx(PyExc_RuntimeWarning, py_warning_message, 1);
    }
    
    // Numpy double type specific - no need to do this if we use numpy integer array and
    // put the data directly into it

    // copy requested samples from last block to output buffer
    for (i=0;i<num_samps;i++){
        if (*(decomp_data + i) == RED_NAN)
            *(numpy_arr_data + i) = NPY_NAN;
        else
            *(numpy_arr_data + i) = (sf8) *(decomp_data + i);
    }
    
    // we're done with the compressed data, get rid of it
    free (decomp_data);
    free (temp_data_buf);
    free (compressed_data_buffer);
    free (rps->difference_buffer);
    free (rps);

    return (PyObject *) py_array_out;
}

/************************************************************************************/
/****************************  MEF clean up functions  ******************************/
/************************************************************************************/

static PyObject *clean_mef_session_metadata(PyObject *self, PyObject *args)
{
    SESSION     *session;
    PyObject    *py_session_obj;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"O",
                          &py_session_obj)){
        return NULL;
    }
    session = (SESSION *) PyArray_DATA((PyArrayObject *) py_session_obj);
    free_session(session, MEF_TRUE);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *clean_mef_channel_metadata(PyObject *self, PyObject *args)
{
    CHANNEL     *channel;
    PyObject    *py_channel_obj;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"O",
                          &py_channel_obj)){
        return NULL;
    }
    channel = (CHANNEL *) PyArray_DATA((PyArrayObject *) py_channel_obj);
    free_channel(channel, MEF_TRUE);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *clean_mef_segment_metadata(PyObject *self, PyObject *args)
{
    SEGMENT     *segment;
    PyObject    *py_segment_obj;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"O",
                          &py_segment_obj)){
        return NULL;
    }
    segment = (SEGMENT *) PyArray_DATA((PyArrayObject *) py_segment_obj);
    free_segment(segment, MEF_TRUE);

    Py_INCREF(Py_None);
    return Py_None;
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

void    map_python_tmd2(PyObject *tmd2_arr, TIME_SERIES_METADATA_SECTION_2 *tmd2)
{
    memcpy(tmd2, PyArray_DATA((PyArrayObject *) tmd2_arr), sizeof(TIME_SERIES_METADATA_SECTION_2));
    return;
}

void    map_python_vmd2(PyObject *vmd2_arr, VIDEO_METADATA_SECTION_2 *vmd2)
{
    memcpy(vmd2, PyArray_DATA((PyArrayObject *) vmd2_arr), sizeof(VIDEO_METADATA_SECTION_2));
    return;
}

void    map_python_vi(PyObject *vi_arr, VIDEO_INDEX *vi)
{
    memcpy(vi, PyArray_DATA((PyArrayObject *) vi_arr), sizeof(VIDEO_INDEX)*PyArray_SHAPE((PyArrayObject *) vi_arr)[0]);
    return;
}

void    map_python_md3(PyObject *md3_arr, METADATA_SECTION_3 *md3)
{
    memcpy(md3, PyArray_DATA((PyArrayObject *) md3_arr), sizeof(METADATA_SECTION_3));
    return;
}

// SEGMENT, CHANNEL, SESSION are constructed from the above when reading

/**************************  Python record struct to Mef  ****************************/


void    map_python_rh(PyObject *rh_arr, RECORD_HEADER  *rh)
{
    memcpy(rh, PyArray_DATA((PyArrayObject *) rh_arr), sizeof(RECORD_HEADER));
    return;
}

void    map_python_EDFA_type(PyObject *EDFA_type_arr, MEFREC_EDFA_1_0  *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) EDFA_type_arr), PyArray_ITEMSIZE((PyArrayObject *) EDFA_type_arr));
    return;
}

void    map_python_LNTP_type(PyObject *LNTP_type_arr, MEFREC_LNTP_1_0  *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) LNTP_type_arr), PyArray_ITEMSIZE((PyArrayObject *) LNTP_type_arr));
    return;
}

void    map_python_Siez_type(PyObject *Siez_type_arr, MEFREC_Seiz_1_0  *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) Siez_type_arr), PyArray_ITEMSIZE((PyArrayObject *) Siez_type_arr));
    return;
}

// NOTE: void r_type because we are mapping all entries
void    map_python_Siez_ch_type(PyObject *Siez_ch_type_arr, si1 *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) Siez_ch_type_arr), PyArray_ITEMSIZE((PyArrayObject *) Siez_ch_type_arr) * PyArray_SIZE((PyArrayObject *) Siez_ch_type_arr));
    return;
}   

void    map_python_CSti_type(PyObject *CSti_type_arr, MEFREC_CSti_1_0  *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) CSti_type_arr), PyArray_ITEMSIZE((PyArrayObject *) CSti_type_arr));
    return;
}

void    map_python_ESti_type(PyObject *ESti_type_arr, MEFREC_ESti_1_0  *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) ESti_type_arr), PyArray_ITEMSIZE((PyArrayObject *) ESti_type_arr));
    return;
}

void    map_python_Curs_type(PyObject *Curs_type_arr, MEFREC_Curs_1_0  *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) Curs_type_arr), PyArray_ITEMSIZE((PyArrayObject *) Curs_type_arr));
    return;
}

void    map_python_Epoc_type(PyObject *Epoc_type_arr, MEFREC_Epoc_1_0  *r_type)
{
    memcpy(r_type, PyArray_DATA((PyArrayObject *) Epoc_type_arr), PyArray_ITEMSIZE((PyArrayObject *) Epoc_type_arr));
    return;
}


/*****************************  Mef struct to Python  *******************************/

PyObject *map_mef3_uh(UNIVERSAL_HEADER *uh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {UNIVERSAL_HEADER_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_uh_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) uh, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_md1(METADATA_SECTION_1 *md1)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {METADATA_SECTION_1_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_md1_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) md1, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_tmd2(TIME_SERIES_METADATA_SECTION_2 *tmd2)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {METADATA_SECTION_2_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_tmd2_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) tmd2, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_vmd2(VIDEO_METADATA_SECTION_2 *vmd2)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {METADATA_SECTION_2_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_vmd2_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) vmd2, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_md3(METADATA_SECTION_3 *md3)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {METADATA_SECTION_3_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_md3_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) md3, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

// TODO - call number of entries from the parent dtype
PyObject *map_mef3_ti(TIME_SERIES_INDEX *ti, si8 number_of_entries)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {number_of_entries};
    npy_intp strides[] = {TIME_SERIES_INDEX_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_ti_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) ti, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

// TODO - call number of entries from the parent dtype
PyObject *map_mef3_vi(VIDEO_INDEX *vi, si8 number_of_entries)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {number_of_entries};
    npy_intp strides[] = {VIDEO_INDEX_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_vi_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) vi, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *create_mef3_TOC(SEGMENT *segment)
{
    // Numpy TOC
    PyArrayObject *py_array_out;

    // Helpers
    TIME_SERIES_INDEX     *tsi;

    si8     number_of_entries;
    si8     prev_time, prev_sample, start_time, start_sample, samp_time_diff, seg_start_sample;
    ui4     n_samples;
    sf8     fs;
    si8     *numpy_arr_data;
    npy_intp dims[2];

    si4     i;
    
    // initialize Numpy
    import_array();

    number_of_entries = segment->time_series_indices_fps->universal_header->number_of_entries;
    tsi = segment->time_series_indices_fps->time_series_indices;
    fs = segment->metadata_fps->metadata.time_series_section_2->sampling_frequency;
    prev_time = tsi->start_time;
    prev_sample = tsi->start_sample;
    seg_start_sample = segment->metadata_fps->metadata.time_series_section_2->start_sample;

    // Create NumPy array and get pointer to data
    dims[0] = 4;
    dims[1] = number_of_entries;
    
    py_array_out = (PyArrayObject *) PyArray_SimpleNew(2, dims, NPY_INT64);
    numpy_arr_data = (si8 *) PyArray_GETPTR2(py_array_out, 0, 0);

    for(i = 0; i < number_of_entries; i++){


        start_time = tsi->start_time;
        start_sample = tsi->start_sample;
        start_sample += seg_start_sample;
        n_samples = tsi->number_of_samples;

        // Have we found a discontinuity?
        numpy_arr_data = (si8 *) PyArray_GETPTR2(py_array_out, 0, i);
        samp_time_diff = (si8) (((start_time - prev_time) - (1e6 * (start_sample - prev_sample)) / fs));
        if (samp_time_diff < (si8) (1e6/fs))
            samp_time_diff = 0;
        if  ((samp_time_diff != 0) | (i == 0)) // First entry is dicontinuity by definition
            *numpy_arr_data = 1;
        else
            *numpy_arr_data = 0;

        // Number of samples of the block
        numpy_arr_data = (si8 *) PyArray_GETPTR2(py_array_out, 1, i);
        *numpy_arr_data = (si8) n_samples;

        // Start sample
        numpy_arr_data = (si8 *) PyArray_GETPTR2(py_array_out, 2, i);
        *numpy_arr_data = (si8) start_sample;

        // Start time
        numpy_arr_data = (si8 *) PyArray_GETPTR2(py_array_out, 3, i);
        *numpy_arr_data = (si8) start_time;

        prev_time = start_time;
        prev_sample = start_sample;
        
        tsi++;
    }

    return (PyObject *) py_array_out;
}

PyObject *map_mef3_segment(SEGMENT *segment, si1 map_indices_flag)
{
    // Dictionaries
    PyObject *metadata_dict;
    PyObject *records_dict;
    PyObject *uh_dict;
    PyObject *uhs_dict;
    PyObject *TOC;
    
    // Helper variables
    si8   number_of_entries;
    
    METADATA_SECTION_1      *md1;
        TIME_SERIES_METADATA_SECTION_2  *tmd2;
        VIDEO_METADATA_SECTION_2    *vmd2;
    METADATA_SECTION_3      *md3;
    TIME_SERIES_INDEX       *tsi;
    VIDEO_INDEX       *vi;

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {sizeof(SEGMENT)};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    import_array();

    // Create python output dictionary
    metadata_dict = PyDict_New();
    
    /* SEGMENT SPECIFIC */
    descr = (PyArray_Descr *) create_segment_dtype();
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) segment, NPY_ARRAY_DEFAULT, Py_None);
    
    // Create specific dictionary
    PyDict_SetItemString(metadata_dict, "segment_specific_metadata", py_array_out); 
                            
    // Read segment records if present and add it to metadata
    if ((segment->record_indices_fps != NULL) & (segment->record_data_fps != NULL)){
        records_dict = map_mef3_records(segment->record_indices_fps, segment->record_data_fps);
        PyDict_SetItemString(metadata_dict, "records_info", records_dict);
    }

    // Assign pointers for reading metadata and read universal headers
    md1 = segment->metadata_fps->metadata.section_1;
    tmd2 = NULL;
    vmd2 = NULL;
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
    PyDict_SetItemString(metadata_dict, "section_1", map_mef3_md1(md1));

    // Create section 2 dictionary
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            PyDict_SetItemString(metadata_dict, "section_2", map_mef3_tmd2(tmd2));
            break;
        case VIDEO_CHANNEL_TYPE:
            PyDict_SetItemString(metadata_dict, "section_2", map_mef3_vmd2(vmd2));
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }
    
    // Create section 3 dictionary
    PyDict_SetItemString(metadata_dict, "section_3", map_mef3_md3(md3));

    // TODO - this should be a list - there is more indices in indices file!!!!

    // Set the TOC to NULL so that the logic works
    TOC = NULL;

    // Create indices dictionary
    
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            number_of_entries = segment->time_series_indices_fps->universal_header->number_of_entries;
            tsi = segment->time_series_indices_fps->time_series_indices;

            PyDict_SetItemString(metadata_dict, "indices", map_mef3_ti(tsi, number_of_entries));

            // Create TOC
            TOC = create_mef3_TOC(segment);
            break;
        case VIDEO_CHANNEL_TYPE:
            number_of_entries = segment->video_indices_fps->universal_header->number_of_entries;
            vi = segment->video_indices_fps->video_indices;

            PyDict_SetItemString(metadata_dict, "indices", map_mef3_vi(vi, number_of_entries));

            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }

    
    if (TOC != NULL)
        PyDict_SetItemString(metadata_dict, "TOC", TOC);

    // Get universal headers
    uhs_dict = PyDict_New();

    // Metadata
    uh_dict = map_mef3_uh(segment->metadata_fps->universal_header);
    PyDict_SetItemString(uhs_dict, "metadata", uh_dict);

    // Data an indices universal headers
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            PyDict_SetItemString(uhs_dict, "time_series_data", map_mef3_uh(segment->time_series_data_fps->universal_header));
            PyDict_SetItemString(uhs_dict, "time_series_indices", map_mef3_uh(segment->time_series_indices_fps->universal_header));
            break;
        case VIDEO_CHANNEL_TYPE:
            PyDict_SetItemString(uhs_dict, "time_series_indices", map_mef3_uh(segment->video_indices_fps->universal_header));
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }

    PyDict_SetItemString(metadata_dict, "universal_headers", uhs_dict);

    return metadata_dict;
}

PyObject *map_mef3_channel(CHANNEL *channel, si1 map_indices_flag) // This funtion also loops through segments
{
    // Dictionaries
    PyObject *metadata_dict;
    PyObject *records_dict;
    PyObject *segment_dict;
    PyObject *segments_dict;
    
    // Helper variables
    si4   i;
    
    // Method
    SEGMENT *segment;
    
    METADATA_SECTION_1      *md1;
        TIME_SERIES_METADATA_SECTION_2  *tmd2;
        VIDEO_METADATA_SECTION_2    *vmd2;
    METADATA_SECTION_3      *md3;

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {sizeof(CHANNEL)};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    import_array();

    // Create python output dictionary
    metadata_dict = PyDict_New();
    
    /* CHANNEL SPECIFIC */
    descr = (PyArray_Descr *) create_channel_dtype();
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) channel, NPY_ARRAY_DEFAULT, Py_None);
    
    // Create specific dictionary
    PyDict_SetItemString(metadata_dict, "channel_specific_metadata", py_array_out);
    
    // Read channel records if present and add it to metadata
    if ((channel->record_indices_fps != NULL) & (channel->record_data_fps != NULL)){
        records_dict = map_mef3_records(channel->record_indices_fps, channel->record_data_fps);
        PyDict_SetItemString(metadata_dict, "records_info", records_dict);
    }

    // Assign pointers for reading metadata
    md1 = channel->metadata.section_1;
    tmd2 = NULL;
    vmd2 = NULL;
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
    PyDict_SetItemString(metadata_dict, "section_1", map_mef3_md1(md1));
    
    // Create section 2 dictionary
    switch (channel->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            PyDict_SetItemString(metadata_dict, "section_2", map_mef3_tmd2(tmd2)); 
            break;
        case VIDEO_CHANNEL_TYPE:
            PyDict_SetItemString(metadata_dict, "section_2", map_mef3_vmd2(vmd2)); 
            break;
        default:
            PyErr_SetString(PyExc_RuntimeError, "Unrecognized channel type, exiting...");
            PyErr_Occurred();
            return NULL;
    }

    // Create section 3 dictionary
    PyDict_SetItemString(metadata_dict, "section_3", map_mef3_md3(md3));
    
    // Loop over segments
    for (i = 0; i < channel->number_of_segments; ++i){
        if (i == 0){
            PyDict_SetItemString(metadata_dict, "segments", PyDict_New()); 
            segments_dict = PyDict_GetItemString(metadata_dict, "segments");
        }
        // Get the channel pointer
        segment = channel->segments + i;
        // Map the channel
        segment_dict = map_mef3_segment(segment,map_indices_flag);
        // Put into ditionary
        PyDict_SetItemString(segments_dict, segment->name, segment_dict); 
    }

    return metadata_dict;
}

PyObject *map_mef3_session(SESSION *session, si1 map_indices_flag) // This funtion also loops through channels
{
    // Dictionaries
    PyObject *metadata_dict;
    PyObject *records_dict;

    PyObject *ts_metadata;
    PyObject *ts_dict;

    PyObject *v_metadata;
    PyObject *v_dict;
    
    // Helper variables
    si4   i;
    
    // Method 
    CHANNEL *channel;

    METADATA_SECTION_1      *md1;
        TIME_SERIES_METADATA_SECTION_2  *tmd2;
        VIDEO_METADATA_SECTION_2    *vmd2;
    METADATA_SECTION_3      *md3;

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {sizeof(SESSION)};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    import_array();

    // Create python output dictionary
    metadata_dict = PyDict_New();
    
    /* SESSION SPECIFIC */
    
    descr = (PyArray_Descr *) create_session_dtype();
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) session, NPY_ARRAY_DEFAULT, Py_None);
    
    // Create specific dictionary
    PyDict_SetItemString(metadata_dict, "session_specific_metadata", py_array_out);

    // Read session records if present and add it to metadata
    if ((session->record_indices_fps != NULL) & (session->record_data_fps != NULL)){
        records_dict = map_mef3_records(session->record_indices_fps, session->record_data_fps);
        PyDict_SetItemString(metadata_dict, "records_info", records_dict);
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
        PyDict_SetItemString(ts_metadata, "section_1", map_mef3_md1(md1));
        
        // Create section time series 2 dictionary
        PyDict_SetItemString(ts_metadata, "section_2", map_mef3_tmd2(tmd2)); 
        
        // Create section time series 3 dictionary
        PyDict_SetItemString(ts_metadata, "section_3", map_mef3_md3(md3));
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
        PyDict_SetItemString(v_metadata, "section_1", map_mef3_md1(md1));
        
        // Create section video 2 dictionary
        PyDict_SetItemString(v_metadata, "section_2", map_mef3_vmd2(vmd2)); 
        
        // Create section video 3 dictionary
        PyDict_SetItemString(v_metadata, "section_3", map_mef3_md3(md3));
    }

    // Loop over time series channels         
    for (i = 0; i < session->number_of_time_series_channels; ++i){
        if (i == 0){
            PyDict_SetItemString(metadata_dict, "time_series_channels", PyDict_New()); 
            ts_dict = PyDict_GetItemString(metadata_dict, "time_series_channels");
        }
        // Get the channel pointer
        channel = session->time_series_channels + i;
        // Put into ditionary
        PyDict_SetItemString(ts_dict, channel->name,  map_mef3_channel(channel, map_indices_flag)); 

    }
    // Loop over video channels
    for (i = 0; i < session->number_of_video_channels; ++i){
        if (i == 0){
            PyDict_SetItemString(metadata_dict, "video_channels", PyDict_New()); 
            v_dict = PyDict_GetItemString(metadata_dict, "video_channels");
        }
        // Get the channel pointer
        channel = session->video_channels + i;
        // Put into ditionary
        PyDict_SetItemString(v_dict, channel->name, map_mef3_channel(channel, map_indices_flag)); 
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

    ui1    *rd;
    si4     i;
    si8     number_of_records;

    RECORD_HEADER   *rh;
    //RECORD_INDEX    *ri;
 
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
    //ri = ri_fps->record_indices; // This is unneccesary now but can be used to filter read records. Will be done in python now.
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
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {RECORD_HEADER_BYTES};
    PyObject    *py_array_out;
    PyArray_Descr    *descr;

    PyObject    *rh_dict;
    MEFREC_Seiz_1_0 *sz;
    ui4     *type_str_int, type_code;

    descr = (PyArray_Descr *) create_rh_dtype();

    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) rh, 1, Py_None);
    // Create output dictionary   
    rh_dict = PyDict_New();

    PyDict_SetItemString(rh_dict, "record_header", py_array_out);
    
    type_str_int = (ui4 *) rh->type_string;
    type_code = *type_str_int;
    switch (type_code) {
        case MEFREC_Note_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_Note_type(rh));
            break;
        case MEFREC_EDFA_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_EDFA_type(rh));
            break;
        case MEFREC_LNTP_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_LNTP_type(rh));
            break;
        case MEFREC_Seiz_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_Seiz_type(rh));
            sz = (MEFREC_Seiz_1_0 *) ((ui1 *) rh + MEFREC_Seiz_1_0_OFFSET);
            if (sz->number_of_channels > 0){
                PyDict_SetItemString(rh_dict, "record_subbody", map_mef3_Seiz_ch_type(rh, sz->number_of_channels));
            }
            break;
        case MEFREC_CSti_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_CSti_type(rh));
            break;
        case MEFREC_ESti_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_ESti_type(rh));
            break;
        case MEFREC_SyLg_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_SyLg_type(rh));
            break;
        case MEFREC_Curs_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_Curs_type(rh));
            break;
        case MEFREC_Epoc_TYPE_CODE:
            PyDict_SetItemString(rh_dict, "record_body", map_mef3_Epoc_type(rh));
            break;
        case MEFREC_UnRc_TYPE_CODE:
            //PyErr_SetString(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Format(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Occurred();
            break;
        default:
            //PyErr_SetString(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Format(PyExc_RuntimeError, "\"%s\" (0x%x) is an unrecognized record type\n", rh->type_string, type_code);
            PyErr_Occurred();
            break;
    }

    return rh_dict;
}

PyObject *map_mef3_ri(RECORD_INDEX *ri)
{
    // Dictionaries
    PyObject *ri_dict;
    
    // Helper variables
    si1   temp_str[256];
    // si8   long_file_time;
 
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
        PyDict_SetItemString(ri_dict, "type_string", Py_None);
                             
    if (ri->version_major != RECORD_INDEX_VERSION_MAJOR_NO_ENTRY)
        PyDict_SetItemString(ri_dict, "version_major",
            Py_BuildValue("B", ri->version_major));
    else
        PyDict_SetItemString(ri_dict, "version_major", Py_None);   

    if (ri->version_minor != RECORD_INDEX_VERSION_MINOR_NO_ENTRY)
        PyDict_SetItemString(ri_dict, "version_minor",
            Py_BuildValue("B", ri->version_minor));
    else
        PyDict_SetItemString(ri_dict, "version_minor", Py_None);

    if (ri->encryption != ENCRYPTION_LEVEL_NO_ENTRY)
        PyDict_SetItemString(ri_dict, "encryption",
            Py_BuildValue("b", ri->encryption));
    else
        PyDict_SetItemString(ri_dict, "encryption", Py_None);

    if (ri->file_offset != RECORD_INDEX_FILE_OFFSET_NO_ENTRY)
        PyDict_SetItemString(ri_dict, "file_offset",
            Py_BuildValue("L", ri->file_offset));
    else
        PyDict_SetItemString(ri_dict, "file_offset", Py_None);

    if (ri->time != RECORD_INDEX_TIME_NO_ENTRY){
        // long_file_time = Py_BuildValue("L", ri->time);
        // PyDict_SetItemString(ri_dict, "time",
        //     ABS(long_file_time));
        PyDict_SetItemString(ri_dict, "time",
            Py_BuildValue("L", ri->time));
    }
    else
        PyDict_SetItemString(ri_dict, "time", Py_None);       

    return ri_dict;
}


PyObject *map_mef3_Note_type(RECORD_HEADER *rh)
{

    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {rh->bytes};
    PyObject    *py_array_out;
    si1     *body_p;
    ui4     str_len;
    PyArray_Descr    *descr;

    body_p = (si1 *) rh+MEFREC_Note_1_0_TEXT_OFFSET;
    str_len = (ui4) strlen(body_p);
    descr = (PyArray_Descr *) create_note_dtype_c(str_len);
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_EDFA_type(RECORD_HEADER *rh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {rh->bytes};
    PyObject    *py_array_out;
    si1     *body_p;
    ui4     str_len;
    PyArray_Descr    *descr;

    body_p = (si1 *) rh+MEFREC_EDFA_1_0_OFFSET;
    str_len = (ui4) strlen(body_p+MEFREC_EDFA_1_0_BYTES);
    descr = (PyArray_Descr *) create_edfa_dtype_c(str_len);
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_LNTP_type(RECORD_HEADER *rh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {rh->bytes};
    PyObject    *py_array_out;
    si1     *body_p;
    PyArray_Descr    *descr;
    ui4     template_size;

    template_size = rh->bytes - MEFREC_LNTP_1_0_BYTES;

    descr = (PyArray_Descr *) create_lntp_dtype_c(template_size);
    body_p = (si1 *) rh+MEFREC_LNTP_1_0_OFFSET;
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_Seiz_type(RECORD_HEADER *rh)
{
    
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {MEFREC_Seiz_1_0_BYTES};
    PyObject    *py_array_out;
    si1     *body_p;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_seiz_dtype();
    body_p = (si1 *) rh+MEFREC_Seiz_1_0_OFFSET;
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_Seiz_ch_type(RECORD_HEADER *rh, si4 number_of_channels)
{
    
    import_array();

    // Numpy array out
    npy_intp dims[] = {number_of_channels};
    npy_intp strides[] = {MEFREC_Seiz_1_0_CHANNEL_BYTES};
    PyObject    *py_array_out;
    si1     *body_p;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_seiz_ch_dtype();
    body_p = (si1 *) rh+MEFREC_Seiz_1_0_CHANNELS_OFFSET;
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_CSti_type(RECORD_HEADER *rh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {MEFREC_CSti_1_0_BYTES};
    PyObject    *py_array_out;
    si1     *body_p;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_csti_dtype();
    body_p = (si1 *) rh+MEFREC_CSti_1_0_OFFSET;
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_ESti_type(RECORD_HEADER *rh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {MEFREC_ESti_1_0_BYTES};
    PyObject    *py_array_out;
    si1     *body_p;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_esti_dtype();
    body_p = (si1 *) rh+MEFREC_ESti_1_0_OFFSET;
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_SyLg_type(RECORD_HEADER *rh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {rh->bytes};
    PyObject    *py_array_out;
    si1     *body_p;
    ui4     str_len;
    PyArray_Descr    *descr;

    body_p = (si1 *) rh+MEFREC_SyLg_1_0_TEXT_OFFSET;
    str_len = (ui4) strlen(body_p);
    descr = (PyArray_Descr *) create_sylg_dtype_c(str_len);
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_Curs_type(RECORD_HEADER *rh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {MEFREC_Curs_1_0_BYTES};
    PyObject    *py_array_out;
    si1     *body_p;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_curs_dtype();
    body_p = (si1 *) rh+MEFREC_Curs_1_0_OFFSET;
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

PyObject *map_mef3_Epoc_type(RECORD_HEADER *rh)
{
    import_array();

    // Numpy array out
    npy_intp dims[] = {1};
    npy_intp strides[] = {MEFREC_Epoc_1_0_BYTES};
    PyObject    *py_array_out;
    si1     *body_p;
    PyArray_Descr    *descr;

    descr = (PyArray_Descr *) create_epoc_dtype();
    body_p = (si1 *) rh+MEFREC_Epoc_1_0_OFFSET;
    py_array_out = PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, strides, (void *) body_p, NPY_ARRAY_DEFAULT, Py_None);
    return py_array_out;
}

/**************************  Numpy data types  ****************************/

// Records
static PyObject *create_rh_dtype()
{

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary
    op = Py_BuildValue("[(s, s),\
                         (s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s)]",

                       "record_CRC", "u4",
                       "type_string", "S", TYPE_BYTES,
                       "version_major", "u1",
                       "version_minor", "u1",
                       "encryption", "i1",
                       "bytes", "u4",
                       "time", "i8");

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);
    return (PyObject *) descr;
}

static PyObject *create_ri_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s)]",

                       "type_string", "S", TYPE_BYTES,
                       "version_major", "u1",
                       "version_minor", "u1",
                       "encryption", "i1",
                       "file_offset", "i8",
                       "time", "i8");

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_edfa_dtype(PyObject *self, PyObject *args)
{

    si8     text_len;

    if (!PyArg_ParseTuple(args,"i",
                          &text_len)){
        return NULL;
    }

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s, i)]",

                       "duration", "i8",
                       "text", "S", text_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_edfa_dtype_c(ui4 text_len)
{

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s, i)]",

                       "duration", "i8",
                       "text", "S", text_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

// // NOTE: not sure this is correct
static PyObject *create_lntp_dtype(PyObject *self, PyObject *args)
{
    si8     template_len;

    if (!PyArg_ParseTuple(args,"i",
                          &template_len)){
        return NULL;
    }

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s, i)]",

                       "length", "i8",
                       "template", "i8", template_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_lntp_dtype_c(ui4 template_len)
{

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s, i)]",

                       "length", "i8",
                       "template", "i8", template_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_note_dtype(PyObject *self, PyObject *args)
{
    si8     text_len;

    if (!PyArg_ParseTuple(args,"i",
                          &text_len)){
        return NULL;
    }

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary
    op = Py_BuildValue("[(s, s, i)]",

                       "text", "S", text_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_note_dtype_c(ui4 text_len)
{

  
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary
    op = Py_BuildValue("[(s, s, i)]",

                       "text", "S", text_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_seiz_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i)]",

                       "earliest_onset", "i8",
                       "latest_offset", "i8",
                       "duration", "i8",
                       "number_of_channels", "i4",
                       "onset_code", "i4",
                       "marker_name_1", "S", MEFREC_Seiz_1_0_MARKER_NAME_BYTES,
                       "marker_name_2", "S", MEFREC_Seiz_1_0_MARKER_NAME_BYTES,
                       "annotation", "S", MEFREC_Seiz_1_0_ANNOTATION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_seiz_ch_dtype()
{

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s, i),\
                         (s, s),\
                         (s, s)]",

                       "name", "S", MEF_BASE_FILE_NAME_BYTES,
                       "onset", "i8",
                       "offset", "i8");

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_sylg_dtype(PyObject *self, PyObject *args)
{
    si8     text_len;

    if (!PyArg_ParseTuple(args,"i",
                          &text_len)){
        return NULL;
    }

    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s, i)]",

                       "text", "S", text_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_sylg_dtype_c(ui4 text_len)
{

  
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary
    op = Py_BuildValue("[(s, s, i)]",

                       "text", "S", text_len);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_csti_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s, i),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i)]",

                        "task_type", "S", MEFREC_CSti_1_0_TASK_TYPE_BYTES,
                        "stimulus_duration", "i8",
                        "stimulus_type", "S", MEFREC_CSti_1_0_STIMULUS_TYPE_BYTES,
                        "patient_response", "S", MEFREC_CSti_1_0_PATIENT_RESPONSE_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_esti_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i)]",

                        "amplitude", "f8",
                        "frequency", "f8",
                        "pulse_width", "i8",
                        "ampunit_code", "i4",
                        "mode_code", "i4",
                        "waveform", "S", MEFREC_ESti_1_0_WAVEFORM_BYTES,
                        "anode", "S", MEFREC_ESti_1_0_ANODE_BYTES,
                        "catode", "S", MEFREC_ESti_1_0_CATODE_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_curs_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i)]",

                        "id_number", "i8",
                        "trace_timestamp", "i8",
                        "latency", "i8",
                        "value", "f8",
                        "name", "S", MEFREC_Curs_1_0_NAME_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_epoc_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i)]",

                        "id_number", "i8",
                        "timestamp", "i8",
                        "end_timestamp", "i8",
                        "duration", "i8",
                        "epoch_type", "S", MEFREC_Epoc_1_0_EPOCH_TYPE_BYTES,
                        "text", "S", MEFREC_Epoc_1_0_TEXT_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

// Library
static PyObject *create_uh_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i)]",

                       "header_CRC", "u4",
                       "body_CRC", "u4",
                       "file_type_string", "S", TYPE_BYTES,
                       "mef_version_major", "u1",
                       "mef_version_minor", "u1",
                       "byte_order_code", "u1",

                       "start_time", "i8",
                       "end_time", "i8",
                       "number_of_entries", "i8",
                       "maximum_entry_size", "i8",
                       "segment_number", "i4",

                       "channel_name", "S", MEF_BASE_FILE_NAME_BYTES,
                       "session_name", "S", MEF_BASE_FILE_NAME_BYTES,
                       "anonymized_name", "S", UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES,

                       "level_UUID", "V", UUID_BYTES,
                       "file_UUID", "V", UUID_BYTES,
                       "provenance_UUID", "V", UUID_BYTES,
                       "level_1_password_validation_field", "V", PASSWORD_VALIDATION_FIELD_BYTES,
                       "level_2_password_validation_field", "V", PASSWORD_VALIDATION_FIELD_BYTES,
                       "protected_region", "V", UNIVERSAL_HEADER_PROTECTED_REGION_BYTES,
                       "discretionary_region", "V", UNIVERSAL_HEADER_DISCRETIONARY_REGION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_md1_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i)]",

                       "section_2_encryption", "i1",
                       "section_3_encryption", "i1",
                       "protected_region", "V", METADATA_SECTION_1_PROTECTED_REGION_BYTES,
                       "discretionary_region", "V", METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_tmd2_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s, i),\
                         (s, s, i),\
                         (s, s),\
                         (s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i)]",
                       // type-independent fields
                       "channel_description", "S", METADATA_CHANNEL_DESCRIPTION_BYTES,
                       "session_description", "S", METADATA_SESSION_DESCRIPTION_BYTES,
                       "recording_duration", "i8",
                       // type-specific fields
                       "reference_description", "S", METADATA_CHANNEL_DESCRIPTION_BYTES,
                       "acquisition_channel_number", "i8",
                       "sampling_frequency", "f8",
                       "low_frequency_filter_setting", "f8",
                       "high_frequency_filter_setting", "f8",
                       "notch_filter_frequency_setting", "f8",
                       "AC_line_frequency", "f8",
                       "units_conversion_factor", "f8",
                       "units_description", "S", TIME_SERIES_METADATA_UNITS_DESCRIPTION_BYTES,
                       "maximum_native_sample_value", "f8",
                       "minimum_native_sample_value", "f8",
                       "start_sample", "i8",
                       "number_of_samples", "i8",
                       "number_of_blocks", "i8",
                       "maximum_block_bytes", "i8",
                       "maximum_block_samples", "u4",
                       "maximum_difference_bytes", "u4",
                       "block_interval", "i8",
                       "number_of_discontinuities", "i8",
                       "maximum_contiguous_blocks", "i8",
                       "maximum_contiguous_block_bytes", "i8",
                       "maximum_contiguous_samples", "i8",
                       "protected_region", "V", TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_BYTES,
                       "discretionary_region", "V", TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_vmd2_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s, i),\
                         (s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i)]",
                       // type-independent fields
                       "channel_description", "S", METADATA_CHANNEL_DESCRIPTION_BYTES,
                       "session_description", "S", METADATA_SESSION_DESCRIPTION_BYTES,
                       "recording_duration", "i8",
                       // type-specific fields
                       "horizontal_resolution", "i8",
                       "vertical_resolution", "i8",
                       "frame_rate", "f8",
                       "number_of_clips", "i8",
                       "maximum_clip_bytes", "i8",
                       "video_format", "S", VIDEO_METADATA_VIDEO_FORMAT_BYTES,
                       "video_file_CRC", "u4",
                       "protected_region", "V", VIDEO_METADATA_SECTION_2_PROTECTED_REGION_BYTES,
                       "discretionary_region", "V", VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_md3_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i)]",

                       "recording_time_offset", "i8",
                       "DST_start_time", "i8",
                       "DST_end_time", "i8",
                       "GMT_offset", "i4",
                       "subject_name_1", "S", METADATA_SUBJECT_NAME_BYTES,
                       "subject_name_2", "S", METADATA_SUBJECT_NAME_BYTES,
                       "subject_ID", "S", METADATA_SUBJECT_ID_BYTES,
                       "recording_location", "S", METADATA_RECORDING_LOCATION_BYTES,
                       "protected_region", "V", METADATA_SECTION_3_PROTECTED_REGION_BYTES,
                       "discretionary_region", "V", METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_ti_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i)]",
                       "file_offset", "i8",
                       "start_time", "i8",
                       "start_sample", "i8",
                       "number_of_samples", "u4",
                       "block_bytes", "u4",
                       "maximum_sample_value", "i4",
                       "minimum_sample_value", "i4",

                       "protected_region", "V", TIME_SERIES_INDEX_PROTECTED_REGION_BYTES,
                       "RED_block_flags", "V",
                       "RED_block_protected_region", "V", RED_BLOCK_PROTECTED_REGION_BYTES,
                       "RED_block_discretionary_region", "V", RED_BLOCK_DISCRETIONARY_REGION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_vi_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i)]",
                       "start_time", "i8",
                       "end_time", "i8",
                       "start_frame", "u4",
                       "end_frame", "u4",
                       "file_offset", "i8",
                       "clip_bytes", "i8",
                       "RED_block_protected_region", "V", VIDEO_INDEX_PROTECTED_REGION_BYTES,
                       "RED_block_discretionary_region", "V", VIDEO_INDEX_DISCRETIONARY_REGION_BYTES);

    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_segment_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i)]",
                       "channel_type", "i4",
                       "metadata_fps", "V", sizeof (void *),
                       "time_series_data_fps", "V", sizeof (void *),
                       "time_series_indices_fps", "V", sizeof (void *),
                       "video_indices_fps", "V", sizeof (void *),
                       "record_data_fps", "V", sizeof (void *),
                       "record_indices_fps", "V", sizeof (void *),
                       "name", "S", MEF_SEGMENT_BASE_FILE_NAME_BYTES,
                       "path", "S", MEF_FULL_FILE_NAME_BYTES,
                       "channel_name", "S", MEF_BASE_FILE_NAME_BYTES,
                       "level_UUID", "i1", UUID_BYTES);


    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_channel_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s)]",
                       "channel_type", "i4",
                       "metadata", "V", sizeof (METADATA),
                       "record_data_fps", "V", sizeof (void *),
                       "record_indices_fps", "V", sizeof (void *),
                       "number_of_segments", "i8",
                       "segments", "V", sizeof (void *),
                       "path", "S", MEF_FULL_FILE_NAME_BYTES,
                       "name", "S", MEF_BASE_FILE_NAME_BYTES,
                       "extension", "S", TYPE_BYTES,
                       "session_name", "S", MEF_BASE_FILE_NAME_BYTES,
                       "level_UUID", "i1", UUID_BYTES,
                       "anonymized_name", "S", UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES,
                       "maximum_number_of_records", "i8",
                       "maximum_record_bytes", "i8",
                       "earliest_start_time", "i8",
                       "latest_end_time", "i8");


    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

static PyObject *create_session_dtype()
{
    import_array();

    // Numpy array out
    PyObject    *op;
    PyArray_Descr    *descr;

    // Build dictionary

    op = Py_BuildValue("[(s, s, i),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s, i),\
                         (s, s),\
                         (s, s),\
                         (s, s),\
                         (s, s)]",
                       "time_series_metadata", "V", sizeof (METADATA),
                       "number_of_time_series_channels", "i4",
                       "time_series_channels", "V", sizeof (void *),

                       "video_metadata", "V", sizeof (METADATA),
                       "number_of_video_channels", "i4",
                       "video_channels", "V", sizeof (void *),

                       "record_data_fps", "V", sizeof (void *),
                       "record_indices_fps", "V", sizeof (void *),

                       "name", "S", MEF_BASE_FILE_NAME_BYTES,
                       "path", "S", MEF_FULL_FILE_NAME_BYTES,
                       "anonymized_name", "S", UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES,
                       "level_UUID", "i1", UUID_BYTES,

                       "maximum_number_of_records", "i8",
                       "maximum_record_bytes", "i8",
                       "earliest_start_time", "i8",
                       "latest_end_time", "i8");


    PyArray_DescrConverter(op, &descr);
    Py_DECREF(op);

    return (PyObject *) descr;
}

/**************************  Other helper functions  ****************************/

si4 check_block_crc(ui1* block_hdr_ptr, ui4 max_samps, ui1* total_data_ptr, ui8 total_data_bytes)
{
    ui8 offset_into_data, remaining_buf_size;
    si1 CRC_valid;
    RED_BLOCK_HEADER* block_header;
    
    offset_into_data = block_hdr_ptr - total_data_ptr;
    remaining_buf_size = total_data_bytes - offset_into_data;
    
    // check if remaining buffer at least contains the RED block header
    if (remaining_buf_size < RED_BLOCK_HEADER_BYTES)
        return 0;
    
    block_header = (RED_BLOCK_HEADER*) block_hdr_ptr;
    
    // check if entire block, based on size specified in header, can possibly fit in the remaining buffer
    if (block_header->block_bytes > remaining_buf_size)
        return 0;
    
    // check if size specified in header is absurdly large
    if (block_header->block_bytes > RED_MAX_COMPRESSED_BYTES(max_samps, 1))
        return 0;
	
    // check if size specified in header is too small to be valid
    if (block_header->block_bytes < RED_BLOCK_HEADER_BYTES)
        return 0;
    
    // at this point we know we have enough data to actually run the CRC calculation, so do it
    CRC_valid = CRC_validate((ui1*) block_header + CRC_BYTES, block_header->block_bytes - CRC_BYTES, block_header->block_CRC);
    
    // return output of CRC heck
    if (CRC_valid == MEF_TRUE)
        return 1;
    else
        return 0;
}

si4 extract_segment_number(si1 *segment_name)
{
    si1     *c;
    si4     segment_number;

    // move pointer to the end of the string
    c = segment_name + strlen(segment_name) - 1;

    // Get to the dash
    while(*--c == '-'){
        if (*c == '/'){
            PyErr_SetString(PyExc_RuntimeError, "Segment name not in valid form XXX-000000");
            PyErr_Occurred();
            return -1;
        }
    }
    c++;

    segment_number = (si4) strtol(c, NULL, 10);

    return segment_number;
}

si8 sample_for_uutc_c(si8 uutc, CHANNEL *channel)
{
    ui8 i, j, sample;
    sf8 native_samp_freq;
    ui8 prev_sample_number;
    si8 prev_time, seg_start_sample;
    si8 next_sample_number;
    
    native_samp_freq = channel->metadata.time_series_section_2->sampling_frequency;
    prev_sample_number = channel->segments[0].metadata_fps->metadata.time_series_section_2->start_sample;
    prev_time = channel->segments[0].time_series_indices_fps->time_series_indices[0].start_time;
    
    for (j = 0; j < (ui8) channel->number_of_segments; j++)
    {
        seg_start_sample = channel->segments[j].metadata_fps->metadata.time_series_section_2->start_sample;
        
        // initialize next_sample_number to end of current segment, in case we're on the last segment and we
        // go all the way to the end of the segment.
        // Otherwise this value will get overridden later on
        next_sample_number = seg_start_sample + channel->segments[j].metadata_fps->metadata.time_series_section_2->number_of_samples;
        
        for (i = 0; i < (ui8) channel->segments[j].metadata_fps->metadata.time_series_section_2->number_of_blocks; ++i) {
            if (channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time > uutc)
            {
                next_sample_number = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample;
                goto done;
            }
            prev_sample_number = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample;
            prev_time = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time;
        }
    }
    
done:
    
    sample = prev_sample_number + (ui8) (((((sf8) (uutc - prev_time)) / 1000000.0) * native_samp_freq) + 0.5);
    if (sample > next_sample_number)
        sample = next_sample_number;  // prevent it from going too far
    
    return(sample);
}

si8 uutc_for_sample_c(si8 sample, CHANNEL *channel)
{
    ui8 i, j, uutc;
    sf8 native_samp_freq;
    ui8 prev_sample_number; 
    si8 prev_time, seg_start_sample;


    native_samp_freq = channel->metadata.time_series_section_2->sampling_frequency;
    prev_sample_number = channel->segments[0].metadata_fps->metadata.time_series_section_2->start_sample;
    prev_time = channel->segments[0].time_series_indices_fps->time_series_indices[0].start_time;

    for (j = 0; j < (ui8) channel->number_of_segments; j++)
    {
        seg_start_sample = channel->segments[j].metadata_fps->metadata.time_series_section_2->start_sample;
        for (i = 0; i < (ui8) channel->segments[j].metadata_fps->metadata.time_series_section_2->number_of_blocks; ++i){
            if (channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample > sample)
                goto done;
            prev_sample_number = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample;
            prev_time = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time;
        }
    }

    done:
        uutc = prev_time + (ui8) ((((sf8) (sample - prev_sample_number) / native_samp_freq) * 1000000.0) + 0.5);
    
    return(uutc);
}

void memset_int(si4 *ptr, si4 value, size_t num)
{
    si4 *temp_ptr;
    int i;
    
    if (num < 1)
        return;
    
    temp_ptr = ptr;
    for (i=0;i<num;i++)
    {
        memcpy(temp_ptr, &value, sizeof(si4));
        temp_ptr++;
    }
}

static PyObject *check_mef_password(PyObject *self, PyObject *args)
{
    // We need to dive into session - get the first
    // Specified by user
    si1    *py_mef_file_path;
    PyObject    *py_password_obj;


    si1     password_arr[PASSWORD_BYTES] = {0};
    si1     *temp_str_bytes;
    si1     *password;
    PyObject    *temp_UTF_str;
        

    si1         password_bytes[PASSWORD_BYTES];
    ui1         sha[SHA256_OUTPUT_SIZE];
    si1         level_1_cumsum,level_2_cumsum;
    si1         putative_level_1_password_bytes[PASSWORD_BYTES];
    si4         i;

    size_t      nb;

    FILE *fp;
    UNIVERSAL_HEADER *uh;
 
    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sO",
                          &py_mef_file_path,
                          &py_password_obj)){
        return NULL;
    }
    
    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    if (PyUnicode_Check(py_password_obj)){
        temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
        if (!*temp_str_bytes){
            password = NULL;
        }else{
            password = strcpy(password_arr,temp_str_bytes);
        }
    }else{
        password = NULL;
    }
    

    // Allocate universal header
    uh = (UNIVERSAL_HEADER *) calloc(1, sizeof(UNIVERSAL_HEADER));
    
    // Read file universal header
    fp = fopen(py_mef_file_path,"rb");
    nb = fread((void *) uh, sizeof(UNIVERSAL_HEADER), 1, fp);
    fclose(fp);
    if (nb != 1){
        PyErr_SetString(PyExc_RuntimeError, "Error reading file, exiting...");
        PyErr_Occurred();
        free(uh);
        return NULL;
    }

    // If password is NULL check if the file is not encrypted
    if (password == NULL){
        level_1_cumsum = 0;
        for (i = 0; i < PASSWORD_VALIDATION_FIELD_BYTES; ++i){  // compare with stored level 1 hash
            level_1_cumsum += uh->level_1_password_validation_field[i];
        }
        level_2_cumsum = 0;
        for (i = 0; i < PASSWORD_VALIDATION_FIELD_BYTES; ++i){  // compare with stored level 1 hash
            level_2_cumsum += uh->level_1_password_validation_field[i];
        }

        // Clean up
        free(uh);
        
        if (level_1_cumsum | level_2_cumsum){
            return PyLong_FromLong(-1); // Wrong password
        }else{
            return PyLong_FromLong(0); // Data not encrypted
        } 
    }

    // Check the password - extracted from process_pasword_data
    extract_terminal_password_bytes(password, password_bytes);
    
    // check for level 1 access
    sha256((ui1 *) password_bytes, PASSWORD_BYTES, sha);  // generate SHA-256 hash of password bytes
    
    level_1_cumsum = 0;
    for (i = 0; i < PASSWORD_VALIDATION_FIELD_BYTES; ++i){  // compare with stored level 1 hash
        level_1_cumsum += uh->level_1_password_validation_field[i];
        if (sha[i] != uh->level_1_password_validation_field[i])
            break;
    }
    if (i == PASSWORD_BYTES) {  // Level 1 password valid - cannot be level 2 password
        // Clean up
        free(uh);
        return PyLong_FromLong(1);
    }
    
    // invalid level 1 => check if level 2 password
    for (i = 0; i < PASSWORD_BYTES; ++i)  // xor with level 2 password validation field
        putative_level_1_password_bytes[i] = sha[i] ^ uh->level_2_password_validation_field[i];
    
    sha256((ui1 *) putative_level_1_password_bytes, PASSWORD_BYTES, sha); // generate SHA-256 hash of putative level 1 password
    
    level_2_cumsum = 0;
    for (i = 0; i < PASSWORD_VALIDATION_FIELD_BYTES; ++i){  // compare with stored level 1 hash
        level_2_cumsum += uh->level_1_password_validation_field[i];
        if (sha[i] != uh->level_1_password_validation_field[i])
            break;
    }
    if (i == PASSWORD_VALIDATION_FIELD_BYTES){ // Level 2 password valid
        // Clean up
        free(uh);
        return PyLong_FromLong(2);
    }

    // Clean up
    free(uh);
    
    if (level_1_cumsum | level_2_cumsum){
        return PyLong_FromLong(-1); // Wrong password
    }else{
        return PyLong_FromLong(0); // Data not encrypted
    } 
}
