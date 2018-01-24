

/************************************************************************************/
/****************************  MEF 3.0 Library Python Wrapper ***********************/
/************************************************************************************/


// Python wrapper for Multiscale Electrophysiology Format (MEF) version 3.0 library
// Copyright 2017, Mayo Foundation, Rochester MN. All rights reserved.
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
    
    PyObject    *py_record_list, *py_record_dict, *py_seiz_chan_list, *py_pass_1_obj, *py_pass_2_obj;
    PyObject    *temp_o, *temp_UTF_str;
    Py_ssize_t  annot_bytes;

    // Method specific
    FILE_PROCESSING_STRUCT *gen_fps, *rec_data_fps, *rec_idx_fps;
    si8     bytes, rb_bytes, max_rec_bytes, file_offset;
    ui4     type_code, *type_str_int, n_records, li, lj;
    ui1     *rd;
    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], record_file_name[MEF_BASE_FILE_NAME_BYTES];
    si1     path_processed;
    
    UNIVERSAL_HEADER        *uh;
    RECORD_HEADER           *rh;
    RECORD_INDEX            *ri;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOlllO!",
                          &py_file_path,
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &recording_time_offset,
                          &PyList_Type, &py_record_list)){
        return NULL;
    }

    /// initialize MEF library
    (void) initialize_meflib();  
    // Apply recording offset
    MEF_globals->recording_time_offset = recording_time_offset;

    // tak care of password entries
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_pass_1_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyUnicode_Check(py_pass_2_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #else
        if (PyString_Check(py_pass_1_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_1_obj);

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyString_Check(py_pass_2_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_2_obj);

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #endif

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
    for (li = 0; li<n_records; li++){

        py_record_dict = PyList_GetItem(py_record_list, li);
        temp_o = PyDict_GetItemString(py_record_dict,"type_string");
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        type_str_int = (ui4 *) temp_str_bytes;
        type_code = *type_str_int;

        // Fork for different record types
        switch (type_code) {
            case MEFREC_EDFA_TYPE_CODE:
                bytes += RECORD_HEADER_BYTES; // + optional annotation
                rb_bytes = MEFREC_EDFA_1_0_BYTES;
                temp_o = PyDict_GetItemString(py_record_dict,"annotation");
                if (temp_o != NULL){
                    #if PY_MAJOR_VERSION >= 3
                        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                        annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    #else
                        temp_str_bytes = PyString_AS_STRING(temp_o);
                        annot_bytes = PyString_GET_SIZE(temp_o);
                    #endif
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
                    #if PY_MAJOR_VERSION >= 3
                        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                        annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    #else
                        temp_str_bytes = PyString_AS_STRING(temp_o);
                        annot_bytes = PyString_GET_SIZE(temp_o);
                    #endif
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
                    #if PY_MAJOR_VERSION >= 3
                        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                        annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    #else
                        temp_str_bytes = PyString_AS_STRING(temp_o);
                        annot_bytes = PyString_GET_SIZE(temp_o);
                    #endif
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
                // ASK should there by types created by this function and passed to subfunctinos?
                map_python_EDFA_type(py_record_dict, (MEFREC_EDFA_1_0 *) rd);

                // Type strings
                rh->bytes = MEFREC_EDFA_1_0_BYTES;
                MEF_strncpy(ri->type_string, MEFREC_EDFA_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_EDFA_TYPE_STRING, TYPE_BYTES);

                // Annotation
                temp_o = PyDict_GetItemString(py_record_dict,"annotation");
                if (temp_o != NULL){
                    #if PY_MAJOR_VERSION >= 3
                        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
                        annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    #else
                        temp_str_bytes = PyString_AS_STRING(temp_o);
                        annot_bytes = PyString_GET_SIZE(temp_o);
                    #endif
                    rh->bytes += MEF_strcpy((si1 *) rd + MEFREC_EDFA_1_0_BYTES, temp_str_bytes);
                }

                // Pad to 16
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_LNTP_TYPE_CODE:
                map_python_LNTP_type(py_record_dict, (MEFREC_LNTP_1_0 *) rd);
                rh->bytes = MEFREC_LNTP_1_0_BYTES;
                MEF_strncpy(ri->type_string, MEFREC_LNTP_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_LNTP_TYPE_STRING, TYPE_BYTES);
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_Seiz_TYPE_CODE:
                map_python_Siez_type(py_record_dict, (MEFREC_Seiz_1_0 *) rd);
                rh->bytes = MEFREC_Seiz_1_0_BYTES;
                MEF_strncpy(ri->type_string, MEFREC_Seiz_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Seiz_TYPE_STRING, TYPE_BYTES);
                // Inidividual channels
                py_seiz_chan_list = PyDict_GetItemString(py_record_dict,"channels");
                for (lj = 0; lj < PyList_Size(py_seiz_chan_list); lj++){
                    py_record_dict = PyList_GetItem(py_seiz_chan_list, lj);
                    map_python_Siez_type_channel(py_record_dict, (MEFREC_Seiz_1_0_CHANNEL *) (rd+MEFREC_Seiz_1_0_BYTES+(lj*MEFREC_Seiz_1_0_CHANNEL_BYTES)));
                    rh->bytes += MEFREC_Seiz_1_0_CHANNEL_BYTES;
                }
                // No need to pad, seizure structs are 16 safe
                break;

            case MEFREC_Note_TYPE_CODE:
                MEF_strncpy(ri->type_string, MEFREC_Note_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_Note_TYPE_STRING, TYPE_BYTES);
                temp_o = PyDict_GetItemString(py_record_dict,"note");
                if (temp_o != NULL){
                    #if PY_MAJOR_VERSION >= 3
                        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
                        annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    #else
                        temp_str_bytes = PyString_AS_STRING(temp_o);
                        annot_bytes = PyString_GET_SIZE(temp_o);
                    #endif
                    rh->bytes += MEF_strcpy((si1 *) rd, temp_str_bytes);
                }
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_CSti_TYPE_CODE:
                map_python_CSti_type(py_record_dict, (MEFREC_CSti_1_0 *) rd);
                MEF_strncpy(ri->type_string, MEFREC_CSti_TYPE_STRING, TYPE_BYTES);
                MEF_strncpy(rh->type_string, MEFREC_CSti_TYPE_STRING, TYPE_BYTES);
                rh->bytes = MEFREC_CSti_1_0_BYTES;
                rh->bytes = MEF_pad(rd, rh->bytes, 16);
                break;

            case MEFREC_ESti_TYPE_CODE:
                map_python_ESti_type(py_record_dict, (MEFREC_ESti_1_0 *) rd);
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
                    #if PY_MAJOR_VERSION >= 3
                        temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
                        temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
                        annot_bytes = PyBytes_GET_SIZE(temp_UTF_str);
                    #else
                        temp_str_bytes = PyString_AS_STRING(temp_o);
                        annot_bytes = PyString_GET_SIZE(temp_o);
                    #endif
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

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_ts_metadata(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_file_path;
    PyObject    *py_pass_1_obj, *py_pass_2_obj;
    
    si8     recording_start_uutc_time, recording_stop_uutc_time, recording_time_offset;
    
    PyObject    *py_tmd2_dict, *py_md3_dict, *temp_o, *temp_UTF_str;

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
    if (!PyArg_ParseTuple(args,"sOOllO!O!",
                          &py_file_path,
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &PyDict_Type, &py_tmd2_dict,
                          &PyDict_Type, &py_md3_dict)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // Apply recording offset
    temp_o = PyDict_GetItemString(py_md3_dict,"recording_time_offset");
    if (temp_o != NULL)
        recording_time_offset = PyLong_AsLong(temp_o);
    else
        recording_time_offset = METADATA_RECORDING_TIME_OFFSET_NO_ENTRY;
    MEF_globals->recording_time_offset = recording_time_offset;

    // tak care of password entries
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_pass_1_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyUnicode_Check(py_pass_2_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #else
        if (PyString_Check(py_pass_1_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_1_obj);

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyString_Check(py_pass_2_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_2_obj);

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #endif

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

    // Get time series metadata section 2 from python dict
    map_python_md3(py_md3_dict, metadata_fps->metadata.section_3);

    write_MEF_file(metadata_fps);
    free_file_processing_struct(metadata_fps); // TODO Commented for now because it is likly to close Python objects - solve later

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
    if (!PyArg_ParseTuple(args,"sOOllO!O!",
                          &py_file_path,
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &PyDict_Type, &py_vmd2_dict,
                          &PyDict_Type, &py_md3_dict)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // Apply recording offset
    MEF_globals->recording_time_offset = recording_start_uutc_time;

    // tak care of password entries
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_pass_1_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyUnicode_Check(py_pass_2_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #else
        if (PyString_Check(py_pass_1_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_1_obj);

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyString_Check(py_pass_2_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_2_obj);

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #endif

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
    
    PyObject *temp_UTF_str;
    si8    samps_per_mef_block;

    // Method specific
    PASSWORD_DATA           *pwd;
    UNIVERSAL_HEADER    *uh;
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
    si8     start_sample, ts_indices_file_bytes, samps_remaining, block_samps, file_offset;
    si8     curr_time, time_inc;

    // Optional arguments
    lossy_flag = 0; // default - no lossy compression

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOlO|i",
                          &py_file_path, // full path including segment
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &samps_per_mef_block,
                          &raw_data,
                          &lossy_flag)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // tak care of password entries
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_pass_1_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyUnicode_Check(py_pass_2_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #else
        if (PyString_Check(py_pass_1_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_1_obj);

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyString_Check(py_pass_2_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_2_obj);

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #endif

    if ((level_1_password == NULL) && (level_2_password != NULL)){
        PyErr_SetString(PyExc_RuntimeError, "Level 2 password cannot be set without level 1 password.");
        PyErr_Occurred();
        return NULL;
    }

    // NOTE: gen_fps is unecessart here if the metadata file with the universal header already exists, or is it?

    // set up a generic mef3 fps for universal header and password data
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
        // Segment - OK - extract segment number and check for time series
        uh->segment_number = extract_segment_number(&name[0]);

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

    // Get the metadata and update some fields - update the rest when processing RED
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
    metadata_fps = read_MEF_file(NULL, full_file_name, level_1_password, pwd, NULL, USE_GLOBAL_BEHAVIOR);

    MEF_globals->recording_time_offset = metadata_fps->metadata.section_3->recording_time_offset;

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
    min_samp = RED_POSITIVE_INFINITY;
    max_samp = RED_NEGATIVE_INFINITY;
    block_samps = samps_per_mef_block;
    file_offset = UNIVERSAL_HEADER_BYTES;

    start_sample = 0;

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

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *write_mef_v_indices(PyObject *self, PyObject *args)
{
    // Specified by user
    PyObject    *vi_list, *vi_dict; 
    si1    *py_file_path;
    PyObject    *py_pass_1_obj, *py_pass_2_obj;
    
    PyObject *temp_UTF_str;
    si8     recording_start_uutc_time, recording_stop_uutc_time;

    // Method specific
    UNIVERSAL_HEADER    *uh;
    FILE_PROCESSING_STRUCT  *gen_fps, *v_idx_fps;
    VIDEO_INDEX   *vi;

    si1     level_1_password_arr[PASSWORD_BYTES] = {0};
    si1     level_2_password_arr[PASSWORD_BYTES] = {0};
    si1    *level_1_password;
    si1    *level_2_password;
    si1    *temp_str_bytes;

    si1     path_in[MEF_FULL_FILE_NAME_BYTES], path_out[MEF_FULL_FILE_NAME_BYTES], name[MEF_BASE_FILE_NAME_BYTES], type[TYPE_BYTES];
    si1     file_path[MEF_FULL_FILE_NAME_BYTES], segment_name[MEF_BASE_FILE_NAME_BYTES];
    si4     li;
    si8     v_indices_file_bytes;

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOllO!",
                          &py_file_path, // full path including segment
                          &py_pass_1_obj,
                          &py_pass_2_obj,
                          &recording_start_uutc_time,
                          &recording_stop_uutc_time,
                          &PyList_Type, &vi_list)){
        return NULL;
    }

    // initialize MEF library
    (void) initialize_meflib();

    // Apply recording offset
    MEF_globals->recording_time_offset = recording_start_uutc_time;

    // NOTE: gen_fps is unecessart here if the metadata file with the universal header already exists, or is it?
    // tak care of password entries
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_pass_1_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyUnicode_Check(py_pass_2_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #else
        if (PyString_Check(py_pass_1_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_1_obj);

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyString_Check(py_pass_2_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_2_obj);

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #endif

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
    li = 0;
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
    si8     start_sample, samps_remaining, block_samps, file_offset;
    sf8     curr_time, time_inc;

    // Optional arguments
    discontinuity_flag = 1; // default - appended samples are discontinuity
    lossy_flag = 0; // default - no lossy compression

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOlllO|ii",
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
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_pass_1_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_1_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyUnicode_Check(py_pass_2_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_pass_2_obj, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #else
        if (PyString_Check(py_pass_1_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_1_obj);

            level_1_password = strcpy(level_1_password_arr, temp_str_bytes);
        }else{
            level_1_password = NULL;
        }

        if (PyString_Check(py_pass_2_obj)){
            temp_str_bytes = PyString_AS_STRING(py_pass_2_obj);

            level_2_password = strcpy(level_2_password_arr, temp_str_bytes);
        }else{
            level_2_password = NULL;
        }
    #endif

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
    metadata_fps = read_MEF_file(NULL, full_file_name, level_1_password, pwd, NULL, USE_GLOBAL_BEHAVIOR);
    tmd2 = metadata_fps->metadata.time_series_section_2;
    // We are appending so get only the end time
    metadata_fps->universal_header->end_time = recording_stop_uutc_time;

    MEF_globals->recording_time_offset = metadata_fps->metadata.section_3->recording_time_offset;

    tmd2->number_of_blocks +=  (si8) ceil((sf8) PyArray_SHAPE(raw_data)[0] / (sf8) samps_per_mef_block);
    if (samps_per_mef_block > tmd2->maximum_block_samples)
        tmd2->maximum_block_samples = samps_per_mef_block;

    // Read in the indices file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
    ts_idx_fps = read_MEF_file(NULL, full_file_name, level_1_password, pwd, gen_directives, USE_GLOBAL_BEHAVIOR);

    // Read in the time series data file
    MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s", file_path, segment_name, TIME_SERIES_DATA_FILE_TYPE_STRING);
    ts_data_fps = read_MEF_file(NULL, full_file_name, level_1_password, pwd, gen_directives, USE_GLOBAL_BEHAVIOR);
    
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
    if (discontinuity_flag != 1)
        block_header->flags = 0;
    start_sample = tmd2->number_of_samples;
    tmd2->number_of_samples = tmd2->number_of_samples + (si8) PyArray_SHAPE(raw_data)[0];
    min_samp = RED_POSITIVE_INFINITY;
    max_samp = RED_NEGATIVE_INFINITY;
    block_samps = samps_per_mef_block; 

    //Move file_offset to the end of RED blocks
    file_offset = ts_data_fps->raw_data_bytes;

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
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_password_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #else
        if (PyString_Check(py_password_obj)){
            temp_str_bytes = PyString_AS_STRING(py_password_obj);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #endif

    MEF_strncpy(session_path, py_session_path, MEF_FULL_FILE_NAME_BYTES);
    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    session = read_MEF_session(NULL, py_session_path, password, NULL, MEF_FALSE, MEF_TRUE);    
    MEF_globals->behavior_on_fail = EXIT_ON_FAIL;

    // Session info
    ses_metadata_dict = map_mef3_session(session, map_indices_flag);

    // clean up
    free_session(session, MEF_TRUE); 
    
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
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_password_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #else
        if (PyString_Check(py_password_obj)){
            temp_str_bytes = PyString_AS_STRING(py_password_obj);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #endif

    
    MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
    channel = read_MEF_channel(NULL, py_channel_dir, UNKNOWN_CHANNEL_TYPE, password, NULL, MEF_FALSE, MEF_TRUE);    

    // map the channel info
    ch_metadata_dict = map_mef3_channel(channel, map_indices_flag);

    // clean up
    free_channel(channel, MEF_TRUE);
    
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
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_password_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #else
        if (PyString_Check(py_password_obj)){
            temp_str_bytes = PyString_AS_STRING(py_password_obj);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #endif

    segment = read_MEF_segment(NULL, py_segment_dir, UNKNOWN_CHANNEL_TYPE, password, NULL, MEF_FALSE, MEF_TRUE);    

    // map the segment info
    seg_metadata_dict = map_mef3_segment(segment,map_indices_flag);

    // clean up
    free_segment(segment, MEF_TRUE);
    
    return seg_metadata_dict; 
}

static PyObject *read_mef_ts_data(PyObject *self, PyObject *args)
{
    // Specified by user
    si1    *py_channel_path;
    PyObject    *py_password_obj;
    PyObject    *ostart, *oend;
    si8     start_time, end_time;
    si8     start_samp, end_samp;
    si4     times_specified;
 
    // Python variables
    PyArrayObject    *py_array_out;

    // Method specific variables
    si4     i, j;
    // si4     offset_to_start_samp;
    //ui8     data_len;
    ui4 n_segments;
    CHANNEL    *channel;
    si4 start_segment, end_segment;
    
    si1  channel_path[MEF_FULL_FILE_NAME_BYTES];
    si8  total_samps;//, samp_counter_base;
    ui8  total_data_bytes;
    ui8 start_idx, end_idx, num_blocks;
    ui1 *compressed_data_buffer, *cdp;
    si8  segment_start_sample, segment_end_sample;
    si8  segment_start_time, segment_end_time;
    si8  block_start_time;
    si4 num_block_in_segment;
    FILE *fp;
    ui8 n_read, bytes_to_read;
    RED_PROCESSING_STRUCT   *rps;
    si4 sample_counter;
    ui4 max_samps;
    si4 *temp_data_buf;
    si4 num_samps;

    si4 *decomp_data;
    sf8 *numpy_arr_data;
    
    si4 offset_into_output_buffer;
    si8 block_start_time_offset;

    si1     password_arr[PASSWORD_BYTES] = {0};
    si1     *temp_str_bytes;
    si1     *password;
    PyObject    *temp_UTF_str;
    npy_intp dims[1];
    
    // Optional arguments
    times_specified = 0; // default behavior - read samples

    // --- Parse the input --- 
    if (!PyArg_ParseTuple(args,"sOOO|i",
                          &py_channel_path,
                          &py_password_obj,
                          &ostart,
                          &oend,
                          &times_specified)){
        return NULL;
    }
        
    // set up mef 3 library
    (void) initialize_meflib();
    MEF_globals->behavior_on_fail = RETURN_ON_FAIL;
    
    // initialize Numpy
    #if PY_MAJOR_VERSION >= 3
        import_array();
    #else
        init_numpy();
    #endif

    // tak care of password entries
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_password_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #else
        if (PyString_Check(py_password_obj)){
            temp_str_bytes = PyString_AS_STRING(py_password_obj);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #endif

    MEF_strncpy(channel_path, py_channel_path, MEF_FULL_FILE_NAME_BYTES); // might be unnecesasry
    channel = read_MEF_channel(NULL, channel_path, TIME_SERIES_CHANNEL_TYPE, password, NULL, MEF_FALSE, MEF_FALSE);
    
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
        start_samp = start_time = PyLong_AsLong(ostart);
    
    if (oend != Py_None)
        end_samp = end_time = PyLong_AsLong(oend);

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
        num_samps = (si4)(((end_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency);
    else
        num_samps = end_samp - start_samp;
        
    // Allocate numpy array
    dims[0] = num_samps;

    // Integers represent the "real" data but cannot use NaNs. this way can put data directly into numpy array
    // when decompressing - cannot do this with floats - have to be copied
    //py_array_out = PyArray_SimpleNew(1, dims, NPY_INT); // Integers for now but might have to convert to floats

    // Usnig doubles so we can use NaN values for discontinuities
    py_array_out = (PyArrayObject *) PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    numpy_arr_data = (sf8 *) PyArray_GETPTR1(py_array_out, 0);
    
    // Iterate through segments, looking for data that matches our criteria
    n_segments = channel->number_of_segments;
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
        if (end_idx < (channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            total_samps += channel->segments[start_segment].time_series_indices_fps->time_series_indices[end_idx+1].start_sample -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
            //fprintf(stderr, "total_samps = %d\n", total_samps);
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
        num_block_in_segment = channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        total_samps += channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
        total_data_bytes +=  channel->segments[start_segment].time_series_data_fps->file_length -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        num_blocks = num_block_in_segment - start_idx;
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            num_block_in_segment = channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_blocks;
            total_samps += channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_samples;
            total_data_bytes += channel->segments[i].time_series_data_fps->file_length -
            channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += num_block_in_segment;
        }
        
        // then last segment
        num_block_in_segment = channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        if (end_idx < (channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
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
        fp = channel->segments[start_segment].time_series_data_fps->fp;
        fseek(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        n_read = fread(cdp, sizeof(si1), (size_t) total_data_bytes, fp);
        if (n_read != total_data_bytes)
            printf("Error reading file\n");
    }
    // spans across segments
    else {
        // start with first segment
        fp = channel->segments[start_segment].time_series_data_fps->fp;
        fseek(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        bytes_to_read = channel->segments[start_segment].time_series_data_fps->file_length -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
        if (n_read != bytes_to_read)
            printf("Error reading file\n");
        cdp += bytes_to_read;
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            fp = channel->segments[i].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[i].time_series_data_fps->file_length - 
            channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read)
                printf("Error reading file\n");
            cdp += bytes_to_read;
        }
        
        // then last segment
        num_block_in_segment = channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        if (end_idx < (channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            fp = channel->segments[end_segment].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read)
                printf("Error reading file\n");
            cdp += bytes_to_read;
        }
        else {
            // case where end_idx is last block in segment
            fp = channel->segments[end_segment].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[end_segment].time_series_data_fps->file_length -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read)
                printf("Error reading file\n");
            cdp += bytes_to_read;
        }

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
    rps->difference_buffer = (si1 *) e_calloc((size_t) RED_MAX_DIFFERENCE_BYTES(max_samps), sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // offset_to_start_samp = start_samp - channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
    
    sample_counter = 0;
    
    cdp = compressed_data_buffer;
    
    // TBD use real max block length
    temp_data_buf = (int *) malloc(33000 * 4);
    rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
    rps->compressed_data = cdp;
    rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
    RED_decode(rps);
    cdp += rps->block_header->block_bytes;
    

    if (times_specified)
        offset_into_output_buffer = (int)((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
    else
        offset_into_output_buffer = channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample - start_samp;

    // copy requested samples from first block to output buffer
    // TBD this loop could be optimized
    for (i=0;i<rps->block_header->number_of_samples;i++)
    {
        if (offset_into_output_buffer < 0)
        {
            offset_into_output_buffer++;
            continue;
        }
        
        if (offset_into_output_buffer >= num_samps)
            break;
        
        *(decomp_data + offset_into_output_buffer) = temp_data_buf[i];
        
        offset_into_output_buffer++;
    }


    // decode bytes to samples
    sample_counter = offset_into_output_buffer;
    for (i=1;i<num_blocks-1;i++) {
        rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        // check that block fits fully within output array
        // this should be true, but it's possible a stray block exists out-of-order, or with a bad timestamp
        
        // we need to manually remove offset, since we are using the time value of the block bevore decoding the block
        // (normally the offset is removed during the decoding process)

        if (times_specified){

            block_start_time_offset = rps->block_header->start_time;
            remove_recording_time_offset( &block_start_time_offset );
            
            if (block_start_time_offset < start_time)
                continue;
            if (block_start_time_offset + ((rps->block_header->number_of_samples / channel->metadata.time_series_section_2->sampling_frequency) * 1e6) >= end_time)
                continue;
            rps->decompressed_ptr = rps->decompressed_data = decomp_data + (int)((((block_start_time_offset - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
        
        }else{
            rps->decompressed_ptr = rps->decompressed_data = decomp_data + sample_counter;
        }

        RED_decode(rps);
        cdp += rps->block_header->block_bytes;
        sample_counter += rps->block_header->number_of_samples;
    }
    
    if (num_blocks > 1)
    {
        // decode last block to temp array
        rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
        RED_decode(rps);
        
        offset_into_output_buffer = (int)((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
        
        if (times_specified)
            offset_into_output_buffer = (int)((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
        else
            offset_into_output_buffer = sample_counter;
        
        // copy requested samples from last block to output buffer
        // TBD this loop could be optimized
        for (i=0;i<rps->block_header->number_of_samples;i++)
        {
            if (offset_into_output_buffer < 0)
            {
                offset_into_output_buffer++;
                continue;
            }
            
            if (offset_into_output_buffer >= num_samps)
                break;
            
            *(decomp_data + offset_into_output_buffer) = temp_data_buf[i];
            
            offset_into_output_buffer++;
        }
    }
    
    // Numpy double type specific - no need to do this if we use numpy integer array and
    // put the data directly into it
    // DAN_CONSULT - is this efficient enough? I am basically reading the data into si4 array and copying it here
    // it i the NaN causing trouble. Let me know if you can think of anything to make it mor efficient.
    for (i=0;i<num_samps;i++){
        if (*(decomp_data + i) == RED_NAN)
            *(numpy_arr_data + i) = NPY_NAN;
        else
            *(numpy_arr_data + i) = (sf8) *(decomp_data + i);
    }
    free(decomp_data);


    // copy requested samples from last block to output buffer
    // we're done with the compressed data, get rid of it
    free (temp_data_buf);
    free (compressed_data_buffer);
    free (rps->difference_buffer);
    free (rps);
    
    if (channel->number_of_segments > 0)
        channel->segments[0].metadata_fps->directives.free_password_data = MEF_TRUE;
    free_channel(channel, MEF_TRUE);


    return (PyObject *) py_array_out;
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
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(tmd2->channel_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(tmd2_dict,"session_description");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(tmd2->session_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(tmd2_dict,"recording_duration");
    if (temp_o != NULL)
        tmd2->recording_duration = PyLong_AsLong(temp_o);

    // Time series specific fields
    temp_o = PyDict_GetItemString(tmd2_dict,"reference_description");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(vmd2->channel_description, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(vmd2_dict,"session_description");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
    PyObject    *temp_o;

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
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(md3->subject_name_1, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(md3_dict,"subject_name_2");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(md3->subject_name_2, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(md3_dict,"subject_ID");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(md3->subject_ID, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(md3_dict,"recording_location");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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

    // Assign from dict to struct
    temp_o = PyDict_GetItemString(rh_dict,"type_string");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(r_type->marker_name_1, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(Siez_type_dict,"marker_name_2");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(r_type->marker_name_2, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(Siez_type_dict,"annotation");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
    temp_UTF_str = NULL;
    if ((temp_o = PyDict_GetItemString(Siez_ch_type_dict,"name"))){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(r_type->name, temp_str_bytes);
    }

    if ((temp_o = PyDict_GetItemString(Siez_ch_type_dict,"onset")))
        r_type->onset = PyLong_AsLong(temp_o);

    if ((temp_o = PyDict_GetItemString(Siez_ch_type_dict,"offset")))
        r_type->offset = PyLong_AsLong(temp_o);

    return;
}   

void    map_python_CSti_type(PyObject *CSti_type_dict, MEFREC_CSti_1_0  *r_type)
{
    // Helpers
    PyObject    *temp_o, *temp_UTF_str;

    si1     *temp_str_bytes;

    // Assign from dict to struct
    temp_UTF_str = NULL;
    temp_o = PyDict_GetItemString(CSti_type_dict,"task_type");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(r_type->task_type, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(CSti_type_dict,"stimulus_duration");
    if (temp_o != NULL)
        r_type->stimulus_duration = PyLong_AsLong(temp_o);

    temp_o = PyDict_GetItemString(CSti_type_dict,"stimulus_type");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(r_type->stimulus_type, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(CSti_type_dict,"patient_response");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
    temp_UTF_str = NULL;
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
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(r_type->waveform, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(ESti_type_dict,"anode");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
        MEF_strcpy(r_type->anode, temp_str_bytes);
    }

    temp_o = PyDict_GetItemString(ESti_type_dict,"catode");
    if (temp_o != NULL){
        #if PY_MAJOR_VERSION >= 3
            temp_UTF_str = PyUnicode_AsEncodedString(temp_o, "utf-8","strict"); // Encode to UTF-8 python objects
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str); // Get the *char 
        #else
            temp_str_bytes = PyString_AS_STRING(temp_o);
        #endif
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
    si1   temp_str[256];
 
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
    si1   temp_str[256];
 
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
            Py_BuildValue("l", vmd->recording_duration));
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
    si1   temp_str[256];
    
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

PyObject *create_mef3_TOC(SEGMENT *segment)
{
    // Numpy TOC
    PyArrayObject *py_array_out;

    // Helpers
    TIME_SERIES_INDEX     *tsi;

    si8     number_of_entries;
    si8     prev_time, prev_sample, start_time, start_sample, samp_time_diff, seg_start_sample;
    sf8     fs;
    si8     *numpy_arr_data;
    npy_intp dims[2];

    si4     i;
    
 
    // initialize Numpy
    #if PY_MAJOR_VERSION >= 3
        import_array();
    #else
        init_numpy();
    #endif

    

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

        // Have we found a discontinuity?
        numpy_arr_data = (si8 *) PyArray_GETPTR2(py_array_out, 0, i);
        samp_time_diff = (si8) (((start_time - prev_time) - (1e6 * (start_sample - prev_sample)) / fs));
        if (samp_time_diff < (si8) (1e6/fs))
            samp_time_diff = 0;
        if  ((samp_time_diff != 0) | (i == 0)) // First entry is dicontinuity by definition
            *numpy_arr_data = 1;
        else
            *numpy_arr_data = 0;

        // Discontinuity duration
        numpy_arr_data = (si8 *) PyArray_GETPTR2(py_array_out, 1, i);
        *numpy_arr_data = (si8) samp_time_diff;

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
    PyObject *spec_dict;
    PyObject *uh_dict;
    PyObject *uhs_dict;
    PyObject *s1_dict;
    PyObject *s2_dict;
    PyObject *s3_dict;
    PyObject *sdi_dict;
    PyObject *idx_list;
    PyObject *TOC;
    
    // Helper variables
    si1   temp_str[256];
    ui4   i;
    si8   number_of_entries;
    
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

    // Set the TOC to NULL so that the logic works
    TOC = NULL;

    // Create indices dictionary
    
    switch (segment->channel_type){
        case TIME_SERIES_CHANNEL_TYPE:
            number_of_entries = segment->time_series_indices_fps->universal_header->number_of_entries;
            tsi = segment->time_series_indices_fps->time_series_indices;
            idx_list = PyList_New(number_of_entries);

            // Map time series indices - this takes long if indices are many
            if (map_indices_flag != 0){
                for(i = 0; i < number_of_entries; i++){
                    sdi_dict = map_mef3_ti(tsi);
                    PyList_SET_ITEM(idx_list, i, sdi_dict);
                    tsi++;
                }
                PyDict_SetItemString(metadata_dict, "indices", idx_list);
            }
            // Create TOC
            TOC = create_mef3_TOC(segment);
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
            PyDict_SetItemString(metadata_dict, "indices", idx_list);
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

PyObject *map_mef3_channel(CHANNEL *channel, si1 map_indices_flag) // This funtion also loops through segments
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
    
    METADATA_SECTION_1      *md1;
        TIME_SERIES_METADATA_SECTION_2  *tmd2;
        VIDEO_METADATA_SECTION_2    *vmd2;
    METADATA_SECTION_3      *md3;

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
    si1   temp_str[256];
    si4   i;
    
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
            PyDict_SetItemString(metadata_dict, "time_series_channels", PyDict_New()); 
            ts_dict = PyDict_GetItemString(metadata_dict, "time_series_channels");
        }
        // Get the channel pointer
        channel = session->time_series_channels + i;
        // Map the channel
        channel_dict = map_mef3_channel(channel, map_indices_flag);
        // Put into ditionary
        PyDict_SetItemString(ts_dict, channel->name, channel_dict); 

    }
    // Loop over video channels
    for (i = 0; i < session->number_of_video_channels; ++i){
        if (i == 0){
            PyDict_SetItemString(metadata_dict, "video_channels", PyDict_New()); 
            v_dict = PyDict_GetItemString(metadata_dict, "video_channels");
        }
        // Get the channel pointer
        channel = session->video_channels + i;
        // Map the channel
        channel_dict = map_mef3_channel(channel, map_indices_flag);
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
    // Dictionaries
    PyObject *rh_dict;
    PyObject *type_dict;

    // Helper variables
    si1   temp_str[256];
    ui4   type_code, *type_str_int;

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

    if (rh->encryption != ENCRYPTION_LEVEL_NO_ENTRY)
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

    type_str_int = (ui4 *) rh->type_string;
    type_code = *type_str_int;

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

    if (ri->encryption != ENCRYPTION_LEVEL_NO_ENTRY)
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
        // long_file_time = Py_BuildValue("l", ri->time);
        // PyDict_SetItemString(ri_dict, "time",
        //     ABS(long_file_time));
        PyDict_SetItemString(ri_dict, "time",
            Py_BuildValue("l", ri->time));
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

    MEFREC_CSti_1_0     *cog_stim;
        
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

    MEFREC_ESti_1_0     *el_stim;
        
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
    
    native_samp_freq = channel->metadata.time_series_section_2->sampling_frequency;
    prev_sample_number = channel->segments[0].metadata_fps->metadata.time_series_section_2->start_sample;
    prev_time = channel->segments[0].time_series_indices_fps->time_series_indices[0].start_time;
    
    for (j = 0; j < channel->number_of_segments; j++)
    {
        seg_start_sample = channel->segments[j].metadata_fps->metadata.time_series_section_2->start_sample;
        for (i = 0; i < channel->segments[j].metadata_fps->metadata.time_series_section_2->number_of_blocks; ++i) {
            if (channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time > uutc)
                goto done;
            prev_sample_number = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample;
            prev_time = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time;
        }
    }
    
    done:
        sample = prev_sample_number + (ui8) (((((sf8) (uutc - prev_time)) / 1000000.0) * native_samp_freq) + 0.5);
    
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

    for (j = 0; j < channel->number_of_segments; j++)
    {
        seg_start_sample = channel->segments[j].metadata_fps->metadata.time_series_section_2->start_sample;
        for (i = 0; i < channel->segments[j].metadata_fps->metadata.time_series_section_2->number_of_blocks; ++i){
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

void init_numpy(void)
{
    //Py_Initialize;
    #if PY_MAJOR_VERSION >= 3
        return;
    #else
        import_array();
        return;
    #endif
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
    #if PY_MAJOR_VERSION >= 3
        if (PyUnicode_Check(py_password_obj)){
            temp_UTF_str = PyUnicode_AsEncodedString(py_password_obj, "utf-8","strict");
            temp_str_bytes = PyBytes_AS_STRING(temp_UTF_str);
            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #else
        if (PyString_Check(py_password_obj)){
            temp_str_bytes = PyString_AS_STRING(py_password_obj);

            password = strcpy(password_arr,temp_str_bytes);
        }else{
            password = NULL;
        }
    #endif

    // Allocate universal header
    uh = (UNIVERSAL_HEADER *) calloc(1, sizeof(UNIVERSAL_HEADER));
    
    // Read file universal header
    fp = fopen(py_mef_file_path,"r");
    nb = fread((void *) uh, sizeof(UNIVERSAL_HEADER), 1, fp);
    if (nb != 1)
        printf("Error reading file");
    fclose(fp);

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
