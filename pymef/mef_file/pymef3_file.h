// Written by Jan Cimbalnik

#include <Python.h>

#include "meflib.h"

#define EPSILON 0.0001
#define FLOAT_EQUAL(x,y) ( ((y - EPSILON) < x) && (x <( y + EPSILON)) )
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


/* Python methods definitions and help */

static char pymef3_file_docstring[] =
    "This submodule provides a wrapper around Mayo Electrophysiology Format (MEF) version 3.0 library.";

/* Documentation to be read in Python - write functions*/
static char write_mef_data_records_docstring[] =
    "Function for writing MEF3 records at any level specified by path parameter.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path where record files will be written (str)\n \
     password_1 - level 1 password (str)\n \
     password_2 - level 2 password (str)\n \
     start_time - uUTC start time for universal header (int)\n \
     end_time - uUTC end time for universal header (int)\n \
     recording_offset - offset for uUTC times (int)\n \
     record_list - list of record dictionaries (list)\n";

static char write_mef_ts_metadata_docstring[] =
    "Function to write MEF3 time series metadata file.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to segment being written (str)\n \
     password_1 - level 1 password (str)\n \
     password_2 - level 2 password (str)\n \
     start_time - uUTC recording start time - also used as recording time offset (int)\n \
     end_time - uUTC recording end time (int)\n \
     section_2_dictionary - dictionary with section 2 metadata values (dictionary)\n \
     section_3_dictionary - dictionary with section 3 metadata values (dictionary)\n";

static char write_mef_v_metadata_docstring[] =
    "Function to write MEF3 video metadata file.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to segment being written (str)\n \
     password_1 - level 1 password (str)\n \
     password_2 - level 2 password (str)\n \
     start_time - uUTC recording start time - also used as recording time offset (int)\n \
     end_time - uUTC recording end time (int)\n \
     section_2_dictionary - dictionary with section 2 metadata values (dictionary)\n \
     section_3_dictionary - dictionary with section 3 metadata values (dictionary)\n";

static char write_mef_ts_data_and_indices_docstring[] =
    "Function to write MEF3 time series data and indices file.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to segment being written (str)\n \
     password_1 - level 1 password (str)\n \
     password_2 - level 2 password (str)\n \
     samples_per_mef_block - number of samples in one MEF RED block (int)\n \
     raw_data - numpy 1D array with raw data (numpy.array, int32)\n \
     lossy_flag (optional) - flag for optional lossy compression (bool, default=False)\n";

static char write_mef_v_indices_docstring[] =
    "Function to write MEF3 video indices file.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to segment being written (str)\n \
     password_1 - level 1 password (str)\n \
     password_2 - level 2 password (str)\n \
     start_time - uUTC recording start time - minimum value of index entries (int)\n \
     end_time - uUTC recording end time - maximum valu of index entries (int)\n \
     index_entries - list of dictionaries with index entries (list)\n";

/* Documentation to be read in Python - append functions*/
static char append_ts_data_and_indices_docstring[] =
    "Function to append MEF3 time series data and indices to existing segment files.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to segment being appended (str)\n \
     password_1 - level 1 password (str)\n \
     password_2 - level 2 password (str)\n \
     samples_per_mef_block - number of samples in one MEF RED block (int)\n \
     raw_data - numpy 1D array with raw data (numpy.array, int32)\n \
     discontinuity_flag (optional) - flag to mark discontinuity at the start of appended data (bool, default=True)\n \
     lossy_flag (optional) - flag for optional lossy compression (bool, default=False)\n";

/* Documentation to be read in Python - read functions*/
static char read_mef_ts_data_docstring[] =
    "Function to read MEF3 time series data.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to time series channel being read (str)\n \
     password - level 1 or level 2 password (str)\n \
     start_sample(time) - start sample or uUTC time to be read (int)\n \
     end_sample(time) - end sample or uUTC time to be read (int)\n \
     time_flag (optional) - flag to indicate if user is reading by samples or uUTC times (bool, default=False - reading by samples)\n\n \
     Returns:\n \
     --------\n \
     data - 1D numpy array (float) with data. If the data is read by uUTC and a gap is present the missing values are filled with NaNs\n";

static char read_mef_session_metadata_docstring[] =
    "Function to read MEF3 session metadata.\n\n \
     Parameters: \
     ----------- \
     target_path - path to MEF3 session being read (str)\n \
     password - level 1 or level 2 password (str)\n \
     Returns:\n \
     --------\n \
     session_metadata - dictionary with session metadata and all channels and segments metadata and records\n";

static char read_mef_channel_metadata_docstring[] =
    "Function to read MEF3 channel metadata.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to MEF3 channel being read (str)\n \
     password - level 1 or level 2 password (str)\n \
     Returns:\n \
     --------\n \
     channel_metadata - dictionary with channel metadata and all segments metadata and records\n";

static char read_mef_segment_metadata_docstring[] =
    "Function to read MEF3 segment metadata.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to MEF3 segment being read (str)\n \
     password - level 1 or level 2 password (str)\n \
     Returns:\n \
     --------\n \
     segment_metadata - dictionary with segment metadata and records\n";

/* Documentation to be read in Python - helper functions*/
static char check_mef_password_docstring[] =
    "Function to check MEF3 password validity.\n\n \
     Parameters:\n \
     -----------\n \
     target_path - path to MEF3 metadata file (str)\n \
     password - level 1 or level 2 password (str)\n \
     Returns:\n \
     --------\n \
     password_type - 0 - incorrect password, 1 - level 1 password, 2 - level 2 password\n";

/* Pyhon object declaration - write functions*/
static PyObject *write_mef_data_records(PyObject *self, PyObject *args);
static PyObject *write_mef_ts_metadata(PyObject *self, PyObject *args);
static PyObject *write_mef_v_metadata(PyObject *self, PyObject *args);
static PyObject *write_mef_ts_data_and_indices(PyObject *self, PyObject *args);
static PyObject *write_mef_v_indices(PyObject *self, PyObject *args);

/* Pyhon object declaration - append functions*/
static PyObject *append_ts_data_and_indices(PyObject *self, PyObject *args);

/* Pyhon object declaration - read functions*/
static PyObject *read_mef_ts_data(PyObject *self, PyObject *args);

static PyObject *read_mef_session_metadata(PyObject *self, PyObject *args);
static PyObject *read_mef_channel_metadata(PyObject *self, PyObject *args);
static PyObject *read_mef_segment_metadata(PyObject *self, PyObject *args);

/* Python object declaration - helper functions */
static PyObject *check_mef_password(PyObject *self, PyObject *args);

/* Specification of the members of the module */
static PyMethodDef module_methods[] = {
    {"write_mef_data_records", write_mef_data_records, METH_VARARGS, write_mef_data_records_docstring},
    {"write_mef_ts_metadata", write_mef_ts_metadata, METH_VARARGS, write_mef_ts_metadata_docstring},
    {"write_mef_v_metadata", write_mef_v_metadata, METH_VARARGS, write_mef_v_metadata_docstring},
    {"write_mef_ts_data_and_indices", write_mef_ts_data_and_indices, METH_VARARGS, write_mef_ts_data_and_indices_docstring},
    {"write_mef_v_indices", write_mef_v_indices, METH_VARARGS, write_mef_v_indices_docstring},
    {"append_ts_data_and_indices", append_ts_data_and_indices, METH_VARARGS, append_ts_data_and_indices_docstring},
    {"read_mef_ts_data", read_mef_ts_data, METH_VARARGS, read_mef_ts_data_docstring},
    {"read_mef_session_metadata", read_mef_session_metadata, METH_VARARGS, read_mef_session_metadata_docstring},
    {"read_mef_channel_metadata", read_mef_channel_metadata, METH_VARARGS, read_mef_channel_metadata_docstring},
    {"read_mef_segment_metadata", read_mef_segment_metadata, METH_VARARGS, read_mef_segment_metadata_docstring},
    {"check_mef_password", check_mef_password, METH_VARARGS, check_mef_password_docstring},
    {NULL, NULL, 0, NULL}
};

//  Fork for python 3 and python 2
#if PY_MAJOR_VERSION >= 3

    /* Definition of struct for python 3 */
    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "pymef.mef_file.pymef3_file",     /* m_name */
        pymef3_file_docstring,  /* m_doc */
        -1,                  /* m_size */
        module_methods,    /* m_methods */
        NULL,                /* m_reload */
        NULL,                /* m_traverse */
        NULL,                /* m_clear */
        NULL,                /* m_free */
    };

    /* Module initialisation */
    PyObject * PyInit_pymef3_file(void)
    {
        PyObject *m = PyModule_Create(&moduledef);

        if (m == NULL)
            return NULL;

        return m;
    }
#else
    /* Module initialisation */
    PyMODINIT_FUNC initpymef3_file(void)
    {
        PyObject *m = Py_InitModule3("pymef3_file", module_methods, pymef3_file_docstring);
        if (m == NULL)
            return;

    }
#endif

/* Function declarations */

// ---------- Python dictionaries to mef3 -----------
// UNIVERSAL_HEADER *map_python_uh(PyObject *);

// METADATA_SECTION_1 *map_python_md1(PyObject *);
void    map_python_tmd2(PyObject *tmd2_dict, TIME_SERIES_METADATA_SECTION_2 *tmd2);
void    map_python_vmd2(PyObject *vmd2_dict, VIDEO_METADATA_SECTION_2 *vmd2);
void    map_python_vi(PyObject *vi_dict, VIDEO_INDEX *vi);
void    map_python_md3(PyObject *md3_dict, METADATA_SECTION_3 *md3);

// Mef record structures
void    map_python_rh(PyObject *rh_dict, RECORD_HEADER  *rh);
void    map_python_EDFA_type(PyObject *EDFA_type_dict, MEFREC_EDFA_1_0  *r_type);
void    map_python_LNTP_type(PyObject *LNTP_type_dict, MEFREC_LNTP_1_0  *r_type);
void    map_python_Siez_type(PyObject *Siez_type_dict, MEFREC_Seiz_1_0  *r_type);
void    map_python_Siez_type_channel(PyObject *Siez_ch_type_dict, MEFREC_Seiz_1_0_CHANNEL *r_type);
void    map_python_Siez_type(PyObject *Siez_type_dict, MEFREC_Seiz_1_0  *r_type);
void    map_python_CSti_type(PyObject *CSti_type_dict, MEFREC_CSti_1_0  *r_type);
void    map_python_ESti_type(PyObject *ESti_type_dict, MEFREC_ESti_1_0  *r_type);


// ---------- Mef3 to python dictionaries -----------
PyObject *map_mef3_uh(UNIVERSAL_HEADER *uh);

PyObject *map_mef3_md1(METADATA_SECTION_1 *md1);
PyObject *map_mef3_tmd2(TIME_SERIES_METADATA_SECTION_2 *tmd);
PyObject *map_mef3_vmd2(VIDEO_METADATA_SECTION_2 *vmd);
PyObject *map_mef3_md3(METADATA_SECTION_3 *md3);

PyObject *map_mef3_ti(TIME_SERIES_INDEX *ti);
PyObject *map_mef3_vi(VIDEO_INDEX *vi);
PyObject *create_mef3_TOC(SEGMENT *segment);

PyObject *map_mef3_segment(SEGMENT *segment, si1 map_indices_flag);
PyObject *map_mef3_channel(CHANNEL *channel, si1 map_indices_flag);
PyObject *map_mef3_session(SESSION *session, si1 map_indices_flag);

// Mef record structures
PyObject *map_mef3_records(FILE_PROCESSING_STRUCT *ri_fps, FILE_PROCESSING_STRUCT *rd_fps);

PyObject *map_mef3_rh(RECORD_HEADER *rh);
PyObject *map_mef3_ri(RECORD_INDEX *ri);

PyObject *map_mef3_Note_type(RECORD_HEADER *rh);
PyObject *map_mef3_EDFA_type(RECORD_HEADER *rh);
PyObject *map_mef3_LNTP_type(RECORD_HEADER *rh);
PyObject *map_mef3_Seiz_type(RECORD_HEADER *rh);
PyObject *map_mef3_CSti_type(RECORD_HEADER *rh);
PyObject *map_mef3_ESti_type(RECORD_HEADER *rh);
PyObject *map_mef3_SyLg_type(RECORD_HEADER *rh);

// Helper functions
si4 extract_segment_number(si1 *segment_name);
si8 sample_for_uutc_c(si8 uutc, CHANNEL *channel);
si8 uutc_for_sample_c(si8 sample, CHANNEL *channel);
void memset_int(si4 *ptr, si4 value, size_t num);
void init_numpy(void);


