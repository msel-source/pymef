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
    "Function for writing MEF3 records at any level specified by path parameter.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path where record files will be written.\n\
     password_1: str\n\
        Level 1 password.\n\
     password_2: str\n\
        Level 2 password.\n\
     start_time: int\n\
        uUTC start time for universal header.\n\
     end_time: int\n\
        uUTC end time for universal header.\n\
     recording_offset: int\n\
        Offset for uUTC times.\n\
     record_list: list\n\
        List of record dictionaries consisting of numpy arrays.";

static char write_mef_ts_metadata_docstring[] =
    "Function to write MEF3 time series metadata file.\n\n\
     Parameters\n\
     ---------- \n\
     target_path: str\n\
        Path to segment being written.\n\
     password_1: str\n\
        Level 1 password.\n\
     password_2: str\n\
        Level 2 password.\n\
     start_time: int\n\
        uUTC recording start time - also used as recording time offset.\n\
     end_time: int\n\
        uUTC recording end time.\n\
     section_2_arr: np.array\n\
        Numpy array of time series metadata section 2 dtype.\n\
     section_3_arr: np.array\n\
        Numpy array of metadata section 3 dtype.";

static char write_mef_v_metadata_docstring[] =
    "Function to write MEF3 video metadata file.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to segment being written.\n\
     password_1: str\n\
        Level 1 password.\n\
     password_2: str\n\
        Level 2 password.\n\
     start_time: int\n\
        uUTC recording start time - also used as recording time offset.\n\
     end_time: int\n\
        uUTC recording end time.\n\
     section_2_arr: np.array\n\
        Numpy array of video metadata section 2 dtype.\n\
     section_3_arr: np.array\n\
        Numpy array of metadata section 3 dtype.";

static char write_mef_ts_data_and_indices_docstring[] =
    "Function to write MEF3 time series data and indices file.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to segment being written.\n\
     password_1: str\n\
        Level 1 password.\n\
     password_2: str\n\
        Level 2 password.\n\
     samples_per_mef_block: int\n\
        Number of samples in one MEF RED block.\n\
     raw_data: np.array\n\
        Numpy 1D array with raw data of dtype int32.\n\
     lossy_flag: bool\n\
        Flag for optional lossy compression (default=False).";

static char write_mef_v_indices_docstring[] =
    "Function to write MEF3 video indices file.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to segment being written.\n\
     password_1: str\n\
        Level 1 password.\n\
     password_2: str\n\
        Level 2 password.\n\
     start_time: int\n\
        uUTC recording start time - minimum value of index entries.\n\
     end_time: int\n\
        uUTC recording end time - maximum value of index entries.\n\
     index_entries: list\n\
        List of numpy arrays with index entries.";

/* Documentation to be read in Python - append functions*/
static char append_ts_data_and_indices_docstring[] =
    "Function to append MEF3 time series data and indices to existing segment files.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to segment being appended\n\
     password_1: str\n\
        Level 1 password.\n\
     password_2: str\n\
        Level 2 password.\n\
     samples_per_mef_block: int\n\
        Number of samples in one MEF RED block.\n\
     raw_data: np.array\n\
        Numpy 1D array with raw data of dtype int32.\n\
     discontinuity_flag: bool\n\
        Flag to mark discontinuity at the start of appended data (default=True)\n\
     lossy_flag: bool\n\
        Flag for optional lossy compression (default=False).";

/* Documentation to be read in Python - read functions*/
static char read_mef_ts_data_docstring[] =
    "Function to read MEF3 time series data.\n\n\
     Parameters\n\
     ----------\n\
     channel_specific_metadata: np.ndarray\n\
        Channel metadata\n\
     start: int\n\
        Start sample or uUTC time to be read.\n\
     end: int\n\
        End sample or uUTC time to be read.\n\
     time_flag: bool\n\
        Flag to indicate if user is reading by samples or uUTC times (default=False - reading by sample)\n\n\
     Returns\n\
     -------\n\
     data: np.array\n\
        1D numpy array (dtype=float) with data. If the data is read by uUTC and a gap is present the missing values are filled with NaNs";

static char read_mef_session_metadata_docstring[] =
    "Function to read MEF3 session metadata.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to MEF3 session being read\n\
     password: str\n\
        Level 1 or level 2 password.\n\
     map_indices_flag: bool\n\
        Flag to enable the mapping of the time-series and video indices (default=True, map indices)\n\
     copy_metadata_to_dict: bool\n\
        Flag to copy metadata into a python dictionary structure (True), instead of returning the metadata by reference in Numpy structured datatypes (Default=False)\n\n\
     Returns\n\
     -------\n\
     session_metadata: dict\n\
        Dictionary with session metadata and all channels and segments metadata and records.";

static char read_mef_channel_metadata_docstring[] =
    "Function to read MEF3 channel metadata.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to MEF3 channel being read.\n\
     password: str\n\
        Level 1 or level 2 password.\n\
     map_indices_flag: bool\n\
        Flag to enable the mapping of the time-series and video indices (default=True, map indices)\n\
     copy_metadata_to_dict: bool\n\
        Flag to copy metadata into a python dictionary structure (True), instead of returning the metadata by reference in Numpy structured datatypes (Default=False)\n\n\
     Returns\n\
     -------\n\
     channel_metadata: dict\n\
        Dictionary with channel metadata and all segments metadata and records.";

static char read_mef_segment_metadata_docstring[] =
    "Function to read MEF3 segment metadata.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to MEF3 segment being read.\n\
     password: str\n\
        Level 1 or level 2 password.\n\
     map_indices_flag: bool\n\
        Flag to enable the mapping of the time-series and video indices (default=True, map indices)\n\
     copy_metadata_to_dict: bool\n\
        Flag to copy metadata into a python dictionary structure (True), instead of returning the metadata by reference in Numpy structured datatypes (Default=False)\n\n\
     Returns\n\
     -------\n\
     segment_metadata: dict\n\
        Dictionary with segment metadata and records.";

/* Documentation to be read in Python - helper functions*/
static char check_mef_password_docstring[] =
    "Function to check MEF3 password validity.\n\n\
     Parameters\n\
     ----------\n\
     target_path: str\n\
        Path to MEF3 metadata file.\n\
     password: str\n\
        Level 1 or level 2 password.\n\
     Returns\n\
     -------\n\
     password_type: int\n\
        - 0 - incorrect password\n\
        - 1 - level 1 password\n\
        - 2 - level 2 password\n";

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
static PyObject *read_mef_session_metadata(PyObject *self, PyObject *args, PyObject* kwargs);
static PyObject *read_mef_channel_metadata(PyObject *self, PyObject *args, PyObject* kwargs);
static PyObject *read_mef_segment_metadata(PyObject *self, PyObject *args, PyObject* kwargs);

/* Pyhon object declaration - clean functions*/
static PyObject *clean_mef_session_metadata(PyObject *self, PyObject *args);
static PyObject *clean_mef_channel_metadata(PyObject *self, PyObject *args);
static PyObject *clean_mef_segment_metadata(PyObject *self, PyObject *args);

/* Python object declaration - helper functions */
static PyObject *check_mef_password(PyObject *self, PyObject *args);

/* Python object declaration - numpy data types */
static PyObject *create_rh_dtype(void);
static PyObject *create_ri_dtype(void);
static PyObject *create_edfa_dtype(PyObject *self, PyObject *args);
static PyObject *create_edfa_dtype_c(ui4 text_len);
static PyObject *create_lntp_dtype(PyObject *self, PyObject *args);
static PyObject *create_lntp_dtype_c(ui4 template_len);
static PyObject *create_note_dtype(PyObject *self, PyObject *args);
static PyObject *create_note_dtype_c(ui4 text_len);
static PyObject *create_seiz_dtype(void);
static PyObject *create_seiz_ch_dtype(void);
static PyObject *create_sylg_dtype(PyObject *self, PyObject *args);
static PyObject *create_sylg_dtype_c(ui4 text_len);
static PyObject *create_csti_dtype(void);
static PyObject *create_esti_dtype(void);
static PyObject *create_curs_dtype(void);
static PyObject *create_epoc_dtype(void);

static PyObject *create_uh_dtype(void);
static PyObject *create_md1_dtype(void);
static PyObject *create_tmd2_dtype(void);
static PyObject *create_vmd2_dtype(void);
static PyObject *create_md3_dtype(void);
static PyObject *create_ti_dtype(void);
static PyObject *create_vi_dtype(void);

static PyObject *create_segment_dtype(void);
static PyObject *create_channel_dtype(void);
static PyObject *create_session_dtype(void);

/* Thread safety for free-threaded Python.
 *
 * meflib keeps its state in the process-wide MEF_globals struct, which every
 * entry point re-initializes, modifies, and frees. With the GIL these calls are
 * already serialized; without it, every function that touches meflib must hold
 * this lock. PyMutex_Lock detaches the thread state while waiting, so it cannot
 * deadlock with stop-the-world pauses. */
#ifdef Py_GIL_DISABLED
static PyMutex pymef_lock = {0};

#define PYMEF_LOCKED(func) \
    static PyObject *func##_locked(PyObject *self, PyObject *args) { \
        PyObject *result; \
        PyMutex_Lock(&pymef_lock); \
        result = func(self, args); \
        PyMutex_Unlock(&pymef_lock); \
        return result; \
    }
#define PYMEF_LOCKED_KW(func) \
    static PyObject *func##_locked(PyObject *self, PyObject *args, PyObject *kwargs) { \
        PyObject *result; \
        PyMutex_Lock(&pymef_lock); \
        result = func(self, args, kwargs); \
        PyMutex_Unlock(&pymef_lock); \
        return result; \
    }
#define PYMEF_FUNC(func) func##_locked
#else
#define PYMEF_LOCKED(func)
#define PYMEF_LOCKED_KW(func)
#define PYMEF_FUNC(func) func
#endif

PYMEF_LOCKED(write_mef_data_records)
PYMEF_LOCKED(write_mef_ts_metadata)
PYMEF_LOCKED(write_mef_v_metadata)
PYMEF_LOCKED(write_mef_ts_data_and_indices)
PYMEF_LOCKED(write_mef_v_indices)
PYMEF_LOCKED(append_ts_data_and_indices)
PYMEF_LOCKED(read_mef_ts_data)
PYMEF_LOCKED_KW(read_mef_session_metadata)
PYMEF_LOCKED_KW(read_mef_channel_metadata)
PYMEF_LOCKED_KW(read_mef_segment_metadata)
PYMEF_LOCKED(clean_mef_session_metadata)
PYMEF_LOCKED(clean_mef_channel_metadata)
PYMEF_LOCKED(clean_mef_segment_metadata)
PYMEF_LOCKED(check_mef_password)

/* Python-callable adapters for the argument-less dtype functions (METH_NOARGS) */
#define PYMEF_NOARGS(func) \
    static PyObject *func##_noargs(PyObject *self, PyObject *Py_UNUSED(args)) { \
        return func(); \
    }

PYMEF_NOARGS(create_rh_dtype)
PYMEF_NOARGS(create_ri_dtype)
PYMEF_NOARGS(create_seiz_dtype)
PYMEF_NOARGS(create_seiz_ch_dtype)
PYMEF_NOARGS(create_csti_dtype)
PYMEF_NOARGS(create_esti_dtype)
PYMEF_NOARGS(create_curs_dtype)
PYMEF_NOARGS(create_epoc_dtype)
PYMEF_NOARGS(create_uh_dtype)
PYMEF_NOARGS(create_md1_dtype)
PYMEF_NOARGS(create_tmd2_dtype)
PYMEF_NOARGS(create_vmd2_dtype)
PYMEF_NOARGS(create_md3_dtype)
PYMEF_NOARGS(create_ti_dtype)
PYMEF_NOARGS(create_vi_dtype)
PYMEF_NOARGS(create_segment_dtype)
PYMEF_NOARGS(create_channel_dtype)
PYMEF_NOARGS(create_session_dtype)

/* Specification of the members of the module */
static PyMethodDef module_methods[] = {
    {"write_mef_data_records", PYMEF_FUNC(write_mef_data_records), METH_VARARGS, write_mef_data_records_docstring},
    {"write_mef_ts_metadata", PYMEF_FUNC(write_mef_ts_metadata), METH_VARARGS, write_mef_ts_metadata_docstring},
    {"write_mef_v_metadata", PYMEF_FUNC(write_mef_v_metadata), METH_VARARGS, write_mef_v_metadata_docstring},
    {"write_mef_ts_data_and_indices", PYMEF_FUNC(write_mef_ts_data_and_indices), METH_VARARGS, write_mef_ts_data_and_indices_docstring},
    {"write_mef_v_indices", PYMEF_FUNC(write_mef_v_indices), METH_VARARGS, write_mef_v_indices_docstring},
    {"append_ts_data_and_indices", PYMEF_FUNC(append_ts_data_and_indices), METH_VARARGS, append_ts_data_and_indices_docstring},
    {"read_mef_ts_data", PYMEF_FUNC(read_mef_ts_data), METH_VARARGS, read_mef_ts_data_docstring},
    {"read_mef_session_metadata", (PyCFunction)PYMEF_FUNC(read_mef_session_metadata), METH_VARARGS | METH_KEYWORDS, read_mef_session_metadata_docstring},
    {"read_mef_channel_metadata", (PyCFunction)PYMEF_FUNC(read_mef_channel_metadata), METH_VARARGS | METH_KEYWORDS, read_mef_channel_metadata_docstring},
    {"read_mef_segment_metadata", (PyCFunction)PYMEF_FUNC(read_mef_segment_metadata), METH_VARARGS | METH_KEYWORDS, read_mef_segment_metadata_docstring},
    {"clean_mef_session_metadata", PYMEF_FUNC(clean_mef_session_metadata), METH_VARARGS, NULL},
    {"clean_mef_channel_metadata", PYMEF_FUNC(clean_mef_channel_metadata), METH_VARARGS, NULL},
    {"clean_mef_segment_metadata", PYMEF_FUNC(clean_mef_segment_metadata), METH_VARARGS, NULL},
    {"check_mef_password", PYMEF_FUNC(check_mef_password), METH_VARARGS, check_mef_password_docstring},

    // New numpy stuff
    {"create_rh_dtype", create_rh_dtype_noargs, METH_NOARGS, NULL},
    {"create_ri_dtype", create_ri_dtype_noargs, METH_NOARGS, NULL},
    {"create_edfa_dtype", create_edfa_dtype, METH_VARARGS, NULL},
    {"create_lntp_dtype", create_lntp_dtype, METH_VARARGS, NULL},
    {"create_note_dtype", create_note_dtype, METH_VARARGS, NULL},
    {"create_seiz_dtype", create_seiz_dtype_noargs, METH_NOARGS, NULL},
    {"create_seiz_ch_dtype", create_seiz_ch_dtype_noargs, METH_NOARGS, NULL},
    {"create_sylg_dtype", create_sylg_dtype, METH_VARARGS, NULL},
    {"create_csti_dtype", create_csti_dtype_noargs, METH_NOARGS, NULL},
    {"create_esti_dtype", create_esti_dtype_noargs, METH_NOARGS, NULL},
    {"create_curs_dtype", create_curs_dtype_noargs, METH_NOARGS, NULL},
    {"create_epoc_dtype", create_epoc_dtype_noargs, METH_NOARGS, NULL},


    {"create_uh_dtype", create_uh_dtype_noargs, METH_NOARGS, NULL},
    {"create_md1_dtype", create_md1_dtype_noargs, METH_NOARGS, NULL},
    {"create_tmd2_dtype", create_tmd2_dtype_noargs, METH_NOARGS, NULL},
    {"create_vmd2_dtype", create_vmd2_dtype_noargs, METH_NOARGS, NULL},
    {"create_md3_dtype", create_md3_dtype_noargs, METH_NOARGS, NULL},
    {"create_ti_dtype", create_ti_dtype_noargs, METH_NOARGS, NULL},
    {"create_vi_dtype", create_vi_dtype_noargs, METH_NOARGS, NULL},

    {"create_segment_dtype", create_segment_dtype_noargs, METH_NOARGS, NULL},
    {"create_channel_dtype", create_channel_dtype_noargs, METH_NOARGS, NULL},
    {"create_session_dtype", create_session_dtype_noargs, METH_NOARGS, NULL},

    {NULL, NULL, 0, NULL}
};

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

#ifdef Py_GIL_DISABLED
    /* Safe without the GIL: all meflib access is serialized by pymef_lock */
    if (PyUnstable_Module_SetGIL(m, Py_MOD_GIL_NOT_USED) < 0) {
        Py_DECREF(m);
        return NULL;
    }
#endif

    return m;
}

/* Function declarations */

// ---------- Python dictionaries to mef3 -----------
// UNIVERSAL_HEADER *map_python_uh(PyObject *);

// METADATA_SECTION_1 *map_python_md1(PyObject *);
void    map_python_tmd2(PyObject *tmd2_arr, TIME_SERIES_METADATA_SECTION_2 *tmd2);
void    map_python_vmd2(PyObject *vmd2_arr, VIDEO_METADATA_SECTION_2 *vmd2);
void    map_python_vi(PyObject *vi_arr, VIDEO_INDEX *vi);
void    map_python_md3(PyObject *md3_arr, METADATA_SECTION_3 *md3);


// Mef record structures
void    map_python_rh(PyObject *rh_dict, RECORD_HEADER  *rh);
void    map_python_EDFA_type(PyObject *EDFA_type_dict, MEFREC_EDFA_1_0  *r_type);
void    map_python_LNTP_type(PyObject *LNTP_type_dict, MEFREC_LNTP_1_0  *r_type);
void    map_python_Siez_type(PyObject *Siez_type_dict, MEFREC_Seiz_1_0  *r_type);
void    map_python_Siez_ch_type(PyObject *Siez_ch_type_dict, si1 *r_type);
void    map_python_Siez_type(PyObject *Siez_type_dict, MEFREC_Seiz_1_0  *r_type);
void    map_python_CSti_type(PyObject *CSti_type_dict, MEFREC_CSti_1_0  *r_type);
void    map_python_ESti_type(PyObject *ESti_type_dict, MEFREC_ESti_1_0  *r_type);
void    map_python_Curs_type(PyObject *Curs_type_dict, MEFREC_Curs_1_0  *r_type);
void    map_python_Epoc_type(PyObject *Epoc_type_dict, MEFREC_Epoc_1_0  *r_type);


// ---------- Mef3 to python -----------

PyObject *map_mef3_decode_maxbytes_to_string(const char *s, size_t max_size);
PyObject *map_mef3_decode_sizebytes_to_string(const char *s, size_t size);

PyObject *map_mef3_uh(UNIVERSAL_HEADER *uh, si1 copy_metadata_to_dict);

PyObject *map_mef3_md1(METADATA_SECTION_1 *md1, si1 copy_metadata_to_dict);
PyObject *map_mef3_tmd2(TIME_SERIES_METADATA_SECTION_2 *tmd, si1 copy_metadata_to_dict);
PyObject *map_mef3_vmd2(VIDEO_METADATA_SECTION_2 *vmd, si1 copy_metadata_to_dict);
PyObject *map_mef3_md3(METADATA_SECTION_3 *md3, si1 copy_metadata_to_dict);


PyObject *map_mef3_ti(TIME_SERIES_INDEX *ti, si8 number_of_entries, si1 copy_metadata_to_dict);
PyObject *map_mef3_vi(VIDEO_INDEX *vi, si8 number_of_entries, si1 copy_metadata_to_dict);
PyObject *create_mef3_TOC(SEGMENT *segment);

PyObject *map_mef3_segment(SEGMENT *segment, si1 map_indices_flag, si1 copy_metadata_to_dict);
PyObject *map_mef3_channel(CHANNEL *channel, si1 map_indices_flag, si1 copy_metadata_to_dict);
PyObject *map_mef3_session(SESSION *session, si1 map_indices_flag, si1 copy_metadata_to_dict);

// Mef record structures
PyObject *map_mef3_records(FILE_PROCESSING_STRUCT *ri_fps, FILE_PROCESSING_STRUCT *rd_fps, si1 copy_metadata_to_dict);

PyObject *map_mef3_rh(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_ri(RECORD_INDEX *ri);

PyObject *map_mef3_Note_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_EDFA_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_LNTP_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_Seiz_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_Seiz_ch_type(RECORD_HEADER *rh, si4 number_of_channels, si1 copy_metadata_to_dict);
PyObject *map_mef3_CSti_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_ESti_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_SyLg_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_Curs_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);
PyObject *map_mef3_Epoc_type(RECORD_HEADER *rh, si1 copy_metadata_to_dict);

// Helper functions
si4 check_block_crc(ui1* block_hdr_ptr, ui4 max_samps, ui1* total_data_ptr, ui8 total_data_bytes);
si4 extract_segment_number(si1 *segment_name);
si8 sample_for_uutc_c(si8 uutc, CHANNEL *channel);
si8 uutc_for_sample_c(si8 sample, CHANNEL *channel);
void memset_int(si4 *ptr, si4 value, size_t num);
void init_numpy(void);
