// Written by Jan Cimbalnik

#include <Python.h>

#include "meflib.h"

#define EPSILON 0.0001
#define FLOAT_EQUAL(x,y) ( ((y - EPSILON) < x) && (x <( y + EPSILON)) )


/* Python methods definitions and help */

static char pymef3_file_docstring[] =
    "This module provides an interface for reading .mef (v 3.x) files.";

/* Documentation to be read in Python - write functions*/
static char write_mef_data_records_docstring[] =
    "Writes .mefd session directory + data and indices files.\n write_mef_session(path_to_session, session_name, level_1_password, level_2_password, uutc_rec_start, uutc_rec_stop, recording_note)";
static char write_mef_ts_metadata_docstring[] =
    "Writes .timd time series directory and .segd segment directory along with the data and indices files. Help to be written";
static char write_mef_v_metadata_docstring[] =
    "Writes .timd time series directory and .segd segment directory along with the data and indices files. Help to be written";
static char write_mef_ts_data_and_indices_docstring[] =
    "Writes .timd time series directory and .segd segment directory along with the data and indices files. Help to be written";
static char write_mef_v_indices_docstring[] =
    "Writes .timd time series directory and .segd segment directory along with the data and indices files. Help to be written";
    
/* Documentation to be read in Python - append functions*/
static char append_ts_data_and_indices_docstring[] =
    "Appends ts data";

/* Documentation to be read in Python - read functions*/
static char read_mef_ts_data_docstring[] =
    "Reads .timd time series directory";
static char read_mef_session_metadata_docstring[] =
    "Reads metadata of a mef session";
static char read_mef_channel_metadata_docstring[] =
    "Reads metadata of a mef channel";
static char read_mef_segment_metadata_docstring[] =
    "Reads metadata of a mef segment";

/* Documentation to be read in Python - helper functions*/
static char check_mef_password_docstring[] =
    "Checks the the mef password";

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


/* Definition of struct for python 3 */
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "pymef.mef_file.pymef3_file",     /* m_name */
    "This module provides an interface operations with MEF3 file format",  /* m_doc */
    -1,                  /* m_size */
    module_methods,    /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};

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

PyObject *map_mef3_segment(SEGMENT *segment);
PyObject *map_mef3_channel(CHANNEL *channel);
PyObject *map_mef3_session(SESSION *session);

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
void memset_int(si4 *ptr, si4 value, size_t num);
si8 uutc_for_sample_c(si8 sample, CHANNEL *channel);

