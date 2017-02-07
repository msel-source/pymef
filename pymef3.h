// Python extension for opeartions with Multiscale Electrophysiology Format (MEF) version 3.0
// For licences please see meflib.c

// Written by Jan Cimbalnik

#include <Python.h>

#include "meflib.h"

#define EPSILON 0.0001
#define FLOAT_EQUAL(x,y) ( ((y - EPSILON) < x) && (x <( y + EPSILON)) )


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


// ---------- Mef3 to python dictionaries -----------
PyObject *map_mef3_universal_header(UNIVERSAL_HEADER *uh);

PyObject *map_mef3_md1(METADATA_SECTION_1 *md1);
PyObject *map_mef3_tmd2(TIME_SERIES_METADATA_SECTION_2 *tmd);
PyObject *map_mef3_vmd2(VIDEO_METADATA_SECTION_2 *vmd);
PyObject *map_mef3_md3(METADATA_SECTION_3 *md3);

PyObject *map_mef3_rh(RECORD_HEADER *rh);
PyObject *map_mef3_ri(RECORD_INDEX *ri);
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
PyObject *map_mef3_SyLg_type(RECORD_HEADER *rh);
