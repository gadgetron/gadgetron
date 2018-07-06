#ifndef GADGETRON_PYTHON_NUMPY_WRAPPERS_H
#define GADGETRON_PYTHON_NUMPY_WRAPPERS_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <ismrmrd/waveform.h>
#include "numpy/ndarraytypes.h"

namespace Gadgetron {

/// Wrappers for NumPy C-API functions. These functions must be used
/// in the same C++ source file as the call to `import_array()`. In this
/// case, that is Python.cpp. The simplest solution is to lightly wrap the API.
EXPORTPYTHON int NumPyArray_NDIM(PyObject* obj);
EXPORTPYTHON npy_intp NumPyArray_DIM(PyObject* obj, int i);
EXPORTPYTHON void *NumPyArray_DATA(PyObject* obj);
EXPORTPYTHON npy_intp *NumPyArray_STRIDES(PyObject* obj);
EXPORTPYTHON int NumPyArray_ITEMSIZE(PyObject* obj);
EXPORTPYTHON npy_intp NumPyArray_SIZE(PyObject* obj);
EXPORTPYTHON PyObject *NumPyArray_SimpleNew(int nd, npy_intp* dims, int typenum);
EXPORTPYTHON PyObject *NumPyArray_EMPTY(int nd, npy_intp* dims, int typenum, int fortran);
EXPORTPYTHON PyObject* NumPyArray_FromAny(PyObject* op, PyArray_Descr* dtype, int min_depth, int max_depth, int requirements, PyObject* context);
/// return the enumerated numpy type for a given C++ type
template <typename T> int get_numpy_type() { return NPY_VOID; }
template <> inline int get_numpy_type< bool >() { return NPY_BOOL; }
template <> inline int get_numpy_type< char >() { return NPY_INT8; }
template <> inline int get_numpy_type< unsigned char >() { return NPY_UINT8; }
template <> inline int get_numpy_type< short >() { return NPY_INT16; }
template <> inline int get_numpy_type< unsigned short >() { return NPY_UINT16; }
template <> inline int get_numpy_type< int >() { return NPY_INT32; }
template <> inline int get_numpy_type< unsigned int >() { return NPY_UINT32; }
template <> inline int get_numpy_type< long >() { return NPY_INT64; }
template <> inline int get_numpy_type< unsigned long >() { return NPY_UINT64; }
template <> inline int get_numpy_type< float >() { return NPY_FLOAT32; }
template <> inline int get_numpy_type< double >() { return NPY_FLOAT64; }
template <> inline int get_numpy_type< std::complex<float> >() { return NPY_COMPLEX64; }
template <> inline int get_numpy_type< std::complex<double> >() { return NPY_COMPLEX128; }
template <> inline int get_numpy_type< ISMRMRD::AcquisitionHeader>() { return NPY_OBJECT;}
template <> inline int get_numpy_type< ISMRMRD::ImageHeader>() { return NPY_OBJECT; }
template <> inline int get_numpy_type< ISMRMRD::ISMRMRD_WaveformHeader>() { return NPY_OBJECT; }
template <> inline int get_numpy_type< ISMRMRD::Waveform>() { return NPY_OBJECT; }
/* Don't define these for now */
/* template <> inline int get_numpy_type< char* >() { return NPY_STRING; } */
/* template <> inline int get_numpy_type< std::string >() { return NPY_STRING; } */

}

#endif // GADGETRON_PYTHON_NUMPY_WRAPPERS_H
