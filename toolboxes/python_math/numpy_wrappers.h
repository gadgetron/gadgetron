#ifndef GADGETRON_PYTHON_MATH_NUMPY_WRAPPERS_H
#define GADGETRON_PYTHON_MATH_NUMPY_WRAPPERS_H

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/ndarraytypes.h"

namespace Gadgetron {

/// return the enumerated numpy type for a given C++ type
template <typename T> int get_numpy_type();
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

/// Wrappers for NumPy C-API functions. These functions must be used
/// in the same C++ source file as the call to `import_array()`. In this
/// case, that is PythonMath.cpp. The simplest solution is to lightly wrap the API.
int NumPyArray_NDIM(PyObject* obj);
long NumPyArray_DIM(PyObject* obj, int i);
void* NumPyArray_DATA(PyObject* obj);
int NumPyArray_ITEMSIZE(PyObject* obj);
PyObject* NumPyArray_SimpleNew(int nd, long* dims, int typenum);

}

#endif // GADGETRON_PYTHON_MATH_NUMPY_WRAPPERS_H
