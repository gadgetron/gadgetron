#ifndef GADGETRON_PYTHON_MATH_CONVERSIONS_H
#define GADGETRON_PYTHON_MATH_CONVERSIONS_H

#include "ismrmrd/ismrmrd.h"

namespace Gadgetron {

/// Interface for registering C++ <-> NumPy type converters.
/// A static function on the `python_converter` struct allows
/// for partial template specialization.
template <typename T>
struct python_converter {
    static void create() { }
};

/// Convenience wrapper for `python_converter<TS...>::create()`
template <typename ...TS>
void register_converter(void) {
    // Parameter packs can only be expanded in specific semantic situations.
    // This creates a fake array to expand and create converters for each
    // variadic type.
    using expander = int[];
    (void) expander {0, (python_converter<TS>::create(), 0)...};
}

}

#include "patchlevel.h"
#include "python_tuple_converter.h"
#include "python_hoNDArray_converter.h"
#include "python_ismrmrd_converter.h"
#include "python_vector_converter.h"
#include "python_IsmrmrdReconData_converter.h"
#include "python_IsmrmrdImageArray_converter.h"

#endif // GADGETRON_PYTHON_MATH_CONVERSIONS_H
