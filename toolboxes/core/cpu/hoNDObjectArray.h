/** \file   hoNDObjectArray.h
\brief  N-dimensional array for objects
        The stored objects should support read/write interfaces
\author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{
    template<class T> using hoNDObjectArray = hoNDArray<T>;
}
