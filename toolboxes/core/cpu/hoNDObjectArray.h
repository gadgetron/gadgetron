/** \file   hoNDObjectArray.h
\brief  CPU-based N-dimensional array for object pointers
if delete_data_on_destruct == true, the object will be released; otherwise, only the object array memory is released
\author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{
    template<class T> using hoNDObjectArray = hoNDArray<T>;

}
