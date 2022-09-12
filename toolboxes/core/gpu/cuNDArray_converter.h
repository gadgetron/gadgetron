#pragma once

#include "hoNDArray.h"
#include "cuNDArray.h"
namespace Gadgetron {
    template<class T>
    boost::shared_ptr <hoNDArray<T>> as_hoNDArray(boost::shared_ptr <cuNDArray<T>> array) {
        return array->to_host();
    }
}
