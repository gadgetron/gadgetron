#pragma once

#include "hoNDArray.h"

namespace Gadgetron {
    template<class T> boost::shared_ptr<hoNDArray<T>> as_hoNDArray(boost::shared_ptr<hoNDArray<T>> array){
        return array;
    }
}
