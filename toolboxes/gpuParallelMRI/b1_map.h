#pragma once

#include "vector_td.h"
#include "cuNDArray.h"
#include <boost/smart_ptr.hpp>

//
// Estimate b1 map
//

template<class REAL, unsigned int D> boost::shared_ptr< cuNDArray<typename complext<REAL>::Type> >
estimate_b1_map( cuNDArray<typename complext<REAL>::Type> *data );
