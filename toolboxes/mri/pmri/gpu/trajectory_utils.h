#pragma once
#include "cuNDArray.h"
#include <boost/shared_ptr.hpp>
#include "vector_td.h"

namespace Gadgetron {
/**
 * Creates a new set of dcw weights, where weights for all points outside the limit are 0. Useful for estimating coil
 * sensitivity maps from k-space centers. Maybe.
 * @param traj
 * @param dcw
 * @param limit
 * @return
 */
boost::shared_ptr< cuNDArray<float> > filter_dcw(cuNDArray<floatd2>* traj, cuNDArray<float>* dcw,float limit);
}
