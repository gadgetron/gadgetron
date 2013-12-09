
#include "gtPlusISMRMRDReconUtil.h"

namespace Gadgetron { namespace gtPlus {

//
// Instantiation
//

template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<float>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<double>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<GT_Complex8>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtil<GT_Complex16>;

template EXPORTGTPLUS class gtPlusISMRMRDReconUtilComplex<GT_Complex8>;
template EXPORTGTPLUS class gtPlusISMRMRDReconUtilComplex<GT_Complex16>;

}}
