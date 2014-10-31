/** \file   GtPlusGadgetOpenMP.h
    \brief  Pack up the OpenMP support in the GtPlus
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "GadgetronTimer.h"
#include "gtPlusISMRMRDReconUtil.h"

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

namespace Gadgetron
{

bool EXPORTGTPLUSGADGET prepOpenMP();

}
