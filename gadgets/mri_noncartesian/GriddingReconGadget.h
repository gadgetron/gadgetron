#pragma once

#include "gadgetron_mri_noncartesian_export.h"
#include "GriddingReconGadgetBase.h"
#include "cuNDArray.h"
#include "cuNDArray_math.h"

namespace Gadgetron {
    EXPORTGADGETSMRINONCARTESIAN class GriddingReconGadget : public GriddingReconGadgetBase<cuNDArray> {
    public:
        GADGET_DECLARE(GriddingReconGadget);
        GriddingReconGadget();
        ~GriddingReconGadget();


    };


}


