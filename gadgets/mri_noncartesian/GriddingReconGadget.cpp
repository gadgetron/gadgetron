//
// Created by dchansen on 10/25/18.
//

#include "GriddingReconGadget.h"
#include "cuNFFT.h"
#include "cuNDArray_converter.h"
    #include "b1_map.h"
#include "GriddingReconGadgetBase.hpp"

namespace Gadgetron {
    GriddingReconGadget::GriddingReconGadget() : GriddingReconGadgetBase() {

    }

    GriddingReconGadget::~GriddingReconGadget() {

    }

    GADGET_FACTORY_DECLARE(GriddingReconGadget);
}