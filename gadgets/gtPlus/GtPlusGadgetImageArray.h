/** \file   GtPlusGadgetImageArray.h
    \brief  The GtPlusGadgetImageArray is used by the triggering gadget to store the ISMRMRD ImageHeader information
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/ismrmrd.h"
#include "GadgetIsmrmrdReadWrite.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"

// GtPlusGadgetImageArray stores the ISMRMRD image header info for every
// 2D kspace

namespace Gadgetron
{

struct  EXPORTGTPLUSGADGET GtPlusGadgetImageExt : public ISMRMRD::ImageHeader
{
    // fields added to store the time_stamp and pmu_time_stamp for every incoming read-out line
    // if one line is not acquried, the corresponding time is -1
    std::vector<int>     time_stamps;
    std::vector<int>     pmu_time_stamps;

    GtPlusGadgetImageExt();
    ~GtPlusGadgetImageExt();

    void copy(GtPlusGadgetImageExt& aMessageImage);
    void set_matrix_size(size_t index, size_t size);

    // interpolation is performed
    // this = weight * this + (1-weight)*aMessageImage
    void recomputeHeader(const GtPlusGadgetImageExt& aMessageImage, double weight);
    void dump();
}; 

// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg Ave]
//  0  1  2   3     4  5   6     7   8   9   10
// store a scan with 11 dimensions
#define GT_DIM_NUM 11

struct  EXPORTGTPLUSGADGET GtPlusGadgetImageArray
{
    // size of the image array
    size_t matrix_size[GT_DIM_NUM];

    size_t max_num_of_images_;

    // message information for every 2D image [RO E1 Cha Slice E2 Contrast Phase Rep Set Seg]
    GtPlusGadgetImageExt* imageArray_;

    GtPlusGadgetImageArray();
    GtPlusGadgetImageArray(const GtPlusGadgetImageArray& imArray);
    GtPlusGadgetImageArray(size_t aSize[GT_DIM_NUM]);
    ~GtPlusGadgetImageArray();

    void findDimIndex(Gadgetron::ISMRMRDDIM& dim, size_t& ind);
    bool getSubImageArray(size_t* startInd, size_t* endInd, GtPlusGadgetImageArray& imageArray);
    void resize(size_t aSize[GT_DIM_NUM]);
    bool copy(const GtPlusGadgetImageArray& imageArray);
    size_t get_offset(size_t slc, size_t e2, size_t con, size_t phs, size_t rep, size_t set, size_t seg, size_t ave);
    bool extractGadgetImageArrayEqual(Gadgetron::ISMRMRDDIM& dim, size_t value, GtPlusGadgetImageArray& imageArray);
    bool extractGadgetImageArrayEqual(Gadgetron::ISMRMRDDIM& dim1, size_t value1, Gadgetron::ISMRMRDDIM& dim2, size_t value2, GtPlusGadgetImageArray& imageArray);
    bool extractGadgetImageArrayLessEqual(Gadgetron::ISMRMRDDIM& dim, size_t value, GtPlusGadgetImageArray& imageArray);
    bool extractGadgetImageArray_Dim1LessEqual_Dim2Equal(Gadgetron::ISMRMRDDIM& dim1, size_t value1, Gadgetron::ISMRMRDDIM& dim2, size_t value2, GtPlusGadgetImageArray& imageArray);

    void dump();
};

}
