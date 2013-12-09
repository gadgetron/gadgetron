/** \file   GtPlusGadgetImageArray.h
    \brief  The GtPlusGadgetImageArray is used by the triggering gadget to store the ISMRMRD ImageHeader information
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd.h"
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
    void set_matrix_size(unsigned long long index, ACE_UINT16 size);

    // interpolation is performed
    // this = weight * this + (1-weight)*aMessageImage
    void recomputeHeader(const GtPlusGadgetImageExt& aMessageImage, double weight);
    void dump();
}; 

// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
//  0  1  2   3     4  5   6     7   8   9
// store a scan with 10 dimensions
struct  EXPORTGTPLUSGADGET GtPlusGadgetImageArray
{
    // size of the image array
    ACE_UINT16 matrix_size[10];

    // message information for every 2D image [RO E1 Cha Slice E2 Contrast Phase Rep Set Seg]
    GtPlusGadgetImageExt* imageArray_;

    GtPlusGadgetImageArray();
    GtPlusGadgetImageArray(const GtPlusGadgetImageArray& imArray);
    GtPlusGadgetImageArray(int aSize[10]);
    ~GtPlusGadgetImageArray();

    void findDimIndex(Gadgetron::gtPlus::ISMRMRDDIM& dim, int& ind);
    bool getSubImageArray(unsigned long long* startInd, unsigned long long* endInd, GtPlusGadgetImageArray& imageArray);
    void resize(int aSize[10]);
    bool copy(const GtPlusGadgetImageArray& imageArray);
    int get_offset(int slc, int e2, int con, int phs, int rep, int set, int seg);
    bool extractGadgetImageArrayEqual(Gadgetron::gtPlus::ISMRMRDDIM& dim, unsigned long long value, GtPlusGadgetImageArray& imageArray);
    bool extractGadgetImageArrayEqual(Gadgetron::gtPlus::ISMRMRDDIM& dim1, unsigned long long value1, Gadgetron::gtPlus::ISMRMRDDIM& dim2, unsigned long long value2, GtPlusGadgetImageArray& imageArray);
    bool extractGadgetImageArrayLessEqual(Gadgetron::gtPlus::ISMRMRDDIM& dim, unsigned long long value, GtPlusGadgetImageArray& imageArray);
    bool extractGadgetImageArray_Dim1LessEqual_Dim2Equal(Gadgetron::gtPlus::ISMRMRDDIM& dim1, unsigned long long value1, Gadgetron::gtPlus::ISMRMRDDIM& dim2, unsigned long long value2, GtPlusGadgetImageArray& imageArray);

    void dump();
};

}
