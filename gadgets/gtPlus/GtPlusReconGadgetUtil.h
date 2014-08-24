/** \file   GtPlusReconGadgetUtil.h
    \brief  Store some utilities functions for reconstruction
    \author Hui Xue
*/

#pragma once

#include <complex>
#include "GtPlusGadgetExport.h"
#include "hoNDArray.h"
#include "GtPlusDefinition.h"
#include "GadgetIsmrmrdReadWrite.h"

#include "ismrmrd.h"
#include "ismrmrd_xml.h"

namespace Gadgetron
{

// [Ro E1 Cha Slice E2 Con Phase Rep Set Seg AVE]
//   0  1  2   3    4   5    6     7  8   9  10

    // find the calibration mode from protocol
    bool EXPORTGTPLUSGADGET findCalibMode(ISMRMRD::IsmrmrdHeader& h, Gadgetron::gtPlus::ISMRMRDCALIBMODE& CalibMode, Gadgetron::gtPlus::ISMRMRDDIM& InterleaveDim, double& acceFactorE1, double& acceFactorE2, bool verbose=false);

    // find the encoding limits from protocol
    bool EXPORTGTPLUSGADGET findEncodingLimits(ISMRMRD::IsmrmrdHeader& h, ISMRMRD::EncodingCounters& meas_max_idx, bool verbose=false);

    // find encoding matrix size and FOV
    void EXPORTGTPLUSGADGET findMatrixSizeEncoding(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_encoding[3]);
    void EXPORTGTPLUSGADGET findFOVEncoding(ISMRMRD::IsmrmrdHeader& h, float field_of_view_encoding[3]);

    // find recon matrix size and FOV
    void EXPORTGTPLUSGADGET findMatrixSizeRecon(ISMRMRD::IsmrmrdHeader& h, size_t matrix_size_recon[3]);
    void EXPORTGTPLUSGADGET findFOVRecon(ISMRMRD::IsmrmrdHeader& h, float field_of_view_recon[3]);

    // find the status of a readout line
    bool EXPORTGTPLUSGADGET checkReadoutStatus(uint64_t flag, int samples, Gadgetron::gtPlus::ISMRMRDCALIBMODE& CalibMode, int roLen, 
                        bool& bIsKSpace, bool& bIsRef, bool& bIsNoise, 
                        bool& bIsPhaseCorr, bool& bIsReflect, bool& bIsOther, 
                        bool& bIsNavigator, bool& bIsRTFeedback, bool& bIsHPFeedback, 
                        bool& bIsDummyScan);

    // estimate the max SEG for a segmented acquisition (number of total segments is segment+1)
    // retro_gated_segment_size : number of readout lines acquired in one segment
    // E1, embedded_ref_lines_E1: number of lines measured along E1 and number of reference lines for embedded mode
    bool EXPORTGTPLUSGADGET estimateMaxSEGForRetroGating(Gadgetron::gtPlus::ISMRMRDCALIBMODE CalibMode, 
                                                      double acceFactorE1, double acceFactorE2, 
                                                      size_t retro_gated_segment_size, 
                                                      uint16_t E1, uint16_t embedded_ref_lines_E1, 
                                                      uint16_t E2, uint16_t embedded_ref_lines_E2, 
                                                      uint16_t& segment, bool verbose=false);


    // get debug folder full path
    void EXPORTGTPLUSGADGET getDebugFolderPath(const std::string& debugFolder, std::string& debugFolderPath, bool verbose=false);

    // create a folder with all permissions for all users
    bool EXPORTGTPLUSGADGET createFolderWithAllPermissions(const std::string& workingdirectory);

}
