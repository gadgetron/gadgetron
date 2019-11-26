#ifndef MRI_CORE_ACQUISITION_BUCKET_H
#define MRI_CORE_ACQUISITION_BUCKET_H

#include "GadgetContainerMessage.h"
#include "mri_core_data.h"

namespace Gadgetron 
{

  /** 
      This class functions as a storage unit for statistics related to
      the @IsmrmrdAcquisitionData objects.

   */
  struct AcquisitionBucketStats {

      // Set of labels found in the data or ref part of a bucket
      //11D, fixed order [RO, E1, E2, CHA, SLC, PHS, CON, REP, SET, SEG, AVE]
      std::set<uint16_t> kspace_encode_step_1;
      std::set<uint16_t> kspace_encode_step_2;
      std::set<uint16_t> slice;
      std::set<uint16_t> phase;
      std::set<uint16_t> contrast;
      std::set<uint16_t> repetition;
      std::set<uint16_t> set;
      std::set<uint16_t> segment;
      std::set<uint16_t> average;
  };


  /**

     This class serves as the storage unit for buffered data. 
     The @IsmrmrdAcquisitionData structure contains pointers 
     to the GadgetContainerMessages with the data. 

     Data stored in these buckets will automatically get deleted when the object is
     destroyed. 

   */ 
  struct AcquisitionBucket {
    std::vector< Core::Acquisition > data_;
    std::vector< Core::Acquisition > ref_;
    std::vector<AcquisitionBucketStats> datastats_;
    std::vector<AcquisitionBucketStats> refstats_;
    std::vector< Core::Waveform > waveform_;
  };
  
}
#endif //MRI_CORE_ACQUISITION_BUCKET_H
