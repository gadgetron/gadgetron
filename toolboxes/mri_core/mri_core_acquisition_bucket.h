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
      void add_stats( const ISMRMRD::AcquisitionHeader& header) {
          average.insert(header.idx.average);
          kspace_encode_step_1.insert(header.idx.kspace_encode_step_1);
          kspace_encode_step_2.insert(header.idx.kspace_encode_step_2);
          slice.insert(header.idx.slice);
          contrast.insert(header.idx.contrast);
          phase.insert(header.idx.phase);
          repetition.insert(header.idx.repetition);
          set.insert(header.idx.set);
          segment.insert(header.idx.segment);
      }

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


      void add_acquisition(Core::Acquisition acq) {
          auto& head  = std::get<ISMRMRD::AcquisitionHeader>(acq);
          auto espace = size_t{head.encoding_space_ref};

          if (ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(head.flags)
              || ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING).isSet(head.flags)) {
              ref_.push_back(acq);
              if (refstats_.size() < (espace + 1)) {
                  refstats_.resize(espace + 1);
              }
              refstats_[espace].add_stats(head);
          }
          if (!(ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PARALLEL_CALIBRATION).isSet(head.flags)
              || ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_PHASECORR_DATA).isSet(head.flags))) {
              if (datastats_.size() < (espace + 1)) {
                  datastats_.resize(espace + 1);
              }
              datastats_[espace].add_stats(head);
              data_.emplace_back(std::move(acq));
          }
      }

  };
  
}
#endif //MRI_CORE_ACQUISITION_BUCKET_H
