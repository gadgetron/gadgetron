#ifndef BUCKETTOBUFFER_H
#define BUCKETTOBUFFER_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include <complex>
#include <map>
#include "mri_core_data.h"
#include "mri_core_acquisition_bucket.h"

namespace Gadgetron{

    // TODO the ignore_segment_ flag is a hack for some EPI sequences
    // should be fixed on the converter side.

    // This gadget fills the IsmrmrdReconData structures with kspace readouts and sets up the sampling limits
    // For the cartesian sampling, the filled kspace ensures its center (N/2) is aligned with the specified center in the encoding limits
    // For the non-cartesian sampling, this "center alignment" constraint is not applied and kspace lines are filled as their E1 and E2 indexes

  class EXPORTGADGETSMRICORE BucketToBufferGadget : 
  public Gadget1<IsmrmrdAcquisitionBucket>
    {
    public:
      GADGET_DECLARE(BucketToBufferGadget);

      BucketToBufferGadget();
      virtual ~BucketToBufferGadget();

      int close(unsigned long flags);
      
    protected:
      GADGET_PROPERTY_LIMITS(N_dimension, std::string, "N-Dimensions", "", 
                 GadgetPropertyLimitsEnumeration,
                 "average",
                 "contrast",
                 "phase",
                 "repetition",
                 "set",
                 "segment",
                 "slice",
                 "");

      GADGET_PROPERTY_LIMITS(S_dimension, std::string, "S-Dimensions", "", 
                 GadgetPropertyLimitsEnumeration,
                 "average",
                 "contrast",
                 "phase",
                 "repetition",
                 "set",
                 "segment",
                 "slice",
                 "");

      GADGET_PROPERTY(split_slices, bool, "Split slices", false);
      GADGET_PROPERTY(ignore_segment, bool, "Ignore segment", false);
      GADGET_PROPERTY(verbose, bool, "Whether to print more information", false);

      IsmrmrdCONDITION N_;
      IsmrmrdCONDITION S_;
      bool split_slices_;
      bool ignore_segment_;
      ISMRMRD::IsmrmrdHeader hdr_;
      
      virtual int process_config(ACE_Message_Block* mb);
      virtual int process(GadgetContainerMessage<IsmrmrdAcquisitionBucket>* m1);
      size_t getKey(ISMRMRD::ISMRMRD_EncodingCounters idx);
      size_t getSlice(ISMRMRD::ISMRMRD_EncodingCounters idx);
      size_t getN(ISMRMRD::ISMRMRD_EncodingCounters idx);
      size_t getS(ISMRMRD::ISMRMRD_EncodingCounters idx);

      IsmrmrdReconBit & getRBit(std::map<size_t, GadgetContainerMessage<IsmrmrdReconData>* > & recon_data_buffers, size_t key, uint16_t espace);
      virtual void allocateDataArrays(IsmrmrdDataBuffered &  dataBuffer, ISMRMRD::AcquisitionHeader & acqhdr, ISMRMRD::Encoding encoding, IsmrmrdAcquisitionBucketStats & stats, bool forref);
      virtual void fillSamplingDescription(SamplingDescription & sampling, ISMRMRD::Encoding & encoding, IsmrmrdAcquisitionBucketStats & stats, ISMRMRD::AcquisitionHeader & acqhdr, bool forref);
      virtual void stuff(std::vector<IsmrmrdAcquisitionData>::iterator it, IsmrmrdDataBuffered & dataBuffer, ISMRMRD::Encoding encoding, bool forref);
    };
}
#endif //BUCKETTOBUFFER_H
