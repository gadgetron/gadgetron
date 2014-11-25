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

namespace Gadgetron{

    // TODO the ignore_segment_ flag is a hack for some EPI sequences
    // should be fixed on the converter side.

  class EXPORTGADGETSMRICORE BucketToBufferGadget : 
  public Gadget1<IsmrmrdAcquisitionBucket>
    {
    public:
      GADGET_DECLARE(BucketToBufferGadget);

      virtual ~BucketToBufferGadget();

      int close(unsigned long flags);
      
    protected:
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
    };

  
}
#endif //BUCKETTOBUFFER_H
