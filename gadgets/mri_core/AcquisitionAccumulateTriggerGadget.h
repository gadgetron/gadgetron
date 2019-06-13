#ifndef ACQUISITIONACCUMULATETRIGGERGADGET_H
#define ACQUISITIONACCUMULATETRIGGERGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <map>
#include "mri_core_acquisition_bucket.h"

namespace Gadgetron{


  class EXPORTGADGETSMRICORE AcquisitionAccumulateTriggerGadget : 
  public Gadgetron::Gadget1Of2<ISMRMRD::AcquisitionHeader, ISMRMRD::ISMRMRD_WaveformHeader >
    {
    public:
      GADGET_DECLARE(AcquisitionAccumulateTriggerGadget);

      typedef std::map< unsigned short int, GadgetContainerMessage<IsmrmrdAcquisitionBucket>* > map_type_;

      virtual ~AcquisitionAccumulateTriggerGadget();

      int close(unsigned long flags);


    protected:
      GADGET_PROPERTY_LIMITS(trigger_dimension, std::string, "Dimension to trigger on", "",
			     GadgetPropertyLimitsEnumeration, 
			     "kspace_encode_step_1",
			     "kspace_encode_step_2",
			     "average",
			     "slice",
			     "contrast",
			     "phase",
			     "repetition",
			     "set",
			     "segment",
			     "user_0",
			     "user_1",
			     "user_2",
			     "user_3",
			     "user_4",
			     "user_5",
			     "user_6",
			     "user_7",
                 	     "n_acquisitions",
			     "");

      GADGET_PROPERTY_LIMITS(sorting_dimension, std::string, "Dimension to sort by", "", 
			     GadgetPropertyLimitsEnumeration, 
			     "kspace_encode_step_1",
			     "kspace_encode_step_2",
			     "average",
			     "slice",
			     "contrast",
			     "phase",
			     "repetition",
			     "set",
			     "segment",
			     "user_0",
			     "user_1",
			     "user_2",
			     "user_3",
			     "user_4",
			     "user_5",
			     "user_6",
			     "user_7",
			     "n_acquisitions",
			     "");
      
      GADGET_PROPERTY(n_acquisitions_before_trigger, unsigned long, "Number of acquisition before first trigger", 40);
      GADGET_PROPERTY(n_acquisitions_before_ongoing_trigger, unsigned long, "Number of acquisition before ongoing triggers", 40);
      
      IsmrmrdCONDITION trigger_;
      IsmrmrdCONDITION sort_;
      map_type_  buckets_;
      IsmrmrdAcquisitionData prev_;
      unsigned long trigger_events_;
      
      unsigned long n_acq_since_trigger_;
      unsigned long n_acquisitions_before_trigger_;
      unsigned long n_acquisitions_before_ongoing_trigger_;
      std::vector<ISMRMRD::Waveform> wav_buf_;

      virtual int process_config(ACE_Message_Block* mb);

      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1);
      virtual int process(GadgetContainerMessage<ISMRMRD::ISMRMRD_WaveformHeader>* m1);

      virtual int trigger();

    };

  
}
#endif //ACQUISITIONACCUMULATETRIGGERGADGET_H
