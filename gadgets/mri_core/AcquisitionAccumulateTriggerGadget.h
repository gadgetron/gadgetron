#ifndef ACQUISITIONACCUMULATETRIGGERGADGET_H
#define ACQUISITIONACCUMULATETRIGGERGADGET_H

#include "Gadget.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd.h>
#include <complex>
#include <map>
#include "mri_core_data.h"

namespace Gadgetron{


  class EXPORTGADGETSMRICORE AcquisitionAccumulateTriggerGadget : 
  public Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >
    {
    public:
      GADGET_DECLARE(AcquisitionAccumulateTriggerGadget);

      typedef std::map< unsigned short int, GadgetContainerMessage<IsmrmrdAcquisitionBucket>* > map_type_;

      enum CONDITION {
	KSPACE_ENCODE_STEP_1,
	KSPACE_ENCODE_STEP_2,
	AVERAGE,
	SLICE,
	CONTRAST,
	PHASE,
	REPETITION,
	SET,
	SEGMENT,
	USER_0,
	USER_1,
	USER_2,
	USER_3,
	USER_4,
	USER_5,
	USER_6,
	USER_7,
	NONE
      };
      
      virtual ~AcquisitionAccumulateTriggerGadget();

      int close(unsigned long flags);


    protected:
      CONDITION trigger_;
      CONDITION sort_;
      map_type_  buckets_;
      IsmrmrdAcquisitionData prev_;
      unsigned long trigger_events_;

      virtual int process_config(ACE_Message_Block* mb);

      virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
			  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2);

      virtual int trigger();

    };

  
}
#endif //ACQUISITIONACCUMULATETRIGGERGADGET_H
