#include "CplxDumpGadget.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"

namespace Gadgetron{

  CplxDumpGadget::CplxDumpGadget() 
    : Gadget2<ISMRMRD::AcquisitionHeader,hoNDArray< std::complex<float> > >()
    , buffer_(ACE_Message_Queue_Base::DEFAULT_HWM * 10000, ACE_Message_Queue_Base::DEFAULT_LWM * 10000)
  {
  }

  CplxDumpGadget::~CplxDumpGadget() {}

  int CplxDumpGadget::process_config(ACE_Message_Block* mb)
  {
    filename_ = filename.value();
    return GADGET_OK;
  }

  int CplxDumpGadget::close(unsigned long flags) {
    
    GDEBUG("CplxDumpGadget::close...\n");
    GDEBUG("Number of items on Q: %d\n", buffer_.message_count());

    int ret = Gadget::close(flags);
    unsigned int readouts_buffered = buffer_.message_count();

    if( readouts_buffered == 0 )
      return GADGET_OK;
    
    // Get the array size from the dimensions of the first buffer entry
    //

    ACE_Message_Block* mbq;
    if (buffer_.dequeue_head(mbq) < 0) {
      GDEBUG("Message dequeue failed\n");
      return GADGET_FAIL;
    }

    GadgetContainerMessage< hoNDArray< std::complex<float> > > *daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
    
    if (!daq) {
      GDEBUG("Unable to interpret data on message queue\n");
      return GADGET_FAIL;
    }

    hoNDArray< std::complex<float> > *entry = daq->getObjectPtr();
    std::vector<size_t> dims_profile = *entry->get_dimensions();
    std::vector<size_t> dims = dims_profile;
    dims.push_back(readouts_buffered);

    // Allocate array for result
    //

    hoNDArray< std::complex<float> > result( dims );

    // And copy over the first profile
    //

    {
      hoNDArray< std::complex<float> > tmp( dims_profile, result.get_data_ptr() );
      tmp = *entry;
    }

    mbq->release();
    
    // Copy the remaining profiles to the array
    //
    
    for (unsigned int i = 1; i < readouts_buffered; i++) {
      
      if (buffer_.dequeue_head(mbq) < 0) {
        GDEBUG("Message dequeue failed\n");
        return GADGET_FAIL;
      }
      
      daq = AsContainerMessage<hoNDArray< std::complex<float> > >(mbq);
      
      if (!daq) {
        GDEBUG("Unable to interpret data on message queue\n");
        return GADGET_FAIL;
      }
      
      entry = daq->getObjectPtr();
      hoNDArray< std::complex<float> > tmp( &dims_profile, result.get_data_ptr()+i*entry->get_number_of_elements() );
      tmp = *entry;
      mbq->release();
    }      
  
    // Reshape to get the coil dimension as the last
    //
  
    std::vector<size_t> order= {0, 2, 1};
    result = permute( result, order);

    // Write out the result
    //
  
    write_nd_array< std::complex<float> >( &result, filename_.c_str() );
  
    return GADGET_OK;
  }
  
  int CplxDumpGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
          GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
    
    // Noise should have been consumed by the noise adjust, but just in case...
    //
    
    bool is_noise = ISMRMRD::FlagBit(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT).isSet(m1->getObjectPtr()->flags);
    if (is_noise) {
      m1->release();
      return GADGET_OK;
    }
    
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* copy = new GadgetContainerMessage< hoNDArray< std::complex<float> > >;
    *copy->getObjectPtr() = *m2->getObjectPtr();
    
    if (buffer_.enqueue_tail(copy) < 0) {
      GDEBUG("Failed to add profile to buffer\n");
      copy->release();
      return GADGET_FAIL;
    }
    
    if (this->next()->putq(m1) < 0) {
      GDEBUG("Unable to put data on queue\n");
      return GADGET_FAIL;
    }
    
    return GADGET_OK;
  }
  
  GADGET_FACTORY_DECLARE(CplxDumpGadget)
}
