#include "FlowPhaseSubtractionGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"

#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

  FlowPhaseSubtractionGadget::FlowPhaseSubtractionGadget() {}

  FlowPhaseSubtractionGadget::~FlowPhaseSubtractionGadget() {}

  int FlowPhaseSubtractionGadget::process_config(ACE_Message_Block* mb)
  {
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    std::vector<long> dims;
    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    sets_ = e_limits.set().present() ? e_limits.set().get().maximum() + 1 : 1;

    if (sets_ > 2) {
      GADGET_DEBUG1("Phase subtraction only implemented for two sets for now\n");
      GADGET_DEBUG2("Number of sets detected: %d, bailing out.\n", sets_);
      return GADGET_FAIL;
    }

    buffer_ = boost::shared_array< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>[sets_]); 

    size_t bsize = sizeof(GadgetContainerMessage< GadgetContainerMessage<ISMRMRD::ImageHeader> >)*10000;

    for( unsigned int i=0; i<sets_; i++ ){
      buffer_[i].high_water_mark(bsize);
      buffer_[i].low_water_mark(bsize);
    }

    return GADGET_OK;
  }

  int FlowPhaseSubtractionGadget::
  process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {

    // We need two sets to make a phase subtraction
    if (sets_ < 2) {
      if (this->next()->putq(m1) < 0) {
	return GADGET_FAIL;
      }
      return GADGET_OK;
    }

    unsigned int set = m1->getObjectPtr()->set;

    // Enqueue until we have images from both sets
    //

    if( buffer_[set].enqueue_tail(m1) < 0 ){
      GADGET_DEBUG1("Message enqueue failed\n");
      return GADGET_FAIL;
    };

    // Phase subtract 
    //

    while( buffer_[0].message_count()>0 && buffer_[1].message_count()>0 ) {

      ACE_Message_Block *mbq1, *mbq2;

      if( buffer_[0].dequeue_head(mbq1) < 0 || buffer_[1].dequeue_head(mbq2) < 0 ) {
	GADGET_DEBUG1("Message dequeue failed\n");
	if( buffer_[set].message_count() > 0 ) 
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	return GADGET_FAIL;
      }
	
      GadgetContainerMessage<ISMRMRD::ImageHeader> *pm1 = 
	AsContainerMessage<ISMRMRD::ImageHeader>(mbq1);

      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cpm1 = 
	AsContainerMessage<hoNDArray< std::complex<float> > >(mbq1->cont());

      GadgetContainerMessage<ISMRMRD::ImageHeader> *pm2 = 
	AsContainerMessage<ISMRMRD::ImageHeader>(mbq2);

      GadgetContainerMessage< hoNDArray< std::complex<float> > > *cpm2 = 
	AsContainerMessage<hoNDArray< std::complex<float> > >(mbq2->cont());
	
      // Some validity checks
      //

      if( pm1->getObjectPtr()->image_index != pm2->getObjectPtr()->image_index ) {
	GADGET_DEBUG2("Mismatch in image indices detected (%d, %d). Bailing out.\n", 
		      pm1->getObjectPtr()->image_index, pm2->getObjectPtr()->image_index);
	pm1->release();
	if( buffer_[set].message_count() > 0 ){
	  pm2->release();		
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	}
	return GADGET_FAIL;
      }
      
      if (cpm1->getObjectPtr()->get_number_of_elements() != cpm2->getObjectPtr()->get_number_of_elements()) {
	GADGET_DEBUG1("Mismatch in number of elements detected. Bailing out.\n");
	pm1->release();
	if( buffer_[set].message_count() > 0 ){
	  pm2->release();
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	}
	return GADGET_FAIL;
      }

      std::complex<float> *p1 = cpm1->getObjectPtr()->get_data_ptr();
      std::complex<float> *p2 = cpm2->getObjectPtr()->get_data_ptr();

#ifdef USE_OMP
#pragma omp parallel for
#endif
      for( long i = 0; i < (long)m2->getObjectPtr()->get_number_of_elements(); i++ ) {
	std::complex<float> tmp = std::polar((std::abs(p1[i])+std::abs(p2[i]))/2.0f, std::arg(p2[i])-std::arg(p1[i]));
	p2[i] = tmp;
      }
      
      pm1->release();	
      pm2->getObjectPtr()->set = 0;

      if (this->next()->putq(pm2) < 0) {
	if( buffer_[set].message_count() > 0 ) {
	  pm2->release();
	  buffer_[set].dequeue_tail(mbq1); // or m1 will be attempted deleted twice
	}
	return GADGET_FAIL;
      }
    }
    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(FlowPhaseSubtractionGadget)
}
