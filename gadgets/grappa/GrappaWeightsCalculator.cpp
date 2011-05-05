#include "GrappaWeightsCalculator.h"

#include "GadgetContainerMessage.h"

template <class T> class GrappaWeightsDescription
{

public:
  std::vector< std::pair<unsigned int, unsigned int> > sampled_region;
  unsigned int acceleration_factor;
  GrappaWeights<T>* destination;
  std::vector<unsigned int> uncombined_channel_weights;
  bool include_uncombined_channels_in_combined_weights;
};

template <class T> int GrappaWeightsCalculator<T>::svc(void) 
{
   ACE_TRACE(( ACE_TEXT("GrappaWeightsCalculator::svc") ));

   ACE_Message_Block *mb;
    
   while (this->getq(mb) >= 0) {   
     if (mb->msg_type() == ACE_Message_Block::MB_HANGUP) {
       if (this->putq(mb) == -1) {
	 ACE_ERROR_RETURN( (LM_ERROR,
			    ACE_TEXT("%p\n"),
			    ACE_TEXT("GrappaWeightsCalculator::svc, putq")),
			   -1);
       }
       break;
     }

     GadgetContainerMessage< GrappaWeightsDescription<T> >* mb1 
       = AsContainerMessage< GrappaWeightsDescription<T> >(mb);

     if (!mb1) {
       mb->release();
       return -2;
     }

     GadgetContainerMessage< hoNDArray< std::complex<T> > >* mb2 
       = AsContainerMessage< hoNDArray< std::complex<T> > >(mb1->cont());

     if (!mb2) {
       mb->release();
       return -3;
     }

     //TODO: call htgrappa weights GPU function


     mb->release();
   }

   return 0;
}

template <class T> int GrappaWeightsCalculator<T>::close(unsigned long flags)
{
  ACE_TRACE(( ACE_TEXT("GrappaWeightsCalculator::close") ));
  
  int rval = 0;
  if (flags == 1) {
    ACE_Message_Block *hangup = new ACE_Message_Block();
    hangup->msg_type( ACE_Message_Block::MB_HANGUP );
    if (this->putq(hangup) == -1) {
	hangup->release();
	ACE_ERROR_RETURN( (LM_ERROR,
			   ACE_TEXT("%p\n"),
			   ACE_TEXT("GrappaWeightsCalculator::close, putq")),
			  -1);
    }
    rval = this->wait();
  }
  return rval;
}


template <class T> int GrappaWeightsCalculator<T>::
add_job( hoNDArray< std::complex<T> >* ref_data,
	 std::vector< std::pair<unsigned int, unsigned int> > sampled_region,
	 unsigned int acceleration_factor,
	 GrappaWeights<T>* destination,
	 std::vector<unsigned int> uncombined_channel_weights,
	 bool include_uncombined_channels_in_combined_weights)
{
   
  GadgetContainerMessage< GrappaWeightsDescription<T> >* mb1 = 
    new GadgetContainerMessage< GrappaWeightsDescription<T> >();

  if (!mb1) {
    return -1;
  }

  mb1->getObjectPtr()->sampled_region = sampled_region;
  mb1->getObjectPtr()->acceleration_factor = acceleration_factor;
  mb1->getObjectPtr()->destination = destination;
  mb1->getObjectPtr()->uncombined_channel_weights = uncombined_channel_weights;
  mb1->getObjectPtr()->include_uncombined_channels_in_combined_weights = 
    include_uncombined_channels_in_combined_weights;


  GadgetContainerMessage< hoNDArray< std::complex<T> > >* mb2 = 
    new GadgetContainerMessage< hoNDArray< std::complex<T> > >();

  if (!mb2) {
    mb1->release();
    return -2;
  }

  mb1->cont(mb2);

  if (!mb2->getObjectPtr()->create(ref_data->get_dimensions())) {
    mb1->release();
    return -3;
  }

  memcpy(mb2->getObjectPtr()->get_data_ptr(), ref_data->get_data_ptr(),
	 ref_data->get_number_of_elements()*sizeof(T)*2);
  
  this->putq(mb1);

  return 0;
}

template class GrappaWeightsCalculator<float>;

