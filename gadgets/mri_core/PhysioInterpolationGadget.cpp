#include "PhysioInterpolationGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include <numeric>

#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

  PhysioInterpolationGadget::PhysioInterpolationGadget() 
    : phys_time_index_(0)
    , phases_to_reconstruct_(30)
    , buffer_(ACE_Message_Queue_Base::DEFAULT_HWM * 10, ACE_Message_Queue_Base::DEFAULT_LWM * 10)
  {}

  PhysioInterpolationGadget::~PhysioInterpolationGadget() {}

  int PhysioInterpolationGadget::process_config(ACE_Message_Block* mb)
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

    return GADGET_OK;
  }

  int PhysioInterpolationGadget::close(unsigned long flags) {
    
    int ret = Gadget::close(flags);

    GADGET_DEBUG1("PhysioInterpolationGadget::close...\n");

    GADGET_DEBUG2("Number of items on Q: %d\n", buffer_.message_count());

    if (time_stamps_.size() != buffer_.message_count()) {
      GADGET_DEBUG1("Inconsistent number of messages and time stamps\n");
      buffer_.flush();
      return GADGET_FAIL;
    }
    
    float previous = -100.0;
    float sum_int  = 0.0; 
    std::vector<float> intervals;
    float int_count = 0.0;
    std::vector<unsigned int> cycle_starts;
    for (unsigned int i = 0; i < time_stamps_.size(); i++) {
      GADGET_DEBUG2("Time %d, %f\n", i, time_stamps_[i]);
      if (time_stamps_[i] < previous) {
	cycle_starts.push_back(i);
      } else if (i > 0 ) {
	sum_int += time_stamps_[i]-time_stamps_[i-1];
	intervals.push_back(time_stamps_[i]-time_stamps_[i-1]);
	int_count += 1.0;
      }
      previous = time_stamps_[i];
    }

    std::sort(intervals.begin(),intervals.end());

    float mean_interval = sum_int/int_count;
    float median_interval = intervals[(intervals.size()>>1)];

    float average_cycle_length = 0.0;
    std::vector<float> cycle_lengths;
    float count = 0;
    for (unsigned int i = 1; i < cycle_starts.size(); i++) {
      float clength = time_stamps_[cycle_starts[i]-1] + median_interval - time_stamps_[cycle_starts[i]];
      cycle_lengths.push_back(clength);
    }

    std::sort(cycle_lengths.begin(),cycle_lengths.end());
    float mean_cycle_length = std::accumulate(cycle_lengths.begin(), cycle_lengths.end(), 0.0)/cycle_lengths.size();
    float median_cycle_length = cycle_lengths[(cycle_lengths.size()>>1)];

    GADGET_DEBUG2("We have %d full cyles, first one starting at %d\n", cycle_starts.size()-1, cycle_starts[0]);
    GADGET_DEBUG2("Mean/Median frame width %f/%f\n", mean_interval,median_interval);
    GADGET_DEBUG2("Mean/Median cycle_length %f/%f\n", mean_cycle_length,median_cycle_length);

    //Correct the first cycle assuming it is of median length:
    float first_cycle_offset = (median_cycle_length-median_interval)+time_stamps_[cycle_starts[0]]-time_stamps_[cycle_starts[0]-1];
    for (unsigned int i = 0; i < cycle_starts[0]; i++) {
      time_stamps_[i] += first_cycle_offset;
    }

    //Calculate relative time stamps
    unsigned int current_cycle = 0;
    std::vector<float> relative_cycle_time;
    for (unsigned int i = 0; i < time_stamps_.size(); i++) {
      if (i >= cycle_starts[current_cycle]) {
	current_cycle++;
      }

      if (current_cycle == 0 || (current_cycle > cycle_lengths.size())) {
	//First and last cycle, we have to assume the cycle length, we will use the median cycle.
	relative_cycle_time.push_back(time_stamps_[i]/median_cycle_length + current_cycle);
	GADGET_DEBUG2("Corrected time stamps: %d, %f * (%d)\n",i,relative_cycle_time[i],current_cycle);
      } else {
	//Here we actually know the cycle length
      	relative_cycle_time.push_back(time_stamps_[i]/cycle_lengths[current_cycle-1] + current_cycle);
	GADGET_DEBUG2("Corrected time stamps: %d, %f  (%d)\n",i,relative_cycle_time[i],current_cycle);
      }
    }
    
    //Make a temporary list of all the data pointers from the Q
    std::vector< ISMRMRD::ImageHeader* > hptrs;
    std::vector< hoNDArray< std::complex<float> > * > aptrs;
    
    ACE_Message_Queue<ACE_MT_SYNCH>::ITERATOR it(buffer_);
    for (ACE_Message_Block* entry = 0;
	 it.next (entry) != 0;
         it.advance ()) 
      {
	GadgetContainerMessage< ISMRMRD::ImageHeader >* tmpm1 =
	  AsContainerMessage< ISMRMRD::ImageHeader >(entry);

	GadgetContainerMessage< hoNDArray< std::complex<float> > > * tmpm2 = 
	  AsContainerMessage< hoNDArray< std::complex<float> >  >(entry->cont());
	
	if (!tmpm1 || !tmpm2) {
	  GADGET_DEBUG1("Failed to cast data on Q, bailing out\n");
	  buffer_.flush();
	  return GADGET_FAIL;
	}
	hptrs.push_back(tmpm1->getObjectPtr());
	aptrs.push_back(tmpm2->getObjectPtr());	
      }

    //Let's figure out which time points we would like to interpolate on:
    ///TODO: Deal with mode 1 and other future modes, we are only implementing mode 0 at the moment
    float phase_interval = 1.0f/static_cast<float>(phases_to_reconstruct_);
    float max_time = floor(relative_cycle_time[relative_cycle_time.size()-1]);
    std::vector<float> recon_cycle_time;
    for (float t=1.0;t<(max_time-0.001);t+=phase_interval) {
      recon_cycle_time.push_back(t);
    }
    

    //Now we can loop over each pixel and estimate the new frames, but first we have to have somewhere to put the data
    std::vector< GadgetContainerMessage< ISMRMRD::ImageHeader >* > out_heads;
    std::vector< GadgetContainerMessage< hoNDArray< std::complex<float> > > * > out_data;
    
    for (unsigned int i = 0; i < recon_cycle_time.size(); i++) {
      

    }



    //We can get rid of the data now
    buffer_.flush();

    return ret;
  }

  int PhysioInterpolationGadget::
  process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {
        
    GadgetContainerMessage<ISMRMRD::ImageHeader>* m3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>;
    GadgetContainerMessage< hoNDArray< std::complex<float> > >* m4 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >;

    
    (*m3->getObjectPtr()) = (*m1->getObjectPtr());
    (*m4->getObjectPtr()) = (*m2->getObjectPtr());
    m3->cont(m4);

    if (buffer_.enqueue_tail(m3) < 0) {
      GADGET_DEBUG1("Failed to add image to buffer\n");
      m3->release();
      return GADGET_FAIL;
    }

    time_stamps_.push_back(m1->getObjectPtr()->physiology_time_stamp[phys_time_index_]);

    if (this->next()->putq(m1) < 0) {
      GADGET_DEBUG1("Unable to put data on next Gadgets Q\n");
      return GADGET_FAIL;
    }

    return GADGET_OK;
  }

  GADGET_FACTORY_DECLARE(PhysioInterpolationGadget)
}
