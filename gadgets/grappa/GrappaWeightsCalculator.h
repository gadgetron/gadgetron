#pragma once

#include <ace/Task.h>

#include "GrappaWeights.h"

template <class T> class GrappaWeightsCalculator : public ACE_Task<ACE_MT_SYNCH>
{
  typedef ACE_Task<ACE_MT_SYNCH> inherited;

 public:
  GrappaWeightsCalculator() 
    : inherited()
   {
    ACE_TRACE(( ACE_TEXT("GrappaWeightsCalculator::GrappaWeightsCalculator") ));
  }

  virtual ~GrappaWeightsCalculator() { }

  virtual int init(void)
  {
    ACE_TRACE(( ACE_TEXT("GrappaWeightsCalculator::init") ));
    return 0;
  }

  virtual int open(void* = 0) 
  {
    ACE_TRACE(( ACE_TEXT("GrappaWeightsCalculator::open") ));
    return this->activate( THR_NEW_LWP | THR_JOINABLE, 1 );
  }

  virtual int close(unsigned long flags);
  virtual int svc(void);

  virtual int add_job( hoNDArray< std::complex<T> >* ref_data,
		       std::vector< std::pair<unsigned int, unsigned int> > sampled_region,
		       unsigned int acceleration_factor,
		       GrappaWeights<T>* destination,
		       std::vector<unsigned int> uncombined_channel_weights,
		       bool include_uncombined_channels_in_combined_weights = true);
};
