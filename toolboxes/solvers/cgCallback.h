#pragma once

#include "real_utilities.h"

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class cgSolver;

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class cgTerminationCallback
{
public:
  cgTerminationCallback() {}
  virtual ~cgTerminationCallback() {}
  
  virtual bool initialize( cgSolver<REAL,ELEMENT_TYPE,ARRAY_TYPE> *cg ) = 0;
  virtual bool iterate( cgSolver<REAL,ELEMENT_TYPE,ARRAY_TYPE> *cg, 
  			unsigned int iteration, REAL *tc_metric, bool *tc_terminate ) = 0;
};

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class relativeResidualTCB 
  : public cgTerminationCallback<REAL,ELEMENT_TYPE,ARRAY_TYPE>
{
public:

  relativeResidualTCB() : cgTerminationCallback<REAL,ELEMENT_TYPE,ARRAY_TYPE>() { 
    rq_0_ = REAL(0); 
    tc_last_ = get_max<REAL>();
  }
  
  virtual ~relativeResidualTCB() {}
  
  virtual bool initialize( cgSolver<REAL,ELEMENT_TYPE,ARRAY_TYPE> *cg )
  {
    tc_last_ = get_max<REAL>();
    rq_0_ = real( cg->solver_dot( cg->get_r().get(), cg->get_p().get() ));
    return true;
  }
  
  virtual bool iterate( cgSolver<REAL,ELEMENT_TYPE,ARRAY_TYPE> *cg, 
			unsigned int iteration, REAL *tc_metric, bool *tc_terminate )
  {
    *tc_metric = cg->get_rq()/rq_0_;
    
    if( cg->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_WARNINGS ) {
      if( cg->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	std::cout << "Iteration " << iteration << ". rq/rq_0 = " << *tc_metric << std::endl;
      }
      if( (tc_last_-(*tc_metric)) < REAL(0) ){
	std::cout << "----- Warning: CG residual increase. Stability problem! -----" << std::endl;
      }
    }
    
    *tc_terminate = ( *tc_metric < cg->get_tc_tolerance() );    
    tc_last_ = *tc_metric;
    return true;
  }
  
protected:
  REAL rq_0_;
  REAL tc_last_;
};
