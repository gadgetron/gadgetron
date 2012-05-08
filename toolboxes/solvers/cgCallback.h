#pragma once

#include "real_utilities.h"
#include "cgSolver.h"

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class cgSolver;

template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class cgTerminationCallback
{
public:
  cgTerminationCallback() {}
  virtual ~cgTerminationCallback() {}
  
  virtual bool initialize( cgSolver<REAL,ELEMENT_TYPE,ARRAY_TYPE> *cg ){cg_ = cg; return true;}
  virtual bool iterate( unsigned int iteration, REAL *tc_metric, bool *tc_terminate ) = 0;


protected:
  cgSolver<REAL,ELEMENT_TYPE,ARRAY_TYPE> *cg_;
  REAL get_rq(){
	  return cg_->rq_;
  }
  REAL get_rq0(){
  	  return cg_->rq0_;
  }
  REAL get_alpha(){
    	  return cg_->alpha_;
  }
  boost::shared_ptr<ARRAY_TYPE> get_x(){
	  return cg_->x_;

  }
  boost::shared_ptr<ARRAY_TYPE> get_p(){
  	  return cg_->p_;
  }
  boost::shared_ptr<ARRAY_TYPE> get_r(){
	  return cg_->r_;
  }
};



template <class REAL, class ELEMENT_TYPE, class ARRAY_TYPE> class relativeResidualTCB 
  : public cgTerminationCallback<REAL,ELEMENT_TYPE,ARRAY_TYPE>
{
	protected:
	  typedef cgTerminationCallback<REAL,ELEMENT_TYPE,ARRAY_TYPE> cgTC;
public:

  relativeResidualTCB() : cgTerminationCallback<REAL,ELEMENT_TYPE,ARRAY_TYPE>() { 
    rq_0_ = REAL(0); 
    tc_last_ = get_max<REAL>();
  }
  
  virtual ~relativeResidualTCB() {}
  
  virtual bool initialize( cgSolver<REAL,ELEMENT_TYPE,ARRAY_TYPE> *cg )
  {
	cgTC::initialize(cg);
    tc_last_ = get_max<REAL>();
    rq_0_ = cgTC::get_rq0();
    return true;
  }
  
  virtual bool iterate( unsigned int iteration, REAL *tc_metric, bool *tc_terminate )
  {
    *tc_metric = cgTC::get_rq()/rq_0_;
    
    if( cgTC::cg_->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_WARNINGS ) {
      if( cgTC::cg_->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	std::cout << "Iteration " << iteration << ". rq/rq_0 = " << *tc_metric << std::endl;
      }
      if( (tc_last_-(*tc_metric)) < REAL(0) ){
	std::cout << "----- Warning: CG residual increase. Stability problem! -----" << std::endl;
      }
    }
    
    *tc_terminate = ( *tc_metric < cgTC::cg_->get_tc_tolerance() );
    tc_last_ = *tc_metric;
    return true;
  }
  
protected:
  REAL rq_0_;
  REAL tc_last_;
};
