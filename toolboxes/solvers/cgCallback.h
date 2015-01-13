/** \file cgCallback.h
    \brief Class to specify the termination criteria for the conjugate gradient solver through a callback mechanism.
*/

#pragma once

#include "real_utilities.h"
#include "cgSolver.h"

namespace Gadgetron{

  template <class ARRAY_TYPE> class cgSolver;

  template <class ARRAY_TYPE> class cgTerminationCallback
  {

  public:

    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;
    cgTerminationCallback() {}
    virtual ~cgTerminationCallback() {}
  
    virtual bool initialize( cgSolver<ARRAY_TYPE> *cg ){cg_ = cg; return true;}
    virtual bool iterate( unsigned int iteration, REAL *tc_metric, bool *tc_terminate ) = 0;

  protected:

    cgSolver<ARRAY_TYPE> *cg_;

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

  template <class ARRAY_TYPE> class relativeResidualTCB
    : public cgTerminationCallback<ARRAY_TYPE>
  {

  protected:
    typedef cgTerminationCallback<ARRAY_TYPE> cgTC;

  public:

    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

    relativeResidualTCB() : cgTerminationCallback<ARRAY_TYPE>() {
      rq_0_ = REAL(0); 
      tc_last_ = get_max<REAL>();
    }
  
    virtual ~relativeResidualTCB() {}
  
    virtual bool initialize( cgSolver<ARRAY_TYPE> *cg )
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
	  GDEBUG_STREAM("Iteration " << iteration << ". rq/rq_0 = " << *tc_metric << std::endl);
	}
	if( (tc_last_-(*tc_metric)) < REAL(0) ){
	  GDEBUG_STREAM("Warning: conjugate gradient residual increase." << std::endl);
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

  template <class ARRAY_TYPE> class residualTCB
    : public cgTerminationCallback<ARRAY_TYPE>
  {

  protected:

    typedef cgTerminationCallback<ARRAY_TYPE> cgTC;

  public:

    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

    residualTCB() : cgTerminationCallback<ARRAY_TYPE>() {
      tc_last_ = get_max<REAL>();
    }

    virtual ~residualTCB() {}

    virtual bool initialize( cgSolver<ARRAY_TYPE> *cg )
    {
      cgTC::initialize(cg);
      tc_last_ = get_max<REAL>();
      return true;
    }

    virtual bool iterate( unsigned int iteration, REAL *tc_metric, bool *tc_terminate )
    {
      *tc_metric = cgTC::get_rq();
      if( cgTC::cg_->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_WARNINGS ) {
        if( cgTC::cg_->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	  GDEBUG_STREAM("Iteration " << iteration << ". rq/rq_0 = " << *tc_metric << std::endl);
        }
        if( (tc_last_-(*tc_metric)) < REAL(0) ){
	  GDEBUG_STREAM("----- Warning: CG residual increase. Stability problem! -----" << std::endl);
        }
      }
      *tc_terminate = ( *tc_metric < cgTC::cg_->get_tc_tolerance() );
      tc_last_ = *tc_metric;
      return true;
    }

  protected:

    REAL tc_last_;
  };

  template <class ARRAY_TYPE> class updateTCB
    : public cgTerminationCallback<ARRAY_TYPE>
  {

  protected:
    typedef cgTerminationCallback<ARRAY_TYPE> cgTC;

  public:
    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;

    updateTCB() : cgTerminationCallback<ARRAY_TYPE>() {

      tc_last_ = get_max<REAL>();
    }

    virtual ~updateTCB() {}

    virtual bool initialize( cgSolver<ARRAY_TYPE> *cg )
    {
      cgTC::initialize(cg);
      tc_last_ = get_max<REAL>();
      return true;
    }

    virtual bool iterate( unsigned int iteration, REAL *tc_metric, bool *tc_terminate )
    {
      *tc_metric = cgTC::cg_->solver_dot(cgTC::get_p().get(),cgTC::get_p().get());
      if( cgTC::cg_->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_WARNINGS ) {
	if( cgTC::cg_->get_output_mode() >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
	  GDEBUG_STREAM("Iteration " << iteration << ". rq/rq_0 = " << *tc_metric << std::endl);
	}
	if( (tc_last_-(*tc_metric)) < REAL(0) ){
	  GDEBUG_STREAM("----- Warning: CG residual increase. Stability problem! -----" << std::endl);
	}
      }
      *tc_terminate = ( *tc_metric < cgTC::cg_->get_tc_tolerance() );
      tc_last_ = *tc_metric;
      return true;
    }

  protected:

    REAL tc_last_;
  };
}
