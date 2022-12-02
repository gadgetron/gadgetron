/** \file cgSolver.h
    \brief Base class for the conjugate gradient solvers.

    The file cgSolver.h is a device independent implementation of the conjugate gradient solver.
    To simplify the actual instantiation we refer to 
    - the class(/file) hoCgSolver(/.h) for a cpu instantiated solver using the hoNDArray class
    - the class(/file) cuCgSolver(/.h) for a gpu instantiated solver using the cuNDArray class
    - the class(/file) hoCuCgSolver(/.h) for a gpu based solver using a host memory interface. 

    The latter version is intended for large reconstructions in which device memory cannot hold 
    the entire data from the image and encoded image domains. 
    In the "hoCu" scenario, suitable encoding and regularization operators
    capable of batching their mult_M and mult_MHM functions should be chosen.

    In all cases, the encoding and regularization operators added to the solver 
    must adhere to the underlying instantiation of the NDArray data type.
*/

#pragma once

#include "linearOperatorSolver.h"
#include "cgCallback.h"
#include "cgPreconditioner.h"
#include "real_utilities.h"
#include "complext.h"

#include <vector>
#include <iostream>
#include <limits>

namespace Gadgetron{

  template <class ARRAY_TYPE> class cgSolver : public linearOperatorSolver<ARRAY_TYPE>
  {
  
  public:

    // Convienient typedefs
    //

    typedef typename ARRAY_TYPE::element_type ELEMENT_TYPE;
    typedef typename realType<ELEMENT_TYPE>::Type REAL;
    friend class cgTerminationCallback<ARRAY_TYPE>;


    // Constructor
    //

    cgSolver() : linearOperatorSolver<ARRAY_TYPE>() {
      alpha_ = std::numeric_limits<ELEMENT_TYPE>::quiet_NaN();
      iterations_ = 10;
      tc_tolerance_ = (REAL)1e-3;
      cb_ = boost::shared_ptr< relativeResidualTCB<ARRAY_TYPE> >( new relativeResidualTCB<ARRAY_TYPE>() );
    }
  

    // Destructor
    //

    virtual ~cgSolver() {}


    // Set preconditioner
    //

    virtual void set_preconditioner( boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond ) {
      precond_ = precond;
    }
  

    // Set termination callback
    //

    virtual void set_termination_callback( boost::shared_ptr< cgTerminationCallback<ARRAY_TYPE> > cb ){
      cb_ = cb;
    }
  

    // Set/get maximally allowed number of iterations
    //

    virtual void set_max_iterations( unsigned int iterations ) { iterations_ = iterations; }
    virtual unsigned int get_max_iterations() { return iterations_; }  


    // Set/get tolerance threshold for termination criterium
    //

    virtual void set_tc_tolerance( REAL tolerance ) { tc_tolerance_ = tolerance; }
    virtual REAL get_tc_tolerance() { return tc_tolerance_; }
  

    // Virtual function that is provided with the intermediate solution at each solver iteration.
    // The default behaviour is to do nothing with this array,
    // but customization is possible by specialization of the virtual function in a derived class.
    //

    virtual void solver_dump( ARRAY_TYPE* ) {}


    //
    // Main solver interface
    //

    virtual boost::shared_ptr<ARRAY_TYPE> solve( ARRAY_TYPE *d )
    {
    
      // Compute right hand side...
      //
      
      boost::shared_ptr<ARRAY_TYPE> rhs = compute_rhs( d );

      // ... and the result
      //

      boost::shared_ptr<ARRAY_TYPE> result =  solve_from_rhs( rhs.get() );
      return result;
    }


    // Alternative solver interface when given the right hand side
    //

    virtual boost::shared_ptr<ARRAY_TYPE> solve_from_rhs( ARRAY_TYPE *rhs ) 
    {
      // For zero iterations we have computed / return the right hand side
      //

      if( iterations_ == 0 ){
        return boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE(*rhs) );
      }

      // Initialize
      //

      initialize(rhs);

      // Iterate
      //

      if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ){
        GDEBUG_STREAM("Iterating..." << std::endl);
      }
    
      for( unsigned int it=0; it<iterations_; it++ ){

        REAL tc_metric;
        bool tc_terminate;
      
        this->iterate( it, &tc_metric, &tc_terminate );

        solver_dump( x_.get());
      
        if( tc_terminate )
          break;
      }
    
      // Clean up and we are done
      //

      boost::shared_ptr<ARRAY_TYPE> tmpx = x_;
      deinitialize();
      return tmpx;
    }


    // Compute right hand side
    //

    virtual boost::shared_ptr<ARRAY_TYPE> compute_rhs( ARRAY_TYPE *d )
    {
    
      if( this->encoding_operator_.get() == 0 ){
      	throw std::runtime_error( "Error: cgSolver::compute_rhs : no encoding operator is set" );
      } 
        
      // Get image space dimensions from the encoding operator
      //

      boost::shared_ptr< std::vector<size_t> > image_dims = this->encoding_operator_->get_domain_dimensions();
      if( image_dims->size() == 0 ){
      	throw std::runtime_error( "Error: cgSolver::compute_rhs : encoding operator has not set domain dimension" );
      }

      // Create result array and clear
      //

      boost::shared_ptr<ARRAY_TYPE> result = boost::shared_ptr<ARRAY_TYPE>(new ARRAY_TYPE(*image_dims));
      clear(result.get());
    
      // Create temporary array
      //

      ARRAY_TYPE tmp( *image_dims );

      // Compute operator adjoint
      //

      this->encoding_operator_->mult_MH( d, &tmp );
    
      // Apply weight
      //

      axpy(ELEMENT_TYPE(this->encoding_operator_->get_weight()), &tmp, result.get() );
    
      return result;
    }

  protected:
  
    //
    // Everything beyond this point is internal to the implementation
    // and not intended to be exposed as a public interface
    //

    // Initialize solver
    //

    virtual void initialize( ARRAY_TYPE *rhs )
    {
      // Input validity test
      //

      if( !rhs || rhs->get_number_of_elements() == 0 ){
      	throw std::runtime_error( "Error: cgSolver::initialize : empty or NULL rhs provided" );
      }
    
      // Result, x
      //

      x_ = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE(rhs->dimensions()) );
    
    
      // Initialize r,p,x
      //

      r_ = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE(*rhs) );
      p_ = boost::shared_ptr<ARRAY_TYPE>( new ARRAY_TYPE(*r_) );
    
      if( !this->get_x0().get() ){ // no starting image provided      
	clear(x_.get());
      }

      // Apply preconditioning, twice (should change preconditioners to do this)
      //
      
      if( precond_.get() ) {	
        precond_->apply( p_.get(), p_.get() );
        precond_->apply( p_.get(), p_.get() );
      }

      rq0_ = real(dot( r_.get(), p_.get() ));

      if (this->get_x0().get()){
	
        if( !this->get_x0()->dimensions_equal( rhs )){
          throw std::runtime_error( "Error: cgSolver::initialize : RHS and initial guess must have same dimensions" );
        }
	
        *x_ = *(this->get_x0());
        
        ARRAY_TYPE mhmX(rhs->dimensions());

        if( this->output_mode_ >= solver<ARRAY_TYPE,ARRAY_TYPE>::OUTPUT_VERBOSE ) {
          GDEBUG_STREAM("Preparing guess..." << std::endl);
        }
        
        mult_MH_M( this->get_x0().get(), &mhmX );
        
        *r_ -= mhmX;
        *p_ = *r_;
        
        // Apply preconditioning, twice (should change preconditioners to do this)
        //
        
        if( precond_.get() ){
          precond_->apply( p_.get(), p_.get() );
          precond_->apply( p_.get(), p_.get() );
        }
      }
      
      rq_ = real( dot( r_.get(), p_.get() ));
      
      // Invoke termination callback initialization
      //
    
      cb_->initialize(this);
    }
  
    // Clean up
    //

    virtual void deinitialize()
    {
      p_.reset();
      r_.reset();
      x_.reset();
    }

    // Perform full cg iteration
    //

    virtual void iterate( unsigned int iteration, REAL *tc_metric, bool *tc_terminate )
    {
      ARRAY_TYPE q = ARRAY_TYPE(x_->dimensions());

      // Perform one iteration of the solver
      //

      mult_MH_M( p_.get(), &q );
    
      // Update solution
      //

      alpha_ = rq_/dot( p_.get(), &q );
      axpy( alpha_, p_.get(), x_.get());

      // Update residual
      //

      axpy( -alpha_, &q, r_.get());

      // Apply preconditioning
      //

      if( precond_.get() ){

        precond_->apply( r_.get(), &q );
        precond_->apply( &q, &q );
        
        REAL tmp_rq = real(dot( r_.get(), &q ));      
        *p_ *= ELEMENT_TYPE((tmp_rq/rq_));
        axpy( ELEMENT_TYPE(1), &q, p_.get() );
        rq_ = tmp_rq;
      } 
      else{
        
        REAL tmp_rq = real(dot( r_.get(), r_.get()) );
        *p_ *= ELEMENT_TYPE((tmp_rq/rq_));           
        axpy( ELEMENT_TYPE(1), r_.get(), p_.get() );
        rq_ = tmp_rq;      
      }
      
      // Invoke termination callback iteration
      //

      if( !cb_->iterate( iteration, tc_metric, tc_terminate ) ){
        throw std::runtime_error( "Error: cgSolver::iterate : termination callback iteration failed" );
      }    
    }
    
    // Perform mult_MH_M of the encoding and regularization matrices
    //

    void mult_MH_M( ARRAY_TYPE *in, ARRAY_TYPE *out )
    {
      // Basic validity checks
      //

      if( !in || !out ){
        throw std::runtime_error( "Error: cgSolver::mult_MH_M : invalid input pointer(s)" );
      }

      if( in->get_number_of_elements() != out->get_number_of_elements() ){
        throw std::runtime_error( "Error: cgSolver::mult_MH_M : array dimensionality mismatch" );
      }
    
      // Intermediate storage
      //

      ARRAY_TYPE q = ARRAY_TYPE(in->dimensions());

      // Start by clearing the output
      //
      clear(out);

      // Apply encoding operator
      //

      this->encoding_operator_->mult_MH_M( in, &q, false );
      axpy( ELEMENT_TYPE (this->encoding_operator_->get_weight()), &q, out );

      // Iterate over regularization operators
      //

      for( unsigned int i=0; i<this->regularization_operators_.size(); i++ ){      
        this->regularization_operators_[i]->mult_MH_M( in, &q, false );
        axpy( ELEMENT_TYPE(this->regularization_operators_[i]->get_weight()), &q, out );
      }      
    }
    
  protected:

    // Preconditioner
    boost::shared_ptr< cgPreconditioner<ARRAY_TYPE> > precond_;

    // Termination criterium callback
    boost::shared_ptr< cgTerminationCallback<ARRAY_TYPE> > cb_;

    // Termination criterium threshold
    REAL tc_tolerance_;

    // Maximum number of iterations
    unsigned int iterations_;

    // Internal variables. 
    REAL rq_;
    REAL rq0_;
    ELEMENT_TYPE alpha_;
    boost::shared_ptr<ARRAY_TYPE> x_, p_, r_;
  };
}
