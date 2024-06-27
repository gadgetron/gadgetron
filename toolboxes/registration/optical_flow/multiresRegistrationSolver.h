/** \file multiresRegistrationSolver.h
    Abstract class defining a multiresolution registration solver.
    Pure virtual functions are expected to do the actual work.
*/

#pragma once

#include "registrationSolver.h"
#include "vector_td_utilities.h"
#include "vector_td_operators.h"
#include <boost/make_shared.hpp>

namespace Gadgetron{

  template<class ARRAY_TYPE_REAL, unsigned int D> class multiresRegistrationSolver
    : public registrationSolver<ARRAY_TYPE_REAL>
  {
  protected:
    typedef typename ARRAY_TYPE_REAL::element_type REAL;

  public:

    multiresRegistrationSolver() : registrationSolver<ARRAY_TYPE_REAL>(){
      num_multires_levels_ = 0;
      max_num_iterations_per_level_ = 500;
    }

    virtual ~multiresRegistrationSolver() {}

    // Utilities to specify the registration settings
    //

    virtual void set_num_multires_levels( unsigned int levels ) {
      num_multires_levels_ = levels; }

    virtual void set_max_num_iterations_per_level( unsigned int iterations ) {
      max_num_iterations_per_level_ = iterations; }

    //
    // The main solver interface
    //

    virtual boost::shared_ptr<ARRAY_TYPE_REAL> solve( registrationData<ARRAY_TYPE_REAL> *rd )
    {
      return registrationSolver<ARRAY_TYPE_REAL>::solve(rd);
    }
  
    virtual boost::shared_ptr<ARRAY_TYPE_REAL> solve(
                                                     ARRAY_TYPE_REAL *fixed_image,
                                                     ARRAY_TYPE_REAL *moving_image,
                                                     bool input_normalization_allowed = false  )
    {
      // Some initial validity tests
      //

      if( !fixed_image || !moving_image ){
        throw std::runtime_error("multiresRegistrationSolver::solve : invalid input pointer.");
      }

      if( !this->interpolator_.get() ){
        throw std::runtime_error("multiresRegistrationSolver::solve : interpolator not set.");
      }

      typename uint64d<D>::Type fixed_dims = from_std_vector<size_t,D>(moving_image->get_dimensions());
      typename uint64d<D>::Type moving_dims = from_std_vector<size_t,D>(fixed_image->get_dimensions());

      if(!(fixed_dims == moving_dims)){
        throw std::runtime_error("multiresRegistrationSolver::solve : fixed/moving image base dimensions mismatch.");
      }

      if( weak_less_equal(fixed_dims>>num_multires_levels_, vector_td<size_t, D>(1)) ){
        throw std::runtime_error("multiresRegistrationSolver::solve : too many multiresolution levels for image dimensionality.");
      }

      // Normalize the input
      //

      ARRAY_TYPE_REAL *normalized_fixed;
      ARRAY_TYPE_REAL *normalized_moving;

      boost::shared_ptr<ARRAY_TYPE_REAL> garbage_collector_fixed, garbage_collector_moving;
      bool use_padding = padding_required(fixed_dims);

      if( input_normalization_allowed ){
        if( use_padding ){
          throw std::runtime_error("multiresRegistrationSolver::solve : input normalization not possible as image padding is required.");
        }
        else{
          normalized_fixed = fixed_image;
          normalized_moving = moving_image;
        }
      }
      else{
        if( use_padding ){
          garbage_collector_fixed = boost::make_shared<ARRAY_TYPE_REAL>(pad<REAL,D>(round_pow2(fixed_dims), *fixed_image));
          garbage_collector_moving = boost::make_shared<ARRAY_TYPE_REAL>(pad<REAL,D>(round_pow2(moving_dims), *moving_image));
          normalized_fixed = garbage_collector_fixed.get();
          normalized_moving = garbage_collector_moving.get();
        }
        else{
          normalized_fixed = new ARRAY_TYPE_REAL(*fixed_image);
          normalized_moving = new ARRAY_TYPE_REAL(*moving_image);
          garbage_collector_fixed = boost::shared_ptr<ARRAY_TYPE_REAL>(normalized_fixed);
          garbage_collector_moving = boost::shared_ptr<ARRAY_TYPE_REAL>(normalized_moving);
        }
      }

      normalize(normalized_fixed, REAL(1));
      normalize(normalized_moving, REAL(1));

      // Invoke multi-resolution solver
      //

      if( this->output_mode_ >= registrationSolver<ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
        GDEBUG_STREAM(std::endl << "Starting multiresolution registration " <<  std::endl);
      }

      boost::shared_ptr<ARRAY_TYPE_REAL> result =
        solveMultiRes( num_multires_levels_, normalized_fixed, normalized_moving, this->stencil_.get() );

      if( use_padding ){
        result = boost::make_shared<ARRAY_TYPE_REAL>(crop<REAL,D>( (round_pow2(fixed_dims)-fixed_dims)>>2, fixed_dims, *result));
      }

      return result;
    }

  protected:

    // Pure virtual fuctions to be implemented in a subclass
    //

    virtual void compute( ARRAY_TYPE_REAL *fixed_image, ARRAY_TYPE_REAL *moving_image, ARRAY_TYPE_REAL *stencil_image, 
                          boost::shared_ptr<ARRAY_TYPE_REAL> &result ) = 0;

    // The recursive multi-resolution solver
    //

    virtual boost::shared_ptr<ARRAY_TYPE_REAL> solveMultiRes(
                                                             unsigned int res_level,
                                                             ARRAY_TYPE_REAL *fixed_image,
                                                             ARRAY_TYPE_REAL *moving_image,
                                                             ARRAY_TYPE_REAL *stencil_image )
    {
      boost::shared_ptr<ARRAY_TYPE_REAL> result;

      if (res_level>0){

        //
        // We are not yet at the end of the multi-resolution chain
        //

        // Downsample input images (and stencil if provided)
        //

        boost::shared_ptr<ARRAY_TYPE_REAL> fixed_image_lowres = boost::make_shared<ARRAY_TYPE_REAL>(downsample<REAL,D>(*fixed_image));
        boost::shared_ptr<ARRAY_TYPE_REAL> moving_image_lowres = boost::make_shared<ARRAY_TYPE_REAL>(downsample<REAL,D>(*moving_image));
        boost::shared_ptr<ARRAY_TYPE_REAL> stencil_image_lowres =
          ((stencil_image) ? boost::make_shared<ARRAY_TYPE_REAL>(downsample<REAL,D>(*stencil_image)) : boost::shared_ptr<ARRAY_TYPE_REAL>());

        // Compute displacement field at the downsampled resolution
        //

        boost::shared_ptr<ARRAY_TYPE_REAL> result_lowres =
          solveMultiRes( res_level-1, fixed_image_lowres.get(), moving_image_lowres.get(), stencil_image_lowres.get() );

        // Clean up low resolution image data
        //

        fixed_image_lowres.reset();
        moving_image_lowres.reset();
        stencil_image_lowres.reset();

        // Upsample lowres results to current resolution
        //

        result = boost::make_shared<ARRAY_TYPE_REAL>(upsample<REAL,D>(*result_lowres));
        *result *= REAL(2); // To adjust the flow vectors to the fact that the resolution is now twice as high

        // Clean up low resolution result
        //

        result_lowres.reset();

        // Some output to track our progress at runtime
        //

        if( this->output_mode_ >= registrationSolver<ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
          GDEBUG_STREAM(std::endl << "Multiresolution level " << res_level);
        }

        // Use estimated (lowres) motion to compute displacements at the current resolution
        //

        boost::shared_ptr<ARRAY_TYPE_REAL> def_moving_image = this->deform( moving_image, result );
      
        // Compute registationnat the current multiresolution level
        //

        compute( fixed_image, def_moving_image.get(), stencil_image, result );
      }	
      else{

        //
        // We are now at the end of the multi-resolution chain
        //

        // Some output to track our progress at runtime
        //

        if( this->output_mode_ >= registrationSolver<ARRAY_TYPE_REAL>::OUTPUT_VERBOSE ) {
          GDEBUG_STREAM(std::endl << "Multiresolution level " << res_level << " (lowest)");
        }

        // Compute displacements at the current resolution (no estimate can be provided)
        //

        compute( fixed_image, moving_image, stencil_image, result );
      }

      return result;
    }

    virtual bool padding_required( typename uint64d<D>::Type dims )
    {
      bool padding_required = false;
      typename uint64d<D>::Type ones(1);
      typename uint64d<D>::Type twos(2);

      for( unsigned int i=0; i<num_multires_levels_; i++ ){

        dims /= (size_t)2;

        if( weak_less( dims, (size_t)12*ones ) ){
          throw std::runtime_error("multiresRegistrationSolver::padding_required : resolution too low. Too many multiresolution levels specified?");
        }

        if( weak_equal(dims%twos, ones) ){
          padding_required = true;
        }
      }
      return padding_required;
    }

  protected:
    unsigned int num_multires_levels_;
    unsigned int max_num_iterations_per_level_;

  private:
    typename uint64d<D>::Type round_pow2(typename uint64d<D>::Type v)
    {
      typename uint64d<D>::Type ones(1);
      typename uint64d<D>::Type out = v-ones;
      for( unsigned int d=0; d<D; d++ ){
        out[d] |= out[d] >> 1;
        out[d] |= out[d] >> 2;
        out[d] |= out[d] >> 4;
        out[d] |= out[d] >> 8;
        out[d] |= out[d] >> 16;
      }
      return out+ones;
    }
  };
}
