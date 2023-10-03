#pragma once

#include "solver.h"
#include "resampleOperator.h"

namespace Gadgetron{
  
  template <class ARRAY_TYPE> class registrationData
  {
  public:
    registrationData( ARRAY_TYPE *fixed_image, ARRAY_TYPE *moving_image )
    {
      fixed_image_ = fixed_image;
      moving_image_ = moving_image;
    }

    virtual ~registrationData() {}
  
    inline ARRAY_TYPE* get_fixed_image () { return fixed_image_; }
    inline ARRAY_TYPE* get_moving_image () { return moving_image_; }
  
  protected:
    ARRAY_TYPE *fixed_image_;
    ARRAY_TYPE *moving_image_;
  };

  template <class ARRAY_TYPE> class registrationSolver 
    : public solver<registrationData<ARRAY_TYPE>, ARRAY_TYPE >
  {
  public:

    // Constructor/destructor
    //

    registrationSolver() : solver<registrationData<ARRAY_TYPE>,ARRAY_TYPE>() {}
    virtual ~registrationSolver() {}

    // Set interpolator for resampling
    //
  
    inline void set_interpolator( boost::shared_ptr< resampleOperator<ARRAY_TYPE,ARRAY_TYPE> > interpolator )
    {
      interpolator_ = interpolator;
    }
  
    // Set zero deformation boundary condition as a stencil image
    //
  
    inline void set_stencil( boost::shared_ptr<ARRAY_TYPE> stencil )
    {
      stencil_ = stencil;
    }
  
    //
    // The solver adds a dimension to ARRAY_TYPE to hold the vector result.
    // I.e. the vector field dimension is the slowest varying.
    //
  
    virtual boost::shared_ptr<ARRAY_TYPE> 
    solve( ARRAY_TYPE *fixed_image, ARRAY_TYPE *moving_image, bool input_normalization_allowed = false ) = 0;
  
    virtual boost::shared_ptr<ARRAY_TYPE> 
    solve( registrationData< ARRAY_TYPE> *rd )
    {
      return solve( rd->get_fixed_image(), rd->get_moving_image() );
    }
  
    // Deform image based on displacement field
    //

    virtual boost::shared_ptr<ARRAY_TYPE> 
    deform( ARRAY_TYPE *image, boost::shared_ptr<ARRAY_TYPE> displacements )
    {
      if( !interpolator_.get() ){
	    throw std::runtime_error("registrationSolver::deform() : interpolator not set");;
      }
    
      boost::shared_ptr<ARRAY_TYPE> out(new ARRAY_TYPE);
      std::vector<size_t> out_dims = displacements->dimensions(); out_dims.pop_back();    
      out->create(out_dims);
    
      interpolator_->set_displacement_field( displacements );
      interpolator_->mult_M( image, out.get() );
      interpolator_->reset();
    
      return out;
    }

    virtual boost::shared_ptr<ARRAY_TYPE> 
    deform_adj( ARRAY_TYPE *image, boost::shared_ptr<ARRAY_TYPE> displacements )
    {
      if( !interpolator_.get() ){
	    throw std::runtime_error("registrationSolver::deform() : interpolator not set");;
      }
    
      boost::shared_ptr<ARRAY_TYPE> out(new ARRAY_TYPE);
      std::vector<size_t> out_dims = displacements->dimensions(); out_dims.pop_back();    
      out->create(out_dims);
    
      interpolator_->set_displacement_field( displacements );
      interpolator_->mult_MH( image, out.get() );
      interpolator_->reset();
    
      return out;
    }
  
    // Deform image based on an invocation of the registration solver
    //
  
    virtual boost::shared_ptr<ARRAY_TYPE> 
    deform( ARRAY_TYPE *fixed_image, ARRAY_TYPE *moving_image )
    {
      boost::shared_ptr<ARRAY_TYPE> displacements = solve( fixed_image, moving_image );
      return deform( moving_image, displacements );
    }
  
  protected:
    boost::shared_ptr< resampleOperator<ARRAY_TYPE,ARRAY_TYPE> > interpolator_;
    boost::shared_ptr<ARRAY_TYPE> stencil_;
  };
}
