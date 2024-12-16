#include "cuSenseBufferCg.h"
#include "vector_td_utilities.h"
#include "cuNDArray_utils.h"
#include "cuNDArray_reductions.h"
#include "cuNDArray_elemwise.h"

namespace Gadgetron{

  template<class REAL, unsigned int D>
  void cuSenseBufferCg<REAL,D>::
  setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W,
	 unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles )
  {
    cuSenseBuffer<REAL,D>::setup( matrix_size, matrix_size_os, W, num_coils, num_cycles, num_sub_cycles );

    D_ = boost::shared_ptr< cuCgPreconditioner<_complext> >( new cuCgPreconditioner<_complext>() );

    cg_.set_encoding_operator( this->E_ );
    cg_.set_preconditioner( D_ );
    cg_.set_max_iterations( 2 );
    cg_.set_tc_tolerance( 1e-6 );
    cg_.set_output_mode( cuCgSolver<_complext>::OUTPUT_SILENT);
  }

  template<class REAL, unsigned int D>
  void cuSenseBufferCg<REAL,D>::preprocess( cuNDArray<_reald> *traj ) {
    this->E_->preprocess(traj);
    std::vector<size_t> dims = traj->get_dimensions();
    dims.push_back(this->num_coils_);
    this->E_->set_codomain_dimensions(dims);
  }

  template<class REAL, unsigned int D>
  boost::shared_ptr< cuNDArray<complext<REAL> > > cuSenseBufferCg<REAL,D>::get_combined_coil_image()
  {
    // Some validity checks
    //

    if( this->csm_.get() == 0x0 ){
      throw std::runtime_error("cuSenseBufferCg::get_combined_coil_image: csm not set");
    }

    if( !this->E_->is_preprocessed() ){
      throw std::runtime_error("cuSenseBufferCg::get_combined_coil_image: preprocessing not performed");
    }

    // Compute (and scale) rhs
    //

    boost::shared_ptr< cuNDArray<_complext> > rhs = cuSenseBuffer<REAL,D>::get_combined_coil_image();

    if( rhs.get() == 0x0 ){
      throw std::runtime_error("cuSenseBufferCg::get_combined_coil_image: failed to compute rhs");
    }

    *rhs *= this->get_normalization_factor();

    // Define preconditioning weights
    //

    boost::shared_ptr< cuNDArray<REAL> > _precon_weights = sum(abs_square(this->csm_.get()).get(), D);
    reciprocal_sqrt_inplace(_precon_weights.get());
    boost::shared_ptr< cuNDArray<_complext> > precon_weights = real_to_complex<_complext>( _precon_weights.get() );
    _precon_weights.reset();
    D_->set_weights( precon_weights );

    // Solve
    //

    return cg_.solve_from_rhs(rhs.get());
  }

  //
  // Instantiations
  //

  template class cuSenseBufferCg<float,2>;
  template class cuSenseBufferCg<float,3>;
  template class cuSenseBufferCg<float,4>;

  template class cuSenseBufferCg<double,2>;
  template class cuSenseBufferCg<double,3>;
  template class cuSenseBufferCg<double,4>;
}
