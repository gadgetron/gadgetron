#include "cuSpiritBuffer.h"
#include "cuCgSolver.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_reductions.h"
#include "hoNDArray_fileio.h"

namespace Gadgetron {

  template<class REAL, unsigned int D>
  void cuSpiritBuffer<REAL,D>::
  setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W, 
	 unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles )
  {      
    cuBuffer<REAL,D>::setup( matrix_size, matrix_size_os, W, num_coils, num_cycles, num_sub_cycles );
    
    E_->setup( this->matrix_size_, this->matrix_size_os_, W );

    cg_.set_encoding_operator( this->E_ );
    cg_.set_max_iterations( 5 );
    cg_.set_tc_tolerance( 1e-8 );
    cg_.set_output_mode( cuCgSolver<_complext>::OUTPUT_VERBOSE);
  }

  template<class REAL, unsigned int D>
  void cuSpiritBuffer<REAL,D>::preprocess( cuNDArray<_reald> *traj ) {
    E_->preprocess(*traj);
    std::vector<size_t> dims = *traj->get_dimensions();
    dims.push_back(this->num_coils_);
    E_->set_codomain_dimensions(&dims);
  }

  template<class REAL, unsigned int D>
  boost::shared_ptr< cuNDArray<complext<REAL> > > cuSpiritBuffer<REAL,D>::get_accumulated_coil_images()
  {
    // Apply adjoint operator to get the rhs
    //

    boost::shared_ptr< cuNDArray<_complext> > rhs = cuBuffer<REAL,D>::get_accumulated_coil_images();

    // Invert by cg solver
    //

    *rhs *= this->get_normalization_factor();
    this->acc_image_ = cg_.solve_from_rhs(rhs.get());

    static int counter = 0;
    char filename[256];
    snprintf((char*)filename, 256, "_coil_images_%d.real", counter);
    write_nd_array<REAL>( abs(this->acc_image_.get())->to_host().get(), filename );
    counter++;

    return this->acc_image_;
  }

  template<class REAL, unsigned int D>
  boost::shared_ptr< cuNDArray<complext<REAL> > > cuSpiritBuffer<REAL,D>::get_combined_coil_image()
  {
    // Get the individual coil images
    //

    if( this->acc_image_.get() == 0x0 ){
      if( this->get_accumulated_coil_images().get() == 0x0 ){ // This updates acc_image_
        throw std::runtime_error("cuSpiritBuffer::get_combined_coil_image: unable to acquire accumulated coil images");
      }
    }
    
    // Compute RSS
    //

    return real_to_complex< complext<REAL> >(sqrt(sum(abs_square(this->acc_image_.get()).get(), 2).get()).get());
  }
  
  //
  // Instantiations
  //

  template class EXPORTGPUPMRI cuSpiritBuffer<float,2>;
  template class EXPORTGPUPMRI cuSpiritBuffer<float,3>;
  template class EXPORTGPUPMRI cuSpiritBuffer<float,4>;

  template class EXPORTGPUPMRI cuSpiritBuffer<double,2>;
  template class EXPORTGPUPMRI cuSpiritBuffer<double,3>;
  template class EXPORTGPUPMRI cuSpiritBuffer<double,4>;
}
