#include "cuSenseBuffer.h"

namespace Gadgetron {

  template<class REAL, unsigned int D, bool ATOMICS>
  void cuSenseBuffer<REAL,D,ATOMICS>
  ::setup( _uint64d matrix_size, _uint64d matrix_size_os, REAL W, 
           unsigned int num_coils, unsigned int num_cycles, unsigned int num_sub_cycles )
  {      
    cuBuffer<REAL,D,ATOMICS>::setup(matrix_size, matrix_size_os, W, num_coils, num_cycles, num_sub_cycles );
    
    if( E_.get() == 0x0 ){   
      std::vector<size_t> dims = to_std_vector(this->matrix_size_);    
      E_ = boost::shared_ptr< cuNonCartesianSenseOperator<REAL,D,ATOMICS> >(new cuNonCartesianSenseOperator<REAL,D,ATOMICS>);      
      E_->set_domain_dimensions(&dims);
      E_->setup( this->matrix_size_, this->matrix_size_os_, W );
    }    
  }
  
  template<class REAL, unsigned int D, bool ATOMICS>
  boost::shared_ptr< cuNDArray<complext<REAL> > > cuSenseBuffer<REAL,D,ATOMICS>::get_combined_coil_image()
  {
    if( this->csm_.get() == 0x0 ){
      throw std::runtime_error("cuSenseBuffer::get_combined_coil_image: csm not set");
    }
    
    if( this->acc_image_.get() == 0x0 ){
      if( this->get_accumulated_coil_images().get() == 0x0 ){ // This updates acc_image_
        throw std::runtime_error("cuSenseBuffer::get_combined_coil_image: unable to acquire accumulated coil images");
      }
    }
    
    std::vector<size_t> dims = to_std_vector(this->matrix_size_);
    boost::shared_ptr< cuNDArray<_complext> > image( new cuNDArray<_complext>(&dims) );

    E_->set_csm(this->csm_);
    E_->mult_csm_conj_sum( this->acc_image_.get(), image.get() );

    return image;
  }
  
  //
  // Instantiations
  //

  template class EXPORTGPUPMRI cuSenseBuffer<float,2,true>;
  template class EXPORTGPUPMRI cuSenseBuffer<float,2,false>;

  template class EXPORTGPUPMRI cuSenseBuffer<float,3,true>;
  template class EXPORTGPUPMRI cuSenseBuffer<float,3,false>;

  template class EXPORTGPUPMRI cuSenseBuffer<float,4,true>;
  template class EXPORTGPUPMRI cuSenseBuffer<float,4,false>;

  template class EXPORTGPUPMRI cuSenseBuffer<double,2,false>;
  template class EXPORTGPUPMRI cuSenseBuffer<double,3,false>;
  template class EXPORTGPUPMRI cuSenseBuffer<double,4,false>;
}
