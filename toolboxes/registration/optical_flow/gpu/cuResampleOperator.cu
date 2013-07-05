#include "cuResampleOperator.h"

#include <thrust/host_vector.h>
#include <thrust/generate.h>
#include <thrust/pair.h>
#include <thrust/sort.h> 
#include <thrust/binary_search.h>
#include <thrust/iterator/counting_iterator.h>

using namespace thrust;

namespace Gadgetron{

  template<class T, unsigned int D> void 
  cuResampleOperator<T,D>::mult_MH_preprocess()
  {
    this->preprocessed_ = false;
  
    // Check if a displacement field has been provided
    //
  
    if( !this->offsets_.get() ){
      throw cuda_error("cuResampleOperator::mult_MH_preprocess(): displacement field not set.");
    }

    // Make a device vector wrap of the displacement field
    //

    std::vector<unsigned int> _dims_disp = *this->offsets_->get_dimensions(); _dims_disp.pop_back(); 
    unsigned int num_elements_disp = D;
    while(!_dims_disp.empty()){
      num_elements_disp *= _dims_disp.back();
      _dims_disp.pop_back();
    }
  
    device_vector<REAL> displacements
      ( device_pointer_cast<REAL>(this->offsets_->get_data_ptr()), 
	device_pointer_cast<REAL>(this->offsets_->get_data_ptr()+num_elements_disp) );
  
    // Make sort keys/values array from the deformation field
    //

    unsigned int num_elements_sort = num_elements_disp/D;
  
    this->lower_bounds_ = device_vector<unsigned int>(num_elements_sort);
    this->upper_bounds_ = device_vector<unsigned int>(num_elements_sort);
  
    this->indices_ = device_vector<unsigned int>(get_num_neighbors()*num_elements_sort);
    this->weights_ = device_vector<REAL>(get_num_neighbors()*num_elements_sort);

    device_vector<unsigned int> sort_keys = device_vector<unsigned int>
      (get_num_neighbors()*num_elements_sort);
  
    // Fill arrays
    //

    write_sort_arrays(sort_keys);
    
    // Make copy of sort_keys before the sort modifies it
    //

    device_vector<unsigned int> sort_keys_copy(sort_keys);
  
    // Sort (twice since we have two value arrays)
    //

    sort_by_key(sort_keys.begin(), sort_keys.end(), this->indices_.begin() );
    sort_by_key(sort_keys_copy.begin(), sort_keys_copy.end(), this->weights_.begin() );
  
    // Find start/end indices (buckets) in the two values arrays
    //
  
    counting_iterator<unsigned int> search_begin(0);
    
    lower_bound( sort_keys.begin(), sort_keys.end(), 
		 search_begin, search_begin + num_elements_sort, this->lower_bounds_.begin() );
  
    upper_bound( sort_keys.begin(), sort_keys.end(), 
		 search_begin, search_begin + num_elements_sort, this->upper_bounds_.begin() );
    
    this->preprocessed_ = true;
  }

  template class EXPORTGPUREG cuResampleOperator<float,1>;
  template class EXPORTGPUREG cuResampleOperator<float_complext,1>;

  template class EXPORTGPUREG cuResampleOperator<float,2>;
  template class EXPORTGPUREG cuResampleOperator<float_complext,2>;

  template class EXPORTGPUREG cuResampleOperator<float,3>;
  template class EXPORTGPUREG cuResampleOperator<float_complext,3>;

  template class EXPORTGPUREG cuResampleOperator<float,4>;
  template class EXPORTGPUREG cuResampleOperator<float_complext,4>;

  template class EXPORTGPUREG cuResampleOperator<double,1>;
  template class EXPORTGPUREG cuResampleOperator<double_complext,1>;

  template class EXPORTGPUREG cuResampleOperator<double,2>;
  template class EXPORTGPUREG cuResampleOperator<double_complext,2>;

  template class EXPORTGPUREG cuResampleOperator<double,3>;
  template class EXPORTGPUREG cuResampleOperator<double_complext,3>;

  template class EXPORTGPUREG cuResampleOperator<double,4>;
  template class EXPORTGPUREG cuResampleOperator<double_complext,4>;
}
