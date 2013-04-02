#pragma once

#include "hoNDArray.h"

namespace Gadgetron {

  class ArrayIterator
  {
  public:
    
    ArrayIterator(std::vector<unsigned int> *dimensions, std::vector<unsigned int> *order)
    {
      dimensions_  = boost::shared_ptr< std::vector<unsigned int> >      (new std::vector<unsigned int>);
      order_       = boost::shared_ptr< std::vector<unsigned int> >      (new std::vector<unsigned int>);
      current_     = boost::shared_ptr< std::vector<unsigned int> >      (new std::vector<unsigned int>);
      block_sizes_ = boost::shared_ptr< std::vector<unsigned long int> > (new std::vector<unsigned long int>);
      
      block_sizes_->push_back(1);
      for (unsigned int i = 0; i < order->size(); i++) {
	dimensions_->push_back((*dimensions)[i]);
	order_->push_back((*order)[i]);
	current_->push_back(0);
	if (i > 0) {
	  block_sizes_->push_back((*block_sizes_)[i-1]*(*dimensions_)[i-1]);
	}
      }
      current_idx_ = 0;
    }
    
    inline unsigned long int advance()
    {
      unsigned int order_index = 0;
	(*current_)[(*order_)[order_index]]++;
	while ((*current_)[(*order_)[order_index]] >= (*dimensions_)[(*order_)[order_index]]) {
	  (*current_)[(*order_)[order_index]] = 0;
	  order_index = (order_index+1)%dimensions_->size();
	  (*current_)[(*order_)[order_index]]++;
	}
	
	current_idx_ = 0;
	for (unsigned int i = 0; i < dimensions_->size(); i++) {
	  current_idx_ += (*current_)[i]*(*block_sizes_)[i];
	}	
	return current_idx_;
    }
    
    inline unsigned long int get_current_idx() {
      return current_idx_;
    }
    
    boost::shared_ptr< std::vector<unsigned int> > get_current_sub() {
      return current_;
    }
    
  protected:
    boost::shared_ptr< std::vector<unsigned int> > dimensions_;
    boost::shared_ptr< std::vector<unsigned int> > order_;
    boost::shared_ptr< std::vector<unsigned int> > current_;
    boost::shared_ptr< std::vector<unsigned long int> > block_sizes_;
    unsigned long int current_idx_;
  };
    
  template<class T> boost::shared_ptr< hoNDArray<T> > shift_dim( hoNDArray<T> *in, int shift )  
  {
    if( in == 0x0 ) {
      BOOST_THROW_EXCEPTION(runtime_error("shift_dim(): invalid input pointer provided"));
    }    
    std::vector<unsigned int> order;
    for (unsigned int i = 0; i < in->get_number_of_dimensions(); i++) {
      order.push_back(static_cast<unsigned int>((i+shift)%in->get_number_of_dimensions()));
    }
    return permute(in,&order);
  }

  template<class T> void shift_dim( hoNDArray<T> *in, hoNDArray<T> *out, int shift )
  {
    if( in == 0x0 || out == 0x0 ) {
      BOOST_THROW_EXCEPTION(runtime_error("shift_dim(): invalid pointer provided"));
    }    
    std::vector<unsigned int> order;
    for (unsigned int i = 0; i < in->get_number_of_dimensions(); i++) {
      order.push_back(static_cast<unsigned int>((i+shift)%in->get_number_of_dimensions()));
    }
    permute(in,out,&order);
  }
  
  template<class T> boost::shared_ptr< hoNDArray<T> > 
  permute( hoNDArray<T> *in, std::vector<unsigned int> *dim_order, int shift_mode = 0) 
  {
    if( in == 0x0 || dim_order == 0x0 ) {
      BOOST_THROW_EXCEPTION(runtime_error("permute(): invalid pointer provided"));
    }    

    std::vector<unsigned int> dims;
    for (unsigned int i = 0; i < dim_order->size(); i++)
      dims.push_back(in->get_dimensions()->at(dim_order->at(i)));
    boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>() );    
    out->create(&dims);
    permute( in, out.get(), dim_order, shift_mode );
    return out;
  }
  
  template<class T> void 
  permute( hoNDArray<T> *in, hoNDArray<T> *out, std::vector<unsigned int> *dim_order, int shift_mode = 0) 
  {
    if( in == 0x0 || out == 0x0 || dim_order == 0x0 ) {
      BOOST_THROW_EXCEPTION(runtime_error("permute(): invalid pointer provided"));
    }    
    
    // Check ordering array
    if (dim_order->size() > in->get_number_of_dimensions()) {
      BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute - Invalid length of dimension ordering array"));
    }
    
    std::vector<unsigned int> dim_count(in->get_number_of_dimensions(),0);
    for (unsigned int i = 0; i < dim_order->size(); i++) {
      if ((*dim_order)[i] >= in->get_number_of_dimensions()) {
	BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute - Invalid dimension order array"));      
      }
      dim_count[(*dim_order)[i]]++;
    }
      
    // Create an internal array to store the dimensions
    std::vector<unsigned int> dim_order_int;
  
    // Check that there are no duplicate dimensions
    for (unsigned int i = 0; i < dim_order->size(); i++) {
      if (dim_count[(*dim_order)[i]] != 1) {
	BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute - Invalid dimension order array (duplicates)"));
      
      }
      dim_order_int.push_back((*dim_order)[i]);
    }
  
    for (unsigned int i = 0; i < dim_order_int.size(); i++) {
      if ((*in->get_dimensions())[dim_order_int[i]] != out->get_size(i)) {
	BOOST_THROW_EXCEPTION(runtime_error("permute(): dimensions of output array do not match the input array"));
      }
    }
    
    // Pad dimension order array with dimension not mentioned in order array
    if (dim_order_int.size() < in->get_number_of_dimensions()) {
      for (unsigned int i = 0; i < dim_count.size(); i++) {
	if (dim_count[i] == 0) {
	  dim_order_int.push_back(i);
	}
      }
    }
    
    T* o = out->get_data_ptr();
    
    ArrayIterator it(in->get_dimensions().get(),&dim_order_int);
    for (unsigned long int i = 0; i < in->get_number_of_elements(); i++) {
      o[i] = in->get_data_ptr()[it.get_current_idx()];
      it.advance();
    }
  }

  // Expand array to new dimension
  template<class T> boost::shared_ptr<hoNDArray<T> > expand(hoNDArray<T> *in, unsigned int new_dim_size )
  {
    if( in == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("expand(): illegal input pointer."));
    }
    
    std::vector<unsigned int> dims = *in->get_dimensions(); dims.push_back(new_dim_size);
    boost::shared_ptr< hoNDArray<T> > out( new hoNDArray<T>()); out->create(&dims);    
    const unsigned int number_of_elements_in = in->get_number_of_elements();    

    for( unsigned int idx=0; idx<number_of_elements_in*new_dim_size; idx++ ){
      out[idx] = in[idx%number_of_elements_in];
    }
    return out;
  }
  
  // Sum over dimension
  template<class T> boost::shared_ptr<hoNDArray<T> > sum(hoNDArray<T> *in, unsigned int dim )
  {
    if( in == 0x0 ){
      BOOST_THROW_EXCEPTION(runtime_error("sum(): illegal input pointer."));
    }

    if( !(in->get_number_of_dimensions()>1) ){
      BOOST_THROW_EXCEPTION(runtime_error("sum(): underdimensioned."));
    }
 
    if( dim > in->get_number_of_dimensions()-1 ){
      BOOST_THROW_EXCEPTION(runtime_error( "sum(): dimension out of range."));
    }

    unsigned int number_of_batches = in->get_size(dim);
    unsigned int number_of_elements = in->get_number_of_elements()/number_of_batches;
    std::vector<unsigned int> dims = *in->get_dimensions(); dims.pop_back();

    boost::shared_ptr< hoNDArray<T> > out(new hoNDArray<T>());
    out->create(&dims);
        
    for( unsigned int idx=0; idx<number_of_elements; idx++ ){
      T val(0);
      for( unsigned int j=0; j<number_of_batches; j++ ){
	unsigned int in_idx = j*number_of_elements+idx;
	val += in->get_data_ptr()[in_idx];      
      }
      out->get_data_ptr()[idx] = val;       
    }
    return out;
  } 
}
