
virtual void permute(std::vector<unsigned int> *dim_order, NDArray<T>* out = 0, int shift_mode = 0) 
{
  hoNDArray<T>* out_int = 0;
  
  // Check ordering array
  if (dim_order->size() > this->dimensions_->size()) {
    BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute - Invalid length of dimension ordering array"));
  }
  
  std::vector<unsigned int> dim_count(this->dimensions_->size(),0);
  for (unsigned int i = 0; i < dim_order->size(); i++) {
    if ((*dim_order)[i] >= this->dimensions_->size()) {
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
  
  // Pad dimension order array with dimension not mentioned in order array
  if (dim_order_int.size() < this->dimensions_->size()) {
    for (unsigned int i = 0; i < dim_count.size(); i++) {
      if (dim_count[i] == 0) {
	dim_order_int.push_back(i);
      }
    }
  }
  
  if (out) {
    out_int = dynamic_cast< hoNDArray<T>* >(out);
    if (!out_int) {
      BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute: failed to dynamic cast out array pointer"));
      
    }
    for (unsigned int i = 0; i < dim_order_int.size(); i++) {
      if ((*this->dimensions_)[dim_order_int[i]] != out_int->get_size(i)) {
	BOOST_THROW_EXCEPTION(runtime_error("hoNDArray::permute: Dimensions of output array do not match the input array"));
      }
    }
  }
  
  T* o = 0;
  if (out_int) {
    o = out_int->get_data_ptr();
  } else {
    o = new T[this->elements_];
  }
  
  ArrayIterator it(this->dimensions_.get(),&dim_order_int);
  for (unsigned long int i = 0; i < this->elements_; i++) {
    o[i] = this->data_[it.get_current_idx()];
    it.advance();
  }
  
  if (!out_int) {
    boost::shared_ptr<std::vector<unsigned int> > tmp_dims(new std::vector<unsigned int>);
    for (unsigned int i = 0; i < this->dimensions_->size(); i++) {
      tmp_dims->push_back((*this->dimensions_)[dim_order_int[i]]);
    }
    this->dimensions_ = tmp_dims;
    
    delete [] this->data_;
    this->data_ = o;
  }
}



void shift_dim(int shift, NDArray<T> *out = 0) {
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < dimensions_->size(); i++) {
    order.push_back(static_cast<unsigned int>((i+shift)%dimensions_->size()));
  }
  permute(&order,out);
}
