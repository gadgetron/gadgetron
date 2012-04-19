#pragma once

#include "linearOperator.h"
#include <iostream>
#include <vector>
#include <boost/smart_ptr.hpp>

template <class REAL, class ARRAY_TYPE> class encodingOperatorContainer
  : public linearOperator<REAL, ARRAY_TYPE>
{
public:
  encodingOperatorContainer() : linearOperator<REAL,ARRAY_TYPE>() { num_elements_ = 0; }
  virtual ~encodingOperatorContainer(){}

  // Get domain and codomain dimensions. 
  // The codomain is a concatenation of the indivudial operators' domains.
  //
  using linearOperator<REAL, ARRAY_TYPE>::get_domain_dimensions;
  virtual boost::shared_ptr< std::vector<unsigned int> > get_codomain_dimensions() 
  { 
    std::vector<unsigned int> *dims = new std::vector<unsigned int>();
    dims->push_back( num_elements_ );
    return boost::shared_ptr< std::vector<unsigned int> >(dims);
  }

  // Get domain and codomain for the individual operators
  //
  virtual boost::shared_ptr< std::vector<unsigned int> > get_domain_dimensions(int i) { return operators_[i]->get_domain_dimensions(); }
  virtual boost::shared_ptr< std::vector<unsigned int> > get_codomain_dimensions(int i) { return operators_[i]->get_codomain_dimensions(); }
  
  // Allocate an array of the codomain dimensions
  //
  boost::shared_ptr< ARRAY_TYPE> create_codomain() 
  {
    std::vector<unsigned int> dims = get_codomain_dimensions();
    ARRAY_TYPE* codomain = new ARRAY_TYPE;
    codomain->create(&dims);
    return boost::shared_ptr<ARRAY_TYPE>(codomain);
  }

  // Get individual operators
  //
  boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE > > get_operator(int i){
    return operators_[i];
  }

  // Get pointer offset into codomain for individual operators "sub-codomains"
  //
  unsigned int get_offset(int i){
    return offsets_[i];
  }

  // Add operator to the container
  //
  bool add_operator(boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE > > op)
  {
    int elements = 1;
    boost::shared_ptr< std::vector<unsigned int> > codomain = op->get_codomain_dimensions();
    for (int i =0; i < codomain->size(); i++){
      elements *= codomain->at(i);
    }
    
    if (offsets_.size() == 0){
      offsets_.push_back(0);
    } else{
      offsets_.push_back(num_elements_);
    }

    num_elements_ += elements;
    operators_.push_back(op);

    return true;
  }
  
  virtual int mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false)
  {
    for (int i =0; i < operators_.size(); i++){

      ARRAY_TYPE tmp;
      if( !tmp.create(operators_[i]->get_codomain_dimensions().get(), out->get_data_ptr()+offsets_[i]) ){
	std::cout << std::endl << "Error: encodingOperatorContainer : failed to create working array" << std::endl;
	return -1;
      }

      if( operators_[i]->mult_M( in, &tmp, accumulate ) < 0 ){
	std::cout << std::endl << "Error: encodingOperatorContainer : mult_M failed on sub-operator" << std::endl;
	return -1;
      }
    }
    return 0;
  }

  virtual int mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false)
  {
    for (int i=0; i < operators_.size(); i++){

      ARRAY_TYPE tmp;
      if( !tmp.create(operators_[i]->get_codomain_dimensions().get(), in->get_data_ptr()+offsets_[i]) ){
	std::cout << std::endl << "Error: encodingOperatorContainer : failed to create working array" << std::endl;
	return -1;
      }
      
      int res;
      if (!accumulate && i == 0){
	res=operators_[i]->mult_MH( &tmp, out, false );
      } else {
	res=operators_[i]->mult_MH( &tmp, out, true );
      }
      if (res < 0){
	std::cout << std::endl << "Error: encodingOperatorContainer : mult_MH failed on sub-operator" << std::endl;
	return -1;
      }
    }
    return 0;
  }
  
  virtual int mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false)
  {
    for (int i=0; i < operators_.size(); i++){
      int res;
      if (!accumulate && i == 0){
	res = operators_[i]->mult_MH_M( in, out, false );
      } else {
	res= operators_[i]->mult_MH_M( in, out, true );
      }
      if (res < 0){
	std::cout << std::endl << "Error: encodingOperatorContainer : mult_MH_M failed on sub-operator" << std::endl;
	return -1;
      }
    }
    return 0;
  }
  
protected:
  std::vector< boost::shared_ptr<linearOperator<REAL, ARRAY_TYPE> > > operators_;
  std::vector<unsigned int> offsets_;
  unsigned int num_elements_;
};
