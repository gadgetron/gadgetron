/** \file encodingOperatorContainer.h
    \brief Encoding operator that can contain multiple other encoding operators. Use when more than one encoding operator is required in a solver.
*/

#pragma once

#include "linearOperator.h"
#include <iostream>
#include <vector>
#include <boost/smart_ptr.hpp>
#include <sstream>
#include <stdexcept>

namespace Gadgetron{
template <class ARRAY_TYPE> class encodingOperatorContainer
  : public linearOperator<ARRAY_TYPE>
{
public:
  encodingOperatorContainer() : linearOperator<ARRAY_TYPE>() { num_elements_ = 0; }
  virtual ~encodingOperatorContainer(){}

  // Get domain and codomain dimensions. 
  // The codomain is a concatenation of the indivudial operators' domains.
  //
  using linearOperator< ARRAY_TYPE>::get_domain_dimensions;
  virtual boost::shared_ptr< std::vector<unsigned int> > get_codomain_dimensions() 
  { 
    std::vector<unsigned int> *dims = new std::vector<unsigned int>();
    dims->push_back( num_elements_ );
    return boost::shared_ptr< std::vector<unsigned int> >(dims);
  }

  // Get domain and codomain for the individual operators
  //
  virtual boost::shared_ptr< std::vector<unsigned int> > get_domain_dimensions(int i) { 
    return operators_[i]->get_domain_dimensions(); }
  
  virtual boost::shared_ptr< std::vector<unsigned int> > get_codomain_dimensions(int i) { 
    return operators_[i]->get_codomain_dimensions(); }
  
  // Allocate an array of the codomain dimensions
  //
  boost::shared_ptr< ARRAY_TYPE> create_codomain() 
  {
    std::vector<unsigned int> dims = get_codomain_dimensions();
    ARRAY_TYPE* codomain = new ARRAY_TYPE;
    codomain->create(&dims);
    return boost::shared_ptr<ARRAY_TYPE>(codomain);
  }

  // Concatenate a vector of codomains into a single array
  //
  boost::shared_ptr< ARRAY_TYPE> create_codomain( std::vector<ARRAY_TYPE*> codoms )
  {
    if (codoms.size() != operators_.size()){
      std::cerr << "encodingOperatorContainter::create_codomain: number of operators and number of codomains do no match" << std::endl;
      return boost::shared_ptr<ARRAY_TYPE>();
    }
    
    std::vector<unsigned int> dims = *get_codomain_dimensions();
    boost::shared_ptr<ARRAY_TYPE> codomain(new ARRAY_TYPE(&dims));


    int offset = 0;

    for (int i = 0; i < operators_.size(); i++){

      if (!codoms[i]->dimensions_equal(get_codomain_dimensions(i).get())){
      	 std::stringstream ss;
      	 ss << "encodingOperatorContainter::create_codomain: input codomain " << i << " does not match corresponding operator codomain" << std::endl;
      	 ss << "Input codomain: ";
      	 std::vector<unsigned int> ico = *codoms[i]->get_dimensions();
      	 for (int k = 0; k < ico.size(); k++) ss << ico[k] << " ";
      	 ss << std::endl;
      	 ss << "Operator codomain: ";
      	 ico = *get_codomain_dimensions(i);
      	 std::cout << "SIZE: " << ico.size() << std::endl;
      	 for (int k = 0; k < ico.size(); k++) ss << ico[k] << " ";
      	 ss << std::endl;
      	 BOOST_THROW_EXCEPTION(runtime_error(ss.str()));

      }

      ARRAY_TYPE slice;
      slice.create(codoms[i]->get_dimensions().get(),codomain->get_data_ptr()+offset);
      slice = *codoms[i];
      offset += slice.get_number_of_elements();
    }

    return codomain;    
  }

  // Get individual operators
  //
  boost::shared_ptr< linearOperator<ARRAY_TYPE> > get_operator(int i){
    return operators_[i];
  }

  // Get pointer offset into codomain for individual operators "sub-codomains"
  //
  unsigned int get_offset(int i){
    return offsets_[i];
  }

  // Add operator to the container
  //
  void add_operator( boost::shared_ptr< linearOperator<ARRAY_TYPE> > op )
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

  }
  
  virtual void mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
  {
    for (int i =0; i < operators_.size(); i++){

      ARRAY_TYPE tmp_data(operators_[i]->get_codomain_dimensions(),out->get_data_ptr()+offsets_[i]);
      operators_[i]->mult_M( in, &tmp_data, accumulate );
    }

  }

  virtual void mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
  {
    
    ARRAY_TYPE tmp_image(get_domain_dimensions());
        
    for (int i=0; i < operators_.size(); i++){
      
      boost::shared_ptr< linearOperator<ARRAY_TYPE> > op = operators_[i];

      ARRAY_TYPE tmp_data(op->get_codomain_dimensions(),in->get_data_ptr()+offsets_[i]);
      
      // This operator is special in that it needs to apply the "internal" operator weights
      //

      op->mult_MH( &tmp_data, &tmp_image );

      if( i == 0 && !accumulate ){
      	*out = tmp_image;
      	*out *= op->get_weight();
      }
      else {
      	axpy( op->get_weight(), &tmp_image, out );
      }
    }

  }
  
  virtual void mult_MH_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
  {

    ARRAY_TYPE tmp_image(get_domain_dimensions());
    
    for (int i=0; i < operators_.size(); i++){
      
      boost::shared_ptr< linearOperator<ARRAY_TYPE> > op = operators_[i];
      
      // This operator is special in that it needs to apply the "internal" operator weights
      //
      
      op->mult_MH_M( in, &tmp_image );
      if( i == 0 && !accumulate ){
      	*out = tmp_image;
      	*out *= op->get_weight();
      }
      else {
      	axpy( op->get_weight(), &tmp_image, out ) ;
      }
    }

  }
  virtual boost::shared_ptr< linearOperator< ARRAY_TYPE> > clone()
	{
      return linearOperator< ARRAY_TYPE >::clone(this);
	}

  
protected:
  std::vector< boost::shared_ptr< linearOperator<ARRAY_TYPE> > > operators_;
  std::vector<unsigned int> offsets_;
  unsigned int num_elements_;
};
}
