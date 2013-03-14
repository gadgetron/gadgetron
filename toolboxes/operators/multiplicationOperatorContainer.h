#pragma once

#include "linearOperator.h"
#include <iostream>
#include <vector>

template <class REAL, class ARRAY_TYPE> class multiplicationOperatorContainer
  : public linearOperator<REAL, ARRAY_TYPE>
{
public:
  multiplicationOperatorContainer() : linearOperator<REAL,ARRAY_TYPE>() {}
  virtual ~multiplicationOperatorContainer(){}

  // Set/get domain and codomain dimensions. 
  //

  virtual bool set_domain_dimensions( std::vector<unsigned int>* )
  { 
    std::cerr << "Warning: multiplicationOperatorContainer::set_domain_dimensions : dimensions ignored, using dimensions of the individual operators instead" << std::endl;
    return false;
  }  

  virtual bool set_codomain_dimensions( std::vector<unsigned int> *dims ) 
  { 
    std::cerr << "Warning: multiplicationOperatorContainer::set_codomain_dimensions : dimensions ignored, using dimensions of the individual operators instead" << std::endl;
    return false;
  }

  virtual boost::shared_ptr< std::vector<unsigned int> > get_domain_dimensions() 
  { 
    if( operators_.size() == 0 )
      return boost::shared_ptr< std::vector<unsigned int> >();
    else
      return operators_[0]->get_domain_dimensions();
  }

  virtual boost::shared_ptr< std::vector<unsigned int> > get_codomain_dimensions() 
  { 
    if( operators_.size() == 0 )
      return boost::shared_ptr< std::vector<unsigned int> >();
    else
      return operators_[operators_.size()-1]->get_codomain_dimensions();
  }
  
  virtual void set_weight( REAL weight ){ 
    REAL op_weight = REAL(1);
    for( int i=0; i<operators_.size(); i++ )
      op_weight *= operators_[i]->get_weight();
    this->weight_ = weight*op_weight;
  }

  // Add operator to the container
  //
  bool add_operator( boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > op )
  {
    if( op.get() == 0x0 ){
      std::cerr << "Error: multiplicationOperatorContainer::add_operator : illegal operator" << std::endl;
      return false;
    } 

    // All operators needs the domain and codomain dimensions set
    //
    if( op->get_domain_dimensions()->size() == 0 ){
      std::cerr << "Error: multiplicationOperatorContainer::add_operator : domain dimensions not set on operator" << std::endl;
      return false;
    }
    if( op->get_codomain_dimensions()->size() == 0 ){
      std::cerr << "Error: multiplicationOperatorContainer::add_operator : codomain dimensions not set on operator" << std::endl;
      return false;
    }

    if( operators_.size() == 0 && !_set_domain_dimensions( op->get_domain_dimensions().get() ) ){
      std::cerr << "Error: multiplicationOperatorContainer::add_operator : failed to set domain dimensions on container" << std::endl;
      return false;
    }

    if( !_set_codomain_dimensions( op->get_codomain_dimensions().get() ) ){
      std::cerr << "Error: multiplicationOperatorContainer::add_operator : failed to set codomain dimensions on container" << std::endl;
      return false;
    }
    
    operators_.push_back( op );
    this->weight_ *= op->get_weight();

    return true;
  }
  
  virtual int mult_M( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
  {
    if( operators_.size() == 0 ){
      std::cerr << "Error: multiplicationOperatorContainer::mult_M : no operators added" << std::endl;
      return -1;
    }
    
    ARRAY_TYPE *tmp_in = in, *tmp_out = 0x0;
    ARRAY_TYPE ping, pong;

    if( operators_.size() > 1 ){
      if( !ping.create( operators_[0]->get_codomain_dimensions().get() )){
	std::cerr << "Error: multiplicationOperatorContainer::mult_M : failed to create intermediate array (1)" << std::endl;
	return -1;
      }
      tmp_out = &ping;
    }
    else{
      tmp_out = out;
    }
    
    // Loop over operators
    //
    for( int i=0; i < operators_.size(); i++ ){
      
      if( operators_[i]->mult_M( tmp_in, tmp_out, (i==operators_.size()-1) ? accumulate : false ) < 0 ){
	std::cerr << "Error: multiplicationOperatorContainer : mult_M failed on sub-operator" << std::endl;
	return -1;
      }
      
      ARRAY_TYPE *tmp_tmp_out = (i==0) ? &pong : tmp_in;
      tmp_in = tmp_out;

      if( operators_.size() > 2 && i < operators_.size()-2 ){
	if( !tmp_tmp_out->create( operators_[i+1]->get_codomain_dimensions().get() )){
	  std::cerr << "Error: multiplicationOperatorContainer::mult_M : failed to create intermediate array (2)" << std::endl;
	  return -1;
	}
	tmp_out = tmp_tmp_out;
      }
      else if( i == operators_.size()-2 ){
	tmp_out = out;
      }      
    }
    return 0;
  }

  virtual int mult_MH( ARRAY_TYPE* in, ARRAY_TYPE* out, bool accumulate = false )
  {
    if( operators_.size() == 0 ){
      std::cerr << "Error: multiplicationOperatorContainer::mult_MH : no operators added" << std::endl;
      return -1;
    }
    
    ARRAY_TYPE *tmp_in = in, *tmp_out = 0x0;
    ARRAY_TYPE ping, pong;
    
    if( operators_.size() > 1 ){
      if( !ping.create( operators_[operators_.size()-1]->get_domain_dimensions().get() )){
	std::cerr << "Error: multiplicationOperatorContainer::mult_MH : failed to create intermediate array (1)" << std::endl;
	return -1;
      }
      tmp_out = &ping;
    }
    else{
      tmp_out = out;
    }
    
    // Loop over operators
    //
    for( int i=operators_.size()-1; i>=0; i-- ){
      
      if( operators_[i]->mult_MH( tmp_in, tmp_out, (i==0) ? accumulate : false ) < 0 ){
	std::cerr << "Error: multiplicationOperatorContainer : mult_MH failed on sub-operator" << std::endl;
	return -1;
      }
      
      ARRAY_TYPE *tmp_tmp_out = (i==operators_.size()-1) ? &pong : tmp_in;
      tmp_in = tmp_out;
      
      if( i > 1 ){
	if( !tmp_tmp_out->create( operators_[i-1]->get_domain_dimensions().get() )){
	  std::cerr << "Error: multiplicationOperatorContainer::mult_MH : failed to create intermediate array (2)" << std::endl;
	  return -1;
	}
	tmp_out = tmp_tmp_out;
      }
      else if( i == 1 ){
	tmp_out = out;
      }      
    }
    return 0;
  }

protected:

  virtual bool _set_domain_dimensions( std::vector<unsigned int> *dims )
  { 
    return linearOperator<REAL, ARRAY_TYPE>::set_domain_dimensions( dims );
  }  

  virtual bool _set_codomain_dimensions( std::vector<unsigned int> *dims )
  { 
    return linearOperator<REAL, ARRAY_TYPE>::set_codomain_dimensions( dims );
  }  
  
protected:
  std::vector< boost::shared_ptr< linearOperator<REAL, ARRAY_TYPE> > > operators_;
};
