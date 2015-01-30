/** \file GadgetronException.h
    \brief An interface to the exception handling used in the Gadgetron to indicate runtime errors.
*/

#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>

namespace Gadgetron{

  class runtime_error: virtual public std::exception 
  {
  public:
    runtime_error() : std::exception(), msg(0){}
    runtime_error(std::string _msg) : std::exception(), msg(_msg.c_str()){
    }
    virtual const  char * what() const throw(){
      if (msg) return msg;
      else return std::exception::what();
    }
  protected:
    const char * msg;
  };
  
  class bad_alloc : public runtime_error 
  {
  public:
    bad_alloc(std::string msg) : runtime_error(msg){}
    bad_alloc() : runtime_error(){}
  };
}
