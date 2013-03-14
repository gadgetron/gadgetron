#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/exception/diagnostic_information.hpp>

namespace Gadgetron{
class runtime_error: virtual public boost::exception, virtual public std::exception {
public:
	runtime_error() : boost::exception(), std::exception(), msg(0){}
	runtime_error(std::string _msg) : boost::exception(), std::exception(), msg(_msg.c_str()){
	}
 virtual const  char * what() const throw(){
	 if (msg) return msg;
	 else return std::exception::what();
 }
protected:
 const char * msg;
}; //(2)

class bad_alloc : public runtime_error {
	public:
		bad_alloc(std::string msg) : runtime_error(msg){}
		bad_alloc() : runtime_error(){}
};
}
