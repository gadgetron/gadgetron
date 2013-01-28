#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/exception/diagnostic_information.hpp>


class gt_runtime_error: virtual public boost::exception, virtual public std::exception {
public:
	gt_runtime_error() : boost::exception(), std::exception(), msg(0){}
	gt_runtime_error(std::string _msg) : boost::exception(), std::exception(), msg(_msg.c_str()){
	}
 virtual const  char * what() const throw(){
	 if (msg) return msg;
	 else return std::exception::what();
 }
protected:
 const char * msg;
}; //(2)

class gt_bad_alloc : public gt_runtime_error {
	public:
		gt_bad_alloc(std::string msg) : gt_runtime_error(msg){}
		gt_bad_alloc() : gt_runtime_error(){}
};
