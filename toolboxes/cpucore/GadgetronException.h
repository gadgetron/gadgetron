#pragma once


#include <iostream>
#include <exception>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <boost/exception/exception.hpp>
#include <boost/exception/info.hpp>
#include <boost/exception/diagnostic_information.hpp>


typedef boost::error_info<struct err_message,std::string> error_message;



class gt_runtime_error: virtual public boost::exception, virtual public std::exception {
public:
	gt_runtime_error() : boost::exception(), std::exception(){}
	gt_runtime_error(std::string msg) : boost::exception(), std::exception(){
		(*this) << error_message(msg);
	}
}; //(2)

class gt_bad_alloc : virtual public gt_runtime_error {
	public:
		gt_bad_alloc(std::string msg) : gt_runtime_error(msg){}
		gt_bad_alloc() : gt_runtime_error(){}
};
