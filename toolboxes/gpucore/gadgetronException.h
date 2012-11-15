#pragma once


#include <iostream>
#include <exception>
#include <stdexcept>

class cuda_error : public std::runtime_error
{
public:
	explicit cuda_error(const std::string& message) : std::runtime_error(message){};
};
