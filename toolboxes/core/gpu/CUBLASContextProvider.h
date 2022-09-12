/*
 * CUBLASContextProvider.h
 *
 *  Created on: Mar 22, 2012
 *      Author: Michael S. Hansen
 */

#ifndef CUBLASCONTEXTPROVIDER_H_
#define CUBLASCONTEXTPROVIDER_H_
#pragma once

#include <cublas_v2.h>
#include <map>
#include <iostream>

class CUBLASContextProvider
{

public:
	static CUBLASContextProvider* instance();

	cublasHandle_t* getCublasHandle(int device_no = 0);

private:
	CUBLASContextProvider() {}
	virtual ~CUBLASContextProvider();

	static CUBLASContextProvider* instance_;

	std::map<int, cublasHandle_t> handles_;
};

#endif /* CUBLASCONTEXTPROVIDER_H_ */
