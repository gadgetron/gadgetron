/*
 * cudaDeviceManager.h
 *
 *  Created on: Jan 25, 2013
 *      Author: Dae
 */

#pragma once
#include <vector>
#include "cublas_v2.h"
class cudaDeviceManager {
private:
	cudaDeviceManager();
	~cudaDeviceManager();
	static void CleanUp();

	int num_devices;
	std::vector<int> _warp_size;
	std::vector<int> _max_blockdim;
	std::vector<int> _max_griddim;
	std::vector<cublasHandle_t>  handle;
	static cudaDeviceManager * _instance;
	cudaDeviceManager(cudaDeviceManager const&);
	cudaDeviceManager& operator=(cudaDeviceManager const&);

public:
	static cudaDeviceManager * Instance();
	int warp_size(int device){return _warp_size[device];};
	int max_blockdim(int device){return _max_blockdim[device];};
	int max_griddim(int device){return _max_griddim[device];};
	int warp_size();
	int max_blockdim();
	int max_griddim();

	size_t getFreeMemory();
	size_t getFreeMemory(int device);
	size_t getTotalMemory();
	size_t getTotalMemory(int device);
	cublasHandle_t getHandle();
	cublasHandle_t getHandle(int device){return handle[device];};

};


