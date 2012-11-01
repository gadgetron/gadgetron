/*
 * cunfft_test.cpp
 *
 *  Created on: Oct 19, 2012
 *      Author: Michael S. Hansen
 */

#include <iostream>
#include "cuNDArray.h"
#include "cuNDFFT.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include <string>
#include "FileInfo.h"
#include <complex>

int main(int argc, char** argv)
{
	if (argc < 2) {
		std::cout << "Usage: " << std::endl
				  << "    " << argv[0] << " <FILENAME.cplx>" << std::endl;
		return -1;
	}

	std::cout << "****** Simple test of the CUDA FFT wrappers ********" << std::endl;
	std::string filename(argv[1]);

	if (!FileInfo(filename).exists()) {
		std::cout << "Unable to find file " << filename << std::endl;
		return -1;
	}

	boost::shared_ptr<hoNDArray<float_complext> > data = read_nd_array<float_complext>(filename.c_str());

	std::cout << "Read array with dimensions [";
	for (unsigned int i = 0; i < data->get_number_of_dimensions(); i++) {
		std::cout << " " << data->get_size(i);
	}
	std::cout << " ]" << std::endl;

	cuNDArray<float_complext> deviceData(data.get());

	cuNDFFT<float_complext> ft;
	std::vector<unsigned int> ftdims(2,0); ftdims[1] = 1;
	ft.fft(&deviceData, &ftdims);
	ft.ifft(&deviceData, &ftdims);

	boost::shared_ptr<hoNDArray<float_complext> > data2 = deviceData.to_host();

	//Calculate RMS difference
	long double sq = 0.0;
	long double sum = 0.0;
	for (unsigned long int i = 0; i < data->get_number_of_elements(); i++) {
		double diff = abs(data->get_data_ptr()[i]-data2->get_data_ptr()[i]);
		sum += abs(data->get_data_ptr()[i]);
		sq += diff*diff;
	}
	sq /= sum;
	sq /= data->get_number_of_elements();
	std::cout << "RMS Difference: " << std::sqrt(sq) << std::endl;

	write_nd_array<float_complext>(data2.get(), "out.cplx");


	return 0;
}




