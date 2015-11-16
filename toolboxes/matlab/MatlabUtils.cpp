#include "MatlabUtils.h"

#include "hoNDArray_math.h"

using namespace Gadgetron;

template<class T> struct isComplex { static constexpr mxComplexity value = mxREAL;};
template<class REAL> struct isComplex<complext<REAL>> { static constexpr mxComplexity value = mxCOMPLEX;};
template<class REAL> struct isComplex<std::complex<REAL>>{ static constexpr mxComplexity value = mxCOMPLEX;};



template<class T>  struct  MatlabClassID {};

template<>	struct MatlabClassID<double>{ static constexpr mxClassID value =  mxDOUBLE_CLASS;};
template<>	struct MatlabClassID<float>{ static constexpr mxClassID value =  mxSINGLE_CLASS;};

template<class REAL> struct MatlabClassID<complext<REAL>>{ static constexpr mxClassID value =  MatlabClassID<REAL>::value;};
template<class REAL> struct MatlabClassID<std::complex<REAL>>{ static constexpr mxClassID value =  MatlabClassID<REAL>::value;};

template<>	struct MatlabClassID<int8_t>{ static constexpr mxClassID value =  mxINT8_CLASS;};
template<>	struct MatlabClassID<uint8_t>{ static constexpr mxClassID value =  mxUINT8_CLASS;};
template<>	struct MatlabClassID<int16_t>{ static constexpr mxClassID value =  mxINT16_CLASS;};
template<>	struct MatlabClassID<uint16_t>{ static constexpr mxClassID value =  mxUINT16_CLASS;};
template<>	struct MatlabClassID<int32_t>{ static constexpr mxClassID value =  mxINT32_CLASS;};
template<>	struct MatlabClassID<uint32_t>{ static constexpr mxClassID value =  mxUINT32_CLASS;};
template<>	struct MatlabClassID<int64_t>{ static constexpr mxClassID value =  mxINT64_CLASS;};
template<>	struct MatlabClassID<uint64_t>{ static constexpr mxClassID value =  mxUINT64_CLASS;};
template<>	struct MatlabClassID<size_t>{ static constexpr mxClassID value =  mxUINT64_CLASS;};

template<class T> struct MatlabConverter {
	static mxArray* convert(hoNDArray<T>* input){

		mwSize ndim = input->get_number_of_dimensions();
		mwSize* dims = new mwSize[ndim];
		for (size_t i = 0; i < ndim; i++)
			dims[i] = input->get_size(i);

		T* raw_data = (T*) mxCalloc(input->get_number_of_elements(),sizeof(T));
		memcpy(raw_data,input->get_data_ptr(),input->get_number_of_bytes());
		auto result =  mxCreateNumericArray(ndim,dims,MatlabClassID<T>::value,isComplex<T>::value);
		mxSetData(result,raw_data);
		return result;

	}

	static hoNDArray<T> convert(mxArray* input) {
		auto ndims = mxGetNumberOfDimensions(input);
		auto dims = mxGetDimensions(input);
		std::vector<size_t> dimensions(ndims);
		for (size_t i = 0; i <ndims; i++) dimensions[i] = dims[i];

		auto result =  hoNDArray<T>(dimensions);

		if (mxGetImagData(input)) //This is for REAL data only
			throw std::runtime_error("Trying to convert complex matlab data to non-complex c++ type");
		if (mxGetClassID(input) == MatlabClassID<T>::value ){ //Same type, so we can just memcpy
			T* raw_data = (T*) mxGetData(input);
			memcpy(result.get_data_ptr(),raw_data,result.get_number_of_elements()*sizeof(T));
		} else {
			switch (mxGetClassID(input)){ // Have to do runtime type conversion, which means cases en-masse.
			case MatlabClassID<double>::value :
			copyMatlabdata<double>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<float>::value:
			copyMatlabdata<float>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<int8_t>::value:
			copyMatlabdata<int8_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<uint8_t>::value:
			copyMatlabdata<uint8_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<int16_t>::value:
			copyMatlabdata<int16_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<uint16_t>::value:
			copyMatlabdata<uint16_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<int32_t>::value:
			copyMatlabdata<int32_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<uint32_t>::value:
			copyMatlabdata<uint32_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;
			case MatlabClassID<int64_t>::value:
			copyMatlabdata<int64_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			case MatlabClassID<uint64_t>::value:
			copyMatlabdata<uint64_t>(input,result.get_data_ptr(),result.get_number_of_elements());
			break;

			default:
				throw std::runtime_error("Trying to convert from unsupported data type");
				break;


			}

		}
		return result;
	}

	template<class R> static void copyMatlabdata(mxArray* input, T* output,size_t len){
		R* raw_ptr = (R*) mxGetData(input);
		for (size_t i = 0; i < len; i++)
			output[i] = T(raw_ptr[i]);
	}
};

template<class REAL,unsigned int N> struct MatlabConverter<vector_td<REAL,N>> {
  static mxArray* convert(hoNDArray<vector_td<REAL,N>>* input){
    std::vector<size_t> dims = *input->get_dimensions();
    std::vector<size_t> dims2;
    dims2.push_back(N);
    for (auto s : dims) dims2.push_back(s);

    hoNDArray<REAL> tmp(dims2,(REAL*)input->get_data_ptr());
    return MatlabConverter<REAL>::convert(&tmp);
  }
  static hoNDArray<vector_td<REAL,N>> convert(mxArray* matarray){

    auto tmp = MatlabConverter<REAL>::convert(matarray);
    auto dims = *tmp.get_dimensions();
    if (dims[0] != N)
      throw std::runtime_error("Converting from Matlab array to hoNDArray with vector_td, but sizes don't match");
    std::vector<size_t> dims2(dims.begin()+1,dims.end());
    auto result = hoNDArray<vector_td<REAL,N>>(dims2);
    memcpy(result.get_data_ptr(),tmp.get_data_ptr(),result.get_number_of_bytes());
    return result;

  }
};
template<class REAL> struct MatlabConverter<complext<REAL>> {
	static mxArray* convert(hoNDArray<complext<REAL>>* input){

		size_t ndim = input->get_number_of_dimensions();

		//Matlab does not support to creation of 7D arrays, but 8,6 and 9 works just fine.
		//If you're on a train that's running Matlab as its control system, you should be very very scared.

		mwSize* dims = new mwSize[ndim];
		for (size_t i = 0; i < ndim; i++)
			dims[i] = input->get_size(i);

		REAL* real_data = (REAL*) mxCalloc(input->get_number_of_elements(),sizeof(REAL));
		REAL* imag_data = (REAL*) mxCalloc(input->get_number_of_elements(),sizeof(REAL));

		complext<REAL>* raw_data = input->get_data_ptr();
		for (size_t i = 0; i < input->get_number_of_elements(); i++){
			real_data[i] = real(raw_data[i]);
			imag_data[i] = imag(raw_data[i]);
		}

		auto result  =  mxCreateNumericArray(ndim,dims,MatlabClassID<REAL>::value,isComplex<complext<REAL>>::value);
		mxSetData(result,real_data);
		mxSetImagData(result,imag_data);

		auto ndims_test = mxGetNumberOfDimensions(result);

		return result;
	}
	static hoNDArray<complext<REAL> > convert(mxArray* input) {
		auto ndims = mxGetNumberOfDimensions(input);
		auto dims = mxGetDimensions(input);
		std::vector<size_t> dimensions(ndims);
		for (size_t i = 0; i <ndims; i++) dimensions[i] = dims[i];
		auto result = hoNDArray<complext<REAL>>(dimensions);
		switch (mxGetClassID(input)){ // Have to do runtime type conversion, which means cases en-masse.
		case MatlabClassID<double>::value :
		copyMatlabdata<double>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<float>::value:
		copyMatlabdata<float>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<int8_t>::value:
		copyMatlabdata<int8_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<uint8_t>::value:
		copyMatlabdata<uint8_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<int16_t>::value:
		copyMatlabdata<int16_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<uint16_t>::value:
		copyMatlabdata<uint16_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<int32_t>::value:
		copyMatlabdata<int32_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<uint32_t>::value:
		copyMatlabdata<uint32_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;
		case MatlabClassID<int64_t>::value:
		copyMatlabdata<int64_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		case MatlabClassID<uint64_t>::value:
		copyMatlabdata<uint64_t>(input,result.get_data_ptr(),result.get_number_of_elements());
		break;

		default:
			throw std::runtime_error("Trying to convert from unsupported data type");
			break;

		}

		return result;

	}

	template<class R> static void copyMatlabdata(mxArray* input, complext<REAL>* output,size_t len){
		R* real_ptr = (R*) mxGetData(input);
		R* imag_ptr = (R*) mxGetImagData(input);
		if (imag_ptr) {
			for (size_t i = 0; i < len; i++)
				output[i] = complext<REAL>(REAL(real_ptr[i]),REAL(imag_ptr[i]));
		} else{
			for (size_t i = 0; i < len; i++)
				output[i] = complext<REAL>(REAL(real_ptr[i]),0);
		}
	}

};

template<class REAL> struct MatlabConverter<std::complex<REAL>> {
	static mxArray* convert(hoNDArray<std::complex<REAL>>* input){
		return MatlabConverter<complext<REAL>>::convert((hoNDArray<complext<REAL>>*) input);
	}

	static hoNDArray<std::complex<REAL>> convert(mxArray* input){
		auto tmp = MatlabConverter<complext<REAL>>::convert(input);
		return std::move(*((hoNDArray<std::complex<REAL>>*)&tmp));
	}
};

template<class T> mxArray* Gadgetron::hoNDArrayToMatlab(hoNDArray<T> * input){
	return MatlabConverter<T>::convert(input);

}


template<class T> hoNDArray<T> Gadgetron::MatlabToHoNDArray(mxArray* data){
	return MatlabConverter<T>::convert(data);
}

mxArray* Gadgetron::BufferToMatlabStruct(IsmrmrdDataBuffered* buffer){

	const char * field_names[] = {"data","trajectory","headers","samplingdescription"};
	mwSize one = 1;
	auto mxstruct = mxCreateStructArray(1,&one,4,field_names);


	if (!mxstruct) throw std::runtime_error("Failed to allocate Matlab struct");

	auto mxdata = hoNDArrayToMatlab(&buffer->data_);
	mxSetField(mxstruct,0,"data",mxdata);
	//Add trajectory if available
	if (buffer->trajectory_){
		auto & trajectory = *buffer->trajectory_;
		int traj_fieldnumber = mxAddField(mxstruct,"trajectory");
		auto mxtraj = hoNDArrayToMatlab(&trajectory);
		mxSetFieldByNumber(mxstruct,0,traj_fieldnumber,mxtraj);
	}

	//Add headers
	std::cout << "Adding headers " << std::endl;
	mwSize num_headers = buffer->headers_.get_number_of_elements();
	auto mxheaders = mxCreateNumericMatrix(sizeof(ISMRMRD::AcquisitionHeader),num_headers,mxUINT8_CLASS,mxREAL);
	memcpy(mxGetData(mxheaders),buffer->headers_.get_data_ptr(),sizeof(ISMRMRD::AcquisitionHeader)*num_headers);
	mxSetField(mxstruct,0,"headers",mxheaders);

	auto samplingdescription = samplingdescriptionToMatlabStruct(&buffer->sampling_);
	mxSetField(mxstruct,0,"samplingdescription",samplingdescription);

	return mxstruct;


}
static SamplingDescription MatlabStructToSamplingdescription(mxArray* mxstruct){

	SamplingDescription samp;
	auto encFOV = mxGetField(mxstruct,0,"encoded_FOV");
	std::cout << "FOV PTR " << encFOV << " " << mxGetData(encFOV)<< std::endl;
	memcpy(samp.encoded_FOV_,mxGetData(encFOV),sizeof(samp.encoded_FOV_));
	auto recFOV = mxGetField(mxstruct,0,"recon_FOV");
	memcpy(samp.recon_FOV_,mxGetData(recFOV),sizeof(samp.recon_FOV_));
	auto encoded_matrix = mxGetField(mxstruct,0,"encoded_matrix");
	memcpy(samp.encoded_matrix_,mxGetData(encoded_matrix),sizeof(samp.encoded_matrix_));
	auto recon_matrix = mxGetField(mxstruct,0,"recon_matrix");
	memcpy(samp.recon_matrix_,mxGetData(recon_matrix),sizeof(samp.recon_matrix_));
	auto sampling_limit = mxGetField(mxstruct,0,"sampling_limits");
	memcpy(samp.sampling_limits_,mxGetData(sampling_limit),sizeof(samp.sampling_limits_));

	return samp;

}


IsmrmrdDataBuffered Gadgetron::MatlabStructToBuffer(mxArray* mxstruct){
	IsmrmrdDataBuffered buffer;

	auto data = mxGetField(mxstruct,0,"data");
	buffer.data_ = MatlabToHoNDArray<std::complex<float>>(data);
	if (buffer.data_.get_number_of_dimensions() != 7){ //Someone (Matlab) got rid of our dimensions. Ghee thanks;
		std::vector<size_t> newdims = *buffer.data_.get_dimensions();
		for (int i = buffer.data_.get_number_of_dimensions(); i<7; i++)
			newdims.push_back(1);
		buffer.data_.reshape(&newdims);
	}
	auto traj = mxGetField(mxstruct,0,"trajectory");
	if (traj){
		buffer.trajectory_ = MatlabToHoNDArray<float>(traj);
		auto & trajectory = *buffer.trajectory_;
		if (trajectory.get_number_of_dimensions() != 7){
			std::vector<size_t> newdims = *trajectory.get_dimensions();
			for (int i = trajectory.get_number_of_dimensions(); i<7; i++)
				newdims.push_back(1);
			trajectory.reshape(&newdims);
		}
	}
	auto headers = mxGetField(mxstruct,0,"headers");

	auto nmat_headers = mxGetN(headers);
	std::vector<size_t> header_dim = {buffer.data_.get_size(1),buffer.data_.get_size(2),buffer.data_.get_size(4),buffer.data_.get_size(5),buffer.data_.get_size(6)};

	buffer.headers_ = hoNDArray<ISMRMRD::AcquisitionHeader>(header_dim);

	std::cout << "Number of headers: " << nmat_headers << " Expected: " << buffer.headers_.get_number_of_elements() << std::endl;
	if (nmat_headers != buffer.headers_.get_number_of_elements())
		throw std::runtime_error("Number of headers does not match number of kspace acquisitions");

	memcpy(buffer.headers_.get_data_ptr(),mxGetData(headers),sizeof(ISMRMRD::AcquisitionHeader)*buffer.headers_.get_number_of_elements());

	auto samplingdescription = mxGetField(mxstruct,0,"samplingdescription");
	buffer.sampling_ = MatlabStructToSamplingdescription(samplingdescription);
	return buffer;


}


mxArray* Gadgetron::samplingdescriptionToMatlabStruct(SamplingDescription* samp){

	const char* fieldnames[5] = {"encoded_FOV","recon_FOV","encoded_matrix","recon_matrix","sampling_limits"};
	mwSize one_dim  = 1;
	auto sampStruct = mxCreateStructArray(1,&one_dim,5,fieldnames);
	//Encoded FOV
	mwSize dims = 3;
	auto encFOV = mxCreateNumericArray(1,&dims,MatlabClassID<float>::value,mxComplexity(0));
	memcpy(mxGetData(encFOV),samp->encoded_FOV_,sizeof(samp->encoded_FOV_));
	mxSetField(sampStruct,0,"encoded_FOV",encFOV);
	//Recon FOV
	auto recFOV = mxCreateNumericArray(1,&dims,MatlabClassID<float>::value,mxComplexity(0));
	memcpy(mxGetData(recFOV),samp->recon_FOV_,sizeof(samp->recon_FOV_));
	mxSetField(sampStruct,0,"recon_FOV",recFOV);
	//Encoded Matrix
	auto encoded_matrix = mxCreateNumericArray(1,&dims,MatlabClassID<uint16_t>::value,mxComplexity(0));
	memcpy(mxGetData(encoded_matrix),samp->encoded_matrix_,sizeof(samp->encoded_matrix_));
	mxSetField(sampStruct,0,"encoded_matrix",encoded_matrix);
	//Recon matrix
	auto recon_matrix = mxCreateNumericArray(1,&dims,MatlabClassID<uint16_t>::value,mxComplexity(0));
	memcpy(mxGetData(recon_matrix),samp->recon_matrix_,sizeof(samp->recon_matrix_));
	mxSetField(sampStruct,0,"recon_matrix",recon_matrix);
	//Sampling Limit
	mwSize twodims[] = {3,3};
	auto sampling_limit = mxCreateNumericArray(2,twodims,MatlabClassID<uint16_t>::value,mxComplexity(0));
	memcpy(mxGetData(sampling_limit),samp->sampling_limits_,sizeof(samp->sampling_limits_));
	mxSetField(sampStruct,0,"sampling_limits",sampling_limit);
	return sampStruct;
}



template mxArray* Gadgetron::hoNDArrayToMatlab<float>(hoNDArray<float> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<double>(hoNDArray<double> *);

template mxArray* Gadgetron::hoNDArrayToMatlab<size_t>(hoNDArray<size_t> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<float_complext>(hoNDArray<float_complext> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<double_complext>(hoNDArray<double_complext> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<std::complex<double>>(hoNDArray<std::complex<double>> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<std::complex<float>>(hoNDArray<std::complex<float>> *);


template hoNDArray<float> Gadgetron::MatlabToHoNDArray<float>(mxArray *);
template hoNDArray<double> Gadgetron::MatlabToHoNDArray<double>(mxArray *);
template hoNDArray<size_t> Gadgetron::MatlabToHoNDArray<size_t>(mxArray *);
template hoNDArray<float_complext> Gadgetron::MatlabToHoNDArray<float_complext>(mxArray *);
template hoNDArray<double_complext> Gadgetron::MatlabToHoNDArray<double_complext>(mxArray *);

template hoNDArray<std::complex<double>> Gadgetron::MatlabToHoNDArray<std::complex<double>>(mxArray *);
template hoNDArray<std::complex<float>> Gadgetron::MatlabToHoNDArray<std::complex<float>>(mxArray *);

template mxArray* Gadgetron::hoNDArrayToMatlab<vector_td<float,1>>(hoNDArray<vector_td<float,1>> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<vector_td<float,2>>(hoNDArray<vector_td<float,2>> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<vector_td<float,3>>(hoNDArray<vector_td<float,3>> *);
template mxArray* Gadgetron::hoNDArrayToMatlab<vector_td<float,4>>(hoNDArray<vector_td<float,4>> *);


template hoNDArray<vector_td<float,1>> Gadgetron::MatlabToHoNDArray<vector_td<float,1>>(mxArray *);
template hoNDArray<vector_td<float,2>> Gadgetron::MatlabToHoNDArray<vector_td<float,2>>(mxArray *);
template hoNDArray<vector_td<float,3>> Gadgetron::MatlabToHoNDArray<vector_td<float,3>>(mxArray *);
template hoNDArray<vector_td<float,4>> Gadgetron::MatlabToHoNDArray<vector_td<float,4>>(mxArray *);
