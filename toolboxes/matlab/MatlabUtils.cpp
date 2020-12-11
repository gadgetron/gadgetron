#include "MatlabUtils.h"

#include "hoNDArray_math.h"

namespace Gadgetron{

template<class T> struct isComplex { static constexpr mxComplexity value = mxREAL;};
template<class REAL> struct isComplex<complext<REAL> > { static constexpr mxComplexity value = mxCOMPLEX;};
template<class REAL> struct isComplex<std::complex<REAL> >{ static constexpr mxComplexity value = mxCOMPLEX;};



template<class T>  struct  MatlabClassID {};

template<>	struct MatlabClassID<double>{ static constexpr mxClassID value =  mxDOUBLE_CLASS;};
template<>	struct MatlabClassID<float>{ static constexpr mxClassID value =  mxSINGLE_CLASS;};

template<class REAL> struct MatlabClassID<complext<REAL> >{ static constexpr mxClassID value =  MatlabClassID<REAL>::value;};
template<class REAL> struct MatlabClassID<std::complex<REAL> >{ static constexpr mxClassID value =  MatlabClassID<REAL>::value;};

template<>	struct MatlabClassID<int8_t>{ static constexpr mxClassID value =  mxINT8_CLASS;};
template<>	struct MatlabClassID<uint8_t>{ static constexpr mxClassID value =  mxUINT8_CLASS;};
template<>	struct MatlabClassID<int16_t>{ static constexpr mxClassID value =  mxINT16_CLASS;};
template<>	struct MatlabClassID<uint16_t>{ static constexpr mxClassID value =  mxUINT16_CLASS;};
template<>	struct MatlabClassID<int32_t>{ static constexpr mxClassID value =  mxINT32_CLASS;};
template<>	struct MatlabClassID<uint32_t>{ static constexpr mxClassID value =  mxUINT32_CLASS;};
template<>	struct MatlabClassID<int64_t>{ static constexpr mxClassID value =  mxINT64_CLASS;};
template<>	struct MatlabClassID<unsigned long>{ static constexpr mxClassID value =  mxUINT64_CLASS;};
template<>	struct MatlabClassID<unsigned long long>{ static constexpr mxClassID value =  mxUINT64_CLASS;};

#ifdef MX_HAS_INTERLEAVED_COMPLEX
template <typename T> struct matlab_type_traits {};
template<> struct matlab_type_traits<double> {
     static constexpr mxDouble* (*get_real)(const mxArray*) = mxGetDoubles;
     static constexpr mxComplexDouble* (*get_cmplx)(const mxArray*) = mxGetComplexDoubles;
};
template<> struct matlab_type_traits<float> {
     static constexpr mxSingle* (*get_real)(const mxArray*) = mxGetSingles;
     static constexpr mxComplexSingle* (*get_cmplx)(const mxArray*) = mxGetComplexSingles;
};
template<> struct matlab_type_traits<int8_t> {
     static constexpr mxInt8* (*get_real)(const mxArray*) = mxGetInt8s;
     static constexpr mxComplexInt8* (*get_cmplx)(const mxArray*) = mxGetComplexInt8s;
};
template<> struct matlab_type_traits<uint8_t> {
     static constexpr mxUint8* (*get_real)(const mxArray*) = mxGetUint8s;
     static constexpr mxComplexUint8* (*get_cmplx)(const mxArray*) = mxGetComplexUint8s;
};
template<> struct matlab_type_traits<int16_t> {
     static constexpr mxInt16* (*get_real)(const mxArray*) = mxGetInt16s;
     static constexpr mxComplexInt16* (*get_cmplx)(const mxArray*) = mxGetComplexInt16s;
};
template<> struct matlab_type_traits<uint16_t> {
     static constexpr mxUint16* (*get_real)(const mxArray*) = mxGetUint16s;
     static constexpr mxComplexUint16* (*get_cmplx)(const mxArray*) = mxGetComplexUint16s;
};
template<> struct matlab_type_traits<int32_t> {
     static constexpr mxInt32* (*get_real)(const mxArray*) = mxGetInt32s;
     static constexpr mxComplexInt32* (*get_cmplx)(const mxArray*) = mxGetComplexInt32s;
};
template<> struct matlab_type_traits<uint32_t> {
     static constexpr mxUint32* (*get_real)(const mxArray*) = mxGetUint32s;
     static constexpr mxComplexUint32* (*get_cmplx)(const mxArray*) = mxGetComplexUint32s;
};
template<> struct matlab_type_traits<int64_t> {
     static constexpr mxInt64* (*get_real)(const mxArray*) = mxGetInt64s;
     static constexpr mxComplexInt64* (*get_cmplx)(const mxArray*) = mxGetComplexInt64s;
};
template<> struct matlab_type_traits<unsigned long> {
     static constexpr mxUint64* (*get_real)(const mxArray*) = mxGetUint64s;
     static constexpr mxComplexUint64* (*get_cmplx)(const mxArray*) = mxGetComplexUint64s;
};
template<> struct matlab_type_traits<unsigned long long> {
     static constexpr mxUint64* (*get_real)(const mxArray*) = mxGetUint64s;
     static constexpr mxComplexUint64* (*get_cmplx)(const mxArray*) = mxGetComplexUint64s;
};
#endif /* MX_HAS_INTERLEAVED_COMPLEX */

// ------------------------
// hoNDArray
// ------------------------

template<class T> struct MatlabConverter {
	static mxArray* convert(hoNDArray<T>* input){

		mwSize ndim = input->get_number_of_dimensions();
		mwSize* dims = new mwSize[ndim];
		for (size_t i = 0; i < ndim; i++)
			dims[i] = input->get_size(i);

		T* raw_data = (T*) mxCalloc(input->get_number_of_elements(),sizeof(T));
		memcpy(raw_data,input->get_data_ptr(),input->get_number_of_bytes());
		auto result =  mxCreateNumericMatrix(0,0,MatlabClassID<T>::value,isComplex<T>::value);
		mxSetDimensions(result,dims,ndim);
		mxSetData(result,raw_data);
		return result;

	}

	static hoNDArray<T> convert(mxArray* input) {
		auto ndims = mxGetNumberOfDimensions(input);
		auto dims = mxGetDimensions(input);
		std::vector<size_t> dimensions(ndims);
		for (size_t i = 0; i <ndims; i++) dimensions[i] = dims[i];

		auto result =  hoNDArray<T>(dimensions);

		if (mxIsComplex(input)) //This is for REAL data only
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

		complext<REAL>* raw_data = input->get_data_ptr();
#ifdef MX_HAS_INTERLEAVED_COMPLEX
		auto result = mxCreateNumericArray(ndim,dims,MatlabClassID<REAL>::value,isComplex<complext<REAL>>::value);

		auto data = matlab_type_traits<REAL>::get_cmplx(result);
		for (size_t i = 0; i < input->get_number_of_elements(); i++){
			data[i].real = real(raw_data[i]);
			data[i].imag = imag(raw_data[i]);
		}
#else
		auto result = mxCreateNumericMatrix(0,0,MatlabClassID<REAL>::value,isComplex<complext<REAL>>::value);

		REAL* real_data = (REAL*) mxCalloc(input->get_number_of_elements(),sizeof(REAL));
		REAL* imag_data = (REAL*) mxCalloc(input->get_number_of_elements(),sizeof(REAL));

		for (size_t i = 0; i < input->get_number_of_elements(); i++){
			real_data[i] = real(raw_data[i]);
			imag_data[i] = imag(raw_data[i]);
		}
		mxSetDimensions(result,dims,ndim);
		mxSetData(result,real_data);
		mxSetImagData(result,imag_data);
#endif /* MX_HAS_INTERLEAVED_COMPLEX */

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

#ifdef MX_HAS_INTERLEAVED_COMPLEX
	template<class R> static void copyMatlabdata(mxArray* input, complext<REAL>* output,size_t len){
		if (mxIsComplex(input)) {
			auto* data = matlab_type_traits<R>::get_cmplx(input);
			for (size_t i = 0; i < len; i++) {
			        output[i]._real = data[i].real;
			        output[i]._imag = data[i].imag;
			}
		} else{
			auto* data = matlab_type_traits<R>::get_real(input);
			for (size_t i = 0; i < len; i++)
				output[i] = data[i];
		}
	}
#else
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
#endif /* MX_HAS_INTERLEAVED_COMPLEX */

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

template<class T> mxArray* hoNDArrayToMatlab(hoNDArray<T> * input){
    return MatlabConverter<T>::convert(input);

}


template<class T> hoNDArray<T> MatlabToHoNDArray(mxArray* data){
    return MatlabConverter<T>::convert(data);
}

template<class T> void MatlabToHoNDArray(mxArray* data, hoNDArray<T>& a)
{
    a = MatlabConverter<T>::convert(data);
}

// ------------------------
// for hoNDImage
// ------------------------

template <typename T, unsigned int D>
class matlab_hoImage_header
{
public:

    typedef hoNDImage<T, D> ImageType;

    typedef typename ImageType::value_type value_type;
    typedef typename ImageType::coord_type coord_type;
    typedef typename ImageType::a_axis_type a_axis_type;
    typedef typename ImageType::axis_type axis_type;

    coord_type pixelSize_[D];
    coord_type origin_[D];
    hoNDPoint<coord_type, D> axis_[D];

    matlab_hoImage_header();
    matlab_hoImage_header(const ImageType& im);
    virtual ~matlab_hoImage_header();

    /// for the axis, it will be a D*D rotation matrix
    /// every column is a oritentation vector for a dimension
    virtual void toMatlab(mxArray*& header);
    virtual void fromMatlab(const mxArray* header);

protected:

    // the header field names
    std::vector<std::string> header_fields_;

    void set_header_fields()
    {
        size_t num = 3; // origin, pixelSize, axis
        header_fields_.resize(3);
        header_fields_[0] = "origin";
        header_fields_[1] = "pixelSize";
        header_fields_[2] = "axis";
    }
};

template <typename T, unsigned int D>
matlab_hoImage_header<T, D>::matlab_hoImage_header()
{
    unsigned int ii;
    for (ii=0;ii<D; ii++)
    {
        pixelSize_[ii] = 1;
        origin_[ii] = 0;
        axis_[ii].fill(0);
        axis_[ii][ii] = coord_type(1.0);
    }

    this->set_header_fields();
}

template <typename T, unsigned int D>
matlab_hoImage_header<T, D>::matlab_hoImage_header(const ImageType& im)
{
    std::vector<coord_type> pixelSize;
    im.get_pixel_size(pixelSize);

    std::vector<coord_type> origin;
    im.get_origin(origin);

    axis_type axis;
    im.get_axis(axis);

    unsigned int ii;
    for (ii=0;ii<D; ii++)
    {
        pixelSize_[ii] = pixelSize[ii];
        origin_[ii] = origin[ii];
        axis_[ii] = axis[ii];
    }

    this->set_header_fields();
}

template <typename T, unsigned int D>
matlab_hoImage_header<T, D>::~matlab_hoImage_header()
{

}

template <typename T, unsigned int D>
void matlab_hoImage_header<T, D>::toMatlab(mxArray*& header)
{
    try
    {
        unsigned int ii, jj;

        mwSize num[2] = {1, 1};

        std::vector<const char*> ptr(header_fields_.size());
        for( ii=0; ii< header_fields_.size(); ii++)
        {
            ptr[ii] = header_fields_[ii].c_str();
        }
        header = mxCreateStructArray(2, num, (int)header_fields_.size(), const_cast<const char**>(&ptr[0]));

        mwSize dims[1];
        dims[0] = D;

        mxArray* aMx = mxCreateNumericArray(1, dims, mxSINGLE_CLASS, mxREAL);
        float* pr = static_cast<float*>(mxGetData(aMx));
        for ( ii=0; ii<D; ii++ )
        {
            pr[ii] = origin_[ii];
        }

        mxSetField(header, 0, header_fields_[0].c_str(), aMx);

        aMx = mxCreateNumericArray(1, dims, mxSINGLE_CLASS, mxREAL);
        pr = static_cast<float*>(mxGetData(aMx));
        for ( ii=0; ii<D; ii++ )
        {
            pr[ii] = pixelSize_[ii];
        }

        mxSetField(header, 0, header_fields_[1].c_str(), aMx);

        mwSize dimsAxis[2];
        dimsAxis[0] = D;
        dimsAxis[1] = D;

        aMx = mxCreateNumericMatrix(D, D, mxSINGLE_CLASS, mxREAL);
        pr = static_cast<float*>(mxGetData(aMx));
        for ( jj=0; jj<D; jj++ )
        {
            for ( ii=0; ii<D; ii++ )
            {
                pr[jj + ii*D] = axis_[jj][ii]; // stored in column-wise
            }
        }

        mxSetField(header, 0, header_fields_[2].c_str(), aMx);
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in matlab_hoImage_header<T, D>::toMatlab(mxArray*& header) ... ");
    }
}

template <typename T, unsigned int D>
void matlab_hoImage_header<T, D>::fromMatlab(const mxArray* header)
{
    try
    {
        GADGET_CHECK_THROW(mxIsStruct(header));

        size_t ii, jj;

        mxArray* aMx = mxGetField(header, 0, header_fields_[0].c_str());
        size_t N = mxGetNumberOfElements(aMx);

        size_t minDN = ( (D<N) ? D : N );

        if ( mxIsSingle(aMx) )
        {
            float* pr = static_cast<float*>(mxGetData(aMx));

            for ( ii=0; ii<minDN; ii++ )
            {
                origin_[ii] = (coord_type)pr[ii];
            }
        }
        else
        {
            double* pr = static_cast<double*>(mxGetData(aMx));

            for ( ii=0; ii<minDN; ii++ )
            {
                origin_[ii] = (coord_type)pr[ii];
            }
        }

        aMx = mxGetField(header, 0, header_fields_[1].c_str());
        N = mxGetNumberOfElements(aMx);

        if ( mxIsSingle(aMx) )
        {
            float* pr = static_cast<float*>(mxGetData(aMx));

            for ( ii=0; ii<minDN; ii++ )
            {
                pixelSize_[ii] = (coord_type)pr[ii];
            }
        }
        else
        {
            double* pr = static_cast<double*>(mxGetData(aMx));

            for ( ii=0; ii<minDN; ii++ )
            {
                pixelSize_[ii] = (coord_type)pr[ii];
            }
        }

        aMx = mxGetField(header, 0, header_fields_[2].c_str());

        if ( mxIsSingle(aMx) )
        {
            float* pr = static_cast<float*>(mxGetData(aMx));

            for ( jj=0; jj<minDN; jj++ )
            {
                for ( ii=0; ii<minDN; ii++ )
                {
                    axis_[jj][ii] = (coord_type)pr[jj + ii*D];
                }
            }
        }
        else
        {
            double* pr = static_cast<double*>(mxGetData(aMx));

            for ( jj=0; jj<minDN; jj++ )
            {
                for ( ii=0; ii<minDN; ii++ )
                {
                    axis_[jj][ii] = (coord_type)pr[jj + ii*D];
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Error happened in matlab_hoImage_header<T, D>::fromMatlab(const mxArray* header) ... ");
    }
}

// ------------------------


template<class T, unsigned int D> 
mxArray* hoNDImageToMatlab(const hoNDImage<T, D>* a, mxArray*& h )
{
    std::vector<size_t> dim(D);
    a->get_dimensions(dim);

    hoNDArray<T> buf(dim, const_cast<T*>(a->get_data_ptr()), false);
    mxArray* m = MatlabConverter<T>::convert(&buf);

    matlab_hoImage_header<T, D> header(*a);
    header.toMatlab(h);

    return m;
}

template<class T, unsigned int D> 
void MatlabToHoNDImage(const mxArray* m, const mxArray* h, hoNDImage<T, D>& a)
{
    mwSize ndim = mxGetNumberOfDimensions(m);
    GADGET_CHECK_THROW( ndim == D );

    hoNDArray<T> buf = MatlabConverter<T>::convert(const_cast<mxArray*>(m));
    a.from_NDArray(buf);

    matlab_hoImage_header<T, D> header;
    header.fromMatlab(h);

    unsigned int ii;
    for ( ii=0; ii<D; ii++ )
    {
        a.set_pixel_size(ii, header.pixelSize_[ii]);
        a.set_origin(ii, header.origin_[ii]);
        a.set_axis(ii, header.axis_[ii]);
    }
}

// ------------------------
// IsmrmrdDataBuffered
// ------------------------

mxArray* BufferToMatlabStruct(IsmrmrdDataBuffered* buffer, bool omitData){

    const char * field_names[] = {"data","trajectory","headers","samplingdescription"};
	mwSize one = 1;
	auto mxstruct = mxCreateStructArray(1,&one,4,field_names);


	if (!mxstruct) throw std::runtime_error("Failed to allocate Matlab struct");

    if(omitData) {
        //auto mxdata = hoNDArrayToMatlab(&buffer->data_);
        
        using namespace std;
        cout << "Compressing data... ";
        clock_t b = clock();
        
        std::complex<float>* raw_data = buffer->data_.get_data_ptr();
        

        size_t nelem  = buffer->data_.get_number_of_elements();
        size_t h_nelem  = buffer->headers_.get_number_of_elements();
        
        size_t nRO    = buffer->data_.get_size(0);
        size_t nPE    = buffer->data_.get_size(1);
        size_t n3D    = buffer->data_.get_size(2);
        size_t nCH    = buffer->data_.get_size(3);
        size_t N      = buffer->data_.get_size(4);
        size_t S      = buffer->data_.get_size(5);
        size_t wtf      = buffer->data_.get_size(6);
        
        // count the number of non-nul RO lines in this buffer (there's probably a more elegant built-in method)
        size_t RO_counter = 0;
        for (size_t l = 0; l < h_nelem; ++l)
            if((bool) buffer->headers_[l].read_dir[2])
                RO_counter += nCH;
        /*
        
        
        RO_counter = 0;
        for (size_t l = 0; l < nelem; l += nRO)
            if(real(raw_data[l]) != 0.0f)
                ++RO_counter;
        std::cout << "RO_counter: " << RO_counter << std::endl;
        
        
        
        for(size_t l=0; l<buffer->headers_.get_number_of_elements(); ++l)
        {
            if(l%64==0)
                cout << "\n";
            
            cout << buffer->headers_[l].read_dir[2];

        }
        */
        /*
        std::cout << "N elem: " << buffer->data_.get_number_of_elements() << std::endl;
        std::cout << "N phase: " << buffer->data_.get_number_of_elements()/buffer->data_.get_size(0) << std::endl;
        std::cout << "RO_counter: " << RO_counter << std::endl;
        std::cout << "data dims: " << buffer->data_.get_size(0) << "," <<
                                      buffer->data_.get_size(1) << "," <<
                                      buffer->data_.get_size(2) << "," <<
                                      buffer->data_.get_size(3) << "," <<
                                      buffer->data_.get_size(4)  << std::endl;
        */
        
        
        // create the packet. A copy of the data is being done here,
        // which overall increase the RAM usage if packets are needed.
        // There may be a more efficient way to do this.
        size_t packet_n_elem = RO_counter * buffer->data_.get_size(0);
        size_t packet_ndim = 2;//buffer->data_.get_number_of_dimensions();
        mwSize* packet_dims = new mwSize[packet_ndim];
        
        packet_dims[0] = buffer->data_.get_size(0);
        packet_dims[1] = RO_counter;

        float* real_data = (float*) mxCalloc(packet_n_elem, sizeof(float));
        float* imag_data = (float*) mxCalloc(packet_n_elem, sizeof(float));

        /*
        size_t counter = 0;
        for (size_t l = 0; l < nelem; l += nRO ){
            if(real(raw_data[l]) != 0.0f) { // need to find a more proper test, e.g. using idx as look up table for getting directly the acquired RO line
                for (size_t j = 0; j < nRO; j++){

                    real_data[counter] = real(raw_data[l + j]);
                    imag_data[counter] = imag(raw_data[l + j]);
                    ++counter;
                }
            }
        }
        */
        
        /*
        size_t counter = 0;
        for (size_t l = 0; l < buffer->headers_.get_number_of_elements(); ++l) {
            
            if((bool) buffer->headers_[l].read_dir[2])
            {
                //for (size_t ch = 0; ch < nCH; ch++){
                    for (size_t j = 0; j < nRO*nCH; j++){

                        real_data[counter] = real(raw_data[l*nRO*nCH + j]);
                        imag_data[counter] = imag(raw_data[l*nRO*nCH + j]);
                        ++counter;
                    }
                //}
            }
        }
        */
        size_t counter = 0;
        size_t h_idx = 0;
        for (size_t ch = 0; ch < nCH; ++ch){
            for (size_t l = 0; l < h_nelem; ++l) {
                if((bool) buffer->headers_[l].read_dir[2])
                {
                    for (size_t r = 0; r < nRO; ++r){
                            h_idx = ch*nRO*nPE*n3D + l*nRO + r;
                            real_data[counter] = real(raw_data[h_idx]);
                            imag_data[counter] = imag(raw_data[h_idx]);
                            ++counter;
                    }
                }
            }
        }
        
        cout << "done (" << (double) (clock() - b)/CLOCKS_PER_SEC << ")\n";
        
        

        auto mxdata =  mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);
        mxSetDimensions(mxdata, packet_dims, packet_ndim);
        mxSetData      (mxdata, real_data);
//        mxSetImagData  (mxdata, imag_data);
        
        
        
        mxSetField(mxstruct,0,"data",mxdata);
        
    }
    else // don't omit data
    {
	auto mxdata = hoNDArrayToMatlab(&buffer->data_);
	mxSetField(mxstruct,0,"data",mxdata);
    }
    
	//Add trajectory if available
	if (buffer->trajectory_){
		auto & trajectory = *buffer->trajectory_;
		int traj_fieldnumber = mxAddField(mxstruct,"trajectory");
		auto mxtraj = hoNDArrayToMatlab(&trajectory);
		mxSetFieldByNumber(mxstruct,0,traj_fieldnumber,mxtraj);
	}

	//Add headers
	std::cout << "Adding headers...";
	mwSize num_headers = buffer->headers_.get_number_of_elements();
	auto mxheaders = mxCreateNumericMatrix(sizeof(ISMRMRD::AcquisitionHeader),num_headers,mxUINT8_CLASS,mxREAL);
	memcpy(mxGetData(mxheaders),buffer->headers_.get_data_ptr(),sizeof(ISMRMRD::AcquisitionHeader)*num_headers);
	mxSetField(mxstruct,0,"headers",mxheaders);

	auto samplingdescription = samplingdescriptionToMatlabStruct(&buffer->sampling_);
	mxSetField(mxstruct,0,"samplingdescription",samplingdescription);
    std::cout << " done." << std::endl;
	return mxstruct;
    
    /*
	const char * field_names[] = {"data","trajectory","headers","samplingdescription"};
	mwSize one = 1;
	auto mxstruct = mxCreateStructArray(1,&one,4,field_names);


	if (!mxstruct) throw std::runtime_error("Failed to allocate Matlab struct");

    if(!omitData) {
        auto mxdata = hoNDArrayToMatlab(&buffer->data_);
        mxSetField(mxstruct,0,"data",mxdata);
    }
    
	//Add trajectory if available
	if (buffer->trajectory_){
		auto & trajectory = *buffer->trajectory_;
		int traj_fieldnumber = mxAddField(mxstruct,"trajectory");
		auto mxtraj = hoNDArrayToMatlab(&trajectory);
		mxSetFieldByNumber(mxstruct,0,traj_fieldnumber,mxtraj);
	}

	//Add headers
	std::cout << "Adding headers...";
	mwSize num_headers = buffer->headers_.get_number_of_elements();
	auto mxheaders = mxCreateNumericMatrix(sizeof(ISMRMRD::AcquisitionHeader),num_headers,mxUINT8_CLASS,mxREAL);
	memcpy(mxGetData(mxheaders),buffer->headers_.get_data_ptr(),sizeof(ISMRMRD::AcquisitionHeader)*num_headers);
	mxSetField(mxstruct,0,"headers",mxheaders);

	auto samplingdescription = samplingdescriptionToMatlabStruct(&buffer->sampling_);
	mxSetField(mxstruct,0,"samplingdescription",samplingdescription);
    std::cout << " done." << std::endl;
	return mxstruct;
    */

}


mxArray* GetSplitReconData(IsmrmrdDataBuffered* buffer, size_t index_begin, size_t index_end) {
    
    // create the packet. A copy of the data is being done here,
    // which overall increase the RAM usage if packets are needed.
    // There may be a more efficient way to do this.
    size_t packet_n_elem = (index_end-index_begin+1) * buffer->data_.get_number_of_elements()/buffer->data_.get_size(0);

    size_t packet_ndim = buffer->data_.get_number_of_dimensions();
    mwSize* packet_dims = new mwSize[packet_ndim];
    packet_dims[0] = index_end-index_begin+1;
    for (size_t j = 1; j < packet_ndim; j++)
        packet_dims[j] = buffer->data_.get_size(j);

    float* real_data = (float*) mxCalloc(packet_n_elem, sizeof(float));
    float* imag_data = (float*) mxCalloc(packet_n_elem, sizeof(float));

    std::complex<float>* raw_data = buffer->data_.get_data_ptr();

    // It appears that the data is not stored as I would have expected it.
    // index 1,2,3,... actually follow the RO dimension, even though RO
    // is the first dimension of the data. In MATLAB this is the other way around.
    /*
    size_t start = index_begin*dim_1_n_elem;
    for (size_t j = 0; j < packet_n_elem; j++){
        real_data[j] = real(raw_data[start + j]);
        imag_data[j] = imag(raw_data[start + j]);
    }
    GDEBUG("Index: start %lu, index_end: %lu\n", start, start + packet_n_elem - 1);
     */

    size_t counter = 0;
    for (size_t l = 0; l < buffer->data_.get_number_of_elements(); l += buffer->data_.get_size(0) ){
        for (size_t j = 0; j < index_end-index_begin+1; j++){

            real_data[counter] = real(raw_data[index_begin + l + j]);
            imag_data[counter] = imag(raw_data[index_begin + l + j]);
            ++counter;
        }
    }

    auto mxdata =  mxCreateNumericMatrix(0, 0, mxSINGLE_CLASS, mxCOMPLEX);
    mxSetDimensions(mxdata, packet_dims, packet_ndim);
    mxSetData      (mxdata, real_data);
//    mxSetImagData  (mxdata, imag_data);
    
    return mxdata;
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


IsmrmrdDataBuffered MatlabStructToBuffer(mxArray* mxstruct){
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


mxArray* samplingdescriptionToMatlabStruct(SamplingDescription* samp){

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

// ------------------------
// std vector
// ------------------------

template<class T> 
mxArray* StdVecToMatlab(const std::vector<T>* a)
{
    hoNDArray<T> t(a->size(), const_cast<T*>(&(*a)[0]));
    return MatlabConverter<T>::convert(&t);
}

template<class T> 
void MatlabToStdVec(const mxArray* m, std::vector<T>& a)
{
    hoNDArray<T> t = MatlabConverter<T>::convert(const_cast<mxArray*>(m));
    a.resize(t.get_number_of_elements());
    memcpy(&(a[0]), t.get_data_ptr(), t.get_number_of_bytes());
}

// ------------------------
// std string
// ------------------------

mxArray* StdStringToMatlab(const std::string* a)
{
    return mxCreateString(a->c_str());
}

std::string MatlabToStdString(const mxArray* a)
{
    mwSize N = mxGetNumberOfElements(a) + 1;

    std::vector<char> buf(N, '\0');
    mxGetString(a, &buf[0], N);
    return std::string(&buf[0]);
}

// ------------------------
// Instantiation
// ------------------------

template EXPORTMATLAB mxArray* hoNDArrayToMatlab<float>(hoNDArray<float> *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<double>(hoNDArray<double> *);

template EXPORTMATLAB mxArray* hoNDArrayToMatlab<size_t>(hoNDArray<size_t> *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<float_complext>(hoNDArray<float_complext> *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<double_complext>(hoNDArray<double_complext> *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<std::complex<double> >(hoNDArray<std::complex<double> > *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<std::complex<float> >(hoNDArray<std::complex<float> > *);


template EXPORTMATLAB hoNDArray<float> MatlabToHoNDArray<float>(mxArray *);
template EXPORTMATLAB hoNDArray<double> MatlabToHoNDArray<double>(mxArray *);
template EXPORTMATLAB hoNDArray<size_t> MatlabToHoNDArray<size_t>(mxArray *);
template EXPORTMATLAB hoNDArray<float_complext> MatlabToHoNDArray<float_complext>(mxArray *);
template EXPORTMATLAB hoNDArray<double_complext> MatlabToHoNDArray<double_complext>(mxArray *);

template EXPORTMATLAB hoNDArray<std::complex<double> > MatlabToHoNDArray<std::complex<double> >(mxArray *);
template EXPORTMATLAB hoNDArray<std::complex<float> > MatlabToHoNDArray<std::complex<float> >(mxArray *);

template EXPORTMATLAB void MatlabToHoNDArray(mxArray * data, hoNDArray<float>& a);
template EXPORTMATLAB void MatlabToHoNDArray(mxArray * data, hoNDArray<double>& a);
template EXPORTMATLAB void MatlabToHoNDArray(mxArray * data, hoNDArray<float_complext>& a);
template EXPORTMATLAB void MatlabToHoNDArray(mxArray * data, hoNDArray<double_complext>& a);
template EXPORTMATLAB void MatlabToHoNDArray(mxArray * data, hoNDArray<std::complex<float> >& a);
template EXPORTMATLAB void MatlabToHoNDArray(mxArray * data, hoNDArray<std::complex<double> >& a);

template EXPORTMATLAB hoNDArray<vector_td<float,1> > MatlabToHoNDArray<vector_td<float,1> >(mxArray *);
template EXPORTMATLAB hoNDArray<vector_td<float,2> > MatlabToHoNDArray<vector_td<float,2> >(mxArray *);
template EXPORTMATLAB hoNDArray<vector_td<float,3> > MatlabToHoNDArray<vector_td<float,3> >(mxArray *);
template EXPORTMATLAB hoNDArray<vector_td<float,4> > MatlabToHoNDArray<vector_td<float,4> >(mxArray *);

template EXPORTMATLAB mxArray* hoNDArrayToMatlab<vector_td<float,1> >(hoNDArray<vector_td<float,1> > *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<vector_td<float,2> >(hoNDArray<vector_td<float,2> > *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<vector_td<float,3> >(hoNDArray<vector_td<float,3> > *);
template EXPORTMATLAB mxArray* hoNDArrayToMatlab<vector_td<float,4> >(hoNDArray<vector_td<float,4> > *);


template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<short, 1>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<short, 2>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<short, 3>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<short, 4>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<short, 5>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<short, 6>* a, mxArray*& header );

template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float, 1>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float, 2>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float, 3>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float, 4>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float, 5>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float, 6>* a, mxArray*& header );

template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double, 1>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double, 2>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double, 3>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double, 4>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double, 5>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double, 6>* a, mxArray*& header );

template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float_complext, 1>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float_complext, 2>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float_complext, 3>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float_complext, 4>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float_complext, 5>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<float_complext, 6>* a, mxArray*& header );

template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double_complext, 1>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double_complext, 2>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double_complext, 3>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double_complext, 4>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double_complext, 5>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<double_complext, 6>* a, mxArray*& header );

template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<float>, 1>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<float>, 2>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<float>, 3>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<float>, 4>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<float>, 5>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<float>, 6>* a, mxArray*& header );

template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<double>, 1>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<double>, 2>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<double>, 3>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<double>, 4>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<double>, 5>* a, mxArray*& header );
template EXPORTMATLAB mxArray* hoNDImageToMatlab(const hoNDImage<std::complex<double>, 6>* a, mxArray*& header );

template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<short, 1>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<short, 2>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<short, 3>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<short, 4>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<short, 5>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<short, 6>& a);

template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float, 1>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float, 2>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float, 3>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float, 4>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float, 5>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float, 6>& a);

template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double, 1>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double, 2>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double, 3>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double, 4>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double, 5>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double, 6>& a);

template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float_complext, 1>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float_complext, 2>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float_complext, 3>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float_complext, 4>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float_complext, 5>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<float_complext, 6>& a);

template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double_complext, 1>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double_complext, 2>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double_complext, 3>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double_complext, 4>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double_complext, 5>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<double_complext, 6>& a);

template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<float>, 1>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<float>, 2>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<float>, 3>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<float>, 4>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<float>, 5>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<float>, 6>& a);

template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<double>, 1>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<double>, 2>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<double>, 3>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<double>, 4>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<double>, 5>& a);
template EXPORTMATLAB void MatlabToHoNDImage(const mxArray* m, const mxArray* header, hoNDImage<std::complex<double>, 6>& a);

template EXPORTMATLAB mxArray* StdVecToMatlab(const std::vector<float>* a);
template EXPORTMATLAB mxArray* StdVecToMatlab(const std::vector<double>* a);
template EXPORTMATLAB mxArray* StdVecToMatlab(const std::vector<float_complext>* a);
template EXPORTMATLAB mxArray* StdVecToMatlab(const std::vector<double_complext>* a);
template EXPORTMATLAB mxArray* StdVecToMatlab(const std::vector<std::complex<float> >* a);
template EXPORTMATLAB mxArray* StdVecToMatlab(const std::vector<std::complex<double> >* a);

template EXPORTMATLAB void MatlabToStdVec(const mxArray* m, std::vector<float>& a);
template EXPORTMATLAB void MatlabToStdVec(const mxArray* m, std::vector<double>& a);
template EXPORTMATLAB void MatlabToStdVec(const mxArray* m, std::vector<float_complext>& a);
template EXPORTMATLAB void MatlabToStdVec(const mxArray* m, std::vector<double_complext>& a);
template EXPORTMATLAB void MatlabToStdVec(const mxArray* m, std::vector<std::complex<float> >& a);
template EXPORTMATLAB void MatlabToStdVec(const mxArray* m, std::vector<std::complex<double> >& a);

}
