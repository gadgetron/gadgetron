
//Include for Visual studio, 'cos reasons.
#define _USE_MATH_DEFINES
#include <cmath>

#include "hoNDFFT.h"
#include "hoMatrix.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_math.h"

namespace Gadgetron{

template<typename T> hoNDFFT<T>* hoNDFFT<T>::instance()
    																{
	if (!instance_) instance_ = new hoNDFFT<T>();
	return instance_;
    																}

template<class T> hoNDFFT<T>* hoNDFFT<T>::instance_ = NULL;


template<class T> void hoNDFFT<T>::fft_int(hoNDArray< ComplexType >* input, size_t dim_to_transform, int sign)	{
	if (sign != -1 && sign != 1) throw std::runtime_error("hoNDFFT::fft_int: illegal sign provided");
	if (dim_to_transform >= input->get_number_of_dimensions()) throw std::runtime_error("hoNDFFT::fft_int: ransform dimension larger than dimension of input array ");

	//Only works for even dimensions. Fall back to slow version
	if (input->get_size(dim_to_transform)%2 == 1){
	    throw std::runtime_error("fft_int only works for even input");
	}
	int stride     = 1;           //Distance between points in transform
	int dist       = 1;           //Distance between vectors
	int trafos     = 1;           //Transformations per chunk
	int chunks     = 1;           //Number of chunks
	int chunk_size = 1;           //Points per chunk
	int length     = 1;           //Length of each transform
	int total_dist = 1;

	typename fftw_types<T>::plan * fft_plan        = 0;


	//Set sizes
	length = (int)input->get_size(dim_to_transform);

	T scale = 1/std::sqrt((T)length);
	if (dim_to_transform != 0)
	{
		for (size_t i = 0; i < dim_to_transform; i++)
		{
			chunk_size *= (int)input->get_size(i);
		}
		stride = chunk_size;
		trafos = chunk_size;
		chunk_size *= length;

		for (size_t i = dim_to_transform+1; i < input->get_number_of_dimensions(); i++)
		{
			chunks *= (int)input->get_size(i);
		}
	}
	else
	{
		for (size_t i = 1; i < input->get_number_of_dimensions(); i++)
		{
			trafos *= (int)input->get_size(i);
		}
		chunk_size = trafos*length;

		dist = length;
	}


	total_dist = trafos*dist;

	//Flip frequencies to center image
    timeswitch(input,dim_to_transform);
//Grab address of data
	ComplexType* data_ptr = input->get_data_ptr();

	//Allocate storage and make plan
	{
            std::lock_guard<std::mutex> guard(mutex_);
            unsigned planner_flags = FFTW_ESTIMATE;
            fft_plan = fftw_plan_many_dft_(1,&length,trafos,data_ptr,&length,stride,dist,data_ptr,&length,stride,dist,sign,planner_flags);
            //fftw_print_plan_(fft_plan);
            if (fft_plan == NULL)
	    {
                throw std::runtime_error("hoNDFFT: failed to create fft plan");
            }
	}

#pragma omp parallel for
	for (int k = 0; k < chunks; k++)
		fftw_execute_dft_(fft_plan,data_ptr+k*chunk_size,data_ptr+k*chunk_size);

//Flip frequencies to center DC freq
    timeswitch(input,dim_to_transform);


	//clean up
	{
            std::lock_guard<std::mutex> guard(mutex_);
            if (fft_plan != 0)
	    {
                fftw_destroy_plan_(fft_plan);
            }
	}


	*input *= scale;
}
template<typename T>
inline size_t hoNDFFT<T>::fftshiftPivot(size_t x)
{
	return (size_t)(ceil(x*0.5));
}

template<typename T>
inline size_t hoNDFFT<T>::ifftshiftPivot(size_t x)
{
	return (size_t)(floor(x*0.5));
}

template<typename T>
inline void hoNDFFT<T>::fftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot)
{
	memcpy(r, a+pivot, sizeof(ComplexType)*(x-pivot));
	memcpy(r+x-pivot, a, sizeof(ComplexType)*pivot);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift1D(const ComplexType* a, ComplexType* r, size_t x, size_t pivot)
{
	return fftshift1D(a, r, x, pivot);
}

template<typename T>
void hoNDFFT<T>::fftshiftPivot1D(ComplexType* a, size_t x, size_t n, size_t pivot)
{
	long long counter;

#pragma omp parallel private(counter) shared(n, x, pivot, a) if ( n > 256 )
	{
		hoNDArray< ComplexType > aTmp(x);

#pragma omp for
		for ( counter=0; counter<(long long)n; counter++ )
		{
			fftshift1D(a+counter*x, aTmp.begin(), x, pivot);
			memcpy(a+counter*x, aTmp.begin(), sizeof(ComplexType)*x);
		}
	}
}

template<typename T>
void hoNDFFT<T>::fftshiftPivot1D(const ComplexType* a, ComplexType* r, size_t x, size_t n, size_t pivot)
{
	long long counter;

#pragma omp parallel for private(counter) shared(n, x, pivot, a, r) if ( n > 256 )
	for ( counter=0; counter<(long long)n; counter++ )
	{
		fftshift1D(a+counter*x, r+counter*x, x, pivot);
	}
}


template<typename T>
void hoNDFFT<T>::fftshift1D(hoNDArray< ComplexType >& a)
{
	size_t x = a.get_size(0);
	size_t pivot = fftshiftPivot(x);
	size_t numOfShifts = a.get_number_of_elements()/x;
	fftshiftPivot1D(a.begin(), x, numOfShifts, pivot);
}

template<typename T>
void hoNDFFT<T>::fftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r = a;
	}

	size_t x = a.get_size(0);
	size_t pivot = fftshiftPivot(x);
	size_t numOfShifts = a.get_number_of_elements()/x;
	fftshiftPivot1D(a.begin(), r.begin(), x, numOfShifts, pivot);
}

template<typename T>
void hoNDFFT<T>::ifftshift1D(hoNDArray< ComplexType >& a)
{
	size_t x = a.get_size(0);
	size_t pivot = ifftshiftPivot(x);
	size_t numOfShifts = a.get_number_of_elements()/x;

	fftshiftPivot1D(a.begin(), x, numOfShifts, pivot);
}

template<typename T>
void hoNDFFT<T>::ifftshift1D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r = a;
	}

	size_t x = a.get_size(0);
	size_t pivot = ifftshiftPivot(x);
	size_t numOfShifts = a.get_number_of_elements()/x;

	fftshiftPivot1D(a.begin(), r.begin(), x, numOfShifts, pivot);
}

template<typename T>
void hoNDFFT<T>::fftshiftPivot2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty)
{
	if (a==NULL || r == NULL) throw std::runtime_error("hoNDFFT::fftshiftPivot2D: void ptr provided");

	long long tt;

#pragma omp parallel for private(tt) shared(a, r, x, y, n, pivotx, pivoty) if (n>16)
	for ( tt=0; tt<(long long)n; tt++ )
	{
		const ComplexType* ac = a + tt*x*y;
		ComplexType* rc = r + tt*x*y;

		size_t ay, ry;

		for ( ay=pivoty; ay<y; ay++ )
		{
			ry = ay - pivoty;
			memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
			memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
		}

		for ( ay=0; ay<pivoty; ay++ )
		{
			ry = ay + y - pivoty;
			memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
			memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
		}
	}
}

template<typename T>
void hoNDFFT<T>::fftshiftPivot2D(ComplexType* a, size_t x, size_t y, size_t n, size_t pivotx, size_t pivoty)
{

	if (a==NULL ) throw std::runtime_error("hoNDFFT::fftshiftPivot2D: void ptr provided");

	long long tt;

#pragma omp parallel private(tt) shared(a, x, y, n, pivotx, pivoty) if (n>16)
	{
		hoNDArray< ComplexType > aTmp(x*y);
		ComplexType* rc = aTmp.begin();

#pragma omp for
		for ( tt=0; tt<(long long)n; tt++ )
		{
			ComplexType* ac = a + tt*x*y;

			size_t ay, ry;

			for ( ay=pivoty; ay<y; ay++ )
			{
				ry = ay - pivoty;
				memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
				memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
			}

			for ( ay=0; ay<pivoty; ay++ )
			{
				ry = ay + y - pivoty;
				memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
				memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
			}

			memcpy(ac, rc, sizeof(ComplexType)*x*y);
		}
	}
}

template<typename T>
inline void hoNDFFT<T>::fftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n)
{

	if (a==NULL || r == NULL ) throw std::runtime_error("hoNDFFT::fftshift2D: void ptr provided");

	size_t pivotx = fftshiftPivot(x);
	size_t pivoty = fftshiftPivot(y);

	fftshiftPivot2D(a, r, x, y, n, pivotx, pivoty);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift2D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t n)
{

	if (a==NULL || r == NULL ) throw std::runtime_error("hoNDFFT::ifftshift2D: void ptr provided");

	size_t pivotx = ifftshiftPivot(x);
	size_t pivoty = ifftshiftPivot(y);

	fftshiftPivot2D(a, r, x, y, n, pivotx, pivoty);
}

template<typename T>
inline void hoNDFFT<T>::fftshift2D(ComplexType* a, size_t x, size_t y, size_t n)
{

	if (a==NULL ) throw std::runtime_error("hoNDFFT::fftshift2D: void ptr provided");

	size_t pivotx = fftshiftPivot(x);
	size_t pivoty = fftshiftPivot(y);

	fftshiftPivot2D(a, x, y, n, pivotx, pivoty);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift2D(ComplexType* a, size_t x, size_t y, size_t n)
{

	size_t pivotx = ifftshiftPivot(x);
	if (a==NULL ) throw std::runtime_error("hoNDFFT::ifftshift2D: void ptr provided");
	size_t pivoty = ifftshiftPivot(y);

	fftshiftPivot2D(a, x, y, n, pivotx, pivoty);
}

template<typename T>
inline void hoNDFFT<T>::fftshift2D(hoNDArray< ComplexType >& a)
{
	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
	return fftshift2D(a.begin(), a.get_size(0), a.get_size(1), n);
}

template<typename T>
inline void hoNDFFT<T>::fftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r = a;
	}

	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
	return fftshift2D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), n);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift2D(hoNDArray< ComplexType >& a)
{
	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
	return ifftshift2D(a.begin(), a.get_size(0), a.get_size(1), n);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift2D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r = a;
	}

	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1));
	return ifftshift2D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), n);
}

template<typename T>
void hoNDFFT<T>::fftshiftPivot3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty,  size_t pivotz)
{

	if (a==NULL  || r == NULL) throw std::runtime_error("hoNDFFT::fftshift2D: void ptr provided");

	long long tt;

#pragma omp parallel for private(tt) shared(a, r, x, y, z, n, pivotx, pivoty, pivotz) if (n>16)
	for ( tt=0; tt<(long long)n; tt++ )
	{
		size_t ay, ry, az, rz;

		for ( az=pivotz; az<z; az++ )
		{
			rz = az - pivotz;

			const ComplexType* ac = a + tt*x*y*z + az*x*y;
			ComplexType* rc = r + tt*x*y*z + rz*x*y;

			for ( ay=pivoty; ay<y; ay++ )
			{
				ry = ay - pivoty;
				memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
				memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
			}

			for ( ay=0; ay<pivoty; ay++ )
			{
				ry = ay + y - pivoty;
				memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
				memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
			}
		}

		for ( az=0; az<pivotz; az++ )
		{
			rz = az + z - pivotz;

			const ComplexType* ac = a + tt*x*y*z + az*x*y;
			ComplexType* rc = r + tt*x*y*z + rz*x*y;

			for ( ay=pivoty; ay<y; ay++ )
			{
				ry = ay - pivoty;
				memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
				memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
			}

			for ( ay=0; ay<pivoty; ay++ )
			{
				ry = ay + y - pivoty;
				memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
				memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
			}
		}
	}
}

template<typename T>
void hoNDFFT<T>::fftshiftPivot3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n, size_t pivotx, size_t pivoty,  size_t pivotz)
{

	if (a==NULL  ) throw std::runtime_error("hoNDFFT::fftshiftPivot3D: void ptr provided");

	long long tt;

#pragma omp parallel private(tt) shared(a, x, y, z, n, pivotx, pivoty, pivotz) if (n>16)
	{
		hoNDArray< ComplexType > aTmp(x*y*z);

#pragma omp for
		for ( tt=0; tt<(long long)n; tt++ )
		{
			size_t ay, ry, az, rz;

			for ( az=pivotz; az<z; az++ )
			{
				rz = az - pivotz;

				const ComplexType* ac = a + tt*x*y*z + az*x*y;
				ComplexType* rc = aTmp.begin() + rz*x*y;

				for ( ay=pivoty; ay<y; ay++ )
				{
					ry = ay - pivoty;
					memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
					memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
				}

				for ( ay=0; ay<pivoty; ay++ )
				{
					ry = ay + y - pivoty;
					memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
					memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
				}
			}

			for ( az=0; az<pivotz; az++ )
			{
				rz = az + z - pivotz;

				const ComplexType* ac = a + tt*x*y*z + az*x*y;
				ComplexType* rc = aTmp.begin() + rz*x*y;

				for ( ay=pivoty; ay<y; ay++ )
				{
					ry = ay - pivoty;
					memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
					memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
				}

				for ( ay=0; ay<pivoty; ay++ )
				{
					ry = ay + y - pivoty;
					memcpy(rc+ry*x, ac+ay*x+pivotx, sizeof(ComplexType)*(x-pivotx));
					memcpy(rc+ry*x+x-pivotx, ac+ay*x, sizeof(ComplexType)*pivotx);
				}
			}

			memcpy(a+tt*x*y*z, aTmp.begin(), sizeof(ComplexType)*x*y*z);
		}
	}
}

template<typename T>
inline void hoNDFFT<T>::fftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n)
{

	if (a==NULL  || r==NULL ) throw std::runtime_error("hoNDFFT::fftshift3D: void ptr provided");

	size_t pivotx = fftshiftPivot(x);
	size_t pivoty = fftshiftPivot(y);
	size_t pivotz = fftshiftPivot(z);

	fftshiftPivot3D(a, r, x, y, z, n, pivotx, pivoty, pivotz);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift3D(const ComplexType* a, ComplexType* r, size_t x, size_t y, size_t z, size_t n)
{

	if (a==NULL  || r==NULL ) throw std::runtime_error("hoNDFFT::ifftshift3D: void ptr provided");

	size_t pivotx = ifftshiftPivot(x);
	size_t pivoty = ifftshiftPivot(y);
	size_t pivotz = ifftshiftPivot(z);

	fftshiftPivot3D(a, r, x, y, z, n, pivotx, pivoty, pivotz);
}

template<typename T>
inline void hoNDFFT<T>::fftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n)
{
	if (a==NULL   ) throw std::runtime_error("hoNDFFT::fftshift3D: void ptr provided");
	size_t pivotx = fftshiftPivot(x);
	size_t pivoty = fftshiftPivot(y);
	size_t pivotz = fftshiftPivot(z);
	fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift3D(ComplexType* a, size_t x, size_t y, size_t z, size_t n)
{
	if (a==NULL   ) throw std::runtime_error("hoNDFFT::ifftshift3D: void ptr provided");

	size_t pivotx = ifftshiftPivot(x);
	size_t pivoty = ifftshiftPivot(y);
	size_t pivotz = ifftshiftPivot(z);

	fftshiftPivot3D(a, x, y, z, n, pivotx, pivoty, pivotz);
}

template<typename T>
inline void hoNDFFT<T>::fftshift3D(hoNDArray< ComplexType >& a)
{
	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
	return fftshift3D(a.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
}

template<typename T>
inline void hoNDFFT<T>::fftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r = a;
	}

	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
	return fftshift3D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift3D(hoNDArray< ComplexType >& a)
{
	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
	return ifftshift3D(a.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
}

template<typename T>
inline void hoNDFFT<T>::ifftshift3D(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r = a;
	}

	size_t n = a.get_number_of_elements()/(a.get_size(0)*a.get_size(1)*a.get_size(2));
	return ifftshift3D(a.begin(), r.begin(), a.get_size(0), a.get_size(1), a.get_size(2), n);
}

// -----------------------------------------------------------------------------------------

template<typename T>
inline void hoNDFFT<T>::fft1(hoNDArray< ComplexType >& a)
{
	return fft1(a, true);
}

template<typename T>
inline void hoNDFFT<T>::ifft1(hoNDArray< ComplexType >& a)
{
	return fft1(a, false);
}

template<typename T>
inline void hoNDFFT<T>::fft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r.create(a.get_dimensions());
	}

	return fft1(const_cast<hoNDArray< ComplexType >&>(a), r, true);
}

template<typename T>
inline void hoNDFFT<T>::ifft1(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	if ( !r.dimensions_equal(&a) )
	{
		r.create(a.get_dimensions());
	}

	return fft1(const_cast<hoNDArray< ComplexType >&>(a), r, false);
}

template<typename T>
inline void hoNDFFT<T>::fft1c(hoNDArray< ComplexType >& a)
{
	ifftshift1D(a);
	fft1(a);
	fftshift1D(a);
}

template<typename T>
inline void hoNDFFT<T>::ifft1c(hoNDArray< ComplexType >& a)
{
	ifftshift1D(a);
	ifft1(a);
	fftshift1D(a);
}

template<typename T>
inline void hoNDFFT<T>::fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	ifftshift1D(a, r);
	fft1(r);
	fftshift1D(r);

}

template<typename T>
inline void hoNDFFT<T>::ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	ifftshift1D(a, r);
	ifft1(r);
	fftshift1D(r);
}

template<typename T>
inline void hoNDFFT<T>::fft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
{
	ifftshift1D(a, r);
	fft1(r, buf);
	fftshift1D(buf, r);
}

template<typename T>
inline void hoNDFFT<T>::ifft1c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
{
	ifftshift1D(a, r);
	ifft1(r, buf);
	fftshift1D(buf, r);
}

// -----------------------------------------------------------------------------------------

template<typename T>
inline void hoNDFFT<T>::fft2(hoNDArray< ComplexType >& a)
{
	fft2(a, true);
}

template<typename T>
inline void hoNDFFT<T>::ifft2(hoNDArray< ComplexType >& a)
{
	fft2(a, false);
}

template<typename T>
inline void hoNDFFT<T>::fft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	//r = a;
	//return fft2(r);
	if ( !r.dimensions_equal(&a) )
	{
		r.create(a.get_dimensions());
	}

	fft2(const_cast<hoNDArray< ComplexType >&>(a), r, true);
}

template<typename T>
inline void hoNDFFT<T>::ifft2(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	/*r = a;
        return ifft2(r);*/

	if ( !r.dimensions_equal(&a) )
	{
		r.create(a.get_dimensions());
	}

	fft2(const_cast<hoNDArray< ComplexType >&>(a), r, false);
}

template<typename T>
inline void hoNDFFT<T>::fft2c(hoNDArray< ComplexType >& a)
{
	ifftshift2D(a);
	fft2(a);
	fftshift2D(a);

}

template<typename T>
inline void hoNDFFT<T>::ifft2c(hoNDArray< ComplexType >& a)
{
	ifftshift2D(a);
	ifft2(a);
	fftshift2D(a);
}

template<typename T>
inline void hoNDFFT<T>::fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	ifftshift2D(a, r);
	fft2(r);
	fftshift2D(r);
}

template<typename T>
inline void hoNDFFT<T>::ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	ifftshift2D(a, r);
	ifft2(r);
	fftshift2D(r);
}

template<typename T>
inline void hoNDFFT<T>::fft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
{
	ifftshift2D(a, r);
	fft2(r, buf);
	fftshift2D(buf, r);
}

template<typename T>
inline void hoNDFFT<T>::ifft2c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
{
	ifftshift2D(a, r);
	ifft2(r, buf);
	fftshift2D(buf, r);
}

// -----------------------------------------------------------------------------------------

template<typename T>
inline void hoNDFFT<T>::fft3(hoNDArray< ComplexType >& a)
{
	fft3(a, true);
}

template<typename T>
inline void hoNDFFT<T>::ifft3(hoNDArray< ComplexType >& a)
{
	fft3(a, false);
}

template<typename T>
inline void hoNDFFT<T>::fft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	/*r = a;
        return fft3(r);*/
	if ( !r.dimensions_equal(&a) )
	{
		r.create(a.get_dimensions());
	}

	fft3(const_cast<hoNDArray< ComplexType >&>(a), r, true);
}

template<typename T>
inline void hoNDFFT<T>::ifft3(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	/*r = a;
        return ifft3(r);*/
	if ( !r.dimensions_equal(&a) )
	{
		r.create(a.get_dimensions());
	}

	fft3(const_cast<hoNDArray< ComplexType >&>(a), r, false);
}

template<typename T>
inline void hoNDFFT<T>::fft3c(hoNDArray< ComplexType >& a)
{
	ifftshift3D(a);
	fft3(a);
	fftshift3D(a);
}

template<typename T>
inline void hoNDFFT<T>::ifft3c(hoNDArray< ComplexType >& a)
{
	ifftshift3D(a);
	ifft3(a);
	fftshift3D(a);
}

template<typename T>
inline void hoNDFFT<T>::fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	ifftshift3D(a, r);
	fft3(r);
	fftshift3D(r);
}

template<typename T>
inline void hoNDFFT<T>::ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r)
{
	ifftshift3D(a, r);
	ifft3(r);
	fftshift3D(r);
}

template<typename T>
inline void hoNDFFT<T>::fft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
{
	ifftshift3D(a, r);
	fft3(r, buf);
	fftshift3D(buf, r);
}

template<typename T>
inline void hoNDFFT<T>::ifft3c(const hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, hoNDArray< ComplexType >& buf)
{
	ifftshift3D(a, r);
	ifft3(r, buf);
	fftshift3D(buf, r);
}

template<typename T>
void hoNDFFT<T>::fft1(hoNDArray< ComplexType >& a, bool forward)
{
	hoNDArray< ComplexType > res(a);
	fft1(res, a, forward);
}

template<typename T>
void hoNDFFT<T>::fft2(hoNDArray< ComplexType >& a, bool forward)
{
	hoNDArray< ComplexType > res(a);
	fft2(res, a, forward);
}

template<typename T>
void hoNDFFT<T>::fft3(hoNDArray< ComplexType >& a, bool forward)
{
	hoNDArray< ComplexType > res(a);
	fft3(res, a, forward);
}

template<typename T>
void hoNDFFT<T>::fft1(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
{
	r = a;

	int n0 = (int)a.get_size(0);
	T fftRatio = T(1.0/std::sqrt( T(n0) ));


	int num = (int)(a.get_number_of_elements()/n0);
	int num_thr = get_num_threads_fft1(n0, num);

	int n;


	typename fftw_types<T>::plan * p;

	if( num_thr > 1 )
	{
		{
                        std::lock_guard<std::mutex> guard(mutex_);
			if ( forward )
			{
				p = fftw_plan_dft_1d_(n0,a.get_data_ptr(),r.get_data_ptr(),FFTW_FORWARD,FFTW_ESTIMATE);
			}
			else
			{
				p = fftw_plan_dft_1d_(n0,a.get_data_ptr(),r.get_data_ptr(),FFTW_BACKWARD,FFTW_ESTIMATE);
			}
		}

#pragma omp parallel for private(n) shared(num, p, a, n0, r) num_threads(num_thr)
		for ( n=0; n<num; n++ )
		{
			fftw_execute_dft_(p, a.get_data_ptr()+n*n0,
					r.get_data_ptr()+n*n0);
		}

		{
                        std::lock_guard<std::mutex> guard(mutex_);
			fftw_destroy_plan_(p);
		}
	}
	else
	{
		// multiple fft interface
		{
                        std::lock_guard<std::mutex> guard(mutex_);
			if ( forward )
			{
				p = fftw_plan_many_dft_(1, &n0, num,
						a.get_data_ptr(), NULL,
						1, n0,
						r.get_data_ptr(), NULL,
						1, n0,
						FFTW_FORWARD, FFTW_ESTIMATE);
			}
			else
			{
				p = fftw_plan_many_dft_(1, &n0, num,
						a.get_data_ptr(),NULL,
						1, n0,
						r.get_data_ptr(), NULL,
						1, n0,
						FFTW_BACKWARD, FFTW_ESTIMATE);
			}
		}

		fftw_execute_(p);

		{
                        std::lock_guard<std::mutex> guard(mutex_);
			fftw_destroy_plan_(p);
		}
	}

	r *= fftRatio;
}

template<typename T>
void hoNDFFT<T>::fft2(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
{
	r = a;

	int n0 = (int)a.get_size(1);
	int n1 = (int)a.get_size(0);

	T fftRatio = T(1.0/std::sqrt( T(n0*n1) ));

	int num = (int)(a.get_number_of_elements()/(n0*n1));
	int num_thr = get_num_threads_fft2(n0, n1, num);

	int n;


	typename fftw_types<T>::plan * p;

	if ( num_thr > 1 )
	{
		{
                        std::lock_guard<std::mutex> guard(mutex_);
			p = fftw_plan_dft_2d_(n0, n1,
					a.begin(),
					r.begin(),
					forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);
		}

#pragma omp parallel for private(n) shared(num, p, a, n0, n1, r) num_threads(num_thr)
		for ( n=0; n<num; n++ )
		{
			fftw_execute_dft_(p, a.begin()+n*n0*n1,
					r.begin()+n*n0*n1);
		}

		{
                        std::lock_guard<std::mutex> guard(mutex_);
			fftw_destroy_plan_(p);
		}
	}
	else
	{
		// multiple fft interface

		int n[] = {n0, n1};
		int idist = n0*n1;
		int odist = n0*n1;

		{
                        std::lock_guard<std::mutex> guard(mutex_);
			p = fftw_plan_many_dft_(2, n, num,
					a.begin(), NULL,
					1, idist,
					r.begin(), NULL,
					1, odist,
					forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

		}

		fftw_execute_(p);

		{
                        std::lock_guard<std::mutex> guard(mutex_);
			fftw_destroy_plan_(p);
		}
	}

	r *= fftRatio;

}

template<typename T>
void hoNDFFT<T>::fft3(hoNDArray< ComplexType >& a, hoNDArray< ComplexType >& r, bool forward)
{
	r = a;

	int n2 = (int)a.get_size(0);
	int n1 = (int)a.get_size(1);
	int n0 = (int)a.get_size(2);

	T fftRatio = T(1.0/std::sqrt( T(n0*n1*n2) ));

	int num = (int)(a.get_number_of_elements()/(n0*n1*n2));
	int num_thr = get_num_threads_fft3(n0, n1, n2, num);

	long long n;

	typename fftw_types<T>::plan * p;

	{
                std::lock_guard<std::mutex> guard(mutex_);
		p = fftw_plan_dft_3d_(n0, n1, n2,
				a.get_data_ptr(),
				r.get_data_ptr(),
				forward ? FFTW_FORWARD : FFTW_BACKWARD, FFTW_ESTIMATE);

	}

#pragma omp parallel for private(n) shared(num, p, a, n0, n1, n2, r) if (num_thr > 1) num_threads(num_thr)
	for ( n=0; n<num; n++ )
	{
		fftw_execute_dft_(p, a.begin()+n*n0*n1*n2,
				r.begin()+n*n0*n1*n2);
	}

	{
                std::lock_guard<std::mutex> guard(mutex_);
		fftw_destroy_plan_(p);
	}

	r *= fftRatio;

}
// TODO: implement more optimized threading strategy
template<typename T>
inline int hoNDFFT<T>::get_num_threads_fft1(size_t n0, size_t num)
{
	if ( num_of_max_threads_ == 1 ) return 1;

	if ( n0*num>1024*128 )
	{
		return num_of_max_threads_;
	}
	else if ( n0*num>512*128 )
	{
		return ( (num_of_max_threads_>8) ? 8 : num_of_max_threads_);
	}
	else if ( n0*num>256*128 )
	{
		return ( (num_of_max_threads_>4) ? 4 : num_of_max_threads_);
	}
	else if ( n0*num>128*128 )
	{
		return 2;
	}

	return 1;
}

template<typename T>
inline int hoNDFFT<T>::get_num_threads_fft2(size_t n0, size_t n1, size_t num)
{
	if ( num_of_max_threads_ == 1 ) return 1;

	if ( n0*n1*num>128*128*64 )
	{
		return num_of_max_threads_;
	}
	else if ( n0*n1*num>128*128*32 )
	{
		return ( (num_of_max_threads_>8) ? 8 : num_of_max_threads_);
	}
	else if ( n0*n1*num>128*128*16 )
	{
		return ( (num_of_max_threads_>4) ? 4 : num_of_max_threads_);
	}
	else if ( n0*n1*num>128*128*8 )
	{
		return 2;
	}

	return 1;
}

template<typename T>
inline int hoNDFFT<T>::get_num_threads_fft3(size_t n0, size_t n1, size_t n2, size_t num)
{
	if ( num_of_max_threads_ == 1 ) return 1;

	if ( num >= num_of_max_threads_ )
	{
		return num_of_max_threads_;
	}

	return 1;
}


template<typename T>
void hoNDFFT<T>::timeswitch(hoNDArray<ComplexType>* inout, int dim_to_transform){
	size_t batchsize = 1;
	for (int i = 0; i < dim_to_transform; i++)
		batchsize *= inout->get_size(i);

	size_t dimsize = inout->get_size(dim_to_transform);


	ComplexType* data = inout->get_data_ptr();
	size_t num_elements = inout->get_number_of_elements();
#pragma omp parallel for
	for (long int k = 0; k < num_elements; k++){
		size_t index = (k/batchsize)%dimsize;
		if (index%2 == 1)
			data[k] *= -1;
	}
}

template<typename T>
void hoNDFFT<T>::phaseshift(hoNDArray<ComplexType>* inout, T phase, int dim_to_transform){
	size_t batchsize = 1;
	for (int i = 0; i < dim_to_transform; i++)
		batchsize *= inout->get_size(i);

	size_t dimsize = inout->get_size(dim_to_transform);


	ComplexType* data = inout->get_data_ptr();
	size_t num_elements = inout->get_number_of_elements();
#pragma omp parallel for
	for (long int k = 0; k < num_elements; k++){
		float index = (k/batchsize)%dimsize-dimsize/2.0;
			data[k] *= ComplexType(0,2.0*M_PI*index*phase);

	}
}
template<> int hoNDFFT<float>::fftw_import_wisdom_from_file_(FILE* file){
	return fftwf_import_wisdom_from_file(file);
}

template<> int hoNDFFT<double>::fftw_import_wisdom_from_file_(FILE* file){
	return fftw_import_wisdom_from_file(file);
}

template<> void hoNDFFT<float>::fftw_export_wisdom_to_file_(FILE* file){
	return fftwf_export_wisdom_to_file(file);
}


template<> void hoNDFFT<double>::fftw_export_wisdom_to_file_(FILE* file){
	return fftw_export_wisdom_to_file(file);
}

template<> void hoNDFFT<float>::fftw_cleanup_(){
	fftwf_cleanup();
}

template<> void hoNDFFT<double>::fftw_cleanup_(){
	fftw_cleanup();
}

template<> void* hoNDFFT<float>::fftw_malloc_(size_t n){
	return fftwf_malloc(n);
}

template<> void* hoNDFFT<double>::fftw_malloc_(size_t n){
	return fftw_malloc(n);
}

template<> void hoNDFFT<float>::fftw_free_(void* ptr){
	fftwf_free(ptr);
}

template<> void hoNDFFT<double>::fftw_free_(void* ptr){
	fftw_free(ptr);
}

template<> void hoNDFFT<float>::fftw_execute_dft_( fftwf_plan_s * ptr, ComplexType* in, ComplexType* out){
	fftwf_execute_dft(ptr,(fftwf_complex*)in, (fftwf_complex*) out);
}


template<> void hoNDFFT<double>::fftw_execute_dft_(fftw_plan_s * ptr, ComplexType* in, ComplexType* out){
	fftw_execute_dft(ptr, (fftw_complex*)in, (fftw_complex*)out);
}

template<> void hoNDFFT<float>::fftw_execute_( fftwf_plan_s * ptr){
	fftwf_execute(ptr);
}


template<> void hoNDFFT<double>::fftw_execute_(fftw_plan_s * ptr){
	fftw_execute(ptr);
}

template<> typename fftw_types<float>::plan * hoNDFFT<float>::fftw_plan_dft_1d_(int rank,  ComplexType* in, ComplexType * out,int sign, unsigned int flags){
	return fftwf_plan_dft_1d(rank,(fftwf_complex*)in,(fftwf_complex*)out,sign,flags);
}
template<> typename fftw_types<double>::plan * hoNDFFT<double>::fftw_plan_dft_1d_(int rank, ComplexType* in, ComplexType * out,int sign, unsigned int flags){
	return fftw_plan_dft_1d(rank,(fftw_complex*)in,(fftw_complex*)out,sign,flags);
}

template<> typename fftw_types<float>::plan * hoNDFFT<float>::fftw_plan_dft_2d_(int n0, int n1,  ComplexType* in, ComplexType * out,int sign, unsigned int flags){
	return fftwf_plan_dft_2d(n0,n1,(fftwf_complex*)in,(fftwf_complex*)out,sign,flags);
}
template<> typename fftw_types<double>::plan * hoNDFFT<double>::fftw_plan_dft_2d_(int n0, int n1, ComplexType* in, ComplexType * out,int sign, unsigned int flags){
	return fftw_plan_dft_2d(n0,n1,(fftw_complex*)in,(fftw_complex*)out,sign,flags);
}


template<> typename fftw_types<float>::plan * hoNDFFT<float>::fftw_plan_dft_3d_(int n0, int n1, int n2, ComplexType* in, ComplexType * out,int sign, unsigned int flags){
	return fftwf_plan_dft_3d(n0,n1,n2,(fftwf_complex*)in,(fftwf_complex*)out,sign,flags);
}
template<> typename fftw_types<double>::plan * hoNDFFT<double>::fftw_plan_dft_3d_(int n0, int n1, int n2, ComplexType* in, ComplexType * out,int sign, unsigned int flags){
	return fftw_plan_dft_3d(n0,n1,n2,(fftw_complex*)in,(fftw_complex*)out,sign,flags);
}

template<> typename fftw_types<float>::plan * hoNDFFT<float>::fftw_plan_many_dft_(int rank, const int *n, int howmany,
		ComplexType *in, const int *inembed,
		int istride, int idist,
		ComplexType *out, const int *onembed,
		int ostride, int odist,
		int sign, unsigned flags){
	return fftwf_plan_many_dft(rank,n,howmany,(fftwf_complex*)in,inembed,istride,idist,(fftwf_complex*)out,onembed,ostride,odist,sign,flags);
}

template<> typename fftw_types<double>::plan * hoNDFFT<double>::fftw_plan_many_dft_(int rank, const int *n, int howmany,
		ComplexType *in, const int *inembed,
		int istride, int idist,
		ComplexType *out, const int *onembed,
		int ostride, int odist,
		int sign, unsigned flags){
	return fftw_plan_many_dft(rank,n,howmany,(fftw_complex*)in,inembed,istride,idist,(fftw_complex*)out,onembed,ostride,odist,sign,flags);
}

template<> void hoNDFFT<float>::fftw_destroy_plan_( typename fftw_types<float>::plan * p ){
	fftwf_destroy_plan(p);
}

template<> void hoNDFFT<double>::fftw_destroy_plan_( typename fftw_types<double>::plan * p ){
	fftw_destroy_plan(p);
}

template<> void hoNDFFT<double>::fftw_print_plan_( typename fftw_types<double>::plan * p ){
	fftw_print_plan(p);
}


template<> void hoNDFFT<float>::fftw_print_plan_( typename fftw_types<float>::plan * p ){
	fftwf_print_plan(p);
}

// -----------------------------------------------------------------------------------------


//
// Instantiation
//

template class EXPORTCPUFFT hoNDFFT<float>;
template class EXPORTCPUFFT hoNDFFT<double>;
}
