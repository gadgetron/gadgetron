#ifndef _NDARRAY_HPP_
#define _NDARRAY_HPP_

#include "export.h"

#ifdef CUDA

namespace std{
#if defined (WIN32) || defined (_WIN32)
	template<class _Ty> class allocator;
	template<class _Ty, class _Ax = allocator<_Ty> > class vector;
#else
	template<typename _Tp> class allocator;
	template<typename _Tp, typename _Alloc = allocator<_Tp> > class vector;
#endif
}

#else // !CUDA

#include <vector>

#endif

namespace mr_recon
{
	template <class T> class DLLEXPORT NDArray
	{

	public:
		NDArray ();
		NDArray(int number_of_dimensions, int* dimensions);
		NDArray(std::vector<int> *dimensions);
		NDArray(int dim1);
		NDArray(int dim1, int dim2);
		NDArray(int dim1, int dim2, int dim3);
		NDArray(int dim1, int dim2, int dim3, int dim4);

		NDArray(const NDArray<T>& a);
		~NDArray();

		T& operator[] (int i);
		NDArray<T>& operator=(const NDArray<T> &a);
		NDArray<T> operator*(const NDArray<T> &a) const;
		NDArray<T> operator*(const T a) const;
		NDArray<T>& operator*=(const T s);

		NDArray<T> dot_multiply(const NDArray<T> &a);
		NDArray<T> dot_divide(const NDArray<T> &a);

		int get_number_of_dimensions();
		int get_size(unsigned int dimension);
		std::vector<int>* get_dimensions();
		unsigned long int get_number_of_elements();

		NDArray<T> get(std::vector<int> *lower_limits, std::vector<int> *upper_limits);
		NDArray<T>& set(const NDArray<T> &a, std::vector<int> *lower_limits, std::vector<int> *upper_limits);

		void flipdim (int dimension);

		T* get_data_ptr();

	private:
		std::vector<int> *m_dimensions;
		int ndim;
		unsigned long int elements;
		int* dims;
		T* data;

		void allocate_memory(std::vector<int> *dimensions);
		void deallocate_memory();
	};
}

#endif //NDARRAY_HPP
