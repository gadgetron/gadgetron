#include "ndarray.hpp"

#include "types.hpp"
#include "types.hcu"

#include <string.h>
#include <vector>
#include <iostream>
#include <complex>

#include <string.h>

using namespace std;
using namespace mr_recon;
using namespace mr_recon_vectors;

template <class T> NDArray<T>::NDArray ()
{
	ndim = 0;
	dims = 0; data = 0;
	m_dimensions = new std::vector<int>;
}

template <class T> NDArray<T>::NDArray (vector<int> *dimensions)
{
	m_dimensions = new std::vector<int>;
	*m_dimensions = *dimensions;
	dims = 0; data = 0;
	allocate_memory(dimensions);
}

template <class T> NDArray<T>::NDArray (int number_of_dimensions, int* dimensions)
{
	m_dimensions = new std::vector<int>;
	for (int i = 0; i < number_of_dimensions; i++) 
		m_dimensions->push_back(dimensions[i]);
	dims = 0; data = 0;
	allocate_memory(m_dimensions);
}

template <class T> NDArray<T>::NDArray (int dim1)
{
	m_dimensions = new std::vector<int>;
	m_dimensions->push_back(dim1);
	dims = 0; data = 0;
	allocate_memory(m_dimensions);
}

template <class T> NDArray<T>::NDArray (int dim1, int dim2)
{
	m_dimensions = new std::vector<int>;
	m_dimensions->push_back(dim1);
	m_dimensions->push_back(dim2);
	dims = 0; data = 0;
	allocate_memory(m_dimensions);
}

template <class T> NDArray<T>::NDArray (int dim1, int dim2, int dim3)
{
	m_dimensions = new std::vector<int>;
	m_dimensions->push_back(dim1);
	m_dimensions->push_back(dim2);
	m_dimensions->push_back(dim3);
	dims = 0; data = 0;
	allocate_memory(m_dimensions);
}

template <class T> NDArray<T>::NDArray (int dim1, int dim2, int dim3, int dim4)
{
	m_dimensions = new std::vector<int>;
	m_dimensions->push_back(dim1);
	m_dimensions->push_back(dim2);
	m_dimensions->push_back(dim3);
	m_dimensions->push_back(dim4);
	dims = 0; data = 0;
	allocate_memory(m_dimensions);
}

template <class T> NDArray<T>::NDArray (const NDArray<T>& a)
{
	dims = 0; data = 0;
	m_dimensions = new std::vector<int>;
	*m_dimensions = *a.m_dimensions;
	allocate_memory(m_dimensions);
	memcpy( data, a.data, elements*sizeof(T));

//	for (unsigned int i = 0; i < elements; i++)
//      {
//	   data[i] = a.data[i];
//      }

}

template <class T> NDArray<T>::~NDArray ()
{
	deallocate_memory();
	delete( m_dimensions );
	m_dimensions = 0x0;
}

template<class T>  NDArray<T> NDArray<T>::get(vector<int> *lower_limits, vector<int> *upper_limits)
{
	if ((lower_limits->size() != upper_limits->size()) || (lower_limits->size() != m_dimensions->size()))
	{
		cout << "NDArray: Invalid Index Parameters. To many indices." << endl;
		cout << "NDArray: Dimensions in array: " << m_dimensions->size() << endl;
		cout << "NDArray: Dimensions in lower limits: " << lower_limits->size() << endl;
		cout << "NDArray: Dimensions in upper limits: " << lower_limits->size() << endl;

		return NDArray();
	}

	vector<int> indexes;
	vector<int> sizes;
	vector<int> sub_matrices;

	for (unsigned int i = 0; i < lower_limits->size(); i++)
	{
		if ((*lower_limits)[i] > (*upper_limits)[i])
		{
			cout << "NDArray: Invalid Index Parameters. Lower limits larger than upper limit." << endl;
			return NDArray();
		}
		if ((*lower_limits)[i] == -1 && (*upper_limits)[i] == -1)
		{
			(*lower_limits)[i] = 0; (*upper_limits)[i] = get_size(i)-1;
		}
		indexes.push_back((*lower_limits)[i]);
		sizes.push_back((*upper_limits)[i]-(*lower_limits)[i]+1);
		if (i == 0)
		{
			sub_matrices.push_back(1);
		}
		else
		{
			sub_matrices.push_back(sub_matrices[i-1]*(*m_dimensions)[i-1]);
		}
	}

	NDArray<T> out(&sizes);

	int index_i = 0;
	int index_o = 0;
	unsigned int p;
	for (p = 0; p < indexes.size(); p++) index_i += sub_matrices[p]*indexes[p];
	while (indexes[indexes.size()-1] <= (*upper_limits)[indexes.size()-1])
	{
		out[index_o] = data[index_i];

		indexes[0]++;
		index_i++;
		index_o++;
		unsigned int current_dim = 0;
		while(indexes[current_dim] > (*upper_limits)[current_dim] && current_dim < indexes.size())
		{
			indexes[current_dim] = (*lower_limits)[current_dim];
			current_dim++;
			if( current_dim == indexes.size() )
				break;
			else
				indexes[current_dim]++;
			index_i = 0;
			for (p = 0; p < indexes.size(); p++) index_i += sub_matrices[p]*indexes[p];
		}
		if (current_dim == indexes.size())
		{
			break;
		}
	}
	return out;
}


template<class T>  NDArray<T>& NDArray<T>::set(const NDArray<T>& a, vector<int> *lower_limits, vector<int> *upper_limits)
{
	if ((lower_limits->size() != upper_limits->size()) || (lower_limits->size() != m_dimensions->size()))
	{
		cout << "NDArray: Invalid Index Parameters. To many indices." << endl;
		cout << "NDArray: Dimensions in array: " << m_dimensions->size() << endl;
		cout << "NDArray: Dimensions in lower limits: " << lower_limits->size() << endl;
		cout << "NDArray: Dimensions in upper limits: " << lower_limits->size() << endl;
		cout << "NDArray: Invalid Index Parameters" << endl;
		return *this;
	}

	vector<int> indexes;
	vector<int> sizes;
	vector<int> sub_matrices;

	for (unsigned int i = 0; i < lower_limits->size(); i++)
	{
		if ((*lower_limits)[i] > (*upper_limits)[i] || (*upper_limits)[i] >= get_size(i) )
		{
			cout << "NDArray: Invalid Index Parameters" << endl;
			return *this;
		}
		if ((*lower_limits)[i] == -1 && (*upper_limits)[i] == -1)
		{
			(*lower_limits)[i] = 0; (*upper_limits)[i] = get_size(i)-1;
		}
		indexes.push_back((*lower_limits)[i]);
		sizes.push_back((*upper_limits)[i]-(*lower_limits)[i]+1);
		if (i == 0)
		{
			sub_matrices.push_back(1);
		}
		else
		{
			sub_matrices.push_back(sub_matrices[i-1]*(*m_dimensions)[i-1]);
		}
	}

	int index_i = 0;
	int index_o = 0;
	unsigned int p;
	for (p = 0; p < indexes.size(); p++) index_i += sub_matrices[p]*indexes[p];
	while (indexes[indexes.size()-1] <= (*upper_limits)[indexes.size()-1])
	{
		/* CHANGE THIS */
		data[index_i] = a.data[index_o];

		indexes[0]++;
		index_i++;
		index_o++;
		unsigned int current_dim = 0;
		while(indexes[current_dim] > (*upper_limits)[current_dim] && current_dim < indexes.size())
		{
			indexes[current_dim] = (*lower_limits)[current_dim];
			current_dim++;
			if( current_dim == indexes.size() )
				break;
			else
				indexes[current_dim]++;
			index_i = 0;
			for (p = 0; p < indexes.size(); p++) index_i += sub_matrices[p]*indexes[p];
		}
		if (current_dim == indexes.size())
		{
			break;
		}
	}
	return *this;
}

template <class T> int NDArray<T>::get_number_of_dimensions()
{
	return m_dimensions->size();
}

template <class T> int NDArray<T>::get_size(unsigned int dim)
{
	if (dim >= m_dimensions->size())
	{
		return 1;
	}
	return (*m_dimensions)[dim];
}


template < class T > vector<int>* NDArray<T>::get_dimensions()
{
	return m_dimensions;
}


template <class T> unsigned long int NDArray<T>::get_number_of_elements()
{
	return elements;
}

template <class T> T* NDArray<T>::get_data_ptr()
{
	return data;
}

template <class T> T& NDArray<T>::operator[] (int i)
{
	if (i < 0 || i >= (int)elements)
	{
		cout << "NDArray: Index out of range: i=" << i << ", elements=" << elements << endl;
		//return data[0];
	}
	return data[i];
}

template <class T> NDArray<T>& NDArray<T>::operator= (const NDArray<T>& a)
{
	if (this != &a)
	{
		if (*m_dimensions != *a.m_dimensions)
		{
			deallocate_memory();
			delete m_dimensions;
			m_dimensions = new vector<int>;
			*m_dimensions = *a.m_dimensions;
			allocate_memory(m_dimensions);
		}
		for (unsigned int i = 0; i < elements; i++) data[i] = a.data[i];
	}
	return *this;
}

template <class T> NDArray<T>& NDArray<T>::operator*= (const T s)
{
	for (unsigned int i = 0; i < elements; i++) data[i] *= s;
	return *this;
}

template <class T> NDArray<T> NDArray<T>::operator* (const NDArray<T>& a) const
{

	NDArray<T> out;

	/* This operation is only really well defined if inputs are matrices (or vectors) and dimensions match */
	if (m_dimensions->size() > 2 || a.m_dimensions->size() > 2)
	{
		cout << "Array multiplication only defined for vectors and matrices." << endl;
		return out;
	} 

	if ((*m_dimensions)[1] != (*(a.m_dimensions))[0])
	{
		cout << "Array dimensions do not match for array multiplication." << endl;
		return out;
	}

	int rows = (*m_dimensions)[0];
	int cols = (*(a.m_dimensions))[1];

	int m = (*m_dimensions)[1];

	out = NDArray<T>(rows,cols);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{

			for (int k = 0; k < m; k++)
			{
				out.data[j*rows + i] += data[k*rows + i]*a.data[j*m + k];
			} 
		} 
	}

	return out;
}

template <class T> NDArray<T> NDArray<T>::operator* (const T scalar) const
{

	NDArray<T> out = *this;

	for (int i=0; i<(int)elements; i++)
	{
		data[i] *= scalar;
	}

	return out;
}

template <class T>  NDArray<T> NDArray<T>::dot_multiply(const NDArray<T>& a)
{
	NDArray<T> out(m_dimensions);
	if (elements != a.elements)
	{
		cout << "Dot multiplication can only be performed on arrays of equal dimensions" << endl;
		return out;
	}

	for (unsigned int i = 0; i < elements; i++)
	{
		out.data[i] = data[i]*a.data[i];
	}

	return out;
}

template <class T>  NDArray<T> NDArray<T>::dot_divide(const NDArray<T>& a)
{
	NDArray<T> out(m_dimensions);
	if (elements != a.elements)
	{
		cout << "Dot multiplication can only be performed on arrays of equal dimensions" << endl;
		return out;
	}

	for (unsigned int i = 0; i < elements; i++)
	{
		out.data[i] = data[i]/a.data[i];
	}

	return out;
}

template <class T> void NDArray<T>::flipdim (int dimension)
{
	if (dimension > get_number_of_dimensions()-1)
	{
		cout << "Error: invalid dimension specified for flipping" << endl;
		return;
	}

	long elements_before = 1;
	for (int i = 0; i < dimension; i++)
	{
		elements_before *= get_size(i);
	}

	long elements_after = 1;
	for (int i = dimension+1; i < get_number_of_dimensions(); i++)
	{
		elements_after *= get_size(i);
	}

	long dimension_length = get_size(dimension);

	T* temp = NULL; 
	try
	{
		temp = new T[dimension_length];
	} 
	catch (...) 
	{
		cout << "Unable to allocate memory for flipdim" << endl;
		return;

	}

	int chunk_size = elements_before*dimension_length;
	for (int j = 0; j < elements_after; j++)
	{
		for (int i = 0; i < elements_before; i++)
		{
			for (int k = 0; k < dimension_length; k++)
			{
				temp[k] = data[j*chunk_size+i*dimension_length + k];
			}
			for (int k = 0; k < dimension_length; k++)
			{
				data[j*chunk_size+i*dimension_length + k] = temp[dimension_length-1-k];
			}
		}
	}

	if (temp)
	{
		delete [] temp;
	}
}

template <class T> void NDArray<T>::allocate_memory (vector<int> *dimensions)
{

	if (dims != 0) delete [] dims;
	if (data != 0) delete [] data;

	ndim = dimensions->size();

	try
	{
		dims = new int[ndim];
	}
	catch (bad_alloc&)
	{
		cout << "Error allocating memory for dimensions array" << endl;
	}

	elements = 1;
	for (int i = 0; i < ndim; i++)
	{
		dims[i] = (*dimensions)[i];
		elements *= dims[i];
	}

	try
	{
		data = new T[elements];
	}
	catch (bad_alloc&)
	{
		cout << "Error allocating memory for data" << endl;
	}
	for (unsigned int i = 0 ; i < elements; i++) data[i] = 0; 

}

template <class T> void NDArray<T>::deallocate_memory ()
{
	delete [] dims; dims = 0;
	delete [] data; data = 0;
	ndim = 0;
	elements = 0;
}


/* The declarations below are necessary to make sure that specific functions are included in the library */

template class DLLEXPORT NDArray< complex<double> >;
template class DLLEXPORT NDArray< complex<float> >;
template class DLLEXPORT NDArray< double >;
template class DLLEXPORT NDArray< float >;
template class DLLEXPORT NDArray< int >;
template class DLLEXPORT NDArray< unsigned int >;
template class DLLEXPORT NDArray< char >;
template class DLLEXPORT NDArray< float2 >;
template class DLLEXPORT NDArray< float3 >;
template class DLLEXPORT NDArray< float4 >;
template class DLLEXPORT NDArray< double2 >;
template class DLLEXPORT NDArray< double3 >;
template class DLLEXPORT NDArray< double4 >;
template class DLLEXPORT NDArray< int2 >;
template class DLLEXPORT NDArray< int3 >;
template class DLLEXPORT NDArray< int4 >;
template class DLLEXPORT NDArray< uint2 >;
template class DLLEXPORT NDArray< uint3 >;
template class DLLEXPORT NDArray< uint4 >;
template class DLLEXPORT NDArray< char2 >;
template class DLLEXPORT NDArray< char3 >;
template class DLLEXPORT NDArray< char4 >;
