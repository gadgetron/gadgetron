/** \file hoNDArray.h
    \brief CPU-based N-dimensional array (data container)
*/

#pragma once

#include "NDArray.h"
#include "complext.h"
#include "vector_td.h"
#include "GadgetronCommon.h"
#include "SerializableObject.h"

#include "cpucore_export.h"
#include <string.h>
#include <float.h>
#include <boost/shared_ptr.hpp>
#include <stdexcept>
#ifdef USE_MKL
    #include "mkl.h"
#endif // USE_MKL

namespace Gadgetron{

  template <typename T> class hoNDArray : public NDArray<T>, public SerializableObject
  {
  public:

    typedef NDArray<T> BaseClass;

    hoNDArray();
    hoNDArray(std::vector<unsigned long long> *dimensions);

    hoNDArray(unsigned long long len);
    hoNDArray(unsigned long long sx, unsigned long long sy);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss);

    hoNDArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long len, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, T* data, bool delete_data_on_destruct = false);
    hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, T* data, bool delete_data_on_destruct = false);

    hoNDArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    hoNDArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~hoNDArray();

    // Copy constructor
    hoNDArray(const hoNDArray<T>& a);

    // Assignment operator
    hoNDArray& operator=(const hoNDArray& rhs);
    void fill(T value);

    virtual void create(std::vector<unsigned long long>& dimensions);

    virtual void create(std::vector<unsigned int> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);

    virtual void create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual void create(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    virtual void create(boost::shared_ptr<std::vector<unsigned long long>  > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual void create(unsigned long long len);
    virtual void create(unsigned long long sx, unsigned long long sy);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su);

    virtual void create(unsigned long long len, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, T* data, bool delete_data_on_destruct = false);
    virtual void create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su, T* data, bool delete_data_on_destruct = false);

    T* begin();
    const T* begin() const;

    T* end();
    const T* end() const;

    T& at( unsigned long long idx );
    const T& at( unsigned long long idx ) const;

    T operator[]( unsigned long long idx );

    T& operator()( unsigned long long idx );
    const T& operator()( unsigned long long idx ) const;

    T& operator()( const std::vector<unsigned long long>& ind );
    const T& operator()( const std::vector<unsigned long long>& ind ) const;

    template<typename T2> 
    bool copyFrom(const hoNDArray<T2>& aArray)
    {
        try
        {
            if ( !this->dimensions_equal(&aArray) )
            {
                this->create(aArray.get_dimensions());
            }

            for ( unsigned long long ii=0; ii<elements_; ii++ )
            {
                data_[ii] = static_cast<T>(aArray(ii));
            }
        }
        catch(...)
        {
            return false;
        }

        return true;
    }

    void get_sub_array(const std::vector<unsigned long long>& start, std::vector<unsigned long long>& size, hoNDArray<T>& out);

    virtual void print(std::ostream& os) const;
    virtual void printContent(std::ostream& os) const;

    virtual bool serialize(char*& buf, unsigned long long& len) const;
    virtual bool deserialize(char* buf, unsigned long long& len);

  protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    virtual void allocate_memory();
    virtual void deallocate_memory();

    // Generic allocator / deallocator
    //

    template<class X> void _allocate_memory( unsigned long long size, X** data )
    {
        *data = new (std::nothrow) X[size];
    }

    template<class X> void _deallocate_memory( X* data )
    {
        delete [] data;
    }

    // Overload these instances to avoid invoking the element class constructor/destructor
    //

    virtual void _allocate_memory( unsigned long long size, float** data );
    virtual void _deallocate_memory( float* data );

    virtual void _allocate_memory( unsigned long long size, double** data );
    virtual void _deallocate_memory( double* data );

    virtual void _allocate_memory( unsigned long long size, std::complex<float>** data );
    virtual void _deallocate_memory( std::complex<float>* data );

    virtual void _allocate_memory( unsigned long long size, std::complex<double>** data );
    virtual void _deallocate_memory( std::complex<double>* data );

    virtual void _allocate_memory( unsigned long long size, float_complext** data );
    virtual void _deallocate_memory( float_complext* data );

    virtual void _allocate_memory( unsigned long long size, double_complext** data );
    virtual void _deallocate_memory( double_complext* data );

    template<class TYPE, unsigned long long D> void _allocate_memory( unsigned long long size, vector_td<TYPE,D>** data )
    {
        *data = (vector_td<TYPE,D>*) malloc( size*sizeof(vector_td<TYPE,D>) );
    }

    template<class TYPE, unsigned long long D>  void _deallocate_memory( vector_td<TYPE,D>* data )
    {
        free( data );
    }
  };
}

#include "hoNDArray.hxx"
