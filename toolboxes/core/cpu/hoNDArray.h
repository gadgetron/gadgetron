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
#endif

namespace Gadgetron{

  template <typename T> class hoNDArray : public NDArray<T>, public SerializableObject
  {
  public:

    typedef NDArray<T> BaseClass;

    hoNDArray();
    hoNDArray(std::vector<size_t> *dimensions);

    hoNDArray(size_t len);
    hoNDArray(size_t sx, size_t sy);
    hoNDArray(size_t sx, size_t sy, size_t sz);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss);

    hoNDArray(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t len, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct = false);

    hoNDArray(boost::shared_ptr< std::vector<size_t> > dimensions);
    hoNDArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~hoNDArray();

    // Copy constructor
    hoNDArray(const hoNDArray<T>& a);

    // Assignment operator
    hoNDArray& operator=(const hoNDArray& rhs);
    void fill(T value);

    virtual void create(std::vector<size_t>& dimensions);

    virtual void create(std::vector<unsigned int> *dimensions);
    virtual void create(std::vector<size_t> *dimensions);

    virtual void create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct = false);
    virtual void create(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual void create(boost::shared_ptr< std::vector<size_t> > dimensions);
    virtual void create(boost::shared_ptr<std::vector<size_t>  > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual void create(size_t len);
    virtual void create(size_t sx, size_t sy);
    virtual void create(size_t sx, size_t sy, size_t sz);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su);

    virtual void create(size_t len, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct = false);
    virtual void create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su, T* data, bool delete_data_on_destruct = false);

    T* begin();
    const T* begin() const;

    T* end();
    const T* end() const;

    T& at( size_t idx );
    const T& at( size_t idx ) const;

    T operator[]( size_t idx );

    T& operator()( size_t idx );
    const T& operator()( size_t idx ) const;

    T& operator()( const std::vector<size_t>& ind );
    const T& operator()( const std::vector<size_t>& ind ) const;

    template<typename T2> 
    bool copyFrom(const hoNDArray<T2>& aArray)
    {
        try
        {
            if ( !this->dimensions_equal(&aArray) )
            {
                this->create(aArray.get_dimensions());
            }

            for ( size_t ii=0; ii<elements_; ii++ )
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

    void get_sub_array(const std::vector<size_t>& start, std::vector<size_t>& size, hoNDArray<T>& out);

    virtual void print(std::ostream& os) const;
    virtual void printContent(std::ostream& os) const;

    virtual bool serialize(char*& buf, size_t& len) const;
    virtual bool deserialize(char* buf, size_t& len);

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

    template<class X> void _allocate_memory( size_t size, X** data )
    {
        *data = new (std::nothrow) X[size];
    }

    template<class X> void _deallocate_memory( X* data )
    {
        delete [] data;
    }

    // Overload these instances to avoid invoking the element class constructor/destructor
    //

    virtual void _allocate_memory( size_t size, float** data );
    virtual void _deallocate_memory( float* data );

    virtual void _allocate_memory( size_t size, double** data );
    virtual void _deallocate_memory( double* data );

    virtual void _allocate_memory( size_t size, std::complex<float>** data );
    virtual void _deallocate_memory( std::complex<float>* data );

    virtual void _allocate_memory( size_t size, std::complex<double>** data );
    virtual void _deallocate_memory( std::complex<double>* data );

    virtual void _allocate_memory( size_t size, float_complext** data );
    virtual void _deallocate_memory( float_complext* data );

    virtual void _allocate_memory( size_t size, double_complext** data );
    virtual void _deallocate_memory( double_complext* data );

    template<class TYPE, unsigned int D> void _allocate_memory( size_t size, vector_td<TYPE,D>** data )
    {
        *data = (vector_td<TYPE,D>*) malloc( size*sizeof(vector_td<TYPE,D>) );
    }

    template<class TYPE, unsigned int D>  void _deallocate_memory( vector_td<TYPE,D>* data )
    {
        free( data );
    }
  };
}

#include "hoNDArray.hxx"
