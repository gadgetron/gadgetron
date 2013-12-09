/** \file hoNDArrayMemoryManaged.h
    \brief CPU-based N-dimensional array (data container) using the memory manager to allocate memory
    \author Hui Xue
*/

#pragma once

#include "hoNDArray.h"
#include "gtPlusMemoryManager.h"

namespace Gadgetron{

  template <typename T> class hoNDArrayMemoryManaged : public Gadgetron::hoNDArray<T>
  {
  public:

    typedef Gadgetron::hoNDArray<T> BaseClass;
    typedef boost::shared_ptr<Gadgetron::gtPlus::gtPlusMemoryManager> MemManagerType;

    hoNDArrayMemoryManaged();
    hoNDArrayMemoryManaged(MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(std::vector<unsigned long long> *dimensions, MemManagerType& mem_manager);

    explicit hoNDArrayMemoryManaged(unsigned long long len, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, MemManagerType& mem_manager);

    explicit hoNDArrayMemoryManaged(boost::shared_ptr< std::vector<unsigned long long> > dimensions, MemManagerType& mem_manager);

    virtual ~hoNDArrayMemoryManaged();

    // Copy constructor
    hoNDArrayMemoryManaged(const hoNDArrayMemoryManaged<T>& a);

    // construct from hoNDArray
    hoNDArrayMemoryManaged(const hoNDArray<T>& a, MemManagerType& mem_manager);

    // Assignment operator
    hoNDArrayMemoryManaged& operator=(const hoNDArrayMemoryManaged<T>& rhs);
    hoNDArrayMemoryManaged& operator=(const hoNDArray<T>& rhs);

    // set memory manager
    void setMemoryManager(MemManagerType& mem_manager);

    virtual void print(std::ostream& os) const;

  protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    MemManagerType mem_manager_;

    // Generic allocator / deallocator
    //

    template<class X> void _allocate_memory( unsigned long long size, X** data )
    {
        if ( mem_manager_ )
        {
            void* ptr = mem_manager_->allocate(sizeof(X)*size);
            *data = (X*)ptr;
            for ( unsigned long long ii=0; ii<size; ii++ )
            {
                X* pt = (X*)ptr + ii;
                data[ii] = new(pt) X();
            }
        }
        else
        {
            *data = new X[size];
        }
    }

    template<class X> void _deallocate_memory( X* data )
    {
        for ( unsigned long long ii=0; ii<this->element_; ii++ )
        {
            data[ii].~X();
        }

        if ( mem_manager_ )
        {
            mem_manager_->free( (void*)data );
        }
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
        if ( mem_manager_ )
        {
            void* ptr = mem_manager_->allocate( size*sizeof(vector_td<TYPE,D>) );
            *data = (vector_td<TYPE,D>*)ptr;
        }
        else
        {
            *data = (vector_td<TYPE,D>*) malloc( size*sizeof(vector_td<TYPE,D>) );
        }
    }

    template<class TYPE, unsigned long long D>  void _deallocate_memory( vector_td<TYPE,D>* data )
    {
        if ( mem_manager_ )
        {
            mem_manager_->free( (void*)data );
        }
        else
        {
            free ( (void*)data );
        }
    }
  };
}

#include "hoNDArrayMemoryManaged.hxx"
