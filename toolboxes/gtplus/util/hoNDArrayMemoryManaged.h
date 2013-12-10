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
    explicit hoNDArrayMemoryManaged(std::vector<size_t> *dimensions, MemManagerType& mem_manager);

    explicit hoNDArrayMemoryManaged(size_t len, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(size_t sx, size_t sy, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(size_t sx, size_t sy, size_t sz, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(size_t sx, size_t sy, size_t sz, size_t st, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, MemManagerType& mem_manager);
    explicit hoNDArrayMemoryManaged(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, MemManagerType& mem_manager);

    explicit hoNDArrayMemoryManaged(boost::shared_ptr< std::vector<size_t> > dimensions, MemManagerType& mem_manager);

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

    template<class X> void _allocate_memory( size_t size, X** data )
    {
        if ( mem_manager_ )
        {
            void* ptr = mem_manager_->allocate(sizeof(X)*size);
            *data = (X*)ptr;
            for ( size_t ii=0; ii<size; ii++ )
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
        for ( size_t ii=0; ii<this->element_; ii++ )
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

    template<class TYPE, unsigned int D>  void _deallocate_memory( vector_td<TYPE,D>* data )
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
