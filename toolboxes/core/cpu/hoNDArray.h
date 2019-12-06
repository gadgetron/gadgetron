/** \file hoNDArray.h
    \brief CPU-based N-dimensional array (data container)
*/

#pragma once

#include "NDArray.h"
#include "complext.h"
#include "vector_td.h"

#include "cpucore_export.h"

#include <string.h>
#include <float.h>
#include <boost/shared_ptr.hpp>
#include <stdexcept>
#include "TypeTraits.h"

namespace Gadgetron{

    namespace Indexing {
        class Slice {};
        constexpr auto slice = Slice{};
    }
   template<class... ARGS>
   struct ValidIndex : std::integral_constant<bool, Core::all_of_v<Core::is_convertible_v<ARGS,size_t>...>> {};

   template<> struct ValidIndex<> : std::true_type {};

   template<class... ARGS>
   struct ValidIndex<Indexing::Slice,ARGS...> : ValidIndex<ARGS...> {};


   template<class T> class hoNDArray;


   template<class T, size_t D>
   class hoNDArrayView {
   public:
       hoNDArrayView& operator=(const hoNDArrayView&);
       hoNDArrayView& operator=(const hoNDArray<T>&);

        template<class... INDICES>
       std::enable_if_t<Core::all_of_v<Core::is_convertible_v<INDICES,size_t>...> && (sizeof...(INDICES) == D),T&>
       operator()(INDICES... indices);

       template<class... INDICES>
       std::enable_if_t<Core::all_of_v<Core::is_convertible_v<INDICES,size_t>...> && (sizeof...(INDICES) == D),const T&>
       operator()(INDICES... indices) const;
   private:
       friend class hoNDArray<T>;
       hoNDArrayView(const std::array<size_t,D>& strides, const std::array<size_t,D>& dimensions, T*);

       vector_td<size_t, D> strides;
       vector_td<size_t, D> dimensions;
       T* data;

   };



  template <typename T> class hoNDArray : public NDArray<T>
  {
  public:

    typedef NDArray<T> BaseClass;
    typedef float coord_type;
    typedef T value_type;
    using iterator = T*;
    using const_iterator = const T*;

    hoNDArray();

    explicit hoNDArray(const std::vector<size_t> &dimensions);
    explicit hoNDArray(const std::vector<size_t> *dimensions);
    explicit hoNDArray(boost::shared_ptr< std::vector<size_t> > dimensions);

    hoNDArray(const std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);
    hoNDArray(const std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct = false);
    hoNDArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct = false);

    hoNDArray(std::initializer_list<size_t> dimensions);
    hoNDArray(std::initializer_list<size_t> dimensions,T* data, bool delete_data_on_destruct = false);

    explicit hoNDArray(size_t len);
    hoNDArray(size_t sx, size_t sy);
    hoNDArray(size_t sx, size_t sy, size_t sz);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss);

    hoNDArray(size_t len, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct = false);
    hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct = false);

    virtual ~hoNDArray();

    // Copy constructors
    hoNDArray(const hoNDArray<T> &a);
    explicit hoNDArray(const hoNDArray<T> *a);
    //Move constructors

    //Move constructors
    hoNDArray(hoNDArray<T>&& a) noexcept;
    hoNDArray& operator=(hoNDArray&& rhs) noexcept;


    // Assignment operator
    hoNDArray& operator=(const hoNDArray& rhs);

    bool operator==(const hoNDArray& rhs) const;
    virtual void create(const std::vector<size_t>& dimensions);
    virtual void create(const std::vector<size_t> *dimensions);
    virtual void create(boost::shared_ptr< std::vector<size_t> > dimensions);

    virtual void create(std::initializer_list<size_t> dimensions);
    virtual void create(std::initializer_list<size_t> dimensions,T* data, bool delete_data_on_destruct = false);

    virtual void create(const std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct = false);
    virtual void create(const std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct = false);
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

    T& operator()( const std::vector<size_t>& ind );
    const T& operator()( const std::vector<size_t>& ind ) const;

    T& operator()( size_t x );
    const T& operator()( size_t x ) const;

    T& operator()( size_t x, size_t y );
    const T& operator()( size_t x, size_t y ) const;

    T& operator()( size_t x, size_t y, size_t z );
    const T& operator()( size_t x, size_t y, size_t z ) const;

    T& operator()( size_t x, size_t y, size_t z, size_t s );
    const T& operator()( size_t x, size_t y, size_t z, size_t s ) const;

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p );
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p ) const;

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r );
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r ) const;

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a );
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a ) const;

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q );
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q ) const;

    T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u );
    const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u ) const;

    template<class... INDICES>
    std::enable_if_t<ValidIndex<Indexing::Slice,INDICES...>::value, hoNDArray<T>>
    operator()(const Indexing::Slice&, const INDICES&... );



  template<class... INDICES, class = std::enable_if_t<Core::any_of_v<Core::is_same_v<INDICES,Indexing::Slice>...>> >
  auto operator()(const INDICES&... );

    void fill(T value);

    T* begin();
    const T* begin() const;

    T* end();
    const T* end() const;

    T& at( size_t idx );
    const T& at( size_t idx ) const;

    T& operator[]( size_t idx );
    const T& operator[]( size_t idx ) const;
    //T& operator()( size_t idx );
    //const T& operator()( size_t idx ) const;

    //T& operator()( const std::vector<size_t>& ind );
    //const T& operator()( const std::vector<size_t>& ind ) const;

    template<typename T2> 
    bool copyFrom(const hoNDArray<T2>& aArray)
    {
        try
        {
            if (!this->dimensions_equal(&aArray))
            {
                this->create(aArray.get_dimensions());
            }

            long long i;
#pragma omp parallel for default(none) private(i) shared(aArray)
            for (i = 0; i < (long long)elements_; i++)
            {
                data_[i] = static_cast<T>(aArray(i));
            }
        }
        catch (...)
        {
            GERROR_STREAM("Exceptions happened in hoNDArray::copyFrom(...) ... ");
            return false;
        }
        return true;
    }

    void get_sub_array(const std::vector<size_t>& start, std::vector<size_t>& size, hoNDArray<T>& out) const;

    virtual void print(std::ostream& os) const;
    virtual void printContent(std::ostream& os) const;

  [[deprecated("Use IO::write instead")]]
    virtual bool serialize(char*& buf, size_t& len) const;
    [[deprecated("Use IO::read instead")]]
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
