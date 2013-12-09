
namespace Gadgetron
{

template <typename T> 
hoNDArray<T>::hoNDArray() : NDArray<T>::NDArray() 
{
}

template <typename T> 
hoNDArray<T>::hoNDArray(std::vector<unsigned long long> *dimensions) : NDArray<T>::NDArray()
{
    this->create(dimensions);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long len) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(1);
    dim[0] = len;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(2);
    dim[0] = sx;
    dim[1] = sy;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(3);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(4);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(5);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(6);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(7);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(8);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    dim[7] = ss;
    this->create(dim);
}

template <typename T> 
hoNDArray<T>::hoNDArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    this->create(dimensions,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long len, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(1);
    dim[0] = len;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(2);
    dim[0] = sx;
    dim[1] = sy;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(3);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(4);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(5);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(6);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(7);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    std::vector<unsigned long long> dim(8);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    dim[7] = ss;
    this->create(&dim,data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::hoNDArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions) : NDArray<T>::NDArray()
{
    this->create(dimensions.get());
}

template <typename T> 
hoNDArray<T>::hoNDArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, 
        T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
{
    this->create(dimensions.get(),data,delete_data_on_destruct);
}

template <typename T> 
hoNDArray<T>::~hoNDArray()
{
    if (this->delete_data_on_destruct_)
    {
        deallocate_memory();
    }
}

template <typename T> 
hoNDArray<T>::hoNDArray(const hoNDArray<T>& a)
{
    this->data_ = 0;
    this->dimensions_ = a.dimensions_;
    this->offsetFactors_ = a.offsetFactors_;
    allocate_memory();
    memcpy( this->data_, a.data_, this->elements_*sizeof(T) );
}

template <typename T> 
hoNDArray<T>& hoNDArray<T>::operator=(const hoNDArray<T>& rhs)
{
    if ( &rhs == this ) return *this;

    if ( rhs.get_number_of_elements() == 0 )
    {
        this->clear();
        return *this;
    }

    // Are the dimensions the same? Then we can just memcpy
    if (this->dimensions_equal(&rhs))
    {
        memcpy(this->data_, rhs.data_, this->elements_*sizeof(T));
    }
    else
    {
        deallocate_memory();
        this->data_ = 0;
        this->dimensions_ = rhs.dimensions_;
        this->offsetFactors_ = rhs.offsetFactors_;
        allocate_memory();
        memcpy( this->data_, rhs.data_, this->elements_*sizeof(T) );
    }
    return *this;
}

template <typename T> 
void hoNDArray<T>::create(std::vector<unsigned long long>& dimensions)
{
    if ( this->dimensions_equal(&dimensions) )
    {
        memset(this->data_, 0, sizeof(T)*this->elements_);
        return;
    }

    this->clear();
    BaseClass::create(dimensions);
}

template <typename T> 
void hoNDArray<T>::create(std::vector<unsigned int> *dimensions) 
{
    std::vector<unsigned long long> dims(dimensions->size(), 0);
    for ( unsigned int ii=0; ii<dimensions->size(); ii++ )
    {
        dims[ii] = (*dimensions)[ii];
    }

    this->create(&dims);
}

template <typename T> 
void hoNDArray<T>::create(std::vector<unsigned long long> *dimensions)
{
    if ( this->dimensions_equal(dimensions) )
    {
        memset(this->data_, 0, sizeof(T)*this->elements_);
        return;
    }

    this->clear();
    BaseClass::create(dimensions);
}

template <typename T> 
void hoNDArray<T>::create(std::vector<unsigned int> *dimensions, T* data, bool delete_data_on_destruct) 
{
    std::vector<unsigned long long> dims(dimensions->size(), 0);
    for ( unsigned int ii=0; ii<dimensions->size(); ii++ )
    {
        dims[ii] = (*dimensions)[ii];
    }

    this->create(&dims, data, delete_data_on_destruct);
}

template <typename T> 
void hoNDArray<T>::create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct) 
{
    if (!data)
    {
        GADGET_THROW("NDArray<T>::create: 0x0 pointer provided");
    }

    this->clear();
    BaseClass::create(dimensions, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(boost::shared_ptr< std::vector<unsigned long long> > dimensions)
{
    this->create(dimensions.get());
}

template <typename T> 
inline void hoNDArray<T>::create(boost::shared_ptr<std::vector<unsigned long long>  > dimensions, T* data, bool delete_data_on_destruct)
{
    this->create(dimensions.get(), data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long len)
{
    std::vector<unsigned long long> dim(1);
    dim[0] = len;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy)
{
    std::vector<unsigned long long> dim(2);
    dim[0] = sx;
    dim[1] = sy;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz)
{
    std::vector<unsigned long long> dim(3);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st)
{
    std::vector<unsigned long long> dim(4);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp)
{
    std::vector<unsigned long long> dim(5);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq)
{
    std::vector<unsigned long long> dim(6);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr)
{
    std::vector<unsigned long long> dim(7);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss)
{
    std::vector<unsigned long long> dim(8);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    dim[7] = ss;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su)
{
    std::vector<unsigned long long> dim(9);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    dim[7] = ss;
    dim[8] = su;
    this->create(dim);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long len, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(1);
    dim[0] = len;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(2);
    dim[0] = sx;
    dim[1] = sy;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(3);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(4);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(5);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(6);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(7);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(8);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    dim[7] = ss;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
inline void hoNDArray<T>::create(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long st, unsigned long long sp, unsigned long long sq, unsigned long long sr, unsigned long long ss, unsigned long long su, T* data, bool delete_data_on_destruct)
{
    std::vector<unsigned long long> dim(9);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = st;
    dim[4] = sp;
    dim[5] = sq;
    dim[6] = sr;
    dim[7] = ss;
    dim[8] = su;
    this->create(&dim, data, delete_data_on_destruct);
}

template <typename T> 
void hoNDArray<T>::fill(T value)
{
    std::fill(this->get_data_ptr(), this->get_data_ptr()+this->get_number_of_elements(), value);
}

template <typename T> 
inline T* hoNDArray<T>::begin()
{
    return this->data_;
}

template <typename T> 
inline const T* hoNDArray<T>::begin() const
{
    return this->data_;
}

template <typename T> 
inline T* hoNDArray<T>::end()
{
    return (this->data_+this->elements_);
}

template <typename T> 
inline const T* hoNDArray<T>::end() const
{
    return (this->data_+this->elements_);
}

template <typename T> 
inline T& hoNDArray<T>::at( unsigned long long idx )
{
    /*if( idx >= this->get_number_of_elements() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::at(): index out of range."));
    }*/
    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    return this->get_data_ptr()[idx];
}

template <typename T> 
inline const T& hoNDArray<T>::at( unsigned long long idx ) const
{
    /*if( idx >= this->get_number_of_elements() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::at(): index out of range."));
    }*/
    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    return this->get_data_ptr()[idx];
}

template <typename T> 
inline T hoNDArray<T>::operator[]( unsigned long long idx )
{
    /*if( idx >= this->get_number_of_elements() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator[]: index out of range."));
    }*/
    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    return this->get_data_ptr()[idx];
}

template <typename T> 
inline T& hoNDArray<T>::operator()( unsigned long long idx )
{
    /*if( idx >= this->get_number_of_elements() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    }*/
    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    return this->get_data_ptr()[idx];
}

template <typename T> 
inline const T& hoNDArray<T>::operator()( unsigned long long idx ) const
{
    /*if( idx >= this->get_number_of_elements() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    }*/
    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    return this->get_data_ptr()[idx];
}

template <typename T> 
inline T& hoNDArray<T>::operator()( const std::vector<unsigned long long>& ind )
{
    long long idx = this->calculate_offset(ind);
    /*if( idx >= this->get_number_of_elements() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    }*/
    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    return this->get_data_ptr()[idx];
}

template <typename T> 
inline const T& hoNDArray<T>::operator()( const std::vector<unsigned long long>& ind ) const
{
    long long idx = this->calculate_offset(ind);
    /*if( idx >= this->get_number_of_elements() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    }*/
    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    return this->get_data_ptr()[idx];
}

template <typename T> 
void hoNDArray<T>::get_sub_array(const std::vector<unsigned long long>& start, std::vector<unsigned long long>& size, hoNDArray<T>& out)
{
    if ( start.size() != size.size() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
    }

    if ( start.size() != (*dimensions_).size() )
    {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
    }

    out.create(&size);

    if ( out.get_number_of_elements() == this->get_number_of_elements() )
    {
        out = *this;
        return;
    }

    std::vector<unsigned long long> end(start.size());

    unsigned long long ii;
    for ( ii=0; ii<start.size(); ii++ )
    {
        end[ii] = start[ii] + size[ii] - 1;
        if ( end[ii] >= (*dimensions_)[ii] )
        {
            BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
        }
    }


}

template <typename T> 
void hoNDArray<T>::print(std::ostream& os) const
{
    using namespace std;

    os.unsetf(std::ios::scientific);
    os.setf(ios::fixed);

    unsigned long long i;

    os << "--------------Gagdgetron ND Array -------------" << endl;
    os << "Array dimension is : " << dimensions_->size() << endl;

    os << "Array size is : ";
    for (i=0; i<dimensions_->size(); i++ ) 
        os << (*dimensions_)[i] << " "; 
    os << endl;

    int elemTypeSize = sizeof(T);
    std::string elemTypeName = std::string(typeid(T).name());

    os << "Array data type is : " << elemTypeName << std::endl;
    os << "Byte number for each element is : " << elemTypeSize << std::endl;
    os << "Number of array size in bytes is : ";
    os << elements_*elemTypeSize << std::endl;

    //os << "-------------------------------------------" << std::endl;
    //unsigned long long numOfPrints = 20;
    //if ( this->elements_ < numOfPrints ) numOfPrints = this->elements_;
    //for (i=0; i<numOfPrints; i++) 
    //{
    //    os << i << " = " << (*this)(i) << std::endl;
    //}
    //os << "-------------------------------------------" << std::endl;

    os << std::endl;
}

template <typename T> 
void hoNDArray<T>::printContent(std::ostream& os) const
{
    using namespace std;

    os.unsetf(std::ios::scientific);
    os.setf(ios::fixed);

    print(os);

    //unsigned long long i;

    //os << "-------------------------------------------" << std::endl;
    //unsigned long long numOfPrints = this->elements_;
    //if ( this->elements_ < numOfPrints ) numOfPrints = this->elements_;
    //for (i=0; i<numOfPrints; i++) 
    //{
    //    os << i << " = " << (*this)(i) << std::endl;
    //}
    //os << "-------------------------------------------" << std::endl;
    //os << std::endl;
}

template <typename T> 
void hoNDArray<T>::allocate_memory()
{
    deallocate_memory();

    this->elements_ = (*this->dimensions_)[0];
    for (unsigned long long i = 1; i < this->dimensions_->size(); i++)
    {
        this->elements_ *= (*this->dimensions_)[i];
    }

    if ( this->elements_ > 0 )
    {
        this->_allocate_memory(this->elements_, &this->data_);

        if( this->data_ == 0x0 )
        {
            BOOST_THROW_EXCEPTION( bad_alloc("hoNDArray<>::allocate memory failed"));
        }

        this->delete_data_on_destruct_ = true;

        // memset(this->data_, 0, sizeof(T)*this->elements_);
    }
}

template <typename T> 
void hoNDArray<T>::deallocate_memory()
{
    if( this->data_ )
    {
        this->_deallocate_memory( this->data_ );
        this->data_ = 0x0;
    }
}

template <typename T> 
inline void hoNDArray<T>::_allocate_memory( unsigned long long size, float** data )
{
    #ifdef USE_MKL
        *data = (float*) mkl_malloc(size*sizeof(float), 4);
    #else
        *data = (float*) malloc( size*sizeof(float) );
    #endif // USE_MKL

    // *data = new float[size];
}

template <typename T> 
inline void hoNDArray<T>::_deallocate_memory( float* data )
{
    #ifdef USE_MKL
        mkl_free(data);
    #else
        free(data);
    #endif // USE_MKL

    // delete [] data;
}

template <typename T> 
inline void hoNDArray<T>::_allocate_memory( unsigned long long size, double** data )
{
    #ifdef USE_MKL
        *data = (double*) mkl_malloc(size*sizeof(double), 4);
    #else
        *data = (double*) malloc( size*sizeof(double) );
    #endif // USE_MKL

    //*data = new double[size];
}

template <typename T> 
inline void hoNDArray<T>::_deallocate_memory( double* data )
{
    #ifdef USE_MKL
        mkl_free(data);
    #else
        free(data);
    #endif // USE_MKL

    //delete [] data;
}

template <typename T> 
inline void hoNDArray<T>::_allocate_memory( unsigned long long size, std::complex<float>** data )
{
    //GADGET_MSG("into hoNDArray<T>::_allocate_memory(std::complex<float>) : " << size);

    #ifdef USE_MKL
        *data = (std::complex<float>*) mkl_malloc(size*sizeof(std::complex<float>), 4);
    #else
        *data = (std::complex<float>*) malloc( size*sizeof(std::complex<float>) );
    #endif // USE_MKL

    //*data = new std::complex<float>[size];
}

template <typename T> 
inline void hoNDArray<T>::_deallocate_memory( std::complex<float>* data )
{
    #ifdef USE_MKL
        mkl_free(data);
    #else
        free(data);
    #endif // USE_MKL

    //delete [] data;
}

template <typename T> 
inline void hoNDArray<T>::_allocate_memory( unsigned long long size, std::complex<double>** data )
{
    #ifdef USE_MKL
        *data = (std::complex<double>*) mkl_malloc(size*sizeof(std::complex<double>), 4);
    #else
        *data = (std::complex<double>*) malloc( size*sizeof(std::complex<double>) );
    #endif // USE_MKL

    //*data = new std::complex<double>[size];
}

template <typename T> 
inline void hoNDArray<T>::_deallocate_memory( std::complex<double>* data )
{
    #ifdef USE_MKL
        mkl_free(data);
    #else
        free(data);
    #endif // USE_MKL

    //delete [] data;
}

template <typename T> 
inline void hoNDArray<T>::_allocate_memory( unsigned long long size, float_complext** data )
{
    *data = (float_complext*) malloc( size*sizeof(float_complext) );
}

template <typename T> 
inline void hoNDArray<T>::_deallocate_memory( float_complext* data )
{
    free( data );
}

template <typename T> 
inline void hoNDArray<T>::_allocate_memory( unsigned long long size, double_complext** data )
{
    *data = (double_complext*) malloc( size*sizeof(double_complext) );
}

template <typename T> 
inline void hoNDArray<T>::_deallocate_memory( double_complext* data )
{
    free( data );
}

template <typename T> 
bool hoNDArray<T>::serialize(char*& buf, unsigned long long& len) const
{
    try
    {
        if ( buf != NULL ) delete[] buf;

        unsigned long long NDim = dimensions_->size();

        // number of dimensions + dimension vector + contents
        len = sizeof(unsigned long long) + sizeof(unsigned long long)*NDim + sizeof(T)*elements_;

        buf = new char[len];
        GADGET_CHECK_RETURN_FALSE( buf != NULL );

        memcpy(buf, &NDim, sizeof(unsigned long long));
        if ( NDim > 0 )
        {
            memcpy(buf+sizeof(unsigned long long), &((*dimensions_)[0]), sizeof(unsigned long long)*NDim);
            memcpy(buf+sizeof(unsigned long long)+sizeof(unsigned long long)*NDim, this->data_, sizeof(T)*elements_);
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors happened in hoNDArray<T>::serialize(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool hoNDArray<T>::deserialize(char* buf, unsigned long long& len)
{
    try
    {
        unsigned long long NDim;
        memcpy(&NDim, buf, sizeof(unsigned long long));

        if ( NDim > 0 )
        {
            std::vector<unsigned long long> dimensions(NDim);
            memcpy(&dimensions[0], buf+sizeof(unsigned long long), sizeof(unsigned long long)*NDim);

            // allocate memory
            this->create(&dimensions);

            // copy the content
            memcpy(this->data_, buf+sizeof(unsigned long long)+sizeof(unsigned long long)*NDim, sizeof(T)*elements_);
        }
        else
        {
            this->clear();
        }

        len = sizeof(unsigned long long)+sizeof(unsigned long long)*NDim+sizeof(T)*elements_;
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors happended in hoNDArray<T>::deserialize(char* buf, unsigned long long len) ...");
        return false;
    }

    return true;
}

//template class hoNDArray<short>;
//template class hoNDArray<unsigned short>;
//template class hoNDArray<int>;
//template class hoNDArray<unsigned long long>;
//template class hoNDArray<float>;
//template class hoNDArray<double>;
//template class hoNDArray< std::complex<float> >;
//template class hoNDArray< std::complex<double> >;

//template class hoNDArray< Gadgetron::complext<float> >;
//template class hoNDArray< Gadgetron::complext<double> >;

}

