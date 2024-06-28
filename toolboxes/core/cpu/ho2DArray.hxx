
namespace Gadgetron{

template <typename T> 
ho2DArray<T>::ho2DArray()
: BaseClass(), accesser_(NULL)
{
}

template <typename T> 
ho2DArray<T>::ho2DArray(size_t sx, size_t sy)
: accesser_(NULL)
{
    this->create(sx,sy);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(const std::vector<size_t>& dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions.size()==2);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions.size()==2);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct)
: BaseClass(sx, sy, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::~ho2DArray()
{
    release_accesser();
}

template <typename T> 
ho2DArray<T>::ho2DArray(const ho2DArray<T>& a)
: BaseClass(a), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>& ho2DArray<T>::operator=(const ho2DArray<T>& rhs)
{
    if ( &rhs == this ) return *this;

    if ( rhs.get_number_of_elements() == 0 )
    {
        this->clear();
        GADGET_CHECK_THROW(init_accesser());
        return *this;
    }

    std::vector<size_t> dim;
    rhs.get_dimensions(dim);
    if (this->dimensions_equal(dim)) 
    {
        memcpy(this->data_, rhs.data_, this->elements_*sizeof(T));
    }
    else
    {
        this->deallocate_memory();
        this->data_ = 0;
        this->dimensions_ = rhs.dimensions_;
        this->allocate_memory();
        memcpy( this->data_, rhs.data_, this->elements_*sizeof(T) );

        GADGET_CHECK_THROW(init_accesser());
    }

    return *this;
}

template <typename T> 
void ho2DArray<T>::create(const std::vector<size_t>& dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho2DArray<T>::create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct)
{
    BaseClass::create(dimensions, data, delete_data_on_destruct);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
bool ho2DArray<T>::createArray(size_t sx, size_t sy)
{
    try
    {
        std::vector<size_t> dim(2);
        dim[0] = sx;
        dim[1] = sy;

        if ( !this->dimensions_equal(dim) )
        {
            this->create(dim);
            GADGET_CHECK_RETURN_FALSE(init_accesser());
        }
        else
        {
            memset(this->data_, 0, sizeof(T)*this->elements_);
        }
    }
    catch(...)
    {
        GERROR_STREAM("ho2DArray<T>::createArray(size_t sx, size_t sy) ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho2DArray<T>::createArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct)
{
    try
    {
        std::vector<size_t> dim(2);
        dim[0] = sx;
        dim[1] = sy;

        this->create(dim, data, delete_data_on_destruct);
        GADGET_CHECK_RETURN_FALSE(init_accesser());
    }
    catch(...)
    {
        GERROR_STREAM("ho2DArray<T>::createArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct) ...");
        return false;
    }

    return true;
}

template <typename T> 
inline T& ho2DArray<T>::operator()(size_t x , size_t y)
{
    GADGET_DEBUG_CHECK_THROW(x<dimensions_[0] && y<dimensions_[1]);
    return accesser_[y][x];
}

template <typename T> 
inline const T& ho2DArray<T>::operator()(size_t x , size_t y) const
{
    GADGET_DEBUG_CHECK_THROW(x<dimensions_[0] && y<dimensions_[1]);
    return accesser_[y][x];
}

template <typename T> 
bool ho2DArray<T>::init_accesser()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(release_accesser());

        if ( elements_ > 0 )
        {
            size_t sx = dimensions_[0];
            size_t sy = dimensions_[1];

            accesser_ = new T*[sy];
            if( accesser_ == NULL) return false;

            accesser_[0] = data_;
            for (size_t y=1; y<sy; y++)
            {
                accesser_[y] = accesser_[y-1] + sx;
            }
        }
        else
        {
            accesser_ = NULL;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in ho2DArray<T>::init_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho2DArray<T>::release_accesser()
{
    try
    {
        if (accesser_ != NULL)
        {
            delete [] accesser_;
        }
        accesser_ = NULL;
    }
    catch(...)
    {
        GERROR_STREAM("Errors in ho2DArray<T>::release_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
void ho2DArray<T>::print(std::ostream& os) const
{
    BaseClass::print(os);
    size_t x, y;
    os << "-------------------------------------------" << std::endl;
    for (y=0; y<dimensions_[1]; y++) 
    {
        os << "y " << y << "\t";
        for (x=0; x<dimensions_[0]; x++)
        {
            os << (*this)(x,y) << "\t";
        }
        os << std::endl;
    }
    os << "-------------------------------------------" << std::endl;
}

}
