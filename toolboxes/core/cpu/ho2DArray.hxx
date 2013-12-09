
namespace Gadgetron{

template <typename T> 
ho2DArray<T>::ho2DArray()
: BaseClass(), accesser_(NULL)
{
}

template <typename T> 
ho2DArray<T>::ho2DArray(unsigned long long sx, unsigned long long sy)
: accesser_(NULL)
{
    std::vector<unsigned long long> dim(2);
    dim[0] = sx;
    dim[1] = sy;

    this->create(&dim);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(std::vector<unsigned long long> *dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==2);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==2);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct)
: BaseClass(sx, sy, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==2);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::ho2DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==2);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho2DArray<T>::~ho2DArray()
{
    GADGET_CHECK_THROW(release_accesser());
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

    if (this->dimensions_equal(&rhs)) 
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
void ho2DArray<T>::create(std::vector<unsigned long long>& dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho2DArray<T>::create(std::vector<unsigned long long> *dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho2DArray<T>::create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
{
    BaseClass::create(dimensions, data, delete_data_on_destruct);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
bool ho2DArray<T>::createArray(unsigned long long sx, unsigned long long sy)
{
    try
    {
        std::vector<unsigned long long> dim(2);
        dim[0] = sx;
        dim[1] = sy;

        if ( !this->dimensions_equal(&dim) )
        {
            this->create(&dim);
            GADGET_CHECK_RETURN_FALSE(init_accesser());
        }
        else
        {
            memset(this->data_, 0, sizeof(T)*this->elements_);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("ho2DArray<T>::createArray(unsigned long long sx, unsigned long long sy) ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho2DArray<T>::createArray(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct)
{
    try
    {
        std::vector<unsigned long long> dim(2);
        dim[0] = sx;
        dim[1] = sy;

        this->create(&dim, data, delete_data_on_destruct);
        GADGET_CHECK_RETURN_FALSE(init_accesser());
    }
    catch(...)
    {
        GADGET_ERROR_MSG("ho2DArray<T>::createArray(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct) ...");
        return false;
    }

    return true;
}

template <typename T> 
inline T& ho2DArray<T>::operator()(unsigned long long x , unsigned long long y)
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1]);
    return accesser_[y][x];
}

template <typename T> 
inline const T& ho2DArray<T>::operator()(unsigned long long x , unsigned long long y) const
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1]);
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
            unsigned long long sx = (*dimensions_)[0];
            unsigned long long sy = (*dimensions_)[1];

            accesser_ = new T*[sy];
            if( accesser_ == NULL) return false;

            accesser_[0] = data_;
            for (unsigned long long y=1; y<sy; y++)
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
        GADGET_ERROR_MSG("Errors in ho2DArray<T>::init_accesser() ...");
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
        GADGET_ERROR_MSG("Errors in ho2DArray<T>::release_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
void ho2DArray<T>::print(std::ostream& os) const
{
    BaseClass::print(os);
    unsigned long long x, y;
    os << "-------------------------------------------" << std::endl;
    for (y=0; y<(*dimensions_)[1]; y++) 
    {
        os << "y " << y << "\t";
        for (x=0; x<(*dimensions_)[0]; x++)
        {
            os << (*this)(x,y) << "\t";
        }
        os << std::endl;
    }
    os << "-------------------------------------------" << std::endl;
}

}
