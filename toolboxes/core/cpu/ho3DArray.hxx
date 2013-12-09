
namespace Gadgetron{

template <typename T> 
ho3DArray<T>::ho3DArray()
: BaseClass(), accesser_(NULL)
{
}

template <typename T> 
ho3DArray<T>::ho3DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz)
: accesser_(NULL)
{
    std::vector<unsigned long long> dim(3);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;

    this->create(&dim);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho3DArray<T>::ho3DArray(std::vector<unsigned long long> *dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==3);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho3DArray<T>::ho3DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==3);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho3DArray<T>::ho3DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct)
: BaseClass(sx, sy, sz, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho3DArray<T>::ho3DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==3);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho3DArray<T>::ho3DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==3);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho3DArray<T>::~ho3DArray()
{
    GADGET_CHECK_THROW(release_accesser());
}

template <typename T> 
ho3DArray<T>::ho3DArray(const ho3DArray<T>& a)
: BaseClass(a), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho3DArray<T>& ho3DArray<T>::operator=(const ho3DArray& rhs)
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
void ho3DArray<T>::create(std::vector<unsigned long long>& dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho3DArray<T>::create(std::vector<unsigned long long> *dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho3DArray<T>::create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
{
    BaseClass::create(dimensions, data, delete_data_on_destruct);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
bool ho3DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz)
{
    try
    {
        std::vector<unsigned long long> dim(3);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;

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
        GADGET_THROW("ho3DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz) ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho3DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct)
{
    try
    {
        std::vector<unsigned long long> dim(3);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;

        this->create(&dim, data, delete_data_on_destruct);
        GADGET_CHECK_RETURN_FALSE(init_accesser());
    }
    catch(...)
    {
        GADGET_THROW("ho3DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct) ...");
        return false;
    }

    return true;
}

template <typename T> 
inline T& ho3DArray<T>::operator()(unsigned long long x , unsigned long long y, unsigned long long z)
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1] && z<(*dimensions_)[2]);
    return accesser_[z][y][x];
}

template <typename T> 
inline const T& ho3DArray<T>::operator()(unsigned long long x , unsigned long long y, unsigned long long z) const
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1] && z<(*dimensions_)[[2]);
    return accesser_[z][y][x];
}

template <typename T> 
bool ho3DArray<T>::init_accesser()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(release_accesser());

        if ( elements_ > 0 )
        {
            unsigned long long sx = (*dimensions_)[0];
            unsigned long long sy = (*dimensions_)[1];
            unsigned long long sz = (*dimensions_)[2];

            unsigned long long y, z;

            accesser_ = new T**[sz];
            if( accesser_ == NULL) return false;

            accesser_[0] = new T*[sy*sz];
            if( accesser_[0] == NULL)
            {
                delete [] accesser_;
                return false;
            }
            for (z = 1; z < sz; z++)
            {
                accesser_[z] = accesser_[z-1] + sy;
            }

            accesser_[0][0] = data_;

            for (z=0; z<sz; z++)
            {
                for (y=0; y<sy; y++)
                {
                    accesser_[z][y] = accesser_[0][0] + (z*sy+y)*sx;
                }
            }
        }
        else
        {
            accesser_ = NULL;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in ho3DArray<T>::init_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho3DArray<T>::release_accesser()
{
    try
    {
        if (accesser_ != NULL)
        {
            delete [] accesser_[0];
            delete [] accesser_;
        }
        accesser_ = NULL;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in ho3DArray<T>::release_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
void ho3DArray<T>::print(std::ostream& os) const
{
    BaseClass::print(os);
    unsigned long long x, y, z;
    os << "-------------------------------------------" << std::endl;
    for (z=0; z<(*dimensions_)[2]; z++) 
    {
        os << "Array3D (:, :, " << z << ") = " << std::endl;
        for (y=0; y<(*dimensions_)[1]; y++) 
        {
            os << "y " << y << "\t";
            for (x=0; x<(*dimensions_)[0]; x++)
            {
                os << (*this)(x,y,z) << "\t";
            }
            os << std::endl;
        }
    }
    os << "-------------------------------------------" << std::endl;
}

}
