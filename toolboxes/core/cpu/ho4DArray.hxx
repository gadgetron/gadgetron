
namespace Gadgetron{

template <typename T> 
ho4DArray<T>::ho4DArray()
: BaseClass(), accesser_(NULL)
{
}

template <typename T> 
ho4DArray<T>::ho4DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss)
: accesser_(NULL)
{
    std::vector<unsigned long long> dim(4);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = ss;

    this->create(&dim);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho4DArray<T>::ho4DArray(std::vector<unsigned long long> *dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==4);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho4DArray<T>::ho4DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==4);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho4DArray<T>::ho4DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, T* data, bool delete_data_on_destruct)
: BaseClass(sx, sy, sz, ss, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho4DArray<T>::ho4DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==4);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho4DArray<T>::ho4DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==4);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho4DArray<T>::~ho4DArray()
{
    GADGET_CHECK_THROW(release_accesser());
}

template <typename T> 
ho4DArray<T>::ho4DArray(const ho4DArray<T>& a)
: BaseClass(a), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho4DArray<T>& ho4DArray<T>::operator=(const ho4DArray<T>& rhs)
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
void ho4DArray<T>::create(std::vector<unsigned long long>& dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho4DArray<T>::create(std::vector<unsigned long long> *dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho4DArray<T>::create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
{
    BaseClass::create(dimensions, data, delete_data_on_destruct);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
bool ho4DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss)
{
    try
    {
        std::vector<unsigned long long> dim(4);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = ss;

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
        GADGET_THROW("ho4DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss) ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho4DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, T* data, bool delete_data_on_destruct)
{
    try
    {
        std::vector<unsigned long long> dim(4);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = ss;

        this->create(&dim, data, delete_data_on_destruct);
        GADGET_CHECK_RETURN_FALSE(init_accesser());
    }
    catch(...)
    {
        GADGET_THROW("ho4DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, T* data, bool delete_data_on_destruct) ...");
        return false;
    }

    return true;
}

template <typename T> 
inline T& ho4DArray<T>::operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s)
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1] && z<(*dimensions_)[2] && s<(*dimensions_)[3]);
    return accesser_[s][z][y][x];
}

template <typename T> 
inline const T& ho4DArray<T>::operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s) const
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1] && z<(*dimensions_)[2] && s<(*dimensions_)[3]);
    return accesser_[s][z][y][x];
}

template <typename T> 
bool ho4DArray<T>::init_accesser()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(release_accesser());

        if ( elements_ > 0 )
        {
            unsigned long long sx = (*dimensions_)[0];
            unsigned long long sy = (*dimensions_)[1];
            unsigned long long sz = (*dimensions_)[2];
            unsigned long long ss = (*dimensions_)[3];

            unsigned long long y, z, s;

            accesser_ = new T***[ss];
            if( accesser_ == NULL) return false;

            accesser_[0] = new T**[sz*ss];
            if( accesser_[0] == NULL)
            {
                delete [] accesser_;
                return false;
            }
            for (s=1; s<ss; s++)
            {
                accesser_[s] = accesser_[s-1] + sz;
            }

            accesser_[0][0] = new T*[sy*sz*ss];
            if (accesser_[0][0] == NULL)
            {
                delete [] accesser_[0];
                delete [] accesser_;
                return false;
            }

            for (s=0; s<ss; s++)
            {
                for (z=0; z<sz; z++)
                {
                    accesser_[s][z] = accesser_[0][0] + s*sz*sy + z*sy;
                }
            }

            accesser_[0][0][0] = data_;
            for (s=0; s<ss; s++)
            {
                for (z=0; z<sz; z++)
                {
                    for (y=0; y<sy; y++)
                    {
                        accesser_[s][z][y] = accesser_[0][0][0] + s*sz*sy*sx + z*sy*sx + y*sx;
                    }
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
        GADGET_ERROR_MSG("Errors in ho4DArray<T>::init_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho4DArray<T>::release_accesser()
{
    try
    {
        if (accesser_ != NULL)
        {
            delete [] accesser_[0][0];
            delete [] accesser_[0];
            delete [] accesser_;
        }
        accesser_ = NULL;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in ho4DArray<T>::release_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
void ho4DArray<T>::print(std::ostream& os) const
{
    BaseClass::print(os);
    unsigned long long x, y, z, s;
    os << "-------------------------------------------" << std::endl;
    for (s=0; s<(*dimensions_)[3]; s++) 
    {
        for (z=0; z<(*dimensions_)[2]; z++) 
        {
            os << "ho4DArray (:, :, " << z << ", " << s << ") = " << std::endl;
            for (y=0; y<(*dimensions_)[1]; y++) 
            {
                os << "y " << y << "\t";
                for (x=0; x<(*dimensions_)[0]; x++)
                {
                    os << (*this)(x,y,z,s) << "\t";
                }
                os << std::endl;
            }
        }
    }
    os << "-------------------------------------------" << std::endl;
}

}
