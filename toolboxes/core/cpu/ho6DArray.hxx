
namespace Gadgetron{

template <typename T> 
ho6DArray<T>::ho6DArray()
: BaseClass(), accesser_(NULL)
{
}

template <typename T> 
ho6DArray<T>::ho6DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr)
: accesser_(NULL)
{
    std::vector<unsigned long long> dim(6);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = ss;
    dim[4] = sp;
    dim[5] = sr;

    this->create(&dim);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho6DArray<T>::ho6DArray(std::vector<unsigned long long> *dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==6);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho6DArray<T>::ho6DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==6);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho6DArray<T>::ho6DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, T* data, bool delete_data_on_destruct)
: BaseClass(sx, sy, sz, ss, sp, sr, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions_->size()==6);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho6DArray<T>::ho6DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==6);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho6DArray<T>::ho6DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==6);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho6DArray<T>::~ho6DArray()
{
    GADGET_CHECK_THROW(release_accesser());
}

template <typename T> 
ho6DArray<T>::ho6DArray(const ho6DArray<T>& a)
: BaseClass(a), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho6DArray<T>& ho6DArray<T>::operator=(const ho6DArray<T>& rhs)
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
void ho6DArray<T>::create(std::vector<unsigned long long>& dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho6DArray<T>::create(std::vector<unsigned long long> *dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho6DArray<T>::create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct)
{
    BaseClass::create(dimensions, data, delete_data_on_destruct);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
bool ho6DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr)
{
    try
    {
        std::vector<unsigned long long> dim(6);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = ss;
        dim[4] = sp;
        dim[5] = sr;

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
        GADGET_THROW("ho6DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr) ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho6DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, T* data, bool delete_data_on_destruct)
{
    try
    {
        std::vector<unsigned long long> dim(6);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = ss;
        dim[4] = sp;
        dim[5] = sr;

        this->create(&dim);
        GADGET_CHECK_RETURN_FALSE(init_accesser());
    }
    catch(...)
    {
        GADGET_THROW("ho6DArray<T>::createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, T* data, bool delete_data_on_destruct) ...");
        return false;
    }

    return true;
}

template <typename T> 
inline T& ho6DArray<T>::operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p, unsigned long long r)
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1] && z<(*dimensions_)[2] && s<(*dimensions_)[3] && p<(*dimensions_)[4] && r<(*dimensions_)[5]);
    return accesser_[r][p][s][z][y][x];
}

template <typename T> 
inline const T& ho6DArray<T>::operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p, unsigned long long r) const
{
    GADGET_DEBUG_CHECK_THROW(x<(*dimensions_)[0] && y<(*dimensions_)[1] && z<(*dimensions_)[2] && s<(*dimensions_)[3] && p<(*dimensions_)[4] && r<(*dimensions_)[5]);
    return accesser_[r][p][s][z][y][x];
}

template <typename T> 
bool ho6DArray<T>::init_accesser()
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
            unsigned long long sp = (*dimensions_)[4];
            unsigned long long sr = (*dimensions_)[5];

            unsigned long long y, z, s, p, r;

            accesser_ = new T*****[sr];
            if( accesser_ == NULL) return false;

            accesser_[0] = new T****[sp*sr];
            if( accesser_[0] == NULL)
            {
                delete [] accesser_;
                return false;
            }
            for (r=1; r<sr; r++)
            {
                accesser_[r] = accesser_[r-1] + sp;
            }

            accesser_[0][0] = new T***[ss*sp*sr];
            if (accesser_[0][0] == NULL)
            {
                delete [] accesser_[0];
                delete [] accesser_;
                return false;
            }

            for (r=0; r<sr; r++)
            {
                for (p=0; p<sp; p++)
                {
                    accesser_[r][p] = accesser_[0][0] + r*sp*ss + p*ss;
                }
            }

            accesser_[0][0][0] = new T**[sz*ss*sp*sr];
            if (accesser_[0][0][0] == NULL)
            {
                delete [] accesser_[0][0];
                delete [] accesser_[0];
                delete [] accesser_;
                return false;
            }

            for (r=0; r<sr; r++)
            {
                for (p=0; p<sp; p++)
                {
                    for (s=0; s<ss; s++)
                    {
                        accesser_[r][p][s] = accesser_[0][0][0] 
                                                + r*sp*ss*sz 
                                                + p*ss*sz 
                                                + s*sz;
                    }
                }
            }

            accesser_[0][0][0][0] = new T*[sy*sz*ss*sp*sr];
            if (accesser_[0][0][0][0] == NULL)
            {
                delete [] accesser_[0][0][0];
                delete [] accesser_[0][0];
                delete [] accesser_[0];
                delete [] accesser_;
                return false;
            }

            for (r=0; r<sr; r++)
            {
                for (p=0; p<sp; p++)
                {
                    for (s=0; s<ss; s++)
                    {
                        for (z=0; z<sz; z++)
                        {
                            accesser_[r][p][s][z] = accesser_[0][0][0][0] 
                                                        + r*sp*ss*sz*sy 
                                                        + p*ss*sz*sy 
                                                        + s*sz*sy 
                                                        + z*sy;
                        }
                    }
                }
            }

            accesser_[0][0][0][0][0] = data_;
            for (r=0; r<sr; r++)
            {
                for (p=0; p<sp; p++)
                {
                    for (s=0; s<ss; s++)
                    {
                        for (z=0; z<sz; z++)
                        {
                            for (y=0; y<sy; y++)
                            {
                                accesser_[r][p][s][z][y] = accesser_[0][0][0][0][0] 
                                                                + r*sp*ss*sz*sy*sx 
                                                                + p*ss*sz*sy*sx 
                                                                + s*sz*sy*sx 
                                                                + z*sy*sx 
                                                                + y*sx;
                            }
                        }
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
        GADGET_ERROR_MSG("Errors in ho6DArray<T>::init_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho6DArray<T>::release_accesser()
{
    try
    {
        if (accesser_ != NULL)
        {
            delete [] accesser_[0][0][0][0];
            delete [] accesser_[0][0][0];
            delete [] accesser_[0][0];
            delete [] accesser_[0];
            delete [] accesser_;
        }
        accesser_ = NULL;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in ho6DArray<T>::release_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
void ho6DArray<T>::print(std::ostream& os) const
{
    BaseClass::print(os);
    unsigned long long x, y, z, s, p, r;
    os << "-------------------------------------------" << std::endl;
    for (r=0; r<(*dimensions_)[5]; r++) 
    {
        for (p=0; p<(*dimensions_)[4]; p++) 
        {
            for (s=0; s<(*dimensions_)[3]; s++) 
            {
                for (z=0; z<(*dimensions_)[2]; z++) 
                {
                    os << "ho6DArray (:, :, " << z << ", " << s << ", " << p << ", " << r << ") = " << std::endl;
                    for (y=0; y<(*dimensions_)[1]; y++) 
                    {
                        os << "y " << y << "\t";
                        for (x=0; x<(*dimensions_)[0]; x++)
                        {
                            os << (*this)(x,y,z,s,p,r) << "\t";
                        }
                        os << std::endl;
                    }
                }
            }
        }
    }
    os << "-------------------------------------------" << std::endl;
}

}
