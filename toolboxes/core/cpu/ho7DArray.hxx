
namespace Gadgetron{

template <typename T> 
ho7DArray<T>::ho7DArray()
: BaseClass(), accesser_(NULL)
{
}

template <typename T> 
ho7DArray<T>::ho7DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa)
: accesser_(NULL)
{
    std::vector<size_t> dim(7);
    dim[0] = sx;
    dim[1] = sy;
    dim[2] = sz;
    dim[3] = ss;
    dim[4] = sp;
    dim[5] = sr;
    dim[6] = sa;

    this->create(dim);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho7DArray<T>::ho7DArray(const std::vector<size_t>& dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions.size()==7);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho7DArray<T>::ho7DArray(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions.size()==7);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho7DArray<T>::ho7DArray(boost::shared_ptr< std::vector<size_t> > dimensions)
: BaseClass(dimensions), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==7);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho7DArray<T>::ho7DArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct)
: BaseClass(dimensions, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions->size()==7);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho7DArray<T>::ho7DArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa, T* data, bool delete_data_on_destruct)
: BaseClass(sx, sy, sz, ss, sp, sr, sa, data, delete_data_on_destruct), accesser_(NULL)
{
    GADGET_CHECK_THROW(dimensions_.size()==7);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho7DArray<T>::~ho7DArray()
{
    release_accesser();
}

template <typename T> 
ho7DArray<T>::ho7DArray(const ho7DArray<T>& a)
: BaseClass(a), accesser_(NULL)
{
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
ho7DArray<T>& ho7DArray<T>::operator=(const ho7DArray<T>& rhs)
{
    if ( &rhs == this ) return *this;

    if ( rhs.get_number_of_elements() == 0 )
    {
        this->clear();
        GADGET_CHECK_THROW(init_accesser());
        return *this;
    }

    if (this->dimensions_equal(rhs)) 
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
void ho7DArray<T>::create(const std::vector<size_t>& dimensions)
{
    BaseClass::create(dimensions);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
void ho7DArray<T>::create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct)
{
    BaseClass::create(dimensions, data, delete_data_on_destruct);
    GADGET_CHECK_THROW(init_accesser());
}

template <typename T> 
bool ho7DArray<T>::createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa)
{
    try
    {
        std::vector<size_t> dim(7);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = ss;
        dim[4] = sp;
        dim[5] = sr;
        dim[6] = sa;

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
        GADGET_THROW("ho7DArray<T>::createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa) ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho7DArray<T>::createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa, T* data, bool delete_data_on_destruct)
{
    try
    {
        std::vector<size_t> dim(7);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = ss;
        dim[4] = sp;
        dim[5] = sr;
        dim[6] = sa;

        this->create(dim, data, delete_data_on_destruct);
        GADGET_CHECK_RETURN_FALSE(init_accesser());
    }
    catch(...)
    {
        GADGET_THROW("ho7DArray<T>::createArray(size_t sx, size_t sy, size_t sz, size_t ss, size_t sp, size_t sr, size_t sa, T* data, bool delete_data_on_destruct) ...");
        return false;
    }

    return true;
}

template <typename T> 
inline T& ho7DArray<T>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a)
{
    GADGET_DEBUG_CHECK_THROW(x<dimensions_[0] && y<dimensions_[1] && z<dimensions_[2] && s<dimensions_[3] && p<dimensions_[4] && r<dimensions_[5] && a<dimensions_[6]);
    return accesser_[a][r][p][s][z][y][x];
}

template <typename T> 
inline const T& ho7DArray<T>::operator()(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const
{
    GADGET_DEBUG_CHECK_THROW(x<dimensions_[0] && y<dimensions_[1] && z<dimensions_[2] && s<dimensions_[3] && p<dimensions_[4] && r<dimensions_[5] && a<dimensions_[6]);
    return accesser_[a][r][p][s][z][y][x];
}

template <typename T> 
bool ho7DArray<T>::init_accesser()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(release_accesser());

        if ( elements_ > 0 )
        {
            size_t sx = dimensions_[0];
            size_t sy = dimensions_[1];
            size_t sz = dimensions_[2];
            size_t ss = dimensions_[3];
            size_t sp = dimensions_[4];
            size_t sr = dimensions_[5];
            size_t sa = dimensions_[6];

            size_t y, z, s, p, r, a;

            accesser_ = new T******[sa];
            if( accesser_ == NULL) return false;

            accesser_[0] = new T*****[sr*sa];
            if( accesser_[0] == NULL)
            {
                delete [] accesser_;
                return false;
            }
            for (a=1; a<sa; a++)
            {
                accesser_[a] = accesser_[a-1] + sr;
            }

            accesser_[0][0] = new T****[sp*sr*sa];
            if (accesser_[0][0] == NULL)
            {
                delete [] accesser_[0];
                delete [] accesser_;
                return false;
            }

            for (a=0; a<sa; a++)
            {
                for (r=0; r<sr; r++)
                {
                    accesser_[a][r] = accesser_[0][0] + a*sr*sp + r*sp;
                }
            }

            accesser_[0][0][0] = new T***[ss*sp*sr*sa];
            if (accesser_[0][0][0] == NULL)
            {
                delete [] accesser_[0][0];
                delete [] accesser_[0];
                delete [] accesser_;
                return false;
            }

            for (a=0; a<sa; a++)
            {
                for (r=0; r<sr; r++)
                {
                    for (p=0; p<sp; p++)
                    {
                        accesser_[a][r][p] = accesser_[0][0][0] 
                                                + a*sr*sp*ss 
                                                + r*sp*ss 
                                                + p*ss;
                    }
                }
            }

            accesser_[0][0][0][0] = new T**[sz*ss*sp*sr*sa];
            if (accesser_[0][0][0][0] == NULL)
            {
                delete [] accesser_[0][0][0];
                delete [] accesser_[0][0];
                delete [] accesser_[0];
                delete [] accesser_;
                return false;
            }

            for (a=0; a<sa; a++)
            {
                for (r=0; r<sr; r++)
                {
                    for (p=0; p<sp; p++)
                    {
                        for (s=0; s<ss; s++)
                        {
                            accesser_[a][r][p][s] = accesser_[0][0][0][0] 
                                                        + a*sr*sp*ss*sz 
                                                        + r*sp*ss*sz 
                                                        + p*ss*sz 
                                                        + s*sz;
                        }
                    }
                }
            }

            accesser_[0][0][0][0][0] = new T*[sy*sz*ss*sp*sr*sa];
            for (a=0; a<sa; a++)
            {
                for (r=0; r<sr; r++)
                {
                    for (p=0; p<sp; p++)
                    {
                        for (s=0; s<ss; s++)
                        {
                            for (z=0; z<sz; z++)
                            {
                                accesser_[a][r][p][s][z] = accesser_[0][0][0][0][0] 
                                                                + a*sr*sp*ss*sz*sy 
                                                                + r*sp*ss*sz*sy 
                                                                + p*ss*sz*sy 
                                                                + s*sz*sy 
                                                                + z*sy;
                            }
                        }
                    }
                }
            }

            accesser_[0][0][0][0][0][0] = data_;
            for (a=0; a<sa; a++)
            {
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
                                    accesser_[a][r][p][s][z][y] = accesser_[0][0][0][0][0][0] 
                                                                    + a*sr*sp*ss*sz*sy*sx 
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
        }
        else
        {
            accesser_ = NULL;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in ho7DArray<T>::init_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
bool ho7DArray<T>::release_accesser()
{
    try
    {
        if (accesser_ != NULL)
        {
            delete [] accesser_[0][0][0][0][0];
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
        GERROR_STREAM("Errors in ho7DArray<T>::release_accesser() ...");
        return false;
    }

    return true;
}

template <typename T> 
void ho7DArray<T>::print(std::ostream& os) const
{
    BaseClass::print(os);
    size_t x, y, z, s, p, r, a;
    os << "-------------------------------------------" << std::endl;
    for (a=0; a<dimensions_[6]; a++) 
    {
        for (r=0; r<dimensions_[5]; r++) 
        {
            for (p=0; p<dimensions_[4]; p++) 
            {
                for (s=0; s<dimensions_[3]; s++) 
                {
                    for (z=0; z<dimensions_[2]; z++) 
                    {
                        os << "ho7DArray (:, :, " << z << ", " << s << ", " << p << ", " << r << ", " << a << ") = " << std::endl;
                        for (y=0; y<dimensions_[1]; y++) 
                        {
                            os << "y " << y << "\t";
                            for (x=0; x<dimensions_[0]; x++)
                            {
                                os << (*this)(x,y,z,s,p,r,a) << "\t";
                            }
                            os << std::endl;
                        }
                    }
                }
            }
        }
    }
    os << "-------------------------------------------" << std::endl;
}

}
