// This file is not to be included by anyone else than hoNDArray.h
// Contains the "private" implementation of the container
//

namespace Gadgetron
{
    template <typename T> 
    hoNDArray<T>::hoNDArray() : NDArray<T>::NDArray() 
    {
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(std::vector<size_t> *dimensions) : NDArray<T>::NDArray()
    {
        this->create(dimensions);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(std::vector<size_t> &dimensions) : NDArray<T>::NDArray()
    {
        this->create(dimensions);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(boost::shared_ptr< std::vector<size_t> > dimensions) : NDArray<T>::NDArray()
    {
        this->create(dimensions);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t len) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(1);
        dim[0] = len;
        this->create(dim);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(2);
        dim[0] = sx;
        dim[1] = sy;
        this->create(dim);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(3);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        this->create(dim);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(4);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        this->create(dim);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(5);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        this->create(dim);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(6);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;
        this->create(dim);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(7);
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
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(8);
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
    hoNDArray<T>::hoNDArray(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        this->create(dimensions,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        this->create(dimensions,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        this->create(dimensions,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t len, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(1);
        dim[0] = len;
        this->create(&dim,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(2);
        dim[0] = sx;
        dim[1] = sy;
        this->create(&dim,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(3);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        this->create(&dim,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(4);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        this->create(&dim,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(5);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        this->create(&dim,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(6);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;
        this->create(&dim,data,delete_data_on_destruct);
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(7);
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
    hoNDArray<T>::hoNDArray(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct) : NDArray<T>::NDArray()
    {
        std::vector<size_t> dim(8);
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
    hoNDArray<T>::~hoNDArray()
    {
        if (this->delete_data_on_destruct_){
            deallocate_memory();
        }
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(const hoNDArray<T>  *a)
    {
        if(!a) throw std::runtime_error("hoNDArray<T>::hoNDArray(): 0x0 pointer provided");
        this->data_ = 0;

        std::vector<size_t>* tmp = new std::vector<size_t>;
        this->dimensions_ = boost::shared_ptr< std::vector<size_t> >(tmp);
        *(this->dimensions_) = *(a->dimensions_);

        tmp = new std::vector<size_t>;
        this->offsetFactors_ = boost::shared_ptr< std::vector<size_t> >(tmp);
        *(this->offsetFactors_) = *(a->offsetFactors_);

        if ( !this->dimensions_->empty() )
        {
            allocate_memory();
            memcpy( this->data_, a->data_, this->elements_*sizeof(T) );
        }
        else
        {
            this->elements_ = 0;
        }
    }

    template <typename T> 
    hoNDArray<T>::hoNDArray(const hoNDArray<T> &a)
    {
        this->data_ = 0;

        std::vector<size_t>* tmp = new std::vector<size_t>;
        this->dimensions_ = boost::shared_ptr< std::vector<size_t> >(tmp);
        *(this->dimensions_) = *(a.dimensions_);

        tmp = new std::vector<size_t>;
        this->offsetFactors_ = boost::shared_ptr< std::vector<size_t> >(tmp);
        *(this->offsetFactors_) = *(a.offsetFactors_);

        if ( !this->dimensions_->empty() )
        {
            allocate_memory();
            memcpy( this->data_, a.data_, this->elements_*sizeof(T) );
        }
        else
        {
            this->elements_ = 0;
        }
    }

    template <typename T> 
    hoNDArray<T>& hoNDArray<T>::operator=(const hoNDArray<T>& rhs)
    {
        if ( &rhs == this ) return *this;

        if ( rhs.get_number_of_elements() == 0 ){
            this->clear();
            return *this;
        }

        // Are the dimensions the same? Then we can just memcpy
        if (this->dimensions_equal(&rhs)){
            memcpy(this->data_, rhs.data_, this->elements_*sizeof(T));
        }
        else{
            deallocate_memory();
            this->data_ = 0;
            *(this->dimensions_) = *(rhs.dimensions_);
            *(this->offsetFactors_) = *(rhs.offsetFactors_);
            allocate_memory();
            memcpy( this->data_, rhs.data_, this->elements_*sizeof(T) );
        }
        return *this;
    }

    template <typename T> 
    void hoNDArray<T>::create(std::vector<size_t>& dimensions)
    {
        if ( this->dimensions_equal(&dimensions) )
        {
            return;
        }

        this->clear();
        BaseClass::create(dimensions);
    }

    template <typename T> 
    void hoNDArray<T>::create(std::vector<size_t> *dimensions)
    {
        if ( this->dimensions_equal(dimensions) )
        {
            return;
        }
        this->clear();
        BaseClass::create(dimensions);
    }

    template <typename T> 
    void hoNDArray<T>::create(boost::shared_ptr< std::vector<size_t> > dimensions)
    {
        if ( this->dimensions_equal(dimensions.get()) )
        {
            return;
        }
        this->clear();
        BaseClass::create(dimensions);
    }

    template <typename T> 
    void hoNDArray<T>::create(std::vector<size_t> *dimensions, T* data, bool delete_data_on_destruct) 
    {
        if(!dimensions) throw std::runtime_error("hoNDArray<T>::create(): 0x0 pointer provided");
        if(!data) throw std::runtime_error("hoNDArray<T>::create(): 0x0 pointer provided");

        if ( this->dimensions_equal(dimensions) )
        {
            if ( this->delete_data_on_destruct_ )
            {
                this->deallocate_memory();
            }

            this->data_ = data;
            this->delete_data_on_destruct_ = delete_data_on_destruct;
        }
        else
        {
            if ( this->delete_data_on_destruct_ )
            {
                this->deallocate_memory();
                this->data_ = NULL;
            }

            BaseClass::create(dimensions, data, delete_data_on_destruct);
        }
    }

    template <typename T> 
    void hoNDArray<T>::create(std::vector<size_t> &dimensions, T* data, bool delete_data_on_destruct) 
    {
        if(!data) throw std::runtime_error("hoNDArray<T>::create(): 0x0 pointer provided");

        if ( this->dimensions_equal(&dimensions) )
        {
            if ( this->delete_data_on_destruct_ )
            {
                this->deallocate_memory();
            }

            this->data_ = data;
            this->delete_data_on_destruct_ = delete_data_on_destruct;
        }
        else
        {
            if ( this->delete_data_on_destruct_ )
            {
                this->deallocate_memory();
                this->data_ = NULL;
            }

            BaseClass::create(dimensions, data, delete_data_on_destruct);
        }
    }

    template <typename T> 
    inline void hoNDArray<T>::create(boost::shared_ptr< std::vector<size_t> > dimensions, T* data, bool delete_data_on_destruct)
    {
        this->create(dimensions.get(), data, delete_data_on_destruct);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t len)
    {
        std::vector<size_t> dim(1);
        dim[0] = len;
        this->create(dim);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy)
    {
        std::vector<size_t> dim(2);
        dim[0] = sx;
        dim[1] = sy;
        this->create(dim);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz)
    {
        std::vector<size_t> dim(3);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        this->create(dim);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st)
    {
        std::vector<size_t> dim(4);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        this->create(dim);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp)
    {
        std::vector<size_t> dim(5);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        this->create(dim);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq)
    {
        std::vector<size_t> dim(6);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;
        this->create(dim);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr)
    {
        std::vector<size_t> dim(7);
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
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss)
    {
        std::vector<size_t> dim(8);
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
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su)
    {
        std::vector<size_t> dim(9);
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
    inline void hoNDArray<T>::create(size_t len, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(1);
        dim[0] = len;
        this->create(&dim, data, delete_data_on_destruct);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(2);
        dim[0] = sx;
        dim[1] = sy;
        this->create(&dim, data, delete_data_on_destruct);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(3);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        this->create(&dim, data, delete_data_on_destruct);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(4);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        this->create(&dim, data, delete_data_on_destruct);
    }
    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(5);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        this->create(&dim, data, delete_data_on_destruct);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(6);
        dim[0] = sx;
        dim[1] = sy;
        dim[2] = sz;
        dim[3] = st;
        dim[4] = sp;
        dim[5] = sq;
        this->create(&dim, data, delete_data_on_destruct);
    }

    template <typename T> 
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(7);
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
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(8);
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
    inline void hoNDArray<T>::create(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, size_t su, T* data, bool delete_data_on_destruct)
    {
        std::vector<size_t> dim(9);
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
    inline T& hoNDArray<T>::at( size_t idx )
    {
        /*if( idx >= this->get_number_of_elements() )
        {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::at(): index out of range."));
        }*/
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->get_data_ptr()[idx];
    }

    template <typename T> 
    inline const T& hoNDArray<T>::at( size_t idx ) const
    {
        /*if( idx >= this->get_number_of_elements() )
        {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::at(): index out of range."));
        }*/
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->get_data_ptr()[idx];
    }

    template <typename T> 
    inline T& hoNDArray<T>::operator[]( size_t idx )
    {
        /*if( idx >= this->get_number_of_elements() )
        {
        BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator[]: index out of range."));
        }*/
        GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
        return this->get_data_ptr()[idx];
    }

    //template <typename T> 
    //inline T& hoNDArray<T>::operator()( size_t idx )
    //{
    //    /*if( idx >= this->get_number_of_elements() )
    //    {
    //    BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    //    }*/
    //    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    //    return this->get_data_ptr()[idx];
    //}

    //template <typename T> 
    //inline const T& hoNDArray<T>::operator()( size_t idx ) const
    //{
    //    /*if( idx >= this->get_number_of_elements() )
    //    {
    //    BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    //    }*/
    //    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    //    return this->get_data_ptr()[idx];
    //}

    //template <typename T> 
    //inline T& hoNDArray<T>::operator()( const std::vector<size_t>& ind )
    //{
    //    size_t idx = this->calculate_offset(ind);
    //    /*if( idx >= this->get_number_of_elements() )
    //    {
    //    BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    //    }*/
    //    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    //    return this->get_data_ptr()[idx];
    //}

    //template <typename T> 
    //inline const T& hoNDArray<T>::operator()( const std::vector<size_t>& ind ) const
    //{
    //    size_t idx = this->calculate_offset(ind);
    //    /*if( idx >= this->get_number_of_elements() )
    //    {
    //    BOOST_THROW_EXCEPTION( runtime_error("hoNDArray::operator(): index out of range."));
    //    }*/
    //    GADGET_DEBUG_CHECK_THROW(idx < this->get_number_of_elements());
    //    return this->get_data_ptr()[idx];
    //}

    template <typename T> 
    void hoNDArray<T>::get_sub_array(const std::vector<size_t>& start, std::vector<size_t>& size, hoNDArray<T>& out)
    {
        if ( start.size() != size.size() ){
            BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
        }

        if ( start.size() != (*dimensions_).size() ){
            BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
        }

        out.create(&size);

        if ( out.get_number_of_elements() == this->get_number_of_elements() ){
            out = *this;
            return;
        }

        std::vector<size_t> end(start.size());

        size_t ii;
        for ( ii=0; ii<start.size(); ii++ ){
            end[ii] = start[ii] + size[ii] - 1;
            if ( end[ii] >= (*dimensions_)[ii] ){
                BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
            }
        }
    }

    template <typename T> 
    void hoNDArray<T>::printContent(std::ostream& os) const
    {
        using namespace std;

        os.unsetf(std::ios::scientific);
        os.setf(ios::fixed);

        size_t i;

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
        os << "Delete data on destruction flag is : " << this->delete_data_on_destruct_ << endl;

        //os << "-------------------------------------------" << std::endl;
        //size_t numOfPrints = 20;
        //if ( this->elements_ < numOfPrints ) numOfPrints = this->elements_;
        //for (i=0; i<numOfPrints; i++) 
        //{
        //    os << i << " = " << (*this)(i) << std::endl;
        //}
        //os << "-------------------------------------------" << std::endl;

        os << std::endl;
    }

    template <typename T> 
    void hoNDArray<T>::print(std::ostream& os) const
    {
        using namespace std;

        os.unsetf(std::ios::scientific);
        os.setf(ios::fixed);

        os << "--------------Gagdgetron ND Array -------------" << endl;
        this->printContent(os);
    }

    template <typename T> 
    void hoNDArray<T>::allocate_memory()
    {
        deallocate_memory();

        if ( !this->dimensions_->empty() )
        {
            this->elements_ = (*this->dimensions_)[0];
            for (size_t i = 1; i < this->dimensions_->size(); i++)
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
        else
        {
            this->elements_ = 0;
        }
    }

    template <typename T> 
    void hoNDArray<T>::deallocate_memory()
    {
        if (!(this->delete_data_on_destruct_)) {
             throw std::runtime_error("You don't own this data.  You cannot deallocate its memory.");
        }
        
        if( this->data_ ){
            this->_deallocate_memory( this->data_ );
            this->data_ = 0x0;
        }
    }

    template <typename T> 
    inline void hoNDArray<T>::_allocate_memory( size_t size, float** data )
    {
        *data = (float*) malloc( size*sizeof(float) );
    }

    template <typename T> 
    inline void hoNDArray<T>::_deallocate_memory( float* data )
    {
        free(data);
    }

    template <typename T> 
    inline void hoNDArray<T>::_allocate_memory( size_t size, double** data )
    {
        *data = (double*) malloc( size*sizeof(double) );
    }

    template <typename T> 
    inline void hoNDArray<T>::_deallocate_memory( double* data )
    {
        free(data);
    }

    template <typename T> 
    inline void hoNDArray<T>::_allocate_memory( size_t size, std::complex<float>** data )
    {
        *data = (std::complex<float>*) malloc( size*sizeof(std::complex<float>) );
    }

    template <typename T> 
    inline void hoNDArray<T>::_deallocate_memory( std::complex<float>* data )
    {
        free(data);
    }

    template <typename T> 
    inline void hoNDArray<T>::_allocate_memory( size_t size, std::complex<double>** data )
    {
        *data = (std::complex<double>*) malloc( size*sizeof(std::complex<double>) );
    }

    template <typename T> 
    inline void hoNDArray<T>::_deallocate_memory( std::complex<double>* data )
    {
        free(data);
    }

    template <typename T> 
    inline void hoNDArray<T>::_allocate_memory( size_t size, float_complext** data )
    {
        *data = (float_complext*) malloc( size*sizeof(float_complext) );
    }

    template <typename T> 
    inline void hoNDArray<T>::_deallocate_memory( float_complext* data )
    {
        free( data );
    }

    template <typename T> 
    inline void hoNDArray<T>::_allocate_memory( size_t size, double_complext** data )
    {
        *data = (double_complext*) malloc( size*sizeof(double_complext) );
    }

    template <typename T> 
    inline void hoNDArray<T>::_deallocate_memory( double_complext* data )
    {
        free( data );
    }

    template <typename T> 
    bool hoNDArray<T>::serialize(char*& buf, size_t& len) const 
    {
        if ( buf != NULL ) delete[] buf;

        size_t NDim = dimensions_->size();

        // number of dimensions + dimension vector + contents
        len = sizeof(size_t) + sizeof(size_t)*NDim + sizeof(T)*elements_;

        buf = new char[len];

        memcpy(buf, &NDim, sizeof(size_t));
        if ( NDim > 0 )
        {
            memcpy(buf+sizeof(size_t), &((*dimensions_)[0]), sizeof(size_t)*NDim);
            memcpy(buf+sizeof(size_t)+sizeof(size_t)*NDim, this->data_, sizeof(T)*elements_);
        }

        return true; // Temporary. Should not be a boolean function.
    }

    template <typename T> 
    bool hoNDArray<T>::deserialize(char* buf, size_t& len)
    {
        size_t NDim;
        memcpy(&NDim, buf, sizeof(size_t));

        if ( NDim > 0 )
        {
            std::vector<size_t> dimensions(NDim);
            memcpy(&dimensions[0], buf+sizeof(size_t), sizeof(size_t)*NDim);

            // allocate memory
            this->create(&dimensions);

            // copy the content
            memcpy(this->data_, buf+sizeof(size_t)+sizeof(size_t)*NDim, sizeof(T)*elements_);
        }
        else
        {
            this->clear();
        }

        len = sizeof(size_t)+sizeof(size_t)*NDim+sizeof(T)*elements_;
        return true; // Temporary. Should not be a boolean function.
    }
}
