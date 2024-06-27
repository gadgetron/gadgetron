/** \file       hoNDImage.hxx
    \brief      Implementation of N-dimensional image class for gadgetron
    \author     Hui Xue
*/

namespace Gadgetron
{
    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage () : BaseClass()
    {
        dimensions_.resize(D, 0);

        unsigned int ii;
        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = 1;
            pixelSize_reciprocal_[ii] = 1;
            origin_[ii] = 0;
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (boost::shared_ptr< std::vector<size_t> > dimensions) : BaseClass( dimensions )
    {
        this->create( *dimensions );
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions, 
        const std::vector<coord_type>& pixelSize) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions, pixelSize);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions, 
        const std::vector<coord_type>& pixelSize, 
        const std::vector<coord_type>& origin) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions, pixelSize, origin);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions, 
                                const std::vector<coord_type>& pixelSize, 
                                const std::vector<coord_type>& origin, 
                                const axis_type& axis) : BaseClass( const_cast<std::vector<size_t>& >(dimensions) )
    {
        this->create(dimensions, pixelSize, origin, axis);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t len) : BaseClass(len)
    {
        std::vector<size_t> dimension(1, len);
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy) : BaseClass(sx, sy)
    {
        std::vector<size_t> dimension(2);
        dimension[0] = sx;
        dimension[1] = sy;
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz) : BaseClass(sx, sy, sz)
    {
        std::vector<size_t> dimension(3);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st) : BaseClass(sx, sy, sz, st)
    {
        std::vector<size_t> dimension(4);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp) : BaseClass(sx, sy, sz, st, sp)
    {
        std::vector<size_t> dimension(5);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq) : BaseClass(sx, sy, sz, st, sp, sq)
    {
        std::vector<size_t> dimension(6);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr) : BaseClass(sx, sy, sz, st, sp, sq, sr)
    {
        std::vector<size_t> dimension(7);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss) : BaseClass(sx, sy, sz, st, sp, sq, sr, ss)
    {
        std::vector<size_t> dimension(8);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        dimension[7] = ss;
        this->create(dimension);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct) : BaseClass(const_cast<std::vector<size_t>*>(&dimensions), data, delete_data_on_destruct)
    {
        this->create(dimensions, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct) : BaseClass(const_cast<std::vector<size_t>*>(&dimensions), data, delete_data_on_destruct)
    {
        this->create(dimensions, pixelSize, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct) : BaseClass(const_cast<std::vector<size_t>*>(&dimensions), data, delete_data_on_destruct)
    {
        this->create(dimensions, pixelSize, origin, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct) : BaseClass(const_cast<std::vector<size_t>*>(&dimensions), data, delete_data_on_destruct)
    {
        this->create(dimensions, pixelSize, origin, axis, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t len, T* data, bool delete_data_on_destruct) : BaseClass(len, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(1, len);
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(2);
        dimension[0] = sx;
        dimension[1] = sy;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(3);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(4);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(5);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(6);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, sr, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(7);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, sr, ss, data, delete_data_on_destruct)
    {
        std::vector<size_t> dimension(8);
        dimension[0] = sx;
        dimension[1] = sy;
        dimension[2] = sz;
        dimension[3] = st;
        dimension[4] = sp;
        dimension[5] = sq;
        dimension[6] = sr;
        dimension[7] = ss;
        this->create(dimension, data, delete_data_on_destruct);
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(const hoNDArray<T>& a) : BaseClass(a)
    {
         boost::shared_ptr< std::vector<size_t> > dim = a.get_dimensions();
         this->create(*dim);
         memcpy(this->data_, a.begin(), this->get_number_of_bytes());
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::hoNDImage(const Self& a) : BaseClass()
    {
        *this = a;
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>& hoNDImage<T, D>::operator=(const Self& rhs)
    {
        if ( &rhs == this ) return *this;

        if ( rhs.get_number_of_elements() == 0 )
        {
            this->clear();
            return *this;
        }

        if ( this->dimensions_equal(rhs) && this->data_!=NULL )
        {
            memcpy(this->data_, rhs.data_, rhs.elements_*sizeof(T));
        }
        else
        {
            this->deallocate_memory();
            this->data_ = 0;

            this->dimensions_ = rhs.dimensions_;
            this->allocate_memory();
            this->calculate_offset_factors( this->dimensions_ );
            memcpy( this->data_, rhs.data_, this->elements_*sizeof(T) );
        }

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->pixelSize_[ii] = rhs.pixelSize_[ii];
            this->pixelSize_reciprocal_[ii] = rhs.pixelSize_reciprocal_[ii];
            this->origin_[ii] = rhs.origin_[ii];
            this->axis_[ii] = rhs.axis_[ii];
        }

        this->image_position_patient_ = rhs.image_position_patient_;
        this->image_orientation_patient_[0] = rhs.image_orientation_patient_[0];
        this->image_orientation_patient_[1] = rhs.image_orientation_patient_[1];
        this->image_orientation_patient_[2] = rhs.image_orientation_patient_[2];

        return *this;
    }

    template <typename T, unsigned int D> 
    hoNDImage<T, D>::~hoNDImage()
    {
        if (this->delete_data_on_destruct_)
        {
            this->deallocate_memory();
        }
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::clear()
    {
        if ( this->delete_data_on_destruct_ )
        {
            this->deallocate_memory();
        }
        this->data_ = 0;
        this->elements_ = 0;
        this->delete_data_on_destruct_ = true;

        unsigned int ii;

        dimensions_.clear();

        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = 1;
            pixelSize_reciprocal_[ii] = 1;
            origin_[ii] = 0;
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions)
    {
        if ( !this->dimensions_equal(dimensions) )
        {
            dimensions_ = dimensions;
            this->allocate_memory();
            this->calculate_offset_factors(dimensions);
        }

        unsigned int ii;
        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = 1;
            pixelSize_reciprocal_[ii] = 1;
            origin_[ii] = 0;
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(boost::shared_ptr< std::vector<size_t> > dimensions)
    {
        this->create(*dimensions);
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize)
    {
        if ( !this->dimensions_equal(dimensions) )
        {
            dimensions_ = dimensions;
            this->allocate_memory();
            this->calculate_offset_factors(dimensions);
        }

        unsigned int ii;
        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = pixelSize[ii];
            pixelSize_reciprocal_[ii] = coord_type(1.0)/pixelSize_[ii];
            origin_[ii] = 0;
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin)
    {
        if ( !this->dimensions_equal(dimensions) )
        {
            dimensions_ = dimensions;
            this->allocate_memory();
            this->calculate_offset_factors(dimensions);
        }

        unsigned int ii;
        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = pixelSize[ii];
            pixelSize_reciprocal_[ii] = coord_type(1.0)/pixelSize_[ii];
            origin_[ii] = origin[ii];
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        if ( D==1 )
        {
            image_position_patient_[0] = origin[0];
        }
        else if ( D == 2 )
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
        }
        else
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
            image_position_patient_[2] = origin[2];
        }

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis)
    {
        if ( !this->dimensions_equal(dimensions) )
        {
            dimensions_ = dimensions;
            this->allocate_memory();
            this->calculate_offset_factors(dimensions);
        }

        unsigned int ii;
        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = pixelSize[ii];
            pixelSize_reciprocal_[ii] = coord_type(1.0)/pixelSize_[ii];
            origin_[ii] = origin[ii];
            axis_[ii] = axis[ii];
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;

        if ( D==1 )
        {
            image_position_patient_[0] = origin[0];
        }
        else if ( D == 2 )
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
        }
        else
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
            image_position_patient_[2] = origin[2];

            image_orientation_patient_[0][0] = axis[0][0]; image_orientation_patient_[0][1] = axis[0][1]; image_orientation_patient_[0][2] = axis[0][2];
            image_orientation_patient_[1][0] = axis[1][0]; image_orientation_patient_[1][1] = axis[1][1]; image_orientation_patient_[1][2] = axis[1][2];
            image_orientation_patient_[2][0] = axis[2][0]; image_orientation_patient_[2][1] = axis[2][1]; image_orientation_patient_[2][2] = axis[2][2];
        }
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct)
    {
        if ( this->delete_data_on_destruct_ )
        {
            this->deallocate_memory();
            this->data_ = NULL;
        }

        this->data_ = data;
        this->delete_data_on_destruct_ = delete_data_on_destruct;
        this->dimensions_ = dimensions;

        unsigned int ii;

        this->elements_ = 1;
        for (ii=0; ii<D; ii++)
        {
            this->elements_ *= dimensions_[ii];
        }
        this->calculate_offset_factors(dimensions);

        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = 1.0;
            pixelSize_reciprocal_[ii] = 1.0;
            origin_[ii] = 0;
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct)
    {
        if ( this->delete_data_on_destruct_ )
        {
            this->deallocate_memory();
            this->data_ = NULL;
        }

        this->data_ = data;
        this->delete_data_on_destruct_ = delete_data_on_destruct;
        this->dimensions_ = dimensions;

        unsigned int ii;

        this->elements_ = 1;
        for (ii=0; ii<D; ii++)
        {
            this->elements_ *= dimensions_[ii];
        }
        this->calculate_offset_factors(dimensions);

        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = pixelSize[ii];
            pixelSize_reciprocal_[ii] = coord_type(1.0)/pixelSize_[ii];
            origin_[ii] = 0;
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct)
    {
        if ( this->delete_data_on_destruct_ )
        {
            this->deallocate_memory();
            this->data_ = NULL;
        }

        this->data_ = data;
        this->delete_data_on_destruct_ = delete_data_on_destruct;
        this->dimensions_ = dimensions;

        unsigned int ii;

        this->elements_ = 1;
        for (ii=0; ii<D; ii++)
        {
            this->elements_ *= dimensions_[ii];
        }
        this->calculate_offset_factors(dimensions);

        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = pixelSize[ii];
            pixelSize_reciprocal_[ii] = coord_type(1.0)/pixelSize_[ii];
            origin_[ii] = origin[ii];
            axis_[ii].fill(0);
            axis_[ii][ii] = coord_type(1.0);
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        if ( D==1 )
        {
            image_position_patient_[0] = origin[0];
        }
        else if ( D == 2 )
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
        }
        else
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
            image_position_patient_[2] = origin[2];
        }

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct)
    {
        if ( this->delete_data_on_destruct_ )
        {
            this->deallocate_memory();
            this->data_ = NULL;
        }

        this->data_ = data;
        this->delete_data_on_destruct_ = delete_data_on_destruct;

        dimensions_ = dimensions;

        unsigned int ii;

        this->elements_ = 1;
        for (ii=0; ii<D; ii++)
        {
            this->elements_ *= dimensions_[ii];
        }
        this->calculate_offset_factors(dimensions);

        for (ii=0;ii<D; ii++)
        {
            pixelSize_[ii] = pixelSize[ii];
            pixelSize_reciprocal_[ii] = coord_type(1.0)/pixelSize_[ii];
            origin_[ii] = origin[ii];
            axis_[ii] = axis[ii];
        }

        image_position_patient_[0] = 0;
        image_position_patient_[1] = 0;
        image_position_patient_[2] = 0;

        image_orientation_patient_[0][0] = 1; image_orientation_patient_[0][1] = 0; image_orientation_patient_[0][2] = 0;
        image_orientation_patient_[1][0] = 0; image_orientation_patient_[1][1] = 1; image_orientation_patient_[1][2] = 0;
        image_orientation_patient_[2][0] = 0; image_orientation_patient_[2][1] = 0; image_orientation_patient_[2][2] = 1;

        if ( D==1 )
        {
            image_position_patient_[0] = origin[0];
        }
        else if ( D == 2 )
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
        }
        else
        {
            image_position_patient_[0] = origin[0];
            image_position_patient_[1] = origin[1];
            image_position_patient_[2] = origin[2];

            image_orientation_patient_[0][0] = axis[0][0]; image_orientation_patient_[0][1] = axis[0][1]; image_orientation_patient_[0][2] = axis[0][2];
            image_orientation_patient_[1][0] = axis[1][0]; image_orientation_patient_[1][1] = axis[1][1]; image_orientation_patient_[1][2] = axis[1][2];
            image_orientation_patient_[2][0] = axis[2][0]; image_orientation_patient_[2][1] = axis[2][1]; image_orientation_patient_[2][2] = axis[2][2];
        }
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::from_NDArray(const hoNDArray<T>& a)
    {
        boost::shared_ptr< std::vector<size_t> > dim = a.get_dimensions();

        size_t ii;

        if ( dim->size() < D )
        {
            std::vector<size_t> dimUsed(D, 1);
            for ( ii=0; ii<dim->size(); ii++ )
            {
                dimUsed[ii] = (*dim)[ii];
            }

            if ( !this->dimensions_equal(dimUsed) )
            {
                this->create(dimUsed);
            }
        }
        else if ( dim->size() > D )
        {
            std::vector<size_t> dimUsed(D, 1);
            for ( ii=0; ii<D; ii++ )
            {
                dimUsed[ii] = (*dim)[ii];
            }

            if ( !this->dimensions_equal(dimUsed) )
            {
                this->create(dimUsed);
            }
        }
        else
        {
            if ( !this->dimensions_equal(*dim) )
            {
                this->create(*dim);
            }
        }

        memcpy(this->data_, a.begin(), this->get_number_of_bytes());
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::to_NDArray(hoNDArray<T>& a) const
    {
        std::vector<size_t> dim;
        this->get_dimensions(dim);

        if ( !a.dimensions_equal(dim) )
        {
            a.create(&dim);
        }

        memcpy(a.begin(), this->data_, a.get_number_of_bytes());
    }

    template <typename T, unsigned int D> 
    inline bool hoNDImage<T, D>::dimensions_equal(const std::vector<size_t>& dimensions) const
    {
        if ( (dimensions.size() != D) || ( dimensions_.size() != dimensions.size() ) ) return false;

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            if ( dimensions_[ii] != dimensions[ii] ) return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    inline typename hoNDImage<T, D>::coord_type hoNDImage<T, D>::get_pixel_size(size_t dimension) const
    {
        GADGET_DEBUG_CHECK_THROW(dimension < D);
        return this->pixelSize_[dimension];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_pixel_size(std::vector<coord_type>& pixelSize) const
    {
        pixelSize.resize(D);
        memcpy(&pixelSize[0], this->pixelSize_, sizeof(coord_type)*D);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_pixel_size(size_t dimension, coord_type v)
    {
        GADGET_DEBUG_CHECK_THROW(dimension < D);
        this->pixelSize_[dimension] = v;
        this->pixelSize_reciprocal_[dimension] = coord_type(1.0)/v;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_pixel_size(const std::vector<coord_type>& pixelSize)
    {
        GADGET_DEBUG_CHECK_THROW(pixelSize.size() >= D);
        memcpy(this->pixelSize_, &pixelSize[0], sizeof(coord_type)*D);

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->pixelSize_reciprocal_[ii] = coord_type(1.0)/this->pixelSize_[ii];
        }
    }

    template <typename T, unsigned int D> 
    inline typename hoNDImage<T, D>::coord_type hoNDImage<T, D>::get_origin(size_t dimension) const
    {
        GADGET_DEBUG_CHECK_THROW(dimension < D);
        return this->origin_[dimension];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_origin(std::vector<coord_type>& origin) const
    {
        origin.resize(D);
        memcpy(&origin[0], this->origin_, sizeof(coord_type)*D);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_origin(size_t dimension, coord_type v)
    {
        GADGET_DEBUG_CHECK_THROW(dimension < D);
        this->origin_[dimension] = v;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_origin(const std::vector<coord_type>& origin)
    {
        GADGET_DEBUG_CHECK_THROW(origin.size() >= D);
        memcpy(this->origin_, &origin[0], sizeof(coord_type)*D);
    }

    template <typename T, unsigned int D> 
    inline typename hoNDImage<T, D>::coord_type hoNDImage<T, D>::get_axis(size_t dimension, size_t elem) const
    {
        GADGET_DEBUG_CHECK_THROW(dimension<D && elem<D);
        return this->axis_[dimension][elem];
    }

    template <typename T, unsigned int D> 
    inline typename hoNDImage<T, D>::a_axis_type hoNDImage<T, D>::get_axis(size_t dimension) const
    {
        GADGET_DEBUG_CHECK_THROW(dimension < D);
        return this->axis_[dimension];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_axis(axis_type& axis) const
    {
        axis.resize(D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            axis[ii] = this->axis_[ii];
        }
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_axis(size_t dimension, size_t elem, coord_type v)
    {
        GADGET_DEBUG_CHECK_THROW(dimension<D && elem<D);
        this->axis_[dimension][elem] = v;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_axis(size_t dimension, const a_axis_type& v)
    {
        GADGET_DEBUG_CHECK_THROW(dimension < D);
        this->axis_[dimension] = v;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_axis(const axis_type& axis)
    {
        GADGET_DEBUG_CHECK_THROW(axis.size() >= D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            this->axis_[ii] = axis[ii];
        }
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_image_position(coord_type pos[3]) const 
    {
        pos[0] = image_position_patient_[0];
        pos[1] = image_position_patient_[1];
        pos[2] = image_position_patient_[2];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_image_position(unsigned int d, coord_type& pos) const 
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        pos = image_position_patient_[d];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_image_position(a_axis_image_patient_type& pos) const 
    {
        pos = image_position_patient_;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_position(coord_type pos[3])
    {
        image_position_patient_[0] = pos[0];
        image_position_patient_[1] = pos[1];
        image_position_patient_[2] = pos[2];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_position(unsigned int d, coord_type pos)
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        image_position_patient_[d] = pos;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_position(const a_axis_image_patient_type& pos)
    {
        image_position_patient_ = pos;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_image_orientation(unsigned int d, coord_type ori[3]) const 
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        ori[0] = image_orientation_patient_[d][0];
        ori[1] = image_orientation_patient_[d][1];
        ori[2] = image_orientation_patient_[d][2];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_image_orientation(unsigned int d, a_axis_image_patient_type& ori) const 
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        ori = image_orientation_patient_[d];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_image_orientation(unsigned int d, unsigned int ind, coord_type& ori) const 
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        GADGET_DEBUG_CHECK_THROW(ind<3);
        ori = image_orientation_patient_[d][ind];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::get_image_orientation(coord_type quat[4]) const 
    {
        coord_type r11 = image_orientation_patient_[0][0], r12 = image_orientation_patient_[1][0], r13 = image_orientation_patient_[2][0];
        coord_type r21 = image_orientation_patient_[0][1], r22 = image_orientation_patient_[1][1], r23 = image_orientation_patient_[2][1];
        coord_type r31 = image_orientation_patient_[0][2], r32 = image_orientation_patient_[1][2], r33 = image_orientation_patient_[2][2];

        double a = 1, b = 0, c = 0, d = 0, s = 0;
        double trace = 0;
        double xd, yd, zd;

        /* verify the sign of the rotation*/
        coord_type deti = (r11 * r22 * r33) + (r12 * r23 * r31) + (r21 * r32 * r13) -
            (r13 * r22 * r31) - (r12 * r21 * r33) - (r11 * r23 * r32);

        if (deti < 0)
        {
            /* flip 3rd column */
            r13 = -r13;
            r23 = -r23;
            r33 = -r33;
        }

        /* Compute quaternion parameters */
        /* http://www.cs.princeton.edu/~gewang/projects/darth/stuff/quat_faq.html#Q55 */
        trace = 1.0l + r11 + r22 + r33;
        if (trace > 0.00001l)
        {                /* simplest case */
            s = std::sqrt(trace) * 2;
            a = (r32 - r23) / s;
            b = (r13 - r31) / s;
            c = (r21 - r12) / s;
            d = 0.25l * s;
        }
        else
        {
            /* trickier case...
             * determine which major diagonal element has
             * the greatest value... */
            xd = 1.0 + r11 - (r22 + r33);  /* 4**b**b */
            yd = 1.0 + r22 - (r11 + r33);  /* 4**c**c */
            zd = 1.0 + r33 - (r11 + r22);  /* 4**d**d */
            /* if r11 is the greatest */
            if (xd > 1.0)
            {
                s = 2.0 * std::sqrt(xd);
                a = 0.25l * s;
                b = (r21 + r12) / s;
                c = (r31 + r13) / s;
                d = (r32 - r23) / s;
            }
            /* else if r22 is the greatest */
            else if (yd > 1.0)
            {
                s = 2.0 * std::sqrt(yd);
                a = (r21 + r12) / s;
                b = 0.25l * s;
                c = (r32 + r23) / s;
                d = (r13 - r31) / s;
            }
            /* else, r33 must be the greatest */
            else
            {
                s = 2.0 * std::sqrt(zd);
                a = (r13 + r31) / s;
                b = (r23 + r32) / s;
                c = 0.25l * s;
                d = (r21 - r12) / s;
            }

            if (a < 0.0l)
            {
                b = -b;
                c = -c;
                d = -d;
                a = -a;
            }
        }

        quat[0] = (coord_type)a; 
        quat[1] = (coord_type)b; 
        quat[2] = (coord_type)c; 
        quat[3] = (coord_type)d;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_orientation(unsigned int d, coord_type ori[3])
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        image_orientation_patient_[d][0] = ori[0];
        image_orientation_patient_[d][1] = ori[1];
        image_orientation_patient_[d][2] = ori[2];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_orientation(unsigned int d, const a_axis_image_patient_type& ori)
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        image_orientation_patient_[d] = ori;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_orientation(unsigned int d, unsigned int ind, coord_type ori)
    {
        GADGET_DEBUG_CHECK_THROW(d<3);
        GADGET_DEBUG_CHECK_THROW(ind<3);
        image_orientation_patient_[d][ind] = ori;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_orientation(coord_type quat[4])
    {
        coord_type a = quat[0], b = quat[1], c = quat[2], d = quat[3];

        image_orientation_patient_[0][0] = 1 - 2*( b*b + c*c );
        image_orientation_patient_[1][0] = 2*( a*b - c*d );
        image_orientation_patient_[2][0] = 2*( a*c + b*d );

        image_orientation_patient_[0][1] = 2*( a*b + c*d );
        image_orientation_patient_[1][1] = 1 - 2*( a*a + c*c );
        image_orientation_patient_[2][1] = 2*( b*c - a*d );

        image_orientation_patient_[0][2] = 2*( a*c - b*d );
        image_orientation_patient_[1][2] = 2*( b*c + a*d );
        image_orientation_patient_[2][2] = 1 - 2*( a*a + b*b );
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(const size_t* ind) const
    {
        GADGET_DEBUG_CHECK_THROW(ind!=NULL);

        size_t offset = ind[0];
        for( size_t i = 1; i < D; i++ )
            offset += ind[i] * offsetFactors_[i];
        return offset;
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(const std::vector<size_t>& ind) const
    {
        return this->calculate_offset(&ind[0]);
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y) const
    {
        GADGET_DEBUG_CHECK_THROW(D==2);
        return x + y * offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y, size_t z) const
    {
        GADGET_DEBUG_CHECK_THROW(D==3);
        return x + (y * offsetFactors_[1]) + (z * offsetFactors_[2]);
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y, size_t z, size_t s) const
    {
        GADGET_DEBUG_CHECK_THROW(D==4);
        return x + (y * offsetFactors_[1]) + (z * offsetFactors_[2]) + (s * offsetFactors_[3]);
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p) const
    {
        GADGET_DEBUG_CHECK_THROW(D==5);
        return x + (y * offsetFactors_[1]) + (z * offsetFactors_[2]) + (s * offsetFactors_[3]) + (p * offsetFactors_[4]);
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const
    {
        GADGET_DEBUG_CHECK_THROW(D==6);
        return x + (y * offsetFactors_[1]) + (z * offsetFactors_[2]) + (s * offsetFactors_[3]) + (p * offsetFactors_[4]) + (r * offsetFactors_[5]);
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const
    {
        GADGET_DEBUG_CHECK_THROW(D==7);
        return x + (y * offsetFactors_[1]) + (z * offsetFactors_[2]) + (s * offsetFactors_[3]) + (p * offsetFactors_[4]) + (r * offsetFactors_[5]) + (a * offsetFactors_[6]);
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const
    {
        GADGET_DEBUG_CHECK_THROW(D==8);
        return x + (y * offsetFactors_[1]) + (z * offsetFactors_[2]) + (s * offsetFactors_[3]) + (p * offsetFactors_[4]) + (r * offsetFactors_[5]) + (a * offsetFactors_[6]) + (q * offsetFactors_[7]);
    }

    template <typename T, unsigned int D> 
    inline size_t hoNDImage<T, D>::calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const
    {
        GADGET_DEBUG_CHECK_THROW(D==9);
        return x + (y * offsetFactors_[1]) + (z * offsetFactors_[2]) + (s * offsetFactors_[3]) + (p * offsetFactors_[4]) + (r * offsetFactors_[5]) + (a * offsetFactors_[6]) + (q * offsetFactors_[7]) + (u * offsetFactors_[8]);
    }

    template <typename T, unsigned int D> 
    inline std::vector<size_t> hoNDImage<T, D>::calculate_index( size_t offset ) const
    {
        std::vector<size_t> index(D, 0);
        this->calculate_index(offset, index);
        return index;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t* index ) const
    {
        GADGET_DEBUG_CHECK_THROW(index!=NULL);

        unsigned int i;
        for( i=D-1; i>0; i-- )
        {
            index[i] = offset / offsetFactors_[i];
            offset %= offsetFactors_[i];
        }
        index[0] = offset;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, std::vector<size_t>& index ) const
    {
        index.resize(D, 0);
        this->calculate_index(offset, &index[0]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, coord_type* index ) const
    {
        unsigned int i;
        for( i=D-1; i>0; i-- )
        {
            index[i] =(coord_type)( offset / offsetFactors_[i] );
            offset %= offsetFactors_[i];
        }
        index[0] = (coord_type)offset;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==2);
        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y, size_t& z ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==3);

        z = offset / offsetFactors_[2];
        offset %= offsetFactors_[2];

        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==4);

        s = offset / offsetFactors_[3];
        offset %= offsetFactors_[3];

        z = offset / offsetFactors_[2];
        offset %= offsetFactors_[2];

        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==5);

        p = offset / offsetFactors_[4];
        offset %= offsetFactors_[4];

        s = offset / offsetFactors_[3];
        offset %= offsetFactors_[3];

        z = offset / offsetFactors_[2];
        offset %= offsetFactors_[2];

        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==6);

        r = offset / offsetFactors_[5];
        offset %= offsetFactors_[5];

        p = offset / offsetFactors_[4];
        offset %= offsetFactors_[4];

        s = offset / offsetFactors_[3];
        offset %= offsetFactors_[3];

        z = offset / offsetFactors_[2];
        offset %= offsetFactors_[2];

        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r, size_t& a ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==7);

        a = offset / offsetFactors_[6];
        offset %= offsetFactors_[6];

        r = offset / offsetFactors_[5];
        offset %= offsetFactors_[5];

        p = offset / offsetFactors_[4];
        offset %= offsetFactors_[4];

        s = offset / offsetFactors_[3];
        offset %= offsetFactors_[3];

        z = offset / offsetFactors_[2];
        offset %= offsetFactors_[2];

        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r, size_t& a, size_t& q ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==8);

        q = offset / offsetFactors_[7];
        offset %= offsetFactors_[7];

        a = offset / offsetFactors_[6];
        offset %= offsetFactors_[6];

        r = offset / offsetFactors_[5];
        offset %= offsetFactors_[5];

        p = offset / offsetFactors_[4];
        offset %= offsetFactors_[4];

        s = offset / offsetFactors_[3];
        offset %= offsetFactors_[3];

        z = offset / offsetFactors_[2];
        offset %= offsetFactors_[2];

        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r, size_t& a, size_t& q, size_t& u ) const
    {
        GADGET_DEBUG_CHECK_THROW(D==9);

        u = offset / offsetFactors_[8];
        offset %= offsetFactors_[8];

        q = offset / offsetFactors_[7];
        offset %= offsetFactors_[7];

        a = offset / offsetFactors_[6];
        offset %= offsetFactors_[6];

        r = offset / offsetFactors_[5];
        offset %= offsetFactors_[5];

        p = offset / offsetFactors_[4];
        offset %= offsetFactors_[4];

        s = offset / offsetFactors_[3];
        offset %= offsetFactors_[3];

        z = offset / offsetFactors_[2];
        offset %= offsetFactors_[2];

        y = offset / offsetFactors_[1];
        x = offset % offsetFactors_[1];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( const size_t* ind )
    {
        size_t idx = this->calculate_offset(ind);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( const size_t* ind ) const
    {
        size_t idx = this->calculate_offset(ind);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( const std::vector<size_t>& ind )
    {
        size_t idx = this->calculate_offset(ind);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( const std::vector<size_t>& ind ) const
    {
        size_t idx = this->calculate_offset(ind);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator[]( size_t x )
    {
        GADGET_DEBUG_CHECK_THROW(x < this->elements_);
        return this->data_[x];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator[]( size_t x ) const
    {
        GADGET_DEBUG_CHECK_THROW(x < this->elements_);
        return this->data_[x];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x )
    {
        GADGET_DEBUG_CHECK_THROW(x < this->elements_);
        return this->data_[x];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x ) const
    {
        GADGET_DEBUG_CHECK_THROW(x < this->elements_);
        return this->data_[x];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y )
    {
        size_t idx = this->calculate_offset(x, y);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y ) const
    {
        size_t idx = this->calculate_offset(x, y);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z )
    {
        size_t idx = this->calculate_offset(x, y, z);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z ) const
    {
        size_t idx = this->calculate_offset(x, y, z);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s )
    {
        size_t idx = this->calculate_offset(x, y, z, s);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s ) const
    {
        size_t idx = this->calculate_offset(x, y, z, s);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p )
    {
        size_t idx = this->calculate_offset(x, y, z, s, p);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p ) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r )
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r ) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a )
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a ) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q )
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q ) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u )
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q, u);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    inline const T& hoNDImage<T, D>::operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u ) const
    {
        size_t idx = this->calculate_offset(x, y, z, s, p, r, a, q, u);
        GADGET_DEBUG_CHECK_THROW(idx < this->elements_);
        return this->data_[idx];
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::fill(T value)
    {
        std::fill(this->get_data_ptr(), this->get_data_ptr()+this->get_number_of_elements(), value);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(const coord_type* ind, coord_type* coord) const
    {
        unsigned int ii, jj;
        for(ii=0; ii<D; ii++)
        {
            coord[ii] = 0;

            for(jj=0; jj<D; jj++)
            {
                coord[ii] += this->axis_[jj][ii] * ( ind[jj] * this->pixelSize_[jj] );
            }

            coord[ii] += this->origin_[ii];
        }
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(const std::vector<coord_type>& ind, std::vector<coord_type>& coord) const
    {
        GADGET_DEBUG_CHECK_THROW(ind.size() >= D);

        if ( coord.size() < D ) coord.resize(D);

        this->image_to_world(&ind[0], &coord[0]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type& cx) const
    {
        GADGET_DEBUG_CHECK_THROW(D==1);
        cx = this->axis_[0][0] * ( x * this->pixelSize_[0] ) + this->origin_[0];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type& cx, coord_type& cy) const
    {
        GADGET_DEBUG_CHECK_THROW(D==2);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];

        cx =    this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy 
              + this->origin_[0];

        cy =    this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy 
              + this->origin_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type z, coord_type& cx, coord_type& cy, coord_type& cz) const
    {
        GADGET_DEBUG_CHECK_THROW(D==3);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];
        coord_type sz = z*this->pixelSize_[2];

        cx =    (this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy) 
              + (this->axis_[2][0] * sz 
              + this->origin_[0]);

        cy =    (this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy) 
              + (this->axis_[2][1] * sz 
              + this->origin_[1]);

        cz =    (this->axis_[0][2] * sx 
              + this->axis_[1][2] * sy) 
              + (this->axis_[2][2] * sz 
              + this->origin_[2]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs) const
    {
        GADGET_DEBUG_CHECK_THROW(D==4);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];
        coord_type sz = z*this->pixelSize_[2];
        coord_type ss = s*this->pixelSize_[3];

        cx =    (this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy) 
              + (this->axis_[2][0] * sz 
              + this->axis_[3][0] * ss) 
              + this->origin_[0];

        cy =    (this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy) 
              + (this->axis_[2][1] * sz 
              + this->axis_[3][1] * ss) 
              + this->origin_[1];

        cz =    (this->axis_[0][2] * sx 
              + this->axis_[1][2] * sy) 
              + (this->axis_[2][2] * sz 
              + this->axis_[3][2] * ss) 
              + this->origin_[2];

        cs =    (this->axis_[0][3] * sx 
              + this->axis_[1][3] * sy) 
              + (this->axis_[2][3] * sz 
              + this->axis_[3][3] * ss) 
              + this->origin_[3];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp) const
    {
        GADGET_DEBUG_CHECK_THROW(D==5);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];
        coord_type sz = z*this->pixelSize_[2];
        coord_type ss = s*this->pixelSize_[3];
        coord_type sp = p*this->pixelSize_[4];

        cx =    (this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy) 
              + (this->axis_[2][0] * sz 
              + this->axis_[3][0] * ss) 
              + (this->axis_[4][0] * sp 
              + this->origin_[0]);

        cy =    (this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy) 
              + (this->axis_[2][1] * sz 
              + this->axis_[3][1] * ss) 
              + (this->axis_[4][1] * sp 
              + this->origin_[1]);

        cz =    (this->axis_[0][2] * sx 
              + this->axis_[1][2] * sy) 
              + (this->axis_[2][2] * sz 
              + this->axis_[3][2] * ss) 
              + (this->axis_[4][2] * sp 
              + this->origin_[2]);

        cs =    (this->axis_[0][3] * sx 
              + this->axis_[1][3] * sy) 
              + (this->axis_[2][3] * sz 
              + this->axis_[3][3] * ss) 
              + (this->axis_[4][3] * sp 
              + this->origin_[3]);

        cp =    (this->axis_[0][4] * sx 
              + this->axis_[1][4] * sy) 
              + (this->axis_[2][4] * sz 
              + this->axis_[3][4] * ss) 
              + (this->axis_[4][4] * sp 
              + this->origin_[4]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr) const
    {
        GADGET_DEBUG_CHECK_THROW(D==6);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];
        coord_type sz = z*this->pixelSize_[2];
        coord_type ss = s*this->pixelSize_[3];
        coord_type sp = p*this->pixelSize_[4];
        coord_type sr = r*this->pixelSize_[5];

        cx =    (this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy) 
              + (this->axis_[2][0] * sz 
              + this->axis_[3][0] * ss) 
              + (this->axis_[4][0] * sp 
              + this->axis_[5][0] * sr) 
              + this->origin_[0];

        cy =    (this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy) 
              + (this->axis_[2][1] * sz 
              + this->axis_[3][1] * ss) 
              + (this->axis_[4][1] * sp 
              + this->axis_[5][1] * sr) 
              + this->origin_[1];

        cz =    (this->axis_[0][2] * sx 
              + this->axis_[1][2] * sy) 
              + (this->axis_[2][2] * sz 
              + this->axis_[3][2] * ss) 
              + (this->axis_[4][2] * sp 
              + this->axis_[5][2] * sr) 
              + this->origin_[2];

        cs =    (this->axis_[0][3] * sx 
              + this->axis_[1][3] * sy) 
              + (this->axis_[2][3] * sz 
              + this->axis_[3][3] * ss) 
              + (this->axis_[4][3] * sp 
              + this->axis_[5][3] * sr) 
              + this->origin_[3];

        cp =    (this->axis_[0][4] * sx 
              + this->axis_[1][4] * sy) 
              + (this->axis_[2][4] * sz 
              + this->axis_[3][4] * ss) 
              + (this->axis_[4][4] * sp 
              + this->axis_[5][4] * sr) 
              + this->origin_[4];

        cr =    (this->axis_[0][5] * sx 
              + this->axis_[1][5] * sy) 
              + (this->axis_[2][5] * sz 
              + this->axis_[3][5] * ss) 
              + (this->axis_[4][5] * sp 
              + this->axis_[5][5] * sr) 
              + this->origin_[5];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca) const
    {
        GADGET_DEBUG_CHECK_THROW(D==7);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];
        coord_type sz = z*this->pixelSize_[2];
        coord_type ss = s*this->pixelSize_[3];
        coord_type sp = p*this->pixelSize_[4];
        coord_type sr = r*this->pixelSize_[5];
        coord_type sa = a*this->pixelSize_[6];

        cx =    (this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy) 
              + (this->axis_[2][0] * sz 
              + this->axis_[3][0] * ss) 
              + (this->axis_[4][0] * sp 
              + this->axis_[5][0] * sr) 
              + (this->axis_[6][0] * sa 
              + this->origin_[0]);

        cy =    (this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy) 
              + (this->axis_[2][1] * sz 
              + this->axis_[3][1] * ss) 
              + (this->axis_[4][1] * sp 
              + this->axis_[5][1] * sr) 
              + (this->axis_[6][1] * sa 
              + this->origin_[1]);

        cz =    (this->axis_[0][2] * sx 
              + this->axis_[1][2] * sy) 
              + (this->axis_[2][2] * sz 
              + this->axis_[3][2] * ss) 
              + (this->axis_[4][2] * sp 
              + this->axis_[5][2] * sr) 
              + (this->axis_[6][2] * sa 
              + this->origin_[2]);

        cs =    (this->axis_[0][3] * sx 
              + this->axis_[1][3] * sy) 
              + (this->axis_[2][3] * sz 
              + this->axis_[3][3] * ss) 
              + (this->axis_[4][3] * sp 
              + this->axis_[5][3] * sr) 
              + (this->axis_[6][3] * sa 
              + this->origin_[3]);

        cp =    (this->axis_[0][4] * sx 
              + this->axis_[1][4] * sy) 
              + (this->axis_[2][4] * sz 
              + this->axis_[3][4] * ss) 
              + (this->axis_[4][4] * sp 
              + this->axis_[5][4] * sr) 
              + (this->axis_[6][4] * sa 
              + this->origin_[4]);

        cr =    (this->axis_[0][5] * sx 
              + this->axis_[1][5] * sy) 
              + (this->axis_[2][5] * sz 
              + this->axis_[3][5] * ss) 
              + (this->axis_[4][5] * sp 
              + this->axis_[5][5] * sr) 
              + (this->axis_[6][5] * sa 
              + this->origin_[5]);

        ca =    (this->axis_[0][6] * sx 
              + this->axis_[1][6] * sy) 
              + (this->axis_[2][6] * sz 
              + this->axis_[3][6] * ss) 
              + (this->axis_[4][6] * sp 
              + this->axis_[5][6] * sr) 
              + (this->axis_[6][6] * sa 
              + this->origin_[6]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq) const
    {
        GADGET_DEBUG_CHECK_THROW(D==8);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];
        coord_type sz = z*this->pixelSize_[2];
        coord_type ss = s*this->pixelSize_[3];
        coord_type sp = p*this->pixelSize_[4];
        coord_type sr = r*this->pixelSize_[5];
        coord_type sa = a*this->pixelSize_[6];
        coord_type sq = q*this->pixelSize_[7];

        cx =    (this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy) 
              + (this->axis_[2][0] * sz 
              + this->axis_[3][0] * ss) 
              + (this->axis_[4][0] * sp 
              + this->axis_[5][0] * sr) 
              + (this->axis_[6][0] * sa 
              + this->axis_[7][0] * sq) 
              + this->origin_[0];

        cy =    (this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy) 
              + (this->axis_[2][1] * sz 
              + this->axis_[3][1] * ss) 
              + (this->axis_[4][1] * sp 
              + this->axis_[5][1] * sr) 
              + (this->axis_[6][1] * sa 
              + this->axis_[7][1] * sq) 
              + this->origin_[1];

        cz =    (this->axis_[0][2] * sx 
              + this->axis_[1][2] * sy) 
              + (this->axis_[2][2] * sz 
              + this->axis_[3][2] * ss) 
              + (this->axis_[4][2] * sp 
              + this->axis_[5][2] * sr) 
              + (this->axis_[6][2] * sa 
              + this->axis_[7][2] * sq) 
              + this->origin_[2];

        cs =    (this->axis_[0][3] * sx 
              + this->axis_[1][3] * sy) 
              + (this->axis_[2][3] * sz 
              + this->axis_[3][3] * ss) 
              + (this->axis_[4][3] * sp 
              + this->axis_[5][3] * sr) 
              + (this->axis_[6][3] * sa 
              + this->axis_[7][3] * sq) 
              + this->origin_[3];

        cp =    (this->axis_[0][4] * sx 
              + this->axis_[1][4] * sy) 
              + (this->axis_[2][4] * sz 
              + this->axis_[3][4] * ss) 
              + (this->axis_[4][4] * sp 
              + this->axis_[5][4] * sr) 
              + (this->axis_[6][4] * sa 
              + this->axis_[7][4] * sq) 
              + this->origin_[4];

        cr =    (this->axis_[0][5] * sx 
              + this->axis_[1][5] * sy) 
              + (this->axis_[2][5] * sz 
              + this->axis_[3][5] * ss) 
              + (this->axis_[4][5] * sp 
              + this->axis_[5][5] * sr) 
              + (this->axis_[6][5] * sa 
              + this->axis_[7][5] * sq) 
              + this->origin_[5];

        ca =    (this->axis_[0][6] * sx 
              + this->axis_[1][6] * sy) 
              + (this->axis_[2][6] * sz 
              + this->axis_[3][6] * ss) 
              + (this->axis_[4][6] * sp 
              + this->axis_[5][6] * sr) 
              + (this->axis_[6][6] * sa 
              + this->axis_[7][6] * sq) 
              + this->origin_[6];

        cq =    (this->axis_[0][7] * sx 
              + this->axis_[1][7] * sy) 
              + (this->axis_[2][7] * sz 
              + this->axis_[3][7] * ss) 
              + (this->axis_[4][7] * sp 
              + this->axis_[5][7] * sr) 
              + (this->axis_[6][7] * sa 
              + this->axis_[7][7] * sq) 
              + this->origin_[7];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u, coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq, coord_type& cu) const
    {
        GADGET_DEBUG_CHECK_THROW(D==9);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];
        coord_type sz = z*this->pixelSize_[2];
        coord_type ss = s*this->pixelSize_[3];
        coord_type sp = p*this->pixelSize_[4];
        coord_type sr = r*this->pixelSize_[5];
        coord_type sa = a*this->pixelSize_[6];
        coord_type sq = q*this->pixelSize_[7];
        coord_type su = u*this->pixelSize_[8];

        cx =    (this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy) 
              + (this->axis_[2][0] * sz 
              + this->axis_[3][0] * ss) 
              + (this->axis_[4][0] * sp 
              + this->axis_[5][0] * sr) 
              + (this->axis_[6][0] * sa 
              + this->axis_[7][0] * sq) 
              + (this->axis_[8][0] * su 
              + this->origin_[0]);

        cy =    (this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy) 
              + (this->axis_[2][1] * sz 
              + this->axis_[3][1] * ss) 
              + (this->axis_[4][1] * sp 
              + this->axis_[5][1] * sr) 
              + (this->axis_[6][1] * sa 
              + this->axis_[7][1] * sq) 
              + (this->axis_[8][1] * su 
              + this->origin_[1]);

        cz =    (this->axis_[0][2] * sx 
              + this->axis_[1][2] * sy) 
              + (this->axis_[2][2] * sz 
              + this->axis_[3][2] * ss) 
              + (this->axis_[4][2] * sp 
              + this->axis_[5][2] * sr) 
              + (this->axis_[6][2] * sa 
              + this->axis_[7][2] * sq) 
              + (this->axis_[8][2] * su 
              + this->origin_[2]);

        cs =    (this->axis_[0][3] * sx 
              + this->axis_[1][3] * sy) 
              + (this->axis_[2][3] * sz 
              + this->axis_[3][3] * ss) 
              + (this->axis_[4][3] * sp 
              + this->axis_[5][3] * sr) 
              + (this->axis_[6][3] * sa 
              + this->axis_[7][3] * sq) 
              + (this->axis_[8][3] * su 
              + this->origin_[3]);

        cp =    (this->axis_[0][4] * sx 
              + this->axis_[1][4] * sy) 
              + (this->axis_[2][4] * sz 
              + this->axis_[3][4] * ss) 
              + (this->axis_[4][4] * sp 
              + this->axis_[5][4] * sr) 
              + (this->axis_[6][4] * sa 
              + this->axis_[7][4] * sq) 
              + (this->axis_[8][4] * su 
              + this->origin_[4]);

        cr =    (this->axis_[0][5] * sx 
              + this->axis_[1][5] * sy) 
              + (this->axis_[2][5] * sz 
              + this->axis_[3][5] * ss) 
              + (this->axis_[4][5] * sp 
              + this->axis_[5][5] * sr) 
              + (this->axis_[6][5] * sa 
              + this->axis_[7][5] * sq) 
              + (this->axis_[8][5] * su 
              + this->origin_[5]);

        ca =    (this->axis_[0][6] * sx 
              + this->axis_[1][6] * sy) 
              + (this->axis_[2][6] * sz 
              + this->axis_[3][6] * ss) 
              + (this->axis_[4][6] * sp 
              + this->axis_[5][6] * sr) 
              + (this->axis_[6][6] * sa 
              + this->axis_[7][6] * sq) 
              + (this->axis_[8][6] * su 
              + this->origin_[6]);

        cq =    (this->axis_[0][7] * sx 
              + this->axis_[1][7] * sy) 
              + (this->axis_[2][7] * sz 
              + this->axis_[3][7] * ss) 
              + (this->axis_[4][7] * sp 
              + this->axis_[5][7] * sr) 
              + (this->axis_[6][7] * sa 
              + this->axis_[7][7] * sq) 
              + (this->axis_[8][7] * su 
              + this->origin_[7]);

        cu =    (this->axis_[0][8] * sx 
              + this->axis_[1][8] * sy) 
              + (this->axis_[2][8] * sz 
              + this->axis_[3][8] * ss) 
              + (this->axis_[4][8] * sp 
              + this->axis_[5][8] * sr) 
              + (this->axis_[6][8] * sa 
              + this->axis_[7][8] * sq) 
              + (this->axis_[8][8] * su 
              + this->origin_[8]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(const size_t* ind, coord_type* coord) const
    {
        unsigned int ii, jj;
        for(ii=0; ii<D; ii++)
        {
            coord[ii] = 0;

            for(jj=0; jj<D; jj++)
            {
                coord[ii] += this->axis_[jj][ii] * ( ind[jj] * this->pixelSize_[jj] );
            }

            coord[ii] += this->origin_[ii];
        }
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(const std::vector<size_t>& ind, std::vector<coord_type>& coord) const
    {
        GADGET_DEBUG_CHECK_THROW(ind.size() >= D);

        if ( coord.size() < D ) coord.resize(D);

        this->image_to_world(&ind[0], &coord[0]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, coord_type& cx) const
    {
        GADGET_DEBUG_CHECK_THROW(D==1);
        cx = this->axis_[0][0] * ( x * this->pixelSize_[0] ) + this->origin_[0];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, coord_type& cx, coord_type& cy) const
    {
        GADGET_DEBUG_CHECK_THROW(D==2);

        coord_type sx = x*this->pixelSize_[0];
        coord_type sy = y*this->pixelSize_[1];

        cx =    this->axis_[0][0] * sx 
              + this->axis_[1][0] * sy 
              + this->origin_[0];

        cy =    this->axis_[0][1] * sx 
              + this->axis_[1][1] * sy 
              + this->origin_[1];
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, size_t z, coord_type& cx, coord_type& cy, coord_type& cz) const
    {
        this->image_to_world(static_cast<coord_type>(x), static_cast<coord_type>(y), static_cast<coord_type>(z), cx, cy, cz);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, size_t z, size_t s, coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs) const
    {
        this->image_to_world(static_cast<coord_type>(x), static_cast<coord_type>(y), static_cast<coord_type>(z), static_cast<coord_type>(s), cx, cy, cz, cs);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp) const
    {
        this->image_to_world(static_cast<coord_type>(x), static_cast<coord_type>(y), static_cast<coord_type>(z), static_cast<coord_type>(s), static_cast<coord_type>(p), cx, cy, cz, cs, cp);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r,
                        coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr) const
    {
        this->image_to_world(static_cast<coord_type>(x), static_cast<coord_type>(y), static_cast<coord_type>(z), static_cast<coord_type>(s), static_cast<coord_type>(p), static_cast<coord_type>(r), cx, cy, cz, cs, cp, cr);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a,
                        coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca) const
    {
        this->image_to_world(static_cast<coord_type>(x), static_cast<coord_type>(y), static_cast<coord_type>(z), static_cast<coord_type>(s), static_cast<coord_type>(p), static_cast<coord_type>(r), static_cast<coord_type>(a), cx, cy, cz, cs, cp, cr, ca);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q,
                        coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq) const
    {
        this->image_to_world(static_cast<coord_type>(x), static_cast<coord_type>(y), static_cast<coord_type>(z), static_cast<coord_type>(s), static_cast<coord_type>(p), static_cast<coord_type>(r), static_cast<coord_type>(a), static_cast<coord_type>(q), cx, cy, cz, cs, cp, cr, ca, cq);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u,
                        coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq, coord_type& cu) const
    {
        this->image_to_world(static_cast<coord_type>(x), static_cast<coord_type>(y), static_cast<coord_type>(z), static_cast<coord_type>(s), static_cast<coord_type>(p), static_cast<coord_type>(r), static_cast<coord_type>(a), static_cast<coord_type>(q), static_cast<coord_type>(u), cx, cy, cz, cs, cp, cr, ca, cq, cu);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::image_to_world_matrix(hoMatrix<coord_type>& image2world) const
    {
        // image to world matrix = tranlation * rotation * pixelSize_Scaling
        image2world.createMatrix(D+1, D+1);

        // rotation matrix
        hoMatrix<coord_type> rotation(D+1, D+1);
        rotation.setIdentity();

        unsigned int ii, jj;
        for ( jj=0; jj<D; jj++ )
        {
            for ( ii=0; ii<D; ii++ )
            {
                rotation(ii, jj) = this->axis_[jj][ii];
            }
        }

        // pixel scaling matrix
        hoMatrix<coord_type> scaling(D+1, D+1);
        scaling.setIdentity();
        for ( ii=0; ii<D; ii++ )
        {
            scaling(ii, ii) = this->pixelSize_[ii];
        }

        // translation matrix
        hoMatrix<coord_type> translation(D+1, D+1);
        translation.setIdentity();
        for ( ii=0; ii<D; ii++ )
        {
            translation(ii, D) = this->origin_[ii];
        }
        Gadgetron::GeneralMatrixProduct(image2world, rotation, false, scaling, false);
        Gadgetron::GeneralMatrixProduct(rotation, translation, false, image2world, false);
        image2world = rotation;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_image_to_world_matrix(const hoMatrix<coord_type>& image2world)
    {
        GADGET_DEBUG_CHECK_THROW(D+1==image2world.rows());
        GADGET_DEBUG_CHECK_THROW(D+1==image2world.cols());

        // origin
        hoMatrix<coord_type> pt(D+1, 1);
        pt(D, 0) = 1.0;

        hoMatrix<coord_type> res(D+1, 1);

        Gadgetron::GeneralMatrixProduct(res, image2world, false, pt, false);
        memcpy(this->origin_, res.begin(), sizeof(coord_type)*D);

        // rotation
        unsigned int ii, jj;
        for ( ii=0; ii<D; ii++ )
        {
            memset(pt.get_data_ptr(), 0, sizeof(coord_type)*(D+1));
            pt(D, 0) = 1.0;
            pt(ii, 0) = 1.0;

            Gadgetron::GeneralMatrixProduct(res, image2world, false, pt, false);
            for ( jj=0; jj<D; jj++ )
            {
                this->axis_[ii][jj] = res(jj, 0) - this->origin_[jj];
            }

            this->pixelSize_[ii] = this->axis_[ii].abs();
            this->axis_[ii].normalize();
        }
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(const coord_type* coord, coord_type* ind) const
    {
        unsigned int ii, jj;
        for(ii=0; ii<D; ii++)
        {
            ind[ii] = 0;
            for(jj=0; jj<D; jj++)
            {
                ind[ii] += this->axis_[ii][jj] * ( coord[jj] - this->origin_[jj] );
            }

            ind[ii] *= this->pixelSize_reciprocal_[ii];
        }
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(const std::vector<coord_type>& coord, std::vector<coord_type>& ind) const
    {
        GADGET_DEBUG_CHECK_THROW(coord.size()>=D);

        if ( ind.size() < D ) ind.resize(D);

        this->world_to_image(&coord[0], &ind[0]);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type& x) const
    {
        GADGET_DEBUG_CHECK_THROW(D==1);
        x = this->pixelSize_reciprocal_[0] * this->axis_[0][0] * ( cx - this->origin_[0] );
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type& x, coord_type& y) const
    {
        GADGET_DEBUG_CHECK_THROW(D==2);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];

        x = this->pixelSize_reciprocal_[0] * (this->axis_[0][0]*sx + this->axis_[0][1]*sy);
        y = this->pixelSize_reciprocal_[1] * (this->axis_[1][0]*sx + this->axis_[1][1]*sy);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type& x, coord_type& y, coord_type& z) const
    {
        GADGET_DEBUG_CHECK_THROW(D==3);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];
        coord_type sz = cz - this->origin_[2];

        x = this->pixelSize_reciprocal_[0] * (this->axis_[0][0]*sx + this->axis_[0][1]*sy + this->axis_[0][2]*sz);
        y = this->pixelSize_reciprocal_[1] * (this->axis_[1][0]*sx + this->axis_[1][1]*sy + this->axis_[1][2]*sz);
        z = this->pixelSize_reciprocal_[2] * (this->axis_[2][0]*sx + this->axis_[2][1]*sy + this->axis_[2][2]*sz);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type& x, coord_type& y, coord_type& z, coord_type& s) const
    {
        GADGET_DEBUG_CHECK_THROW(D==4);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];
        coord_type sz = cz - this->origin_[2];
        coord_type ss = cs - this->origin_[3];

        x = this->pixelSize_reciprocal_[0] * ((this->axis_[0][0]*sx + this->axis_[0][1]*sy) + (this->axis_[0][2]*sz + this->axis_[0][3]*ss));
        y = this->pixelSize_reciprocal_[1] * ((this->axis_[1][0]*sx + this->axis_[1][1]*sy) + (this->axis_[1][2]*sz + this->axis_[1][3]*ss));
        z = this->pixelSize_reciprocal_[2] * ((this->axis_[2][0]*sx + this->axis_[2][1]*sy) + (this->axis_[2][2]*sz + this->axis_[2][3]*ss));
        s = this->pixelSize_reciprocal_[3] * ((this->axis_[3][0]*sx + this->axis_[3][1]*sy) + (this->axis_[3][2]*sz + this->axis_[3][3]*ss));
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p) const
    {
        GADGET_DEBUG_CHECK_THROW(D==5);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];
        coord_type sz = cz - this->origin_[2];
        coord_type ss = cs - this->origin_[3];
        coord_type sp = cp - this->origin_[4];

        x = this->pixelSize_reciprocal_[0] * ((this->axis_[0][0]*sx + this->axis_[0][1]*sy) + (this->axis_[0][2]*sz + this->axis_[0][3]*ss) + this->axis_[0][4]*sp);
        y = this->pixelSize_reciprocal_[1] * ((this->axis_[1][0]*sx + this->axis_[1][1]*sy) + (this->axis_[1][2]*sz + this->axis_[1][3]*ss) + this->axis_[1][4]*sp);
        z = this->pixelSize_reciprocal_[2] * ((this->axis_[2][0]*sx + this->axis_[2][1]*sy) + (this->axis_[2][2]*sz + this->axis_[2][3]*ss) + this->axis_[2][4]*sp);
        s = this->pixelSize_reciprocal_[3] * ((this->axis_[3][0]*sx + this->axis_[3][1]*sy) + (this->axis_[3][2]*sz + this->axis_[3][3]*ss) + this->axis_[3][4]*sp);
        p = this->pixelSize_reciprocal_[4] * ((this->axis_[4][0]*sx + this->axis_[4][1]*sy) + (this->axis_[4][2]*sz + this->axis_[4][3]*ss) + this->axis_[4][4]*sp);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr, coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r) const
    {
        GADGET_DEBUG_CHECK_THROW(D==6);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];
        coord_type sz = cz - this->origin_[2];
        coord_type ss = cs - this->origin_[3];
        coord_type sp = cp - this->origin_[4];
        coord_type sr = cr - this->origin_[5];

        x = this->pixelSize_reciprocal_[0] * ((this->axis_[0][0]*sx + this->axis_[0][1]*sy) + (this->axis_[0][2]*sz + this->axis_[0][3]*ss) + (this->axis_[0][4]*sp + this->axis_[0][5]*sr));
        y = this->pixelSize_reciprocal_[1] * ((this->axis_[1][0]*sx + this->axis_[1][1]*sy) + (this->axis_[1][2]*sz + this->axis_[1][3]*ss) + (this->axis_[1][4]*sp + this->axis_[1][5]*sr));
        z = this->pixelSize_reciprocal_[2] * ((this->axis_[2][0]*sx + this->axis_[2][1]*sy) + (this->axis_[2][2]*sz + this->axis_[2][3]*ss) + (this->axis_[2][4]*sp + this->axis_[2][5]*sr));
        s = this->pixelSize_reciprocal_[3] * ((this->axis_[3][0]*sx + this->axis_[3][1]*sy) + (this->axis_[3][2]*sz + this->axis_[3][3]*ss) + (this->axis_[3][4]*sp + this->axis_[3][5]*sr));
        p = this->pixelSize_reciprocal_[4] * ((this->axis_[4][0]*sx + this->axis_[4][1]*sy) + (this->axis_[4][2]*sz + this->axis_[4][3]*ss) + (this->axis_[4][4]*sp + this->axis_[4][5]*sr));
        r = this->pixelSize_reciprocal_[5] * ((this->axis_[5][0]*sx + this->axis_[5][1]*sy) + (this->axis_[5][2]*sz + this->axis_[5][3]*ss) + (this->axis_[5][4]*sp + this->axis_[5][5]*sr));
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr, coord_type ca, coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r, coord_type& a) const
    {
        GADGET_DEBUG_CHECK_THROW(D==7);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];
        coord_type sz = cz - this->origin_[2];
        coord_type ss = cs - this->origin_[3];
        coord_type sp = cp - this->origin_[4];
        coord_type sr = cr - this->origin_[5];
        coord_type sa = ca - this->origin_[6];

        x = this->pixelSize_reciprocal_[0] * ((this->axis_[0][0]*sx + this->axis_[0][1]*sy) + (this->axis_[0][2]*sz + this->axis_[0][3]*ss) + (this->axis_[0][4]*sp + this->axis_[0][5]*sr) + this->axis_[0][6]*sa);
        y = this->pixelSize_reciprocal_[1] * ((this->axis_[1][0]*sx + this->axis_[1][1]*sy) + (this->axis_[1][2]*sz + this->axis_[1][3]*ss) + (this->axis_[1][4]*sp + this->axis_[1][5]*sr) + this->axis_[1][6]*sa);
        z = this->pixelSize_reciprocal_[2] * ((this->axis_[2][0]*sx + this->axis_[2][1]*sy) + (this->axis_[2][2]*sz + this->axis_[2][3]*ss) + (this->axis_[2][4]*sp + this->axis_[2][5]*sr) + this->axis_[2][6]*sa);
        s = this->pixelSize_reciprocal_[3] * ((this->axis_[3][0]*sx + this->axis_[3][1]*sy) + (this->axis_[3][2]*sz + this->axis_[3][3]*ss) + (this->axis_[3][4]*sp + this->axis_[3][5]*sr) + this->axis_[3][6]*sa);
        p = this->pixelSize_reciprocal_[4] * ((this->axis_[4][0]*sx + this->axis_[4][1]*sy) + (this->axis_[4][2]*sz + this->axis_[4][3]*ss) + (this->axis_[4][4]*sp + this->axis_[4][5]*sr) + this->axis_[4][6]*sa);
        r = this->pixelSize_reciprocal_[5] * ((this->axis_[5][0]*sx + this->axis_[5][1]*sy) + (this->axis_[5][2]*sz + this->axis_[5][3]*ss) + (this->axis_[5][4]*sp + this->axis_[5][5]*sr) + this->axis_[5][6]*sa);
        a = this->pixelSize_reciprocal_[6] * ((this->axis_[6][0]*sx + this->axis_[6][1]*sy) + (this->axis_[6][2]*sz + this->axis_[6][3]*ss) + (this->axis_[6][4]*sp + this->axis_[6][5]*sr) + this->axis_[6][6]*sa);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr, coord_type ca, coord_type cq, coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r, coord_type& a, coord_type& q) const
    {
        GADGET_DEBUG_CHECK_THROW(D==8);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];
        coord_type sz = cz - this->origin_[2];
        coord_type ss = cs - this->origin_[3];
        coord_type sp = cp - this->origin_[4];
        coord_type sr = cr - this->origin_[5];
        coord_type sa = ca - this->origin_[6];
        coord_type sq = cq - this->origin_[7];

        x = this->pixelSize_reciprocal_[0] * ((this->axis_[0][0]*sx + this->axis_[0][1]*sy) + (this->axis_[0][2]*sz + this->axis_[0][3]*ss) + (this->axis_[0][4]*sp + this->axis_[0][5]*sr) + (this->axis_[0][6]*sa + this->axis_[0][7]*sq));
        y = this->pixelSize_reciprocal_[1] * ((this->axis_[1][0]*sx + this->axis_[1][1]*sy) + (this->axis_[1][2]*sz + this->axis_[1][3]*ss) + (this->axis_[1][4]*sp + this->axis_[1][5]*sr) + (this->axis_[1][6]*sa + this->axis_[1][7]*sq));
        z = this->pixelSize_reciprocal_[2] * ((this->axis_[2][0]*sx + this->axis_[2][1]*sy) + (this->axis_[2][2]*sz + this->axis_[2][3]*ss) + (this->axis_[2][4]*sp + this->axis_[2][5]*sr) + (this->axis_[2][6]*sa + this->axis_[2][7]*sq));
        s = this->pixelSize_reciprocal_[3] * ((this->axis_[3][0]*sx + this->axis_[3][1]*sy) + (this->axis_[3][2]*sz + this->axis_[3][3]*ss) + (this->axis_[3][4]*sp + this->axis_[3][5]*sr) + (this->axis_[3][6]*sa + this->axis_[3][7]*sq));
        p = this->pixelSize_reciprocal_[4] * ((this->axis_[4][0]*sx + this->axis_[4][1]*sy) + (this->axis_[4][2]*sz + this->axis_[4][3]*ss) + (this->axis_[4][4]*sp + this->axis_[4][5]*sr) + (this->axis_[4][6]*sa + this->axis_[4][7]*sq));
        r = this->pixelSize_reciprocal_[5] * ((this->axis_[5][0]*sx + this->axis_[5][1]*sy) + (this->axis_[5][2]*sz + this->axis_[5][3]*ss) + (this->axis_[5][4]*sp + this->axis_[5][5]*sr) + (this->axis_[5][6]*sa + this->axis_[5][7]*sq));
        a = this->pixelSize_reciprocal_[6] * ((this->axis_[6][0]*sx + this->axis_[6][1]*sy) + (this->axis_[6][2]*sz + this->axis_[6][3]*ss) + (this->axis_[6][4]*sp + this->axis_[6][5]*sr) + (this->axis_[6][6]*sa + this->axis_[6][7]*sq));
        q = this->pixelSize_reciprocal_[7] * ((this->axis_[7][0]*sx + this->axis_[7][1]*sy) + (this->axis_[7][2]*sz + this->axis_[7][3]*ss) + (this->axis_[7][4]*sp + this->axis_[7][5]*sr) + (this->axis_[7][6]*sa + this->axis_[7][7]*sq));
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr, coord_type ca, coord_type cq, coord_type cu, coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r, coord_type& a, coord_type& q, coord_type& u) const
    {
        GADGET_DEBUG_CHECK_THROW(D==9);

        coord_type sx = cx - this->origin_[0];
        coord_type sy = cy - this->origin_[1];
        coord_type sz = cz - this->origin_[2];
        coord_type ss = cs - this->origin_[3];
        coord_type sp = cp - this->origin_[4];
        coord_type sr = cr - this->origin_[5];
        coord_type sa = ca - this->origin_[6];
        coord_type sq = cq - this->origin_[7];
        coord_type su = cu - this->origin_[8];

        x = this->pixelSize_reciprocal_[0] * ((this->axis_[0][0]*sx + this->axis_[0][1]*sy) + (this->axis_[0][2]*sz + this->axis_[0][3]*ss) + (this->axis_[0][4]*sp + this->axis_[0][5]*sr) + (this->axis_[0][6]*sa + this->axis_[0][7]*sq) + this->axis_[0][8]*su);
        y = this->pixelSize_reciprocal_[1] * ((this->axis_[1][0]*sx + this->axis_[1][1]*sy) + (this->axis_[1][2]*sz + this->axis_[1][3]*ss) + (this->axis_[1][4]*sp + this->axis_[1][5]*sr) + (this->axis_[1][6]*sa + this->axis_[1][7]*sq) + this->axis_[1][8]*su);
        z = this->pixelSize_reciprocal_[2] * ((this->axis_[2][0]*sx + this->axis_[2][1]*sy) + (this->axis_[2][2]*sz + this->axis_[2][3]*ss) + (this->axis_[2][4]*sp + this->axis_[2][5]*sr) + (this->axis_[2][6]*sa + this->axis_[2][7]*sq) + this->axis_[2][8]*su);
        s = this->pixelSize_reciprocal_[3] * ((this->axis_[3][0]*sx + this->axis_[3][1]*sy) + (this->axis_[3][2]*sz + this->axis_[3][3]*ss) + (this->axis_[3][4]*sp + this->axis_[3][5]*sr) + (this->axis_[3][6]*sa + this->axis_[3][7]*sq) + this->axis_[3][8]*su);
        p = this->pixelSize_reciprocal_[4] * ((this->axis_[4][0]*sx + this->axis_[4][1]*sy) + (this->axis_[4][2]*sz + this->axis_[4][3]*ss) + (this->axis_[4][4]*sp + this->axis_[4][5]*sr) + (this->axis_[4][6]*sa + this->axis_[4][7]*sq) + this->axis_[4][8]*su);
        r = this->pixelSize_reciprocal_[5] * ((this->axis_[5][0]*sx + this->axis_[5][1]*sy) + (this->axis_[5][2]*sz + this->axis_[5][3]*ss) + (this->axis_[5][4]*sp + this->axis_[5][5]*sr) + (this->axis_[5][6]*sa + this->axis_[5][7]*sq) + this->axis_[5][8]*su);
        a = this->pixelSize_reciprocal_[6] * ((this->axis_[6][0]*sx + this->axis_[6][1]*sy) + (this->axis_[6][2]*sz + this->axis_[6][3]*ss) + (this->axis_[6][4]*sp + this->axis_[6][5]*sr) + (this->axis_[6][6]*sa + this->axis_[6][7]*sq) + this->axis_[6][8]*su);
        q = this->pixelSize_reciprocal_[7] * ((this->axis_[7][0]*sx + this->axis_[7][1]*sy) + (this->axis_[7][2]*sz + this->axis_[7][3]*ss) + (this->axis_[7][4]*sp + this->axis_[7][5]*sr) + (this->axis_[7][6]*sa + this->axis_[7][7]*sq) + this->axis_[7][8]*su);
        u = this->pixelSize_reciprocal_[8] * ((this->axis_[8][0]*sx + this->axis_[8][1]*sy) + (this->axis_[8][2]*sz + this->axis_[8][3]*ss) + (this->axis_[8][4]*sp + this->axis_[8][5]*sr) + (this->axis_[8][6]*sa + this->axis_[8][7]*sq) + this->axis_[8][8]*su);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::world_to_image_matrix(hoMatrix<coord_type>& world2image) const
    {
        // world to image matrix = inv(pixelSize_Scaling) * inv(rotation) * inv(tranlation)
        world2image.createMatrix(D+1, D+1);

        // rotation matrix
        hoMatrix<coord_type> rotation(D+1, D+1);
        rotation.setIdentity();

        unsigned int ii, jj;
        for ( jj=0; jj<D; jj++ )
        {
            for ( ii=0; ii<D; ii++ )
            {
                rotation(jj, ii) = this->axis_[jj][ii];
            }
        }

        // pixel scaling matrix
        hoMatrix<coord_type> scaling(D+1, D+1);
        scaling.setIdentity();
        for ( ii=0; ii<D; ii++ )
        {
            scaling(ii, ii) = this->pixelSize_reciprocal_[ii];
        }

        // translation matrix
        hoMatrix<coord_type> translation(D+1, D+1);
        translation.setIdentity();
        for ( ii=0; ii<D; ii++ )
        {
            translation(ii, D) = -this->origin_[ii];
        }

        Gadgetron::GeneralMatrixProduct(world2image, rotation, false, translation, false);
        Gadgetron::GeneralMatrixProduct(rotation, scaling, false, world2image, false);

        world2image = rotation;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::set_world_to_image_matrix(const hoMatrix<coord_type>& world2image)
    {
        GADGET_THROW("This function is not implemented ... ");
    }

    template <typename T, unsigned int D> 
    inline bool hoNDImage<T, D>::in_image_region(const std::vector<size_t>& start, std::vector<size_t>& size)
    {
        GADGET_DEBUG_CHECK_THROW(start.size()>=D);
        GADGET_DEBUG_CHECK_THROW(size.size()>=D);

        if ( !this->point_in_range(start) ) return false;

        std::vector<size_t> end(D);
        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            end[ii] = start[ii]+size[ii]-1;
        }

        if ( !this->point_in_range(end) ) return false;

        return true;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::get_sub_image(const std::vector<size_t>& start, std::vector<size_t>& size, Self& out)
    {
        GADGET_DEBUG_CHECK_THROW(start.size()>=D);
        GADGET_DEBUG_CHECK_THROW(size.size()>=D);

        if ( !this->in_image_region(start, size) )
        {
        	GWARN_STREAM("Sub-image regin is not in the image ... ");
            return;
        }

        out.create(size);

        memcpy(out.pixelSize_, this->pixelSize_, sizeof(coord_type)*D);
        memcpy(out.pixelSize_reciprocal_, this->pixelSize_reciprocal_, sizeof(coord_type)*D);

        unsigned int ii;
        for ( ii=0; ii<D; ii++ )
        {
            out.axis_[ii] = this->axis_[ii];
        }

        size_t N = out.get_number_of_elements() / size[0];

        long long t;

        #pragma omp parallel default(none) private(t) shared(N, size, out)
        {
            std::vector<size_t> indOut(D), ind(D);

            #pragma omp for
            for ( t=0; t<(long long)N; t++ )
            {
                out.calculate_index( (size_t)(t*size[0]), indOut);

                unsigned int ii;
                for ( ii=0; ii<D; ii++ )
                {
                    ind[ii] = indOut[ii]+start[ii];
                }

                size_t offset = this->calculate_offset(ind);

                memcpy(out.begin()+t*size[0], this->data_+offset, sizeof(T)*size[0]);
            }
        }

        std::vector<coord_type> origin_out(D);
        this->image_to_world(start, origin_out);

        memcpy(out.origin_, &origin_out[0], sizeof(coord_type)*D);
    }

    template <typename T, unsigned int D> 
    bool hoNDImage<T, D>::serialize(char*& buf, size_t& len) const 
    {
        try
        {
            if ( buf != NULL ) delete[] buf;

            // number of dimensions + dimension vector + pixel size + origin + axis + contents
            len = sizeof(unsigned int) + sizeof(size_t)*D 
                + sizeof(coord_type)*D + sizeof(coord_type)*D + sizeof(coord_type)*D*D 
                + sizeof(T)*this->elements_;

            buf = new char[len];
            GADGET_CHECK_RETURN_FALSE(buf!=NULL);

            unsigned int NDim=D;

            size_t offset = 0;
            memcpy(buf, &NDim, sizeof(unsigned int));
            offset += sizeof(unsigned int);

            if ( NDim > 0 )
            {
                memcpy(buf+offset, &(dimensions_[0]), sizeof(size_t)*D);
                offset += sizeof(size_t)*D;

                memcpy(buf+offset, this->pixelSize_, sizeof(coord_type)*D);
                offset += sizeof(coord_type)*D;

                memcpy(buf+offset, this->origin_, sizeof(coord_type)*D);
                offset += sizeof(coord_type)*D;

                unsigned int ii;
                for ( ii=0; ii<D; ii++ )
                {
                    memcpy(buf+offset, this->axis_[ii].begin(), sizeof(coord_type)*D);
                    offset += sizeof(coord_type)*D;
                }

                memcpy(buf+offset, this->data_, sizeof(T)*elements_);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImage<T, D>::serialize(...) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool hoNDImage<T, D>::deserialize(char* buf, size_t& len)
    {
        try
        {
            unsigned int NDim;
            memcpy(&NDim, buf, sizeof(unsigned int));
            if ( NDim != D )
            {
                GERROR_STREAM("hoNDImage<T, D>::deserialize(...) : number of image dimensions does not match ... ");
                return false;
            }

            size_t offset = sizeof(unsigned int);

            unsigned int ii;

            if ( NDim > 0 )
            {
                std::vector<size_t> dimensions(NDim);

                memcpy(&dimensions[0], buf+offset, sizeof(size_t)*D);
                offset += sizeof(size_t)*D;

                this->create(dimensions);

                memcpy(this->pixelSize_, buf+offset, sizeof(coord_type)*D);
                offset += sizeof(coord_type)*D;

                for ( ii=0; ii<D; ii++ )
                {
                    this->pixelSize_reciprocal_[ii] = coord_type(1.0)/this->pixelSize_[ii];
                }

                memcpy(this->origin_, buf+offset, sizeof(coord_type)*D);
                offset += sizeof(coord_type)*D;

                for ( ii=0; ii<D; ii++ )
                {
                    memcpy(this->axis_[ii].begin(), buf+offset, sizeof(coord_type)*D);
                    offset += sizeof(coord_type)*D;
                }

                // copy the content
                memcpy(this->data_, buf+offset, sizeof(T)*elements_);
                offset += sizeof(T)*elements_;
            }
            else
            {
                this->clear();
            }

            len = offset;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImage<T, D>::deserialize(...) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    void hoNDImage<T, D>::print(std::ostream& os) const
    {
        using namespace std;
        os << "-------------- Gagdgetron ND Image -------------" << endl;
        this->printContent(os);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImage<T, D>::printContent(std::ostream& os) const
    {
        using namespace std;

        os.unsetf(std::ios::scientific);
        os.setf(ios::fixed);

        size_t i, j;

        os << "Image dimension is : " << D << endl;

        os << "Image size is : ";
        for (i=0; i<D; i++ ) 
            os << dimensions_[i] << " "; 
        os << endl;

        int elemTypeSize = sizeof(T);
        std::string elemTypeName = std::string(typeid(T).name());

        os << "Image data type is : " << elemTypeName << std::endl;
        os << "Byte number for each element is : " << elemTypeSize << std::endl;
        os << "Number of array size in bytes is : ";
        os << elements_*elemTypeSize << std::endl;

        os << "Pixel size is : ";
        for (i=0; i<D; i++ ) 
            os << this->pixelSize_[i] << " "; 
        os << endl;

        os << "Origin is : ";
        for (i=0; i<D; i++ ) 
            os << this->origin_[i] << " "; 
        os << endl;

        for (i=0; i<D; i++ )
        {
            os << "Axis " << i << " : [ ";
            for (j=0; j<D; j++ )
            {
                os << this->axis_[i][j] << " "; 
            }
            os << "] " << endl;
        }
        os << endl << ends;
    }
}
