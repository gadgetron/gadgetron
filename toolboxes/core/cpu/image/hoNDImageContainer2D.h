/** \file       hoNDImageContainer2D.h
    \brief      a container class to store a matrix of hoNDImages

                This name "container2D" does not mean the 2D images. It means the container is a 2D array in its storage logic.

                The points of images are stored in this container. However, the images can be deleted if delete_data_on_destruct_==true
                The images are stored as 2D arrays. But every row can have differnet number of images (or columns). So it is not exactly an 
                image matrix.

    \author     Hui Xue
*/

#pragma once

#include "hoNDImage.h"
#include "hoNDArray_elemwise.h"

namespace Gadgetron
{

    template <typename ImageType>
    class hoNDImageContainer2D
    {
    public:

        typedef hoNDImageContainer2D<ImageType> Self;

        typedef typename ImageType::value_type value_type;
        typedef typename ImageType::coord_type coord_type;
        typedef typename ImageType::a_axis_type a_axis_type;
        typedef typename ImageType::axis_type axis_type;

        /// constructors
        hoNDImageContainer2D(bool delete_data_on_destruct=true);
        hoNDImageContainer2D(const hoNDImageContainer2D<ImageType>& a);

        Self& operator=(const Self& rhs);

        virtual ~hoNDImageContainer2D();

        /// create a container with images
        bool create(const std::vector<size_t>& col, bool createImage=true);

        /// create a container with images at certain sizes/pixel sizes/axis
        /// the image will not be filled with zeros
        bool create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions);
        bool create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize);
        bool create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin);
        bool create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis);

        /// create a container with every row being its own image size
        bool create(const std::vector<size_t>& col, const std::vector< std::vector<size_t> >& dimensions);

        /// create a container from a chunk of memory
        /// the dim.size() = ImageType.get_number_of_dimensions()+1
        /// e.g., a 3D memory chunk [RO E1 N] is used to allocate N [RO E1] images
        /// the container will have 1 row and N columns
        bool create(value_type* buf, const std::vector<size_t>& dim);

        /// clear the matrix, if delete_data_on_destruct_==true, delete all stored images
        bool clear();

        /// copy from a container 2D, deep copy image content
        bool copyFrom(const Self& a);

        /// fill all images with zeros
        bool fillWithZeros();

        /// whether two containers have the same size
        template <typename ImageType2> 
        bool dimensions_equal_container(const hoNDImageContainer2D<ImageType2>& a) const
        {
            if ( this->rows() != a.rows() ) return false;

            unsigned int row;
            for ( row=0; row<this->rows(); row++ )
            {
                if ( this->cols(row) != a.cols(row) )
                {
                    return false;
                }
            }

            return true;
        }

        /// add one image to a row at end
        bool push_back(ImageType& im, size_t row);

        /// add one image to a row at head
        bool push_front(ImageType& im, size_t row);

        /// add one image to a row
        bool insert(ImageType& im, size_t row, size_t col);

        /// pop an image from a row end
        bool pop_back(ImageType*& im, size_t row);

        /// pop an image from a row head
        bool pop_front(ImageType*& im, size_t row);

        /// remove an image from the storage
        bool remove(ImageType*& im, size_t row, size_t col);

        /// if delete_data_on_destruct_==true, the image will be deleted
        bool remove(size_t row, size_t col);

        /// expand the container by certain number of rows
        bool expand(size_t newRows);

        /// insert one row
        bool insert(std::vector<ImageType*>& im_array, size_t row);

        /// remove one row
        bool remove(std::vector<ImageType*>& im_array, size_t row);
        /// if delete_data_on_destruct_==true, the image will be deleted
        bool remove(size_t row);

        /// get image pointers
        ImageType& get(size_t row, size_t col);
        const ImageType& get(size_t row, size_t col) const;

        ImageType& operator() (size_t row, size_t col);
        const ImageType& operator() (size_t row, size_t col) const;

        /// get one row
        bool get(std::vector<ImageType*>& im_array, size_t row) const;

        /// get number of all images in the container
        size_t get_number_of_all_images();

        /// get all images
        bool get_all_images(std::vector<ImageType*>& im_array);

        /// set image pointer
        bool set(ImageType* pImage, size_t row, size_t col);

        /// convert one row to a hoNDArray
        /// all images in this row should have the same dimensions; if not, return false
        bool to_NDArray(size_t row, hoNDArray<value_type>& a) const;

        /// whether to delete the memory on destruction
        bool delete_data_on_destruct() const;
        void delete_data_on_destruct(bool d);

        /// get number of row and column
        size_t rows() const;
        size_t cols(size_t row) const;
        std::vector<size_t> cols() const;

        /// check whether all images in a row have the same dimensions/pixelSizes/axises
        bool has_identical_dimensions(unsigned int row) const;
        bool has_identical_pixel_size(unsigned int row) const;
        bool has_identical_axis(unsigned int row) const;
        bool has_identical_image_geometry(unsigned int row) const;

        /// serialize/deserialize
        virtual bool serialize(char*& buf, size_t& len) const;
        virtual bool deserialize(char* buf, size_t& len);

        /// print out the image container information
        virtual void print(std::ostream& os) const;

    protected:

        std::vector< std::vector<ImageType*> > image_container_;

        bool delete_data_on_destruct_;
    };

    template <typename ImageType> 
    hoNDImageContainer2D<ImageType>::hoNDImageContainer2D(bool delete_data_on_destruct) : delete_data_on_destruct_(delete_data_on_destruct)
    {
    }

    template <typename ImageType> 
    hoNDImageContainer2D<ImageType>::hoNDImageContainer2D(const hoNDImageContainer2D<ImageType>& a) : delete_data_on_destruct_(false)
    {
        *this = a;
    }

    template <typename ImageType> 
    hoNDImageContainer2D<ImageType>& hoNDImageContainer2D<ImageType>::operator=(const Self& rhs)
    {
        if ( this == &rhs ) return *this;

        this->clear();
        size_t row = rhs.rows();

        size_t ii;
        for ( ii=0; ii<row; ii++ )
        {
            std::vector<ImageType*> a_row;
            rhs.get(a_row, ii);
            this->image_container_.push_back(a_row);
        }

        this->delete_data_on_destruct_ = false;

        return *this;
    }

    template <typename ImageType> 
    hoNDImageContainer2D<ImageType>::~hoNDImageContainer2D()
    {
        this->clear();
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, bool createImage)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->clear());
            if ( createImage )
            {
                this->delete_data_on_destruct(true);
            }
            else
            {
                this->delete_data_on_destruct(false);
            }

            size_t row = col.size();
            image_container_.resize(row);

            unsigned int r, c;
            for ( r=0; r<row; r++ )
            {
                image_container_[r].resize(col[r], NULL);

                if ( createImage )
                {
                    for ( c=0; c<col[r]; c++ )
                    {
                        image_container_[r][c] = new ImageType();
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->clear());
            this->delete_data_on_destruct(true);

            size_t row = col.size();
            image_container_.resize(row);

            unsigned int r, c;
            for ( r=0; r<row; r++ )
            {
                image_container_[r].resize(col[r], NULL);

                for ( c=0; c<col[r]; c++ )
                {
                    image_container_[r][c] = new ImageType(dimensions);
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->clear());
            this->delete_data_on_destruct(true);

            size_t row = col.size();
            image_container_.resize(row);

            unsigned int r, c;
            for ( r=0; r<row; r++ )
            {
                image_container_[r].resize(col[r], NULL);

                for ( c=0; c<col[r]; c++ )
                {
                    image_container_[r][c] = new ImageType(dimensions, pixelSize);
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->clear());
            this->delete_data_on_destruct(true);

            size_t row = col.size();
            image_container_.resize(row);

            unsigned int r, c;
            for ( r=0; r<row; r++ )
            {
                image_container_[r].resize(col[r], NULL);

                for ( c=0; c<col[r]; c++ )
                {
                    image_container_[r][c] = new ImageType(dimensions, pixelSize, origin);
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->clear());
            this->delete_data_on_destruct(true);

            size_t row = col.size();
            image_container_.resize(row);

            unsigned int r, c;
            for ( r=0; r<row; r++ )
            {
                image_container_[r].resize(col[r], NULL);

                for ( c=0; c<col[r]; c++ )
                {
                    image_container_[r][c] = new ImageType(dimensions, pixelSize, origin, axis);
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType>
    inline bool hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector< std::vector<size_t> >& dimensions)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(col.size()==dimensions.size());
            GADGET_CHECK_RETURN_FALSE(this->clear());
            this->delete_data_on_destruct(true);

            size_t row = col.size();
            image_container_.resize(row);

            unsigned int r, c;
            for (r = 0; r<row; r++)
            {
                image_container_[r].resize(col[r], NULL);

                for (c = 0; c<col[r]; c++)
                {
                    image_container_[r][c] = new ImageType(dimensions[r]);
                }
            }
        }
        catch (...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::create(const std::vector<size_t>& col, const std::vector< std::vector<size_t> >& dimensions ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::create(value_type* buf, const std::vector<size_t>& dim)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE( (dim.size()==ImageType::NDIM) || (dim.size()==ImageType::NDIM+1) );

            GADGET_CHECK_RETURN_FALSE(this->clear());
            this->delete_data_on_destruct(true);

            unsigned int ii;
            size_t col;
            std::vector<size_t> dim_im;
            if ( dim.size()==ImageType::NDIM )
            {
                dim_im = dim;
                col = 1;
            }
            else
            {
                dim_im.resize(ImageType::NDIM);
                memcpy(&dim_im[0], &dim[0], sizeof(size_t)*ImageType::NDIM);
                col = dim[ImageType::NDIM];
            }

            size_t row = 1;
            image_container_.resize(row);
            image_container_[0].resize(col);

            size_t numOfPixels = 1;
            for ( ii=0; ii<dim_im.size(); ii++ )
            {
                numOfPixels *= dim_im[ii];
            }

            unsigned int c;
            for ( c=0; c<col; c++ )
            {
                image_container_[0][c] = new ImageType();
                image_container_[0][c]->create(dim_im, buf+c*numOfPixels, false);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::create(value_type* buf, const std::vector<size_t>& dim) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::clear()
    {
        try
        {
            if ( delete_data_on_destruct_ )
            {
                size_t row = this->rows();

                unsigned int ii, jj;
                for ( ii=0; ii<row; ii++ )
                {
                    size_t col = this->cols(ii);
                    for ( jj=0; jj<col; jj++ )
                    {
                        ImageType* pImg = image_container_[ii][jj];
                        if ( pImg != NULL )
                        {
                            delete pImg;
                            image_container_[ii][jj] = NULL;
                        }
                    }
                }
            }

            image_container_.clear();
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::clear() ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::copyFrom(const Self& a)
    {
        try
        {
            if ( !this->dimensions_equal_container(a) )
            {
                GADGET_CHECK_RETURN_FALSE(this->clear());
                this->delete_data_on_destruct(true);

                GADGET_CHECK_RETURN_FALSE(this->create(a.cols()));
            }

            size_t row = this->rows();

            unsigned int ii, jj;
            for ( ii=0; ii<row; ii++ )
            {
                size_t col = this->cols(ii);
                for ( jj=0; jj<col; jj++ )
                {
                    image_container_[ii][jj]->copyImageInfoAndContent(a(ii, jj));
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::copyFrom(const Self& a) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::fillWithZeros()
    {
        try
        {
            size_t row = this->rows();

            unsigned int ii, jj;
            for ( ii=0; ii<row; ii++ )
            {
                size_t col = this->cols(ii);
                for ( jj=0; jj<col; jj++ )
                {
                    memset(image_container_[ii][jj]->begin(), 0, image_container_[ii][jj]->get_number_of_bytes());
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::fillWithZeros() ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::push_back(ImageType& im, size_t row)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            image_container_[row].push_back(&im);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::push_back(ImageType& im, size_t row) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::push_front(ImageType& im, size_t row)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            image_container_[row].insert(image_container_[row].begin(), &im);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::push_front(ImageType& im, size_t row) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::insert(ImageType& im, size_t row, size_t col)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            GADGET_CHECK_RETURN_FALSE(col<this->cols(row));

            image_container_[row].insert(image_container_[row].begin()+col, &im);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::insert(ImageType& im, size_t row, size_t col) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::pop_back(ImageType*& im, size_t row)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            im = image_container_[row][this->cols(row)-1];
            image_container_[row].pop_back();
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::pop_back(ImageType*& im, size_t row) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::pop_front(ImageType*& im, size_t row)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());

            if ( this->cols(row) == 0 )
            {
                im = NULL;
                return true;
            }

            image_container_[row].erase(image_container_[row].begin());
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::pop_front(ImageType*& im, size_t row) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::remove(ImageType*& im, size_t row, size_t col)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            GADGET_CHECK_RETURN_FALSE(col<this->cols(row));

            im = image_container_[row][col];

            image_container_[row].erase(image_container_[row].begin()+col);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::remove(ImageType*& im, size_t row, size_t col) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::remove(size_t row, size_t col)
    {
        try
        {
            ImageType* im = NULL;
            GADGET_CHECK_RETURN_FALSE(this->remove(im, row, col));
            if( delete_data_on_destruct_ ) delete im;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::remove(size_t row, size_t col) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::expand(size_t newRows)
    {
        try
        {
            size_t row = this->rows();
            if ( newRows > 0 )
            {
                image_container_.resize(row+newRows);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::expand(size_t newRows) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::insert(std::vector<ImageType*>& im_array, size_t row)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            image_container_.insert(image_container_.begin()+row, im_array);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::insert(std::vector<ImageType*>& im_array, size_t row) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::remove(std::vector<ImageType*>& im_array, size_t row)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            im_array = image_container_[row];
            image_container_.erase(image_container_.begin()+row);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::remove(std::vector<ImageType*>& im_array, size_t row) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::remove(size_t row)
    {
        try
        {
            std::vector<ImageType*> im_array;
            GADGET_CHECK_RETURN_FALSE(this->remove(im_array, row));

            if( delete_data_on_destruct_ )
            {
                size_t N = im_array.size();
                unsigned int ii;
                for ( ii=0; ii<N; ii++ )
                {
                    if ( im_array[ii] != NULL )
                    {
                        delete im_array[ii];
                        im_array[ii] = NULL;
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::remove(size_t row) ... ");
            return false;
        }

        return true;
    }

    /// get image pointers
    template <typename ImageType> 
    inline ImageType& hoNDImageContainer2D<ImageType>::get(size_t row, size_t col)
    {
        GADGET_DEBUG_CHECK_THROW(row<this->rows());
        GADGET_DEBUG_CHECK_THROW(col<this->cols(row));

        return *(image_container_[row][col]);
    }

    template <typename ImageType> 
    inline const ImageType& hoNDImageContainer2D<ImageType>::get(size_t row, size_t col) const
    {
        GADGET_DEBUG_CHECK_THROW(row<this->rows());
        GADGET_DEBUG_CHECK_THROW(col<this->cols(row));

        return *(image_container_[row][col]);
    }

    template <typename ImageType> 
    inline ImageType& hoNDImageContainer2D<ImageType>::operator() (size_t row, size_t col)
    {
        GADGET_DEBUG_CHECK_THROW(row<this->rows());
        GADGET_DEBUG_CHECK_THROW(col<this->cols(row));

        return *(image_container_[row][col]);
    }

    template <typename ImageType> 
    inline const ImageType& hoNDImageContainer2D<ImageType>::operator() (size_t row, size_t col) const
    {
        GADGET_DEBUG_CHECK_THROW(row<this->rows());
        GADGET_DEBUG_CHECK_THROW(col<this->cols(row));

        return *(image_container_[row][col]);
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::get(std::vector<ImageType*>& im_array, size_t row) const
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());
            im_array = image_container_[row];
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::get(std::vector<ImageType*>& im_array, size_t row) const ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline size_t hoNDImageContainer2D<ImageType>::get_number_of_all_images()
    {
        try
        {
            size_t num = 0;

            size_t row = this->rows();
            if ( row == 0 ) return num;

            unsigned int r;
            for ( r=0; r<row; r++ )
            {
                num += this->cols(r);
            }

            return num;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::get(std::vector<ImageType*>& im_array, size_t row) ... ");
            return false;
        }
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::get_all_images(std::vector<ImageType*>& im_array)
    {
        try
        {
            im_array.clear();

            size_t row = this->rows();
            if ( row == 0 ) return true;

            size_t num = this->get_number_of_all_images();

            im_array.resize(num, NULL);

            unsigned int r, c, ind(0);
            for ( r=0; r<row; r++ )
            {
                for ( c=0; c<this->cols(r); c++ )
                {
                    im_array[ind] = image_container_[r][c];
                    ind++;
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::get(std::vector<ImageType*>& im_array, size_t row) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::set(ImageType* pImage, size_t row, size_t col)
    {
        GADGET_DEBUG_CHECK_RETURN_FALSE(row<this->rows());
        GADGET_DEBUG_CHECK_RETURN_FALSE(col<this->cols(row));

        if ( image_container_[row][col] != NULL )
        {
            if ( this->delete_data_on_destruct() ) delete image_container_[row][col];
        }

        image_container_[row][col] = pImage;

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::to_NDArray(size_t row, hoNDArray<value_type>& a) const
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(row<this->rows());

            size_t col = this->cols(row);
            if ( col == 0 ) return true;

            GADGET_CHECK_RETURN_FALSE(this->has_identical_dimensions( (unsigned int)row));

            std::vector<size_t> dim;
            image_container_[row][0]->get_dimensions(dim);

            size_t numOfElements = image_container_[row][0]->get_number_of_elements();
            size_t numOfBytes = image_container_[row][0]->get_number_of_bytes();

            std::vector<size_t> dim_out(dim.size()+1);
            memcpy(&dim_out[0], &dim[0], sizeof(size_t)*dim.size());
            dim_out[ dim.size() ] = col;

            a.create(dim_out);

            unsigned int c;
            for ( c=0; c<col; c++ )
            {
                memcpy(a.begin()+c*numOfElements, image_container_[row][c]->begin(), numOfBytes);
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::to_NDArray(size_t row, hoNDArray<value_type>& a) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::delete_data_on_destruct() const
    {
        return this->delete_data_on_destruct_;
    }

    template <typename ImageType> 
    inline void hoNDImageContainer2D<ImageType>::delete_data_on_destruct(bool d)
    {
        this->delete_data_on_destruct_ = d;
    }

    template <typename ImageType> 
    inline size_t hoNDImageContainer2D<ImageType>::rows() const
    {
        return image_container_.size();
    }

    template <typename ImageType> 
    inline size_t hoNDImageContainer2D<ImageType>::cols(size_t row) const
    {
        GADGET_DEBUG_CHECK_THROW(row<this->rows());
        return image_container_[row].size();
    }

    template <typename ImageType> 
    inline std::vector<size_t> hoNDImageContainer2D<ImageType>::cols() const
    {
        std::vector<size_t> col(this->rows(), 0);
        unsigned int row;
        for ( row=0; row<this->rows(); row++ )
        {
            col[row] = this->cols(row);
        }
        return col;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::has_identical_dimensions(unsigned int row) const
    {
        GADGET_CHECK_RETURN_FALSE(row<this->rows());

        size_t col = this->cols(row);
        if ( col == 0 ) return true;

        unsigned int c;
        for ( c=1; c<col; c++ )
        {
            if ( !image_container_[row][0]->dimensions_equal( *image_container_[row][c] ) )
            {
                return false;
            }
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::has_identical_pixel_size(unsigned int row) const
    {
        GADGET_CHECK_RETURN_FALSE(row<this->rows());

        size_t col = this->cols(row);
        if ( col == 0 ) return true;

        unsigned int c;
        for ( c=1; c<col; c++ )
        {
            if ( !image_container_[row][0]->pixel_size_equal( *image_container_[row][c] ) )
            {
                return false;
            }
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::has_identical_axis(unsigned int row) const
    {
        GADGET_CHECK_RETURN_FALSE(row<this->rows());

        size_t col = this->cols(row);
        if ( col == 0 ) return true;

        unsigned int c;
        for ( c=1; c<col; c++ )
        {
            if ( !image_container_[row][0]->axis_equal( *image_container_[row][c] ) )
            {
                return false;
            }
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::has_identical_image_geometry(unsigned int row) const
    {
        GADGET_CHECK_RETURN_FALSE(row<this->rows());

        if ( !this->has_identical_dimensions() ) return false;
        if ( !this->has_identical_pixel_size() ) return false;
        if ( !this->has_identical_axis() ) return false;

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::serialize(char*& buf, size_t& totalLen) const 
    {
        try
        {
            // memory layout
            // number of row, number of col for row 1 (col1), number of col for row 2, ..., number of col for row n
            // offset for image[0][0], len of buffer for image[0][0], offset for image[0][1], len of buffer for image[0][1], ..., offset for image[0][0], len of buffer for image[0][col1-1], 
            // ...
            // offset for image[row-1][0], len of buffer for image[row-1][0], offset for image[row-1][1], len of buffer for image[row-1][1], ..., offset for image[row-1][0], len of buffer for image[row-1][col1-1], 
            // content for image[0][0], ..., image[0][col1-1], image[1][0], ..., image[row-1][col1-1]

            // starting for image content
            size_t offsetImage = sizeof(size_t) + this->rows()*sizeof(size_t);

            std::vector<size_t> col(this->rows());

            std::vector< std::vector<size_t> > offset(this->rows());
            std::vector< std::vector<size_t> > len(this->rows());
            std::vector< std::vector<char*> > bufIm(this->rows());

            size_t row, c;
            for ( row=0; row<this->rows(); row++ )
            {
                col[row] = this->cols(row);
                offset[row].resize(col[row], 0);
                len[row].resize(col[row], 0);

                offsetImage += sizeof(size_t)*col[row]*2;

                bufIm[row].resize(col[row], NULL);
            }

            totalLen = offsetImage;
            offset[0][0] = offsetImage;

            for ( row=0; row<this->rows(); row++ )
            {
                for ( c=0; c<col[row]; c++ )
                {
                    ImageType* im = image_container_[row][c];
                    if ( im != NULL )
                    {
                        char* bufImCurr=NULL;
                        size_t lenIm;

                        im->serialize(bufImCurr, lenIm);

                        bufIm[row][c] = bufImCurr;
                        len[row][c] = lenIm;
                    }
                    else
                    {
                        len[row][c] = 0;
                    }

                    totalLen += len[row][c];

                    if ( row==0 && c== 0 ) continue;

                    offset[row][c] = offset[row][c-1] + len[row][c-1];
                }
            }

            buf = new char[totalLen];
            GADGET_CHECK_RETURN_FALSE(buf!=NULL);

            size_t offsetBuf = 0;

            size_t numOfRows = this->rows();
            memcpy(buf+offsetBuf, &numOfRows, sizeof(size_t));
            offsetBuf += sizeof(size_t);

            memcpy(buf+offsetBuf, &col[0], sizeof(size_t)*numOfRows);
            offsetBuf += sizeof(size_t)*numOfRows;

            for ( row=0; row<this->rows(); row++ )
            {
                for ( c=0; c<col[row]; c++ )
                {
                    size_t v = offset[row][c];
                    size_t lv = len[row][c];

                    memcpy(buf+offsetBuf, &v, sizeof(size_t));
                    offsetBuf += sizeof(size_t);

                    memcpy(buf+offsetBuf, &lv, sizeof(size_t));
                    offsetBuf += sizeof(size_t);
                }
            }

            for ( row=0; row<this->rows(); row++ )
            {
                for ( c=0; c<col[row]; c++ )
                {
                    if ( bufIm[row][c] != NULL )
                    {
                        memcpy(buf+offsetBuf, bufIm[row][c], len[row][c]);
                        offsetBuf += len[row][c];

                        delete [] bufIm[row][c];
                        bufIm[row][c] = NULL;
                    }
                }
            }

            GADGET_CHECK_RETURN_FALSE(totalLen == offsetBuf);
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::serialize(char*& buf, size_t& len) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline bool hoNDImageContainer2D<ImageType>::deserialize(char* buf, size_t& len)
    {
        try
        {
            GADGET_CHECK_RETURN_FALSE(this->clear());
            this->delete_data_on_destruct(true);

            size_t offsetBuf = 0;

            size_t numOfRows(0);
            memcpy(&numOfRows, buf+offsetBuf, sizeof(size_t));
            offsetBuf += sizeof(size_t);

            if ( numOfRows == 0 ) return true;

            image_container_.resize(numOfRows);

            std::vector<size_t> col(numOfRows);

            memcpy(&col[0], buf+offsetBuf, sizeof(size_t)*numOfRows);
            offsetBuf += sizeof(size_t)*numOfRows;

            size_t row, c;
            for ( row=0; row<this->rows(); row++ )
            {
                image_container_[row].resize(col[row], NULL);

                for ( c=0; c<col[row]; c++ )
                {
                    size_t offsetCurr, lenCurr;

                    memcpy(&offsetCurr, buf+offsetBuf, sizeof(size_t));
                    offsetBuf += sizeof(size_t);

                    memcpy(&lenCurr, buf+offsetBuf, sizeof(size_t));
                    offsetBuf += sizeof(size_t);

                    image_container_[row][c] = new ImageType();
                    GADGET_CHECK_RETURN_FALSE(image_container_[row][c]!=NULL);

                    if ( lenCurr > 0 )
                    {
                        GADGET_CHECK_RETURN_FALSE(image_container_[row][c]->deserialize(buf+offsetCurr, lenCurr));
                    }
                }
            }
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageContainer2D<ImageType>::deserialize(char* buf, size_t& len) ... ");
            return false;
        }

        return true;
    }

    template <typename ImageType> 
    inline void hoNDImageContainer2D<ImageType>::print(std::ostream& os) const
    {
        using namespace std;

        os.unsetf(std::ios::scientific);
        os.setf(ios::fixed);

        size_t r, c;

        os << "--------------Gagdgetron Image Container 2D -------------" << endl;
        os << "Image type is : " << std::string(typeid(ImageType).name()) << endl;
        os << "Number of stored image rows is : " << this->rows() << endl;
        for ( r=0; r<this->rows(); r++ )
        {
            os << "Row " << r << " has " << this->cols(r) << " images " << endl;
        }
        os << "---------------------------------------------------------" << endl;
        for ( r=0; r<this->rows(); r++ )
        {
            os << "Row " << r << " : "<< endl;
            os << "=========================================================" << endl;
            for ( c=0; c<this->cols(r); c++ )
            {
                if ( c > 2 ) break;

                if ( image_container_[r][c] != NULL )
                {
                    os << "--> Image " << c << " : "<< endl;
                    image_container_[r][c]->print(os);
                    os << "=========================================================" << endl;
                }
            }
        }
    }
}
