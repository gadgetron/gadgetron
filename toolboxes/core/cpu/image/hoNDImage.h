/** \file       hoNDImage.h
    \brief      N-dimensional image class for gadgetron

                The default N-dimensional image is defined by the origin (the first pixel indexed by [0 0 0 ...]),
                the pixel size and the axis for every coordinate. This defines an Euclidean space.

                If this N-dimensional image is used with other coordinate systems, e.g. polar coordinate system, then the axis 
                should not be used to compute the image-to-world transformation.

    \author     Hui Xue
*/

#pragma once

#include "hoNDPoint.h"
#include "hoMatrix.h"

namespace Gadgetron
{

    template <typename T, unsigned int D>
    class hoNDImage : public hoNDArray<T>
    {
    public:

        typedef hoNDArray<T> BaseClass;
        typedef hoNDImage<T, D> Self;

        typedef T element_type;
        typedef T value_type;
        typedef double coord_type;

        typedef hoNDPoint<coord_type, D> a_axis_type;
        typedef std::vector<a_axis_type> axis_type;

        typedef hoNDPoint<coord_type, 3> a_axis_image_patient_type;

        enum { NDIM = D };

        void* operator new (size_t bytes)
        {
            return ::new char[bytes];
        }

        void operator delete (void *ptr)
        {
            delete [] static_cast <char *> (ptr);
        } 

        void * operator new(size_t s, void * p)
        {
            return p;
        }

        /// constructors
        hoNDImage ();
        hoNDImage (const std::vector<size_t>& dimensions);
        hoNDImage (boost::shared_ptr< std::vector<size_t> > dimensions);
        hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize);
        hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin);
        hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis);

        hoNDImage(size_t len);
        hoNDImage(size_t sx, size_t sy);
        hoNDImage(size_t sx, size_t sy, size_t sz);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss);

        /// attach memory constructors
        hoNDImage (const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);
        hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct = false);
        hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct = false);
        hoNDImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct = false);

        hoNDImage(size_t len, T* data, bool delete_data_on_destruct = false);
        hoNDImage(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);
        hoNDImage(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct = false);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct = false);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct = false);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct = false);
        hoNDImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct = false);

        hoNDImage(const hoNDArray<T>& a);

        hoNDImage(const Self& a);
        Self& operator=(const Self& rhs);

        virtual ~hoNDImage();

        /// clear the images, release all memory it holds, set pixelsize/axis/origin to zero-status
        void clear();

        /// create the image, called by constructors
        virtual void create(const std::vector<size_t>& dimensions);
        virtual void create(boost::shared_ptr< std::vector<size_t> > dimensions);
        virtual void create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize);
        virtual void create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin);
        virtual void create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis);

        /// create the image from another image
        /// not copy its content
        template<typename T2> 
        void createFrom(const hoNDImage<T2, D>& im)
        {
            this->clear();

            std::vector<size_t> dim;
            im.get_dimensions(dim);

            std::vector<coord_type> pixelSize;
            im.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            im.get_origin(origin);

            axis_type axis;
            im.get_axis(axis);

            this->create(dim, pixelSize, origin, axis);
        }

        /// create the image from another image
        /// copy its content
        template<typename T2> 
        void create(const hoNDImage<T2, D>& im)
        {
            this->createFrom(im);

            size_t ii;
            size_t N = this->get_number_of_elements();
            for ( ii=0; ii<N; ii++ )
            {
                this->data_[ii] = static_cast<T>(im.get_data_ptr()[ii]);
            }
        }

        template<typename T2> 
        inline void copyImageInfo(const hoNDImage<T2, D>& im)
        {
            this->createFrom(im);
        }

        template<typename T2> 
        inline void copyImageInfoAndContent(const hoNDImage<T2, D>& im)
        {
            this->create(im);
        }

        template<typename T2> 
        inline void copyImageInfoWithoutImageSize(const hoNDImage<T2, D>& im)
        {
            std::vector<coord_type> pixelSize;
            im.get_pixel_size(pixelSize);

            std::vector<coord_type> origin;
            im.get_origin(origin);

            axis_type axis;
            im.get_axis(axis);

            this->set_pixel_size(pixelSize);
            this->set_origin(origin);
            this->set_axis(axis);
        }

        virtual void create(const std::vector<size_t>& dimensions,
                            T* data, 
                            bool delete_data_on_destruct = false);

        virtual void create(const std::vector<size_t>& dimensions, 
                            const std::vector<coord_type>& pixelSize, 
                            T* data, 
                            bool delete_data_on_destruct = false);

        virtual void create(const std::vector<size_t>& dimensions, 
                            const std::vector<coord_type>& pixelSize, 
                            const std::vector<coord_type>& origin, 
                            T* data, 
                            bool delete_data_on_destruct = false);

        virtual void create(const std::vector<size_t>& dimensions, 
                            const std::vector<coord_type>& pixelSize, 
                            const std::vector<coord_type>& origin, 
                            const axis_type& axis, 
                            T* data, 
                            bool delete_data_on_destruct = false);

        /// convert from/to hoNDArray
        void from_NDArray(const hoNDArray<T>& a);
        void to_NDArray(hoNDArray<T>& a) const;

        /// whether two images have the same size
        bool dimensions_equal(const std::vector<size_t>& dimensions) const;

        template<class S> 
        bool dimensions_equal(const hoNDArray<S>& im) const
        {
            std::vector<size_t> dim;
            im.get_dimensions(dim);

            return this->dimensions_equal(dim);
        }

        template<class S> 
        bool dimensions_equal(const hoNDImage<S, D>& im) const
        {
            std::vector<size_t> dim;
            im.get_dimensions(dim);

            return this->dimensions_equal(dim);
        }

        template<class S> 
        bool dimensions_equal(const hoNDImage<S, D>* im) const
        {
            return this->dimensions_equal(*im);
        }

        template<class S> 
        bool pixel_size_equal(const hoNDImage<S, D>& im) const
        {
            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                if ( std::abs(this->pixelSize_[ii] - im.pixelSize_[ii]) > FLT_EPSILON ) return false;
            }

            return true;
        }

        template<class S> 
        bool axis_equal(const hoNDImage<S, D>& im) const
        {
            unsigned int ii;
            for ( ii=0; ii<D; ii++ )
            {
                if ( this->axis_[ii] != im.axis_[ii] ) return false;
            }

            return true;
        }

        /// get the pixel size
        coord_type get_pixel_size(size_t dimension) const;
        void get_pixel_size(std::vector<coord_type>& pixelSize) const;

        void set_pixel_size(size_t dimension, coord_type v);
        void set_pixel_size(const std::vector<coord_type>& pixelSize);

        /// get origin
        coord_type get_origin(size_t dimension) const;
        void get_origin(std::vector<coord_type>& origin) const;

        void set_origin(size_t dimension, coord_type v);
        void set_origin(const std::vector<coord_type>& origin);

        /// get axis
        coord_type get_axis(size_t dimension, size_t elem) const;
        a_axis_type get_axis(size_t dimension) const;
        void get_axis(axis_type& axis) const;

        void set_axis(size_t dimension, size_t elem, coord_type v);
        void set_axis(size_t dimension, const a_axis_type& v);
        void set_axis(const axis_type& axis);

        /// get image position patient
        void get_image_position(coord_type pos[3]) const;
        void get_image_position(unsigned int d, coord_type& pos) const;
        void get_image_position(a_axis_image_patient_type& pos) const;

        void set_image_position(coord_type pos[3]);
        void set_image_position(unsigned int d, coord_type pos);
        void set_image_position(const a_axis_image_patient_type& pos);

        /// get image orientation patient
        void get_image_orientation(unsigned int d, coord_type ori[3]) const;
        void get_image_orientation(unsigned int d, a_axis_image_patient_type& ori) const;
        
        /// for dimension d and index ind
        void get_image_orientation(unsigned int d, unsigned int ind, coord_type& ori) const;
        /// get image orientation as a quaternion
        void get_image_orientation(coord_type quat[4]) const;

        void set_image_orientation(unsigned int d, coord_type ori[3]);
        void set_image_orientation(unsigned int d, const a_axis_image_patient_type& ori);
        void set_image_orientation(unsigned int d, unsigned int ind, coord_type ori);
        void set_image_orientation(coord_type quat[4]);

        size_t get_number_of_dimensions() const { return D; }

        size_t calculate_offset(const size_t* ind) const;
        size_t calculate_offset(const std::vector<size_t>& ind) const;

        size_t calculate_offset(size_t x, size_t y) const;
        size_t calculate_offset(size_t x, size_t y, size_t z) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q) const;
        size_t calculate_offset(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u) const;

        /// given the 1D offset, compute the corresponding indexes
        std::vector<size_t> calculate_index( size_t offset ) const;
        void calculate_index( size_t offset, size_t* index ) const;
        void calculate_index( size_t offset, std::vector<size_t>& index ) const;
        void calculate_index( size_t offset, coord_type* index ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r, size_t& a ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r, size_t& a, size_t& q ) const;
        void calculate_index( size_t offset, size_t& x, size_t& y, size_t& z, size_t& s, size_t& p, size_t& r, size_t& a, size_t& q, size_t& u ) const;

        /// access the pixel value
        T& operator()( const size_t* ind );
        const T& operator()( const size_t* ind ) const;

        T& operator()( const std::vector<size_t>& ind );
        const T& operator()( const std::vector<size_t>& ind ) const;

        T& operator[]( size_t x );
        const T& operator[]( size_t x ) const;

        T& operator()( size_t x );
        const T& operator()( size_t x ) const;

        T& operator()( size_t x, size_t y );
        const T& operator()( size_t x, size_t y ) const;

        T& operator()( size_t x, size_t y, size_t z );
        const T& operator()( size_t x, size_t y, size_t z ) const;

        T& operator()( size_t x, size_t y, size_t z, size_t s );
        const T& operator()( size_t x, size_t y, size_t z, size_t s ) const;

        T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p );
        const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p ) const;

        T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r );
        const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r ) const;

        T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a );
        const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a ) const;

        T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q );
        const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q ) const;

        T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u );
        const T& operator()( size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u ) const;

        /// fill the image with a value
        void fill(T value);

        template<typename T2> 
        void copyFrom(const hoNDImage<T2, D>& aIm)
        {
            this->create(aIm);
        }

        /// image pixel index to world coordinate
        void image_to_world(const coord_type* ind, coord_type* coord) const;
        void image_to_world(const std::vector<coord_type>& ind, std::vector<coord_type>& coord) const;

        void image_to_world(coord_type x, coord_type& cx) const;

        void image_to_world(coord_type x, coord_type y, 
                            coord_type& cx, coord_type& cy) const;

        void image_to_world(coord_type x, coord_type y, coord_type z,
                            coord_type& cx, coord_type& cy, coord_type& cz) const;

        void image_to_world(coord_type x, coord_type y, coord_type z, coord_type s,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs) const;

        void image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp) const;

        void image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr) const;

        void image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca) const;

        void image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq) const;

        void image_to_world(coord_type x, coord_type y, coord_type z, coord_type s, coord_type p, coord_type r, coord_type a, coord_type q, coord_type u,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq, coord_type& cu) const;

        /// for integer pixel indexes
        void image_to_world(const size_t* ind, coord_type* coord) const;
        void image_to_world(const std::vector<size_t>& ind, std::vector<coord_type>& coord) const;

        void image_to_world(size_t x, coord_type& cx) const;

        void image_to_world(size_t x, size_t y, 
                            coord_type& cx, coord_type& cy) const;

        void image_to_world(size_t x, size_t y, size_t z,
                            coord_type& cx, coord_type& cy, coord_type& cz) const;

        void image_to_world(size_t x, size_t y, size_t z, size_t s,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs) const;

        void image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp) const;

        void image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr) const;

        void image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca) const;

        void image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq) const;

        void image_to_world(size_t x, size_t y, size_t z, size_t s, size_t p, size_t r, size_t a, size_t q, size_t u,
                            coord_type& cx, coord_type& cy, coord_type& cz, coord_type& cs, coord_type& cp, coord_type& cr, coord_type& ca, coord_type& cq, coord_type& cu) const;

        /// get the image-to-world transformation matrix
        /// the Homogeneous coordinate transformation matrix is computed
        void image_to_world_matrix(hoMatrix<coord_type>& image2world) const;
        void set_image_to_world_matrix(const hoMatrix<coord_type>& image2world);

        /// world coordinate to image pixel index
        void world_to_image(const coord_type* coord, coord_type* ind) const;
        void world_to_image(const std::vector<coord_type>& coord, std::vector<coord_type>& ind) const;

        void world_to_image(coord_type cx, coord_type& x) const;

        void world_to_image(coord_type cx, coord_type cy, 
                            coord_type& x, coord_type& y) const;

        void world_to_image(coord_type cx, coord_type cy, coord_type cz,
                            coord_type& x, coord_type& y, coord_type& z) const;

        void world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs,
                            coord_type& x, coord_type& y, coord_type& z, coord_type& s) const;

        void world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp,
                            coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p) const;

        void world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr,
                            coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r) const;

        void world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr, coord_type ca,
                            coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r, coord_type& a) const;

        void world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr, coord_type ca, coord_type cq,
                            coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r, coord_type& a, coord_type& q) const;

        void world_to_image(coord_type cx, coord_type cy, coord_type cz, coord_type cs, coord_type cp, coord_type cr, coord_type ca, coord_type cq, coord_type cu,
                            coord_type& x, coord_type& y, coord_type& z, coord_type& s, coord_type& p, coord_type& r, coord_type& a, coord_type& q, coord_type& u) const;

        /// get the world_to_image transformation matrix
        /// the Homogeneous coordinate transformation matrix is computed
        void world_to_image_matrix(hoMatrix<coord_type>& world2image) const;
        void set_world_to_image_matrix(const hoMatrix<coord_type>& world2image);

        /// is the sub region in the image
        bool in_image_region(const std::vector<size_t>& start, std::vector<size_t>& size);

        /// get the sub image
        void get_sub_image(const std::vector<size_t>& start, std::vector<size_t>& size, Self& out);

        /// serialize/deserialize image content
        virtual bool serialize(char*& buf, size_t& len) const;
        virtual bool deserialize(char* buf, size_t& len);

        /// print out the image information
        virtual void print(std::ostream& os) const;
        virtual void printContent(std::ostream& os) const;

    protected:

        using BaseClass::dimensions_;
        using BaseClass::offsetFactors_;
        using BaseClass::data_;
        using BaseClass::elements_;
        using BaseClass::delete_data_on_destruct_;

        coord_type pixelSize_[D];
        coord_type pixelSize_reciprocal_[D];
        coord_type origin_[D];
        hoNDPoint<coord_type, D> axis_[D];

        /// for the dicom coordinate system
        a_axis_image_patient_type image_position_patient_;
        /// image orientation for row/column/slice directions
        a_axis_image_patient_type image_orientation_patient_[3];
    };
}

#include "hoNDImage.hxx"
