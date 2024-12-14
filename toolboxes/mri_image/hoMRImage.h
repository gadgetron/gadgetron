/** \file       hoMRImage.h
    \brief      N-dimensional image class for gadgetron MRI

                This image class is derived from the hoMRImage and attached with MRD header and attributes
                The purpose of this class is to provide an easy-to-use image class for MR recon and image processing

    \author     Hui Xue
*/

#pragma once

#include "hoNDImage.h"
#include "io/primitives.h"
#include "mrd/types.h"

namespace Gadgetron
{

    template <typename T, unsigned int D>
    class hoMRImage : public hoNDImage<T, D>
    {
    public:

        typedef hoNDImage<T, D> BaseClass;
        typedef hoMRImage<T, D> Self;

        typedef T element_type;
        typedef T value_type;
        typedef typename BaseClass::coord_type coord_type;
        typedef typename BaseClass::a_axis_type a_axis_type;
        typedef typename BaseClass::axis_type axis_type;
        typedef typename BaseClass::a_axis_image_patient_type a_axis_image_patient_type;

        enum { NDIM = D };

        /// constructors
        hoMRImage ();
        hoMRImage (const std::vector<size_t>& dimensions);
        hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize);
        hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin);
        hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis);

        hoMRImage(size_t len);
        hoMRImage(size_t sx, size_t sy);
        hoMRImage(size_t sx, size_t sy, size_t sz);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss);

        /// attach memory constructors
        hoMRImage (const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);
        hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct = false);
        hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct = false);
        hoMRImage (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct = false);

        hoMRImage(size_t len, T* data, bool delete_data_on_destruct = false);
        hoMRImage(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);
        hoMRImage(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct = false);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct = false);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct = false);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct = false);
        hoMRImage(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct = false);

        hoMRImage(const hoNDArray<T>& a);

        hoMRImage(const Self& a);
        Self& operator=(const Self& rhs);

        virtual ~hoMRImage();

        /// clear the images, release all memory it holds, set pixelsize/axis/origin to zero-status
        void clear();

        /// create the image, called by constructors
        virtual void create(const std::vector<size_t>& dimensions);
        virtual void create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize);
        virtual void create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin);
        virtual void create(const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis);

        /// create the image from another image
        /// not copy its content
        template<typename T2>
        void createFrom(const hoMRImage<T2, D>& im)
        {
            BaseClass::createFrom(im);

            this->header_ = im.header_;
            this->attrib_ = im.attrib_;
        }

        /// create the image from another image
        /// copy its content
        template<typename T2>
        void create(const hoMRImage<T2, D>& im)
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
        inline void copyImageInfo(const hoMRImage<T2, D>& im)
        {
            this->createFrom(im);
        }

        template<typename T2>
        inline void copyImageInfoAndContent(const hoMRImage<T2, D>& im)
        {
            this->create(im);
        }

        template<typename T2>
        inline void copyImageInfoWithoutImageSize(const hoMRImage<T2, D>& im)
        {
            BaseClass::copyImageInfoWithoutImageSize(im);

            this->header_ = im.header_;
            this->attrib_ = im.attrib_;
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

        template<typename T2>
        void copyFrom(const hoMRImage<T2, D>& aIm)
        {
            this->create(aIm);
        }

        /// get the sub image
        void get_sub_image(const std::vector<size_t>& start, std::vector<size_t>& size, Self& out);

        /// mrd image header structure
        mrd::ImageHeader header_;

        /// meta attributes
        mrd::ImageMeta attrib_;

        /// print out the image information
        virtual void printContent(std::ostream& os) const;

    protected:

        using BaseClass::dimensions_;
        using BaseClass::offsetFactors_;
        using BaseClass::data_;
        using BaseClass::elements_;
        using BaseClass::delete_data_on_destruct_;

        using BaseClass::pixelSize_;
        using BaseClass::pixelSize_reciprocal_;
        using BaseClass::origin_;
        using BaseClass::axis_;

        using BaseClass::image_position_patient_;
        using BaseClass::image_orientation_patient_;
    };
}

#include "hoMRImage.hxx"
