/** \file       hoNDImageAttrib.h
    \brief      N-dimensional image class for gadgetron with meta attributes

                The serialize and deserialize function includes the meta attribute structure as well
                The image data are first serialized, followed by the xml meta attribute representation

    \author     Hui Xue
*/

#pragma once

#include "hoNDImage.h"
#include "hoNDMetaAttributes.h"

namespace Gadgetron
{
    template <typename T, unsigned int D>
    class hoNDImageAttrib : public hoNDImage<T, D>
    {
    public:

        typedef hoNDImage<T, D> BaseClass;
        typedef hoNDImageAttrib<T, D> Self;

        typedef T element_type;
        typedef T value_type;
        typedef float coord_type;

        typedef typename BaseClass::a_axis_type a_axis_type;
        typedef typename BaseClass::axis_type axis_type;

        /// constructors
        hoNDImageAttrib ();
        hoNDImageAttrib (const std::vector<size_t>& dimensions);
        hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize);
        hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin);
        hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis);

        hoNDImageAttrib(size_t len);
        hoNDImageAttrib(size_t sx, size_t sy);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss);

        /// attach memory constructors
        hoNDImageAttrib (const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct = false);

        hoNDImageAttrib(size_t len, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib(size_t sx, size_t sy, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct = false);
        hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct = false);

        hoNDImageAttrib(const hoNDArray<T>& a);
        hoNDImageAttrib(const Self& a);

        virtual ~hoNDImageAttrib();

        /// meta attributes
        GtImageAttribType attrib_;

        /// serialize/deserialize
        virtual bool serialize(char*& buf, size_t& len);
        virtual bool deserialize(char* buf, size_t& len);

        /// print out the image information
        virtual void print(std::ostream& os) const;
        virtual void printContent(std::ostream& os) const;

    protected:

        using BaseClass::dimensions_;
        using BaseClass::offsetFactors_;
        using BaseClass::pixelSize_;
        using BaseClass::pixelSize_reciprocal_;
        using BaseClass::origin_;
        using BaseClass::axis_;
        using BaseClass::data_;
        using BaseClass::elements_;
        using BaseClass::delete_data_on_destruct_;
    };

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib () : BaseClass()
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions) : BaseClass(dimensions)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize) : BaseClass(dimensions, pixelSize)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin) : BaseClass(dimensions, pixelSize, origin)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis) : BaseClass(dimensions, pixelSize, origin, axis)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t len) : BaseClass(len)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy) : BaseClass(sx, sy)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz) : BaseClass(sx, sy, sz)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st) : BaseClass(sx, sy, sz, st)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp) : BaseClass(sx, sy, sz, st, sp)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq) : BaseClass(sx, sy, sz, st, sp, sq)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr) : BaseClass(sx, sy, sz, st, sp, sq, sr)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss) : BaseClass(sx, sy, sz, st, sp, sq, sr, ss)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, pixelSize, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, pixelSize, origin, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib (const std::vector<size_t>& dimensions, const std::vector<coord_type>& pixelSize, const std::vector<coord_type>& origin, const axis_type& axis, T* data, bool delete_data_on_destruct) : BaseClass(dimensions, pixelSize, origin, axis, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t len, T* data, bool delete_data_on_destruct) : BaseClass(len, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, sr, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(size_t sx, size_t sy, size_t sz, size_t st, size_t sp, size_t sq, size_t sr, size_t ss, T* data, bool delete_data_on_destruct) : BaseClass(sx, sy, sz, st, sp, sq, sr, ss, data, delete_data_on_destruct)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(const hoNDArray<T>& a) : BaseClass(a)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::hoNDImageAttrib(const Self& a) : BaseClass(a)
    {
    }

    template <typename T, unsigned int D> 
    hoNDImageAttrib<T, D>::~hoNDImageAttrib()
    {
    }

    template <typename T, unsigned int D> 
    bool hoNDImageAttrib<T, D>::serialize(char*& buf, size_t& len) 
    {
        char* bufImage = NULL;
        char* bufAttrib = NULL;

        try
        {
            size_t lenImage(0);
            GADGET_CHECK_THROW(BaseClass::serialize(bufImage, lenImage));

            size_t lenAttrib(0);
            GADGET_CHECK_THROW(attrib_.serialize(bufAttrib, lenAttrib));

            len = sizeof(unsigned long long) + lenImage + sizeof(unsigned long long) + lenAttrib;

            if ( buf != NULL )
            {
                delete [] buf;
                buf = NULL;
            }

            buf = new char[len];
            GADGET_CHECK_THROW(buf != NULL);

            size_t offset = 0;
            memcpy(buf, &lenImage, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(buf+offset, bufImage, lenImage);
            offset += lenImage;

            memcpy(buf+offset, &lenAttrib, sizeof(size_t));
            offset += sizeof(size_t);

            memcpy(buf+offset, bufAttrib, lenAttrib);
            offset += lenAttrib;

            if ( bufImage != NULL ) delete [] bufImage;
            if ( bufAttrib != NULL ) delete [] bufAttrib;
        }
        catch(...)
        {
            if ( bufImage != NULL ) delete [] bufImage;
            if ( bufAttrib != NULL ) delete [] bufAttrib;

            GERROR_STREAM("Errors happened in hoNDImageAttrib<T, D>::serialize(char*& buf, size_t& len) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    bool hoNDImageAttrib<T, D>::deserialize(char* buf, size_t& len)
    {
        try
        {
            size_t lenImage(0);
            size_t lenAttrib(0);

            size_t offset = 0;
            memcpy(&lenImage, buf, sizeof(size_t));
            offset += sizeof(size_t);

            GADGET_CHECK_RETURN_FALSE(BaseClass::deserialize(buf+offset, lenImage));
            offset += lenImage;

            memcpy(&lenAttrib, buf+offset, sizeof(size_t));
            offset += sizeof(size_t);

            GADGET_CHECK_RETURN_FALSE(attrib_.deserialize(buf+offset, lenAttrib));
            offset += lenAttrib;
        }
        catch(...)
        {
            GERROR_STREAM("Errors happened in hoNDImageAttrib<T, D>::deserialize(char* buf, size_t& len) ... ");
            return false;
        }

        return true;
    }

    template <typename T, unsigned int D> 
    inline void hoNDImageAttrib<T, D>::print(std::ostream& os) const
    {
        using namespace std;
        os << "-------------- Gagdgetron ND Image with meta attributes -------------" << endl;
        this->printContent(os);
    }

    template <typename T, unsigned int D> 
    inline void hoNDImageAttrib<T, D>::printContent(std::ostream& os) const
    {
        BaseClass::printContent(os);
        attrib_.print(os);
    }
}
