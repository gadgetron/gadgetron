/** \file       hoNDMetaAttributes.h
    \brief      meta attributes class for the gadgetron

                The field key is a string. The key-value pairs are stored in the std::map.
                The sereialization and deserialization of meta info is using the xml format.

    \author     Hui Xue
*/

#pragma once

#include "cpucore_export.h"
#include "hoNDMetaAttributesBase.h"

//#include <ace/Task.h>
//#include <ace/SOCK_Stream.h>
#include <string>
#include <sstream>
#include <vector>
#include <complex>

namespace Gadgetron
{
    template <typename TInteger, typename TFloat, typename TComplexFloat>
    class EXPORTCPUCORE hoNDMetaAttributes
    {
    public:

        typedef hoNDMetaAttributes<TInteger, TFloat, TComplexFloat> Self;

        typedef typename TComplexFloat::value_type cx_value_type;

        typedef hoNDMetaAttributesBase<TInteger> AttributeType1;
        typedef typename AttributeType1::AttributeStorageType AttributeStorageType1;

        typedef hoNDMetaAttributesBase<TFloat> AttributeType2;
        typedef typename AttributeType2::AttributeStorageType AttributeStorageType2;

        typedef hoNDMetaAttributesBase<TComplexFloat> AttributeType3;
        typedef typename AttributeType3::AttributeStorageType AttributeStorageType3;

        typedef hoNDMetaAttributesBase<std::string> AttributeType4;
        typedef typename AttributeType4::AttributeStorageType AttributeStorageType4;

        typedef typename AttributeType1::size_t_type size_t_type;

        typedef TInteger AttribIntegerType;
        typedef TFloat AttribFloatType;
        typedef TComplexFloat AttribComplexType;

        hoNDMetaAttributes();
        virtual ~hoNDMetaAttributes();

        hoNDMetaAttributes(const Self& attrib);

        Self& operator=(const Self& attrib)
        {
            if ( this == &attrib ) return *this;
            attribute1_ = attrib.attribute1_;
            attribute2_ = attrib.attribute2_;
            attribute3_ = attrib.attribute3_;
            attribute4_ = attrib.attribute4_;
            return *this;
        }

        bool clear()
        {
            GADGET_CHECK_RETURN_FALSE(attribute1_.clear());
            GADGET_CHECK_RETURN_FALSE(attribute2_.clear());
            GADGET_CHECK_RETURN_FALSE(attribute3_.clear());
            GADGET_CHECK_RETURN_FALSE(attribute4_.clear());
            return true;
        }

        /// swap the attributes
        bool swap(Self& attrib)
        {
            GADGET_CHECK_RETURN_FALSE(attribute1_.swap(attrib.attribute1_));
            GADGET_CHECK_RETURN_FALSE(attribute2_.swap(attrib.attribute2_));
            GADGET_CHECK_RETURN_FALSE(attribute3_.swap(attrib.attribute3_));
            GADGET_CHECK_RETURN_FALSE(attribute4_.swap(attrib.attribute4_));
            return true;
        }

        // serialize and deserialize to/from the buffer
        virtual bool serialize(char*& buf, size_t_type& len);
        virtual bool deserialize(char* buf, size_t_type& len);
        // deserialize the content of attributes
        virtual bool deserializeContent(char* buf, size_t_type& xmlLen);

        /// print out the image information
        virtual void print(std::ostream& os) const;

        /// all other operations should be performed directly on attribute itself
        AttributeType1 attribute1_;
        AttributeType2 attribute2_;
        AttributeType3 attribute3_;
        AttributeType4 attribute4_;
    private:

        // ACE_Thread_Mutex mtx_;
    };

    typedef hoNDMetaAttributes<long long, float, std::complex<float> > GtImageAttribType;
}
