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

        typedef hoNDMetaAttributesBase<TInteger> AttributeTypeInteger;
        typedef typename AttributeTypeInteger::AttributeStorageType AttributeTypeStringAttributeTypeInteger;

        typedef hoNDMetaAttributesBase<TFloat> AttributeTypeFloat;
        typedef typename AttributeTypeFloat::AttributeStorageType AttributeTypeStringAttributeTypeFloat;

        typedef hoNDMetaAttributesBase<TComplexFloat> AttributeTypeComplexFloat;
        typedef typename AttributeTypeComplexFloat::AttributeStorageType AttributeTypeStringAttributeTypeComplexFloat;

        typedef hoNDMetaAttributesBase<std::string> AttributeTypeString;
        typedef typename AttributeTypeString::AttributeStorageType AttributeTypeStringAttributeTypeString;

        typedef typename AttributeTypeInteger::size_t_type size_t_type;

        typedef TInteger AttribIntegerType;
        typedef TFloat AttribFloatType;
        typedef TComplexFloat AttribComplexType;

        hoNDMetaAttributes();
        virtual ~hoNDMetaAttributes();

        hoNDMetaAttributes(const Self& attrib);

        Self& operator=(const Self& attrib)
        {
            if ( this == &attrib ) return *this;
            attributeInteger_ = attrib.attributeInteger_;
            attributeFloat_ = attrib.attributeFloat_;
            attributeComplexFloat_ = attrib.attributeComplexFloat_;
            attributeString_ = attrib.attributeString_;
            return *this;
        }

        bool clear()
        {
            GADGET_CHECK_RETURN_FALSE(attributeInteger_.clear());
            GADGET_CHECK_RETURN_FALSE(attributeFloat_.clear());
            GADGET_CHECK_RETURN_FALSE(attributeComplexFloat_.clear());
            GADGET_CHECK_RETURN_FALSE(attributeString_.clear());
            return true;
        }

        /// swap the attributes
        bool swap(Self& attrib)
        {
            GADGET_CHECK_RETURN_FALSE(attributeInteger_.swap(attrib.attributeInteger_));
            GADGET_CHECK_RETURN_FALSE(attributeFloat_.swap(attrib.attributeFloat_));
            GADGET_CHECK_RETURN_FALSE(attributeComplexFloat_.swap(attrib.attributeComplexFloat_));
            GADGET_CHECK_RETURN_FALSE(attributeString_.swap(attrib.attributeString_));
            return true;
        }

        // serialize and deserialize to/from the buffer
        virtual bool serialize(char*& buf, size_t_type& len) const;
        virtual bool deserialize(char* buf, size_t_type& len);
        // deserialize the content of attributes
        virtual bool deserializeContent(char* buf, size_t_type& xmlLen);

        /// print out the image information
        virtual void print(std::ostream& os) const;

        /// all other operations should be performed directly on attribute itself
        AttributeTypeInteger attributeInteger_;
        AttributeTypeFloat attributeFloat_;
        AttributeTypeComplexFloat attributeComplexFloat_;
        AttributeTypeString attributeString_;
    private:

        // ACE_Thread_Mutex mtx_;
    };

    typedef hoNDMetaAttributes<long long, float, std::complex<float> > GtImageAttribType;
}
