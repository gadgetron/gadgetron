/** \file       hoNDMetaAttributes.cpp
    \brief      meta attributes class for the gadgetron

                The field key is a string. The key-value pairs are stored in the std::map.
                The sereialization and deserialization of meta info is using the xml format.

    \author     Hui Xue
*/

#include "hoNDMetaAttributes.h"

#include "url_encode.h"
#include "gadgetronMetaAttributes.hxx"

namespace Gadgetron
{
    template <typename TInteger, typename TFloat, typename TComplexFloat> 
    hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::hoNDMetaAttributes()
    {
    }

    template <typename TInteger, typename TFloat, typename TComplexFloat> 
    hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::~hoNDMetaAttributes()
    {
    }

    template <typename TInteger, typename TFloat, typename TComplexFloat> 
    hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::hoNDMetaAttributes(const Self& attrib)
    {
        *this = attrib;
    }

    template <typename TInteger, typename TFloat, typename TComplexFloat>
    bool hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::serialize(char*& buf, size_t_type& len) const 
    {
        // ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, mtx_, false);

        try
        {
            //mutex_.lock();

            gadgetron::gadgetronMetaAttributes attrib;

            // attribute1
            const AttributeTypeStringAttributeTypeInteger& stroage1 = attributeInteger_.get_attrib_storage();
            typename AttributeTypeStringAttributeTypeInteger::const_iterator iter1;
            for ( iter1=stroage1.begin(); iter1!=stroage1.end(); iter1++ )
            {
                gadgetron::AttributeInteger item(iter1->first);

                gadgetron::AttributeInteger::value_sequence s;

                size_t_type n;
                for ( n=0; n<iter1->second.size(); n++ )
                {
                    s.push_back(iter1->second[n]);
                }

                item.value(s);

                attrib.AttributeInteger().push_back(item);
            }

            // attribute2
            const AttributeTypeStringAttributeTypeFloat& stroage2 = attributeFloat_.get_attrib_storage();
            typename AttributeTypeStringAttributeTypeFloat::const_iterator iter2;
            for ( iter2=stroage2.begin(); iter2!=stroage2.end(); iter2++ )
            {
                gadgetron::AttributeFloat item(iter2->first);

                gadgetron::AttributeFloat::value_sequence s;

                size_t_type n;
                for ( n=0; n<iter2->second.size(); n++ )
                {
                    s.push_back(iter2->second[n]);
                }

                item.value(s);

                attrib.AttributeFloat().push_back(item);
            }

            // attribute3
            const AttributeTypeStringAttributeTypeComplexFloat& stroage3 = attributeComplexFloat_.get_attrib_storage();
            typename AttributeTypeStringAttributeTypeComplexFloat::const_iterator iter3;
            for ( iter3=stroage3.begin(); iter3!=stroage3.end(); iter3++ )
            {
                gadgetron::AttributeComplex item(iter3->first);

                gadgetron::AttributeComplex::value_real_sequence s;
                gadgetron::AttributeComplex::value_imag_sequence t;

                size_t_type n;
                for ( n=0; n<iter3->second.size(); n++ )
                {
                    s.push_back(iter3->second[n].real());
                    t.push_back(iter3->second[n].imag());
                }

                item.value_real(s);
                item.value_imag(t);

                attrib.AttributeComplex().push_back(item);
            }

            // attribute4
            const AttributeTypeStringAttributeTypeString& stroage4 = attributeString_.get_attrib_storage();
            typename AttributeTypeStringAttributeTypeString::const_iterator iter4;
            for ( iter4=stroage4.begin(); iter4!=stroage4.end(); iter4++ )
            {
                gadgetron::AttributeString item(iter4->first);

                gadgetron::AttributeString::value_sequence s;

                size_t_type n;
                for ( n=0; n<iter4->second.size(); n++ )
                {
                    s.push_back(iter4->second[n]);
                }

                item.value(s);

                attrib.AttributeString().push_back(item);
            }

            // get the xml
            xml_schema::namespace_infomap map;
            map[""].name = "http://gadgetron.sf.net/gadgetron";
            map[""].schema = "gadgetronMetaAttributes.xsd";
            std::stringstream str;
            gadgetron::gadgetronMetaAttributes_(str, attrib, map);
            std::string xmltext = str.str();

            // fill the buffer
            // the first size_t_type is the length of xml document
            size_t_type xmlLen = xmltext.length();
            len = xmlLen + 1 + sizeof(size_t_type);
            buf = new char[len];
            memset(buf, '\0', len);
            memcpy(buf, &len, sizeof(size_t_type));
            memcpy(buf+sizeof(size_t_type), xmltext.c_str(), xmlLen);

            // mutex_.unlock();
        }
        catch(...)
        {
            // mutex_.unlock();
            GADGET_ERROR_MSG("Errors happened in hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::serialize(char*& buf, size_t_type& len) ... ");
            return false;
        }

        return true;
    }

    template <typename TInteger, typename TFloat, typename TComplexFloat>
    bool hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::deserialize(char* buf, size_t_type& len)
    {
        try
        {
            size_t_type xmlLen;
            memcpy(&xmlLen, buf, sizeof(size_t_type) );
            xmlLen -= sizeof(size_t_type);
            GADGET_CHECK_RETURN_FALSE(this->deserializeContent( buf+sizeof(size_t_type), xmlLen ) );
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::deserialize(char* buf, size_t_type& len) ... ");
            return false;
        }

        return true;
    }

    template <typename TInteger, typename TFloat, typename TComplexFloat>
    bool hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::deserializeContent(char* buf, size_t_type& xmlLen)
    {
        // ACE_GUARD_RETURN(ACE_Thread_Mutex, guard, mtx_, false);

        try
        {
            // clear the attributes
            GADGET_CHECK_RETURN_FALSE(this->clear());

            char * gadgetron_home = std::getenv("GADGETRON_HOME");
            std::string schema_file_name = std::string(gadgetron_home);
            schema_file_name.append("/schema/gadgetronMetaAttributes.xsd");

            std::string tmp(schema_file_name);
            tmp = url_encode(tmp);
            schema_file_name = tmp;

            xml_schema::properties props;
            props.schema_location (
                "http://gadgetron.sf.net/gadgetron",
                std::string (schema_file_name));

            std::string xml(buf);
            std::istringstream str_stream(xml, std::stringstream::in);

            boost::shared_ptr<gadgetron::gadgetronMetaAttributes> cfg;

            try
            {
                cfg = boost::shared_ptr<gadgetron::gadgetronMetaAttributes>(gadgetron::gadgetronMetaAttributes_ (str_stream,0,props));
            }
            catch (const xml_schema::exception& e)
            {
                GADGET_ERROR_MSG("Failed to parse gadgetron meta attributes Parameters: " << e.what());
            }

            gadgetron::gadgetronMetaAttributes::AttributeInteger_sequence& integer = cfg->AttributeInteger();
            gadgetron::gadgetronMetaAttributes::AttributeInteger_sequence::const_iterator iter = integer.begin();
            for ( ; iter!=integer.end(); iter++ )
            {
                std::string name = iter->name();
                std::vector<TInteger> item;

                const gadgetron::AttributeInteger::value_sequence& s = iter->value();
                gadgetron::AttributeInteger::value_sequence::const_iterator iter_s = s.begin();

                for ( ; iter_s!= s.end(); iter_s++ )
                {
                    item.push_back( (TInteger)(*iter_s) );
                }

                attributeInteger_.set(name, item);
            }

            gadgetron::gadgetronMetaAttributes::AttributeFloat_sequence& floatSeq = cfg->AttributeFloat();
            gadgetron::gadgetronMetaAttributes::AttributeFloat_sequence::const_iterator iterFloat = floatSeq.begin();
            for ( ; iterFloat!=floatSeq.end(); iterFloat++ )
            {
                std::string name = iterFloat->name();
                std::vector<TFloat> item;

                const gadgetron::AttributeFloat::value_sequence& s = iterFloat->value();
                gadgetron::AttributeFloat::value_sequence::const_iterator iter_s = s.begin();

                for ( ; iter_s!= s.end(); iter_s++ )
                {
                    item.push_back( (TFloat)(*iter_s) );
                }

                attributeFloat_.set(name, item);
            }

            gadgetron::gadgetronMetaAttributes::AttributeComplex_sequence& complexSeq = cfg->AttributeComplex();
            gadgetron::gadgetronMetaAttributes::AttributeComplex_sequence::const_iterator iterComplex = complexSeq.begin();
            for ( ; iterComplex!=complexSeq.end(); iterComplex++ )
            {
                std::string name = iterComplex->name();
                std::vector<TComplexFloat> item;

                const gadgetron::AttributeComplex::value_real_sequence& r = iterComplex->value_real();
                gadgetron::AttributeComplex::value_real_sequence::const_iterator iter_r = r.begin();

                const gadgetron::AttributeComplex::value_imag_sequence& i = iterComplex->value_imag();
                gadgetron::AttributeComplex::value_imag_sequence::const_iterator iter_i = i.begin();

                TComplexFloat v;
                cx_value_type rv, iv;
                for ( ; (iter_r!=r.end() && iter_i!=i.end()); iter_r++,iter_i++ )
                {
                    rv = (cx_value_type)(*iter_r);
                    iv = (cx_value_type)(*iter_i);

                    v = TComplexFloat(rv, iv);

                    item.push_back( v );
                }

                attributeComplexFloat_.set(name, item);
            }

            gadgetron::gadgetronMetaAttributes::AttributeString_sequence& StringSeq = cfg->AttributeString();
            gadgetron::gadgetronMetaAttributes::AttributeString_sequence::const_iterator iterString = StringSeq.begin();
            for ( ; iterString!=StringSeq.end(); iterString++ )
            {
                std::string name = iterString->name();
                std::vector<std::string> item;

                const gadgetron::AttributeString::value_sequence& s = iterString->value();
                gadgetron::AttributeString::value_sequence::const_iterator iter_s = s.begin();

                for ( ; iter_s!= s.end(); iter_s++ )
                {
                    item.push_back( *iter_s );
                }

                attributeString_.set(name, item);
            }

            //// parse the xml
            //TiXmlDocument doc;
            //doc.Parse(buf);

            //TiXmlHandle docHandle(&doc);

            //TiXmlElement* pElem=docHandle.FirstChildElement().Element()->FirstChild()->ToElement();
            //for ( pElem; pElem!=NULL; pElem=pElem->NextSiblingElement() )
            //{
            //    std::string elemStr = std::string(pElem->Value());

            //    TiXmlElement* nameElem = pElem->FirstChildElement("name")->ToElement();
            //    GADGET_CHECK_RETURN_FALSE(nameElem!=NULL);
            //    std::string name = nameElem->FirstChild()->ToText()->ValueStr();

            //    size_t_type n;

            //    if ( elemStr == "AttributeInteger" )
            //    {
            //        std::vector<TInteger> item;
            //        TiXmlElement* value = pElem->FirstChildElement("value")->ToElement();

            //        if ( value != NULL )
            //        {
            //            do
            //            {
            //                std::string valueStr = value->FirstChild()->ToText()->ValueStr();
            //                item.push_back( (TInteger)std::atoi( valueStr.c_str() ) );
            //                value = value->NextSiblingElement( "value" );
            //            }
            //            while ( value != NULL );
            //        }

            //        attributeInteger_.set(name, item);
            //    }
            //    else if ( elemStr == "AttributeFloat" )
            //    {
            //        std::vector<TFloat> item;
            //        TiXmlElement* value = pElem->FirstChildElement("value")->ToElement();

            //        if ( value != NULL )
            //        {
            //            do
            //            {
            //                std::string valueStr = value->FirstChild()->ToText()->ValueStr();
            //                item.push_back( (TFloat)std::atof( valueStr.c_str() ) );
            //                value = value->NextSiblingElement( "value" );
            //            }
            //            while ( value != NULL );
            //        }

            //        attributeFloat_.set(name, item);
            //    }
            //    else if ( elemStr == "AttributeComplex" )
            //    {
            //        std::vector<TComplexFloat> item;
            //        TiXmlElement* value_real = pElem->FirstChildElement("value_real")->ToElement();
            //        TiXmlElement* value_imag = pElem->FirstChildElement("value_imag")->ToElement();

            //        cx_value_type r_v, i_v;

            //        if ( value_real!=NULL && value_imag!=NULL )
            //        {
            //            do
            //            {
            //                std::string valueStr = value_real->FirstChild()->ToText()->ValueStr();
            //                r_v = (cx_value_type)std::atof( valueStr.c_str() );

            //                valueStr = value_imag->FirstChild()->ToText()->ValueStr();
            //                i_v = (cx_value_type)std::atof( valueStr.c_str() );

            //                item.push_back( TComplexFloat(r_v, i_v) );

            //                value_real = value_real->NextSiblingElement( "value_real" );
            //                value_imag = value_imag->NextSiblingElement( "value_imag" );
            //            }
            //            while ( value_real!=NULL && value_imag!=NULL );
            //        }

            //        attributeComplexFloat_.set(name, item);
            //    }
            //    else if ( elemStr == "AttributeString" )
            //    {
            //        std::vector< std::string > item;
            //        TiXmlElement* value = pElem->FirstChildElement("value")->ToElement();

            //        if ( value != NULL )
            //        {
            //            std::string valueStr;

            //            do
            //            {
            //                valueStr = value->FirstChild()->ToText()->ValueStr();
            //                item.push_back( valueStr );
            //                value = value->NextSiblingElement( "value" );
            //            }
            //            while ( value != NULL );
            //        }

            //        attributeString_.set(name, item);
            //    }
            //}
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::deserializeContent(char* buf, size_t_type& xmlLen) ... ");
            return false;
        }

        return true;
    }

    template <typename TInteger, typename TFloat, typename TComplexFloat>
    void hoNDMetaAttributes<TInteger, TFloat, TComplexFloat>::print(std::ostream& os) const
    {
        using namespace std;
        os << "-------------- Gagdgetron attributes -------------" << endl;
        os << "Attribute 1 : " << std::endl;
        attributeInteger_.printContent(os);
        os << "--------------------------------------------------" << endl;
        os << "Attribute 2 : " << std::endl;
        attributeFloat_.printContent(os);
        os << "--------------------------------------------------" << endl;
        os << "Attribute 3 : " << std::endl;
        attributeComplexFloat_.printContent(os);
        os << "--------------------------------------------------" << endl;
        os << "Attribute 4 : " << std::endl;
        attributeString_.printContent(os);
        os << "--------------------------------------------------" << endl;
    }

    template EXPORTCPUCORE class hoNDMetaAttributes<int, float, std::complex<float> >;
    template EXPORTCPUCORE class hoNDMetaAttributes<int, double, std::complex<double> >;
    template EXPORTCPUCORE class hoNDMetaAttributes<int, float, std::complex<double> >;
    template EXPORTCPUCORE class hoNDMetaAttributes<int, double, std::complex<float> >;
    template EXPORTCPUCORE class hoNDMetaAttributes<long long, float, std::complex<float> >;
    template EXPORTCPUCORE class hoNDMetaAttributes<long long, double, std::complex<double> >;
    template EXPORTCPUCORE class hoNDMetaAttributes<long long, float, std::complex<double> >;
    template EXPORTCPUCORE class hoNDMetaAttributes<long long, double, std::complex<float> >;
}
