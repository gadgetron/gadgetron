/** \file       hoNDMetaAttributesBase.h
    \brief      meta attributes class for the gadgetron

                The field key is a string. The key-value pairs are stored in the std::map.
                The sereialization and deserialization of meta info is using the xml format.

    \author     Hui Xue
*/

#pragma once

#include "GadgetronException.h"
#include "GadgetronCommon.h"

#include <map>
#include <vector>
#include <iostream>
#include <strstream>

namespace Gadgetron
{
    template <typename T>
    class hoNDMetaAttributesBase
    {
    public:

        typedef hoNDMetaAttributesBase<T> Self;

        typedef T element_type;
        typedef T value_type;
        typedef float coord_type;

        typedef unsigned long long size_t_type;

        typedef std::map< std::string, std::vector<T> > AttributeStorageType;

        hoNDMetaAttributesBase();
        virtual ~hoNDMetaAttributesBase();

        hoNDMetaAttributesBase(const Self& attrib);

        Self& operator=(const Self& attrib);

        /// clear the attributes
        bool clear();

        /// swap the attributes
        bool swap(Self& attrib);

        /// whether an field exists
        bool exist(const std::string& key);
        /// whether the ind element of an field exists
        bool exist(const std::string& key, size_t_type ind);

        /// get an element or an item
        bool get(const std::string& key, size_t_type ind, T& v) const;
        bool get(const std::string& key, std::vector<T>& elem) const;

        /// set an element or an item
        /// if the item or element exists, replace it
        /// if the item does not exist, add this item to the attributes
        bool set(const std::string& key, size_t_type ind, const T& v);
        bool set(const std::string& key, const std::vector<T>& elem);

        /// if the item exists, add this element to the end of vector
        /// if not, create a new item
        bool set(const std::string& key, const T& v);

        /// erase an element or an item
        bool erase(const std::string& key, size_t_type ind, T& v);
        bool erase(const std::string& key, size_t_type ind);
        bool erase(const std::string& key, std::vector<T>& elem);
        bool erase(const std::string& key);

        /// get the attribute storage
        AttributeStorageType& get_attrib_storage() { return meta_attribute_; }
        const AttributeStorageType& get_attrib_storage() const { return meta_attribute_; }

        /// print out the image information
        virtual void print(std::ostream& os) const;
        virtual void printContent(std::ostream& os) const;

    protected:

        std::map< std::string, std::vector<T> > meta_attribute_;
    };

    template <typename T>
    hoNDMetaAttributesBase<T>::hoNDMetaAttributesBase()
    {
    }

    template <typename T>
    hoNDMetaAttributesBase<T>::~hoNDMetaAttributesBase()
    {
        this->clear();
    }

    template <typename T>
    hoNDMetaAttributesBase<T>::hoNDMetaAttributesBase(const Self& attrib)
    {
        if ( this != &attrib )
        {
            meta_attribute_ = attrib.meta_attribute_;
        }
    }

    template <typename T>
    inline hoNDMetaAttributesBase<T>& hoNDMetaAttributesBase<T>::operator=(const Self& attrib)
    {
        if ( this == &attrib ) return *this;
        meta_attribute_ = attrib.meta_attribute_;
        return *this;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::clear()
    {
        try
        {
            meta_attribute_.clear();
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in hoNDMetaAttributesBase<T>::clear() ... ");
            return false;
        }

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::swap(Self& attrib)
    {
        try
        {
            meta_attribute_.swap(attrib.meta_attribute_);
        }
        catch(...)
        {
            GADGET_ERROR_MSG("Errors happened in hoNDMetaAttributesBase<T>::swap(Self& attrib) ... ");
            return false;
        }

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::exist(const std::string& key)
    {
        typename AttributeStorageType::const_iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return false;
        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::exist(const std::string& key, size_t_type ind)
    {
        typename AttributeStorageType::const_iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return false;
        if ( ind >= iter->second.size() ) return false;

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::get(const std::string& key, size_t_type ind, T& v) const
    {
        typename AttributeStorageType::const_iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return false;

        if ( ind >= iter->second.size() ) return false;

        v = iter->second[ind];

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::get(const std::string& key, std::vector<T>& elem) const
    {
        typename AttributeStorageType::const_iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return false;
        elem = iter->second;

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::set(const std::string& key, size_t_type ind, const T& v)
    {
        typename AttributeStorageType::iterator iter = meta_attribute_.find(key);

        if( iter == meta_attribute_.end() )
        {
            // insert a new item
            std::vector<T> item(1);
            item[0] = v;

            meta_attribute_[key] = item;
        }
        else
        {
            if ( ind > iter->second.size() ) return false;

            if ( ind == iter->second.size() )
            {
                iter->second.push_back(v);
            }
            else
            {
                iter->second[ind] = v;
            }
        }

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::set(const std::string& key, const std::vector<T>& item)
    {
        meta_attribute_[key] = item;
        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::set(const std::string& key, const T& v)
    {
        typename AttributeStorageType::iterator iter = meta_attribute_.find(key);

        if( iter == meta_attribute_.end() )
        {
            // insert a new item
            std::vector<T> item(1);
            item[0] = v;

            meta_attribute_[key] = item;
        }
        else
        {
            iter->second.push_back(v);
        }

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::erase(const std::string& key, size_t_type ind, T& v)
    {
        typename AttributeStorageType::iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return true; // no need to erase anything

        size_t_type nElem = iter->second.size();
        if ( ind >= nElem ) return true;

        v = iter->second[ind];

        if ( nElem == 1 )
        {
            meta_attribute_.erase(iter);
        }
        else
        {
            iter->second.erase(iter->second.begin()+ind);
        }

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::erase(const std::string& key, size_t_type ind)
    {
        typename AttributeStorageType::iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return true; // no need to erase anything

        size_t_type nElem = iter->second.size();
        if ( ind >= nElem ) return true;

        if ( nElem == 1 )
        {
            meta_attribute_.erase(iter);
        }
        else
        {
            iter->second.erase(iter->second.begin()+ind);
        }

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::erase(const std::string& key, std::vector<T>& elem)
    {
        typename AttributeStorageType::iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return true;

        elem = iter->second;
        meta_attribute_.erase(iter);

        return true;
    }

    template <typename T>
    inline bool hoNDMetaAttributesBase<T>::erase(const std::string& key)
    {
        typename AttributeStorageType::iterator iter = meta_attribute_.find(key);
        if( iter == meta_attribute_.end() ) return true;

        meta_attribute_.erase(iter);

        return true;
    }

    template <typename T>
    inline void hoNDMetaAttributesBase<T>::printContent(std::ostream& os) const
    {
        using namespace std;
        os << "Attributes data type is : " << std::string(typeid(T).name()) << std::endl;
        os << "Number of items in attributes is : " << meta_attribute_.size() << std::endl;

        typename AttributeStorageType::const_iterator iter;
        for ( iter=meta_attribute_.begin(); iter!=meta_attribute_.end(); iter++ )
        {
            os << " [" << iter->first << "] : \t";

            size_t_type nElem = iter->second.size();

            unsigned int n;
            for ( n=0; n<GT_MIN(nElem, 20); n++ )
            {
                os << iter->second[n] << " ";
            }
            if ( n != nElem ) os << " ... ";
            os << endl;
        }
    }

    template <typename T>
    inline void hoNDMetaAttributesBase<T>::print(std::ostream& os) const
    {
        using namespace std;
        os << "-------------- Gagdgetron attributes -------------" << endl;
        this->printContent(os);
    }
}
