/** \file   hoNDObjectArray.h
\brief  N-dimensional array for objects
        The stored objects should support read/write interfaces
\author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{
    template <typename TObjectType> class hoNDObjectArray : public hoNDArray<TObjectType>
    {
    public:

        typedef hoNDArray<TObjectType> BaseClass;

        hoNDObjectArray();
        virtual ~hoNDObjectArray();

        hoNDObjectArray(const std::vector<size_t> *dimensions);
        hoNDObjectArray(const std::vector<size_t> &dimensions);
        hoNDObjectArray(boost::shared_ptr< std::vector<size_t> > dimensions);
        hoNDObjectArray(const hoNDObjectArray<TObjectType>* a);
        hoNDObjectArray(const hoNDObjectArray<TObjectType> &a);
    };

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray() : BaseClass()
    {
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const std::vector<size_t> *dimensions) : BaseClass(dimensions)
    {
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const std::vector<size_t> &dimensions) : BaseClass(dimensions)
    {
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(boost::shared_ptr< std::vector<size_t> > dimensions) : BaseClass(dimensions)
    {
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::~hoNDObjectArray()
    {
        this->clear();
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const hoNDObjectArray<TObjectType>* a)
    {
        this->copyFrom(*a);
    }

    template <typename TObjectType>
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const hoNDObjectArray<TObjectType> &a)
    {
        this->copyFrom(a);
    }
}