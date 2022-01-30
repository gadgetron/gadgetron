/** \file   hoNDObjectArray.h
\brief  store the objects as a ND array
        the stored object should support read/write interface
\author Hui Xue
*/

#pragma once

#include "hoNDArray.h"
#include "io/primitives.h"

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

        void write(std::ostream &stream) const;
        void read(std::istream& stream);
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

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::write(std::ostream &stream) const
    {
        Gadgetron::Core::IO::write(stream, this->dimensions_);
        size_t N = this->get_number_of_elements();
        if(N>0)
        {
            for (auto ii = 0; ii < N; ii++)
            {
                Gadgetron::Core::IO::write(stream, this->data_[ii]);
            }
        }
    }

    template <typename TObjectType>
    void hoNDObjectArray<TObjectType>::read(std::istream& stream)
    {
        std::vector<size_t> dimensions;
        Gadgetron::Core::IO::read(stream, dimensions);

        this->create(dimensions);

        size_t N = this->get_number_of_elements();
        for (auto ii = 0; ii < N; ii++)
        {
            Gadgetron::Core::IO::read(stream, this->data_[ii]);
        }
    }
}
