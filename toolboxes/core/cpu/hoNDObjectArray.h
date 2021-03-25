/** \file   hoNDObjectArray.h
\brief  CPU-based N-dimensional array for object pointers
        if delete_data_on_destruct == true, the object will be released; otherwise, only the object array memory is released
        the stored object should implement the serialize/deserialize functions
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

        virtual bool serialize(char*& buf, size_t& len) const;
        virtual bool deserialize(char* buf, size_t& len);
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
    bool hoNDObjectArray<TObjectType>::serialize(char*& buf, size_t& len) const
    {
        if (buf != NULL) delete[] buf;

        size_t NDim = this->dimensions_.size();

        // number of dimensions + dimension vector + contents
        len = sizeof(size_t) + sizeof(size_t) * NDim;

        // find len for every element
        size_t N = this->get_number_of_elements();

        std::vector<char*> buf_obj;
        std::vector<size_t> buf_len;
        size_t ii;

        if (N > 0)
        {
            buf_obj.resize(N, NULL);
            buf_len.resize(N, 0);

            for (ii = 0; ii < N; ii++)
            {
                this->data_[ii].serialize(buf_obj[ii], buf_len[ii]);

                len += sizeof(size_t); // buf_len[ii]
                len += buf_len[ii]; // buf[ii]
            }
        }

        buf = new char[len];

        size_t ind = 0;
        memcpy(buf, &NDim, sizeof(size_t));
        ind += sizeof(size_t);

        if (NDim > 0)
        {
            memcpy(buf + ind, &(this->dimensions_[0]), sizeof(size_t) * NDim);
            ind += sizeof(size_t) * NDim;

            if(N>0)
            {
                for (ii = 0; ii < N; ii++)
                {
                    memcpy(buf + ind, &(buf_len[ii]), sizeof(size_t)); ind += sizeof(size_t);
                    memcpy(buf + ind, buf_obj[ii], buf_len[ii]); ind += buf_len[ii];
                }

                for (ii = 0; ii < N; ii++)
                {
                    delete[] buf_obj[ii];
                    buf_obj[ii] = NULL;
                }
            }
        }

        return true;
    }

    template <typename TObjectType>
    bool hoNDObjectArray<TObjectType>::deserialize(char* buf, size_t& len)
    {
        size_t ind = 0;

        size_t NDim;
        memcpy(&NDim, buf, sizeof(size_t));
        ind += sizeof(size_t);

        if (NDim > 0)
        {
            std::vector<size_t> dimensions(NDim);
            memcpy(&dimensions[0], buf + ind, sizeof(size_t) * NDim);
            ind += sizeof(size_t) * NDim;

            // allocate memory
            this->create(&dimensions);

            // deserialize the content
            size_t N = this->get_number_of_elements();
            size_t ii;
            for (ii = 0; ii < N; ii++)
            {
                TObjectType a;

                size_t len_obj(0);
                memcpy(&len_obj, buf + ind, sizeof(size_t));
                ind += sizeof(size_t);

                a.deserialize(buf + ind, len_obj);
                ind += len_obj;

                this->data_[ii] = a;
            }
        }
        else
        {
            this->clear();
        }

        return true;
    }
}