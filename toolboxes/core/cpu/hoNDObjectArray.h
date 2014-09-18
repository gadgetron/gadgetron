/** \file   hoNDObjectArray.h
\brief  CPU-based N-dimensional array for object pointers
if delete_data_on_destruct == true, the object will be released; otherwise, only the object array memory is released
\author Hui Xue
*/

#pragma once

#include "hoNDArray.h"

namespace Gadgetron
{

    template <typename TObjectType> class hoNDObjectArray : public hoNDArray<TObjectType*>
    {
    public:

        typedef hoNDArray<TObjectType*> BaseClass;
        typedef float coord_type;
        typedef typename BaseClass::value_type value_type;

        hoNDObjectArray();

        explicit hoNDObjectArray(std::vector<size_t> &dimensions);
        explicit hoNDObjectArray(std::vector<size_t> *dimensions);
        explicit hoNDObjectArray(boost::shared_ptr< std::vector<size_t> > dimensions);

        virtual ~hoNDObjectArray();

        // Copy constructors
        hoNDObjectArray(const hoNDObjectArray<TObjectType> &a);
        explicit hoNDObjectArray(const hoNDObjectArray<TObjectType> *a);

        // Assignment operator
        hoNDObjectArray& operator=(const hoNDObjectArray& rhs);

        virtual void create(std::vector<size_t>& dimensions);
        virtual void create(std::vector<size_t> *dimensions);
        virtual void create(boost::shared_ptr< std::vector<size_t> > dimensions);

        void get_sub_array(const std::vector<size_t>& start, std::vector<size_t>& size, hoNDObjectArray<TObjectType>& out);

        virtual void print(std::ostream& os) const;

    protected:

        using BaseClass::dimensions_;
        using BaseClass::offsetFactors_;
        using BaseClass::data_;
        using BaseClass::elements_;
        using BaseClass::delete_data_on_destruct_;
    };

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>::hoNDObjectArray() : BaseClass() 
    {
    }

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>::hoNDObjectArray(std::vector<size_t> *dimensions) : BaseClass(dimensions)
    {
        this->create(dimensions);
    }

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>::hoNDObjectArray(std::vector<size_t> &dimensions) : BaseClass(dimensions)
    {
        this->create(dimensions);
    }

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>::hoNDObjectArray(boost::shared_ptr< std::vector<size_t> > dimensions) : BaseClass(dimensions)
    {
        this->create(dimensions);
    }

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>::~hoNDObjectArray()
    {
        if (this->delete_data_on_destruct_)
        {
            size_t n;
            for ( n=0; n<this->elements_; n++ )
            {
                if ( this->data_[n] != NULL )
                {
                    delete this->data_[n];
                    this->data_[n] = NULL;
                }
            }

            this->deallocate_memory();
        }
    }

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const hoNDObjectArray<TObjectType>  *a) : BaseClass(a)
    {
        this->delete_data_on_destruct_ = false;
    }

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>::hoNDObjectArray(const hoNDObjectArray<TObjectType> &a) : BaseClass(a)
    {
        this->delete_data_on_destruct_ = false;
    }

    template <typename TObjectType> 
    hoNDObjectArray<TObjectType>& hoNDObjectArray<TObjectType>::operator=(const hoNDObjectArray<TObjectType>& rhs)
    {
        if ( &rhs == this ) return *this;

        BaseClass::operator=(rhs);

        this->delete_data_on_destruct_ = false;

        return *this;
    }

    template <typename TObjectType> 
    void hoNDObjectArray<TObjectType>::create(std::vector<size_t>& dimensions)
    {
        BaseClass::create(dimensions);

        for ( size_t n=0; n<this->elements_; n++ )
        {
            this->data_[n] = NULL;
        }
    }

    template <typename TObjectType> 
    void hoNDObjectArray<TObjectType>::create(std::vector<size_t> *dimensions)
    {
        BaseClass::create(dimensions);

        for ( size_t n=0; n<this->elements_; n++ )
        {
            this->data_[n] = NULL;
        }
    }

    template <typename TObjectType> 
    void hoNDObjectArray<TObjectType>::create(boost::shared_ptr< std::vector<size_t> > dimensions)
    {
        BaseClass::create(dimensions);

        for ( size_t n=0; n<this->elements_; n++ )
        {
            this->data_[n] = NULL;
        }
    }

    template <typename TObjectType> 
    void hoNDObjectArray<TObjectType>::get_sub_array(const std::vector<size_t>& start, std::vector<size_t>& size, hoNDObjectArray<TObjectType>& out)
    {
        if ( start.size() != size.size() )
        {
            BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
        }

        if ( start.size() != (*dimensions_).size() )
        {
            BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
        }

        out.create(&size);

        if ( out.get_number_of_elements() == this->get_number_of_elements() )
        {
            out = *this;
            return;
        }

        std::vector<size_t> end(start.size());

        size_t ii;
        for ( ii=0; ii<start.size(); ii++ )
        {
            end[ii] = start[ii] + size[ii] - 1;
            if ( end[ii] >= (*dimensions_)[ii] )
            {
                BOOST_THROW_EXCEPTION( runtime_error("hoNDArray<>::get_sub_array failed"));
            }
        }

        out.delete_data_on_destruct(false);
    }

    template <typename TObjectType> 
    void hoNDObjectArray<TObjectType>::print(std::ostream& os) const
    {
        using namespace std;

        os.unsetf(std::ios::scientific);
        os.setf(ios::fixed);

        os << "-------------- Gagdgetron ND Object Array -------------" << endl;
        this->printContent(os);
    }
}
