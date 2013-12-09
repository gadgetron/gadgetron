#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho2DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho2DArray();
    ho2DArray(unsigned long long sx, unsigned long long sy);
    ho2DArray(std::vector<unsigned long long> *dimensions);
    ho2DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho2DArray(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct = false);
    ho2DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    ho2DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho2DArray();

    ho2DArray(const ho2DArray<T>& a);
    ho2DArray<T>& operator=(const ho2DArray<T>& rhs);

    virtual void create(std::vector<unsigned long long>& dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(unsigned long long sx, unsigned long long sy);
    virtual bool createArray(unsigned long long sx, unsigned long long sy, T* data, bool delete_data_on_destruct = false);

    T& operator()(unsigned long long x , unsigned long long y);
    const T& operator()(unsigned long long x , unsigned long long y) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T** accesser_;
};

}

#include <ho2DArray.hxx>
