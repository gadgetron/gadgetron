#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho3DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho3DArray();
    ho3DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz);
    ho3DArray(std::vector<unsigned long long> *dimensions);
    ho3DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho3DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct = false);
    ho3DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    ho3DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho3DArray();

    ho3DArray(const ho3DArray<T>& a);
    ho3DArray<T>& operator=(const ho3DArray<T>& rhs);

    virtual void create(std::vector<unsigned long long>& dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz);
    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, T* data, bool delete_data_on_destruct = false);

    T& operator()(unsigned long long x, unsigned long long y, unsigned long long z);
    const T& operator()(unsigned long long x, unsigned long long y, unsigned long long z) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T*** accesser_;
};

}

#include <ho3DArray.hxx>
