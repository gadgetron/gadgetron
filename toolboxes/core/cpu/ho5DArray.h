#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho5DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho5DArray();
    ho5DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp);
    ho5DArray(std::vector<unsigned long long> *dimensions);
    ho5DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho5DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, T* data, bool delete_data_on_destruct = false);
    ho5DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    ho5DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho5DArray();

    ho5DArray(const ho5DArray<T>& a);
    ho5DArray<T>& operator=(const ho5DArray<T>& rhs);

    virtual void create(std::vector<unsigned long long>& dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp);
    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, T* data, bool delete_data_on_destruct = false);

    T& operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p);
    const T& operator()(unsigned long long x , unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T***** accesser_;
};

}

#include <ho5DArray.hxx>
