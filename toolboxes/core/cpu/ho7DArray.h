#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho7DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho7DArray();
    ho7DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, unsigned long long sa);
    ho7DArray(std::vector<unsigned long long> *dimensions);
    ho7DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho7DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, unsigned long long sa, T* data, bool delete_data_on_destruct = false);
    ho7DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    ho7DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho7DArray();

    ho7DArray(const ho7DArray<T>& a);
    ho7DArray<T>& operator=(const ho7DArray<T>& rhs);

    virtual void create(std::vector<unsigned long long>& dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, unsigned long long sa);
    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, unsigned long long sa, T* data, bool delete_data_on_destruct = false);

    T& operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p, unsigned long long r, unsigned long long a);
    const T& operator()(unsigned long long x , unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p, unsigned long long r, unsigned long long a) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T******* accesser_;
};

}

#include <ho7DArray.hxx>
