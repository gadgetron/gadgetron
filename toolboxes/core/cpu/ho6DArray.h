#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho6DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho6DArray();
    ho6DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr);
    ho6DArray(std::vector<unsigned long long> *dimensions);
    ho6DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho6DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, T* data, bool delete_data_on_destruct = false);
    ho6DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    ho6DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho6DArray();

    ho6DArray(const ho6DArray<T>& a);
    ho6DArray<T>& operator=(const ho6DArray<T>& rhs);

    virtual void create(std::vector<unsigned long long>& dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr);
    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, unsigned long long sp, unsigned long long sr, T* data, bool delete_data_on_destruct = false);

    T& operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p, unsigned long long r);
    const T& operator()(unsigned long long x , unsigned long long y, unsigned long long z, unsigned long long s, unsigned long long p, unsigned long long r) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T****** accesser_;
};

}

#include <ho6DArray.hxx>
