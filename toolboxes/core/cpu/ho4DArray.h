#pragma once

#include "hoNDArray.h"

namespace Gadgetron{

template <class T> class ho4DArray : public hoNDArray<T>
{
public:

    typedef hoNDArray<T> BaseClass;

    ho4DArray();
    ho4DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss);
    ho4DArray(std::vector<unsigned long long> *dimensions);
    ho4DArray(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);
    ho4DArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, T* data, bool delete_data_on_destruct = false);
    ho4DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions);
    ho4DArray(boost::shared_ptr< std::vector<unsigned long long> > dimensions, T* data, bool delete_data_on_destruct = false);

    virtual ~ho4DArray();

    ho4DArray(const ho4DArray<T>& a);
    ho4DArray<T>& operator=(const ho4DArray<T>& rhs);

    virtual void create(std::vector<unsigned long long>& dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions);
    virtual void create(std::vector<unsigned long long> *dimensions, T* data, bool delete_data_on_destruct = false);

    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss);
    virtual bool createArray(unsigned long long sx, unsigned long long sy, unsigned long long sz, unsigned long long ss, T* data, bool delete_data_on_destruct = false);

    T& operator()(unsigned long long x, unsigned long long y, unsigned long long z, unsigned long long s);
    const T& operator()(unsigned long long x , unsigned long long y, unsigned long long z, unsigned long long s) const;

    virtual void print(std::ostream& os) const;

protected:

    using BaseClass::dimensions_;
    using BaseClass::offsetFactors_;
    using BaseClass::data_;
    using BaseClass::elements_;
    using BaseClass::delete_data_on_destruct_;

    bool init_accesser();
    bool release_accesser();

    T**** accesser_;
};

}

#include <ho4DArray.hxx>
