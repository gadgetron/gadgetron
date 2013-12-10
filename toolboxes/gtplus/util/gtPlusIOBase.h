/** \file       gtPlusIOBase.h
    \brief      Define the base IO funcatinality for GtPlus toolbox
    \author     Hui Xue
*/

#pragma once

#include <iostream>
#include <typeinfo>

#include "GtPlusExport.h"
#include "NDArray.h"
#include "complext.h"
#include "vector_td.h"
#include "GadgetronException.h"
#include "GadgetronCommon.h"

#include <mkl.h>

#include "hoNDArray.h"
#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"

#include "hoNDArray_fileio.h"
#include "hoNDArray_elemwise.h"

// the file input/output utility functions

#ifdef GT_Complex8
    #undef GT_Complex8
#endif // GT_Complex8
typedef std::complex<float> GT_Complex8;

#ifdef GT_Complex16
    #undef GT_Complex16
#endif // GT_Complex16
typedef std::complex<double> GT_Complex16;

namespace Gadgetron { namespace gtPlus {

class EXPORTGTPLUS gtPlusIOWorker
{
public:

    gtPlusIOWorker(const std::string& ioTag, bool readFlag=true);
    virtual ~gtPlusIOWorker();

    // open the file stream
    // readFlag: true, read mode; false, write mode
    virtual bool open();

    // close the file stream
    virtual bool close();

    // the current file offset
    long tell();

    // set the file offset
    bool seek(long long offset);

    // reset the file to the beginning
    bool reset() { return (this->seek(0)); }

    // check the status of i/o operations
    bool IOinError();

    // read/write
    // len: number of bytes in data
    bool read(char* data, long long len);
    bool write(const char* data, long long len);

protected:

    std::string ioTag_;
    std::fstream fid_;
    bool readFlag_;
};

class EXPORTGTPLUS gtPlusIOBase
{
public:

    gtPlusIOBase() {}
    virtual ~gtPlusIOBase() {}

public:

    void printInfo(std::ostream& os);

    // buffer read/write functions
    // length: number of bytes
    bool readFromFile(const std::string& filename, char*& data, long long& length);
    bool writeToFile(const std::string& filename, char* data, long long length);

    // general export/input for ND array
    //template <typename T> bool exportNDArray(const hoNDArray<T>& a, const std::string& filename) const;
    //template <typename T> bool importNDArray(hoNDArray<T>& a, std::string& filename) const;
};

/*template <typename T>
bool gtPlusIOBase::exportNDArray(const hoNDArray<T>& a, const std::string& filename) const
{
    GADGET_CHECK_RETURN_FALSE( Gadgetron::write_nd_array(const_cast<hoNDArray<T>* >(&a), filename) == 0 );
    return true;
}

template <typename T> 
bool gtPlusIOBase::importNDArray(hoNDArray<T>& a, std::string& filename) const
{
    try
    {
        boost::shared_ptr< hoNDArray<T> > aRead;
        aRead = Gadgetron::read_nd_array(filename);
        a = *aRead;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOBase::importNDArray(hoNDArray<T>& a, const std::string& filename) ... ");
        return false;
    }

    return true;
}*/

}}
