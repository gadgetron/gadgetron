/** \file       gtPlusIOBase.cpp
    \brief      Define the base IO funcatinality for GtPlus toolbox
    \author     Hui Xue
*/

#include <gtPlusIOBase.h>

namespace Gadgetron { namespace gtPlus {

gtPlusIOWorker::gtPlusIOWorker(const std::string& ioTag, bool readFlag) : ioTag_(ioTag), readFlag_(readFlag)
{
}

gtPlusIOWorker::~gtPlusIOWorker()
{
    if ( !close() )
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOWorker::~gtPlusIOWorker() ... ");
    }
}

bool gtPlusIOWorker::open()
{
    try
    {
        if ( fid_.is_open() )
        {
            fid_.close();
        }

        if ( readFlag_ )
        {
            fid_.open(ioTag_.c_str(), std::ios::in | std::ios::binary);
        }
        else
        {
            fid_.open(ioTag_.c_str(), std::ios::out | std::ios::binary);
        }

        if ( !fid_ )
        {
            GADGET_ERROR_MSG("gtPlusIOWorker::open() cannot open file stream : " << ioTag_);
            return false;
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOWorker::open() ... ");
        return false;
    }

    return true;
}

bool gtPlusIOWorker::close()
{
    try
    {
        if ( fid_.is_open() )
        {
            fid_.close();
        }
    }
    catch (...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusIOWorker::close() ... ");
        return false;
    }

    return true;
}

long gtPlusIOWorker::tell()
{
    if ( !fid_.is_open() ) return -1;

    if ( readFlag_ )
    {
        return (long)fid_.tellg();
    }

    return (long)fid_.tellp();
}

bool gtPlusIOWorker::seek(long long offset)
{
    if ( !fid_.is_open() ) return false;

    if ( readFlag_ )
    {
        fid_.seekg(offset, std::ios::beg);
        return this->IOinError();
    }

    fid_.seekp(offset, std::ios::beg);
    return this->IOinError();
}

bool gtPlusIOWorker::IOinError()
{
    std::ios::iostate s;
    s = fid_.rdstate();

    if ( (s&std::ios::failbit) || (s&std::ios::badbit) )
    {
        return false;
    }

    return true;
}

bool gtPlusIOWorker::read(char* data, long long len)
{
    if ( !fid_.is_open() ) return false;
    fid_.read(data, len*sizeof(char));
    return IOinError();
}

bool gtPlusIOWorker::write(const char* data, long long len)
{
    if ( !fid_.is_open() ) return false;
    fid_.write(data, len*sizeof(char));
    return IOinError();
}

// --------------------------------------------------------------------------

//void gtPlusIOBase::printInfo(std::ostream& os)
//{
//    using namespace std;
//
//    os << "-------------- GTPlus IO Util ---------------" << endl;
//    os << "Implementation of file input/output operations" << endl;
//    os << "---------------------------------------------" << endl;
//}
//
//bool gtPlusIOBase::readFromFile(const std::string& filename, char*& data, long long& length)
//{
//    try
//    {
//        if (data!=NULL) delete [] data;
//
//        gtPlusIOWorker ioworker_(filename, true);
//
//        GADGET_CHECK_RETURN_FALSE(ioworker_.open());
//
//        // read the total length
//        long long totalLen;
//        GADGET_CHECK_RETURN_FALSE(ioworker_.read(reinterpret_cast<char*>(&totalLen), sizeof(long long)));
//
//        length = totalLen - sizeof(long long);
//
//        data = new char[length];
//        GADGET_CHECK_RETURN_FALSE(data!=NULL);
//
//        GADGET_CHECK_RETURN_FALSE(ioworker_.read(data, length));
//
//        GADGET_CHECK_RETURN_FALSE(ioworker_.close());
//    }
//    catch (...)
//    {
//        GADGET_ERROR_MSG("Errors in gtPlusIOBase::readFromFile(const std::string& filename, char*& data, long long& length) ... ");
//        return false;
//    }
//
//    return true;
//}
//
//bool gtPlusIOBase::writeToFile(const std::string& filename, char* data, long long length)
//{
//    try
//    {
//        if ( length == 0 ) return true;
//
//        GADGET_CHECK_RETURN_FALSE(data!=NULL);
//
//        gtPlusIOWorker ioworker_(filename, false);
//
//        GADGET_CHECK_RETURN_FALSE(ioworker_.open());
//
//        // write the total lengh
//        const long long totalLen = length+sizeof(long long);
//        GADGET_CHECK_RETURN_FALSE(ioworker_.write(reinterpret_cast<const char*>(&totalLen), sizeof(long long)));
//
//        // write the data
//        GADGET_CHECK_RETURN_FALSE(ioworker_.write(data, length));
//
//        // close the file
//        GADGET_CHECK_RETURN_FALSE(ioworker_.close());
//    }
//    catch (...)
//    {
//        GADGET_ERROR_MSG("Errors in gtPlusIOBase::writeToFile(const std::string& filename, char* data, long long length) ... ");
//        return false;
//    }
//
//    return true;
//}

}}
