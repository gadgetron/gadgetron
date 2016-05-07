/** \file       ImageIOBase.cpp
    \brief      Define the base IO funcatinality for image_io toolbox
    \author     Hui Xue
*/

#include "ImageIOBase.h"

namespace Gadgetron { 

ImageIOWorker::ImageIOWorker(const std::string& ioTag, bool readFlag) : ioTag_(ioTag), readFlag_(readFlag)
{
}

ImageIOWorker::~ImageIOWorker()
{
    if ( !close() )
    {
        GERROR_STREAM("Errors in ImageIOWorker::~ImageIOWorker() ... ");
    }
}

bool ImageIOWorker::open()
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
            GERROR_STREAM("ImageIOWorker::open() cannot open file stream : " << ioTag_);
            return false;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in ImageIOWorker::open() ... ");
        return false;
    }

    return true;
}

bool ImageIOWorker::close()
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
        GERROR_STREAM("Errors in ImageIOWorker::close() ... ");
        return false;
    }

    return true;
}

long ImageIOWorker::tell()
{
    if ( !fid_.is_open() ) return -1;

    if ( readFlag_ )
    {
        return (long)fid_.tellg();
    }

    return (long)fid_.tellp();
}

bool ImageIOWorker::seek(long long offset)
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

bool ImageIOWorker::IOinError()
{
    std::ios::iostate s;
    s = fid_.rdstate();

    if ( (s&std::ios::failbit) || (s&std::ios::badbit) )
    {
        return false;
    }

    return true;
}

bool ImageIOWorker::read(char* data, long long len)
{
    if ( !fid_.is_open() ) return false;
    fid_.read(data, len*sizeof(char));
    return IOinError();
}

bool ImageIOWorker::write(const char* data, long long len)
{
    if ( !fid_.is_open() ) return false;
    fid_.write(data, len*sizeof(char));
    return IOinError();
}

}
