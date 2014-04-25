
#include "gtPlusUtil.h"

namespace Gadgetron { namespace gtPlus {

// ------------------------------------------------------------------------
// random generator
// ------------------------------------------------------------------------

#ifdef _WIN32

template <typename T> 
gtPlusRandNorm<T>::gtPlusRandNorm()
{
    rng_.seed();
    this->setPara(0, 1);
}

template <typename T> 
gtPlusRandNorm<T>::gtPlusRandNorm(long long s, T mean, T sigma)
{
    this->seed(s);
    this->setPara(mean, sigma);
}

template <typename T> 
gtPlusRandNorm<T>::~gtPlusRandNorm()
{
}

template <typename T> 
void gtPlusRandNorm<T>::seed(unsigned long s)
{
    rng_.seed(s);
}

template <typename T> 
void gtPlusRandNorm<T>::setPara(T mean, T sigma)
{
    typename std::normal_distribution<T>::param_type para(mean, sigma);
    dist_norm_.param(para);
}

template <typename T> 
inline bool gtPlusRandNorm<T>::gen(hoNDArray<T>& randNum)
{
    try
    {
        size_t N = randNum.get_number_of_elements();
        size_t n;
        for ( n=0; n<N; n++ )
        {
            randNum(n) = dist_norm_(rng_);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusRandNorm<T>::gen(hoNDArray<T>& randNum) ... ");
        return false;
    }

    return true;
}

template <typename T> 
inline bool gtPlusRandNorm<T>::gen(hoNDArray< std::complex<T> >& randNum)
{
    try
    {
        size_t N = randNum.get_number_of_elements();
        size_t n;

        T real, imag;
        for ( n=0; n<N; n++ )
        {
            real = dist_norm_(rng_);
            imag = dist_norm_(rng_);

            randNum(n) = std::complex<T>(real, imag);
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusRandNorm<T>::gen(hoNDArray< std::complex<T> >& randNum) ... ");
        return false;
    }

    return true;
}
#endif // _WIN32

// ------------------------------------------------------------------------
// utility functions for various things
// ------------------------------------------------------------------------

template <typename T> 
bool gtPlusUtil<T>::getCurrentTime(size_t time[6])
{
    try
    {
        time_t rawtime;
        struct tm* timeinfo;

        std::time(&rawtime);
        timeinfo = std::gmtime (&rawtime);

        time[0] = timeinfo->tm_year+1900;
        time[1] = timeinfo->tm_mon+1;
        time[2] = timeinfo->tm_mday;
        time[3] = timeinfo->tm_hour;
        time[4] = timeinfo->tm_min;
        time[5] = timeinfo->tm_sec;
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Error happened in gtPlusUtil<T>::getCurrentTime(size_t time[6]) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusUtil<T>::convertTimeToUTC(size_t time[6], double& tmUTC)
{
    try
    {
        struct tm timeinfo;

        timeinfo.tm_year   = time[0]-1900;
        timeinfo.tm_mon    = time[1] - 1;
        timeinfo.tm_mday   = time[2];
        timeinfo.tm_hour   = time[3];
        timeinfo.tm_min    = time[4];
        timeinfo.tm_sec    = time[5];
        timeinfo.tm_isdst  = 0;

        tmUTC = (double)mktime(&timeinfo);
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Error happened in gtPlusUtil<T>::convertTimeToUTC(size_t time[6], double& tmUTC) ... ");
        return false;
    }

    return true;
}

}}
