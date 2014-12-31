/** \file   gtPlusUtil.h
    \brief  Define the symbols and implement common functionalities for GtPlus toolbox

    \author Hui Xue
*/

#pragma once

#include "GtPlusExport.h"

#include "ho2DArray.h"
#include "ho3DArray.h"
#include "ho4DArray.h"
#include "ho5DArray.h"
#include "ho6DArray.h"
#include "ho7DArray.h"
#include "hoMatrix.h"
#include "hoNDFFT.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDImage_util.h"
#include "gtPlusIOAnalyze.h"
#include "hoNDArrayMemoryManaged.h"
#include "GadgetronTimer.h"

#ifdef _WIN32
    #include <random>
    #include <array>
#endif // _WIN32

#ifdef USE_OMP
    #include <omp.h>
#endif // USE_OMP

#include "mri_core_definition.h"

namespace Gadgetron { namespace gtPlus {

// ------------------------------------------------------------------------
// random generator
// ------------------------------------------------------------------------

#ifdef _WIN32

/// norm distribution random number generator
template <typename T> 
class gtPlusRandNorm
{
public:

    typedef std::mt19937 RandomGeneratorType;

    gtPlusRandNorm();
    gtPlusRandNorm(long long seed, T mean=0, T sigma=1);
    ~gtPlusRandNorm();

    void seed(unsigned long seed);
    void setPara(T mean=0, T sigma=1);

    RandomGeneratorType& getRandomer() { return rng_; }
    const RandomGeneratorType& getRandomer() const { return rng_; }

    bool gen(hoNDArray<T>& randNum);
    bool gen(hoNDArray< std::complex<T> >& randNum);

protected:

    RandomGeneratorType rng_;
    std::normal_distribution<T> dist_norm_;
};

#endif // _WIN32

template <typename T> 
class gtPlusUtil
{
public:

    gtPlusUtil() {}
    ~gtPlusUtil() {}

    // ------------------------------------------------------------------------
    // utility functions for various things
    // ------------------------------------------------------------------------

    /// get the current time in system
    /// time stores year, month, date, hour, minute and second
    bool getCurrentTime(size_t time[6]);

    /// get UTC (Coordinated Universal Time) time from current time
    bool convertTimeToUTC(size_t time[6], double& tmUTC);
};

}}

#include "gtPlusUtil.hxx"
