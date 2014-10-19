
#ifdef USE_OMP
    #include "omp.h"
#endif // USE_OMP

#include "iostream"
#include "sstream"
#include "fstream"

#include <gtest/gtest.h>

#include "hoNDArray_elemwise.h"
#include "hoNDMath_util.h"
#include "GadgetronTimer.h"
#include "ho2DArray.h"

#include <boost/thread/mutex.hpp>

#ifdef USE_MKL
    #include "mkl.h"

    #ifdef MKL_ILP64
        #define ILP_MODE_ON 1
    #endif
#endif // USE_MKL

#ifndef lapack_int
    #ifdef USE_MKL
        #ifdef MKL_ILP64
            #define lapack_int __int64
        #else
            #define lapack_int int
        #endif // MKL_ILP64
    #else
        #define lapack_int int
    #endif // USE_MKL
#endif // lapack_int

using namespace Gadgetron;
using testing::Types;

const size_t REP = 1000;

/// Computes the sum of magnitudes of the vector elements.
extern "C" {

    float snrm2_(lapack_int* N, float* x, lapack_int* incx);
    float scnrm2_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
    double dnrm2_(lapack_int* N, double* x, lapack_int* incx);
    double dznrm2_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);

    float sasum_(lapack_int* N, float* x, lapack_int* incx);
    float scasum_(lapack_int* N, lapack_complex_float* x, lapack_int* incx);
    double dasum_(lapack_int* N, double* x, lapack_int* incx);
    double dzasum_(lapack_int* N, lapack_complex_double* x, lapack_int* incx);

    void cdotc_(void* r, lapack_int* N, lapack_complex_float* x, lapack_int* incx, lapack_complex_float* y, lapack_int* incy);
    void zdotc_(void* r, lapack_int* N, lapack_complex_double* x, lapack_int* incx, lapack_complex_double* y, lapack_int* incy);
}

template <typename T> class math_speed_test : public ::testing::Test 
{
protected:
    virtual void SetUp()
    {
        GADGET_MSG("=============================================================================================");
        GADGET_MSG("Unit Test for speed of math utils");
        GADGET_MSG("=============================================================================================");

        this->timer_.set_timing_in_destruction(false);

    #ifdef USE_OMP
        this->max_num_thread_ = omp_get_num_procs();
    #else
        this->max_num_thread_ = 1;
    #endif // USE_OMP

        GADGET_MSG("max_num_thread_ is " << this->max_num_thread_);

        GADGET_MSG("Initializing the testing data ... ");

        size_t N = 12;

        num_of_elements_.resize(N);
        num_of_elements_[0] = 128;

        size_t n;
        for ( n=1; n<N; n++ )
        {
            num_of_elements_[n] = 2*num_of_elements_[n-1];
        }

        A_.resize(N);
        B_.resize(N);
        res_.resize(N);

        for ( n=0; n<N; n++ )
        {
            GADGET_MSG("Allocate " << num_of_elements_[n] << " elements - " << num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb ...");

            GADGET_START_TIMING(this->timer_, "Allocate ... ");
            A_[n].create(num_of_elements_[n]); Gadgetron::clear(A_[n]);
            B_[n].create(num_of_elements_[n]); Gadgetron::math::fill( B_[n].get_number_of_elements(), B_[n].begin(), std::complex<float>(2, 4) );
            res_[n].create(num_of_elements_[n]);
            GADGET_STOP_TIMING(this->timer_);
        }
    }

    std::vector< size_t > num_of_elements_;
    std::vector< hoNDArray<std::complex<float> > > A_;
    std::vector< hoNDArray<std::complex<float> > > B_;
    std::vector< hoNDArray<std::complex<float> > > res_;

    int max_num_thread_;
    GadgetronTimer timer_;
};

typedef Types<float, double> realImplementations;
typedef Types< std::complex<float> > cpfloatImplementations;

TYPED_TEST_CASE(math_speed_test, cpfloatImplementations);

#define NumElementsUseThreading 128*128*6

void save_timing(const std::string& testName, const std::vector< size_t >& num_of_elements, const ho2DArray<double>& timing)
{
    std::stringstream str;

    str << testName<< std::endl;
    str << "===================================================================" << std::endl;

    size_t N = timing.get_size(0);
    size_t P = timing.get_size(1);

    str << "Number of arrays : " << N << std::endl;
    str << "Number of processors : " << P << std::endl;
    str << "===================================================================" << std::endl;

    size_t n, p;

    for ( n=0; n<N; n++ )
    {
        str << "---------------------------------------------------------------" << std::endl;
        str << "number of array elements " << num_of_elements[n] << std::endl;
        str << "--------------------------" << std::endl;
        for ( p=0; p<P; p++ )
        {
            str << "number of cores " << p+1 << "  -  " << timing(n, p)/1000.0 << " ms " << std::endl;
        }
    }
    str << std::endl;

    std::string filename = testName + ".txt";
    std::ofstream f;
    f.open (filename.c_str(), std::ios::out|std::ios::binary);
    std::string content = str.str();
    f.write(content.c_str(), content.size());
    f.close();
}

/// ============================================================================================================

void add(size_t N, const std::complex<float>* x, const std::complex<float>* y, std::complex<float>* r, int numOfThreads)
{
    long long n;

    if ( numOfThreads > 1 )
    {
        #pragma omp parallel for private(n) shared(N, r, x, y) num_threads(numOfThreads) 
        for ( n=0; n<(long long)N; ++n)
        {
            const std::complex<float>& vx = x[n];
            const float re1 = vx.real();
            const float im1 = vx.imag();

            const std::complex<float>& vy = y[n];
            const float re2 = vy.real();
            const float im2 = vy.imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = re1 + re2;
            reinterpret_cast<float(&)[2]>(r[n])[1] = im1 + im2;
        }
    }
    else
    {
        for ( n=0; n<(long long)N; ++n)
        {
            const std::complex<float>& vx = x[n];
            const float re1 = vx.real();
            const float im1 = vx.imag();

            const std::complex<float>& vy = y[n];
            const float re2 = vy.real();
            const float im2 = vy.imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = re1 + re2;
            reinterpret_cast<float(&)[2]>(r[n])[1] = im1 + im2;
        }
    }
}

void add_thread(size_t N, const std::complex<float>* x, const std::complex<float>* y, std::complex<float>* r, int numOfThreads)
{
    long long n;

    #pragma omp parallel for private(n) shared(N, r, x, y) if(numOfThreads>1) num_threads(numOfThreads) 
    for ( n=0; n<(long long)N; ++n)
    {
        const std::complex<float>& vx = x[n];
        const float re1 = vx.real();
        const float im1 = vx.imag();

        const std::complex<float>& vy = y[n];
        const float re2 = vy.real();
        const float im2 = vy.imag();

        reinterpret_cast<float(&)[2]>(r[n])[0] = re1 + re2;
        reinterpret_cast<float(&)[2]>(r[n])[1] = im1 + im2;
    }
}

#ifdef USE_MKL

void add_mkl(size_t N, const std::complex<float>* x, const std::complex<float>* y, std::complex<float>* r)
{
    vcAdd((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
    GADGET_CHECK_THROW(vmlGetErrStatus()==0);
}

#endif // USE_MKL

TYPED_TEST(math_speed_test, add)
{
    size_t N = this->A_.size();

    size_t n, p, r;

    GadgetronTimer timer;
    timer.set_timing_in_destruction(false);

    ho2DArray<double> time_used(N, this->max_num_thread_);

    GADGET_MSG("-----------------> Loop add <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        add(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), 2);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);

            timer.start();

            for ( r=0; r<REP; r++ )
            {
                add(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("add", this->num_of_elements_, time_used);

    GADGET_MSG("-----------------> Loop add thread <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        add(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), 2);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);

            timer.start();

            for ( r=0; r<REP; r++ )
            {
                add_thread(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("add_thread", this->num_of_elements_, time_used);

#ifdef USE_MKL

    GADGET_MSG("-----------------> mkl add <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        add_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin());

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                add_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin());
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("add_mkl", this->num_of_elements_, time_used);

#endif // USE_MKL
}

/// ============================================================================================================

void multiply(size_t N, const std::complex<float>* x, const std::complex<float>* y, std::complex<float>* r, int numOfThreads)
{
    long long i;

    if ( numOfThreads > 1 )
    {
        #pragma omp parallel for private(i) num_threads(numOfThreads)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& a1 = x[i];
            const std::complex<float>& b1 = y[i];
            const float a = a1.real();
            const float b = a1.imag();
            const float c = b1.real();
            const float d = b1.imag();

            reinterpret_cast<float(&)[2]>(r[i])[0] = a*c-b*d;
            reinterpret_cast<float(&)[2]>(r[i])[1] = a*d+b*c;
        }
    }
    else
    {
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& a1 = x[i];
            const std::complex<float>& b1 = y[i];
            const float a = a1.real();
            const float b = a1.imag();
            const float c = b1.real();
            const float d = b1.imag();

            reinterpret_cast<float(&)[2]>(r[i])[0] = a*c-b*d;
            reinterpret_cast<float(&)[2]>(r[i])[1] = a*d+b*c;
        }
    }
}

#ifdef USE_MKL
void multiply_mkl(size_t N, const std::complex<float>* x, const std::complex<float>* y, std::complex<float>* r)
{
    vcMul((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
    GADGET_CHECK_THROW(vmlGetErrStatus()==0);
}
#endif USE_MKL

TYPED_TEST(math_speed_test, multiply)
{
    size_t N = this->A_.size();

    size_t n, p, r;

    GadgetronTimer timer;
    timer.set_timing_in_destruction(false);

    ho2DArray<double> time_used(N, this->max_num_thread_);

    GADGET_MSG("-----------------> Loop multiply <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        multiply(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), 1);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                multiply(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("multiply", this->num_of_elements_, time_used);

#ifdef USE_MKL
    GADGET_MSG("-----------------> mkl multiply <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        multiply_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin());

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                multiply_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin());
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("multiply_mkl", this->num_of_elements_, time_used);
#endif // USE_MKL
}

/// ============================================================================================================

void multiplyConj(size_t N, const std::complex<float>* x, const std::complex<float>* y, std::complex<float>* r, int numOfThreads)
{
    long long n;

    if ( numOfThreads > 1 )
    {
    #pragma omp parallel for private(n) shared(N, x, y, r) num_threads(numOfThreads)
        for ( n=0; n<(long long)N; n++ )
        {
            const float a = x[n].real();
            const float b = x[n].imag();
            const float c = y[n].real();
            const float d = y[n].imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = (a*c + b*d);
            reinterpret_cast<float(&)[2]>(r[n])[1] = (c*b - a*d);
        }
    }
    else
    {
        for ( n=0; n<(long long)N; n++ )
        {
            const float a = x[n].real();
            const float b = x[n].imag();
            const float c = y[n].real();
            const float d = y[n].imag();

            reinterpret_cast<float(&)[2]>(r[n])[0] = (a*c + b*d);
            reinterpret_cast<float(&)[2]>(r[n])[1] = (c*b - a*d);
        }
    }
}

#ifdef USE_MKL
void multiplyConj_mkl(size_t N, const std::complex<float>* x, const std::complex<float>* y, std::complex<float>* r)
{
    vcMulByConj((lapack_int)N, (MKL_Complex8*)(x), (MKL_Complex8*)(y), (MKL_Complex8*)(r));
    GADGET_CHECK_THROW(vmlGetErrStatus()==0);
}
#endif // USE_MKL

TYPED_TEST(math_speed_test, multiplyConj)
{
    size_t N = this->A_.size();

    size_t n, p, r;

    GadgetronTimer timer;
    timer.set_timing_in_destruction(false);

    ho2DArray<double> time_used(N, this->max_num_thread_);

    GADGET_MSG("-----------------> Loop multiplyConj <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        multiplyConj(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), 1);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                multiplyConj(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin(), p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("multiplyConj", this->num_of_elements_, time_used);

#ifdef USE_MKL
    GADGET_MSG("-----------------> mkl multiplyConj <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        multiplyConj_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin());

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                multiplyConj_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), this->res_[n].begin());
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("multiplyConj_mkl", this->num_of_elements_, time_used);
#endif // USE_MKL
}

/// ============================================================================================================

void norm2(size_t N, const std::complex<float>* x, float& r, int numOfThreads)
    {
        long long i;

        float sum(0);

        long long num = (long long)N;

        if( numOfThreads > 1 )
        {
            #pragma omp parallel for reduction(+:sum) num_threads(numOfThreads)
            for (i = 0; i < (long long)N; i++)
            {
                const std::complex<float>& c = x[i];
                const float re = c.real();
                const float im = c.imag();
                sum += ( (re*re) + (im * im) );
            }
        }
        else
        {
            for (i = 0; i < num; i++)
            {
                const std::complex<float>& c = x[i];
                const float re = c.real();
                const float im = c.imag();
                sum += ( (re*re) + (im * im) );
            }
        }

        r = std::sqrt(sum);
    }

void norm2_blas(size_t N, const std::complex<float>* x, float& r)
{
    lapack_int num = (lapack_int)N;
    lapack_int incx = 1;

    r = scnrm2_(&num, (lapack_complex_float*)(x), &incx);
    GADGET_CHECK_THROW(vmlGetErrStatus()==0);
}

#ifdef USE_MKL
void norm2_mkl(size_t N, const std::complex<float>* x, float& r)
{
    int num = (int)N;
    int incx = 1;

    r = scnrm2(&num, (MKL_Complex8*)(x), &incx);
    GADGET_CHECK_THROW(vmlGetErrStatus()==0);
}
#endif // USE_MKL

TYPED_TEST(math_speed_test, norm2)
{
    size_t N = this->A_.size();

    size_t n, p, r;

    GadgetronTimer timer;
    timer.set_timing_in_destruction(false);

    ho2DArray<double> time_used(N, this->max_num_thread_);

    float res;

    GADGET_MSG("-----------------> Loop norm2 <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        norm2(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res, 1);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                norm2(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res, p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("norm2", this->num_of_elements_, time_used);

    GADGET_MSG("-----------------> blas norm2 <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        norm2_blas(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                norm2_blas(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("norm2_blas", this->num_of_elements_, time_used);

#ifdef USE_MKL
    GADGET_MSG("-----------------> mkl norm2 <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        norm2_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                norm2_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("norm2_mkl", this->num_of_elements_, time_used);
#endif // USE_MKL
}

/// ============================================================================================================
void norm1(size_t N, const std::complex<float>* x, float& r, int numOfThreads)
{
    long long i;
    float sum = 0.0f;

    if( numOfThreads > 1 )
    {
        #pragma omp parallel for reduction(+:sum) num_threads(numOfThreads)
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += std::sqrt( (re*re) + (im * im) );
        }
    }
    else
    {
        for (i = 0; i < (long long)N; i++)
        {
            const std::complex<float>& c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += std::sqrt( (re*re) + (im * im) );
        }
    }

    r = sum;
}

TYPED_TEST(math_speed_test, norm1)
{
    size_t N = this->A_.size();

    size_t n, p, r;

    GadgetronTimer timer;
    timer.set_timing_in_destruction(false);

    ho2DArray<double> time_used(N, this->max_num_thread_);

    float res;

    GADGET_MSG("-----------------> Loop norm1 <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        norm1(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res, 2);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                norm1(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res, p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("norm1", this->num_of_elements_, time_used);
}
/// ============================================================================================================
void dotc(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r, int numOfThreads)
    {
        long long n;

        float sum(0);

        float sa(0), sb(0);

        if ( numOfThreads > 1 )
        {
            #pragma omp parallel for private(n) reduction(+:sa) num_threads(numOfThreads)
            for (n = 0; n < (long long)N; n++)
            {
                const float a = x[n].real();
                const float b = x[n].imag();
                const float c = y[n].real();
                const float d = y[n].imag();

                sa += (a*c + b*d);
                sb += (c*b - a*d);
            }
        }
        else
        {
            for (n = 0; n < (long long)N; n++)
            {
                const float a = x[n].real();
                const float b = x[n].imag();
                const float c = y[n].real();
                const float d = y[n].imag();

                sa += (a*c + b*d);
                sb += (c*b - a*d);
            }
        }

        reinterpret_cast<float(&)[2]>(r)[0] = sa;
        reinterpret_cast<float(&)[2]>(r)[1] = sb;
    }

void dotc_blas(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r)
{
    lapack_int num = (lapack_int)N;
    lapack_int incx=1, incy=1;
    cdotc_((lapack_complex_float*)(&r), &num, (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(y), &incy);
}

#ifdef USE_MKL
void dotc_mkl(size_t N, const GT_Complex8* x, const GT_Complex8* y, GT_Complex8& r)
{
    lapack_int num = (lapack_int)N;
    lapack_int incx=1, incy=1;
    cdotc((lapack_complex_float*)(&r), &num, (lapack_complex_float*)(x), &incx, (lapack_complex_float*)(y), &incy);
}
#endif // USE_MKL

TYPED_TEST(math_speed_test, dotc)
{
    size_t N = this->A_.size();

    size_t n, p, r;

    GadgetronTimer timer;
    timer.set_timing_in_destruction(false);

    ho2DArray<double> time_used(N, this->max_num_thread_);

    GT_Complex8 res;

    GADGET_MSG("-----------------> Loop dotc <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        dotc(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), res, 1);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                dotc(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), res, p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("dotc", this->num_of_elements_, time_used);

    GADGET_MSG("-----------------> blas dotc <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        dotc_blas(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), res);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                dotc_blas(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), res);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("dotc_blas", this->num_of_elements_, time_used);

#ifdef USE_MKL
    GADGET_MSG("-----------------> mkl dotc <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        dotc_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), res);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                dotc_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), this->B_[n].begin(), res);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("dotc_mkl", this->num_of_elements_, time_used);
#endif // USE_MKL
}

/// ============================================================================================================
void asum(size_t N, const GT_Complex8* x, float& r, int numOfThreads)
{
    long long i;
    float sum(0);

    if ( numOfThreads > 1 )
    {
        for (i = 0; i < (long long)N; i++)
        {
            const GT_Complex8& c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += ( std::abs(re) + std::abs(im) );
        }
    }
    else
    {
        #pragma omp parallel for private(i) reduction(+:sum) num_threads(numOfThreads)
        for (i = 0; i < (long long)N; i++)
        {
            const GT_Complex8& c = x[i];
            const float re = c.real();
            const float im = c.imag();
            sum += ( std::abs(re) + std::abs(im) );
        }
    }

    r = sum;
}

void asum_blas(size_t N, const GT_Complex8* x, float& r)
{
    lapack_int num = (lapack_int)(N);
    lapack_int incx = 1;
    r = scasum_(&num, (lapack_complex_float*)(x), &incx);
}

#ifdef USE_MKL
void asum_mkl(size_t N, const GT_Complex8* x, float& r)
{
    lapack_int num = (lapack_int)(N);
    lapack_int incx = 1;
    r = scasum(&num, (lapack_complex_float*)(x), &incx);
}
#endif // USE_MKL

TYPED_TEST(math_speed_test, asum)
{
    size_t N = this->A_.size();

    size_t n, p, r;

    GadgetronTimer timer;
    timer.set_timing_in_destruction(false);

    ho2DArray<double> time_used(N, this->max_num_thread_);

    float res;

    GADGET_MSG("-----------------> Loop asum <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        asum(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res, 1);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                asum(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res, p+1);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("asum", this->num_of_elements_, time_used);

    GADGET_MSG("-----------------> blas asum <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        asum_blas(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                asum_blas(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("asum_blas", this->num_of_elements_, time_used);

#ifdef USE_MKL
    GADGET_MSG("-----------------> mkl asum <--------------------------------------");
    for ( n=0; n<N; n++ )
    {
        asum_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);

        for ( p=0; p<this->max_num_thread_; p++ )
        {
            GADGET_MSG("Array size : " << this->num_of_elements_[n]*sizeof(std::complex<float>)/1024.0/1024 << "mb - number of cores : " << p+1);
            timer.start();

            for ( r=0; r<REP; r++ )
            {
                asum_mkl(this->A_[n].get_number_of_elements(), this->A_[n].begin(), res);
            }

            time_used(n, p) = timer.stop();
        }
    }

    save_timing("asum_mkl", this->num_of_elements_, time_used);
#endif // USE_MKL
}

/// ============================================================================================================
