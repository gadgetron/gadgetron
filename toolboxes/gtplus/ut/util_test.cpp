
#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#include "Gadget.h"
#include "Gadgetron.h"
#include "ismrmrd.h"
#include "hoNDArray_elemwise.h"
#include "complext.h"

#include <gtest/gtest.h>

#include "hoNDArray_utils.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
// #include "gtPlusISMRMRDReconWorkOrder.h"
#include "gtPlusISMRMRDReconWorker2DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker2DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorker3DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian2DT.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian3DT.h"
#include "gtPlusMemoryManager.h"
#include "hoNDArrayMemoryManaged.h"
#include "gtPlusSPIRIT2DOperator.h"
#include "gtPlusSPIRIT2DTOperator.h"
#include "gtPlusSPIRIT3DOperator.h"
#include "gtPlusSPIRITNoNullSpace2DOperator.h"
#include "gtPlusSPIRITNoNullSpace2DTOperator.h"
#include "gtPlusSPIRITNoNullSpace3DOperator.h"
#include "gtPlusNCGSolver.h"

#include "GadgetronTimer.h"

#include <boost/thread/mutex.hpp>

#ifdef max
#undef max
#endif // max

using namespace Gadgetron;
using namespace Gadgetron::gtPlus;
using testing::Types;

template <typename T> class gtPlus_IO_Test : public ::testing::Test 
{
protected:
    virtual void SetUp()
    {
        GADGET_MSG("=============================================================================================");
        gtPluse_ut_folder_ = std::string(::getenv("GTPLUS_UNITTEST_DIRECTORY"));
        GADGET_MSG("=============================================================================================");
        GADGET_MSG("Unit Test for GtPlus");
        gtPluse_ut_data_folder_ = gtPluse_ut_folder_ + "/data/";
        gtPluse_ut_res_folder_ = gtPluse_ut_folder_ + "/result/";
        GADGET_MSG("gtPluse_ut_data_folder_ is " << gtPluse_ut_data_folder_);
        GADGET_MSG("gtPluse_ut_res_folder_ is " << gtPluse_ut_res_folder_);

        timer_.set_timing_in_destruction(false);

#ifdef WIN32
    #ifdef USE_OMP
        /// lock the threads
        #pragma omp parallel default(shared)
        {
            int tid = omp_get_thread_num();
            // std::cout << tid << std::endl;
            DWORD_PTR mask = (1 << tid);
            SetThreadAffinityMask( GetCurrentThread(), mask );
        }
    #endif // USE_OMP
#endif // WIN32
    }

    std::string gtPluse_ut_folder_;
    std::string gtPluse_ut_data_folder_;
    std::string gtPluse_ut_res_folder_;

    gtPlusIOAnalyze gt_io_;
    gtPlusISMRMRDReconUtil<T> util_;
    gtPlusISMRMRDReconUtilComplex<T> utilCplx_;
    GadgetronTimer timer_;
};

typedef Types<float, double> realImplementations;

typedef Types< std::complex<float> > cpfloatImplementations;

typedef Types<std::complex<float>, std::complex<double>, float_complext, double_complext> cplxImplementations;
typedef Types<std::complex<float>, std::complex<double> > stdCplxImplementations;
typedef Types<float_complext, double_complext> cplxtImplementations;

TYPED_TEST_CASE(gtPlus_IO_Test, cpfloatImplementations);

TYPED_TEST(gtPlus_IO_Test, recon2DCoilMapGPU)
{
    typedef GT_Complex8 T;

    gtPlusIOAnalyze gt_io;

    float v;

    // image data
    hoNDArray<GT_Complex8> data;
    // gt_io.importArrayComplex(data, this->gtPluse_ut_data_folder_ + "fullkspace__REAL", this->gtPluse_ut_data_folder_ + "fullkspace__IMAG");
    gt_io.importArrayComplex(data, this->gtPluse_ut_data_folder_ + "aveComplexIm_REAL", this->gtPluse_ut_data_folder_ + "aveComplexIm_IMAG");
    data.print(std::cout);

    data.squeeze();

    GadgetronTimer timer(false);

    unsigned int RO = data.get_size(0);
    unsigned int E1 = data.get_size(1);
    unsigned int CHA = data.get_size(2);
    unsigned int N = data.get_size(3);

    Gadgetron::norm2(data, v);
    GADGET_MSG("data = " << v);

    {
        GPUTimer t("all steps");
    }

    hoNDArray<GT_Complex8> data2D(RO, E1, CHA, data.begin());

    hoNDArray<T> CoilMap2D;
    timer.start("coilMap2DNIHGPU 2D");
    gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIHGPU(data2D, CoilMap2D, ISMRMRD_SOUHEIL, 7, 3, 3, 1e-3);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, CoilMap2D, "CoilMap2D_1");

    {
    // call the old coil map code
    timer.start("coilMap2DNIHGPU 2D old");
    hoNDArray<float_complext> host_data(RO, E1, CHA, reinterpret_cast<float_complext*>(data2D.begin()));
    cuNDArray<float_complext> device_data(host_data);
    boost::shared_ptr< cuNDArray<float_complext> > csm = Gadgetron::estimate_b1_map<float, 2>( &device_data, CHA);
    boost::shared_ptr< hoNDArray<float_complext> > csm_host = csm->to_host();
    memcpy(CoilMap2D.begin(), csm_host->begin(), csm_host->get_number_of_bytes());
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, CoilMap2D, "CoilMap2D_1_old");
    }

    hoNDArray<T> CoilMap;
    timer.start("coilMap2DNIHGPU");
    gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIHGPU(data, CoilMap, ISMRMRD_SOUHEIL, 7, 3, 3, 1e-3);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, CoilMap, "CoilMap2D");

    hoNDArray<T> CoilMap2;
    timer.start("coilMap2DNIH");
    gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(data2D, CoilMap2, ISMRMRD_SOUHEIL, 7, 3, 3, 1e-3, false);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, CoilMap2, "CoilMap2D_2");

    hoNDArray<T> combined;
    timer.start("coil combine");
    gtPlusISMRMRDReconUtilComplex<T>().coilCombine(data, CoilMap, combined);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, combined, "combined2D");

    cudaDeviceReset();
}

TYPED_TEST(gtPlus_IO_Test, recon3DCoilMapGPU)
{
    typedef GT_Complex8 T;

    gtPlusIOAnalyze gt_io;

    float v;

    // image data
    hoNDArray<GT_Complex8> data;
    gt_io.importArrayComplex(data, this->gtPluse_ut_data_folder_ + "fullkspace__REAL", this->gtPluse_ut_data_folder_ + "fullkspace__IMAG");
    data.print(std::cout);

    data.squeeze();

    GadgetronTimer timer(false);

    unsigned int RO = data.get_size(0);
    unsigned int E1 = data.get_size(1);
    unsigned int E2 = data.get_size(2);
    unsigned int CHA = data.get_size(3);

    Gadgetron::norm2(data, v);
    GADGET_MSG("data = " << v);

    {
        GPUTimer t("all steps");
    }

    hoNDArray<GT_Complex8> Im2;
    timer.start("ifft3c");
    hoNDFFT<float>::instance()->ifft3c(data, Im2);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, Im2, "Im2");

    hoNDArray<T> CoilMap;
    timer.start("coilMap3DNIHGPU");
    gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIHGPU_FullResMap(Im2, CoilMap, ISMRMRD_SOUHEIL, 7, 3, true);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, CoilMap, "CoilMap");

    omp_set_nested(1);

    hoNDArray<T> CoilMap2;
    timer.start("coilMap3DNIH");
    gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(Im2, CoilMap2, ISMRMRD_SOUHEIL, 7, 3, true);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, CoilMap2, "CoilMap2");

    hoNDArray<T> combined;
    timer.start("coil combine");
    gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(Im2, CoilMap, combined);
    timer.stop();
    GADGET_EXPORT_ARRAY_COMPLEX(this->gtPluse_ut_res_folder_, gt_io, combined, "combined");

    cudaDeviceReset();
}

//TYPED_TEST(gtPlus_IO_Test, reconCoilCompression)
//{
//    typedef float T;
//    typedef std::complex<T> TValueType;
//
//    gtPlusIOAnalyze gt_io;
//    GadgetronTimer timer(false);
//
//    gtPlusISMRMRDReconUtil<TValueType> util;
//    gtPlusISMRMRDReconUtilComplex<TValueType> utilCplx;
//    std::string filename;
//
//    hoNDArray<GT_Complex8> data;
//    gt_io.importArrayComplex(data, this->gtPluse_ut_data_folder_ + "refRecon_REAL", this->gtPluse_ut_data_folder_ + "refRecon_IMAG");
//    data.print(std::cout);
//
//    // export images
//    hoNDArray<GT_Complex8> complexIm;
//    Gadgetron::hoNDFFT<T>::instance()->ifft2c(data, complexIm);
//
//    hoNDArray<TValueType> sos;
//    utilCplx.sumOfSquare(complexIm, sos);
//
//    filename = this->gtPluse_ut_res_folder_ + "refRecon_SoS";
//    gt_io.exportArrayComplex(sos, filename);
//
//    hoMatrix<GT_Complex8> coeff, eigenValues;
//    utilCplx.computeKLCoilCompressionCoeff(data, 1e-3, coeff, eigenValues);
//    eigenValues.print(std::cout);
//}
//
//TYPED_TEST(gtPlus_IO_Test, MatrixComputation)
//{
//    MKL_INT n = 4, nrhs = 2, ldb = 2;
//
//    /* Local arrays */
//    //MKL_Complex8 a[16] = 
//    //{
//    //    { 5.96f,  0.00f}, { 0.40f,  -1.19f}, { -0.83f, -0.48f}, { -0.57f, 0.40f},
//    //    { 0.40f,  1.19f}, { 7.95f,  0.00f}, { 0.33f,  0.09f}, { 0.22f, 0.74f},
//    //    {-0.83f,  0.48f}, { 0.33f, -0.09f}, { 4.43f,  0.00f}, { -1.09f, 0.32f},
//    //    {-0.57f, -0.40f}, { 0.22f, -0.74f}, {-1.09f, -0.32f}, { 3.46f,  0.00f}
//    //};
//
//    //MKL_Complex8 b[8] = 
//    //{
//    //    {-2.94f,  5.79f}, { 8.44f,  3.07f},
//    //    { 8.12f, -9.12f}, { 1.00f, -4.62f},
//    //    { 9.09f, -5.03f}, { 3.64f, -2.33f},
//    //    { 7.36f,  6.77f}, { 8.04f,  2.87f}
//    //};
//
//    MKL_Complex8 a[16] = 
//    {
//        { 5.96f,  0.00f},   { 0.40f,  1.19f},   { -0.83f, 0.48f},   { -0.57f, -0.40f},
//        { 0.40f,  -1.19f},  { 7.95f,  0.00f},   { 0.33f,  -0.09f},  { 0.22f, -0.74f},
//        {-0.83f,  -0.48f},  { 0.33f, 0.09f},    { 4.43f,  0.00f},   { -1.09f, -0.32f},
//        {-0.57f,  0.40f},   { 0.22f, 0.74f},    {-1.09f, 0.32f},    { 3.46f,  0.00f}
//    };
//
//    MKL_Complex8 b[8] = 
//    {
//        {-2.94f,  5.79f}, { 8.12f, -9.12f}, { 9.09f, -5.03f}, { 7.36f,  6.77f}, 
//        { 8.44f,  3.07f}, { 1.00f, -4.62f}, { 3.64f, -2.33f}, { 8.04f,  2.87f}
//    };
//
//    hoMatrix< std::complex<float> > A(n, n, reinterpret_cast<std::complex<float>*>(a));
//    hoMatrix< std::complex<float> > B(n, ldb, reinterpret_cast<std::complex<float>*>(b));
//
//    hoMatrix< std::complex<float> > AB;
//    GeneralMatrixProduct_gemm(AB, A, false, B, false);
//    AB.print(std::cout);
//
//    GeneralMatrixProduct_gemm(AB, A, true, B, false);
//    AB.print(std::cout);
//
//    //A*B
//    //ans =
//    //                   -41.9895 +               20.0944i                    35.3342 +                17.026i
//    //                    56.5497 -               67.5926i                     8.7286 -               19.3177i
//    //                    31.5997 -               37.2643i                    -2.1214 -               10.9889i
//    //                    12.9773 +               15.8586i                    16.3236 +                4.4228i
//
//    hoMatrix< std::complex<float> > A2(A);
//    hoMatrix< std::complex<float> > B2(B);
//
//    SymmetricHermitianPositiveDefiniteLinearSystem_posv(A2, B2);
//
//    A2.print(std::cout);
//    B2.print(std::cout);
//
//    //    Solution
//    //    (  0.80,  1.62) (  2.52,  0.61)
//    //    (  1.26, -1.78) (  0.01, -1.38)
//    //    (  3.38, -0.29) (  2.42, -0.52)
//    //    (  3.46,  2.92) (  3.77,  1.37)
//
//    //    Details of Cholesky factorization
//    //    (  2.44,  0.00) (  0.00,  0.00) (  0.00,  0.00) (  0.00,  0.00)
//    //    (  0.16,  0.49) (  2.77,  0.00) (  0.00,  0.00) (  0.00,  0.00)
//    //    ( -0.34,  0.20) (  0.10, -0.10) (  2.06,  0.00) (  0.00,  0.00)
//    //    ( -0.23, -0.16) (  0.12, -0.30) ( -0.57, -0.20) (  1.71,  0.00)
//
//    A2 = A;
//    CholeskyHermitianPositiveDefinite_potrf(A2, 'L');
//    A2.print(std::cout);
//
//    A2 = A;
//    A2.print(std::cout);
//
//    hoMatrix< std::complex<float> > eigenValue;
//    EigenAnalysis_syev_heev2(A2, eigenValue);
//    A2.print(std::cout);
//    eigenValue.print(std::cout);
//
//    hoMatrix< std::complex<float> > C;
//    GeneralMatrixProduct_gemm(C, A2, false, A2, true);
//    C.print(std::cout);
//
//    A2 = A;
//    B2 = B;
//    hoMatrix< std::complex<float> > x;
//    double lamda = 1e-4;
//    SolveLinearSystem_Tikhonov(A2, B2, x, lamda);
//    x.print(std::cout);
//}
//
//TYPED_TEST(gtPlus_IO_Test, memoryManager)
//{
//    typedef GT_Complex8 T;
//
//    unsigned int RO = 256;
//    unsigned int E1 = 256;
//    unsigned int E2 = 256;
//    unsigned int CHA = 32;
//
//    size_t num = (size_t)RO*E1*E2*CHA*sizeof(T);
//    std::cout << "Allocate " << num/1024/1024 << " MegaBytes ..." << std::endl;
//
//    GadgetronTimer timer(false);
//
//    timer.start("Allocate 2D array...");
//    hoNDArray<T> a2D(RO, E1);
//    timer.stop();
//
//    timer.start("Allocate 3D array...");
//    hoNDArray<T> a3D(RO, E1, E2);
//    timer.stop();
//
//    timer.start("Allocate 3D array...");
//    T* p3D = new T[RO*E1*E2];
//    timer.stop();
//    memset(p3D, 0, sizeof(T)*RO*E1*E2);
//    delete [] p3D;
//
//    timer.start("Allocate 3D array...");
//    p3D = (T*)mkl_malloc(sizeof(T)*RO*E1*E2, 4);
//    timer.stop();
//    p3D[12] = T(2.3);
//    memset(p3D, 0, sizeof(T)*RO*E1*E2);
//
//    timer.start("Allocate 4D array...");
//    hoNDArray<T> a4D(RO, E1, E2, CHA);
//    timer.stop();
//
//    timer.start("Allocate 4D array...");
//    T* p4D = new T[RO*E1*E2*CHA];
//    timer.stop();
//    memset(p4D, 0, sizeof(T)*RO*E1*E2*CHA);
//    delete [] p4D;
//
//    timer.start("Allocate 4D array...");
//    p4D = (T*)mkl_malloc(sizeof(T)*RO*E1*E2*CHA, 4);
//    timer.stop();
//    p4D[12560] = T(2.3);
//    timer.start("Allocate 4D array...");
//    memset(p4D, 0, sizeof(T)*RO*E1*E2*CHA);
//    timer.stop();
//
//    timer.start("Allocate ...");
//    boost::shared_ptr<gtPlusMemoryManager> memMagnager(new gtPlusMemoryManager(4, num));
//    timer.stop();
//
//    timer.start("Allocate 3 pieces ...");
//    void* ptr = memMagnager->allocate(num/2);
//    ptr = memMagnager->allocate(num/4);
//    ptr = memMagnager->allocate(num/8);
//    timer.stop();
//
//    memMagnager->printInfo(std::cout);
//
//    boost::shared_ptr<gtPlusMemoryManager> memMagnager2;
//
//    if ( memMagnager2 )
//    {
//        std::cout << "Test " << std::endl;
//    }
//
//    if ( memMagnager )
//    {
//        std::cout << "Test " << std::endl;
//    }
//
//    boost::mutex mutex_;
//
//    timer.start("mutex cost ...");
//    mutex_.lock();
//    mutex_.unlock();
//    timer.stop();
//
//    std::cout << memMagnager.use_count() << std::endl;
//
//    timer.start("Allocate hoNDArrayMemoryManaged ...");
//    mutex_.lock();
//    Gadgetron::hoNDArrayMemoryManaged<T> a(256, 256, 128, memMagnager);
//    mutex_.unlock();
//    timer.stop();
//
//    std::cout << memMagnager.use_count() << std::endl;
//
//    memMagnager->printInfo(std::cout);
//
//    a.clear();
//
//    memMagnager->printInfo(std::cout);
//
//    int ii;
//    #pragma omp parallel
//    {
//        Gadgetron::hoNDArrayMemoryManaged<T> b(256, 256, memMagnager);
//    }
//
//    //timer.start("Allocate hoNDArrayMemoryManaged 2...");
//    //Gadgetron::hoNDArrayMemoryManaged<T> b(256, 256, 128, *memMagnager);
//    //timer.stop();
//}
//
//TYPED_TEST(gtPlus_IO_Test, kspaceFilter)
//{
//    typedef GT_Complex8 T;
//
//    gtPlusIOAnalyze gt_io;
//
//    gtPlusISMRMRDReconUtil<T> util;
//
//    hoNDArray<T> filter;
//
//    unsigned int len = 12;
//    double sigma = 1.5;
//    unsigned int width = len*0.15;
//
//    ISMRMRDKSPACEFILTER filterType = ISMRMRD_FILTER_NONE;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_GAUSSIAN;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_HANNING;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TUKEY;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//
//    len = 13;
//
//    filterType = ISMRMRD_FILTER_NONE;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_GAUSSIAN;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_HANNING;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TUKEY;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateSymmetricFilter(len, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//
//    len = 13;
//    unsigned int start = 0;
//    unsigned int end = 9;
//
//    filterType = ISMRMRD_FILTER_NONE;
//    util.generateAsymmetricFilter(len, start, end, filter, filterType, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateAsymmetricFilter(len, start, end, filter, filterType, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//
//    start = 4;
//    end = 12;
//
//    filterType = ISMRMRD_FILTER_NONE;
//    util.generateAsymmetricFilter(len, start, end, filter, filterType, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateAsymmetricFilter(len, start, end, filter, filterType, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//
//    len = 12;
//
//    start = 0;
//    end = 9;
//    filterType = ISMRMRD_FILTER_NONE;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_GAUSSIAN;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TUKEY;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//
//    start = 4;
//    end = len-1;
//
//    filterType = ISMRMRD_FILTER_NONE;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_GAUSSIAN;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TUKEY;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//
//    len = 13;
//
//    start = 0;
//    end = 9;
//    filterType = ISMRMRD_FILTER_NONE;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_GAUSSIAN;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TUKEY;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//
//    start = 4;
//    end = len-1;
//
//    filterType = ISMRMRD_FILTER_NONE;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_GAUSSIAN;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TUKEY;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    filterType = ISMRMRD_FILTER_TAPERED_HANNING;
//    util.generateSymmetricFilterForRef(len, start, end, filter, filterType, sigma, width);
//    filter.printContent(std::cout);
//
//    GADGET_MSG("------------------------------------------------");
//}
//
//TYPED_TEST(gtPlus_IO_Test, FFT)
//{
//    {
//        hoNDArray< std::complex<float> > A1D(7);
//        for ( unsigned int ii=0; ii<7; ii++ )
//        {
//            A1D(ii) = ii;
//        }
//
//        A1D.print(std::cout);
//
//        hoNDArray< std::complex<float> > A1Ds;
//        hoNDFFT<float>::instance()->ifftshift1D(A1D, A1Ds);
//        A1Ds.print(std::cout);
//
//        hoNDFFT<float>::instance()->fftshift1D(A1D, A1Ds);
//        A1Ds.print(std::cout);
//
//        hoNDFFT<float>::instance()->ifftshift1D(A1Ds, A1D);
//        A1D.print(std::cout);
//
//        hoNDArray< std::complex<float> > AR(A1D);
//        hoNDFFT<float>::instance()->fft1(A1D, AR);
//        AR.print(std::cout);
//
//        //0 = (7.937254,0.000000)
//        //1 = (-1.322875,2.746980)
//        //2 = (-1.322876,1.054958)
//        //3 = (-1.322875,0.301938)
//        //4 = (-1.322875,-0.301938)
//        //5 = (-1.322876,-1.054958)
//        //6 = (-1.322875,-2.746980)
//
//        hoNDFFT<float>::instance()->ifft1(A1D, AR);
//        AR.print(std::cout);
//
//        //0 = (7.937254,0.000000)
//        //1 = (-1.322875,-2.746980)
//        //2 = (-1.322876,-1.054958)
//        //3 = (-1.322875,-0.301938)
//        //4 = (-1.322875,0.301938)
//        //5 = (-1.322876,1.054958)
//        //6 = (-1.322875,2.746980)
//
//        hoNDFFT<float>::instance()->fft1c(A1D, AR);
//        AR.print(std::cout);
//
//        //0 = (0.000000,1.356896)
//        //1 = (0.000000,-1.692022)
//        //2 = (0.000000,3.048917)
//        //3 = (7.937254,0.000000)
//        //4 = (0.000000,-3.048917)
//        //5 = (0.000000,1.692022)
//        //6 = (0.000000,-1.356896)
//
//        hoNDFFT<float>::instance()->ifft1c(A1D, AR);
//        AR.print(std::cout);
//
//        //0 = (0.000000,-1.356896)
//        //1 = (0.000000,1.692022)
//        //2 = (0.000000,-3.048917)
//        //3 = (7.937254,0.000000)
//        //4 = (0.000000,3.048917)
//        //5 = (0.000000,-1.692022)
//        //6 = (0.000000,1.356896)
//    }
//
//    {
//        int nx = 5, ny = 6, nz = 3;
//
//        ho3DArray< std::complex<float> > A(nx, ny, nz);
//        A.fill(2.0);
//        A(1, 4, 2) = std::complex<float>(12, -5.0);
//        A.print(std::cout);
//
//        ho3DArray< std::complex<float> > AR(A);
//        hoNDFFT<float>::instance()->fft2(AR);
//        AR.print(std::cout);
//
//        //AR(:,:,1) =
//        //    10.9545         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    AR(:,:,2) =
//        //    10.9545         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    0         0         0         0         0
//        //    AR(:,:,3) =
//        //    12.7802 - 0.9129i  -0.3040 - 2.0185i  -2.0136 - 0.3346i  -0.9405 + 1.8117i   1.4324 + 1.4543i
//        //    -0.1223 + 2.0376i   1.9001 + 0.7460i   1.2966 - 1.5765i  -1.0987 - 1.7203i  -1.9756 + 0.5133i
//        //    -1.7034 - 1.1247i  -1.5960 + 1.2725i   0.7170 + 1.9112i   2.0392 - 0.0914i   0.5433 - 1.9676i
//        //    1.8257 - 0.9129i  -0.3040 - 2.0185i  -2.0136 - 0.3346i  -0.9405 + 1.8117i   1.4324 + 1.4543i
//        //    -0.1223 + 2.0376i   1.9001 + 0.7460i   1.2966 - 1.5765i  -1.0987 - 1.7203i  -1.9756 + 0.5133i
//        //    -1.7034 - 1.1247i  -1.5960 + 1.2725i   0.7170 + 1.9112i   2.0392 - 0.0914i   0.5433 - 1.9676i
//
//        ho3DArray< std::complex<float> > AR_I(A);
//        hoNDFFT<float>::instance()->ifft2(AR_I);
//        AR_I.print(std::cout);
//
//        //AR_I(:,:,1) =
//        //10.9545         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //AR_I(:,:,2) =
//        //10.9545         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //AR_I(:,:,3) =
//        //12.7802 - 0.9129i   1.4324 + 1.4543i  -0.9405 + 1.8117i  -2.0136 - 0.3346i  -0.3040 - 2.0185i
//        //-1.7034 - 1.1247i   0.5433 - 1.9676i   2.0392 - 0.0914i   0.7170 + 1.9112i  -1.5960 + 1.2725i
//        //-0.1223 + 2.0376i  -1.9756 + 0.5133i  -1.0987 - 1.7203i   1.2966 - 1.5765i   1.9001 + 0.7460i
//        //1.8257 - 0.9129i   1.4324 + 1.4543i  -0.9405 + 1.8117i  -2.0136 - 0.3346i  -0.3040 - 2.0185i
//        //-1.7034 - 1.1247i   0.5433 - 1.9676i   2.0392 - 0.0914i   0.7170 + 1.9112i  -1.5960 + 1.2725i
//        //-0.1223 + 2.0376i  -1.9756 + 0.5133i  -1.0987 - 1.7203i   1.2966 - 1.5765i   1.9001 + 0.7460i
//
//        ho3DArray< std::complex<float> > ARc(A);
//        hoNDFFT<float>::instance()->fft2c(ARc);
//        ARc.print(std::cout);
//
//        //ARc(:,:,1) =
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0   10.9545         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //ARc(:,:,2) =
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0   10.9545         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //ARc(:,:,3) =
//        //2.0136 + 0.3346i   0.3040 + 2.0185i  -1.8257 + 0.9129i  -1.4324 - 1.4543i   0.9405 - 1.8117i
//        //1.2966 - 1.5765i   1.9001 + 0.7460i  -0.1223 + 2.0376i  -1.9756 + 0.5133i  -1.0987 - 1.7203i
//        //-0.7170 - 1.9112i   1.5960 - 1.2725i   1.7034 + 1.1247i  -0.5433 + 1.9676i  -2.0392 + 0.0914i
//        //-2.0136 - 0.3346i  -0.3040 - 2.0185i  12.7802 - 0.9129i   1.4324 + 1.4543i  -0.9405 + 1.8117i
//        //-1.2966 + 1.5765i  -1.9001 - 0.7460i   0.1223 - 2.0376i   1.9756 - 0.5133i   1.0987 + 1.7203i
//        //0.7170 + 1.9112i  -1.5960 + 1.2725i  -1.7034 - 1.1247i   0.5433 - 1.9676i   2.0392 - 0.0914i
//
//        ARc = A;
//        hoNDFFT<float>::instance()->ifft2c(ARc);
//        ARc.print(std::cout);
//
//        //ARc(:,:,1) =
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0   10.9545         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //ARc(:,:,2) =
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //0         0   10.9545         0         0
//        //0         0         0         0         0
//        //0         0         0         0         0
//        //ARc(:,:,3) =
//        //0.9405 - 1.8117i  -1.4324 - 1.4543i  -1.8257 + 0.9129i   0.3040 + 2.0185i   2.0136 + 0.3346i
//        //2.0392 - 0.0914i   0.5433 - 1.9676i  -1.7034 - 1.1247i  -1.5960 + 1.2725i   0.7170 + 1.9112i
//        //1.0987 + 1.7203i   1.9756 - 0.5133i   0.1223 - 2.0376i  -1.9001 - 0.7460i  -1.2966 + 1.5765i
//        //-0.9405 + 1.8117i   1.4324 + 1.4543i  12.7802 - 0.9129i  -0.3040 - 2.0185i  -2.0136 - 0.3346i
//        //-2.0392 + 0.0914i  -0.5433 + 1.9676i   1.7034 + 1.1247i   1.5960 - 1.2725i  -0.7170 - 1.9112i
//        //-1.0987 - 1.7203i  -1.9756 + 0.5133i  -0.1223 + 2.0376i   1.9001 + 0.7460i   1.2966 - 1.5765i
//    }
//
//    {
//        int nx = 5, ny = 6, nz = 3;
//
//        ho3DArray< std::complex<float> > A(nx, ny, nz);
//        A.fill(2.0);
//        A(1, 4, 2) = std::complex<float>(12, -5.0);
//        A.print(std::cout);
//
//        ho3DArray< std::complex<float> > AR(A);
//        hoNDFFT<float>::instance()->fft3(AR);
//        AR.print(std::cout);
//
//        hoNDFFT<float>::instance()->ifft3(AR);
//        AR.print(std::cout);
//
//        ho3DArray< std::complex<float> > AR_I(A);
//        hoNDFFT<float>::instance()->ifft3(AR_I);
//        AR_I.print(std::cout);
//
//        ho3DArray< std::complex<float> > ARc(A);
//        hoNDFFT<float>::instance()->fft3c(ARc);
//        ARc.print(std::cout);
//
//        ARc = A;
//        hoNDFFT<float>::instance()->ifft3c(ARc);
//        ARc.print(std::cout);
//    }
//}
//
//TYPED_TEST(gtPlus_IO_Test, recon3D)
//{
//    typedef GT_Complex8 T;
//
//    gtPlusIOAnalyze gt_io;
//
//    float v;
//
//    std::string debugFolder;
//
//    // image data
//    hoNDArray<GT_Complex8> data;
//    gt_io.importArrayComplex(data, debugFolder + "data_dst__REAL", 
//        debugFolder + "data_dst__IMAG");
//    data.print(std::cout);
//
//    GadgetronTimer timer(false);
//
//    unsigned int RO = data.get_size(0);
//    unsigned int E1 = data.get_size(1);
//    unsigned int E2 = data.get_size(2);
//    unsigned int CHA = data.get_size(3);
//
//    Gadgetron::norm2(data, v);
//    GADGET_MSG("data = " << v);
//
//    hoNDArray<GT_Complex8> Im;
//    hoNDFFT<float>::instance()->ifft3(data, Im);
//
//    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder, gt_io, Im, "Im");
//
//    hoNDArray<GT_Complex8> Im2;
//    timer.start("ifft3c");
//    hoNDFFT<float>::instance()->ifft3c(data, Im2);
//    timer.stop();
//    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder, gt_io, Im2, "Im2");
//
//    hoNDArray<GT_Complex8> Im3(RO, E1, 4, CHA);
//    memcpy(Im3.begin(), Im2.begin(), Im3.get_number_of_bytes());
//    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder, gt_io, Im3, "Im3");
//
//    hoNDArray<T> CoilMap;
//    timer.start("coilMap3DNIH");
//    gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(Im3, CoilMap, ISMRMRD_SOUHEIL, 7, 3, true);
//    timer.stop();
//    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder, gt_io, CoilMap, "CoilMap");
//}
//
//TYPED_TEST(gtPlus_IO_Test, KLTransform)
//{
//    typedef GT_Complex8 T;
//
//    gtPlusIOAnalyze gt_io;
//
//    float v;
//
//    // image data
//    hoNDArray<float> real_data;
//    std::string filename = this->gtPluse_ut_data_folder_ + "fullkspace_REAL";
//    gt_io.importArray(real_data, filename);
//    real_data.print(std::cout);
//
//    hoNDArray<float> imag_data;
//    filename = this->gtPluse_ut_data_folder_ + "fullkspace_IMAG";
//    gt_io.importArray(imag_data, filename);
//    imag_data.print(std::cout);
//
//    boost::shared_ptr< hoNDArray<GT_Complex8> > tmp = real_imag_to_complex<GT_Complex8>(&real_data, &imag_data);
//
//    unsigned int RO = tmp->get_size(0);
//    unsigned int E1 = tmp->get_size(1);
//    unsigned int CHA = tmp->get_size(2);
//    unsigned int PHS = tmp->get_size(3);
//
//    hoNDArray<GT_Complex8> kspace(RO, E1, CHA, PHS, tmp->begin());
//
//    gtPlusISMRMRDReconUtil<GT_Complex8> util;
//    gtPlusISMRMRDReconUtilComplex<GT_Complex8> utilCplx;
//
//    hoNDArray<GT_Complex8> complexIm;
//    Gadgetron::hoNDFFT<float>::instance()->ifft2c(kspace, complexIm);
//
//    hoNDArray<GT_Complex8> complexImSoS;
//    utilCplx.sumOfSquare(complexIm, complexImSoS);
//
//    gt_io.export3DArrayComplex(complexImSoS, this->gtPluse_ut_res_folder_+"complexImSoS");
//
//    unsigned int numOfModes = 10;
//
//    hoNDArray<GT_Complex8> complexImKLF;
//    util.computeKLFilter(complexIm, numOfModes, complexImKLF);
//
//    utilCplx.sumOfSquare(complexImKLF, complexImSoS);
//
//    gt_io.export3DArrayComplex(complexImSoS, this->gtPluse_ut_res_folder_+"complexImKLFSoS");
//}
//
//TYPED_TEST(gtPlus_IO_Test, reconRemoveROOversampling)
//{
//    typedef float T;
//    typedef std::complex<T> TValueType;
//
//    gtPlusIOAnalyze gt_io;
//    GadgetronTimer timer(false);
//
//    gtPlusISMRMRDReconUtil<TValueType> util;
//    gtPlusISMRMRDReconUtilComplex<TValueType> utilCplx;
//    std::string filename;
//
//    hoNDArray<GT_Complex8> data;
//    gt_io.importArrayComplex(data, this->gtPluse_ut_data_folder_ + "kspace_DownSampleFE_real", this->gtPluse_ut_data_folder_ + "kspace_DownSampleFE_imag");
//    // real_imag_to_complex<GT_Complex8>(real_data, imag_data, data);
//    data.print(std::cout);
//
//    // export images
//    hoNDArray<GT_Complex8> complexIm;
//    Gadgetron::hoNDFFT<T>::instance()->ifft2c(data, complexIm);
//
//    hoNDArray<TValueType> sos;
//    utilCplx.sumOfSquare(complexIm, sos);
//
//    filename = this->gtPluse_ut_res_folder_ + "kspace_DownSampleFE_SoS";
//    gt_io.exportArrayComplex(sos, filename);
//
//    // cut down RO oversampling
//    hoNDArray<TValueType> dataCut;
//    Gadgetron::hoNDFFT<T>::instance()->ifft1c(data);
//    utilCplx.cutpad2D(data, data.get_size(0)/2, data.get_size(1), dataCut);
//    Gadgetron::hoNDFFT<T>::instance()->fft1c(dataCut);
//
//    Gadgetron::hoNDFFT<T>::instance()->ifft2c(dataCut, complexIm);
//    utilCplx.sumOfSquare(complexIm, sos);
//
//    filename = this->gtPluse_ut_res_folder_ + "kspace_DownSampleFE_SoS_CutRO";
//    gt_io.exportArrayComplex(sos, filename);
//}
//
//TYPED_TEST(gtPlus_IO_Test, reconNoisePrewhitening)
//{
//    typedef float T;
//    typedef std::complex<T> TValueType;
//
//    gtPlusIOAnalyze gt_io;
//    GadgetronTimer timer(false);
//
//    hoNDArray<float> real_noise;
//    std::string filename = this->gtPluse_ut_data_folder_ + "Noise_real";
//    gt_io.importArray(real_noise, filename);
//    real_noise.print(std::cout);
//
//    hoNDArray<float> imag_noise;
//    filename = this->gtPluse_ut_data_folder_ + "Noise_imag";
//    gt_io.importArray(imag_noise, filename);
//    imag_noise.print(std::cout);
//
//    ho3DArray<GT_Complex8> noise;
//    real_imag_to_complex<GT_Complex8>(real_noise, imag_noise, noise);
//
//    int COL = noise.get_size(0);
//    int E1 = noise.get_size(1);
//    int CHA = noise.get_size(2);
//
//    GADGET_MSG(noise(12, 0, 10));
//
//    // compute noise prewhitener
//    double rxDwellTimeData = 2100;
//    hoMatrix<TValueType> noisePrewhitener(CHA, CHA);
//
//    gtPlusISMRMRDReconUtilComplex<TValueType> utilCplx;
//
//    double noiseBandWidth = 130;
//    double receiverBWRatio = 0.79;
//    double ADCSamplingTimeinSecond = 2100/1e9;
//
//    hoMatrix<TValueType> prewhiteningMatrix;
//
//    GADGET_START_TIMING(timer, "computeNoisePrewhiteningMatrix");
//    utilCplx.computeNoisePrewhiteningMatrix(noise, noiseBandWidth, receiverBWRatio, ADCSamplingTimeinSecond, prewhiteningMatrix);
//    GADGET_STOP_TIMING(timer);
//    // prewhiteningMatrix.print(std::cout);
//
//    EXPECT_NEAR(prewhiteningMatrix(0, 0).real(), 5.1331672e+004, 0.01);
//    EXPECT_NEAR(prewhiteningMatrix(0, 0).imag(), 0.0, 0.01);
//
//    EXPECT_NEAR(prewhiteningMatrix(1, 0).real(), -5791.2319, 0.01);
//    EXPECT_NEAR(prewhiteningMatrix(1, 0).imag(), -1603.6230, 0.01);
//
//    EXPECT_NEAR(prewhiteningMatrix(2, 1).real(), -9597.3955, 0.01);
//    EXPECT_NEAR(prewhiteningMatrix(2, 1).imag(), 4500.7114, 0.01);
//
//    EXPECT_NEAR(prewhiteningMatrix(4, 3).real(), -7718.3286, 0.01);
//    EXPECT_NEAR(prewhiteningMatrix(4, 3).imag(), -3565.7336, 0.01);
//
//    EXPECT_NEAR(prewhiteningMatrix(31, 31).real(), 60350.840, 0.01);
//    EXPECT_NEAR(prewhiteningMatrix(31, 31).imag(), 0.0, 0.01);
//
//    /// load the data scan
//    hoNDArray<float> real_data;
//    filename = this->gtPluse_ut_data_folder_ + "noisePrewhitening_DataScan_real";
//    gt_io.importArray(real_data, filename);
//    real_data.print(std::cout);
//
//    hoNDArray<float> imag_data;
//    filename = this->gtPluse_ut_data_folder_ + "noisePrewhitening_DataScan_imag";
//    gt_io.importArray(imag_data, filename);
//    imag_data.print(std::cout);
//
//    ho3DArray<GT_Complex8> data;
//    real_imag_to_complex<GT_Complex8>(real_data, imag_data, data);
//
//    GADGET_MSG(data(42, 12, 10));
//
//    // apply the noise matrix
//    GADGET_START_TIMING(timer, "performNoisePrewhitening");
//    utilCplx.performNoisePrewhitening(data, prewhiteningMatrix);
//    GADGET_STOP_TIMING(timer);
//    GADGET_MSG(data(42, 12, 10));
//    EXPECT_LE(std::abs(data(42, 12, 10)-TValueType(-0.068069, -0.185625)), 1e-6);
//}
//
//
//TYPED_TEST(gtPlus_IO_Test, IOTest)
//{
//    typedef GT_Complex8 T;
//
//    gtPlusIOAnalyze gt_io;
//
//    hoNDArray<float> real_Im;
//    std::string filename = this->gtPluse_ut_data_folder_ + "KSpaceBinning_IncomingKSpace_real";
//    gt_io.importArray(real_Im, filename);
//    real_Im.print(std::cout);
//
//    hoNDArray<float> imag_Im;
//    filename = this->gtPluse_ut_data_folder_ + "KSpaceBinning_IncomingKSpace_imag";
//    gt_io.importArray(imag_Im, filename);
//    imag_Im.print(std::cout);
//
//    filename = this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_real2";
//    gt_io.exportArray(real_Im, filename);
//
//    filename = this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_imag2";
//    gt_io.exportArray(imag_Im, filename);
//
//    boost::shared_ptr< hoNDArray<GT_Complex8> > tmp = real_imag_to_complex<GT_Complex8>(&real_Im, &imag_Im);
//
//    unsigned int RO = tmp->get_size(0);
//    unsigned int E1 = tmp->get_size(1);
//    unsigned int CHA = tmp->get_size(2);
//    unsigned int PHS = tmp->get_size(3);
//
//    hoNDArray<GT_Complex8> kspace(RO, E1, CHA, PHS, tmp->begin());
//
//    float nrm2;
//    Gadgetron::norm2(kspace, nrm2);
//    GADGET_MSG("nrm2 = " << nrm2);
//
//    gtPlusISMRMRDReconUtil<GT_Complex8> util;
//    gtPlusISMRMRDReconUtilComplex<GT_Complex8> utilCplx;
//
//    // sum of square
//    hoNDArray<GT_Complex8> complexIm, sosIm;
//
//    GadgetronTimer timer(false);
//    timer.start("ifft2c");
//    hoNDFFT<float>::instance()->ifft2c(kspace, complexIm);
//    timer.stop();
//
//    timer.start("sumOfSquare");
//    utilCplx.sumOfSquare(complexIm, sosIm);
//    timer.stop();
//
//    hoNDArray<float> magSoS;
//    timer.start("absolute");
//    Gadgetron::absolute(sosIm, magSoS);
//    timer.stop();
//
//    filename = this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_SoS";
//    gt_io.exportArray(magSoS, filename);
//
//    // coil map estimation
//
//    hoNDArray<GT_Complex8> meanKSpace;
//    sumOverLastDimension(kspace, meanKSpace);
//
//    filename = this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_mean";
//    gt_io.export3DArrayComplex(meanKSpace, filename);
//
//    Gadgetron::norm2(meanKSpace, nrm2);
//    GADGET_MSG("nrm2 = " << nrm2);
//
//    hoNDArray<GT_Complex8> meanIm;
//    hoNDFFT<float>::instance()->ifft2c(meanKSpace, meanIm);
//    Gadgetron::norm2(meanIm, nrm2);
//    GADGET_MSG("nrm2 = " << nrm2);
//
//    filename = this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_meanIm";
//    gt_io.export3DArrayComplex(meanIm, filename);
//
//    hoNDArray<GT_Complex8> coilMap;
//    timer.start("coilMap2DNIH");
//    utilCplx.coilMap2DNIH(meanIm, coilMap, ISMRMRD_SOUHEIL, 7, 3);
//    timer.stop();
//
//    filename = this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_meanIm_coilMap";
//    gt_io.export3DArrayComplex(coilMap, filename);
//
//    hoNDArray<GT_Complex8> combined;
//    timer.start("coilCombine");
//    utilCplx.coilCombine(meanIm, coilMap, combined);
//    timer.stop();
//
//    gt_io.export3DArrayComplex(combined, this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_meanIm_coilMap_combined");
//
//    // KLT
//    hoMatrix<T> coeff, eigenValues;
//    timer.start("computeKLTCoeff");
//    util.computeKLTCoeff(meanKSpace, coeff, eigenValues);
//    timer.stop();
//    eigenValues.print(std::cout);
//
//    double thres = 0.001;
//    timer.start("computeKLCoilCompressionCoeff, thres");
//    util.computeKLCoilCompressionCoeff(meanKSpace, thres, coeff, eigenValues);
//    timer.stop();
//    eigenValues.print(std::cout);
//
//    hoNDArray<T> dataEigen;
//    int numOfModeKept = 20;
//    util.computeKLCoilCompression(meanKSpace, numOfModeKept, coeff, eigenValues, dataEigen);
//    Gadgetron::norm2(dataEigen, nrm2);
//    GADGET_MSG("nrm2 = " << nrm2);
//
//    hoNDFFT<float>::instance()->ifft2c(dataEigen, meanIm);
//    gt_io.export3DArrayComplex(meanIm, this->gtPluse_ut_res_folder_ + "KSpaceBinning_IncomingKSpace_meanIm_dataEigen");
//}
