
#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#include "Gadget.h"
#include "ismrmrd/ismrmrd.h"
#include "complext.h"

#include <gtest/gtest.h>

#include "hoNDArray_utils.h"

#include "gtPlusIOAnalyze.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "gtPlusISMRMRDReconWorker2DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker2DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorker3DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian2DT.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian3DT.h"

#include "GadgetronTimer.h"

#include <boost/thread/mutex.hpp>

#ifdef max
#undef max
#endif // max

using namespace Gadgetron;
using namespace Gadgetron::gtPlus;
using testing::Types;

template <typename T> class gtPlus_grappa_Test : public ::testing::Test 
{
protected:
    virtual void SetUp()
    {
        GDEBUG_STREAM("=============================================================================================");
        gt_ut_folder_ = std::string(::getenv("GADGETRON_UNITTEST_DIRECTORY"));
        GDEBUG_STREAM("=============================================================================================");
        GDEBUG_STREAM("Unit Test");
        gt_ut_data_folder_ = gt_ut_folder_ + "/data/";
        gt_ut_res_folder_ = gt_ut_folder_ + "/result/";
        GDEBUG_STREAM("gt_ut_data_folder_ is " << gt_ut_data_folder_);
        GDEBUG_STREAM("gt_ut_res_folder_ is " << gt_ut_res_folder_);

        timer_.set_timing_in_destruction(false);

#ifdef WIN32
    #ifdef USE_OMP
        /// lock the threads
        #pragma omp parallel default(shared)
        {
            int tid = omp_get_thread_num();
            // GDEBUG_STREAM(tid << std::endl);
            DWORD_PTR mask = (1 << tid);
            SetThreadAffinityMask( GetCurrentThread(), mask );
        }
    #endif // USE_OMP
#endif // WIN32
    }

    std::string gt_ut_folder_;
    std::string gt_ut_data_folder_;
    std::string gt_ut_res_folder_;

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

TYPED_TEST_CASE(gtPlus_grappa_Test, cpfloatImplementations);

TYPED_TEST(gtPlus_grappa_Test, reconWorker2DTGRAPPA_SNRUnit)
{
    typedef std::complex<float> T;

    gtPlusIOAnalyze gt_io;

    float v;

    // image data
    hoNDArray<std::complex<float> > data;
    gt_io.importArrayComplex(data, this->gt_ut_data_folder_ + "StdandardDataR2_Kspace_real", 
        this->gt_ut_data_folder_ + "StdandardDataR2_Kspace_imag");
    data.print(std::cout);

    unsigned long long RO = data.get_size(0);
    unsigned long long E1 = data.get_size(1);
    unsigned long long CHA = data.get_size(2);
    unsigned long long PHS = data.get_size(3);

    unsigned long long reconE1 = 144;

    // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    unsigned long long SLC = 1;
    unsigned long long E2 = 1;
    unsigned long long CON = 1;
    unsigned long long REP = 1;
    unsigned long long SET = 1;
    unsigned long long SEG = 1;

    hoNDArray<std::complex<float> > kspace(RO, E1, CHA, SLC, E2, CON, PHS);
    memcpy(kspace.begin(), data.begin(), data.get_number_of_bytes());

    Gadgetron::norm2(kspace, v);
    GDEBUG_STREAM("kspace = " << v);

    // ref
    hoNDArray<T> refTmp;
    gt_io.importArrayComplex(refTmp, this->gt_ut_data_folder_ + "StdandardDataR2_Ref_real", 
        this->gt_ut_data_folder_ + "StdandardDataR2_Ref_imag");

    hoNDArray<T> ref(refTmp.get_size(0), refTmp.get_size(1), refTmp.get_size(2), SLC, E2, CON, PHS);
    memcpy(ref.begin(), refTmp.begin(), refTmp.get_number_of_bytes());
    ref.print(std::cout);

    // noise
    hoNDArray<T> noise;
    gt_io.importArrayComplex(noise, this->gt_ut_data_folder_ + "StdandardDataR2_Noise_real", 
        this->gt_ut_data_folder_ + "StdandardDataR2_Noise_imag");
    noise.print(std::cout);

    // call the recon
    typedef std::complex<float> ValueType;
    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder2DT<ValueType> WorkOrderType;
    typedef std::pair<Gadgetron::ISMRMRDDIM, unsigned long long> DimensionRecordType;

    WorkOrderType* workOrder = new WorkOrderType;

    boost::shared_ptr< std::vector<size_t> > dims = kspace.get_dimensions();

    GDEBUG_STREAM("[Ro E1 Cha Slice E2 Con Phase Rep Set Seg] = [" 
        << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4] 
        << " " << (*dims)[5] << " " << (*dims)[6] << " " << 1 << " " << 1 << " " << 1 << "]");

    std::vector<size_t> dimensions_ = *dims;

        // work flow
    Gadgetron::gtPlus::gtPlusISMRMRDReconWorkFlowCartesian2DT<ValueType> workflow_;

    // worker
    Gadgetron::gtPlus::gtPlusReconWorker2DTGRAPPA<ValueType> worker_grappa_;

    // parameters
    Gadgetron::ISMRMRDDIM dim_4th_ = DIM_Phase;
    Gadgetron::ISMRMRDDIM dim_5th_ = DIM_Slice;
    Gadgetron::ISMRMRDDIM workOrder_ShareDim_ = DIM_NONE;

    bool interleaved_same_combinationcoeff_allS_ = false;
    int interleaved_whichS_combinationcoeff_ = 0;

    bool embedded_averageall_ref_ = false;
    bool embedded_fullres_coilmap_ = true;
    bool embedded_same_combinationcoeff_allS_ = false;
    int embedded_whichS_combinationcoeff_ = 0;
    bool embedded_ref_fillback_ = true;

    bool separate_averageall_ref_ = true;
    bool separate_fullres_coilmap_ = false;
    bool separate_same_combinationcoeff_allS_ = false;
    int separate_whichS_combinationcoeff_ = 0;

    bool same_coil_compression_coeff_allS_ = true;
    bool downstream_coil_compression_ = true;
    double coil_compression_thres_ = 1e-3;
    int coil_compression_num_modesKept_ = -1;

    unsigned long long csm_kSize_ = 7;
    unsigned long long csm_powermethod_num_ = 3;

    Gadgetron::ISMRMRDALGO recon_algorithm_ = ISMRMRD_GRAPPA;
    bool recon_kspace_needed_ = true;

    unsigned long long grappa_kSize_RO_ = 5;
    unsigned long long grappa_kSize_E1_ = 4;
    unsigned long long grappa_kSize_E2_ = 4;
    double grappa_reg_lamda_ = 1e-4;

    // recon
    workflow_.setDataArray(kspace);
    workflow_.setRefArray(ref);
    workflow_.noise_ = &noise;

    workflow_.noiseBW_ = 130;
    workflow_.overSamplingRatioRO_ = 2.0;
    workflow_.ADCSamplingTimeinSecond_ = 7800/1e9;

    // for this ut data, the oversampling removal and noise prewhitening on ref are not needed
    workflow_.ref_remove_oversampling_RO_ = false;
    workflow_.ref_apply_noisePreWhitening_ = false;

    workflow_.reconSizeRO_ = RO/2;
    workflow_.reconSizeE1_ = reconE1;
    workflow_.reconSizeE2_ = 1;
    // workflow_.dataDimStartingIndexes_ = workOrder->dataDimStartingIndexes_;
    workflow_.dim4th_ = dim_4th_;
    workflow_.dim5th_ = dim_5th_;

    workOrder->CalibMode_ = ISMRMRD_separate;
    workOrder->acceFactorE1_ = 2;
    workOrder->acceFactorE2_ = 1;

    workOrder->downstream_coil_compression_ = downstream_coil_compression_;
    workOrder->coil_compression_thres_ = coil_compression_thres_;
    workOrder->coil_compression_num_modesKept_ = coil_compression_num_modesKept_;
    workOrder->csm_kSize_ = csm_kSize_;
    workOrder->csm_powermethod_num_ = csm_powermethod_num_;
    workOrder->grappa_kSize_RO_ = grappa_kSize_RO_;
    workOrder->grappa_kSize_E1_ = grappa_kSize_E1_;
    workOrder->grappa_kSize_E2_ = grappa_kSize_E2_;
    workOrder->grappa_reg_lamda_ = grappa_reg_lamda_;
    workOrder->recon_kspace_needed_ = recon_kspace_needed_;

    if ( coil_compression_thres_>0 || coil_compression_num_modesKept_>0 )
    {
        workOrder->coil_compression_ = true;
    }
    else
    {
        workOrder->coil_compression_ = false;
    }

    workOrder->same_coil_compression_coeff_allS_ = same_coil_compression_coeff_allS_;
    workOrder->embedded_averageall_ref_ = embedded_averageall_ref_;
    workOrder->embedded_fullres_coilmap_ = embedded_fullres_coilmap_;
    workOrder->embedded_same_combinationcoeff_allS_ = embedded_same_combinationcoeff_allS_;
    workOrder->embedded_whichS_combinationcoeff_ = embedded_whichS_combinationcoeff_;
    workOrder->embedded_ref_fillback_ = embedded_ref_fillback_;
    workOrder->separate_averageall_ref_ = separate_averageall_ref_;
    workOrder->separate_fullres_coilmap_ = separate_fullres_coilmap_;
    workOrder->separate_same_combinationcoeff_allS_ = separate_same_combinationcoeff_allS_;
    workOrder->separate_whichS_combinationcoeff_ = separate_whichS_combinationcoeff_;
    workOrder->interleaved_same_combinationcoeff_allS_ = interleaved_same_combinationcoeff_allS_;
    workOrder->interleaved_whichS_combinationcoeff_ = interleaved_whichS_combinationcoeff_;

    worker_grappa_.performTiming_ = true;
    worker_grappa_.debugFolder_ = this->gt_ut_res_folder_;

    workflow_.debugFolder_ = this->gt_ut_res_folder_;

    workflow_.worker_ = &worker_grappa_;
    workflow_.workOrder_ = workOrder;

    workflow_.preProcessing();
    workflow_.recon();
    workflow_.postProcessing();

    gt_io.exportArrayComplex(workflow_.res_, this->gt_ut_res_folder_+"StdandardDataR2_res");

    workflow_.res_.squeeze();
    gt_io.export3DArrayComplex(workflow_.res_, this->gt_ut_res_folder_+"StdandardDataR2_res_squeezed");

    hoNDArray<T> std;
    bool NMinusOne = true;
    stdOver3rdDimension(workflow_.res_, std, NMinusOne);
    gt_io.export2DArrayComplex(std, this->gt_ut_res_folder_+"StdandardDataR2_res_squeezed_std");
}

TYPED_TEST(gtPlus_grappa_Test, reconWorker2DTGRAPPA)
{
    typedef std::complex<float> T;

    gtPlusIOAnalyze gt_io;

    float v;

    // image data
    hoNDArray<float> real_data;
    std::string filename = this->gt_ut_data_folder_ + "underSampledKSpace_real";
    gt_io.importArray(real_data, filename);
    real_data.print(std::cout);

    hoNDArray<float> imag_data;
    filename = this->gt_ut_data_folder_ + "underSampledKSpace_imag";
    gt_io.importArray(imag_data, filename);
    imag_data.print(std::cout);

    hoNDArray<std::complex<float> > tmp;
    Gadgetron::real_imag_to_complex<std::complex<float> >(real_data, imag_data, tmp);

    unsigned long long RO = tmp.get_size(0);
    unsigned long long E1 = tmp.get_size(1);
    unsigned long long CHA = tmp.get_size(2);
    unsigned long long PHS = tmp.get_size(3);

    unsigned long long reconE1 = 120;

    // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    unsigned long long SLC = 1;
    unsigned long long E2 = 1;
    unsigned long long CON = 1;
    unsigned long long REP = 1;
    unsigned long long SET = 1;
    unsigned long long SEG = 1;

    hoNDArray<std::complex<float> > kspace(RO, E1, CHA, SLC, E2, CON, PHS, tmp.begin());

    Gadgetron::norm2(kspace, v);
    GDEBUG_STREAM("kspace = " << v);

    // ref
    hoNDArray<float> real_ref;
    filename = this->gt_ut_data_folder_ + "ref_real";
    gt_io.importArray(real_ref, filename);
    real_ref.print(std::cout);

    hoNDArray<float> imag_ref;
    filename = this->gt_ut_data_folder_ + "ref_imag";
    gt_io.importArray(imag_ref, filename);
    imag_ref.print(std::cout);

    hoNDArray<T> ref;
    real_imag_to_complex<std::complex<float> >(real_ref, imag_ref, ref);

    Gadgetron::norm2(ref, v);
    GDEBUG_STREAM("ref = " << v);

    // call the recon
    typedef std::complex<float> ValueType;
    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder2DT<ValueType> WorkOrderType;
    typedef std::pair<Gadgetron::ISMRMRDDIM, unsigned long long> DimensionRecordType;

    WorkOrderType* workOrder = new WorkOrderType;

    workOrder->data_ = kspace;
    workOrder->ref_ = ref;

    boost::shared_ptr< std::vector<size_t> > dims = workOrder->data_.get_dimensions();

    GDEBUG_STREAM("[Ro E1 Cha Slice E2 Con Phase Rep Set Seg] = [" 
        << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4] 
        << " " << (*dims)[5] << " " << (*dims)[6] << " " << 1 << " " << 1 << " " << 1 << "]");

    std::vector<size_t> dimensions_ = *dims;

        // work flow
    Gadgetron::gtPlus::gtPlusISMRMRDReconWorkFlowCartesian2DT<ValueType> workflow_;

    // worker
    Gadgetron::gtPlus::gtPlusReconWorker2DTGRAPPA<ValueType> worker_grappa_;

    // parameters
    Gadgetron::ISMRMRDDIM dim_4th_ = DIM_Phase;
    Gadgetron::ISMRMRDDIM dim_5th_ = DIM_Slice;
    Gadgetron::ISMRMRDDIM workOrder_ShareDim_ = DIM_NONE;

    bool interleaved_same_combinationcoeff_allS_ = false;
    int interleaved_whichS_combinationcoeff_ = 0;

    bool embedded_averageall_ref_ = false;
    bool embedded_fullres_coilmap_ = true;
    bool embedded_same_combinationcoeff_allS_ = false;
    int embedded_whichS_combinationcoeff_ = 0;
    bool embedded_ref_fillback_ = true;

    bool separate_averageall_ref_ = false;
    bool separate_fullres_coilmap_ = true;
    bool separate_same_combinationcoeff_allS_ = false;
    int separate_whichS_combinationcoeff_ = 0;

    bool same_coil_compression_coeff_allS_ = true;
    bool downstream_coil_compression_ = true;
    double coil_compression_thres_ = 1e-3;
    int coil_compression_num_modesKept_ = -1;

    unsigned long long csm_kSize_ = 7;
    unsigned long long csm_powermethod_num_ = 3;

    Gadgetron::ISMRMRDALGO recon_algorithm_ = ISMRMRD_GRAPPA;
    bool recon_kspace_needed_ = true;

    unsigned long long grappa_kSize_RO_ = 5;
    unsigned long long grappa_kSize_E1_ = 4;
    unsigned long long grappa_kSize_E2_ = 4;
    double grappa_reg_lamda_ = 1e-4;

    // recon
    workflow_.setDataArray(kspace);
    workflow_.setRefArray(ref);

    Gadgetron::norm2(workOrder->data_, v); GDEBUG_STREAM("workOrder->data_ = " << v);
    Gadgetron::norm2(workOrder->ref_, v); GDEBUG_STREAM("workOrder->ref_ = " << v);

    workflow_.reconSizeRO_ = RO;
    workflow_.reconSizeE1_ = reconE1;
    workflow_.reconSizeE2_ = 1;
    // workflow_.dataDimStartingIndexes_ = workOrder->dataDimStartingIndexes_;
    workflow_.dim4th_ = dim_4th_;
    workflow_.dim5th_ = dim_5th_;

    workOrder->CalibMode_ = ISMRMRD_separate;
    workOrder->start_RO_ = 34;
    workOrder->end_RO_ = (int)RO-1;
    workOrder->acceFactorE1_ = 4;
    workOrder->acceFactorE2_ = 1;

    workOrder->downstream_coil_compression_ = downstream_coil_compression_;
    workOrder->coil_compression_thres_ = coil_compression_thres_;
    workOrder->coil_compression_num_modesKept_ = coil_compression_num_modesKept_;
    workOrder->csm_kSize_ = csm_kSize_;
    workOrder->csm_powermethod_num_ = csm_powermethod_num_;
    workOrder->grappa_kSize_RO_ = grappa_kSize_RO_;
    workOrder->grappa_kSize_E1_ = grappa_kSize_E1_;
    workOrder->grappa_kSize_E2_ = grappa_kSize_E2_;
    workOrder->grappa_reg_lamda_ = grappa_reg_lamda_;
    workOrder->recon_kspace_needed_ = recon_kspace_needed_;

    if ( coil_compression_thres_>0 || coil_compression_num_modesKept_>0 )
    {
        workOrder->coil_compression_ = true;
    }
    else
    {
        workOrder->coil_compression_ = false;
    }

    workOrder->same_coil_compression_coeff_allS_ = same_coil_compression_coeff_allS_;
    workOrder->embedded_averageall_ref_ = embedded_averageall_ref_;
    workOrder->embedded_fullres_coilmap_ = embedded_fullres_coilmap_;
    workOrder->embedded_same_combinationcoeff_allS_ = embedded_same_combinationcoeff_allS_;
    workOrder->embedded_whichS_combinationcoeff_ = embedded_whichS_combinationcoeff_;
    workOrder->embedded_ref_fillback_ = embedded_ref_fillback_;
    workOrder->separate_averageall_ref_ = separate_averageall_ref_;
    workOrder->separate_fullres_coilmap_ = separate_fullres_coilmap_;
    workOrder->separate_same_combinationcoeff_allS_ = separate_same_combinationcoeff_allS_;
    workOrder->separate_whichS_combinationcoeff_ = separate_whichS_combinationcoeff_;
    workOrder->interleaved_same_combinationcoeff_allS_ = interleaved_same_combinationcoeff_allS_;
    workOrder->interleaved_whichS_combinationcoeff_ = interleaved_whichS_combinationcoeff_;

    worker_grappa_.performTiming_ = true;
    worker_grappa_.debugFolder_ = this->gt_ut_res_folder_;

    workflow_.debugFolder_ = this->gt_ut_res_folder_;
    workflow_.worker_ = &worker_grappa_;
    workflow_.workOrder_ = workOrder;

    gt_io.exportArrayComplex(workflow_.workOrder_->ref_, this->gt_ut_res_folder_+"ref");

    workflow_.preProcessing();
    workflow_.recon();
    workflow_.postProcessing();

    gt_io.exportArrayComplex(workflow_.res_, this->gt_ut_res_folder_+"grappa2D_gtPlus_res");
}

TYPED_TEST(gtPlus_grappa_Test, grappa2D)
{
    typedef std::complex<float> T;

    gtPlusIOAnalyze gt_io;

    float v;

    // image data
    hoNDArray<float> real_data;
    std::string filename = this->gt_ut_data_folder_ + "underSampledKSpace_real";
    gt_io.importArray(real_data, filename);
    real_data.print(std::cout);

    hoNDArray<float> imag_data;
    filename = this->gt_ut_data_folder_ + "underSampledKSpace_imag";
    gt_io.importArray(imag_data, filename);
    imag_data.print(std::cout);

    hoNDArray<std::complex<float> > tmp;
    Gadgetron::real_imag_to_complex<std::complex<float> >(real_data, imag_data, tmp);

    unsigned long long RO = tmp.get_size(0);
    unsigned long long E1 = tmp.get_size(1);
    unsigned long long CHA = tmp.get_size(2);
    unsigned long long PHS = tmp.get_size(3);

    hoNDArray<std::complex<float> > kspace(RO, E1, CHA, PHS, tmp.begin());

    // ref
    hoNDArray<float> real_ref;
    filename = this->gt_ut_data_folder_ + "ref_real";
    gt_io.importArray(real_ref, filename);
    real_ref.print(std::cout);

    hoNDArray<float> imag_ref;
    filename = this->gt_ut_data_folder_ + "ref_imag";
    gt_io.importArray(imag_ref, filename);
    imag_ref.print(std::cout);

    hoNDArray<T> ref;
    real_imag_to_complex<std::complex<float> >(real_ref, imag_ref, ref);

    Gadgetron::norm2(ref, v);
    GDEBUG_STREAM("ref = " << v);

    // recon
    gtPlusISMRMRDReconUtil<std::complex<float> > util;
    gtPlusISMRMRDReconUtilComplex<std::complex<float> > utilCplx;

    // sum of square
    hoNDArray<std::complex<float> > complexIm, sosIm;

    GadgetronTimer timer(false);
    timer.start("ifft2c");
    hoNDFFT<float>::instance()->ifft2c(kspace, complexIm);
    timer.stop();

    timer.start("sumOfSquare");
    utilCplx.sumOfSquare(complexIm, sosIm);
    timer.stop();

    hoNDArray<float> magSoS;
    timer.start("absolute");
    Gadgetron::abs(sosIm, magSoS);
    timer.stop();

    filename = this->gt_ut_res_folder_ + "SoS";
    gt_io.exportArray(magSoS, filename);

    // coil map estimation
    hoNDFFT<float>::instance()->ifft2c(ref, complexIm);

    filename = this->gt_ut_res_folder_ + "complexIm";
    gt_io.export3DArrayComplex(complexIm, filename);

    hoNDArray<std::complex<float> > coilMap;
    timer.start("coilMap2DNIH");
    utilCplx.coilMap2DNIH(complexIm, coilMap, ISMRMRD_SOUHEIL, 7, 3, 3, true);
    timer.stop();

    filename = this->gt_ut_res_folder_ + "coilMap";
    gt_io.export3DArrayComplex(coilMap, filename);

    // grappa kernel estimation
    gtPlusReconWorker2DTGRAPPA<T> grappa;

    unsigned long long kRO = 5;
    unsigned long long kNE1 = 4;
    unsigned long long srcCHA = CHA;
    unsigned long long dstCHA = CHA;

    double grappa_reg_lamda_ = 1e-4;

    ho3DArray<T> acsSrc(RO, E1, CHA, const_cast<T*>(ref.begin()));
    ho3DArray<T> acsDst(RO, E1, CHA, const_cast<T*>(ref.begin()));

    Gadgetron::norm2(acsSrc, v);
    GDEBUG_STREAM("acsSrc = " << v);

    int accelFactor = 4;
    bool fitItself = true;

    ho4DArray<T> convKer;
    timer.start("grappa2d_calib_convolution_kernel");
    Gadgetron::grappa2d_calib_convolution_kernel(acsSrc, acsDst, accelFactor, grappa_reg_lamda_, kRO, kNE1, convKer);
    timer.stop();

    Gadgetron::norm2(convKer, v);
    GDEBUG_STREAM("convKer = " << v);
    gt_io.exportArrayComplex(convKer, this->gt_ut_res_folder_ + "convKer");

    ho4DArray<T> kIm(RO, E1, srcCHA, dstCHA);
    timer.start("grappa2d_image_domain_kernel");
    Gadgetron::grappa2d_image_domain_kernel(convKer, RO, E1, kIm);
    timer.stop();
    gt_io.exportArrayComplex(kIm, this->gt_ut_res_folder_ + "kIm");

    Gadgetron::norm2(kIm, v);
    GDEBUG_STREAM("kIm = " << v);

    ho3DArray<T> unmixC(RO, E1, srcCHA);
    ho2DArray<float> gFactor(RO, E1);

    ho3DArray<T> coilMap2(RO, E1, dstCHA, coilMap.begin());

    Gadgetron::norm2(coilMap2, v);
    GDEBUG_STREAM("coilMap2 = " << v);

    Gadgetron::grappa2d_unmixing_coeff(kIm, coilMap2, accelFactor, unmixC, gFactor);

    Gadgetron::norm2(unmixC, v);
    GDEBUG_STREAM("unmixC = " << v);

    gt_io.export3DArrayComplex(unmixC, this->gt_ut_res_folder_ + "unmixC");
    gt_io.export2DArray(gFactor, this->gt_ut_res_folder_ + "gFactor");

    // unwarpping
    hoNDArray<T> res;
    grappa.applyImageDomainKernel(kspace, kIm, res);
    gt_io.export3DArrayComplex(res, this->gt_ut_res_folder_ + "grappa2D_res");

    Gadgetron::apply_unmix_coeff_kspace(kspace, unmixC, res);
    gt_io.export2DArrayComplex(res, this->gt_ut_res_folder_ + "res_unmixC");
}
