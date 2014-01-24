
#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#include "Gadget.h"
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
#include "gtPlusISMRMRDReconWorker2DTL1SPIRITNCG.h"
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
#include "gtPlusLSQRSolver.h"
#include "gtPlusNCGSolver.h"
#include "gtPlusWavelet2DOperator.h"
#include "gtPlusWavelet3DOperator.h"
#include "gtPlusWaveletNoNullSpace2DOperator.h"
#include "gtPlusWaveletNoNullSpace3DOperator.h"
#include "gtPlusDataFidelityOperator.h"
#include "gtPlusISMRMRDReconWorkFlowCartesian3DT.h"
#include "gtPlusISMRMRDReconWorker3DTGRAPPA.h"
#include "gtPlusISMRMRDReconWorker3DTNoAcceleration.h"
#include "gtPlusISMRMRDReconWorker3DTSPIRIT.h"
#include "gtPlusISMRMRDReconWorker3DTL1SPIRITNCG.h"
#include "gtPlusMemoryManager.h"

#include "GadgetronTimer.h"

#include <boost/thread/mutex.hpp>

#ifdef max
#undef max
#endif // max

using namespace Gadgetron;
using namespace Gadgetron::gtPlus;
using testing::Types;

template <typename T> class gtPlus_spirit_Test : public ::testing::Test 
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

TYPED_TEST_CASE(gtPlus_spirit_Test, cpfloatImplementations);

TYPED_TEST(gtPlus_spirit_Test, reconWorker2DTSPIRIT)
{
    typedef GT_Complex8 T;

    gtPlusIOAnalyze gt_io;

    float v;

    // image data
    hoNDArray<float> real_data;
    std::string filename = this->gtPluse_ut_data_folder_ + "underSampledKSpace_real";
    gt_io.importArray(real_data, filename);
    real_data.print(std::cout);

    hoNDArray<float> imag_data;
    filename = this->gtPluse_ut_data_folder_ + "underSampledKSpace_imag";
    gt_io.importArray(imag_data, filename);
    imag_data.print(std::cout);

    boost::shared_ptr< hoNDArray<GT_Complex8> > tmp = real_imag_to_complex<GT_Complex8>(&real_data, &imag_data);

    unsigned long long RO = tmp->get_size(0);
    unsigned long long E1 = tmp->get_size(1);
    unsigned long long CHA = tmp->get_size(2);
    unsigned long long PHS = tmp->get_size(3);

    unsigned long long reconE1 = 120;

    // [Ro E1 Cha Slice E2 Con Phase Rep Set Seg]
    unsigned long long SLC = 1;
    unsigned long long E2 = 1;
    unsigned long long CON = 1;
    unsigned long long REP = 1;
    unsigned long long SET = 1;
    unsigned long long SEG = 1;

    hoNDArray<GT_Complex8> kspace(RO, E1, CHA, SLC, E2, CON, PHS, tmp->begin());

    Gadgetron::norm2(kspace, v);
    GADGET_MSG("kspace = " << v);

    // ref
    hoNDArray<float> real_ref;
    filename = this->gtPluse_ut_data_folder_ + "ref_real";
    gt_io.importArray(real_ref, filename);
    real_ref.print(std::cout);

    hoNDArray<float> imag_ref;
    filename = this->gtPluse_ut_data_folder_ + "ref_imag";
    gt_io.importArray(imag_ref, filename);
    imag_ref.print(std::cout);

    hoNDArray<T> ref;
    real_imag_to_complex<GT_Complex8>(real_ref, imag_ref, ref);

    Gadgetron::norm2(ref, v);
    GADGET_MSG("ref = " << v);

    // call the recon
    typedef std::complex<float> ValueType;
    typedef Gadgetron::gtPlus::gtPlusReconWorkOrder2DT<ValueType> WorkOrderType;
    typedef std::pair<Gadgetron::gtPlus::ISMRMRDDIM, unsigned long long> DimensionRecordType;

    WorkOrderType* workOrder = new WorkOrderType;

    workOrder->data_ = kspace;
    workOrder->ref_ = ref;

    boost::shared_ptr< std::vector<size_t> > dims = workOrder->data_.get_dimensions();

    GADGET_MSG("[Ro E1 Cha Slice E2 Con Phase Rep Set Seg] = [" 
        << (*dims)[0] << " " << (*dims)[1] << " " << (*dims)[2] << " " << (*dims)[3] << " " << (*dims)[4] 
        << " " << (*dims)[5] << " " << (*dims)[6] << " " << 1 << " " << 1 << " " << 1 << "]");

    std::vector<size_t> dimensions_ = *dims;

        // work flow
    Gadgetron::gtPlus::gtPlusISMRMRDReconWorkFlowCartesian2DT<ValueType> workflow_;

    // worker
    Gadgetron::gtPlus::gtPlusReconWorker2DTSPIRIT<ValueType> worker_spirit_;

    // parameters
    Gadgetron::gtPlus::ISMRMRDDIM dim_4th_ = DIM_Phase;
    Gadgetron::gtPlus::ISMRMRDDIM dim_5th_ = DIM_Slice;
    Gadgetron::gtPlus::ISMRMRDDIM workOrder_ShareDim_ = DIM_NONE;

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

    Gadgetron::gtPlus::ISMRMRDALGO recon_algorithm_ = ISMRMRD_SPIRIT;
    bool recon_kspace_needed_ = true;

    unsigned long long spirit_kSize_RO_ = 5;
    unsigned long long spirit_kSize_E1_ = 5;
    unsigned long long spirit_kSize_E2_ = 5;

    double spirit_reg_lamda_ = 0.005;
    unsigned long long spirit_iter_max_ = 100;
    double spirit_iter_thres_ = 1e-5;

    // recon
    workflow_.setDataArray(kspace);
    workflow_.setRefArray(ref);

    Gadgetron::norm2(workOrder->data_, v); GADGET_MSG("workOrder->data_ = " << v);
    Gadgetron::norm2(workOrder->ref_, v); GADGET_MSG("workOrder->ref_ = " << v);

    workflow_.reconSizeRO_ = RO;
    workflow_.reconSizeE1_ = reconE1;
    workflow_.reconSizeE2_ = 1;
    // workflow_.dataDimStartingIndexes_ = workOrder->dataDimStartingIndexes_;
    workflow_.dim4th_ = dim_4th_;
    workflow_.dim5th_ = dim_5th_;

    workOrder->CalibMode_ = ISMRMRD_separate;
    workOrder->start_RO_ = 34;
    workOrder->end_RO_ = RO-1;
    workOrder->acceFactorE1_ = 4;
    workOrder->acceFactorE2_ = 1;

    workOrder->downstream_coil_compression_ = downstream_coil_compression_;
    workOrder->coil_compression_thres_ = coil_compression_thres_;
    workOrder->coil_compression_num_modesKept_ = coil_compression_num_modesKept_;
    workOrder->csm_kSize_ = csm_kSize_;
    workOrder->csm_powermethod_num_ = csm_powermethod_num_;;

    workOrder->recon_algorithm_ = recon_algorithm_;

    workOrder->spirit_kSize_RO_ = spirit_kSize_RO_;
    workOrder->spirit_kSize_E1_ = spirit_kSize_E1_;
    workOrder->spirit_kSize_E2_ = spirit_kSize_E2_;
    workOrder->spirit_reg_lamda_ = spirit_reg_lamda_;
    workOrder->spirit_iter_max_ = spirit_iter_max_;
    workOrder->spirit_iter_thres_ = spirit_iter_thres_;

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

    worker_spirit_.performTiming_ = true;
    worker_spirit_.debugFolder_ = this->gtPluse_ut_res_folder_;

    workflow_.debugFolder_ = this->gtPluse_ut_res_folder_;
    workflow_.worker_ = &worker_spirit_;
    workflow_.workOrder_ = workOrder;

    gt_io.exportArrayComplex(workflow_.workOrder_->ref_, this->gtPluse_ut_res_folder_+"ref");

    boost::shared_ptr<Gadgetron::gtPlus::gtPlusMemoryManager> mem_manager_(new Gadgetron::gtPlus::gtPlusMemoryManager(4, 640*1024*1024));
    worker_spirit_.gtPlus_mem_manager_ = mem_manager_;

    workflow_.preProcessing();
    workflow_.recon();
    workflow_.postProcessing();

    gt_io.exportArrayComplex(workflow_.res_, this->gtPluse_ut_res_folder_+"spirit2D_gtPlus_res");
}

TYPED_TEST(gtPlus_spirit_Test, testNCGSolver2DTSPIRIT_neuro_3by3)
{
    typedef GT_Complex8 T;

    gtPlusIOAnalyze gt_io;

    float v;

    // image data
    hoNDArray<float> real_data;
    std::string filename = this->gtPluse_ut_data_folder_ + "Job2DT_kspace_ID6_REAL";
    gt_io.importArray(real_data, filename);
    real_data.print(std::cout);

    hoNDArray<float> imag_data;
    filename = this->gtPluse_ut_data_folder_ + "Job2DT_kspace_ID6_IMAG";
    gt_io.importArray(imag_data, filename);
    imag_data.print(std::cout);

    boost::shared_ptr< hoNDArray<GT_Complex8> > kspace = real_imag_to_complex<GT_Complex8>(&real_data, &imag_data);

    hoNDArray<float> real_ker;
    filename = this->gtPluse_ut_data_folder_ + "Job2DT_ker_ID6_REAL";
    gt_io.importArray(real_ker, filename);
    real_ker.print(std::cout);

    hoNDArray<float> imag_ker;
    filename = this->gtPluse_ut_data_folder_ + "Job2DT_ker_ID6_IMAG";
    gt_io.importArray(imag_ker, filename);
    imag_ker.print(std::cout);

    boost::shared_ptr< hoNDArray<GT_Complex8> > ker = real_imag_to_complex<GT_Complex8>(&real_ker, &imag_ker);

    hoNDArray<float> real_kspaceLinear;
    filename = this->gtPluse_ut_data_folder_ + "Job2DT_kspaceLinear_ID6_REAL";
    gt_io.importArray(real_kspaceLinear, filename);
    real_kspaceLinear.print(std::cout);

    hoNDArray<float> imag_kspaceLinear;
    filename = this->gtPluse_ut_data_folder_ + "Job2DT_kspaceLinear_ID6_IMAG";
    gt_io.importArray(imag_kspaceLinear, filename);
    imag_kspaceLinear.print(std::cout);

    boost::shared_ptr< hoNDArray<GT_Complex8> > kspaceLinear = real_imag_to_complex<GT_Complex8>(&real_kspaceLinear, &imag_kspaceLinear);

    Gadgetron::gtPlus::gtPlusReconWorker3DTL1SPIRITNCG<GT_Complex8> worker_spirit_L1_ncg_;
    worker_spirit_L1_ncg_.performTiming_ = true;
    worker_spirit_L1_ncg_.debugFolder_ = this->gtPluse_ut_res_folder_;

    Gadgetron::gtPlus::gtPlusReconJob2DT< std::complex<float> > job;

    job.kspace = *kspace;
    job.ker = *ker;

    job.workOrder2DT.CalibMode_ = ISMRMRD_embedded;
    job.workOrder2DT.InterleaveDim_ = DIM_Phase;

    job.workOrder2DT.acceFactorE1_ = 3;
    job.workOrder2DT.acceFactorE2_ = 3;

    job.workOrder2DT.kSpaceCenterRO_ = 128;
    job.workOrder2DT.kSpaceCenterEncode1_ = 127;
    job.workOrder2DT.kSpaceCenterEncode2_ = 96;

    job.workOrder2DT.kSpaceMaxRO_ = 256;
    job.workOrder2DT.kSpaceMaxEncode1_ = 255;
    job.workOrder2DT.kSpaceMaxEncode2_ = 191;

    job.workOrder2DT.recon_algorithm_ = ISMRMRD_L1SPIRIT;
    job.workOrder2DT.recon_auto_parameters_ = false;

    job.workOrder2DT.spirit_kSize_RO_ = 7;
    job.workOrder2DT.spirit_kSize_E1_ = 7;
    job.workOrder2DT.spirit_kSize_E2_ = 5;

    job.workOrder2DT.spirit_reg_lamda_ = 0.01;
    job.workOrder2DT.spirit_calib_over_determine_ratio_ = 15;

    job.workOrder2DT.spirit_solve_symmetric_ = false;

    job.workOrder2DT.spirit_iter_max_ = 100;
    job.workOrder2DT.spirit_iter_thres_ = 0.005;
    job.workOrder2DT.spirit_print_iter_ = true;

    job.workOrder2DT.spirit_perform_linear_ = true;
    job.workOrder2DT.spirit_perform_nonlinear_ = true;

    job.workOrder2DT.spirit_parallel_imaging_lamda_ = 1;
    job.workOrder2DT.spirit_image_reg_lamda_ = 0.0025;
    job.workOrder2DT.spirit_data_fidelity_lamda_ = 0;

    job.workOrder2DT.spirit_ncg_iter_max_ = 10;
    job.workOrder2DT.spirit_ncg_iter_thres_ = 0.001;
    job.workOrder2DT.spirit_ncg_print_iter_ = true;
    job.workOrder2DT.spirit_ncg_scale_factor_ = 1;

    job.workOrder2DT.spirit_use_coil_sen_map_ = false;
    job.workOrder2DT.spirit_use_moco_enhancement_ = false;
    job.workOrder2DT.spirit_recon_moco_images_ = false;

    job.workOrder2DT.spirit_temporal_enhancement_ratio_ = 5;
    job.workOrder2DT.spirit_2D_scale_per_chunk_ = false;

    job.workOrder2DT.spirit_E2_enhancement_ratio_ = 1.0;
    job.workOrder2DT.spirit_3D_scale_per_chunk_ = false;

    bool succeed = true;
    GADGET_START_TIMING_CONDITION(this->timer_, "Recon 2DT job ... ", true);

    job.res = job.kspace;

    worker_spirit_L1_ncg_.performUnwarppingImplROPermuted(&(job.workOrder2DT), job.kspace, job.ker, *job.workOrder2DT.coilMap_, job.res);
    // worker_spirit_L1_ncg_.performUnwarppingImplROPermuted(&(job.workOrder2DT), job.kspace, job.ker, *job.workOrder2DT.coilMap_, *kspaceLinear, job.res);

    // succeed = worker_spirit_L1_ncg_.performUnwarppingImpl(job);

    GADGET_STOP_TIMING_CONDITION(this->timer_, true);

    gt_io.exportArrayComplex(job.res, this->gtPluse_ut_res_folder_+"NCGSolver2DTSPIRIT_neuro_3by3_res");
}
