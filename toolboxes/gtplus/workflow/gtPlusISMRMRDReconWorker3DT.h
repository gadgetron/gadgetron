/** \file   gtPlusISMRMRDReconWorker3DT.h
    \brief  Define the base class for the GtPlus worker for 3D or 3DT reconstruction cases

            Four different strategies were implemented for partial fourier or asymmetric echo acquisition, including:

            ISMRMRD_PF_ZEROFILLING          : only zero filling the unacquired k-space

            ISMRMRD_PF_ZEROFILLING_FILTER   : zero filling the unacquired k-space and apply a transition filter on the edges between
                                              acquired and unacquired regions

            ISMRMRD_PF_POCS                 : perform the iterative POCS reconstruction
                                              Magnetic Resonance Imaging: Physical Principles and Sequence Design. Page 296-297.
                                              E. Mark Haacke, Robert W. Brown, Michael R. Thompson, Ramesh Venkatesan. 
                                              Wiley-Liss, ISBN-10: 0471351288.

            ISMRMRD_PF_FENGHUANG            : perform a k-space convolution based partial fourier reconstruction. 

                                              Feng Huang, Wei Lin, and Yu Li. 
                                              Partial Fourier Reconstruction Through Data Fitting and Convolution in k-Space.
                                              Magnetic Resonance in Medicine, Vol 62, page 1261�1269, 2009.
    \author Hui Xue
*/

#pragma once

#include "gtPlusISMRMRDReconWorker.h"

#include "mri_core_kspace_filter.h"
#include "mri_core_partial_fourier.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusReconWorker3DT : public gtPlusReconWorker<T>
{
public:

    typedef gtPlusReconWorker<T> BaseClass;
    typedef gtPlusReconWorkOrder3DT<T> WorkOrderType;
    typedef typename BaseClass::value_type value_type;

    gtPlusReconWorker3DT() : BaseClass(), startE1_(0), endE1_(0), startE2_(0), endE2_(0) {}
    virtual ~gtPlusReconWorker3DT() {}

    virtual bool performRecon(gtPlusReconWorkOrder<T>* workOrder)
    {
        // check whether we have all-zeros input
        value_type v(1);
        Gadgetron::norm2(workOrder->data_, v);
        if ( v <= 0 )
        {
            GWARN_STREAM("gtPlusReconWorker2DT, performRecon(workOrder) : incoming data contains all-zeros ... ");

            boost::shared_ptr< std::vector<size_t> > dims = workOrder->data_.get_dimensions();
            (*dims)[3] = workOrder->num_channels_res_;
            workOrder->complexIm_.create(dims);
            Gadgetron::clear(workOrder->complexIm_);

            return true;
        }

        gtPlusReconWorkOrder3DT<T>* workOrder3DT = dynamic_cast<gtPlusReconWorkOrder3DT<T>*>(workOrder);
        if ( workOrder3DT == NULL ) return false;

        if ( workOrder3DT->recon_auto_parameters_ )
        {
            this->autoReconParameter(workOrder3DT);
            GDEBUG_STREAM("Gt Plus 3DT -- automatic paramter selection ---");
            if ( !this->debugFolder_.empty() ) { workOrder3DT->print(std::cout); }
        }

        return this->performRecon(workOrder3DT);
    }

    // the common functionalities are performed here for 3DT recon
    // compute the coil compression coefficients
    // prepare the ref data array
    virtual bool performRecon(gtPlusReconWorkOrder3DT<T>* workOrder3DT);

    virtual bool estimateCoilMap(WorkOrderType* workOrder3DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst);
    virtual bool performCalib(WorkOrderType* workOrder3DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst);
    virtual bool performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT);
    virtual bool performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, WorkOrderType* workOrder3DT, size_t usedN);

    virtual bool performUnwrapping(WorkOrderType* workOrder3DT, const hoNDArray<T>& data);

    // the partial fourier handling for the 3DT reconstruction
    // the computation is performed on the reconstructed full kspace
    virtual bool performPartialFourierHandling(WorkOrderType* workOrder3DT);

    // perform the kspace filter on ref data for coil map estimation
    virtual bool performRefFilter(gtPlusReconWorkOrder3DT<T>* workOrder3DT, 
                                        const hoNDArray<T>& ref, hoNDArray<T>& refFiltered, 
                                        int startRO, int endRO, int startE1, int endE1, int startE2, int endE2);

    // for interleave, compute mean ref
    // for embedded and separate, squeeze out the zero lines
    virtual bool prepRef(WorkOrderType* workOrder3DT, 
                        const hoNDArray<T>& ref, 
                        hoNDArray<T>& refRecon, 
                        hoNDArray<T>& refCoilMap, 
                        int startRO, int endRO, 
                        int startE1, int endE1, 
                        int startE2, int endE2, 
                        size_t dataE1, 
                        size_t dataE2);

    // implement reference data preparation
    virtual bool prepRefByAveragingCrossN(WorkOrderType* workOrder3DT, const hoNDArray<T>& ref, bool averageAllRef, int numOfModes, hoNDArray<T>& refRecon);

    // compute coil compression coefficients
    virtual bool coilCompression(WorkOrderType* workOrder3DT);

    // after unwrapping, for embedded and separate, the full res coil map may be estimated
    // for embedded, the ref may be filled back to fullkspace
    virtual bool afterUnwrapping(WorkOrderType* workOrder3DT);

    // whether to recon kspace, if true, the coil combination may not be performed, only the fullkspace is computed
    virtual bool computeKSpace(gtPlusReconWorkOrder3DT<T>* workOrder3DT) = 0;

    // ----------------------------------------------------
    // common functions for 3DT reconstruction
    // ----------------------------------------------------
    // image domain kernel with coil sensitivity
    // kerIm: [RO E1 E2 srcCHA dstCHA]
    // coilMap: [RO E1 E2 dstCHA]
    // unmixCoeff: [RO E1 E2 srcCHA]
    // gFactor: [RO E1 E2]
    bool unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor);

    // apply image domain kernel
    // kspace: [RO E1 E2 srcCHA ...]
    // complexIm : [RO E1 E2 dstCHA ...]
    bool applyImageDomainKernel(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm);
    // aliasedIm : [RO E1 E2 srcCHA ...]
    bool applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm);
    // for speed, a buffer can be provided
    bool applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& kerImBuffer, hoNDArray<T>& complexIm);

    // apply unmixCoeff
    // kspace: [RO E1 E2 srcCHA ...]
    // unmixCoeff : [RO E1 E2 srcCHA]
    // complexIm : [RO E1 E2 ...]
    bool applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);
    // aliasedIm : [RO E1 E2 srcCHA ...]
    bool applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    // apply coil compression coefficient
    bool appy_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, const std::vector< hoNDKLT<T> >& coeff, hoNDArray<T>& dataEigen);

    // ----------------------------------------------------
    // Partial fourier handling for 3DT reconstruction
    // ----------------------------------------------------
    // apply the partial fourier filer along the edges
    bool performPartialFourierFilter(WorkOrderType& workOrder3DT, hoNDArray<T>& kspace);
    // apply the iterative POCS for partial fourier reconstruction
    bool performPartialFourierPOCSRecon(WorkOrderType& workOrder3DT, hoNDArray<T>& kspace);
    // apply the Feng Huang partial fourier reconstruction
    bool performPartialFourierFengHuangRecon(WorkOrderType& workOrder3DT, hoNDArray<T>& kspace);

    //// compute Feng Huang kernel and perform recon
    bool calibFengHuang(WorkOrderType& workOrder3DT, const hoNDArray<T>& src, const hoNDArray<T>& dst, ho6DArray<T>& kernel);
    bool performReconFangHuang(WorkOrderType& workOrder3DT, const hoNDArray<T>& kspaceConj, hoNDArray<T>& kspace, int startRO, int endRO, int startE1, int endE1, int startE2, int endE2, ho6DArray<T>& kernel);

    // estimate job size for 3DT recon
    virtual bool estimateJobSize(gtPlusReconWorkOrder<T>* workOrder3DT, size_t maxNumOfBytesPerJob, size_t overlapBetweenJobs, size_t numOfNodes, size_t& jobSize);

    using BaseClass::partial_fourier_handling_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::verbose_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;

protected:

    // helper memory for computation
    hoNDArray<T> ref_src_;
    hoNDArray<T> ref_dst_;
    hoNDArray<T> data_dst_;
    hoNDArray<T> ref_coil_map_dst_;

    // sampled region along E1/E2
    size_t startE1_;
    size_t endE1_;

    size_t startE2_;
    size_t endE2_;
};

template <typename T> 
bool gtPlusReconWorker3DT<T>::performRecon(WorkOrderType* workOrder3DT)
{
    // the 3DT recon on 5D array [RO E1 E2 CHA N]
    try
    {
        GADGET_CHECK_RETURN_FALSE(workOrder3DT!=NULL);

        if ( !workOrder3DT->workFlow_use_BufferedKernel_ )
        {
            if ( performTiming_ ) { gt_timer1_.start("prepRef"); }
            GADGET_CHECK_RETURN_FALSE(prepRef(workOrder3DT, workOrder3DT->ref_, 
                                            workOrder3DT->ref_recon_, 
                                            workOrder3DT->ref_coil_map_, 
                                            workOrder3DT->start_RO_, workOrder3DT->end_RO_, 
                                            workOrder3DT->start_E1_, workOrder3DT->end_E1_, 
                                            workOrder3DT->start_E2_, workOrder3DT->end_E2_, 
                                            workOrder3DT->data_.get_size(1), workOrder3DT->data_.get_size(2)));
            if ( performTiming_ ) { gt_timer1_.stop(); }

            if ( performTiming_ ) { gt_timer1_.start("coilCompression"); }
            GADGET_CHECK_RETURN_FALSE(coilCompression(workOrder3DT));
            if ( performTiming_ ) { gt_timer1_.stop(); }
        }

        // apply coil compression coefficients
        if ( workOrder3DT->workFlow_use_BufferedKernel_ )
        {
            if ( workOrder3DT->coil_compression_ 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_SPIRIT 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP )
            {
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->data_, debugFolder_+"data_"); }
                this->appy_KLT_coil_compression_coeff_3D(workOrder3DT->data_, *workOrder3DT->coilCompressionCoef_, data_dst_);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(data_dst_, debugFolder_+"data_dst_"); }
            }
            else
            {
                data_dst_ = workOrder3DT->data_;
            }
        }
        else
        {
            if ( workOrder3DT->coil_compression_ 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_SPIRIT 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP 
                && (workOrder3DT->acceFactorE1_>1 || workOrder3DT->acceFactorE2_>1) )
            {
                ref_src_ = workOrder3DT->ref_recon_;

                if ( performTiming_ ) { gt_timer2_.start("apply coil compression ... "); }

                #pragma omp parallel sections default(shared)
                {
                    #pragma omp section
                    {
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_src_, debugFolder_+"ref_src_"); }
                        this->appy_KLT_coil_compression_coeff_3D(ref_src_, *workOrder3DT->coilCompressionCoef_, ref_dst_);
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_dst_, debugFolder_+"ref_dst_"); }
                    }

                    #pragma omp section
                    {
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->data_, debugFolder_+"data_"); }
                        this->appy_KLT_coil_compression_coeff_3D(workOrder3DT->data_, *workOrder3DT->coilCompressionCoef_, data_dst_);
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(data_dst_, debugFolder_+"data_dst_"); }
                    }

                    #pragma omp section
                    {
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->ref_coil_map_, debugFolder_+"ref_coil_map_"); }
                        this->appy_KLT_coil_compression_coeff_3D(workOrder3DT->ref_coil_map_, *workOrder3DT->coilCompressionCoef_, ref_coil_map_dst_);
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_coil_map_dst_, debugFolder_+"ref_coil_map_dst_"); }
                    }
                }

                if ( performTiming_ ) { gt_timer2_.stop(); }

                if ( !workOrder3DT->downstream_coil_compression_ 
                    || workOrder3DT->recon_algorithm_==ISMRMRD_SPIRIT 
                    || workOrder3DT->recon_algorithm_==ISMRMRD_L1SPIRIT_SLEP 
                    || workOrder3DT->recon_algorithm_==ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP )
                {
                    ref_src_ = ref_dst_;
                }
            }
            else
            {
                ref_src_ = workOrder3DT->ref_recon_;
                ref_dst_ = workOrder3DT->ref_recon_;
                data_dst_ = workOrder3DT->data_;
                ref_coil_map_dst_ = workOrder3DT->ref_coil_map_;
            }

            if ( performTiming_ ) { gt_timer1_.start("estimate coil map"); }
            GADGET_CHECK_RETURN_FALSE(this->estimateCoilMap(workOrder3DT, ref_src_, ref_dst_, ref_coil_map_dst_));
            if ( performTiming_ ) { gt_timer1_.stop(); }

            if ( workOrder3DT->acceFactorE1_>1 || workOrder3DT->acceFactorE2_>1 )
            {
                if ( performTiming_ ) { gt_timer1_.start("performCalib"); }
                GADGET_CHECK_RETURN_FALSE(this->performCalib(workOrder3DT, ref_src_, ref_dst_, ref_coil_map_dst_));
                if ( performTiming_ ) { gt_timer1_.stop(); }
            }
        }

        if ( performTiming_ ) { gt_timer1_.start("performUnwrapping"); }
        GADGET_CHECK_RETURN_FALSE(this->performUnwrapping(workOrder3DT, data_dst_));
        if ( performTiming_ ) { gt_timer1_.stop(); }

        if ( performTiming_ ) { gt_timer1_.start("afterUnwrapping"); }
        GADGET_CHECK_RETURN_FALSE(this->afterUnwrapping(workOrder3DT));
        if ( performTiming_ ) { gt_timer1_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performRecon(WorkOrderType* workOrder3DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::
estimateCoilMap(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst)
{
    try
    {
        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t N = workOrder3DT->data_.get_size(4);
        size_t srcCHA = workOrder3DT->data_.get_size(3);

        size_t refRO = ref_dst.get_size(0);
        size_t refE1 = ref_dst.get_size(1);
        size_t refE2 = ref_dst.get_size(2);
        size_t refN = ref_dst.get_size(4);
        size_t dstCHA = ref_coil_map_dst.get_size(3);

        bool same_combinationcoeff_allN = false;
        size_t whichN_combinationcoeff = 0;
        if ( workOrder3DT->CalibMode_ == ISMRMRD_interleaved )
        {
            same_combinationcoeff_allN = true;
            whichN_combinationcoeff = 0;
        }

        if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
        {
            same_combinationcoeff_allN = workOrder3DT->embedded_same_combinationcoeff_allN_;
            whichN_combinationcoeff = workOrder3DT->embedded_whichN_combinationcoeff_;
        }

        if ( workOrder3DT->CalibMode_ == ISMRMRD_separate )
        {
            same_combinationcoeff_allN = workOrder3DT->separate_same_combinationcoeff_allN_;
            whichN_combinationcoeff = workOrder3DT->separate_whichN_combinationcoeff_;
        }

        if ( whichN_combinationcoeff >= refN ) whichN_combinationcoeff=refN-1;

        bool reconKSpace = this->computeKSpace(workOrder3DT);

        // if the coil map has not been preset
        if ( !reconKSpace )
        {
            if ( (workOrder3DT->coilMap_->get_size(0)!=RO) 
                || (workOrder3DT->coilMap_->get_size(1)!=E1)
                || (workOrder3DT->coilMap_->get_size(2)!=E2) )
            {
                if ( same_combinationcoeff_allN )
                {
                    size_t usedN = whichN_combinationcoeff;

                    hoNDArray<T> refCoilMapN(RO, E1, E2, dstCHA, const_cast<T*>(ref_coil_map_dst.begin()+usedN*RO*E1*E2*dstCHA));

                    workOrder3DT->coilMap_->create(RO, E1, E2, dstCHA, refN);

                    hoNDArray<T> coilMapN(RO, E1, E2, dstCHA);

                    hoNDArray<T> buffer3DT(RO, E1, E2, dstCHA);

                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(refCoilMapN, buffer3DT);

                    if ( performTiming_ ) { gt_timer3_.start("coil map estimation ... "); }
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, 
                            coilMapN, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, 
                            workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    memcpy(workOrder3DT->coilMap_->begin()+usedN*RO*E1*E2*dstCHA, coilMapN.begin(), coilMapN.get_number_of_bytes());
                    GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder3DT->coilMap_, usedN));
                }
                else
                {
                    hoNDArray<T> buffer3DT(ref_coil_map_dst.get_dimensions());

                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(ref_coil_map_dst, buffer3DT);

                    if ( performTiming_ ) { gt_timer3_.start("coil map estimation ... "); }
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, 
                            *workOrder3DT->coilMap_, workOrder3DT->coil_map_algorithm_, 
                            workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, 
                            workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                    if ( performTiming_ ) { gt_timer3_.stop(); }
                }
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->coilMap_, debugFolder_+"coilMap_"); }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::estimateCoilMap(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::
performCalib(gtPlusReconWorkOrder3DT<T>* workOrder3DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst)
{
    try
    {
        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t N = workOrder3DT->data_.get_size(4);
        size_t srcCHA = workOrder3DT->data_.get_size(3);

        size_t refRO = ref_dst.get_size(0);
        size_t refE1 = ref_dst.get_size(1);
        size_t refE2 = ref_dst.get_size(2);
        size_t refN = ref_dst.get_size(4);
        size_t dstCHA = ref_coil_map_dst.get_size(3);

        bool same_combinationcoeff_allN = false;
        size_t whichN_combinationcoeff = 0;
        if ( workOrder3DT->CalibMode_ == ISMRMRD_interleaved )
        {
            same_combinationcoeff_allN = true;
            whichN_combinationcoeff = 0;
        }

        if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
        {
            same_combinationcoeff_allN = workOrder3DT->embedded_same_combinationcoeff_allN_;
            whichN_combinationcoeff = workOrder3DT->embedded_whichN_combinationcoeff_;
        }

        if ( workOrder3DT->CalibMode_ == ISMRMRD_separate )
        {
            same_combinationcoeff_allN = workOrder3DT->separate_same_combinationcoeff_allN_;
            whichN_combinationcoeff = workOrder3DT->separate_whichN_combinationcoeff_;
        }

        if ( whichN_combinationcoeff >= refN ) whichN_combinationcoeff=refN-1;

        bool reconKSpace = this->computeKSpace(workOrder3DT);

        // calibration
        if ((workOrder3DT->kernel_->get_size(3) != srcCHA) || (workOrder3DT->kernel_->get_size(4) != dstCHA))
        {
           GADGET_CHECK_RETURN_FALSE(this->performCalibPrep(ref_src, ref_dst, workOrder3DT));

            // perform calibration
            if ( same_combinationcoeff_allN )
            {
                size_t usedN = whichN_combinationcoeff;

                this->performCalibImpl(ref_src, ref_dst, workOrder3DT, usedN);

                GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder3DT->kernel_, usedN));

                if ( !reconKSpace )
                {
                    if ( workOrder3DT->unmixingCoeffIm_ && (workOrder3DT->unmixingCoeffIm_->get_number_of_elements()>0) )
                    {
                        GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder3DT->unmixingCoeffIm_, usedN));
                        GADGET_CHECK_RETURN_FALSE(repmatLastDimension(workOrder3DT->gfactor_, usedN));

                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->unmixingCoeffIm_, debugFolder_+"unmixingCoeffIm_"); }
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->gfactor_, debugFolder_+"gfactor_"); }
                    }
                }
            }
            else
            {
                int usedN;
                #ifdef USE_OMP
                    omp_set_nested(1);
                #endif // USE_OMP

                #pragma omp parallel for default(none) private(usedN) shared(N, ref_src, ref_dst, workOrder3DT, reconKSpace)
                for ( usedN=0; usedN<(int)N; usedN++ )
                {
                    this->performCalibImpl(ref_src, ref_dst, workOrder3DT, usedN);
                }

                #ifdef USE_OMP
                    omp_set_nested(0);
                #endif // USE_OMP
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performCalib(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::
performCalibPrep(const hoNDArray<T>& , const hoNDArray<T>& , WorkOrderType* /*workOrder3DT*/)
{
    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::
performCalibImpl(const hoNDArray<T>& , const hoNDArray<T>& , WorkOrderType* , size_t )
{
    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::performUnwrapping(WorkOrderType* , const hoNDArray<T>& )
{
    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::performRefFilter(WorkOrderType* workOrder3DT, 
                                            const hoNDArray<T>& ref, 
                                            hoNDArray<T>& refFiltered, 
                                            int startRO, int endRO, 
                                            int startE1, int endE1, 
                                            int startE2, int endE2)
{
    try
    {
        refFiltered = ref;

        size_t RO = ref.get_size(0);
        size_t E1 = ref.get_size(1);
        size_t E2 = ref.get_size(2);

        if ( workOrder3DT->filterROE1E2_ref_.get_size(0)==RO 
            && workOrder3DT->filterROE1E2_ref_.get_size(1)==E1 
            && workOrder3DT->filterROE1E2_ref_.get_size(2)==E2 )
        {
            Gadgetron::apply_kspace_filter_ROE1E2(ref, workOrder3DT->filterROE1E2_ref_, refFiltered);
        }
        else if ( (workOrder3DT->filterRO_ref_.get_number_of_elements()==RO) 
            && (workOrder3DT->filterE1_ref_.get_number_of_elements()==E1) 
            && (workOrder3DT->filterE2_ref_.get_number_of_elements()==E2) )
        {
            Gadgetron::apply_kspace_filter_ROE1E2(ref, workOrder3DT->filterRO_ref_, workOrder3DT->filterE1_ref_, workOrder3DT->filterE2_ref_, refFiltered);
        }
        else
        {
            if ( (workOrder3DT->filterRO_ref_.get_number_of_elements()==RO) 
                && (workOrder3DT->filterE1_ref_.get_number_of_elements()!=E1) 
                && (workOrder3DT->filterE2_ref_.get_number_of_elements()!=E2) )
            {
                Gadgetron::apply_kspace_filter_RO(ref, workOrder3DT->filterRO_ref_, refFiltered);
            }

            if ( (workOrder3DT->filterRO_ref_.get_number_of_elements()!=RO) 
                && (workOrder3DT->filterE1_ref_.get_number_of_elements()==E1) 
                && (workOrder3DT->filterE2_ref_.get_number_of_elements()!=E2) )
            {
                Gadgetron::apply_kspace_filter_E1(ref, workOrder3DT->filterE1_ref_, refFiltered);
            }

            if ( (workOrder3DT->filterRO_ref_.get_number_of_elements()!=RO) 
                && (workOrder3DT->filterE1_ref_.get_number_of_elements()!=E1) 
                && (workOrder3DT->filterE2_ref_.get_number_of_elements()==E2) )
            {
                Gadgetron::apply_kspace_filter_E2(ref, workOrder3DT->filterE2_ref_, refFiltered);
            }

            if ( (workOrder3DT->filterRO_ref_.get_number_of_elements()==RO) 
                && (workOrder3DT->filterE1_ref_.get_number_of_elements()==E1) 
                && (workOrder3DT->filterE2_ref_.get_number_of_elements()!=E2) )
            {
                Gadgetron::apply_kspace_filter_ROE1(ref, workOrder3DT->filterRO_ref_, workOrder3DT->filterE1_ref_, refFiltered);
            }

            if ( (workOrder3DT->filterRO_ref_.get_number_of_elements()==RO) 
                && (workOrder3DT->filterE1_ref_.get_number_of_elements()!=E1) 
                && (workOrder3DT->filterE2_ref_.get_number_of_elements()==E2) )
            {
                Gadgetron::apply_kspace_filter_ROE2(ref, workOrder3DT->filterRO_ref_, workOrder3DT->filterE2_ref_, refFiltered);
            }

            if ( (workOrder3DT->filterRO_ref_.get_number_of_elements()!=RO) 
                && (workOrder3DT->filterE1_ref_.get_number_of_elements()==E1) 
                && (workOrder3DT->filterE2_ref_.get_number_of_elements()==E2) )
            {
                Gadgetron::apply_kspace_filter_E1E2(ref, workOrder3DT->filterE1_ref_, workOrder3DT->filterE2_ref_, refFiltered);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performRefFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::prepRefByAveragingCrossN(WorkOrderType* workOrder3DT, const hoNDArray<T>& ref, bool averageAllRef, int numOfModes, hoNDArray<T>& refRecon)
{
    try
    {
        size_t RO = ref.get_size(0);
        size_t E1 = ref.get_size(1);
        size_t E2 = ref.get_size(2);
        size_t CHA = ref.get_size(3);
        size_t N = ref.get_size(4);

        if ( !averageAllRef && ( (numOfModes<1) || (numOfModes>N-1) ) )
        {
            refRecon = ref;
        }
        else if ( averageAllRef && ( (numOfModes<1) || (numOfModes>N-1) ) )
        {
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace5D(ref, refRecon));
        }
        else if ( averageAllRef && (numOfModes>=1) && (numOfModes<=N-1) )
        {
            hoNDArray<T> refKLF(RO, E1, E2, CHA, N);
            Gadgetron::clear(refKLF);

            hoMatrix<T> A(RO*E1*E2*CHA, N, const_cast<T*>(ref.begin()));
            hoMatrix<T> A_KLF(RO*E1*E2*CHA, N, refKLF.begin());
            hoNDKLT<T> klt;
            klt.prepare(A, 1, (size_t)0);
            klt.KL_filter(A, A_KLF, 1, numOfModes);
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refKLF, debugFolder_+"refKLF"); }

            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace5D(refKLF, refRecon));
        }
        else if ( !averageAllRef && (numOfModes>=1) && (numOfModes<=N-1) )
        {
            refRecon.create(RO, E1, E2, CHA, N);
            Gadgetron::clear(refRecon);

            hoMatrix<T> A(RO*E1*E2*CHA, N, const_cast<T*>(ref.begin()));
            hoMatrix<T> A_KLF(RO*E1*E2*CHA, N, refRecon.begin());
            hoNDKLT<T> klt;
            klt.prepare(A, 1, (size_t)0);
            klt.KL_filter(A, A_KLF, 1, numOfModes);
        }
        else
        {
            refRecon = ref;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::prepRefByAveragingCrossN(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::prepRef(WorkOrderType* workOrder3DT, const hoNDArray<T>& ref, 
                                hoNDArray<T>& refRecon, hoNDArray<T>& refCoilMap, 
                                int startRO, int endRO, 
                                int startE1, int endE1, 
                                int startE2, int endE2, 
                                size_t dataE1, size_t dataE2)
{
    try
    {
        size_t dataRO = workOrder3DT->data_.get_size(0);
        size_t dataN = workOrder3DT->data_.get_size(4);

        size_t RO = ref.get_size(0);
        size_t E1 = ref.get_size(1);
        size_t E2 = ref.get_size(2);
        size_t srcCHA = ref.get_size(3);
        size_t N = ref.get_size(4);

        if ( workOrder3DT->acceFactorE1_ == 1 && workOrder3DT->acceFactorE2_ == 1 )
        {
            if ( workOrder3DT->no_acceleration_averageall_ref_ )
            {
                GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder3DT, ref, workOrder3DT->no_acceleration_averageall_ref_, 0, refRecon));
            }

            GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder3DT, refRecon, refCoilMap, startRO, endRO, startE1, endE1, startE2, endE2));
        }
        else if ( workOrder3DT->CalibMode_ == ISMRMRD_interleaved )
        {
            GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder3DT, ref, true, 0, refRecon));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_+"refRecon_interleaved"); }

            GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder3DT, refRecon, refCoilMap, startRO, endRO, startE1, endE1, startE2, endE2));

            if ( (startRO>=0 && endRO>0 && endRO>startRO) || (startE1>=0 && endE1>0 && endE1>startE1) || (startE2>=0 && endE2>0 && endE2>startE2) )
            {
                std::vector<size_t> crop_offset(5), crop_size(5);

                crop_offset[0] = 0;
                crop_offset[1] = 0;
                crop_offset[2] = 0;
                crop_offset[3] = 0;
                crop_offset[4] = 0;

                crop_size[0] = RO;
                crop_size[1] = E1;
                crop_size[2] = refRecon.get_size(2);
                crop_size[3] = refRecon.get_size(3);
                crop_size[4] = refRecon.get_size(4);

                if (startRO>=0 && endRO>0 && endRO>startRO)
                {
                    crop_offset[0] = startRO;
                    crop_size[0] = endRO-startRO+1;
                }

                if (startE1>=0 && endE1>0 && endE1>startE1)
                {
                    crop_offset[1] = startE1;
                    crop_size[1] = endE1-startE1+1;
                }

                if (startE2>=0 && endE2>0 && endE2>startE2)
                {
                    crop_offset[2] = startE2;
                    crop_size[2] = endE2-startE2+1;
                }

                hoNDArray<T> croppedRef;
                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                refRecon = croppedRef;
            }
        }
        else if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded 
                || workOrder3DT->CalibMode_ == ISMRMRD_separate 
                || workOrder3DT->CalibMode_ == ISMRMRD_external )
        {
            if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
            {
                if (performTiming_) { gt_timer2_.start("detectSampledRegionE1E2 ... "); }
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.detectSampledRegionE1E2(ref, startE1_, endE1_, startE2_, endE2_));
                if (performTiming_) { gt_timer2_.stop(); }
            }

            if ( workOrder3DT->CalibMode_ == ISMRMRD_separate )
            {
                GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder3DT, ref, workOrder3DT->separate_averageall_ref_, 0, refRecon));

                if (performTiming_) { gt_timer2_.start("detectSampledRegionE1E2 ... "); }
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.detectSampledRegionE1E2(refRecon, startE1_, endE1_, startE2_, endE2_));
                if (performTiming_) { gt_timer2_.stop(); }

                if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_ + "refRecon_beforeCrop"); }
            }

            std::vector<size_t> crop_offset(5);
            crop_offset[0] = 0;
            crop_offset[1] = startE1_;
            crop_offset[2] = startE2_;
            crop_offset[3] = 0;
            crop_offset[4] = 0;

            std::vector<size_t> crop_size(5);

            crop_size[1] = endE1_ - startE1_ + 1;
            crop_size[2] = endE2_ - startE2_ + 1;
            crop_size[3] = srcCHA;

            if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
            {
                crop_size[0] = ref.get_size(0);
                crop_size[4] = ref.get_size(4);

                if ( performTiming_ ) { gt_timer2_.start("crop sampled region ... "); }
                hoNDArray<T> croppedRef;
                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(ref, croppedRef, crop_offset, crop_size));
                if ( performTiming_ ) { gt_timer2_.stop(); }
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(croppedRef, debugFolder_+"refRecon_afterCrop"); }

                if ( workOrder3DT->recon_algorithm_ == ISMRMRD_SPIRIT 
                    || workOrder3DT->recon_algorithm_ == ISMRMRD_L1SPIRIT_SLEP 
                    || workOrder3DT->recon_algorithm_ == ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP )
                {
                    // copy the ref into the data
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.copyAlongROE1E2(ref, workOrder3DT->data_, 0, refRecon.get_size(0)-1, startE1_, endE1_, startE2_, endE2_));
                }

                GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder3DT, croppedRef, workOrder3DT->embedded_averageall_ref_, 0, refRecon));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_+"refRecon_afterCrop_prepCrossN"); }

                crop_size[4] = refRecon.get_size(4);

                if ( performTiming_ ) { gt_timer2_.start("set up ref for coil map ... "); }
                refCoilMap.create(RO, E1, E2, srcCHA, refRecon.get_size(4));
                GADGET_CHECK_RETURN_FALSE(setSubArrayUpTo11DArray(refRecon, refCoilMap, crop_offset, crop_size));
                if ( performTiming_ ) { gt_timer2_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"refCoilMap"); }

                if ( performTiming_ ) { gt_timer2_.start("perform ref coil map filter ... "); }
                GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder3DT, refCoilMap, refCoilMap, startRO, endRO, startE1, endE1, startE2, endE2));
                if ( performTiming_ ) { gt_timer2_.stop(); }
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"refCoilMap_filtered"); }

                if ( refRecon.get_size(0) == RO )
                {
                    if ( startRO>=0 && endRO>0 && endRO>startRO && startRO<RO && endRO<RO )
                    {
                        crop_offset[0] = startRO;
                        crop_size[0] = endRO-startRO+1;

                        crop_offset[1] = 0;
                        crop_size[1] = refRecon.get_size(1);

                        crop_offset[2] = 0;
                        crop_size[2] = refRecon.get_size(2);
                    }
                }

                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                refRecon = croppedRef;
            }
            else
            {
                crop_size[0] = refRecon.get_size(0);
                crop_size[4] = refRecon.get_size(4);

                hoNDArray<T> croppedRef;
                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(croppedRef, debugFolder_+"croppedRef"); }

                GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder3DT, croppedRef, refCoilMap, startRO, endRO, startE1, endE1, startE2, endE2));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"croppedRef_filtered"); }

                refRecon = croppedRef;

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(dataRO, dataE1, dataE2, &refCoilMap, &croppedRef));
                refCoilMap = croppedRef;
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"refCoilMap"); }

                if ( refRecon.get_size(0) == RO )
                {
                    if ( startRO>=0 && endRO>0 && endRO>startRO && startRO<RO && endRO<RO )
                    {
                        crop_offset[0] = startRO;
                        crop_size[0] = endRO-startRO+1;

                        crop_offset[1] = 0;
                        crop_size[1] = refRecon.get_size(1);

                        crop_offset[2] = 0;
                        crop_size[2] = refRecon.get_size(2);

                        GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                        refRecon = croppedRef;
                    }
                }
            }

            // if the ref N is smaller than the data N, e.g. in some cases with the separate mode
            // make sure every data N have its ref data
            if ( N < dataN )
            {
                hoNDArray<T> refReconDataN(refRecon.get_size(0), refRecon.get_size(1), refRecon.get_size(2), refRecon.get_size(3), dataN);
                hoNDArray<T> refCoilMapDataN(refCoilMap.get_size(0), refCoilMap.get_size(1), refCoilMap.get_size(2), refCoilMap.get_size(3), dataN);

                memcpy(refReconDataN.begin(), refRecon.begin(), refRecon.get_number_of_bytes());
                memcpy(refCoilMapDataN.begin(), refCoilMap.begin(), refCoilMap.get_number_of_bytes());

                size_t refReconN4D = refRecon.get_size(0)*refRecon.get_size(1)*refRecon.get_size(2)*refRecon.get_size(3);
                size_t refCoilMapN4D = refCoilMap.get_size(0)*refCoilMap.get_size(1)*refCoilMap.get_size(2)*refCoilMap.get_size(3);

                size_t n;
                for ( n=N; n<dataN; n++ )
                {
                    memcpy(refReconDataN.begin()+n*refReconN4D, refRecon.begin()+(N-1)*refReconN4D, sizeof(T)*refReconN4D);
                    memcpy(refCoilMapDataN.begin()+n*refCoilMapN4D, refCoilMap.begin()+(N-1)*refCoilMapN4D, sizeof(T)*refCoilMapN4D);
                }

                refRecon = refReconDataN;
                refCoilMap = refCoilMapDataN;
            }
        }
        else
        {
            GERROR_STREAM("CalibMode is not supported in gtPlusReconWorker3DT<T>::prepRef(...) : " << workOrder3DT->CalibMode_);
            return false;
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_+"refRecon"); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"refCoilMap"); }

        // if the upstream coil compression is needed
        if ( workOrder3DT->upstream_coil_compression_ )
        {
            if ( !debugFolder_.empty() ) { GDEBUG_STREAM("Upstream coil compression ... "); }

            if ( performTiming_ ) { gt_timer2_.start("average along N ... "); }
            hoNDArray<T> aveAll;
            if (refRecon.get_size(4) > 1)
            {
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace5D(refRecon, aveAll));
                aveAll.squeeze();
            }
            else
            {
                aveAll.create(refRecon.get_size(0), refRecon.get_size(1), refRecon.get_size(2), refRecon.get_size(3), refRecon.begin());
            }
            if ( performTiming_ ) { gt_timer2_.stop(); }

            std::vector<hoNDKLT<T> > upstreamKLTfRef(workOrder3DT->ref_.get_size(4)), upstreamKLTfRefRecon(refRecon.get_size(4)), upstreamKLTfData(workOrder3DT->data_.get_size(4));

            if ( performTiming_ ) { gt_timer2_.start("compute coil compression coefficients ... "); }
            hoNDArray<T> eigenValues;
            hoNDKLT<T> klt;
            if ( workOrder3DT->upstream_coil_compression_num_modesKept_ > 0 )
            {
                klt.prepare(aveAll, 3, (size_t)workOrder3DT->upstream_coil_compression_num_modesKept_);
                klt.eigen_value(eigenValues);
            }
            else
            {
                klt.prepare(aveAll, 3, (value_type)workOrder3DT->upstream_coil_compression_thres_);
                klt.eigen_value(eigenValues);
            }
            if ( performTiming_ ) { gt_timer2_.stop(); }

            {
                eigenValues.print(std::cout);
                for (size_t i = 0; i<eigenValues.get_size(0); i++)
                {
                    GDEBUG_STREAM(i << " = " << eigenValues(i));
                }
            }
            GDEBUG_STREAM("Upstream coil compression, number of channel kept is " << klt.output_length());

            size_t n;
            for ( n=0; n<upstreamKLTfRef.size(); n++ )
            {
                upstreamKLTfRef[n] = klt;
            }

            for ( n=0; n<upstreamKLTfRefRecon.size(); n++ )
            {
                upstreamKLTfRefRecon[n] = klt;
            }

            for ( n=0; n<upstreamKLTfData.size(); n++ )
            {
                upstreamKLTfData[n] = klt;
            }

            if (klt.output_length()<srcCHA)
            {
                if ( performTiming_ ) { gt_timer2_.start("apply upstream coil compression ... "); }
                {
                    {
                        if ( performTiming_ ) { gt_timer3_.start("appy_KLT_coil_compression_coeff_3D ... "); }
                        this->appy_KLT_coil_compression_coeff_3D(workOrder3DT->data_, upstreamKLTfData, data_dst_);
                        if ( performTiming_ ) { gt_timer3_.stop(); }

                        if ( performTiming_ ) { gt_timer3_.start("copy data ... "); }
                        workOrder3DT->data_ = data_dst_;
                        if ( performTiming_ ) { gt_timer3_.stop(); }

                        data_dst_.clear();
                    }

                    {
                        this->appy_KLT_coil_compression_coeff_3D(workOrder3DT->ref_, upstreamKLTfRef, ref_dst_);
                        workOrder3DT->ref_ = ref_dst_;

                        ref_dst_.clear();
                    }

                    {
                        hoNDArray<T> refRecon_upstream;
                        this->appy_KLT_coil_compression_coeff_3D(refRecon, upstreamKLTfRefRecon, refRecon_upstream);
                        refRecon = refRecon_upstream;
                        refRecon_upstream.clear();
                    }

                    {
                        hoNDArray<T> refCoilMap_upstream;
                        this->appy_KLT_coil_compression_coeff_3D(refCoilMap, upstreamKLTfRefRecon, refCoilMap_upstream);
                        refCoilMap = refCoilMap_upstream;
                        refCoilMap_upstream.clear();
                    }
                }

                if ( performTiming_ ) { gt_timer2_.stop(); }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::prepRef(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::coilCompression(WorkOrderType* workOrder3DT)
{
    // the 3DT recon on 5D array [RO E1 E2 CHA N]
    try
    {
        size_t RO = workOrder3DT->ref_recon_.get_size(0);
        size_t E1 = workOrder3DT->ref_recon_.get_size(1);
        size_t E2 = workOrder3DT->ref_recon_.get_size(2);
        size_t srcCHA = workOrder3DT->ref_recon_.get_size(3);
        size_t N = workOrder3DT->ref_recon_.get_size(4);

        size_t dataN = workOrder3DT->data_.get_size(4);

        size_t n;

        if (workOrder3DT->CalibMode_ == ISMRMRD_noacceleration) return true;

        // compute coil compression coeff
        if ( workOrder3DT->coil_compression_ 
            && workOrder3DT->recon_algorithm_!=ISMRMRD_SPIRIT 
            && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP 
            && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP )
        {
            // check whether coil compression coeff has been preset
            if ( workOrder3DT->coilCompressionCoef_->size()!=dataN )
            {
                if ( workOrder3DT->same_coil_compression_coeff_allN_ )
                {
                    hoNDArray<T> aveAll;
                    if (workOrder3DT->ref_recon_.get_size(4) > 1)
                    {
                        GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace5D(workOrder3DT->ref_recon_, aveAll));
                        aveAll.squeeze();
                    }
                    else
                    {
                        aveAll.create(RO, E1, E2, srcCHA, workOrder3DT->ref_recon_.begin());
                    }
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(aveAll, debugFolder_+"aveAll"); }

                    Gadgetron::hoNDKLT<T> klt;
                    hoNDArray<T> eigenValues;
                    if ( workOrder3DT->coil_compression_num_modesKept_ > 0 )
                    {
                        klt.prepare(aveAll, 3, (size_t)workOrder3DT->upstream_coil_compression_num_modesKept_);
                        klt.eigen_value(eigenValues);
                    }
                    else
                    {
                        klt.prepare(aveAll, 3, (value_type)workOrder3DT->coil_compression_thres_);
                        klt.eigen_value(eigenValues);
                    }

                    workOrder3DT->coilCompressionCoef_->resize(dataN);

                    for ( n=0; n<dataN; n++ )
                    {
                        (*workOrder3DT->coilCompressionCoef_)[n] = klt;
                    }

                    if (!debugFolder_.empty())
                    {
                        eigenValues.print(std::cout);
                        for (size_t i = 0; i<eigenValues.get_size(0); i++)
                        {
                            GDEBUG_STREAM(i << " = " << eigenValues(i));
                        }
                    }
                    GDEBUG_STREAM("Coil compression, number of channel kept is " << klt.output_length());
                }
                else
                {
                    std::vector<size_t> allNDim(4);
                    allNDim[0] = RO;
                    allNDim[1] = E1;
                    allNDim[2] = E2;
                    allNDim[3] = srcCHA;

                    size_t num_modesKept = srcCHA;

                    workOrder3DT->coilCompressionCoef_->resize(N);

                    for ( n=0; n<N; n++ )
                    {
                        hoNDArray<T> dataCurrN(&allNDim, workOrder3DT->ref_recon_.begin()+n*RO*E1*E2*srcCHA, false);

                        hoNDKLT<T> klt;
                        hoNDArray<T> eigenValues;

                        if ( n == 0 )
                        {
                            if ( workOrder3DT->coil_compression_num_modesKept_ > 0 )
                            {
                                klt.prepare(dataCurrN, 3, (size_t)workOrder3DT->upstream_coil_compression_num_modesKept_);
                                klt.eigen_value(eigenValues);
                            }
                            else
                            {
                                klt.prepare(dataCurrN, 3, (value_type)workOrder3DT->coil_compression_thres_);
                                klt.eigen_value(eigenValues);
                            }

                            num_modesKept = klt.output_length();
                            (*workOrder3DT->coilCompressionCoef_)[0] = klt;
                        }
                        else
                        {
                            klt.prepare(dataCurrN, 3, (size_t)num_modesKept);
                            klt.eigen_value(eigenValues);

                            (*workOrder3DT->coilCompressionCoef_)[n] = klt;
                        }

                        if (!debugFolder_.empty())
                        {
                            eigenValues.print(std::cout);
                            for (size_t i = 0; i<eigenValues.get_size(0); i++)
                            {
                                GDEBUG_STREAM(i << " = " << eigenValues(i));
                            }
                        }
                        GDEBUG_STREAM("Coil compression, number of channel kept is " << klt.output_length());
                    }
                }
            }

            if ( N < dataN )
            {
                std::vector<hoNDKLT<T> > coilCompressionCoef(dataN);
                for ( n=0; n<N; n++ )
                {
                    coilCompressionCoef[n] = (*workOrder3DT->coilCompressionCoef_)[n];
                }

                for ( n=N; n<dataN; n++ )
                {
                    coilCompressionCoef[n] = (*workOrder3DT->coilCompressionCoef_)[N-1];
                }

                *(workOrder3DT->coilCompressionCoef_) = coilCompressionCoef;
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::coilCompression(WorkOrderType* workOrder3DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t E2 = kerIm.get_size(2);
        size_t srcCHA = kerIm.get_size(3);
        size_t dstCHA = kerIm.get_size(4);

        GADGET_CHECK_RETURN_FALSE(coilMap.get_size(0)==RO);
        GADGET_CHECK_RETURN_FALSE(coilMap.get_size(1)==E1);
        GADGET_CHECK_RETURN_FALSE(coilMap.get_size(2)==E2);
        GADGET_CHECK_RETURN_FALSE(coilMap.get_size(3)==dstCHA);

        unmixCoeff.create(RO, E1, E2, srcCHA);
        Gadgetron::clear(&unmixCoeff);

        int src;

        T* pKerIm = const_cast<T*>(kerIm.begin());
        T* pCoilMap = const_cast<T*>(coilMap.begin());
        T* pCoeff = unmixCoeff.begin();

        std::vector<size_t> dim(3);
        dim[0] = RO;
        dim[1] = E1;
        dim[2] = E2;

        #pragma omp parallel default(none) private(src) shared(RO, E1, E2, srcCHA, dstCHA, pKerIm, pCoilMap, pCoeff, dim)
        {
            hoNDArray<T> coeff2D, coeffTmp(&dim);
            hoNDArray<T> coilMap2D;
            hoNDArray<T> kerIm2D;

            #pragma omp for
            for ( src=0; src<(int)srcCHA; src++ )
            {
                coeff2D.create(&dim, pCoeff+src*RO*E1*E2);

                for ( size_t dst=0; dst<dstCHA; dst++ )
                {
                    kerIm2D.create(&dim, pKerIm+src*RO*E1*E2+dst*RO*E1*E2*srcCHA);
                    coilMap2D.create(&dim, pCoilMap+dst*RO*E1*E2);
                    Gadgetron::multiplyConj(kerIm2D, coilMap2D, coeffTmp);
                    Gadgetron::add(coeff2D, coeffTmp, coeff2D);
                }
            }
        }

        hoNDArray<T> conjUnmixCoeff(unmixCoeff);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiplyConj(unmixCoeff, conjUnmixCoeff, conjUnmixCoeff));

        gFactor.create(RO, E1, E2);
        Gadgetron::clear(&gFactor);

        hoNDArray<T> gFactorBuf(RO, E1, E2, 1, gFactor.begin());
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(conjUnmixCoeff, gFactorBuf, 3));
        Gadgetron::sqrt(gFactor, gFactor);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::applyImageDomainKernel(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t E2 = kerIm.get_size(2);
        size_t srcCHA = kerIm.get_size(3);
        size_t dstCHA = kerIm.get_size(4);

        GADGET_CHECK_RETURN_FALSE(kspace.get_size(0)==RO);
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(1)==E1);
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(2)==E2);
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(3)==srcCHA);

        // buffer3DT_unwrapping_ = kspace;

        hoNDArray<T> buffer3DT(kspace.get_dimensions());

        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspace, buffer3DT));

        GADGET_CHECK_RETURN_FALSE(applyImageDomainKernelImage(buffer3DT, kerIm, complexIm));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::applyImageDomainKernel(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm)
{
    return applyImageDomainKernelImage(aliasedIm, kerIm, this->buf4D, complexIm);
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& kerImBuffer, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t E2 = kerIm.get_size(2);
        size_t srcCHA = kerIm.get_size(3);
        size_t dstCHA = kerIm.get_size(4);

        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(0)==RO);
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(1)==E1);
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(2)==E2);
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(3)==srcCHA);

        boost::shared_ptr< std::vector<size_t> > dim = aliasedIm.get_dimensions();
        std::vector<size_t> dimIm(*dim);
        dimIm[3] = dstCHA;

        if ( !complexIm.dimensions_equal(&dimIm) )
        {
            complexIm.create(&dimIm);
        }
        Gadgetron::clear(&complexIm);

        std::vector<size_t> dim4D(4);
        dim4D[0] = RO;
        dim4D[1] = E1;
        dim4D[2] = E2;
        dim4D[3] = srcCHA;

        std::vector<size_t> dimIm4D(4);
        dimIm4D[0] = RO;
        dimIm4D[1] = E1;
        dimIm4D[2] = E2;
        dimIm4D[3] = dstCHA;

        size_t num = aliasedIm.get_number_of_elements()/ (RO*E1*E2*srcCHA);

        int n;

        #pragma omp parallel default(none) private(n) shared(num, dim4D, aliasedIm, RO, E1, E2, srcCHA, dstCHA, kerIm, complexIm) num_threads( (int)((num<16) ? num : 16) )
        {
            hoNDArray<T> unwrapped4D(RO, E1, E2, srcCHA);

            #pragma omp for
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> buf4D(&dim4D, const_cast<T*>(aliasedIm.begin()+n*RO*E1*E2*srcCHA));

                int dCha;

                for ( dCha=0; dCha<(int)dstCHA; dCha++ )
                {
                    hoNDArray<T> kerIm4D(RO, E1, E2, srcCHA, const_cast<T*>(kerIm.begin()+dCha*RO*E1*E2*srcCHA));
                    hoNDArray<T> complexIm3D(RO, E1, E2, 1, complexIm.begin()+n*RO*E1*E2*dstCHA+dCha*RO*E1*E2);
                    Gadgetron::multiply(kerIm4D, buf4D, unwrapped4D);
                    Gadgetron::sum_over_dimension(unwrapped4D, complexIm3D, 3);
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& kerImBuffer, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(2)==unmixCoeff.get_size(2));
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(3)==unmixCoeff.get_size(3));

        // buffer3DT_unwrapping_ = kspace;
        hoNDArray<T> buffer3DT(kspace.get_dimensions());

        GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspace, buffer3DT));
        GADGET_CHECK_RETURN_FALSE(applyUnmixCoeffImage(buffer3DT, unmixCoeff, complexIm));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(2)==unmixCoeff.get_size(2));
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(3)==unmixCoeff.get_size(3));

        boost::shared_ptr< std::vector<size_t> > dim = aliasedIm.get_dimensions();
        std::vector<size_t> dimIm(*dim);
        dimIm[3] = 1;

        if ( !complexIm.dimensions_equal(&dimIm) )
        {
            complexIm.create(&dimIm);
        }
        Gadgetron::clear(&complexIm);

        hoNDArray<T> buffer3DT(aliasedIm.get_dimensions());

        Gadgetron::multiply(aliasedIm, unmixCoeff, buffer3DT);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(buffer3DT, complexIm, 3));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T>
bool gtPlusReconWorker3DT<T>::appy_KLT_coil_compression_coeff_3D(const hoNDArray<T>& data, const std::vector< hoNDKLT<T> >& coeff, hoNDArray<T>& dataEigen)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim >= 4);

        GADGET_CHECK_RETURN_FALSE(coeff.size() >= data.get_size(NDim - 1));

        size_t LastDim = coeff.size();
        size_t dstCHA = coeff[0].output_length();

        size_t n;
        for (n = 1; n<LastDim; n++)
        {
            GADGET_CHECK_RETURN_FALSE(coeff[n].output_length() == dstCHA);
        }

        size_t LastDimData = data.get_size(NDim - 1);
        std::vector<size_t> dim;
        data.get_dimensions(dim);
        long long N = data.get_number_of_elements() / LastDimData;

        std::vector<size_t> dimEigen(dim);
        dimEigen[3] = dstCHA;

        dataEigen.create(&dimEigen);
        long long eigenN = dataEigen.get_number_of_elements() / LastDimData;

        std::vector<size_t> dimLastDim(NDim - 1);
        for (n = 0; n<NDim - 1; n++)
        {
            dimLastDim[n] = dim[n];
        }

        std::vector<size_t> dimEigenLastDim(dimLastDim);
        dimEigenLastDim[3] = dstCHA;

        if (LastDimData>1)
        {
            hoNDArray<T> dataEigenLastDim;
            for (n = 0; n<LastDimData; n++)
            {
                hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin() + n*N));
                coeff[n].transform(dataLastDim, dataEigenLastDim, 3);
                memcpy(dataEigen.begin() + n*eigenN, dataEigenLastDim.begin(), dataEigenLastDim.get_number_of_bytes());
            }
        }
        else
        {
            hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin()));
            hoNDArray<T> dataEigenLastDim(&dimEigenLastDim, dataEigen.begin());
            coeff[0].transform(dataLastDim, dataEigenLastDim, 3);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::appy_KLT_coil_compression_coeff_3D(std::vector<hoNDKLT<T> >& coeff) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::afterUnwrapping(WorkOrderType* workOrder3DT)
{
    try
    {
        bool fullres_coilmap = false;
        bool ref_fillback = false;
        bool averageallN_coilmap = false;
        bool same_coilmap_allN = false;
        size_t whichN_coilmap = 0;

        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t srcCHA = workOrder3DT->kernel_->get_size(3);
        size_t dstCHA = workOrder3DT->kernel_->get_size(4);
        size_t N = workOrder3DT->data_.get_size(4);

        if ( workOrder3DT->CalibMode_ == ISMRMRD_embedded )
        {
            if ( workOrder3DT->embedded_fullres_coilmap_ )
            {
                fullres_coilmap = true;
            }

            if ( workOrder3DT->embedded_ref_fillback_ )
            {
                ref_fillback = true;
            }

            if ( workOrder3DT->embedded_averageall_ref_ )
            {
                averageallN_coilmap = true;
            }

            if ( workOrder3DT->embedded_same_combinationcoeff_allN_ )
            {
                same_coilmap_allN = true;
                whichN_coilmap = workOrder3DT->embedded_whichN_combinationcoeff_;
            }
        }

        if ( workOrder3DT->CalibMode_ == ISMRMRD_separate )
        {
            if ( workOrder3DT->separate_fullres_coilmap_ )
            {
                fullres_coilmap = true;
            }

            if ( workOrder3DT->separate_averageall_ref_ )
            {
                averageallN_coilmap = true;
            }

            if ( workOrder3DT->separate_same_combinationcoeff_allN_ )
            {
                same_coilmap_allN = true;
                whichN_coilmap = workOrder3DT->separate_whichN_combinationcoeff_;
            }
        }

        if ( whichN_coilmap >= N ) whichN_coilmap = N-1;

        if ( ref_fillback )
        {
            if ( performTiming_ ) { gt_timer2_.start("ref fill back ... "); }

            hoNDArray<T> ref_dst;
            if ( workOrder3DT->coil_compression_ 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_SPIRIT 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP 
                && workOrder3DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP )
            {
                this->appy_KLT_coil_compression_coeff_3D(workOrder3DT->ref_, *workOrder3DT->coilCompressionCoef_, ref_dst);
            }
            else
            {
                ref_dst = workOrder3DT->ref_;
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_dst, debugFolder_+"ref_dst"); }

            if ( (ref_dst.get_size(3)==dstCHA) && (ref_dst.get_size(4)==N) )
            {
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->fullkspace_, debugFolder_+"fullkspace_"); }

                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.copyAlongROE1E2(ref_dst, workOrder3DT->fullkspace_, 0, RO-1, startE1_, endE1_, startE2_, endE2_));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->fullkspace_, debugFolder_+"fullkspace_After"); }
            }

            if ( performTiming_ ) { gt_timer2_.stop(); }
        }

        // partial fourier handling
        if ( partial_fourier_handling_ )
        {
            GADGET_CHECK_RETURN_FALSE(this->performPartialFourierHandling(workOrder3DT));
        }

        if ( this->computeKSpace(workOrder3DT) || fullres_coilmap )
        {
            if ( performTiming_ ) { gt_timer2_.start("full res coil map : allocate buffer 3DT ...  "); }
            hoNDArray<T> buffer3DT(workOrder3DT->fullkspace_.get_dimensions());
            hoNDArray<T> buffer3DT_Two(workOrder3DT->fullkspace_.get_dimensions());
            if ( performTiming_ ) { gt_timer2_.stop(); }

            if ( performTiming_ ) { gt_timer2_.start("full res coil map : go to image domain ...  "); }
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->fullkspace_, buffer3DT, buffer3DT_Two);
            if ( performTiming_ ) { gt_timer2_.stop(); }
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer3DT, debugFolder_+"ComplexIm_afterRefFill"); }

            if ( averageallN_coilmap )
            {
                if ( workOrder3DT->workFlow_use_BufferedKernel_ )
                {
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_));
                    Gadgetron::coil_combine(buffer3DT, *workOrder3DT->coilMap_, 3, workOrder3DT->complexIm_);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"complexImCombined"); }
                }
                else
                {
                    if ( performTiming_ ) { gt_timer2_.start("full res coil map : allocate coil map ...  "); }
                    workOrder3DT->coilMap_->create(RO, E1, E2, dstCHA, 1);
                    if ( performTiming_ ) { gt_timer2_.stop(); }

                    if ( N > 1 )
                    {
                        hoNDArray<T> aveComplexIm(RO, E1, E2, dstCHA, 1);
                        GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace5D(buffer3DT, aveComplexIm));

                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(aveComplexIm, debugFolder_+"aveComplexIm"); }

                        if ( performTiming_ ) { gt_timer2_.start("full res coil map : compute 3D coil map ...  "); }
                        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(aveComplexIm, *workOrder3DT->coilMap_, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                        if ( performTiming_ ) { gt_timer2_.stop(); }
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->coilMap_, debugFolder_+"coilMap_fullres"); }
                    }
                    else
                    {
                        if ( performTiming_ ) { gt_timer2_.start("full res coil map : compute 3D coil map ...  "); }
                        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                        if ( performTiming_ ) { gt_timer2_.stop(); }
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->coilMap_, debugFolder_+"coilMap_fullres"); }
                    }

                    if ( performTiming_ ) { gt_timer2_.start("full res coil map : coil combine 3D ...  "); }
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_));
                    Gadgetron::coil_combine(buffer3DT, *workOrder3DT->coilMap_, 3, workOrder3DT->complexIm_);
                    if ( performTiming_ ) { gt_timer2_.stop(); }
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"complexImCombined"); }
                }
            }
            else
            {
                if ( workOrder3DT->workFlow_use_BufferedKernel_ )
                {
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_));
                    Gadgetron::coil_combine(buffer3DT, *workOrder3DT->coilMap_, 3, workOrder3DT->complexIm_);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"complexIm_"); }
                }
                else
                {
                    workOrder3DT->coilMap_->create(RO, E1, E2, dstCHA, N);

                    if ( same_coilmap_allN )
                    {
                        hoNDArray<T> complexImN(RO, E1, E2, dstCHA, buffer3DT.begin()+whichN_coilmap*RO*E1*E2*dstCHA);
                        hoNDArray<T> coilMapN(RO, E1, E2, dstCHA, workOrder3DT->coilMap_->begin()+whichN_coilmap*RO*E1*E2*dstCHA);

                        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(complexImN, coilMapN, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                        GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder3DT->coilMap_, whichN_coilmap));
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->coilMap_, debugFolder_+"coilMap_fullres"); }
                    }
                    else
                    {
                        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap3DNIH(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->coil_map_algorithm_, workOrder3DT->csm_kSize_, workOrder3DT->csm_powermethod_num_, workOrder3DT->csm_iter_num_, (value_type)workOrder3DT->csm_iter_thres_, workOrder3DT->csm_true_3D_));
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder3DT->coilMap_, debugFolder_+"coilMap_fullres"); }
                    }

                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine3D(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_));
                    Gadgetron::coil_combine(buffer3DT, *workOrder3DT->coilMap_, 3, workOrder3DT->complexIm_);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"complexIm_"); }
                }
            }
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(workOrder3DT->complexIm_.get_size(0)==RO);
            GADGET_CHECK_RETURN_FALSE(workOrder3DT->complexIm_.get_size(1)==E1);
            GADGET_CHECK_RETURN_FALSE(workOrder3DT->complexIm_.get_size(2)==E2);

            if ( partial_fourier_handling_ )
            {
                bool partialFourierHandling = true;
                if ( (workOrder3DT->start_RO_<0 || workOrder3DT->end_RO_<0 || (workOrder3DT->end_RO_-workOrder3DT->start_RO_+1==RO) ) 
                        && (workOrder3DT->start_E1_<0 || workOrder3DT->end_E1_<0 || (workOrder3DT->end_E1_-workOrder3DT->start_E1_+1==E1) ) 
                        && (workOrder3DT->start_E2_<0 || workOrder3DT->end_E2_<0 || (workOrder3DT->end_E2_-workOrder3DT->start_E2_+1==E2) ) )
                {
                    partialFourierHandling = false;
                }

                // if the partial fourier handling is used to compute updated full kspace, the coil combination needs to be repeated
                if ( partialFourierHandling )
                {
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"complexIm_origin_noFullResCoilMap_"); }

                    if ( performTiming_ ) { gt_timer2_.start("after partial fourier handling, allocate buffer 3DT ...  "); }
                    hoNDArray<T> buffer3DT(workOrder3DT->fullkspace_.get_dimensions());
                    hoNDArray<T> buffer3DT_Two(workOrder3DT->fullkspace_.get_dimensions());
                    if ( performTiming_ ) { gt_timer2_.stop(); }

                    // if the partial fourier handling is performed on the fullkspace, an extra coil combination is needed
                    if (workOrder3DT->CalibMode_ == ISMRMRD_noacceleration)
                    {
                        hoNDArray<T> buffer3DT_Two(workOrder3DT->data_.get_dimensions());
                        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->data_, buffer3DT, buffer3DT_Two);
                        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_));
                        Gadgetron::coil_combine(buffer3DT, *workOrder3DT->coilMap_, 2, workOrder3DT->complexIm_);
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"complexIm_noFullResCoilMap_"); }
                    }
                    else if ( workOrder3DT->fullkspace_.get_number_of_elements() > 0 )
                    {
                        if (workOrder3DT->fullkspace_.get_size(3) == workOrder3DT->coilMap_->get_size(3))
                        {
                            hoNDArray<T> buffer3DT_Two(workOrder3DT->fullkspace_.get_dimensions());
                            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->fullkspace_, buffer3DT, buffer3DT_Two);
                            // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer3DT, *workOrder3DT->coilMap_, workOrder3DT->complexIm_));
                            Gadgetron::coil_combine(buffer3DT, *workOrder3DT->coilMap_, 2, workOrder3DT->complexIm_);
                            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_+"complexIm_noFullResCoilMap_"); }
                        }
                        else
                        {
                            if (workOrder3DT->fullkspace_.get_size(3) == 1)
                            {
                                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(workOrder3DT->fullkspace_, buffer3DT);
                                memcpy(workOrder3DT->complexIm_.begin(), buffer3DT.begin(), workOrder3DT->complexIm_.get_number_of_bytes());
                                if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(workOrder3DT->complexIm_, debugFolder_ + "complexIm_noFullResCoilMap__noReconKSpace_"); }
                            }
                        }
                    }
                }
            }
        }

        // flip along E2
        if ( performTiming_ ) { gt_timer2_.start("flip along E2 ...  "); }

        size_t imRO = workOrder3DT->complexIm_.get_size(0);
        size_t imE1 = workOrder3DT->complexIm_.get_size(1);
        size_t imE2 = workOrder3DT->complexIm_.get_size(2);
        size_t imCHA = workOrder3DT->complexIm_.get_size(3);

        hoNDArray<T> complexIm(workOrder3DT->complexIm_);

        T* pSrc = workOrder3DT->complexIm_.begin();
        T* pDst = complexIm.begin();

        size_t mid_RO = imRO/2;
        size_t mid_E1 = imE1/2;
        size_t mid_E2 = imE2/2;

        size_t n, cha;
        for ( n=0; n<workOrder3DT->complexIm_.get_size(4); n++ )
        {
            for ( cha=0; cha<imCHA; cha++ )
            {
                size_t offset = n*imRO*imE1*imE2*imCHA+cha*imRO*imE1*imE2;

                for ( size_t e2=0; e2<imE2; e2++ )
                {
                    size_t e2_from = 2*mid_E2-e2;
                    if ( e2_from >= imE2 ) e2_from -= imE2;

                    memcpy(pDst+offset+e2*imRO*imE1, pSrc+offset+e2_from*imRO*imE1, sizeof(T)*imRO*imE1);
                }
            }
        }
        if ( performTiming_ ) { gt_timer2_.stop(); }

        workOrder3DT->complexIm_ = complexIm;
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::afterUnwrapping(WorkOrderType* workOrder3DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::performPartialFourierHandling(WorkOrderType* workOrder3DT)
{
    try
    {
        // value_type partialFourierCompensationFactor = 1;

        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);

        /*if ( !( workOrder3DT->start_RO_<0 || workOrder3DT->end_RO_<0 || (workOrder3DT->end_RO_-workOrder3DT->start_RO_+1==RO) ) )
        {
            partialFourierCompensationFactor *= (value_type)(RO)/(value_type)(workOrder3DT->end_RO_-workOrder3DT->start_RO_+1);
        }

        if ( !( workOrder3DT->start_E1_<0 || workOrder3DT->end_E1_<0 || (workOrder3DT->end_E1_-workOrder3DT->start_E1_+1==E1) ) )
        {
            if ( workOrder3DT->end_E1_-workOrder3DT->start_E1_+1 <= E1 )
            {
                partialFourierCompensationFactor *= (value_type)(E1)/(value_type)(workOrder3DT->end_E1_-workOrder3DT->start_E1_+1);
            }
        }

        if ( !( workOrder3DT->start_E2_<0 || workOrder3DT->end_E2_<0 || (workOrder3DT->end_E2_-workOrder3DT->start_E2_+1==E2) ) )
        {
            if ( workOrder3DT->end_E2_-workOrder3DT->start_E2_+1 <= E2 )
            {
                partialFourierCompensationFactor *= (value_type)(E2)/(value_type)(workOrder3DT->end_E2_-workOrder3DT->start_E2_+1);
            }
        }

        partialFourierCompensationFactor = std::sqrt(partialFourierCompensationFactor);
        if ( performTiming_ ) { GDEBUG_STREAM("Partial fourier scaling factor : " << partialFourierCompensationFactor); }*/

        // if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING ) return true;

        if (workOrder3DT->CalibMode_ == ISMRMRD_noacceleration)
        {
            //if ( (workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING || workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER) && (std::abs(partialFourierCompensationFactor-1)>FLT_EPSILON) )
            //{
            //    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal(partialFourierCompensationFactor, workOrder3DT->data_));
            //}

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFilter(*workOrder3DT, workOrder3DT->data_));
            }

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_POCS )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierPOCSRecon(*workOrder3DT, workOrder3DT->data_));
            }

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_FENGHUANG )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFengHuangRecon(*workOrder3DT, workOrder3DT->data_));
            }
        }
        else if ( workOrder3DT->fullkspace_.get_number_of_elements() > 0 )
        {
            //if ( (workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING || workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER) && (std::abs(partialFourierCompensationFactor-1)>FLT_EPSILON) )
            //{
            //    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal(partialFourierCompensationFactor, workOrder3DT->fullkspace_));
            //}

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFilter(*workOrder3DT, workOrder3DT->fullkspace_));
            }

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_POCS )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierPOCSRecon(*workOrder3DT, workOrder3DT->fullkspace_));
            }

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_FENGHUANG )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFengHuangRecon(*workOrder3DT, workOrder3DT->fullkspace_));
            }
        }
        else
        {
            // perform partial fourier handling on the complex images after coil combination
            hoNDArray<T> kspace(workOrder3DT->complexIm_.get_dimensions());
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(workOrder3DT->complexIm_, kspace);

            //if ( (workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING || workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER) && (std::abs(partialFourierCompensationFactor-1)>FLT_EPSILON) )
            //{
            //    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal(partialFourierCompensationFactor, kspace));
            //}

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFilter(*workOrder3DT, kspace));
            }

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_POCS )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierPOCSRecon(*workOrder3DT, kspace));
            }

            if ( workOrder3DT->partialFourier_algo_ == ISMRMRD_PF_FENGHUANG )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFengHuangRecon(*workOrder3DT, kspace));
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(kspace, workOrder3DT->complexIm_);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performPartialFourierHandling(gtPlusReconworkOrder3DT<T>* workOrder3DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::performPartialFourierFilter(gtPlusReconWorkOrder3DT<T>& workOrder3DT, hoNDArray<T>& kspace)
{
    try
    {
        GDEBUG_STREAM("--> Into gt Plus 3DT partial fourier filter ... ");

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);

        // check whether partial fourier is used
        if ( (workOrder3DT.start_RO_<0 || workOrder3DT.end_RO_<0 || (workOrder3DT.end_RO_-workOrder3DT.start_RO_+1==RO) ) 
            && (workOrder3DT.start_E1_<0 || workOrder3DT.end_E1_<0 || (workOrder3DT.end_E1_-workOrder3DT.start_E1_+1==E1) )
            && (workOrder3DT.start_E2_<0 || workOrder3DT.end_E2_<0 || (workOrder3DT.end_E2_-workOrder3DT.start_E2_+1==E2) ) )
        {
            return true;
        }

        hoNDArray<T> input(kspace);

        if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(kspace, debugFolder_ + "kspace_before_PF_Filter"); }

        double filter_pf_width_RO_ = 0.15;
        double filter_pf_width_E1_ = 0.15;
        double filter_pf_width_E2_ = 0.15;

        bool filter_pf_density_comp_ = false;

        Gadgetron::partial_fourier_filter(input,
            workOrder3DT.start_RO_, workOrder3DT.end_RO_, workOrder3DT.start_E1_, workOrder3DT.end_E1_, workOrder3DT.start_E2_, workOrder3DT.end_E2_,
            filter_pf_width_RO_, filter_pf_width_E1_, filter_pf_width_E2_, filter_pf_density_comp_,
            workOrder3DT.filterRO_partialfourier_, workOrder3DT.filterE1_partialfourier_, workOrder3DT.filterE2_partialfourier_, kspace);

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_PF_Filter"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performPartialFourierFilter(gtPlusReconWorkOrder3DT<T>& workOrder3DT, hoNDArray<T>& kspace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::performPartialFourierPOCSRecon(WorkOrderType& workOrder3DT, hoNDArray<T>& kspace)
{
    try
    {
        GDEBUG_STREAM("--> Into gt Plus 3DT partial fourier POCS ... ");

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);
        size_t CHA = kspace.get_size(3);
        size_t N = kspace.get_size(4);

        // check whether partial fourier is used
        if ( (workOrder3DT.start_RO_<0 || workOrder3DT.end_RO_<0 || (workOrder3DT.end_RO_-workOrder3DT.start_RO_+1==RO) ) 
            && (workOrder3DT.start_E1_<0 || workOrder3DT.end_E1_<0 || (workOrder3DT.end_E1_-workOrder3DT.start_E1_+1==E1) )
            && (workOrder3DT.start_E2_<0 || workOrder3DT.end_E2_<0 || (workOrder3DT.end_E2_-workOrder3DT.start_E2_+1==E2) ) )
        {
            return true;
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_before_POCS"); }

        hoNDArray<T> input(kspace);

        Gadgetron::partial_fourier_POCS(input,
            workOrder3DT.start_RO_, workOrder3DT.end_RO_, workOrder3DT.start_E1_, workOrder3DT.end_E1_, workOrder3DT.start_E2_, workOrder3DT.end_E2_,
            workOrder3DT.partialFourier_POCS_transitBand_, workOrder3DT.partialFourier_POCS_transitBand_, workOrder3DT.partialFourier_POCS_transitBand_E2_,
            workOrder3DT.partialFourier_POCS_iters_, workOrder3DT.partialFourier_POCS_thres_, kspace);

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_POCS"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performPartialFourierPOCSRecon(WorkOrderType& workOrder3DT, hoNDArray<T>& kspace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::performPartialFourierFengHuangRecon(WorkOrderType& workOrder3DT, hoNDArray<T>& kspace)
{
    try
    {
        GDEBUG_STREAM("--> Into gt Plus 3DT partial fourier FengHuang ... ");

        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t E2 = kspace.get_size(2);
        size_t CHA = kspace.get_size(3);
        size_t N = kspace.get_size(4);

        // check whether partial fourier is used
        if ( (workOrder3DT.start_RO_<0 || workOrder3DT.end_RO_<0 || (workOrder3DT.end_RO_-workOrder3DT.start_RO_+1==RO) ) 
            && (workOrder3DT.start_E1_<0 || workOrder3DT.end_E1_<0 || (workOrder3DT.end_E1_-workOrder3DT.start_E1_+1==E1) )
            && (workOrder3DT.start_E2_<0 || workOrder3DT.end_E2_<0 || (workOrder3DT.end_E2_-workOrder3DT.start_E2_+1==E2) ) )
        {
            return true;
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_before_FengHuang"); }

        size_t startRO(0), endRO(RO-1);
        if ( workOrder3DT.start_RO_>=0 && workOrder3DT.end_RO_<RO )
        {
            startRO = workOrder3DT.start_RO_;
            endRO = workOrder3DT.end_RO_;
        }

        size_t startE1(0), endE1(E1-1);
        if ( workOrder3DT.start_E1_>=0 && workOrder3DT.end_E1_<E1 )
        {
            startE1 = workOrder3DT.start_E1_;
            endE1 = workOrder3DT.end_E1_;
        }

        size_t startE2(0), endE2(E2-1);
        if ( workOrder3DT.start_E2_>=0 && workOrder3DT.end_E2_<E2 )
        {
            startE2 = workOrder3DT.start_E2_;
            endE2 = workOrder3DT.end_E2_;
        }

        // compute the conjugate symmetric kspace
        hoNDArray<T> buffer3DT(kspace.get_dimensions());

        if ( performTiming_ ) { gt_timer1_.start("conjugateSymmetry3D"); }
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().conjugateSymmetry3D(kspace, buffer3DT));
        if ( performTiming_ ) { gt_timer1_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer3DT, debugFolder_+"kspaceConj_FengHuang"); }

        // find the symmetric region in the kspace
        size_t startSymRO, endSymRO;
        Gadgetron::find_symmetric_sampled_region(startRO, endRO, RO / 2, startSymRO, endSymRO);

        size_t startSymE1, endSymE1;
        Gadgetron::find_symmetric_sampled_region(startE1, endE1, E1 / 2, startSymE1, endSymE1);

        size_t startSymE2, endSymE2;
        Gadgetron::find_symmetric_sampled_region(startE2, endE2, E2 / 2, startSymE2, endSymE2);

        // the reference kspace for kernel estimation
        hoNDArray<T> src, dst;
        std::vector<size_t> start(5), size(5);

        start[0] = startSymRO;
        start[1] = startSymE1;
        start[2] = startSymE2;
        start[3] = 0;
        start[4] = 0;

        size[0] = endSymRO-startSymRO+1;
        size[1] = endSymE1-startSymE1+1;
        size[2] = endSymE2-startSymE2+1;;
        size[3] = CHA;
        size[4] = N;

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(buffer3DT, src, start, size));
        GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(kspace, dst, start, size));

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(src, debugFolder_+"src_FengHuang"); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(dst, debugFolder_+"dst_FengHuang"); }

        if ( workOrder3DT.partialFourier_FengHuang_sameKernel_allN_ )
        {
            hoNDArray<T> ave4D;
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace5D(src, ave4D));
            src = ave4D;

            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace5D(dst, ave4D));
            dst = ave4D;

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(src, debugFolder_+"src_ave4D_FengHuang"); }
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(dst, debugFolder_+"dst_ave4D_FengHuang"); }
        }

        // estimate the kernels
        ho6DArray<T> kernel; // [RO E1 E2 srcCHA 1 N]
        if ( performTiming_ ) { gt_timer1_.start("calibFengHuang"); }
        GADGET_CHECK_RETURN_FALSE(this->calibFengHuang(workOrder3DT, src, dst, kernel));
        if ( performTiming_ ) { gt_timer1_.stop(); }

        // perform the recon
        if ( workOrder3DT.partialFourier_FengHuang_transitBand_==0 )
        {
            if ( performTiming_ ) { gt_timer1_.start("performReconFangHuang"); }
            GADGET_CHECK_RETURN_FALSE(this->performReconFangHuang(workOrder3DT, buffer3DT, kspace, (int)startRO, (int)endRO, (int)startE1, (int)endE1, (int)startE2, (int)endE2, kernel));
            if ( performTiming_ ) { gt_timer1_.stop(); }
        }
        else
        {
            if ( performTiming_ ) { gt_timer1_.start("performReconFangHuang with transition band"); }

            long long tb =  (long long)workOrder3DT.partialFourier_FengHuang_transitBand_;

            long long sRO(startRO), eRO(endRO), sE1(startE1), eE1(endE1), sE2(startE2), eE2(endE2);

            if ( startRO > 0 )
            {
                startRO += tb;
                if ( startRO > RO ) startRO = 0;
            }

            if ( endRO < RO-1 )
            {
                endRO -= tb;
            }

            if ( startRO > endRO )
            {
                startRO = 0;
                endRO = RO-1;
            }

            if ( startE1 > 0 )
            {
                startE1 += tb;
                if ( startE1 > E1 ) startE1 = 0;
            }

            if ( endE1 < E1-1 )
            {
                endE1 -= tb;
            }

            if ( startE1 > endE1 )
            {
                startE1 = 0;
                endE1 = E1-1;
            }

            if ( startE2 > 0 )
            {
                startE2 += tb;
                if ( startE2 > E2 ) startE2 = 0;
            }

            if ( endE2 < E2-1 )
            {
                endE2 -= tb;
            }

            if ( startE2 > endE2 )
            {
                startE2 = 0;
                endE2 = E2-1;
            }

            hoNDArray<T> buffer3DT_partial_fourier_kspaceIter(kspace.get_dimensions());
            GADGET_CHECK_RETURN_FALSE(this->performReconFangHuang(workOrder3DT, buffer3DT, 
                    buffer3DT_partial_fourier_kspaceIter, (int)startRO, (int)endRO, (int)startE1, (int)endE1, (int)startE2, (int)endE2, kernel));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer3DT_partial_fourier_kspaceIter, debugFolder_+"kspace_FengHuang_recon"); }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_FengHuang_original"); }

            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.copyAlongROE1E2TransitionBand(kspace, buffer3DT_partial_fourier_kspaceIter, 
                    sRO, eRO, sE1, eE1, sE2, eE2, workOrder3DT.partialFourier_FengHuang_transitBand_, 
                    workOrder3DT.partialFourier_FengHuang_transitBand_, workOrder3DT.partialFourier_FengHuang_transitBand_E2_));

            kspace = buffer3DT_partial_fourier_kspaceIter;

            if ( performTiming_ ) { gt_timer1_.stop(); }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_FengHuang"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performPartialFourierFengHuangRecon(WorkOrderType& workOrder3DT, hoNDArray<T>& kspace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::calibFengHuang(WorkOrderType& workOrder3DT, const hoNDArray<T>& src, const hoNDArray<T>& dst, ho6DArray<T>& kernel)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(src.dimensions_equal(&dst));

        long long RO = (long long)src.get_size(0);
        long long E1 = (long long)src.get_size(1);
        long long E2 = (long long)src.get_size(2);
        long long srcCHA = (long long)src.get_size(3);
        long long N = (long long)src.get_size(4);

        long long kx = (long long)workOrder3DT.partialFourier_FengHuang_kSize_RO_;
        long long ky = (long long)workOrder3DT.partialFourier_FengHuang_kSize_E1_;
        long long kz = (long long)workOrder3DT.partialFourier_FengHuang_kSize_E2_;

        if ( kx%2 == 0 ) kx++;
        if ( ky%2 == 0 ) ky++;
        if ( kz%2 == 0 ) kz++;

        long long halfKx = (long long)kx/2;
        long long halfKy = (long long)ky/2;
        long long halfKz = (long long)kz/2;

        // the cross-channel kernel is not estimated
        kernel.createArray(kx, ky, kz, srcCHA, 1, N);

        long long ii=0;
        long long num = N*srcCHA;

        long long startRO = halfKx;
        long long endRO = RO - halfKx - 1;

        long long startE1 = halfKy;
        long long endE1 = E1 - halfKy - 1;

        long long startE2 = halfKz;
        long long endE2 = E2 - halfKz - 1;

        long long rowA, colA, rowB, colB;
        rowA = (endE2-startE2+1)*(endE1-startE1+1)*(endRO-startRO+1); 
        colA = kx*ky*kz;

        rowB = rowA;
        colB = 1;

        double thresReg = workOrder3DT.partialFourier_FengHuang_thresReg_;

        #ifdef USE_OMP
            omp_set_nested(1);
        #endif // USE_OMP

        #pragma omp parallel default(none) private(ii) shared(num, RO, E1, E2, srcCHA, N, kx, ky, kz, src, dst, kernel, rowA, colA, rowB, colB, startRO, endRO, startE1, endE1, startE2, endE2, halfKx, halfKy, halfKz, thresReg) if ( num > 1 ) num_threads( (int)(num<16 ? num : 16) )
        {
            hoNDArray<T> A_mem(rowA, colA);
            hoNDArray<T> B_mem(rowB, colB);
            hoNDArray<T> K_mem(colA, colB);

            hoMatrix<T> A(rowA, colA, A_mem.begin());
            hoMatrix<T> B(rowB, colB, B_mem.begin());
            hoMatrix<T> K(colA, colB, K_mem.begin());

            #pragma omp for
            for ( ii=0; ii<num; ii ++ )
            {
                ho3DArray<T> src3D(RO, E1, E2, const_cast<T*>(src.begin())+ii*RO*E1*E2);
                ho3DArray<T> dst3D(RO, E1, E2, const_cast<T*>(dst.begin())+ii*RO*E1*E2);

                long long ro, e1, e2, row(0);
                long long x, y, z;

                for ( e2=startE2; e2<=endE2; e2++ )
                {
                    for ( e1=startE1; e1<=endE1; e1++ )
                    {
                        for ( ro=startRO; ro<=endRO; ro++ )
                        {

                            size_t colInd(0);
                            for ( z=-halfKz; z<=halfKz; z++ )
                            {
                                for ( y=-halfKy; y<=halfKy; y++ )
                                {
                                    for ( x=-halfKx; x<=halfKx; x++ )
                                    {
                                        A(row, colInd++) = src3D(ro+x, e1+y, e2+z);
                                    }
                                }
                            }

                            B(row, 0) = dst3D(ro, e1, e2);

                            row++;
                        }
                    }
                }

                Gadgetron::SolveLinearSystem_Tikhonov(A, B, K, thresReg);

                memcpy(kernel.begin()+ii*kx*ky*kz, K.begin(), sizeof(T)*kx*ky*kz);
            }
        }

        #ifdef USE_OMP
            omp_set_nested(0);
        #endif // USE_OMP
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::calibFengHuang(WorkOrderType& workOrder3DT, const hoNDArray<T>& src, const hoNDArray<T>& dst, ho6DArray<T>& kernel) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::performReconFangHuang(WorkOrderType& workOrder3DT, 
                                                const hoNDArray<T>& kspaceConj, hoNDArray<T>& kspace, 
                                                int startRO, int endRO, int startE1, int endE1, 
                                                int startE2, int endE2, ho6DArray<T>& kernel)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(kspaceConj.dimensions_equal(&kspace));

        long long RO = (long long)kspace.get_size(0);
        long long E1 = (long long)kspace.get_size(1);
        long long E2 = (long long)kspace.get_size(2);
        long long CHA = (long long)kspace.get_size(3);
        long long N = (long long)kspace.get_size(4);

        long long kx = (long long)kernel.get_size(0);
        long long ky = (long long)kernel.get_size(1);
        long long kz = (long long)kernel.get_size(2);

        long long halfKx = kx/2;
        long long halfKy = ky/2;
        long long halfKz = kz/2;

        long long kerN = kernel.get_size(5);
        GADGET_CHECK_RETURN_FALSE( (kerN==1) || (kerN==N) );

        long long num = CHA*N;

        long long rowD = RO*E1*E2 - ( (endE2-startE2+1) * (endE1-startE1+1) * (endRO-startRO+1) );
        long long colD = kx*ky*kz;

        ho2DArray<long long> coeffX(colD, rowD);
        long long* pCx = coeffX.begin();

        ho2DArray<long long> coeffY(colD, rowD);
        long long* pCy = coeffY.begin();

        ho2DArray<long long> coeffZ(colD, rowD);
        long long* pCz = coeffZ.begin();

        long long ro, e1, e2;
        long long row(0);
        long long x, y, z;

        ho2DArray<long long> rowInd(3, rowD);
        long long* pRowInd = rowInd.begin();

        hoNDArray<long long> offsetX(colD);
        long long* pOffsetX = offsetX.begin();

        hoNDArray<long long> offsetY(colD);
        long long* pOffsetY = offsetY.begin();

        hoNDArray<long long> offsetZ(colD);
        long long* pOffsetZ = offsetZ.begin();

        long long colInd(0);
        for ( z=-halfKz; z<=halfKz; z++ )
        {
            for ( y=-halfKy; y<=halfKy; y++ )
            {
                for ( x=-halfKx; x<=halfKx; x++ )
                {
                    offsetX(colInd) = x;
                    offsetY(colInd) = y;
                    offsetZ(colInd) = z;
                    colInd++;
                }
            }
        }

        if ( performTiming_ ) { gt_timer3_.start("performReconFangHuang - compute coeff array"); }

        if ( performTiming_ ) { gt_timer2_.start("performReconFangHuang - compute coeff array - internal"); }

        long long* pRowIndCurr;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                for ( ro=0; ro<RO; ro++ )
                {
                    if ( (ro>=startRO) && (ro<=endRO) && (e1>=startE1) && (e1<=endE1) && (e2>=startE2) && (e2<=endE2) )
                    {
                        continue;
                    }

                    pRowIndCurr = pRowInd + row*3;

                    pRowIndCurr[0] = ro;
                    pRowIndCurr[1] = e1;
                    pRowIndCurr[2] = e2;

                    row++;
                }
            }
        }

        long long r;
        #pragma omp parallel for default(none) private(r) shared(rowD, colD, pCx, pCy, pCz, pRowInd, pRowIndCurr, pOffsetX, pOffsetY, pOffsetZ)
        for ( r=0; r<rowD; r++ )
        {
            long long offsetC = r*colD;
            pRowIndCurr = pRowInd + r*3;

            for ( int colInd=0; colInd<colD; colInd++ )
            {
                pCx[offsetC+colInd] = pRowIndCurr[0]+pOffsetX[colInd];
                pCy[offsetC+colInd] = pRowIndCurr[1]+pOffsetY[colInd];
                pCz[offsetC+colInd] = pRowIndCurr[2]+pOffsetZ[colInd];
            }
        }

        if ( performTiming_ ) { gt_timer2_.stop(); }

        #pragma omp parallel for default(none) private(r) shared(rowD, colD, pCx, pCy, pCz, RO, E1, E2)
        for ( r=0; r<rowD; r++ )
        {
            for ( int c=0; c<colD; c++ )
            {
                long long offset = c + r*colD;

                //pCx[offset] += pOffsetX[c];

                if ( pCx[offset] < 0 )
                {
                    pCx[offset] += RO;
                }
                else if ( pCx[offset] > RO-1 )
                {
                    pCx[offset] -= RO;
                }

                //pCy[offset] += pOffsetY[c];

                if ( pCy[offset] < 0 )
                {
                    pCy[offset] += E1;
                }
                else if ( pCy[offset] > E1-1 )
                {
                    pCy[offset] -= E1;
                }

                //pCz[offset] += pOffsetZ[c];

                if ( pCz[offset] < 0 )
                {
                    pCz[offset] += E2;
                }
                else if ( pCz[offset] > E2-1 )
                {
                    pCz[offset] -= E2;
                }
            }
        }

        /*row = 0;
        for ( e2=0; e2<E2; e2++ )
        {
            for ( e1=0; e1<E1; e1++ )
            {
                for ( ro=0; ro<RO; ro++ )
                {
                    if ( (ro>=startRO) && (ro<=endRO) && (e1>=startE1) && (e1<=endE1) && (e2>=startE2) && (e2<=endE2) )
                    {
                        continue;
                    }

                    size_t colInd(0);

                    for ( z=-halfKz; z<=halfKz; z++ )
                    {
                        dz = e2 + z;
                        if ( dz < 0 ) dz += E2;
                        if ( dz > E2-1 ) dz -= E2;

                        for ( y=-halfKy; y<=halfKy; y++ )
                        {
                            dy = e1 + y;
                            if ( dy < 0 ) dy += E1;
                            if ( dy > E1-1 ) dy -= E1;

                            for ( x=-halfKx; x<=halfKx; x++ )
                            {
                                dx = ro + x;
                                if ( dx < 0 ) dx += RO;
                                if ( dx > RO-1 ) dx -= RO;

                                coeffX(row, colInd) = dx;
                                coeffY(row, colInd) = dy;
                                coeffZ(row, colInd) = dz;
                                colInd++;
                            }
                        }
                    }

                    row++;
                }
            }
        }*/
        if ( performTiming_ ) { gt_timer3_.stop(); }

        long long ii;
        int numOfThreads = (int)((num>4) ? 4 : num);
        #pragma omp parallel default(none) private(ii) shared(num, RO, E1, E2, CHA, N, kerN, kspaceConj, kspace, kernel, rowD, colD, coeffX, coeffY, coeffZ, pCx, pCy, pCz) if ( num > 1 ) num_threads( numOfThreads )
        {
            hoNDArray<T> D_mem(rowD, colD);

            hoMatrix<T> D(rowD, colD, D_mem.begin());
            T* pD = D.begin();

            hoMatrix<T> K(colD, 1);
            hoMatrix<T> R(rowD, 1);

            Gadgetron::clear(D);
            Gadgetron::clear(K);
            Gadgetron::clear(R);

            #pragma omp for
            for ( ii=0; ii<num; ii ++ )
            {
                ho3DArray<T> src3D(RO, E1, E2, const_cast<T*>(kspaceConj.begin())+ii*RO*E1*E2);
                ho3DArray<T> dst3D(RO, E1, E2, kspace.begin()+ii*RO*E1*E2);

                long long row;

                if ( performTiming_ ) { gt_timer2_.start("fill data matrix ... "); }
                #pragma omp parallel for private(row) shared(colD, rowD, D, src3D, pD)
                for ( row=0; row<rowD; row++ )
                {
                    for ( long long col=0; col<colD; col++ )
                    {
                        long long offset = col + row*colD;
                        pD[offset] = src3D(pCx[offset], pCy[offset], pCz[offset]);
                    }
                }
                if ( performTiming_ ) { gt_timer2_.stop(); }

                if ( kerN == 1 )
                {
                    long long ind = ii;
                    long long currS = ind/(CHA*N);
                    ind %= CHA*N;
                    long long currN = ind/CHA;
                    ind %= CHA;
                    memcpy(K.begin(), kernel.begin()+(ind+currS*CHA)*colD, sizeof(T)*colD);
                }
                else
                {
                    memcpy(K.begin(), kernel.begin()+ii*colD, sizeof(T)*colD);
                }

                // R = D*K
                if ( performTiming_ ) { gt_timer2_.start("matrix multiplication ... "); }
                Gadgetron::gemm(R, D, false, K, false);
                if ( performTiming_ ) { gt_timer2_.stop(); }

                size_t colCenter = colD/2;

                if ( performTiming_ ) { gt_timer2_.start("fill the result array ... "); }
                #pragma omp parallel for private(row) default(none) shared(rowD, dst3D, colCenter, coeffX, coeffY, coeffZ, R)
                for ( row=0; row<rowD; row++ )
                {
                    dst3D( coeffX(colCenter, row), coeffY(colCenter, row), coeffZ(colCenter, row) ) = R(row, 0);
                }
                if ( performTiming_ ) { gt_timer2_.stop(); }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::performReconFangHuang(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker3DT<T>::
estimateJobSize(gtPlusReconWorkOrder<T>* workOrder3DT, size_t maxNumOfBytesPerJob, size_t overlapBetweenJobs, size_t numOfNodes, size_t& jobSize)
{
    try
    {
        size_t nodeN = numOfNodes;
        GADGET_CHECK_RETURN_FALSE(this->computeEffectiveNodeNumberBasedOnComputingPowerIndex(workOrder3DT, nodeN));
        if ( workOrder3DT->job_perform_on_control_node_ ) nodeN++;

        GDEBUG_STREAM("GtPlus Cloud 3DT - job_perform_on_control_node is " << workOrder3DT->job_perform_on_control_node_  << " - nodeN is " << nodeN << " - overlapBetweenJobs is " << overlapBetweenJobs << " ... ");

        // adjust jobN according to cloud size
        size_t RO = workOrder3DT->data_.get_size(0);
        size_t E1 = workOrder3DT->data_.get_size(1);
        size_t E2 = workOrder3DT->data_.get_size(2);
        size_t N = workOrder3DT->data_.get_size(4);

        size_t srcCHA = workOrder3DT->kernel_->get_size(3);
        size_t dstCHA = workOrder3DT->kernel_->get_size(4);

        size_t totalJobNum = RO;
        jobSize = (size_t)std::ceil( (double)(totalJobNum+overlapBetweenJobs*(nodeN-1))/(double)nodeN );

        size_t numOfBytesPerJob = sizeof(T)*( E1*E2*srcCHA*dstCHA*jobSize + 2*E1*E2*srcCHA*jobSize );

        // here a 64Mb graceful size is given to job
        while ( numOfBytesPerJob > maxNumOfBytesPerJob-64.0*1024*1024 )
        {
            nodeN *= 2;
            jobSize = (size_t)std::ceil( (double)(totalJobNum+overlapBetweenJobs*(nodeN-1))/(double)nodeN );
            numOfBytesPerJob = sizeof(T)*( E1*E2*srcCHA*dstCHA*jobSize + 2*E1*E2*srcCHA*jobSize );
        }

        GDEBUG_STREAM("GtPlus Cloud 3DT - jobSize is " << jobSize << "; every job has " << numOfBytesPerJob/1024.0/1024 << " MBytes ... ");
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker3DT<T>::estimateJobSize(...) ... ");
        return false;
    }

    return true;
}

}}
