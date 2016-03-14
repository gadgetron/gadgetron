/** \file   gtPlusISMRMRDReconWorker2DT.h
    \brief  Define the base class for the GtPlus worker for 2DT reconstruction cases

            Five different strategies were implemented for partial fourier or asymmetric echo acquisition, including:

            ISMRMRD_PF_ZEROFILLING          : only zero filling the unacquired k-space

            ISMRMRD_PF_ZEROFILLING_FILTER   : zero filling the unacquired k-space and apply a transition filter on the edges between
                                              acquired and unacquired regions

            ISMRMRD_PF_HOMODYNE             : perform the iterative homodyne filter
                                              Handbook of MRI Pulse Sequences. Page 556.
                                              Matt A. Bernstein, Kevin F. King, Xiaohong Joe Zhou. 
                                              Academic Press, ISBN-10: 0120928612.

            ISMRMRD_PF_POCS                 : perform the iterative POCS reconstruction
                                              Magnetic Resonance Imaging: Physical Principles and Sequence Design. Page 296-297.
                                              E. Mark Haacke, Robert W. Brown, Michael R. Thompson, Ramesh Venkatesan. 
                                              Wiley-Liss, ISBN-10: 0471351288.

            ISMRMRD_PF_FENGHUANG            : perform a k-space convolution based partial fourier reconstruction. 
                                              This is our recommendation for 2D, 2DT cases.

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
class gtPlusReconWorker2DT : public gtPlusReconWorker<T>
{
public:

    typedef gtPlusReconWorker<T> BaseClass;
    typedef typename realType<T>::Type value_type;

    gtPlusReconWorker2DT() : BaseClass(), startE1_(0), endE1_(1024) {}
    virtual ~gtPlusReconWorker2DT() {}

    virtual bool performRecon(gtPlusReconWorkOrder<T>* workOrder)
    {
        // check whether we have all-zeros input
        value_type v(1);
        Gadgetron::norm2(workOrder->data_, v);
        if ( v <= 0 )
        {
            GWARN_STREAM("gtPlusReconWorker2DT, performRecon(workOrder) : incoming data contains all-zeros ... ");

            boost::shared_ptr< std::vector<size_t> > dims = workOrder->data_.get_dimensions();
            (*dims)[2] = workOrder->num_channels_res_;
            workOrder->complexIm_.create(dims);
            Gadgetron::clear(workOrder->complexIm_);

            return true;
        }

        gtPlusReconWorkOrder2DT<T>* workOrder2DT = dynamic_cast<gtPlusReconWorkOrder2DT<T>*>(workOrder);
        if ( workOrder2DT == NULL ) return false;

        if ( workOrder2DT->recon_auto_parameters_ )
        {
            this->autoReconParameter(workOrder2DT);
            GDEBUG_STREAM("Gt Plus 2DT -- automatic paramter selection ---");
            if ( !this->debugFolder_.empty() ) { workOrder2DT->print(std::cout); }
        }

        return this->performRecon(workOrder2DT);
    }

    // the common functionalities are performed here for 2DT recon
    // compute the coil compression coefficients
    // prepare the ref data array
    virtual bool performRecon(gtPlusReconWorkOrder2DT<T>* workOrder2DT);

    virtual bool estimateCoilMap(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst);
    virtual bool performCalib(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst);
    virtual bool performCalibPrep(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT);
    virtual bool performCalibImpl(const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, gtPlusReconWorkOrder2DT<T>* workOrder2DT, size_t n, size_t usedS);

    virtual bool performUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& data);

    // the partial fourier handling for the 2DT reconstruction
    // the computation is performed on the reconstructed full kspace
    virtual bool performPartialFourierHandling(gtPlusReconWorkOrder2DT<T>* workOrder2DT);

    // perform the kspace filter on ref data for coil map estimation
    virtual bool performRefFilter(gtPlusReconWorkOrder2DT<T>* workOrder2DT, 
                                        const hoNDArray<T>& ref, hoNDArray<T>& refFiltered, 
                                        int startRO, int endRO, int startE1, int endE1);

    // for interleave, compute mean ref
    // for embedded and separate, squeeze out the zero lines
    virtual bool prepRef(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref, 
                    hoNDArray<T>& refRecon, hoNDArray<T>& refCoilMap, 
                    int startRO, int endRO, int startE1, int endE1, size_t dataE1);

    // implement reference data preparation
    virtual bool prepRefByAveragingCrossN(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref, bool averageAllRef, int numOfModes, hoNDArray<T>& refRecon);

    // compute coil compression coefficients
    virtual bool coilCompression(gtPlusReconWorkOrder2DT<T>* workOrder2DT);

    // after unwrapping, for embedded and separate, the full res coil map may be estimated
    // for embedded, the ref may be filled back to fullkspace
    virtual bool afterUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT);

    // pick the frame with highest signal from the 2DT buffer
    // data: [RO E1 CHA N S], res: [RO E1 CHA 1 S]
    bool pickHighestSignalForN(const hoNDArray<T>& data, hoNDArray<T>& res);

    // ----------------------------------------------------
    // common functions for 2DT reconstruction
    // ----------------------------------------------------
    // image domain kernel with coil sensitivity
    // kerIm: [RO E1 srcCHA dstCHA]
    // coilMap: [RO E1 dstCHA]
    // unmixCoeff: [RO E1 srcCHA]
    // gFactor: [RO E1]
    bool unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor);

    // apply image domain kernel
    // kspace: [RO E1 srcCHA ...]
    // complexIm : [RO E1 dstCHA ...]
    bool applyImageDomainKernel(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm);
    // aliasedIm : [RO E1 srcCHA ...]
    bool applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm);
    // for speed, a buffer can be provided
    bool applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& kerImBuffer, hoNDArray<T>& complexIm);

    // apply unmixCoeff
    // kspace: [RO E1 srcCHA ...]
    // unmixCoeff : [RO E1 srcCHA]
    // complexIm : [RO E1 ...]
    bool applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);
    // aliasedIm : [RO E1 srcCHA ...]
    bool applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm);

    // apply the KL coil compression coefficients
    bool appy_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, const std::vector< hoNDKLT<T> >& coeff, hoNDArray<T>& dataEigen);

    // ----------------------------------------------------
    // Partial fourier handling for 2DT reconstruction
    // ----------------------------------------------------
    // apply the partial fourier filer along the edges
    bool performPartialFourierFilter(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace);
    // apply the iterative homodyne filter for partial fourier reconstruction
    bool performPartialFourierHomodyneRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace);
    // apply the iterative POCS for partial fourier reconstruction
    bool performPartialFourierPOCSRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace);
    // apply the Feng Huang partial fourier reconstruction
    bool performPartialFourierFengHuangRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace);

    // compute Feng Huang kernel and perform recon
    bool calibFengHuang(gtPlusReconWorkOrder2DT<T>& workOrder2DT, const hoNDArray<T>& src, const hoNDArray<T>& dst, ho6DArray<T>& kernel);
    bool performReconFangHuang(gtPlusReconWorkOrder2DT<T>& workOrder2DT, const hoNDArray<T>& kspaceConj, hoNDArray<T>& kspace, int startRO, int endRO, int startE1, int endE1, ho6DArray<T>& kernel);

    // estimate the job size, given the maximal memory usage for every job
    virtual bool estimateJobSize(gtPlusReconWorkOrder<T>* workOrder, size_t maxNumOfBytesPerJob, size_t overlapBetweenJobs, size_t numOfNodes, size_t& jobSize);

    using BaseClass::partial_fourier_handling_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::verbose_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;
    using BaseClass::gtPlus_util_;
    using BaseClass::gtPlus_util_cplx_;

protected:

    // helper memory for computation
    hoNDArray<T> buffer2DT_;
    hoNDArray<T> buffer2DT_unwrapping_;
    hoNDArray<T> buffer2DT_partial_fourier_;
    hoNDArray<T> buffer2DT_partial_fourier_kspaceIter_;
    hoNDArray<T> ref_src_;
    hoNDArray<T> ref_dst_;
    hoNDArray<T> data_dst_;
    hoNDArray<T> ref_coil_map_dst_;

    // sampled region along E1
    size_t startE1_;
    size_t endE1_;
};

template <typename T> 
bool gtPlusReconWorker2DT<T>::performRefFilter(gtPlusReconWorkOrder2DT<T>* workOrder2DT, 
                                        const hoNDArray<T>& ref, hoNDArray<T>& refFiltered, 
                                        int startRO, int endRO, int startE1, int endE1)
{
    try
    {
        refFiltered = ref;

        size_t RO = ref.get_size(0);
        size_t E1 = ref.get_size(1);

        if ( workOrder2DT->filterROE1_ref_.get_size(0)==RO && workOrder2DT->filterROE1_ref_.get_size(1)==E1 )
        {
            Gadgetron::apply_kspace_filter_ROE1(ref, workOrder2DT->filterROE1_ref_, refFiltered);
        }
        else if ( (workOrder2DT->filterRO_ref_.get_number_of_elements()==RO) && (workOrder2DT->filterE1_ref_.get_number_of_elements()==E1) )
        {
            Gadgetron::apply_kspace_filter_ROE1(ref, workOrder2DT->filterRO_ref_, workOrder2DT->filterE1_ref_, refFiltered);
        }
        else
        {
            if ( (workOrder2DT->filterRO_ref_.get_number_of_elements()==RO) && (workOrder2DT->filterE1_ref_.get_number_of_elements()!=E1) )
            {
                Gadgetron::apply_kspace_filter_RO(ref, workOrder2DT->filterRO_ref_, refFiltered);
            }

            if ( (workOrder2DT->filterRO_ref_.get_number_of_elements()!=RO) && (workOrder2DT->filterE1_ref_.get_number_of_elements()==E1) )
            {
                Gadgetron::apply_kspace_filter_E1(ref, workOrder2DT->filterE1_ref_, refFiltered);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performRefFilter(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::prepRefByAveragingCrossN(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref, bool averageAllRef, int numOfModes, hoNDArray<T>& refRecon)
{
    try
    {
        size_t RO = ref.get_size(0);
        size_t E1 = ref.get_size(1);
        size_t CHA = ref.get_size(2);
        size_t N = ref.get_size(3);
        size_t S = ref.get_size(4);

        std::vector<size_t> sampledTimes;

        if ( !averageAllRef && ( (numOfModes<1) || (numOfModes>N-1) ) )
        {
            refRecon = ref;
        }
        else if ( averageAllRef && ( (numOfModes<1) || (numOfModes>N-1) ) )
        {
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(ref, refRecon, sampledTimes));
        }
        else if ( averageAllRef && (numOfModes>=1) && (numOfModes<=N-1) )
        {
            hoNDArray<T> refKLF(RO, E1, CHA, N, S);

            size_t s;
            for ( s=0; s<S; s++ )
            {
                hoMatrix<T> A(RO*E1*CHA, N, const_cast<T*>(ref.begin()+s*RO*E1*CHA*N));
                hoMatrix<T> A_KLF(RO*E1*CHA, N, refKLF.begin()+s*RO*E1*CHA*N);

                hoNDKLT<T> klt;
                klt.prepare(A, (size_t)1, (size_t)0);
                klt.KL_filter(A, A_KLF, 1, numOfModes);
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refKLF, debugFolder_+"refKLF"); }

            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(refKLF, refRecon, sampledTimes));
        }
        else if ( !averageAllRef && (numOfModes>=1) && (numOfModes<=N-1) )
        {
            refRecon.create(RO, E1, CHA, N, S);

            size_t s;
            for ( s=0; s<S; s++ )
            {
                hoMatrix<T> A(RO*E1*CHA, N, const_cast<T*>(ref.begin()+s*RO*E1*CHA*N));
                hoMatrix<T> A_KLF(RO*E1*CHA, N, refRecon.begin()+s*RO*E1*CHA*N);

                hoNDKLT<T> klt;
                klt.prepare(A, (size_t)1, (size_t)0);
                klt.KL_filter(A, A_KLF, 1, numOfModes);
            }
        }
        else
        {
            refRecon = ref;
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::prepRefByAveragingCrossN(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::prepRef(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref, 
                    hoNDArray<T>& refRecon, hoNDArray<T>& refCoilMap, 
                    int startRO, int endRO, int startE1, int endE1, size_t dataE1)
{
    try
    {
        size_t dataRO = workOrder2DT->data_.get_size(0);
        size_t dataS = workOrder2DT->data_.get_size(4);

        size_t RO = ref.get_size(0);
        size_t E1 = ref.get_size(1);
        size_t srcCHA = ref.get_size(2);
        size_t N = ref.get_size(3);
        size_t S = ref.get_size(4);

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref, debugFolder_+"ref"); }

        if (workOrder2DT->CalibMode_ == ISMRMRD_noacceleration)
        {
            if ( workOrder2DT->no_acceleration_averageall_ref_ )
            {
                GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder2DT, ref, workOrder2DT->no_acceleration_averageall_ref_, workOrder2DT->no_acceleration_ref_numOfModes_, refRecon));
            }

            GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder2DT, refRecon, refCoilMap, startRO, endRO, startE1, endE1));
        }
        else if ( workOrder2DT->CalibMode_ == ISMRMRD_interleaved )
        {
            GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder2DT, ref, true, workOrder2DT->interleaved_ref_numOfModes_, refRecon));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_+"refRecon_interleaved"); }

            GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder2DT, refRecon, refCoilMap, startRO, endRO, startE1, endE1));

            if ( (startRO>=0 && endRO>0 && endRO>startRO) || (startE1>=0 && endE1>0 && endE1>startE1) )
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

                hoNDArray<T> croppedRef;
                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                refRecon = croppedRef;
            }
        }
        else if ( workOrder2DT->CalibMode_ == ISMRMRD_embedded 
                || workOrder2DT->CalibMode_ == ISMRMRD_separate 
                || workOrder2DT->CalibMode_ == ISMRMRD_external )
        {
            if ( workOrder2DT->CalibMode_ == ISMRMRD_embedded )
            {
                refRecon = ref;
            }

            if ( workOrder2DT->CalibMode_ == ISMRMRD_separate )
            {
                GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder2DT, ref, workOrder2DT->separate_averageall_ref_, workOrder2DT->separate_ref_numOfModes_, refRecon));
            }

            hoNDArray<typename realType<T>::Type> refMag(refRecon.get_dimensions()), refMagSum;
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::abs(refRecon, refMag));

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(sum_over_dimension(refMag, refMagSum, refMag.get_number_of_dimensions()-1));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(sum_over_dimension(refMagSum, refMag, refMagSum.get_number_of_dimensions() - 2));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(sum_over_dimension(refMag, refMagSum, refMag.get_number_of_dimensions() - 3));

            refMagSum.squeeze();
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<float>().detectSampledRegionE1(refMagSum, startE1_, endE1_));

            std::vector<size_t> crop_offset(5);
            crop_offset[0] = 0;
            crop_offset[1] = startE1_;
            crop_offset[2] = 0;
            crop_offset[3] = 0;
            crop_offset[4] = 0;

            std::vector<size_t> crop_size(5);
            crop_size[0] = refRecon.get_size(0);
            crop_size[1] = endE1_-startE1_+1;
            crop_size[2] = srcCHA;
            crop_size[3] = refRecon.get_size(3);
            crop_size[4] = refRecon.get_size(4);

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_+"refRecon_beforeCrop"); }

            if ( workOrder2DT->CalibMode_ == ISMRMRD_embedded )
            {
                hoNDArray<T> croppedRef;
                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(croppedRef, debugFolder_+"refRecon_afterCrop"); }

                if ( workOrder2DT->recon_algorithm_ == ISMRMRD_SPIRIT 
                    || workOrder2DT->recon_algorithm_ == ISMRMRD_L1SPIRIT 
                    || workOrder2DT->recon_algorithm_ == ISMRMRD_L1SPIRIT_SLEP 
                    || workOrder2DT->recon_algorithm_ == ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP )
                {
                    // copy the ref into the data
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.copyAlongE1(refRecon, workOrder2DT->data_, startE1_, endE1_));
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->data_, debugFolder_+"data_copyAlongE1"); }
                }

                GADGET_CHECK_RETURN_FALSE(prepRefByAveragingCrossN(workOrder2DT, croppedRef, workOrder2DT->embedded_averageall_ref_, workOrder2DT->embedded_ref_numOfModes_, refRecon));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_+"refRecon_afterCrop_prepCrossN"); }

                crop_size[3] = refRecon.get_size(3);

                refCoilMap.create(RO, E1, srcCHA, refRecon.get_size(3), S);
                GADGET_CHECK_RETURN_FALSE(setSubArrayUpTo11DArray(refRecon, refCoilMap, crop_offset, crop_size));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"refCoilMap"); }

                hoNDArray<T> refCoilMapTmp(refCoilMap);
                GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder2DT, refCoilMapTmp, refCoilMap, startRO, endRO, startE1, endE1));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"refCoilMap_filtered"); }

                if ( refRecon.get_size(0) == RO )
                {
                    if ( startRO>=0 && endRO>0 && endRO>startRO && startRO<RO && endRO<RO )
                    {
                        crop_offset[0] = startRO;
                        crop_size[0] = endRO-startRO+1;

                        crop_offset[1] = 0;
                        crop_size[1] = refRecon.get_size(1);
                    }
                }

                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                refRecon = croppedRef;
            }
            else
            {
                hoNDArray<T> croppedRef;
                GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(croppedRef, debugFolder_+"croppedRef"); }

                GADGET_CHECK_RETURN_FALSE(performRefFilter(workOrder2DT, croppedRef, refCoilMap, startRO, endRO, startE1, endE1));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"croppedRef_filtered"); }

                refRecon = croppedRef;

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::pad(dataRO, dataE1, &refCoilMap, &croppedRef));
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

                        GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(refRecon, croppedRef, crop_offset, crop_size));
                        refRecon = croppedRef;
                    }
                }
            }

            if ( S < dataS )
            {
                hoNDArray<T> refReconDataS(refRecon.get_size(0), refRecon.get_size(1), refRecon.get_size(2), refRecon.get_size(3), dataS);
                hoNDArray<T> refCoilMapDataS(refCoilMap.get_size(0), refCoilMap.get_size(1), refCoilMap.get_size(2), refCoilMap.get_size(3), dataS);

                memcpy(refReconDataS.begin(), refRecon.begin(), refRecon.get_number_of_bytes());
                memcpy(refCoilMapDataS.begin(), refCoilMap.begin(), refCoilMap.get_number_of_bytes());

                size_t refReconN4D = refRecon.get_size(0)*refRecon.get_size(1)*refRecon.get_size(2)*refRecon.get_size(3);
                size_t refCoilMapN4D = refCoilMap.get_size(0)*refCoilMap.get_size(1)*refCoilMap.get_size(2)*refCoilMap.get_size(3);

                size_t s;
                for ( s=S; s<dataS; s++ )
                {
                    memcpy(refReconDataS.begin()+s*refReconN4D, refRecon.begin()+(S-1)*refReconN4D, sizeof(T)*refReconN4D);
                    memcpy(refCoilMapDataS.begin()+s*refCoilMapN4D, refCoilMap.begin()+(S-1)*refCoilMapN4D, sizeof(T)*refCoilMapN4D);
                }

                refRecon = refReconDataS;
                refCoilMap = refCoilMapDataS;
            }
        }
        else
        {
            GERROR_STREAM("CalibMode is not supported in gtPlusReconWorker2DT<T>::prepRef(...) : " << workOrder2DT->CalibMode_);
            return false;
        }

        // if the upstream coil compression is needed
        if ( workOrder2DT->upstream_coil_compression_ )
        {
            if ( !debugFolder_.empty() ) { GDEBUG_STREAM("Upstream coil compression ... "); }

            std::vector<hoNDKLT<T> > upstreamKLTRef(workOrder2DT->ref_.get_size(4)), upstreamKLTRefRecon(refRecon.get_size(4));
            std::vector<hoNDKLT<T> > upstreamKLTData(workOrder2DT->data_.get_size(4));

            if ( workOrder2DT->same_coil_compression_coeff_allS_ )
            {
                hoNDArray<T> aveAllS;

                std::vector<size_t> allSDim(4);
                allSDim[0] = refRecon.get_size(0);
                allSDim[1] = refRecon.get_size(1);
                allSDim[2] = refRecon.get_size(2);
                allSDim[3] = refRecon.get_size(3)*refRecon.get_size(4);

                hoNDArray<T> dataAllS(&allSDim, refRecon.begin(), false);
                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(dataAllS, aveAllS));

                hoNDArray<T> eigenValues;
                if ( workOrder2DT->coil_compression_num_modesKept_ > 0 )
                {
                    upstreamKLTData[0].prepare(aveAllS, 2, (size_t)workOrder2DT->upstream_coil_compression_num_modesKept_);
                    upstreamKLTData[0].eigen_value(eigenValues);
                }
                else
                {
                    upstreamKLTData[0].prepare(aveAllS, 2, (value_type)(workOrder2DT->upstream_coil_compression_thres_));
                    upstreamKLTData[0].eigen_value(eigenValues);
                }

                {
                    eigenValues.print(std::cout);
                    for (size_t i = 0; i < eigenValues.get_size(0); i++)
                    {
                        GDEBUG_STREAM(i << " = " << eigenValues(i));
                    }
                }
                GDEBUG_STREAM("Upstream coil compression, number of channel kept is " << upstreamKLTData[0].output_length());

                size_t n;
                for (n = 1; n<upstreamKLTRef.size(); n++)
                {
                    upstreamKLTRef[n] = upstreamKLTData[0];
                }

                for (n = 1; n<upstreamKLTRefRecon.size(); n++)
                {
                    upstreamKLTRefRecon[n] = upstreamKLTData[0];
                }

                for (n = 1; n<upstreamKLTData.size(); n++)
                {
                    upstreamKLTData[n] = upstreamKLTData[0];
                }
            }
            else
            {
                std::vector<size_t> allSDim(4);
                allSDim[0] = refRecon.get_size(0);
                allSDim[1] = refRecon.get_size(1);
                allSDim[2] = refRecon.get_size(2);
                allSDim[3] = refRecon.get_size(3);

                size_t N_refRecon = allSDim[0]*allSDim[1]*allSDim[2]*allSDim[3];

                size_t num_modesKept = srcCHA;

                size_t s;
                for ( s=0; s<refRecon.get_size(4); s++ )
                {
                    hoNDArray<T> dataCurrS(&allSDim, refRecon.begin()+s*N_refRecon, false);

                    hoNDArray<T> eigenValues;

                    if ( s == 0 )
                    {
                        if ( workOrder2DT->coil_compression_num_modesKept_ > 0 )
                        {
                            upstreamKLTData[0].prepare(dataCurrS, 2, (size_t)workOrder2DT->upstream_coil_compression_num_modesKept_);
                            upstreamKLTData[0].eigen_value(eigenValues);
                        }
                        else
                        {
                            upstreamKLTData[0].prepare(dataCurrS, 2, (value_type)(workOrder2DT->upstream_coil_compression_thres_));
                            upstreamKLTData[0].eigen_value(eigenValues);
                        }

                        num_modesKept = upstreamKLTData[0].output_length();
                    }
                    else
                    {
                        upstreamKLTData[s].prepare(dataCurrS, 2, (size_t)num_modesKept);
                        upstreamKLTData[s].eigen_value(eigenValues);
                    }

                    if ( !debugFolder_.empty() )
                    {
                        eigenValues.print(std::cout);
                        for (size_t i = 0; i<eigenValues.get_size(0); i++)
                        {
                            GDEBUG_STREAM(i << " = " << eigenValues(i));
                        }
                    }
                    GDEBUG_STREAM("Upstream coil compression, number of channel kept is " << num_modesKept);

                    if ( s < upstreamKLTRef.size() )
                    {
                        upstreamKLTRef[s] = upstreamKLTData[s];
                    }

                    upstreamKLTRefRecon[s] = upstreamKLTData[s];
                }
            }

            // apply the coil compression
            #ifdef USE_OMP
                omp_set_nested(1);
            #endif // USE_OMP

            if ( performTiming_ ) { gt_timer2_.start("apply upstream coil compression ... "); }
            #pragma omp parallel sections default(shared)
            {

                #pragma omp section
                {
                    if ( performTiming_ ) { gt_timer3_.start("appy_KLT_coil_compression_coeff_2D ... "); }
                    this->appy_KLT_coil_compression_coeff_2D(workOrder2DT->data_, upstreamKLTData, data_dst_);
                    if ( performTiming_ ) { gt_timer3_.stop(); }

                    if ( performTiming_ ) { gt_timer3_.start("copy data ... "); }
                    workOrder2DT->data_ = data_dst_;
                    if ( performTiming_ ) { gt_timer3_.stop(); }
                }

                #pragma omp section
                {
                    this->appy_KLT_coil_compression_coeff_2D(workOrder2DT->ref_, upstreamKLTRef, ref_dst_);
                    workOrder2DT->ref_ = ref_dst_;
                }

                #pragma omp section
                {
                    hoNDArray<T> refRecon_upstream;
                    this->appy_KLT_coil_compression_coeff_2D(refRecon, upstreamKLTRefRecon, refRecon_upstream);
                    refRecon = refRecon_upstream;
                    refRecon_upstream.clear();
                }

                #pragma omp section
                {
                    hoNDArray<T> refCoilMap_upstream;
                    this->appy_KLT_coil_compression_coeff_2D(refCoilMap, upstreamKLTRefRecon, refCoilMap_upstream);
                    refCoilMap = refCoilMap_upstream;
                    refCoilMap_upstream.clear();
                }
            }

            if ( performTiming_ ) { gt_timer2_.stop(); }

            #ifdef USE_OMP
                omp_set_nested(0);
            #endif // USE_OMP
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refRecon, debugFolder_+"refRecon"); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(refCoilMap, debugFolder_+"refCoilMap"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::prepRef(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::coilCompression(gtPlusReconWorkOrder2DT<T>* workOrder2DT)
{
    // the 2DT recon on 5D array [RO E1 CHA N S]
    try
    {
        size_t RO = workOrder2DT->ref_recon_.get_size(0);
        size_t E1 = workOrder2DT->ref_recon_.get_size(1);
        size_t srcCHA = workOrder2DT->ref_recon_.get_size(2);
        size_t N = workOrder2DT->ref_recon_.get_size(3);
        size_t S = workOrder2DT->ref_recon_.get_size(4);

        size_t dataS = workOrder2DT->data_.get_size(4);

        // if ( workOrder2DT->acceFactorE1_ == 1 ) return true;

        // compute coil compression coeff
        if ( workOrder2DT->coil_compression_ )
        {
            // check whether coil compression coeff has been preset
            if ( workOrder2DT->coilCompressionCoef_->size()!=S )
            {
                if ( workOrder2DT->same_coil_compression_coeff_allS_ )
                {
                    hoNDArray<T> aveAllS;

                    std::vector<size_t> allSDim(4);
                    allSDim[0] = RO;
                    allSDim[1] = E1;
                    allSDim[2] = srcCHA;
                    allSDim[3] = N*S;

                    hoNDArray<T> dataAllS(&allSDim, workOrder2DT->ref_recon_.begin(), false);
                    GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(dataAllS, aveAllS));

                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(aveAllS, debugFolder_+"aveAllS"); }

                    workOrder2DT->coilCompressionCoef_->resize(dataS);

                    hoNDArray<T> eigenValues;
                    if ( workOrder2DT->coil_compression_num_modesKept_ > 0 )
                    {
                        (*workOrder2DT->coilCompressionCoef_)[0].prepare(aveAllS, 2, (size_t)workOrder2DT->coil_compression_num_modesKept_);
                        (*workOrder2DT->coilCompressionCoef_)[0].eigen_value(eigenValues);
                    }
                    else
                    {
                        (*workOrder2DT->coilCompressionCoef_)[0].prepare(aveAllS, 2, (value_type)workOrder2DT->coil_compression_thres_);
                        (*workOrder2DT->coilCompressionCoef_)[0].eigen_value(eigenValues);
                    }

                    size_t s;
                    for ( s=1; s<dataS; s++ )
                    {
                        (*workOrder2DT->coilCompressionCoef_)[s] = (*workOrder2DT->coilCompressionCoef_)[0];
                    }

                    if (!debugFolder_.empty())
                    {
                        eigenValues.print(std::cout);
                        for (size_t i = 0; i<eigenValues.get_size(0); i++)
                        {
                            GDEBUG_STREAM(i << " = " << eigenValues(i));
                        }
                    }
                    GDEBUG_STREAM("Coil compression, number of channel kept is " << (*workOrder2DT->coilCompressionCoef_)[0].output_length());
                }
                else
                {
                    std::vector<size_t> allSDim(4);
                    allSDim[0] = RO;
                    allSDim[1] = E1;
                    allSDim[2] = srcCHA;
                    allSDim[3] = N;

                    size_t num_modesKept = srcCHA;

                    workOrder2DT->coilCompressionCoef_->resize(S);

                    size_t s;
                    for ( s=0; s<S; s++ )
                    {
                        hoNDArray<T> dataCurrS(&allSDim, workOrder2DT->ref_recon_.begin()+s*RO*E1*srcCHA*N, false);

                        hoNDArray<T> eigenValues;

                        if ( s == 0 )
                        {
                            if ( workOrder2DT->coil_compression_num_modesKept_ > 0 )
                            {
                                (*workOrder2DT->coilCompressionCoef_)[0].prepare(dataCurrS, 2, (size_t)workOrder2DT->coil_compression_num_modesKept_);
                                (*workOrder2DT->coilCompressionCoef_)[0].eigen_value(eigenValues);
                            }
                            else
                            {
                                (*workOrder2DT->coilCompressionCoef_)[0].prepare(dataCurrS, 2, (value_type)workOrder2DT->coil_compression_thres_);
                                (*workOrder2DT->coilCompressionCoef_)[0].eigen_value(eigenValues);
                            }

                            num_modesKept = (*workOrder2DT->coilCompressionCoef_)[0].output_length();
                        }
                        else
                        {
                            (*workOrder2DT->coilCompressionCoef_)[s].prepare(dataCurrS, 2, (size_t)num_modesKept);
                            (*workOrder2DT->coilCompressionCoef_)[s].eigen_value(eigenValues);
                        }

                        if (!debugFolder_.empty())
                        {
                            eigenValues.print(std::cout);
                            for (size_t i = 0; i<eigenValues.get_size(0); i++)
                            {
                                GDEBUG_STREAM(i << " = " << eigenValues(i));
                            }
                        }
                        GDEBUG_STREAM("Coil compression, number of channel kept is " << (*workOrder2DT->coilCompressionCoef_)[s].output_length());
                    }

                    if ( S < dataS )
                    {
                        std::vector<hoNDKLT<T> > coilCompressionCoef(dataS);
                        for ( s=0; s<S; s++ )
                        {
                            coilCompressionCoef[s] = (*workOrder2DT->coilCompressionCoef_)[s];
                        }

                        for ( s=S; s<dataS; s++ )
                        {
                            coilCompressionCoef[s] = (*workOrder2DT->coilCompressionCoef_)[S-1];
                        }

                        *(workOrder2DT->coilCompressionCoef_) = coilCompressionCoef;
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::coilCompression(gtPlusReconWorkOrder2DT<T>* workOrder2DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performRecon(gtPlusReconWorkOrder2DT<T>* workOrder2DT)
{
    // the 2DT recon on 5D array [RO E1 CHA N S]
    try
    {
        if ( !workOrder2DT->workFlow_use_BufferedKernel_ )
        {
            if ( performTiming_ ) { gt_timer1_.start("prepRef"); }
            GADGET_CHECK_RETURN_FALSE(prepRef(workOrder2DT, workOrder2DT->ref_, workOrder2DT->ref_recon_, workOrder2DT->ref_coil_map_, 
                        workOrder2DT->start_RO_, workOrder2DT->end_RO_, workOrder2DT->start_E1_, workOrder2DT->end_E1_, workOrder2DT->data_.get_size(1)));
            if ( performTiming_ ) { gt_timer1_.stop(); }

            if ( performTiming_ ) { gt_timer1_.start("coilCompression"); }
            GADGET_CHECK_RETURN_FALSE(coilCompression(workOrder2DT));
            if ( performTiming_ ) { gt_timer1_.stop(); }
        }

         // apply coil compression coefficients
        if ( workOrder2DT->workFlow_use_BufferedKernel_ )
        {
            if ( workOrder2DT->coil_compression_ )
            {
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->data_, debugFolder_+"data_"); }
                this->appy_KLT_coil_compression_coeff_2D(workOrder2DT->data_, *workOrder2DT->coilCompressionCoef_, data_dst_);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(data_dst_, debugFolder_+"data_dst_"); }
            }
            else
            {
                data_dst_ = workOrder2DT->data_;
            }
        }
        else
        {
            if ( workOrder2DT->coil_compression_ )
            {
                ref_src_ = workOrder2DT->ref_recon_;

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_src_, debugFolder_+"ref_src_"); }
                this->appy_KLT_coil_compression_coeff_2D(ref_src_, *workOrder2DT->coilCompressionCoef_, ref_dst_);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_dst_, debugFolder_+"ref_dst_"); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->data_, debugFolder_+"data_"); }
                this->appy_KLT_coil_compression_coeff_2D(workOrder2DT->data_, *workOrder2DT->coilCompressionCoef_, data_dst_);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(data_dst_, debugFolder_+"data_dst_"); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->ref_coil_map_, debugFolder_+"ref_coil_map_"); }
                this->appy_KLT_coil_compression_coeff_2D(workOrder2DT->ref_coil_map_, *workOrder2DT->coilCompressionCoef_, ref_coil_map_dst_);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_coil_map_dst_, debugFolder_+"ref_coil_map_dst_"); }

                if ( !workOrder2DT->downstream_coil_compression_ 
                    || workOrder2DT->recon_algorithm_==ISMRMRD_SPIRIT 
                    || workOrder2DT->recon_algorithm_==ISMRMRD_L1SPIRIT 
                    || workOrder2DT->recon_algorithm_==ISMRMRD_L1SPIRIT_SLEP 
                    || workOrder2DT->recon_algorithm_==ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP )
                {
                    ref_src_ = ref_dst_;
                }
            }
            else
            {
                ref_src_ = workOrder2DT->ref_recon_;
                ref_dst_ = workOrder2DT->ref_recon_;
                data_dst_ = workOrder2DT->data_;
                ref_coil_map_dst_ = workOrder2DT->ref_coil_map_;
            }

            if ( performTiming_ ) { gt_timer1_.start("estimateCoilMap"); }
            GADGET_CHECK_RETURN_FALSE(this->estimateCoilMap(workOrder2DT, ref_src_, ref_dst_, ref_coil_map_dst_));
            if ( performTiming_ ) { gt_timer1_.stop(); }

            if ( performTiming_ ) { gt_timer1_.start("performCalib"); }
            GADGET_CHECK_RETURN_FALSE(this->performCalib(workOrder2DT, ref_src_, ref_dst_, ref_coil_map_dst_));
            if ( performTiming_ ) { gt_timer1_.stop(); }
        }

        if ( performTiming_ ) { gt_timer1_.start("performUnwrapping"); }
        GADGET_CHECK_RETURN_FALSE(this->performUnwrapping(workOrder2DT, data_dst_));
        if ( performTiming_ ) { gt_timer1_.stop(); }

        if ( performTiming_ ) { gt_timer1_.start("afterUnwrapping"); }
        GADGET_CHECK_RETURN_FALSE(this->afterUnwrapping(workOrder2DT));
        if ( performTiming_ ) { gt_timer1_.stop(); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performRecon(gtPlusReconWorkOrder2DT<T>* workOrder2DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::
estimateCoilMap(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst)
{
    try
    {
        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t srcCHA = ref_src.get_size(2);

        size_t refRO = ref_dst.get_size(0);
        size_t refE1 = ref_dst.get_size(1);
        size_t refN = ref_dst.get_size(3);
        size_t dstCHA = ref_coil_map_dst.get_size(2);

        bool same_combinationcoeff_allS = false;
        size_t whichS_combinationcoeff = 0;
        if ( workOrder2DT->CalibMode_ == ISMRMRD_interleaved )
        {
            same_combinationcoeff_allS = workOrder2DT->interleaved_same_combinationcoeff_allS_;
            whichS_combinationcoeff = workOrder2DT->interleaved_whichS_combinationcoeff_;
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_embedded )
        {
            same_combinationcoeff_allS = workOrder2DT->embedded_same_combinationcoeff_allS_;
            whichS_combinationcoeff = workOrder2DT->embedded_whichS_combinationcoeff_;
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_separate )
        {
            same_combinationcoeff_allS = workOrder2DT->separate_same_combinationcoeff_allS_;
            whichS_combinationcoeff = workOrder2DT->separate_whichS_combinationcoeff_;
        }

        if ( whichS_combinationcoeff >= S ) whichS_combinationcoeff=S-1;

        // if the coil map has not been preset
        if ( (workOrder2DT->coilMap_->get_size(0)!=RO) 
            || (workOrder2DT->coilMap_->get_size(1)!=E1)
            || (workOrder2DT->coilMap_->get_size(4)!=S) )
        {
            if ( same_combinationcoeff_allS )
            {
                size_t usedS = whichS_combinationcoeff;

                hoNDArray<T> refCoilMapS(RO, E1, dstCHA, refN, const_cast<T*>(ref_coil_map_dst.begin()+usedS*RO*E1*dstCHA*refN));

                workOrder2DT->coilMap_->create(RO, E1, dstCHA, refN, S);

                hoNDArray<T> coilMapS(RO, E1, dstCHA, refN, workOrder2DT->coilMap_->begin()+usedS*RO*E1*dstCHA*refN);

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(refCoilMapS, buffer2DT_);
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(buffer2DT_, 
                        coilMapS, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, 
                        workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_));

                GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder2DT->coilMap_, usedS));
            }
            else
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(ref_coil_map_dst, buffer2DT_);
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(buffer2DT_, 
                        *workOrder2DT->coilMap_, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, 
                        workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_));
            }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder2DT->coilMap_, debugFolder_+"coilMap_"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::estimateCoilMap(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::
performCalib(gtPlusReconWorkOrder2DT<T>* workOrder2DT, const hoNDArray<T>& ref_src, const hoNDArray<T>& ref_dst, const hoNDArray<T>& ref_coil_map_dst)
{
    try
    {
        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t srcCHA = ref_src.get_size(2);

        size_t refRO = ref_dst.get_size(0);
        size_t refE1 = ref_dst.get_size(1);
        size_t refN = ref_dst.get_size(3);
        size_t dstCHA = ref_coil_map_dst.get_size(2);

        bool same_combinationcoeff_allS = false;
        size_t whichS_combinationcoeff = 0;
        if ( workOrder2DT->CalibMode_ == ISMRMRD_interleaved )
        {
            same_combinationcoeff_allS = workOrder2DT->interleaved_same_combinationcoeff_allS_;
            whichS_combinationcoeff = workOrder2DT->interleaved_whichS_combinationcoeff_;
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_embedded )
        {
            same_combinationcoeff_allS = workOrder2DT->embedded_same_combinationcoeff_allS_;
            whichS_combinationcoeff = workOrder2DT->embedded_whichS_combinationcoeff_;
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_separate )
        {
            same_combinationcoeff_allS = workOrder2DT->separate_same_combinationcoeff_allS_;
            whichS_combinationcoeff = workOrder2DT->separate_whichS_combinationcoeff_;
        }

        if ( whichS_combinationcoeff >= S ) whichS_combinationcoeff=S-1;

        // calibration
        if ( (workOrder2DT->kernelIm_->get_size(0)!=RO) 
                || (workOrder2DT->kernelIm_->get_size(1)!=E1)
                || (workOrder2DT->kernelIm_->get_size(2)!=srcCHA)
                || (workOrder2DT->kernelIm_->get_size(3)!=dstCHA)
                || (workOrder2DT->kernelIm_->get_size(5)!=S) )
        {
            GADGET_CHECK_RETURN_FALSE(this->performCalibPrep(ref_src, ref_dst, workOrder2DT));

            size_t n;

            // perform calibration
            if ( same_combinationcoeff_allS )
            {
                size_t usedS = whichS_combinationcoeff;

                for ( n=0; n<refN; n++ )
                {
                    GADGET_CHECK_RETURN_FALSE(this->performCalibImpl(ref_src, ref_dst, workOrder2DT, n, usedS));
                }

                if ( S > 1 )
                {
                    GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder2DT->kernel_, usedS));
                    GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder2DT->kernelIm_, usedS));
                    if (workOrder2DT->unmixingCoeffIm_ && workOrder2DT->unmixingCoeffIm_->get_size(4) == S) GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder2DT->unmixingCoeffIm_, usedS));
                    if ( workOrder2DT->gfactor_needed_ ) { GADGET_CHECK_RETURN_FALSE(repmatLastDimension(workOrder2DT->gfactor_, usedS)); }
                    if ( workOrder2DT->wrap_around_map_needed_ ) { GADGET_CHECK_RETURN_FALSE(repmatLastDimension(workOrder2DT->wrap_around_map_, usedS)); }
                }

                if (!debugFolder_.empty())
                {
                    gt_exporter_.exportArrayComplex(workOrder2DT->gfactor_, debugFolder_ + "gfactor_after_calib");
                }
            }
            else
            {
                int usedS;
                #ifdef USE_OMP
                    if ( S < omp_get_num_procs()/2 )
                    {
                        omp_set_nested(1);
                        GDEBUG_STREAM("performCalib, nested omp is on ... ");
                    }
                #endif // USE_OMP

                #pragma omp parallel for default(none) private(usedS) shared(S, refN, ref_src, ref_dst, workOrder2DT) if (S>1)
                for ( usedS=0; usedS<(int)S; usedS++ )
                {
                    for ( size_t n=0; n<refN; n++ )
                    {
                        this->performCalibImpl(ref_src, ref_dst, workOrder2DT, n, usedS);
                    }
                }

                #ifdef USE_OMP
                    omp_set_nested(0);
                #endif // USE_OMP

                if (!debugFolder_.empty())
                {
                    gt_exporter_.exportArrayComplex(workOrder2DT->gfactor_, debugFolder_ + "gfactor_after_calib");
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performCalib(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::
performCalibPrep(const hoNDArray<T>& , const hoNDArray<T>& , gtPlusReconWorkOrder2DT<T>* /*workOrder2DT*/)
{
    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::
performCalibImpl(const hoNDArray<T>& , const hoNDArray<T>& , gtPlusReconWorkOrder2DT<T>* , size_t , size_t )
{
    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performUnwrapping(gtPlusReconWorkOrder2DT<T>* , const hoNDArray<T>& )
{
    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t srcCHA = kerIm.get_size(2);
        size_t dstCHA = kerIm.get_size(3);

        GADGET_CHECK_RETURN_FALSE(coilMap.get_size(0)==RO);
        GADGET_CHECK_RETURN_FALSE(coilMap.get_size(1)==E1);
        GADGET_CHECK_RETURN_FALSE(coilMap.get_size(2)==dstCHA);

        unmixCoeff.create(RO, E1, srcCHA);
        Gadgetron::clear(&unmixCoeff);

        gFactor.create(RO, E1);
        Gadgetron::clear(&gFactor);

        int src;

        T* pKerIm = const_cast<T*>(kerIm.begin());
        T* pCoilMap = const_cast<T*>(coilMap.begin());
        T* pCoeff = unmixCoeff.begin();

        std::vector<size_t> dim(2);
        dim[0] = RO;
        dim[1] = E1;

        #pragma omp parallel default(none) private(src) shared(RO, E1, srcCHA, dstCHA, pKerIm, pCoilMap, pCoeff, dim)
        {
            hoNDArray<T> coeff2D, coeffTmp(&dim);
            hoNDArray<T> coilMap2D;
            hoNDArray<T> kerIm2D;

            #pragma omp for
            for ( src=0; src<(int)srcCHA; src++ )
            {
                coeff2D.create(&dim, pCoeff+src*RO*E1);

                for ( size_t dst=0; dst<dstCHA; dst++ )
                {
                    kerIm2D.create(&dim, pKerIm+src*RO*E1+dst*RO*E1*srcCHA);
                    coilMap2D.create(&dim, pCoilMap+dst*RO*E1);
                    Gadgetron::multiplyConj(kerIm2D, coilMap2D, coeffTmp);
                    Gadgetron::add(coeff2D, coeffTmp, coeff2D);
                }
            }
        }

        hoNDArray<T> conjUnmixCoeff(unmixCoeff);
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiplyConj(unmixCoeff, conjUnmixCoeff, conjUnmixCoeff));

        hoNDArray<T> gFactorBuf(RO, E1, 1, gFactor.begin());
        GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(conjUnmixCoeff, gFactorBuf, 2));

        Gadgetron::sqrt(gFactor, gFactor);
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::unmixCoeff(const hoNDArray<T>& kerIm, const hoNDArray<T>& coilMap, hoNDArray<T>& unmixCoeff, hoNDArray<T>& gFactor) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::applyImageDomainKernel(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t srcCHA = kerIm.get_size(2);
        size_t dstCHA = kerIm.get_size(3);

        GADGET_CHECK_RETURN_FALSE(kspace.get_size(0)==RO);
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(1)==E1);
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(2)==srcCHA);

        buffer2DT_unwrapping_ = kspace;

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, buffer2DT_unwrapping_);

        GADGET_CHECK_RETURN_FALSE(applyImageDomainKernelImage(buffer2DT_unwrapping_, kerIm, complexIm));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::applyImageDomainKernel(const hoNDArray<T>& kspace, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& complexIm)
{
    hoNDArray<T> buf4D(kerIm.get_dimensions());
    return applyImageDomainKernelImage(aliasedIm, kerIm, buf4D, complexIm);
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& kerImBuffer, hoNDArray<T>& complexIm)
{
    try
    {
        size_t RO = kerIm.get_size(0);
        size_t E1 = kerIm.get_size(1);
        size_t srcCHA = kerIm.get_size(2);
        size_t dstCHA = kerIm.get_size(3);

        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(0)==RO);
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(1)==E1);
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(2)==srcCHA);

        boost::shared_ptr< std::vector<size_t> > dim = aliasedIm.get_dimensions();

        std::vector<size_t> dimIm(*dim);
        dimIm[2] = dstCHA;

        if ( !complexIm.dimensions_equal(&dimIm) )
        {
            complexIm.create(&dimIm);
        }
        Gadgetron::clear(&complexIm);

        std::vector<size_t> dim3D(3);
        dim3D[0] = RO;
        dim3D[1] = E1;
        dim3D[2] = srcCHA;

        std::vector<size_t> dimIm3D(3);
        dimIm3D[0] = RO;
        dimIm3D[1] = E1;
        dimIm3D[2] = dstCHA;

        size_t num = aliasedIm.get_number_of_elements()/ (RO*E1*srcCHA);

        int n;

        if ( num <= 8 )
        {
            if ( performTiming_ ) { gt_timer3_.start("apply image domain kernel image ... "); }
            for ( n=0; n<(int)num; n++ )
            {
                hoNDArray<T> buf3D(&dim3D, const_cast<T*>(aliasedIm.begin()+n*RO*E1*srcCHA));
                hoNDArray<T> bufIm3D(RO, E1, 1, dstCHA, complexIm.begin() + n*RO*E1*dstCHA);

                Gadgetron::multiply(kerIm, buf3D, kerImBuffer);
                Gadgetron::sum_over_dimension(kerImBuffer, bufIm3D, 2);
            }
            if ( performTiming_ ) { gt_timer3_.stop(); }
        }
        else
        {
            #pragma omp parallel default(none) private(n) shared(kerIm, num, dim3D, aliasedIm, RO, E1, srcCHA, dimIm3D, dstCHA, complexIm) 
            {
                hoNDArray<T> buf3D;
                hoNDArray<T> bufIm3D;
                hoNDArray<T> buf4D(kerIm.get_dimensions());

                #pragma omp for
                for ( n=0; n<(int)num; n++ )
                {
                    buf3D.create(&dim3D, const_cast<T*>(aliasedIm.begin()+n*RO*E1*srcCHA));
                    bufIm3D.create(RO, E1, 1, dstCHA, complexIm.begin() + n*RO*E1*dstCHA);

                    Gadgetron::multiply(kerIm, buf3D, buf4D);
                    Gadgetron::sum_over_dimension(buf4D, bufIm3D, 2);
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::applyImageDomainKernelImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& kerIm, hoNDArray<T>& kerImBuffer, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_RETURN_FALSE(kspace.get_size(2)==unmixCoeff.get_size(2));

        buffer2DT_unwrapping_ = kspace;

        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, buffer2DT_unwrapping_);
        GADGET_CHECK_RETURN_FALSE(applyUnmixCoeffImage(buffer2DT_unwrapping_, unmixCoeff, complexIm));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::applyUnmixCoeff(const hoNDArray<T>& kspace, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(0)==unmixCoeff.get_size(0));
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(1)==unmixCoeff.get_size(1));
        GADGET_CHECK_RETURN_FALSE(aliasedIm.get_size(2)==unmixCoeff.get_size(2));

        boost::shared_ptr< std::vector<size_t> > dim = aliasedIm.get_dimensions();

        std::vector<size_t> dimIm(*dim);
        dimIm[2] = 1;

        if ( !complexIm.dimensions_equal(&dimIm) )
        {
            complexIm.create(&dimIm);
        }
        Gadgetron::clear(&complexIm);

        buffer2DT_unwrapping_ = aliasedIm;

        Gadgetron::multiply(aliasedIm, unmixCoeff, buffer2DT_unwrapping_);
        GADGET_CATCH_THROW(Gadgetron::sum_over_dimension(buffer2DT_unwrapping_, complexIm, 2));
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::applyUnmixCoeffImage(const hoNDArray<T>& aliasedIm, const hoNDArray<T>& unmixCoeff, hoNDArray<T>& complexIm) ... ");
        return false;
    }

    return true;
}

template <typename T>
bool gtPlusReconWorker2DT<T>::appy_KLT_coil_compression_coeff_2D(const hoNDArray<T>& data, const std::vector< hoNDKLT<T> >& coeff, hoNDArray<T>& dataEigen)
{
    try
    {
        size_t NDim = data.get_number_of_dimensions();
        GADGET_CHECK_RETURN_FALSE(NDim >= 3);

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
        dimEigen[2] = dstCHA;

        dataEigen.create(&dimEigen);
        long long eigenN = dataEigen.get_number_of_elements() / LastDimData;

        std::vector<size_t> dimLastDim(NDim - 1);
        for (n = 0; n<NDim - 1; n++)
        {
            dimLastDim[n] = dim[n];
        }

        std::vector<size_t> dimEigenLastDim(dimLastDim);
        dimEigenLastDim[2] = dstCHA;

        if (LastDimData>1)
        {
            hoNDArray<T> dataEigenLastDim;
            for (n = 0; n < LastDimData; n++)
            {
                hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin() + n*N));
                coeff[n].transform(dataLastDim, dataEigenLastDim, 2);
                memcpy(dataEigen.begin() + n*eigenN, dataEigenLastDim.begin(), dataEigenLastDim.get_number_of_bytes());
            }
        }
        else
        {
            hoNDArray<T> dataLastDim(&dimLastDim, const_cast<T*>(data.begin()));
            hoNDArray<T> dataEigenLastDim(&dimEigenLastDim, dataEigen.begin());
            coeff[0].transform(dataLastDim, dataEigenLastDim, 2);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::appy_KLT_coil_compression_coeff_2D(std::vector<hoNDArray<T> >& coeff) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::afterUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT)
{
    try
    {
        bool fullres_coilmap = false;
        bool fullres_coilmap_useHighestSignal = false;
        bool ref_fillback = false;
        bool averageallN_coilmap = false;
        int numOfModesKept = 0;
        bool same_coilmap_allS = false;
        size_t whichS_coilmap = 0;

        size_t RO = workOrder2DT->kernelIm_->get_size(0);
        size_t E1 = workOrder2DT->kernelIm_->get_size(1);
        size_t srcCHA = workOrder2DT->kernelIm_->get_size(2);
        size_t dstCHA = workOrder2DT->kernelIm_->get_size(3);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        if ( workOrder2DT->CalibMode_ == ISMRMRD_noacceleration )
        {
            fullres_coilmap = false;
            ref_fillback = false;
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_embedded )
        {
            if ( workOrder2DT->embedded_fullres_coilmap_ )
            {
                fullres_coilmap = true;
                fullres_coilmap_useHighestSignal = workOrder2DT->embedded_fullres_coilmap_useHighestSignal_;
            }

            if ( workOrder2DT->embedded_ref_fillback_ 
                && (workOrder2DT->recon_algorithm_!=ISMRMRD_SPIRIT) 
                && (workOrder2DT->recon_algorithm_!=ISMRMRD_L1SPIRIT)
                && (workOrder2DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP)
                && (workOrder2DT->recon_algorithm_!=ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP) )
            {
                ref_fillback = true;
            }

            if ( workOrder2DT->embedded_averageall_ref_ )
            {
                averageallN_coilmap = true;
            }

            if ( workOrder2DT->embedded_same_combinationcoeff_allS_ )
            {
                same_coilmap_allS = true;
                whichS_coilmap = workOrder2DT->embedded_whichS_combinationcoeff_;
            }

            numOfModesKept = workOrder2DT->embedded_ref_numOfModes_;
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_separate )
        {
            if ( workOrder2DT->separate_fullres_coilmap_ )
            {
                fullres_coilmap = true;
            }

            if ( workOrder2DT->separate_averageall_ref_ )
            {
                averageallN_coilmap = true;
            }

            if ( workOrder2DT->separate_same_combinationcoeff_allS_ )
            {
                same_coilmap_allS = true;
                whichS_coilmap = workOrder2DT->separate_whichS_combinationcoeff_;
            }

            numOfModesKept = workOrder2DT->separate_ref_numOfModes_;
        }

        if ( whichS_coilmap >= S ) whichS_coilmap = S-1;

        if ( ref_fillback )
        {
            GDEBUG_STREAM("Fill back the reference kspace lines to the reconstruction ");

            hoNDArray<T> ref_dst;
            if ( workOrder2DT->coil_compression_ )
            {
                this->appy_KLT_coil_compression_coeff_2D(workOrder2DT->ref_, *workOrder2DT->coilCompressionCoef_, ref_dst);
            }
            else
            {
                ref_dst = workOrder2DT->ref_;
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(ref_dst, debugFolder_+"ref_dst"); }

            if ( (ref_dst.get_size(2)==dstCHA) && (ref_dst.get_size(3)==N) && (ref_dst.get_size(4)==S) )
            {
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->fullkspace_, debugFolder_+"fullkspace_"); }

                GADGET_CHECK_RETURN_FALSE(gtPlus_util_.copyAlongE1(ref_dst, workOrder2DT->fullkspace_, startE1_, endE1_));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->fullkspace_, debugFolder_+"fullkspace_After"); }
            }
        }

        // partial fourier handling
        if ( partial_fourier_handling_ )
        {
            GADGET_CHECK_RETURN_FALSE(this->performPartialFourierHandling(workOrder2DT));
        }

        if ( fullres_coilmap )
        {
            if ( performTiming_ ) { gt_timer2_.start("full res coil map : allocate buffer 2DT ...  "); }
            hoNDArray<T> buffer2DT_Two(workOrder2DT->fullkspace_.get_dimensions());
            if ( performTiming_ ) { gt_timer2_.stop(); }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(workOrder2DT->fullkspace_, buffer2DT_, buffer2DT_Two);
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_, debugFolder_+"ComplexIm_afterRefFill"); }

            if ( averageallN_coilmap )
            {
                if ( workOrder2DT->workFlow_use_BufferedKernel_ && workOrder2DT->coilMap_->get_size(3)==1 && workOrder2DT->coilMap_->get_size(4)==S )
                {
                    size_t s;
                    for ( s=0; s<S; s++ )
                    {
                        hoNDArray<T> coilMapS(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+s*RO*E1*dstCHA);
                        hoNDArray<T> complexImS(RO, E1, dstCHA, N, buffer2DT_.begin()+s*RO*E1*dstCHA*N);
                        hoNDArray<T> complexImCombinedS(RO, E1, N, workOrder2DT->complexIm_.begin()+s*RO*E1*N);

                        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImS, coilMapS, complexImCombinedS));
                        Gadgetron::coil_combine(complexImS, coilMapS, 2, complexImCombinedS);
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(complexImCombinedS, debugFolder_+"complexImCombinedS"); }
                    }
                }
                else
                {
                    workOrder2DT->coilMap_->create(RO, E1, dstCHA, 1, S);
                    //Gadgetron::clear( *(workOrder2DT->coilMap_) );

                    size_t s;

                    if ( same_coilmap_allS )
                    {
                        hoNDArray<T> aveComplexImS(RO, E1, dstCHA, 1);
                        //Gadgetron::clear(aveComplexImS);

                        buffer2DT_unwrapping_.create(RO, E1, dstCHA, N);
                        //Gadgetron::clear(aveComplexImS);

                        hoMatrix<T> A(RO*E1*dstCHA, N, buffer2DT_.begin()+whichS_coilmap*RO*E1*dstCHA*N);
                        hoMatrix<T> A_KLF(RO*E1*dstCHA, N, buffer2DT_unwrapping_.begin());

                        if ( numOfModesKept>0 && numOfModesKept<dstCHA )
                        {
                            hoNDKLT<T> klt;
                            klt.prepare(A, (size_t)1, (size_t)0);
                            klt.KL_filter(A, A_KLF, 1, numOfModesKept);
                        }
                        else
                        {
                            memcpy(A_KLF.begin(), A.begin(), A_KLF.get_number_of_bytes());
                        }

                        if ( fullres_coilmap_useHighestSignal )
                        {
                            GADGET_CHECK_RETURN_FALSE(pickHighestSignalForN(buffer2DT_unwrapping_, aveComplexImS));
                        }
                        else
                        {
                            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(buffer2DT_unwrapping_, aveComplexImS));
                        }

                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(aveComplexImS, debugFolder_+"aveComplexImS"); }

                        hoNDArray<T> coilMapS(RO, E1, dstCHA, 1, workOrder2DT->coilMap_->begin()+whichS_coilmap*RO*E1*dstCHA);

                        if ( performTiming_ ) { gt_timer2_.start("coilMap2DNIH ...  "); }
                        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(aveComplexImS, coilMapS, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_));
                        if ( performTiming_ ) { gt_timer2_.stop(); }

                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(coilMapS, debugFolder_+"coilMapS"); }

                        GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder2DT->coilMap_, whichS_coilmap));

                        for ( s=0; s<S; s++ )
                        {
                            hoNDArray<T> coilMapS(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+s*RO*E1*dstCHA);
                            hoNDArray<T> complexImS(RO, E1, dstCHA, N, buffer2DT_.begin()+s*RO*E1*dstCHA*N);
                            hoNDArray<T> complexImCombinedS(RO, E1, N, workOrder2DT->complexIm_.begin()+s*RO*E1*N);

                            // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImS, coilMapS, complexImCombinedS));
                            Gadgetron::coil_combine(complexImS, coilMapS, 2, complexImCombinedS);
                            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(complexImCombinedS, debugFolder_+"complexImCombinedS"); }
                        }
                    }
                    else
                    {
                        hoNDArray<T> aveComplexIm(RO, E1, dstCHA, 1, S);
                        //Gadgetron::clear(aveComplexIm);

                        buffer2DT_unwrapping_ = buffer2DT_;
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_unwrapping_, debugFolder_+"buffer2DT_unwrapping"); }

                        if ( numOfModesKept>0 && numOfModesKept<dstCHA )
                        {
                            for ( s=0; s<S; s++ )
                            {
                                hoMatrix<T> A(RO*E1*dstCHA, N, buffer2DT_.begin()+s*RO*E1*dstCHA*N);
                                hoMatrix<T> A_KLF(RO*E1*dstCHA, N, buffer2DT_unwrapping_.begin()+s*RO*E1*dstCHA*N);

                                hoNDKLT<T> klt;
                                klt.prepare(A, (size_t)1, (size_t)0);
                                klt.KL_filter(A, A_KLF, 1, numOfModesKept);
                            }

                            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_unwrapping_, debugFolder_+"ComplexIm_KLF"); }
                        }

                        if ( fullres_coilmap_useHighestSignal )
                        {
                            GADGET_CHECK_RETURN_FALSE(pickHighestSignalForN(buffer2DT_unwrapping_, aveComplexIm));
                        }
                        else
                        {
                            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(buffer2DT_unwrapping_, aveComplexIm));
                        }

                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(aveComplexIm, debugFolder_+"aveComplexIm"); }

                        if ( performTiming_ ) { gt_timer2_.start("coilMap2DNIH ...  "); }

                        gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(aveComplexIm, *workOrder2DT->coilMap_, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_);

                        // gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer2DT_, *workOrder2DT->coilMap_, workOrder2DT->complexIm_);
                        Gadgetron::coil_combine(buffer2DT_, *workOrder2DT->coilMap_, 2, workOrder2DT->complexIm_);

                        //long long ss;
                        //#pragma omp parallel for private(s) if (S>2)
                        //for ( ss=0; ss<S; ss++ )
                        //{
                        //    hoNDArray<T> aveComplexImS(RO, E1, dstCHA, aveComplexIm.begin()+ss*RO*E1*dstCHA);
                        //    hoNDArray<T> coilMapS(RO, E1, dstCHA, workOrder2DT->coilMap_->begin()+ss*RO*E1*dstCHA);

                        //    //GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(aveComplexImS, coilMapS, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, workOrder2DT->csm_iter_thres_));
                        //    gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(aveComplexImS, coilMapS, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, workOrder2DT->csm_iter_thres_);

                        //    hoNDArray<T> complexImS(RO, E1, dstCHA, N, buffer2DT_.begin()+ss*RO*E1*dstCHA*N);
                        //    hoNDArray<T> complexImCombinedS(RO, E1, N, workOrder2DT->complexIm_.begin()+ss*RO*E1*N);

                        //    //GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImS, coilMapS, complexImCombinedS));
                        //    gtPlusISMRMRDReconUtilComplex<T>().coilCombine(complexImS, coilMapS, complexImCombinedS);
                        //}
                        if ( performTiming_ ) { gt_timer2_.stop(); }

                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder2DT->coilMap_, debugFolder_+"coilMap_fullres"); }
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_+"complexImCombined"); }
                    }
                }
            }
            else
            {
                if ( workOrder2DT->workFlow_use_BufferedKernel_ && workOrder2DT->coilMap_->get_size(3)==N && workOrder2DT->coilMap_->get_size(4)==S )
                {
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer2DT_, *workOrder2DT->coilMap_, workOrder2DT->complexIm_));
                    Gadgetron::coil_combine(buffer2DT_, *workOrder2DT->coilMap_, 2, workOrder2DT->complexIm_);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_+"complexIm_"); }
                }
                else
                {
                    workOrder2DT->coilMap_->create(RO, E1, dstCHA, N, S);

                    if ( performTiming_ ) { gt_timer2_.start("coilMap2DNIH ...  "); }
                    if ( same_coilmap_allS )
                    {
                        hoNDArray<T> complexImS(RO, E1, dstCHA, N, buffer2DT_.begin()+whichS_coilmap*RO*E1*dstCHA*N);
                        hoNDArray<T> coilMapS(RO, E1, dstCHA, N, workOrder2DT->coilMap_->begin()+whichS_coilmap*RO*E1*dstCHA*N);

                        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(complexImS, coilMapS, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_));
                        GADGET_CHECK_RETURN_FALSE(repmatLastDimension(*workOrder2DT->coilMap_, whichS_coilmap));
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder2DT->coilMap_, debugFolder_+"coilMap_fullres"); }
                    }
                    else
                    {
                        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilMap2DNIH(buffer2DT_, *workOrder2DT->coilMap_, workOrder2DT->coil_map_algorithm_, workOrder2DT->csm_kSize_, workOrder2DT->csm_powermethod_num_, workOrder2DT->csm_iter_num_, (value_type)workOrder2DT->csm_iter_thres_, workOrder2DT->csm_use_gpu_));
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*workOrder2DT->coilMap_, debugFolder_+"coilMap_fullres"); }
                    }
                    if ( performTiming_ ) { gt_timer2_.stop(); }

                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer2DT_, *workOrder2DT->coilMap_, workOrder2DT->complexIm_));
                    Gadgetron::coil_combine(buffer2DT_, *workOrder2DT->coilMap_, 2, workOrder2DT->complexIm_);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_+"complexIm_"); }
                }
            }
        }
        else
        {
            if ( partial_fourier_handling_ )
            {
                bool partialFourierHandling = true;
                if ( (workOrder2DT->start_RO_<0 || workOrder2DT->end_RO_<0 || (workOrder2DT->end_RO_-workOrder2DT->start_RO_+1==RO) ) 
                        && (workOrder2DT->start_E1_<0 || workOrder2DT->end_E1_<0 || (workOrder2DT->end_E1_-workOrder2DT->start_E1_+1==E1) ) )
                {
                    partialFourierHandling = false;
                }

                // if the partial fourier handling is used to compute updated full kspace, the coil combination needs to be repeated
                if ( partialFourierHandling )
                {
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_+"complexIm_origin_noFullResCoilMap_"); }

                    // if the partial fourier handling is performed on the fullkspace, an extra coil combination is needed
                    if ( workOrder2DT->CalibMode_ == ISMRMRD_noacceleration )
                    {
                        hoNDArray<T> buffer2DT_Two(workOrder2DT->data_.get_dimensions());
                        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(workOrder2DT->data_, buffer2DT_, buffer2DT_Two);
                        // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer2DT_, *workOrder2DT->coilMap_, workOrder2DT->complexIm_));
                        Gadgetron::coil_combine(buffer2DT_, *workOrder2DT->coilMap_, 2, workOrder2DT->complexIm_);
                        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_+"complexIm_noFullResCoilMap_"); }
                    }
                    else if ( workOrder2DT->fullkspace_.get_number_of_elements() > 0 )
                    {
                        if ( workOrder2DT->fullkspace_.get_size(2) == workOrder2DT->coilMap_->get_size(2) )
                        {
                            hoNDArray<T> buffer2DT_Two(workOrder2DT->fullkspace_.get_dimensions());
                            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(workOrder2DT->fullkspace_, buffer2DT_, buffer2DT_Two);
                            // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().coilCombine(buffer2DT_, *workOrder2DT->coilMap_, workOrder2DT->complexIm_));
                            Gadgetron::coil_combine(buffer2DT_, *workOrder2DT->coilMap_, 2, workOrder2DT->complexIm_);
                            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_+"complexIm_noFullResCoilMap_"); }
                        }
                        else if (workOrder2DT->fullkspace_.get_size(2) == 1) // if recon kspace is not required
                        {
                            hoNDArray<T> buffer2DT_ComplexIm(workOrder2DT->fullkspace_.get_dimensions());
                            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(workOrder2DT->fullkspace_, buffer2DT_ComplexIm);
                            memcpy(workOrder2DT->complexIm_.begin(), buffer2DT_ComplexIm.begin(), workOrder2DT->complexIm_.get_number_of_bytes());
                            if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(workOrder2DT->complexIm_, debugFolder_ + "complexIm_noFullResCoilMap_noReconKSpace_"); }
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::afterUnwrapping(gtPlusReconWorkOrder2DT<T>* workOrder2DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::pickHighestSignalForN(const hoNDArray<T>& data, hoNDArray<T>& res)
{
    try
    {
        size_t RO = data.get_size(0);
        size_t E1 = data.get_size(1);
        size_t CHA = data.get_size(2);
        size_t N = data.get_size(3);
        size_t S = data.get_size(4);

        res.create(RO, E1, CHA, 1, S);

        size_t s;
        for ( s=0; s<S; s++ )
        {
            size_t maxInd=0;
            typename realType<T>::Type maxNorm;

            hoNDArray<T> data3D(RO, E1, CHA, const_cast<T*>(data.begin()+s*RO*E1*CHA*N));
            Gadgetron::norm2(data3D, maxNorm);

            size_t n;
            for ( n=1; n<N; n++ )
            {
                data3D.create(RO, E1, CHA, const_cast<T*>(data.begin()+n*RO*E1*CHA+s*RO*E1*CHA*N));

                typename realType<T>::Type currNorm;
                Gadgetron::norm2(data3D, currNorm);

                if ( maxNorm < currNorm )
                {
                    maxNorm = currNorm;
                    maxInd = n;
                }
            }

            memcpy(res.begin()+s*RO*E1*CHA*N, data.begin()+maxInd*RO*E1*CHA+s*RO*E1*CHA*N, sizeof(T)*RO*E1*CHA);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::pickHighestSignalForN() ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performPartialFourierHandling(gtPlusReconWorkOrder2DT<T>* workOrder2DT)
{
    try
    {
        // compensate for the partial fourier to preserve the SNR unit
        // value_type partialFourierCompensationFactor = 1;

        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);

        //if ( !( workOrder2DT->start_RO_<0 || workOrder2DT->end_RO_<0 || (workOrder2DT->end_RO_-workOrder2DT->start_RO_+1==RO) ) )
        //{
        //    partialFourierCompensationFactor *= (value_type)(RO)/(value_type)(workOrder2DT->end_RO_-workOrder2DT->start_RO_+1);
        //}

        /*if ( !( workOrder2DT->start_E1_<0 || workOrder2DT->end_E1_<0 || (workOrder2DT->end_E1_-workOrder2DT->start_E1_+1==E1) ) )
        {
            if ( workOrder2DT->end_E1_-workOrder2DT->start_E1_+1 <= E1 )
            {
                partialFourierCompensationFactor *= (value_type)(E1)/(value_type)(workOrder2DT->end_E1_-workOrder2DT->start_E1_+1);
            }
        }*/

        /*partialFourierCompensationFactor = std::sqrt(partialFourierCompensationFactor);
        if ( performTiming_ ) { GDEBUG_STREAM("Partial fourier scaling factor : " << partialFourierCompensationFactor); }*/

        if ( performTiming_ )
        {
            GDEBUG_STREAM("Partial fourier algorithm : " << gtPlus_util_.getNameFromISMRMRDPartialFourierReconAlgo(workOrder2DT->partialFourier_algo_));
            GDEBUG_STREAM("RO : [" << workOrder2DT->start_RO_ << " " << workOrder2DT->end_RO_ << "]; E1 : [" << workOrder2DT->start_E1_ << " " << workOrder2DT->end_E1_ << "]");
        }

        if ( workOrder2DT->CalibMode_ == ISMRMRD_noacceleration )
        {
            /*if ( (workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING || workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER) && (std::abs(partialFourierCompensationFactor-1)>FLT_EPSILON) )
            {
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal(partialFourierCompensationFactor, workOrder2DT->data_));
            }*/

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFilter(*workOrder2DT, workOrder2DT->data_));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_HOMODYNE )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierHomodyneRecon(*workOrder2DT, workOrder2DT->data_));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_POCS )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierPOCSRecon(*workOrder2DT, workOrder2DT->data_));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_FENGHUANG )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFengHuangRecon(*workOrder2DT, workOrder2DT->data_));
            }
        }
        else if ( workOrder2DT->fullkspace_.get_number_of_elements() > 0 )
        {
            /*if ( (workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING || workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER) && (std::abs(partialFourierCompensationFactor-1)>FLT_EPSILON) )
            {
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal(partialFourierCompensationFactor, workOrder2DT->fullkspace_));
            }*/

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFilter(*workOrder2DT, workOrder2DT->fullkspace_));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_HOMODYNE )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierHomodyneRecon(*workOrder2DT, workOrder2DT->fullkspace_));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_POCS )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierPOCSRecon(*workOrder2DT, workOrder2DT->fullkspace_));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_FENGHUANG )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFengHuangRecon(*workOrder2DT, workOrder2DT->fullkspace_));
            }
        }
        else
        {
            // perform partial fourier handling on the complex images after coil combination
            hoNDArray<T> kspace(workOrder2DT->complexIm_);
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(workOrder2DT->complexIm_, kspace);

            /*if ( (workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING || workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER) && (std::abs(partialFourierCompensationFactor-1)>FLT_EPSILON) )
            {
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal(partialFourierCompensationFactor, kspace));
            }*/

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_ZEROFILLING_FILTER )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFilter(*workOrder2DT, kspace));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_HOMODYNE )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierHomodyneRecon(*workOrder2DT, kspace));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_POCS )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierPOCSRecon(*workOrder2DT, kspace));
            }

            if ( workOrder2DT->partialFourier_algo_ == ISMRMRD_PF_FENGHUANG )
            {
                GADGET_CHECK_RETURN_FALSE(performPartialFourierFengHuangRecon(*workOrder2DT, kspace));
            }

            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, workOrder2DT->complexIm_);
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performPartialFourierHandling(gtPlusReconWorkOrder2DT<T>* workOrder2DT) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performPartialFourierFilter(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        // check whether partial fourier is used
        if ( (workOrder2DT.start_RO_<0 || workOrder2DT.end_RO_<0 || (workOrder2DT.end_RO_-workOrder2DT.start_RO_+1==RO) ) 
            && (workOrder2DT.start_E1_<0 || workOrder2DT.end_E1_<0 || (workOrder2DT.end_E1_-workOrder2DT.start_E1_+1==E1) ) )
        {
            return true;
        }

        hoNDArray<T> input;
        input.create(RO, E1, 1, CHA, N, S, kspace.begin());

        if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(kspace, debugFolder_ + "kspace_before_PF_Filter"); }

        double filter_pf_width_RO_ = 0.15;
        double filter_pf_width_E1_ = 0.15;
        double filter_pf_width_E2_ = 0.15;

        bool filter_pf_density_comp_ = false;

        hoNDArray<T> res;

        Gadgetron::partial_fourier_filter(input,
            workOrder2DT.start_RO_, workOrder2DT.end_RO_, workOrder2DT.start_E1_, workOrder2DT.end_E1_, 0, 0,
            filter_pf_width_RO_, filter_pf_width_E1_, filter_pf_width_E2_, filter_pf_density_comp_,
            workOrder2DT.filterRO_partialfourier_, workOrder2DT.filterE1_partialfourier_, workOrder2DT.filterE2_partialfourier_, res);

        memcpy(kspace.begin(), res.begin(), sizeof(T)*RO*E1*CHA*N*S);

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_PF_Filter"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performPartialFourierFilter(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performPartialFourierHomodyneRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        // check whether partial fourier is used
        if ( (workOrder2DT.start_RO_<0 || workOrder2DT.end_RO_<0 || (workOrder2DT.end_RO_-workOrder2DT.start_RO_+1==RO) ) 
            && (workOrder2DT.start_E1_<0 || workOrder2DT.end_E1_<0 || (workOrder2DT.end_E1_-workOrder2DT.start_E1_+1==E1) ) )
        {
            return true;
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_before_homodyne"); }

        // create kspace filter for homodyne phase estimation
        ISMRMRDKSPACEFILTER filter_ref_type_ = ISMRMRD_FILTER_HANNING;
        double filter_ref_sigma_ = 1.5;
        double filter_ref_width_ = 0.15;

        size_t startRO(0), endRO(RO-1);
        hoNDArray<T> filterRO(RO);
        if ( (workOrder2DT.start_RO_<0 || workOrder2DT.end_RO_<0) )
        {
            Gadgetron::generate_symmetric_filter_ref(RO, 0, RO - 1, filterRO);
        }
        else
        {
            Gadgetron::generate_symmetric_filter_ref(RO, workOrder2DT.start_RO_, workOrder2DT.end_RO_, filterRO);

            startRO = workOrder2DT.start_RO_;
            endRO = workOrder2DT.end_RO_;
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(filterRO, "filterRO_homodyne"); }

        size_t startE1(0), endE1(E1-1);
        hoNDArray<T> filterE1(E1);
        if ( (workOrder2DT.start_E1_<0 || workOrder2DT.end_E1_<0) )
        {
            Gadgetron::generate_symmetric_filter_ref(E1, 0, E1 - 1, filterE1);
        }
        else
        {
            Gadgetron::generate_symmetric_filter_ref(E1, workOrder2DT.start_E1_, workOrder2DT.end_E1_, filterE1);

            startE1 = workOrder2DT.start_E1_;
            endE1 = workOrder2DT.end_E1_;
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(filterE1, debugFolder_+"filterE1_homodyne"); }

        hoNDArray<T> kspaceIter(kspace.get_dimensions());
        kspaceIter = kspace;
        // store the filtered kspace
        buffer2DT_partial_fourier_ = kspace;
        // store the phase images
        buffer2DT_ = kspace;
        // magnitude of complex images
        hoNDArray<typename realType<T>::Type> mag(kspace.get_dimensions());
        hoNDArray<T> magComplex(kspace.get_dimensions());

        // complex images
        hoNDArray<T> complexIm(kspace.get_dimensions());
        Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(kspace, complexIm);

        hoNDArray<T> complexImPrev(complexIm);

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"homodyne_kspace_beforeIteration"); }

        size_t ii;
        for ( ii=0; ii<workOrder2DT.partialFourier_homodyne_iters_; ii++ )
        {
            // kspace filter before phase extraction
            Gadgetron::apply_kspace_filter_ROE1(kspaceIter, filterRO, filterE1, buffer2DT_partial_fourier_);

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_partial_fourier_, debugFolder_+"homodyne_kspaceIter_afterFiltered"); }

            // go to image domain
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(buffer2DT_partial_fourier_);
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_partial_fourier_, debugFolder_+"homodyne_complexIm"); }

            // get the phase
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::abs(buffer2DT_partial_fourier_, mag));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::addEpsilon(mag));
            GADGET_CHECK_RETURN_FALSE(magComplex.copyFrom(mag));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::divide(buffer2DT_partial_fourier_, magComplex, buffer2DT_));
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_, debugFolder_+"homodyne_phase"); }

            // remove the phase from complex images
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::conjugate(buffer2DT_, buffer2DT_));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::multiply(complexIm, buffer2DT_, complexIm));
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(complexIm, debugFolder_+"homodyne_complexIm_removePhase"); }

            // go back to kspace
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(complexIm, kspaceIter);
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIter, debugFolder_+"homodyne_complexIm_removePhase_kspace"); }

            // compute threshold to stop the iteration
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::subtract(complexImPrev, complexIm, buffer2DT_));
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_, debugFolder_+"homodyne_diff_complexIm"); }

            typename realType<T>::Type diff, prev;
            Gadgetron::norm2(complexImPrev, prev);
            Gadgetron::norm2(buffer2DT_, diff);

            typename realType<T>::Type thres = diff/prev;

            if ( !debugFolder_.empty() )
            {
                GDEBUG_STREAM("Homodyne iter : " << ii << " - thres : " << thres << " ... ");
            }

            if ( thres < workOrder2DT.partialFourier_homodyne_thres_ )
            {
                break;
            }

            complexImPrev = complexIm;
        }

        // restore the acquired region
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIter, debugFolder_+"kspaceIter_after_homodyne_beforeCopy"); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_homodyne_beforeCopy"); }

        if ( workOrder2DT.partialFourier_homodyne_densityComp_ )
        {
            size_t width_RO = (size_t)std::floor(0.1*RO);
            size_t width_E1 = (size_t)std::floor(0.1*E1);

            // compute PF filter for RO and E1
            hoNDArray<T> filterPF_RO, filterPF_E1;

            if ( workOrder2DT.start_RO_<0 || workOrder2DT.end_RO_<0 || (workOrder2DT.start_RO_==0 && workOrder2DT.end_RO_==RO-1) )
            {
                Gadgetron::generate_asymmetric_filter(RO, workOrder2DT.start_RO_, workOrder2DT.end_RO_, filterPF_RO, ISMRMRD_FILTER_NONE, width_RO, true);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(RO, workOrder2DT.start_RO_, workOrder2DT.end_RO_, filterPF_RO, ISMRMRD_FILTER_TAPERED_HANNING, width_RO, true);
            }

            if ( workOrder2DT.start_E1_<0 || workOrder2DT.end_E1_<0 || (workOrder2DT.start_E1_==0 && workOrder2DT.end_E1_==E1-1) )
            {
                Gadgetron::generate_asymmetric_filter(E1, workOrder2DT.start_E1_, workOrder2DT.end_E1_, filterPF_E1, ISMRMRD_FILTER_NONE, width_E1, true);
            }
            else
            {
                Gadgetron::generate_asymmetric_filter(E1, workOrder2DT.start_E1_, workOrder2DT.end_E1_, filterPF_E1, ISMRMRD_FILTER_TAPERED_HANNING, width_E1, true);
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(filterPF_RO, debugFolder_+"filterPF_RO_homodyne"); }
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(filterPF_E1, debugFolder_+"filterPF_E1_homodyne"); }

            // compensate filter for homodyne filtered kspace
            hoNDArray<T> filterPF_homodyne_RO(filterPF_RO), filterPF_homodyne_E1(filterPF_E1);

            T midValue = filterPF_RO(RO/2);
            for ( ii=0; ii<RO; ii++ )
            {
                if ( std::abs(filterPF_homodyne_RO(ii)) > std::abs(midValue) )
                {
                    filterPF_homodyne_RO(ii) = T(0.0);
                }
                else
                {
                    filterPF_homodyne_RO(ii) = midValue - filterPF_homodyne_RO(ii);
                }
            }

            midValue = filterPF_E1(E1/2);
            for ( ii=0; ii<E1; ii++ )
            {
                if ( std::abs(filterPF_homodyne_E1(ii)) > std::abs(midValue) )
                {
                    filterPF_homodyne_E1(ii) = T(0.0);
                }
                else
                {
                    filterPF_homodyne_E1(ii) = midValue - filterPF_homodyne_E1(ii);
                }
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(filterPF_homodyne_RO, "filterPF_homodyne_RO_homodyne"); }
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(filterPF_homodyne_E1, "filterPF_homodyne_E1_homodyne"); }

            T scaleFactor(1.0);
            hoNDArray<T> filterPF;

            if ( workOrder2DT.start_RO_<0 || workOrder2DT.end_RO_<0 || (workOrder2DT.start_RO_==0 && workOrder2DT.end_RO_==RO-1) )
            {
                Gadgetron::apply_kspace_filter_E1(kspace, filterPF_E1, kspace);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_homodyne_PF_Filter"); }

                Gadgetron::apply_kspace_filter_E1(kspaceIter, filterPF_homodyne_E1, kspaceIter);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIter, debugFolder_+"kspaceIter_after_homodyne_PF_Filter"); }

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::add(filterPF_E1, filterPF_homodyne_E1, filterPF));
                Gadgetron::compute_filter_SNR_unit_scale_factor(filterPF, scaleFactor);
            }
            else if ( workOrder2DT.start_E1_<0 || workOrder2DT.end_E1_<0 || (workOrder2DT.start_E1_==0 && workOrder2DT.end_E1_==E1-1) )
            {
                Gadgetron::apply_kspace_filter_RO(kspace, filterPF_RO, kspace);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_homodyne_PF_Filter"); }

                Gadgetron::apply_kspace_filter_RO(kspaceIter, filterPF_homodyne_RO, kspaceIter);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIter, debugFolder_+"kspaceIter_after_homodyne_PF_Filter"); }

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::add(filterPF_RO, filterPF_homodyne_RO, filterPF));
                Gadgetron::compute_filter_SNR_unit_scale_factor(filterPF, scaleFactor);
            }
            else
            {
                Gadgetron::apply_kspace_filter_ROE1(kspace, filterPF_RO, filterPF_E1, kspace);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_homodyne_PF_Filter"); }

                Gadgetron::apply_kspace_filter_ROE1(kspaceIter, filterPF_homodyne_RO, filterPF_homodyne_E1, kspaceIter);
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIter, debugFolder_+"kspaceIter_after_homodyne_PF_Filter"); }

                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::add(filterPF_RO, filterPF_homodyne_RO, filterPF));
                Gadgetron::compute_filter_SNR_unit_scale_factor(filterPF, scaleFactor);

                T scaleFactorE1(1.0);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::add(filterPF_E1, filterPF_homodyne_E1, filterPF));
                Gadgetron::compute_filter_SNR_unit_scale_factor(filterPF, scaleFactorE1);

                scaleFactor *= scaleFactorE1;
            }

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::add(kspace, kspaceIter, kspace));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::scal(scaleFactor, kspace));
        }
        else
        {
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.copyAlongROE1(kspace, kspaceIter, startRO, endRO, startE1, endE1));
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspaceIter, debugFolder_+"kspaceIter_after_homodyne_afterCopy"); }
            kspace = kspaceIter;
        }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_homodyne"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performPartialFourierHomodyneRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performPartialFourierPOCSRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        // check whether partial fourier is used
        if ( (workOrder2DT.start_RO_<0 || workOrder2DT.end_RO_<0 || (workOrder2DT.end_RO_-workOrder2DT.start_RO_+1==RO) ) 
            && (workOrder2DT.start_E1_<0 || workOrder2DT.end_E1_<0 || (workOrder2DT.end_E1_-workOrder2DT.start_E1_+1==E1) ) )
        {
            return true;
        }

        hoNDArray<T> input, res;
        input.create(RO, E1, 1, CHA, N, S, kspace.begin());

        if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(kspace, debugFolder_ + "kspace_before_POCS"); }

        Gadgetron::partial_fourier_POCS(input,
            workOrder2DT.start_RO_, workOrder2DT.end_RO_, workOrder2DT.start_E1_, workOrder2DT.end_E1_, 0, 0,
            workOrder2DT.partialFourier_POCS_transitBand_, workOrder2DT.partialFourier_POCS_transitBand_, 0,
            workOrder2DT.partialFourier_POCS_iters_, workOrder2DT.partialFourier_POCS_thres_, res);

        memcpy(kspace.begin(), res.begin(), sizeof(T)*RO*E1*CHA*N*S);

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_POCS"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performPartialFourierPOCSRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performPartialFourierFengHuangRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace)
{
    try
    {
        size_t RO = kspace.get_size(0);
        size_t E1 = kspace.get_size(1);
        size_t CHA = kspace.get_size(2);
        size_t N = kspace.get_size(3);
        size_t S = kspace.get_size(4);

        // check whether partial fourier is used
        if ( (workOrder2DT.start_RO_<0 || workOrder2DT.end_RO_<0 || (workOrder2DT.end_RO_-workOrder2DT.start_RO_+1==RO) ) 
            && (workOrder2DT.start_E1_<0 || workOrder2DT.end_E1_<0 || (workOrder2DT.end_E1_-workOrder2DT.start_E1_+1==E1) ) )
        {
            return true;
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_before_FengHuang"); }

        size_t startRO(0), endRO(RO-1);
        if ( workOrder2DT.start_RO_>=0 && workOrder2DT.end_RO_<RO )
        {
            startRO = workOrder2DT.start_RO_;
            endRO = workOrder2DT.end_RO_;
        }

        size_t startE1(0), endE1(E1-1);
        if ( workOrder2DT.start_E1_>=0 && workOrder2DT.end_E1_<E1 )
        {
            startE1 = workOrder2DT.start_E1_;
            endE1 = workOrder2DT.end_E1_;
        }

        // compute the conjugate symmetric kspace
        if ( performTiming_ ) { gt_timer1_.start("conjugateSymmetry2D"); }
        GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().conjugateSymmetry2D(kspace, buffer2DT_));
        if ( performTiming_ ) { gt_timer1_.stop(); }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_, debugFolder_+"kspaceConj_FengHuang"); }

        // find the symmetric region in the kspace
        size_t startSymRO, endSymRO;
        Gadgetron::find_symmetric_sampled_region(startRO, endRO, RO/2, startSymRO, endSymRO);

        size_t startSymE1, endSymE1;
        Gadgetron::find_symmetric_sampled_region(startE1, endE1, E1/2, startSymE1, endSymE1);

        // the reference kspace for kernel estimation
        hoNDArray<T> src, dst;
        std::vector<size_t> start(5), size(5);

        start[0] = startSymRO;
        start[1] = startSymE1;
        start[2] = 0;
        start[3] = 0;
        start[4] = 0;

        size[0] = endSymRO-startSymRO+1;
        size[1] = endSymE1-startSymE1+1;
        size[2] = CHA;
        size[3] = N;
        size[4] = S;

        GADGET_CHECK_RETURN_FALSE(Gadgetron::cropUpTo11DArray(buffer2DT_, src, start, size));
        GADGET_CHECK_RETURN_FALSE(cropUpTo11DArray(kspace, dst, start, size));

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(src, debugFolder_+"src_FengHuang"); }
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(dst, debugFolder_+"dst_FengHuang"); }

        if ( workOrder2DT.partialFourier_FengHuang_sameKernel_allN_ )
        {
            hoNDArray<T> ave4D;
            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(src, ave4D));
            src = ave4D;

            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.averageKSpace4D(dst, ave4D));
            dst = ave4D;

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(src, debugFolder_+"src_ave4D_FengHuang"); }
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(dst, debugFolder_+"dst_ave4D_FengHuang"); }
        }

        // estimate the kernels
        ho6DArray<T> kernel; // [RO E1 srcCHA dstCHA N S]
        if ( performTiming_ ) { gt_timer1_.start("calibFengHuang"); }
        GADGET_CHECK_RETURN_FALSE(this->calibFengHuang(workOrder2DT, src, dst, kernel));
        if ( performTiming_ ) { gt_timer1_.stop(); }

        // perform the recon
        if ( workOrder2DT.partialFourier_FengHuang_transitBand_==0 )
        {
            if ( performTiming_ ) { gt_timer1_.start("performReconFangHuang"); }
            GADGET_CHECK_RETURN_FALSE(this->performReconFangHuang(workOrder2DT, buffer2DT_, kspace, (int)startRO, (int)endRO, (int)startE1, (int)endE1, kernel));
            if ( performTiming_ ) { gt_timer1_.stop(); }
        }
        else
        {
            if ( performTiming_ ) { gt_timer1_.start("performReconFangHuang with transition band"); }

            size_t tb =  (int)workOrder2DT.partialFourier_FengHuang_transitBand_;

            size_t sRO(startRO), eRO(endRO), sE1(startE1), eE1(endE1);

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

            buffer2DT_partial_fourier_kspaceIter_ = kspace;
            GADGET_CHECK_RETURN_FALSE(this->performReconFangHuang(workOrder2DT, buffer2DT_, 
                    buffer2DT_partial_fourier_kspaceIter_, (int)startRO, (int)endRO, (int)startE1, (int)endE1, kernel));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(buffer2DT_partial_fourier_kspaceIter_, debugFolder_+"kspace_FengHuang_recon"); }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_FengHuang_original"); }

            GADGET_CHECK_RETURN_FALSE(gtPlus_util_.copyAlongROE1TransitionBand(kspace, buffer2DT_partial_fourier_kspaceIter_, 
                    sRO, eRO, sE1, eE1, workOrder2DT.partialFourier_FengHuang_transitBand_, workOrder2DT.partialFourier_FengHuang_transitBand_));

            kspace = buffer2DT_partial_fourier_kspaceIter_;

            if ( performTiming_ ) { gt_timer1_.stop(); }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(kspace, debugFolder_+"kspace_after_FengHuang"); }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performPartialFourierFengHuangRecon(gtPlusReconWorkOrder2DT<T>& workOrder2DT, hoNDArray<T>& kspace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::calibFengHuang(gtPlusReconWorkOrder2DT<T>& workOrder2DT, const hoNDArray<T>& src, const hoNDArray<T>& dst, ho6DArray<T>& kernel)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(src.dimensions_equal(&dst));

        long long RO = (long long)src.get_size(0);
        long long E1 = (long long)src.get_size(1);
        long long srcCHA = (long long)src.get_size(2);
        long long N = (long long)src.get_size(3);
        long long S = (long long)src.get_size(4);

        long long kx = (long long)workOrder2DT.partialFourier_FengHuang_kSize_RO_;
        long long ky = (long long)workOrder2DT.partialFourier_FengHuang_kSize_E1_;

        if ( kx%2 == 0 ) kx++;
        if ( ky%2 == 0 ) ky++;

        long long halfKx = (long long)kx/2;
        long long halfKy = (long long)ky/2;

        // the cross-channel kernel is not estimated
        kernel.createArray(kx, ky, srcCHA, 1, N, S);

        long long ii=0;
        long long num = N*S*srcCHA;

        size_t startRO = halfKx;
        size_t endRO = RO - halfKx - 1;

        size_t startE1 = halfKy;
        size_t endE1 = E1 - halfKy - 1;

        long long rowA, colA, rowB, colB;
        rowA = (endE1-startE1+1)*(endRO-startRO+1); 
        colA = kx*ky;

        rowB = rowA;
        colB = 1;

        double thresReg = workOrder2DT.partialFourier_FengHuang_thresReg_;

        #pragma omp parallel default(none) private(ii) shared(num, RO, E1, srcCHA, N, S, kx, ky, src, dst, kernel, rowA, colA, rowB, colB, startRO, endRO, startE1, endE1, halfKx, halfKy, thresReg)
        {
            hoMatrix<T> A(rowA, colA);
            T* pA = A.begin();

            hoMatrix<T> B(rowB, colB);
            T* pB = B.begin();

            hoMatrix<T> K(colA, colB);

            #pragma omp for
            for ( ii=0; ii<num; ii ++ )
            {
                T* pSrc2D = const_cast<T*>(src.begin())+ii*RO*E1;
                T* pDst2D = const_cast<T*>(dst.begin())+ii*RO*E1;
                //ho2DArray<T> src2D(RO, E1, const_cast<T*>(src.begin())+ii*RO*E1);
                //ho2DArray<T> dst2D(RO, E1, const_cast<T*>(dst.begin())+ii*RO*E1);

                size_t ro, e1, row(0);
                long long x, y;

                for ( e1=startE1; e1<=endE1; e1++ )
                {
                    for ( ro=startRO; ro<=endRO; ro++ )
                    {

                        size_t colInd(0);
                        for ( y=-halfKy; y<=halfKy; y++ )
                        {
                            for ( x=-halfKx; x<=halfKx; x++ )
                            {
                                // A(row, colInd++) = src2D(ro+x, e1+y);
                                pA[row + colInd*rowA] = pSrc2D[ro+x + (e1+y)*RO];
                                colInd++;
                            }
                        }

                        // B(row, 0) = dst2D(ro, e1);
                        pB[row] = pDst2D[ro + e1*RO];

                        row++;
                    }
                }

                Gadgetron::SolveLinearSystem_Tikhonov(A, B, K, thresReg);

                memcpy(kernel.begin()+ii*kx*ky, K.begin(), sizeof(T)*kx*ky);
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::calibFengHuang(gtPlusReconWorkOrder2DT<T>& workOrder2DT, const hoNDArray<T>& src, const hoNDArray<T>& dst, ho6DArray<T>& kernel) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::performReconFangHuang(gtPlusReconWorkOrder2DT<T>& workOrder2DT, 
                                                const hoNDArray<T>& kspaceConj, hoNDArray<T>& kspace, 
                                                int startRO, int endRO, int startE1, int endE1, ho6DArray<T>& kernel)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(kspaceConj.dimensions_equal(&kspace));

        long long RO = (long long)kspace.get_size(0);
        long long E1 = (long long)kspace.get_size(1);
        long long CHA = (long long)kspace.get_size(2);
        long long N = (long long)kspace.get_size(3);
        long long S = (long long)kspace.get_size(4);

        long long kx = (long long)kernel.get_size(0);
        long long ky = (long long)kernel.get_size(1);

        long long halfKx = kx/2;
        long long halfKy = ky/2;
        long long kerN = (long long)kernel.get_size(4);
        GADGET_CHECK_RETURN_FALSE( (kerN==1) || (kerN==N) );

        long long num = CHA*N*S;

        long long rowD = RO*E1 - ( (endE1-startE1+1) * (endRO-startRO+1) );
        long long colD = kx*ky;

        ho2DArray<size_t> coeffX(rowD, colD);
        ho2DArray<size_t> coeffY(rowD, colD);

        long long ro, e1, row(0);
        long long x, y, dx, dy;

        for ( e1=0; e1<E1; e1++ )
        {
            for ( ro=0; ro<RO; ro++ )
            {
                if ( (ro>=startRO) && (ro<=endRO) && (e1>=startE1) && (e1<=endE1) )
                {
                    continue;
                }

                size_t colInd(0);
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
                        colInd++;
                    }
                }

                row++;
            }
        }

        long long ii;
        #pragma omp parallel default(none) private(ii) shared(num, RO, E1, CHA, N, S, kerN, kspaceConj, kspace, kernel, rowD, colD, coeffX, coeffY)
        {
            hoMatrix<T> D(rowD, colD);
            hoMatrix<T> K(colD, 1);
            hoMatrix<T> R(rowD, 1);

            Gadgetron::clear(D);
            Gadgetron::clear(K);
            Gadgetron::clear(R);

            #pragma omp for
            for ( ii=0; ii<num; ii ++ )
            {
                ho2DArray<T> src2D(RO, E1, const_cast<T*>(kspaceConj.begin())+ii*RO*E1);
                ho2DArray<T> dst2D(RO, E1, kspace.begin()+ii*RO*E1);

                long long row, col;
                for ( col=0; col<colD; col++ )
                {
                    for ( row=0; row<rowD; row++ )
                    {
                        D(row, col) = src2D(coeffX(row, col), coeffY(row, col));
                    }
                }

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
                Gadgetron::gemm(R, D, false, K, false);

                for ( row=0; row<rowD; row++ )
                {
                    dst2D( coeffX(row, colD/2), coeffY(row, colD/2) ) = R(row, 0);
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::performReconFangHuang(...) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusReconWorker2DT<T>::
estimateJobSize(gtPlusReconWorkOrder<T>* workOrder2DT, size_t maxNumOfBytesPerJob, size_t overlapBetweenJobs, size_t numOfNodes, size_t& jobSize)
{
    try
    {
        size_t nodeN = numOfNodes;
        GADGET_CHECK_RETURN_FALSE(this->computeEffectiveNodeNumberBasedOnComputingPowerIndex(workOrder2DT, nodeN));
        if ( workOrder2DT->job_perform_on_control_node_ ) nodeN++;

        GDEBUG_STREAM("GtPlus Cloud 2DT - job_perform_on_control_node is " << workOrder2DT->job_perform_on_control_node_  << " - nodeN is " << nodeN << " - overlapBetweenJobs is " << overlapBetweenJobs << " ... ");

        // adjust jobN according to cloud size
        size_t RO = workOrder2DT->data_.get_size(0);
        size_t E1 = workOrder2DT->data_.get_size(1);
        size_t N = workOrder2DT->data_.get_size(3);
        size_t S = workOrder2DT->data_.get_size(4);

        size_t srcCHA = workOrder2DT->kernelIm_->get_size(2);
        size_t dstCHA = workOrder2DT->kernelIm_->get_size(3);

        size_t totalJobNum = N;
        jobSize = (size_t)std::ceil( (double)(totalJobNum+overlapBetweenJobs*(nodeN-1))/(double)nodeN );

        size_t numOfBytesPerJob = sizeof(T)*( RO*E1*srcCHA*dstCHA*jobSize + 2*RO*E1*srcCHA*jobSize );

        // here a 64Mb graceful size is given to job
        while ( numOfBytesPerJob > maxNumOfBytesPerJob*1024*1024*1024-64.0*1024*1024 )
        {
            nodeN *= 2;
            jobSize = (size_t)std::ceil( (double)(totalJobNum+overlapBetweenJobs*(nodeN-1))/(double)nodeN );
            numOfBytesPerJob = sizeof(T)*( RO*E1*srcCHA*dstCHA*jobSize + 2*RO*E1*srcCHA*jobSize );
        }

        GDEBUG_STREAM("GtPlus Cloud 2DT - jobSize is " << jobSize << "; every job has " << numOfBytesPerJob/1024.0/1024 << " MBytes ... ");
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusReconWorker2DT<T>::estimateJobSize(...) ... ");
        return false;
    }

    return true;
}

}}
