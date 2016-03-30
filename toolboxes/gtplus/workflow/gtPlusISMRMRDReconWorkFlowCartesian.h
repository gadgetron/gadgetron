/** \file   gtPlusISMRMRDReconWorkFlowCartesian.h
    \brief  Define the base class for the GtPlus reconstruction workflow for cartesian sampling
    \author Hui Xue
*/

#pragma once

#include "gtPlusISMRMRDReconWorkFlow.h"

namespace Gadgetron { namespace gtPlus {

template <typename T> 
class gtPlusISMRMRDReconWorkFlowCartesian : public gtPlusISMRMRDReconWorkFlow<T>
{
public:

    typedef gtPlusISMRMRDReconWorkFlow<T> BaseClass;
    typedef typename BaseClass::DimensionRecordType DimensionRecordType;
    typedef typename BaseClass::real_value_type real_value_type;

    gtPlusISMRMRDReconWorkFlowCartesian();
    virtual ~gtPlusISMRMRDReconWorkFlowCartesian();

    void printInfo(std::ostream& os);

    virtual bool preProcessing();

    virtual bool postProcessing();
    virtual bool postProcessing(hoNDArray<T>& res, bool process_gfactor=true, bool process_wrap_around_map=true);

    virtual bool configureWorkOrder(const std::vector<ISMRMRDDIM>& dims);

    // resize or cut the reconstruected images to the recon space
    // res_ [RO E1 CHA SLC E2 ...]
    virtual bool convertToReconSpace2D(hoNDArray<T>& input_, hoNDArray<T>& output_, bool isKSpace);
    // res_ [RO E1 E2 CHA ...]
    virtual bool convertToReconSpace3D(hoNDArray<T>& input_, hoNDArray<T>& output_, bool isKSpace);

    // predict the workOrder dimensions
    virtual bool predictDimensions() = 0;

    using BaseClass::data_;
    using BaseClass::time_stamp_;
    using BaseClass::physio_time_stamp_;
    using BaseClass::ref_;
    using BaseClass::gfactor_;
    using BaseClass::wrap_around_map_;
    using BaseClass::noise_;
    using BaseClass::noiseBW_;
    using BaseClass::receriverBWRatio_;
    using BaseClass::ADCSamplingTimeinSecond_;
    using BaseClass::overSamplingRatioRO_;
    using BaseClass::reconSizeRO_;
    using BaseClass::reconSizeE1_;
    using BaseClass::reconSizeE2_;
    using BaseClass::encodingFOV_RO_;
    using BaseClass::encodingFOV_E1_;
    using BaseClass::encodingFOV_E2_;
    using BaseClass::reconFOV_RO_;
    using BaseClass::reconFOV_E1_;
    using BaseClass::reconFOV_E2_;

    using BaseClass::res_;
    using BaseClass::res_second_;
    using BaseClass::res_time_stamp_;
    using BaseClass::res_physio_time_stamp_;
    using BaseClass::res_time_stamp_second_;
    using BaseClass::res_physio_time_stamp_second_;

    using BaseClass::worker_;
    using BaseClass::workOrder_;

    using BaseClass::dimsRes_;

    using BaseClass::dataDimStartingIndexes_;

    using BaseClass::gt_timer1_;
    using BaseClass::gt_timer2_;
    using BaseClass::gt_timer3_;
    using BaseClass::performTiming_;
    using BaseClass::gt_exporter_;
    using BaseClass::debugFolder_;

    // the workOrder can share the kernel computation results
    // if this WorkOrderShareDim_ is not DIM_NONE, then 
    // workOrders will share kernel estimation results along this 
    // dimensions
    ISMRMRDDIM WorkOrderShareDim_;

    // work flow can buffer the kernel computed from previous work order and apply them to other work orders
    // work flow looks at the workFlow_BufferKernel_ and workFlow_use_BufferedKernel_ fields of work order
    // buffered kernels
    boost::shared_ptr< hoNDArray<T> > workFlowBuffer_kernel_;
    boost::shared_ptr< hoNDArray<T> > workFlowBuffer_kernelIm_;
    boost::shared_ptr< hoNDArray<T> > workFlowBuffer_unmixingCoeffIm_;
    boost::shared_ptr< hoNDArray<T> > workFlowBuffer_coilMap_;
    boost::shared_ptr< std::vector<hoMatrix<T> > > workFlowBuffer_coilCompressionCoef_;

    // whether to perform oversampling removal for ref data
    bool ref_remove_oversampling_RO_;
    // whether to apply noise prewhitening on ref data
    bool ref_apply_noisePreWhitening_;

protected:

    using BaseClass::dataCurr_;
    using BaseClass::refCurr_;
    using BaseClass::gfactorCurr_;
    using BaseClass::wrap_around_mapCurr_;

    using BaseClass::RO_;
    using BaseClass::E1_;
    using BaseClass::CHA_;
    using BaseClass::SLC_;
    using BaseClass::E2_;
    using BaseClass::CON_;
    using BaseClass::PHS_;
    using BaseClass::REP_;
    using BaseClass::SET_;
    using BaseClass::SEG_;
    using BaseClass::AVE_;

    using BaseClass::RO_ref_;
    using BaseClass::E1_ref_;
    using BaseClass::CHA_ref_;
    using BaseClass::SLC_ref_;
    using BaseClass::E2_ref_;
    using BaseClass::CON_ref_;
    using BaseClass::PHS_ref_;
    using BaseClass::REP_ref_;
    using BaseClass::SET_ref_;
    using BaseClass::SEG_ref_;
    using BaseClass::AVE_ref_;

    using BaseClass::gtPlus_util_;

    /// permute the array to the fixed order
    template <typename T2> 
    bool permuteArrayOrder(hoNDArray<T2>& data, std::vector<size_t>& order)
    {
        try
        {
            boost::shared_ptr< hoNDArray<T2> > data_permuted = Gadgetron::permute(&data, &order);
            data.reshape(data_permuted->get_dimensions());
            memcpy(data.begin(), data_permuted->begin(), data_permuted->get_number_of_bytes());
        }
        catch(...)
        {
            GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::permuteArrayOrder(hoNDArray<T>& data, const std::vector<int>& order) ... ");
            return false;
        }

        return true;
    }

    /// copy workOrder results to workflow
    bool copyReconResultsSecond(size_t dim5, size_t dim6, size_t dim7, size_t dim8, size_t dim9);
    virtual bool copyGFactor(size_t dim5, size_t dim6, size_t dim7, size_t dim8, size_t dim9, bool gfactor_needed);
    virtual bool copyWrapAroundMap(size_t dim5, size_t dim6, size_t dim7, size_t dim8, size_t dim9, bool wrap_around_map_needed);
};

template <typename T> 
gtPlusISMRMRDReconWorkFlowCartesian<T>::
gtPlusISMRMRDReconWorkFlowCartesian() : BaseClass(), WorkOrderShareDim_(DIM_NONE), ref_remove_oversampling_RO_(true), ref_apply_noisePreWhitening_(true)
{
}

template <typename T> 
gtPlusISMRMRDReconWorkFlowCartesian<T>::~gtPlusISMRMRDReconWorkFlowCartesian() 
{
}

template <typename T> 
void gtPlusISMRMRDReconWorkFlowCartesian<T>::printInfo(std::ostream& os)
{
    using namespace std;

    os << "-------------- GTPlus ISMRMRD Recon workflow Cartesian -------------" << endl;
    os << "Implementation of general reconstruction workflow for cartesian sampling" << endl;
    os << "Typical PreProcessing includes:" << endl;
    os << "a) Combine SEG dimension" << endl;
    os << "b) Remove readout oversampling if any" << endl;
    os << "c) If input noise scan is available, compute and apply the noise prewhitening matrix" << endl;
    os << "d) Apply the kspace filter along the RO direction if required" << endl;
    os << endl;
    os << "Typical PostProcessing includes:" << endl;
    os << "a) Apply the kspace filter along the E1 and E2 directions if required" << endl;
    os << "b) Perform the zero-padding resize if required" << endl;
    os << endl;
    os << "Data buffers are named to reflect the typical nature of MR acquisition" << endl;
    os << "data: image kspace data, 10D array [RO E1 CHA SLC E2 CON PHS REP SET SEG AVE]" << endl;
    os << "ref: calibration data, 10D array [RO E1 CHA SLC E2 CON PHS REP SET SEG AVE]" << endl;
    os << "----------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
preProcessing()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data_!=NULL);

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*data_, debugFolder_+"incomingKSpace"); }

        // combine the segment dimension
        if ( SEG_.second > 1 )
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(*data_, dataCurr_));

            std::vector<size_t> dim, dimCurr;
            data_->get_dimensions(dim);

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(*data_, dataCurr_, data_->get_number_of_dimensions()-1));

            dimCurr.resize(dim.size() - 1);
            memcpy(&dimCurr[0], &dim[0], sizeof(size_t)*dimCurr.size());
            dataCurr_.reshape(dimCurr);

            *data_ = dataCurr_;
            SEG_.second = 1;

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*data_, debugFolder_+"incomingKSpace_SEGCombined"); }
        }

        if ( (ref_ != NULL) && (ref_->get_number_of_elements()>0) )
        {
            if ( !debugFolder_.empty() )
            {
                gt_exporter_.exportArrayComplex(*ref_, debugFolder_+"incomingRef");
            }
        }

        if ( ref_!=NULL && SEG_ref_.second>1 )
        {
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(*ref_, refCurr_));

            std::vector<size_t> dim, dimCurr;
            ref_->get_dimensions(dim);

            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::sum_over_dimension(*ref_, refCurr_, ref_->get_number_of_dimensions() - 1));

            dimCurr.resize(dim.size() - 1);
            memcpy(&dimCurr[0], &dim[0], sizeof(size_t)*dimCurr.size());
            refCurr_.reshape(dimCurr);

            *ref_ = refCurr_;
            SEG_ref_.second = 1;
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*data_, debugFolder_+"incomingRef_SEGCombined"); }
        }

        // if needed, remove the readout oversampling
        if ( overSamplingRatioRO_ > 1.0 )
        {
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(*data_);
            // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().cutpad2D(*data_, (size_t)(data_->get_size(0)/overSamplingRatioRO_), data_->get_size(1), dataCurr_));
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop((size_t)(data_->get_size(0) / overSamplingRatioRO_), data_->get_size(1), data_, &dataCurr_));
            Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(dataCurr_);
            *data_ = dataCurr_;
            RO_.second = data_->get_size(0);

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*data_, debugFolder_+"kspace_oversamplingRORemoved"); }

            if ( ref_ != NULL && ref_remove_oversampling_RO_ )
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(*ref_);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop((size_t)(ref_->get_size(0) / overSamplingRatioRO_), ref_->get_size(1), ref_, &refCurr_));
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(refCurr_);
                *ref_ = refCurr_;
                RO_ref_.second = ref_->get_size(0);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*ref_, debugFolder_+"ref_oversamplingRORemoved"); }
            }

            if ( workOrder_->start_RO_>=0 && workOrder_->end_RO_>=0 )
            {
                workOrder_->start_RO_ = (int)(workOrder_->start_RO_/overSamplingRatioRO_);
                workOrder_->end_RO_ = (int)(workOrder_->end_RO_/overSamplingRatioRO_);
            }
        }

        // if needed, perform the noise prewhitening
        if ( noise_ != NULL )
        {
            hoMatrix<T> prewhiteningMatrix;
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().computeNoisePrewhiteningMatrix(*noise_, noiseBW_, receriverBWRatio_, ADCSamplingTimeinSecond_, prewhiteningMatrix));
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().performNoisePrewhitening(*data_, prewhiteningMatrix));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*data_, debugFolder_+"kspace_noiseprewhitenned"); }

            if ( ref_!=NULL && ref_apply_noisePreWhitening_ )
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().performNoisePrewhitening(*ref_, prewhiteningMatrix));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*ref_, debugFolder_+"ref_noiseprewhitenned"); }
            }
        }

        // if asymmetric echo is used, set the corresponding RO regions as 0
        size_t RO = data_->get_size(0);
        if ( !( workOrder_->start_RO_<0 || workOrder_->end_RO_<0 || (workOrder_->end_RO_-workOrder_->start_RO_+1==RO) ) )
        {
            size_t num = data_->get_number_of_elements() / RO;
            long long n;

            long long startRO = workOrder_->start_RO_;
            long long endRO = workOrder_->end_RO_;
            T* pData = data_->begin();

            #pragma omp parallel for default(none) private(n) shared(num, RO, startRO, endRO, pData)
            for ( n=0; n<(long long)num; n++ )
            {
                if ( startRO > 0 )
                {
                    memset(pData+n*RO, 0, startRO*sizeof(T) );
                }
                else
                {
                    memset(pData+n*RO+endRO+1, 0, (RO-endRO)*sizeof(T) );
                }
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*data_, debugFolder_+"incomingKSpace_RO_setzeros"); }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::preProcessing() ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
convertToReconSpace2D(hoNDArray<T>& input_, hoNDArray<T>& output_, bool isKSpace)
{
    try
    {
        size_t RO = data_->get_size(0);
        size_t E1 = data_->get_size(1);

        size_t inputRO = input_.get_size(0);
        size_t inputE1 = input_.get_size(1);

        output_ = input_;

        // if encoded FOV are the same as recon FOV
        if ( (std::abs(encodingFOV_RO_/2 - reconFOV_RO_)<0.1) && (std::abs(encodingFOV_E1_-reconFOV_E1_)<0.1) )
        {
            hoNDArray<T> buffer2D;

            if ( isKSpace )
            {
                if (reconSizeE1_ >= E1)
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2DOnKSpace(input_, RO, reconSizeE1_, buffer2D));
                }
                else
                {
                    hoNDArray<T> bufferIm2D;
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(input_, bufferIm2D);
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, reconSizeE1_, &bufferIm2D, &buffer2D));
                }

                if (reconSizeRO_ >= RO)
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2D(buffer2D, reconSizeRO_, reconSizeE1_, output_));
                }
                else
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, reconSizeE1_, &buffer2D, &output_));
                }
            }
            else
            {
                if (reconSizeE1_ >= E1)
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2D(input_, RO, reconSizeE1_, buffer2D));
                }
                else
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, reconSizeE1_, &input_, &buffer2D));
                }

                if(reconSizeRO_>=RO)
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2D(buffer2D, reconSizeRO_, reconSizeE1_, output_));
                }
                else
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, reconSizeE1_, &buffer2D, &output_));
                }
            }
        }
        else if (encodingFOV_E1_>=reconFOV_E1_)
        {
            size_t encodingE1 = reconSizeE1_;
            if ( encodingFOV_E1_ > reconFOV_E1_ )
            {
                float spacingE1 = reconFOV_E1_/reconSizeE1_;
                encodingE1 = (size_t)std::floor(encodingFOV_E1_/spacingE1+0.5);
            }

            hoNDArray<T>* pSrc = &input_;
            hoNDArray<T>* pDst = &output_;
            hoNDArray<T>* pTmp;

            hoNDArray<T> buffer2D;

            // adjust E1
            if ( encodingE1>E1 && encodingE1>inputE1 )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2DOnKSpace(*pSrc, RO, encodingE1, *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2D(*pSrc, RO, encodingE1, *pDst));
                }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*pDst, debugFolder_+"complexIm_zpadResize2D_enlarged"); }

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }
            else if ( encodingE1 < E1 )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, encodingE1, pSrc, pDst));
                }
                else
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(*pSrc, buffer2D);
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, encodingE1, &buffer2D, pDst));
                }

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(*pDst);

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*pDst, "complexIm_zpadResize2D_cut"); }

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }

            //adjust RO
            if ( RO < reconSizeRO_ )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2DOnKSpace(*pSrc, reconSizeRO_, pSrc->get_size(1), *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2D(*pSrc, reconSizeRO_, pSrc->get_size(1), *pDst));
                }

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }
            else if ( RO > reconSizeRO_ )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, pSrc->get_size(1), pSrc, pDst));
                }
                else
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(*pSrc, buffer2D);
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, pSrc->get_size(1), &buffer2D, pDst));
                }

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(*pDst);

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }

            // final cut
            if ( isKSpace )
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(*pSrc, buffer2D);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, reconSizeE1_, &buffer2D, pDst));
            }
            else
            {
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, reconSizeE1_, pSrc, pDst));
            }

            if ( pDst != &output_ )
            {
                output_ = *pDst;
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::convertToReconSpace2D(const hoNDArray& input_, hoNDArray& output_, bool isKSpace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
convertToReconSpace3D(hoNDArray<T>& input_, hoNDArray<T>& output_, bool isKSpace)
{
    try
    {
        size_t RO = input_.get_size(0);
        size_t E1 = input_.get_size(1);
        size_t E2 = input_.get_size(2);

        output_ = input_;

        // if encoded FOV are the same as recon FOV
        if ( (std::abs(encodingFOV_RO_/2 - reconFOV_RO_)<0.1) && (std::abs(encodingFOV_E1_-reconFOV_E1_)<0.1) && (std::abs(encodingFOV_E2_-reconFOV_E2_)<0.1) )
        {
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3DOnKSpace(input_, reconSizeRO_, reconSizeE1_, reconSizeE2_, output_));
            }
            else
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3D(input_, reconSizeRO_, reconSizeE1_, reconSizeE2_, output_));
            }
        }
        else if ( (encodingFOV_E1_>=reconFOV_E1_) && (encodingFOV_E2_>=reconFOV_E2_) )
        {
            size_t encodingE1 = reconSizeE1_;
            if ( encodingFOV_E1_ > reconFOV_E1_ )
            {
                float spacingE1 = reconFOV_E1_/reconSizeE1_;
                encodingE1 = (size_t)std::floor(encodingFOV_E1_/spacingE1+0.5);
            }

            size_t encodingE2 = reconSizeE2_;
            if ( encodingFOV_E2_ > reconFOV_E2_ )
            {
                float spacingE2 = reconFOV_E2_/reconSizeE2_;
                encodingE2 = (size_t)std::floor(encodingFOV_E2_/spacingE2+0.5);
            }

            hoNDArray<T>* pSrc = &input_;
            hoNDArray<T>* pDst = &output_;
            hoNDArray<T>* pTmp;

            hoNDArray<T> buffer3D;

            // adjust E1
            if ( encodingE1 >= E1+1 )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3DOnKSpace(*pSrc, RO, encodingE1, E2, *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3D(*pSrc, RO, encodingE1, E2, *pDst));
                }

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }
            else if ( encodingE1 <= E1-1 )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, encodingE1, E2, pSrc, pDst));
                }
                else
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pSrc, buffer3D);
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, encodingE1, E2, &buffer3D, pDst));
                }

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst);

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }

            // adjust E2
            if ( encodingE2 >= E2+1 )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3DOnKSpace(*pSrc, RO, pSrc->get_size(1), encodingE2, *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3D(*pSrc, RO, pSrc->get_size(1), encodingE2, *pDst));
                }

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }
            else if ( encodingE2 <= E2-1 )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, pSrc->get_size(1), encodingE2, pSrc, pDst));
                }
                else
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pSrc, buffer3D);
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(RO, pSrc->get_size(1), encodingE2, &buffer3D, pDst));
                }

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst);

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }

            //adjust RO
            if ( RO < reconSizeRO_ )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3DOnKSpace(*pSrc, reconSizeRO_, pSrc->get_size(1), pSrc->get_size(2), *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize3D(*pSrc, reconSizeRO_, pSrc->get_size(1), pSrc->get_size(2), *pDst));
                }

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }
            else if ( RO > reconSizeRO_ )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, pSrc->get_size(1), pSrc->get_size(2), pSrc, pDst));
                }
                else
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pSrc, buffer3D);
                    GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, pSrc->get_size(1), pSrc->get_size(2), &buffer3D, pDst));
                }

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst);

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*pSrc, debugFolder_ + "res_beforeCut"); }

            // final cut on image
            if ( isKSpace )
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pSrc, buffer3D);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, reconSizeE1_, reconSizeE2_, &buffer3D, pDst));
            }
            else
            {
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::crop(reconSizeRO_, reconSizeE1_, reconSizeE2_, pSrc, pDst));
            }
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(*pDst, debugFolder_+"res_AfterCut"); }

            if ( pDst != &output_ )
            {
                output_ = *pDst;
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::convertToReconSpace3D(const hoNDArray& input_, hoNDArray& output_, bool isKSpace) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
postProcessing(hoNDArray<T>& res, bool process_gfactor, bool process_wrap_around_map)
{
    try
    {
        size_t RO = res.get_size(0);
        size_t E1 = res.get_size(1);
        size_t E2 = res.get_size(4);

        // whether to process gfactor
        bool has_gfactor = false;
        if ( (gfactor_.get_size(0)==RO) && (gfactor_.get_size(1)==E1) )
        {
            has_gfactor = true;
        }

        if ( !process_gfactor )
        {
            has_gfactor = false;
        }

        // whehter to process wrap_around map
        bool has_wrap_around = false;
        if ( (wrap_around_map_.get_size(0)==RO) && (wrap_around_map_.get_size(1)==E1) )
        {
            has_wrap_around = true;
        }

        if ( !process_wrap_around_map )
        {
            has_wrap_around = false;
        }

        if ( E2_.second > 1 )
        {
            if ( performTiming_ ) { gt_timer1_.start("postProcessing - permute res array ... "); }

            // permute E2 dimension from [RO E1 CHA SLC E2 ...] to [RO E1 E2 CHA SLC ...]

            std::vector<size_t> dim_order(5);
            dim_order[0] = 0;
            dim_order[1] = 1;
            dim_order[2] = 4;
            dim_order[3] = 2;
            dim_order[4] = 3;

            std::vector<size_t> dim, dimPermuted;
            res.get_dimensions(dim);
            dimPermuted = dim;
            dimPermuted[2] = dim[4];
            dimPermuted[3] = dim[2];
            dimPermuted[4] = dim[3];

            dataCurr_.create(dimPermuted);
            GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&res, &dataCurr_, &dim_order));

            if ( has_gfactor )
            {
                gfactor_.get_dimensions(dim);
                dimPermuted = dim;
                dimPermuted[2] = dim[4];
                dimPermuted[3] = dim[2];
                dimPermuted[4] = dim[3];

                gfactorCurr_.create(dimPermuted);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&gfactor_, &gfactorCurr_, &dim_order));
            }

            if ( has_wrap_around )
            {
                wrap_around_map_.get_dimensions(dim);
                dimPermuted = dim;
                dimPermuted[2] = dim[4];
                dimPermuted[3] = dim[2];
                dimPermuted[4] = dim[3];

                wrap_around_mapCurr_.create(dimPermuted);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&wrap_around_map_, &wrap_around_mapCurr_, &dim_order));
            }

            if ( performTiming_ ) { gt_timer1_.stop(); }

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(dataCurr_, debugFolder_+"data_permuted"); }

            // dataCurr_ = *data_permuted;

            res.reshape(dataCurr_.get_dimensions());

            bool inKSpace = false;

            if ( workOrder_->filterROE1E2_.get_size(0)==RO 
                    && workOrder_->filterROE1E2_.get_size(1)==E1 
                    && workOrder_->filterROE1E2_.get_size(2)==E2 )
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(dataCurr_, res);
                Gadgetron::apply_kspace_filter_ROE1E2(res, workOrder_->filterROE1E2_, dataCurr_);
                inKSpace = true;
            }
            else if ( (workOrder_->filterRO_.get_number_of_elements() == RO) 
                        && (workOrder_->filterE1_.get_number_of_elements() == E1) 
                        && (workOrder_->filterE2_.get_number_of_elements() == E2) )
            {
                if ( performTiming_ ) { gt_timer1_.start("postProcessing - fft3c ... "); }
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(dataCurr_, res);
                if ( performTiming_ ) { gt_timer1_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"kspace_beforefiltered"); }

                if ( performTiming_ ) { gt_timer1_.start("postProcessing - 3D kspace filter ... "); }
                Gadgetron::apply_kspace_filter_ROE1E2(res, workOrder_->filterRO_, workOrder_->filterE1_, workOrder_->filterE2_, dataCurr_);
                if ( performTiming_ ) { gt_timer1_.stop(); }

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(dataCurr_, debugFolder_+"kspace_afterfiltered"); }
                inKSpace = true;
            }
            else
            {
                hoNDArray<T>* pSrc = &res;
                hoNDArray<T>* pDst = &dataCurr_;

                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pDst, *pSrc);

                bool filterPerformed = false;

                if ( workOrder_->filterRO_.get_number_of_elements() == RO )
                {
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterRO(*pSrc, workOrder_->filterRO_, *pDst));
                    Gadgetron::apply_kspace_filter_RO(*pSrc, workOrder_->filterRO_, *pDst);
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if ( workOrder_->filterE1_.get_number_of_elements() == E1 )
                {
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterE1(*pSrc, workOrder_->filterE1_, *pDst));
                    Gadgetron::apply_kspace_filter_E1(*pSrc, workOrder_->filterE1_, *pDst);
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if ( workOrder_->filterE2_.get_number_of_elements() == E2 )
                {
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspace3DfilterE2(*pSrc, workOrder_->filterE2_, *pDst));
                    Gadgetron::apply_kspace_filter_E2(*pSrc, workOrder_->filterE2_, *pDst);
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if ( filterPerformed )
                {
                    if ( pDst != &dataCurr_ )
                    {
                        dataCurr_ = *pDst;
                    }
                }
                else
                {
                    dataCurr_ = res;
                }

                inKSpace = true;
            }

            if ( inKSpace )
            {
                if ( !debugFolder_.empty() )
                {
                    hoNDArray<T> Im(dataCurr_);
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(Im);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(Im, debugFolder_+"complexIm_filtered"); }
                }
            }
            else
            {
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"complexIm_filtered"); }
            }

            GADGET_CHECK_RETURN_FALSE(convertToReconSpace3D(dataCurr_, res, inKSpace));

            {
                std::vector<size_t> dim_order(5);
                dim_order[0] = 0;
                dim_order[1] = 1;
                dim_order[2] = 3;
                dim_order[3] = 4;
                dim_order[4] = 2;

                std::vector<size_t> dim, dimPermuted;
                res.get_dimensions(dim);
                dimPermuted = dim;
                dimPermuted[2] = dim[3];
                dimPermuted[3] = dim[4];
                dimPermuted[4] = dim[2];

                dataCurr_.create(dimPermuted);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&res, &dataCurr_, &dim_order));
            }

            res.reshape(dataCurr_.get_dimensions());
            memcpy(res.begin(), dataCurr_.begin(), res.get_number_of_bytes());

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"complexIm_zpadResize3D"); }

            if ( has_gfactor )
            {
                GADGET_CHECK_RETURN_FALSE(convertToReconSpace3D(gfactorCurr_, gfactor_, false));

                std::vector<size_t> dim_order(5);
                dim_order[0] = 0;
                dim_order[1] = 1;
                dim_order[2] = 3;
                dim_order[3] = 4;
                dim_order[4] = 2;

                std::vector<size_t> dim, dimPermuted;
                gfactor_.get_dimensions(dim);
                dimPermuted = dim;
                dimPermuted[2] = dim[3];
                dimPermuted[3] = dim[4];
                dimPermuted[4] = dim[2];

                gfactorCurr_.create(dimPermuted);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&gfactor_, &gfactorCurr_, &dim_order));

                gfactor_.reshape(gfactorCurr_.get_dimensions());
                memcpy(gfactor_.begin(), gfactorCurr_.begin(), gfactor_.get_number_of_bytes());

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(gfactor_, debugFolder_+"gfactor_zpadResize3D"); }
            }

            if ( has_wrap_around )
            {
                GADGET_CHECK_RETURN_FALSE(convertToReconSpace3D(wrap_around_mapCurr_, wrap_around_map_, false));

                std::vector<size_t> dim_order(5);
                dim_order[0] = 0;
                dim_order[1] = 1;
                dim_order[2] = 3;
                dim_order[3] = 4;
                dim_order[4] = 2;

                std::vector<size_t> dim, dimPermuted;
                wrap_around_map_.get_dimensions(dim);
                dimPermuted = dim;
                dimPermuted[2] = dim[3];
                dimPermuted[3] = dim[4];
                dimPermuted[4] = dim[2];

                wrap_around_mapCurr_.create(dimPermuted);
                GADGET_CHECK_EXCEPTION_RETURN_FALSE(Gadgetron::permute(&wrap_around_map_, &wrap_around_mapCurr_, &dim_order));

                wrap_around_map_.reshape(wrap_around_mapCurr_.get_dimensions());
                memcpy(wrap_around_map_.begin(), wrap_around_mapCurr_.begin(), wrap_around_map_.get_number_of_bytes());

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(wrap_around_map_, debugFolder_+"wrap_around_map_zpadResize3D"); }
            }
        }
        else
        {
            dataCurr_ = res;
            bool inKSpace = false;

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(dataCurr_, debugFolder_+"complexIm_before_filtered"); }

            if ( workOrder_->filterROE1_.get_size(0)==RO && workOrder_->filterROE1_.get_size(1)==E1 )
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res);
                // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterROE1(res, workOrder_->filterROE1_, dataCurr_));
                Gadgetron::apply_kspace_filter_ROE1(res, workOrder_->filterROE1_, dataCurr_);
                inKSpace = true;
            }
            else if ( (workOrder_->filterRO_.get_number_of_elements() == RO) && (workOrder_->filterE1_.get_number_of_elements() == E1) )
            {
                Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res);
                // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterROE1(res, workOrder_->filterRO_, workOrder_->filterE1_, dataCurr_));
                Gadgetron::apply_kspace_filter_ROE1(res, workOrder_->filterRO_, workOrder_->filterE1_, dataCurr_);
                inKSpace = true;
            }
            else
            {
                if ( (workOrder_->filterRO_.get_number_of_elements() == RO) && (workOrder_->filterE1_.get_number_of_elements() != E1) )
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res);
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterRO(res, workOrder_->filterRO_, dataCurr_));
                    Gadgetron::apply_kspace_filter_RO(res, workOrder_->filterRO_, dataCurr_);
                    inKSpace = true;
                }

                if ( (workOrder_->filterRO_.get_number_of_elements() != RO) && (workOrder_->filterE1_.get_number_of_elements() == E1) )
                {
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res);
                    // GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterE1(res, workOrder_->filterE1_, dataCurr_));
                    Gadgetron::apply_kspace_filter_E1(res, workOrder_->filterE1_, dataCurr_);
                    inKSpace = true;
                }
            }

            if ( inKSpace )
            {
                if ( !debugFolder_.empty() )
                {
                    hoNDArray<T> Im(dataCurr_);
                    Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(Im);
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(Im, debugFolder_+"complexIm_after_filtered"); }
                }
            }
            else
            {
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"complexIm_after_filtered"); }
            }

            GADGET_CHECK_RETURN_FALSE(convertToReconSpace2D(dataCurr_, res, inKSpace));

            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res, debugFolder_+"complexIm_zpadResize2D"); }

            if ( has_gfactor )
            {
                gfactorCurr_ = gfactor_;
                GADGET_CHECK_RETURN_FALSE(convertToReconSpace2D(gfactorCurr_, gfactor_, false));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(gfactor_, debugFolder_+"gfactor_zpadResize2D"); }
            }

            if ( has_wrap_around )
            {
                wrap_around_mapCurr_ = wrap_around_map_;
                GADGET_CHECK_RETURN_FALSE(convertToReconSpace2D(wrap_around_mapCurr_, wrap_around_map_, false));

                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(wrap_around_map_, debugFolder_+"wrap_around_map_zpadResize2D"); }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::postProcessing(res) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
postProcessing()
{
    try
    {
        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res_, debugFolder_+"complexIm_afterRecon"); }
        GADGET_CHECK_RETURN_FALSE(this->postProcessing(res_, true, true));
        if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(res_, debugFolder_ + "complexIm_afterRecon_afterPostProcessing"); }

        if ( this->res_second_.get_number_of_elements() > 0 )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res_second_, debugFolder_+"complexImSecond_afterRecon"); }
            GADGET_CHECK_RETURN_FALSE(this->postProcessing(res_second_, false, false));
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::postProcessing() ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
configureWorkOrder(const std::vector<ISMRMRDDIM>& dims)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data_!=NULL);
        GADGET_CHECK_RETURN_FALSE(worker_!=NULL);
        GADGET_CHECK_RETURN_FALSE(workOrder_!=NULL);

        if ( ref_ == NULL )
        {
            ref_ = data_;
        }

        size_t dd;

        // find the dimension size for data and ref
        size_t E2 = data_->get_size(4);

        std::vector<size_t> dimSize(dims.size());
        std::vector<size_t> dimSizeRef(dims.size(), 1);
        size_t indChannelDim = 2;
        for ( dd=0; dd<dims.size(); dd++ )
        {
            dimSize[dd] = data_->get_size(dims[dd]-DIM_ReadOut);
            if ( ref_ != NULL )
            {
                dimSizeRef[dd] = ref_->get_size(dims[dd]-DIM_ReadOut);
            }

            if ( dims[dd] == DIM_Channel )
            {
                indChannelDim = dd;
            }
        }

        GDEBUG_CONDITION_STREAM(!debugFolder_.empty(), "Recon dimensions : " << this->printISMRMRDDimensions(dims));
        GDEBUG_CONDITION_STREAM(!debugFolder_.empty(), "Recon size       : " << this->printISMRMRDDimensionSize(dimSize));
        GDEBUG_CONDITION_STREAM(!debugFolder_.empty(), "Recon ref size   : " << this->printISMRMRDDimensionSize(dimSizeRef));

        bool gfactor_needed = workOrder_->gfactor_needed_;
        bool wrap_around_map_needed = workOrder_->wrap_around_map_needed_;

        // recon workOrder size
        std::vector<size_t> dimReconSize(5);
        dimReconSize[0] = dimSize[0];
        dimReconSize[1] = dimSize[1];
        dimReconSize[2] = dimSize[2];
        dimReconSize[3] = dimSize[3];
        dimReconSize[4] = dimSize[4];

        std::vector<size_t> dimReconTimeStampSize(dimReconSize);
        dimReconTimeStampSize[0] = 1; // RO = 1
        dimReconTimeStampSize[2] = 1; // CHA = 1

        std::vector<size_t> dimReconSizeRef(5);
        dimReconSizeRef[0] = dimSizeRef[0];
        dimReconSizeRef[1] = dimSizeRef[1];
        dimReconSizeRef[2] = dimSizeRef[2];
        dimReconSizeRef[3] = dimSizeRef[3];
        dimReconSizeRef[4] = dimSizeRef[4];

        // first two dimension are always RO and E1
        size_t N2D = dimReconSize[0]*dimReconSize[1];
        size_t N2DRef = dimReconSizeRef[0]*dimReconSizeRef[1];

        size_t N3D = N2D*dimReconSize[2];
        size_t N3DRef = N2DRef*dimReconSizeRef[2];

        // allocate the results
        size_t num_channels_res = workOrder_->num_channels_res_;

        std::vector<size_t> dimResSize(dimSize);

        if ( gfactor_needed )
        {
            dimResSize[indChannelDim] = 1;
            gfactor_.create(&dimResSize);
        }

        if ( wrap_around_map_needed )
        {
            if ( workOrder_->acceFactorE2_ > 1 ) // 3D acquisition
            {
                dimResSize[indChannelDim] = 3;
            }
            else
            {
                dimResSize[indChannelDim] = 2;
            }

            wrap_around_map_.create(&dimResSize);
        }

        dimResSize[indChannelDim] = num_channels_res;
        res_.create(&dimResSize);
        Gadgetron::clear(res_);

        res_second_.create(&dimResSize);
        Gadgetron::clear(res_second_);

        std::vector<size_t> dimReconResTimeStampSize(dimResSize);
        dimReconResTimeStampSize[0] = 1;
        dimReconResTimeStampSize[1] = 1;
        dimReconResTimeStampSize[2] = 1;

        res_time_stamp_.create(dimReconResTimeStampSize);
        Gadgetron::fill(res_time_stamp_, (real_value_type)(-1) );

        res_physio_time_stamp_.create(dimReconResTimeStampSize);
        Gadgetron::fill(res_physio_time_stamp_, (real_value_type)(-1) );

        res_time_stamp_second_.create(dimReconResTimeStampSize);
        Gadgetron::fill(res_time_stamp_second_, (real_value_type)(-1) );

        res_physio_time_stamp_second_.create(dimReconResTimeStampSize);
        Gadgetron::fill(res_physio_time_stamp_second_, (real_value_type)(-1) );

        std::vector<ISMRMRDDIM> dimsRes(dims);

        GDEBUG_CONDITION_STREAM(!debugFolder_.empty(), "Recon res dimensions : " << this->printISMRMRDDimensions(dimsRes));
        GDEBUG_CONDITION_STREAM(!debugFolder_.empty(), "Recon res size       : " << this->printISMRMRDDimensionSize(dimResSize));

        bool shareAcrossWorkOrders = (WorkOrderShareDim_!=DIM_NONE);

        if ( !debugFolder_.empty() )
        {
            gt_exporter_.exportArrayComplex(*data_, debugFolder_ + "data_");
            gt_exporter_.exportArrayComplex(*ref_, debugFolder_ + "ref_");

            if ( time_stamp_ != NULL )
            {
                gt_exporter_.exportArray(*time_stamp_, debugFolder_ + "time_stamp_");
            }

            if ( physio_time_stamp_ != NULL )
            {
                gt_exporter_.exportArray(*physio_time_stamp_, debugFolder_ + "physio_time_stamp_");
            }
        }

        bool workFlow_use_BufferedKernel_ = workOrder_->workFlow_use_BufferedKernel_;

        bool has_second_res = false;
        bool has_recon_time_stamp = false;
        bool has_recon_physio_time_stamp = false;
        bool has_recon_time_stamp_second = false;
        bool has_recon_physio_time_stamp_second = false;

        size_t numOfRecon = dimSize[9] * dimSize[8] * dimSize[7] * dimSize[6] * dimSize[5];

        // call up the recon
        size_t num_recon = 0;

        size_t dim9, dim8, dim7, dim6, dim5, dim4, dim3, dim2;
        for ( dim9=0; dim9<dimSize[9]; dim9++ )
        {
            for ( dim8=0; dim8<dimSize[8]; dim8++ )
            {
                for ( dim7=0; dim7<dimSize[7]; dim7++ )
                {
                    for ( dim6=0; dim6<dimSize[6]; dim6++ )
                    {
                        for ( dim5=0; dim5<dimSize[5]; dim5++ )
                        {
                            std::vector<size_t> ind(11, 0);
                            this->ismrmrdDimIndex10D(ind, dims[9], dim9);
                            this->ismrmrdDimIndex10D(ind, dims[8], dim8);
                            this->ismrmrdDimIndex10D(ind, dims[7], dim7);
                            this->ismrmrdDimIndex10D(ind, dims[6], dim6);
                            this->ismrmrdDimIndex10D(ind, dims[5], dim5);

                            // ---------------------------
                            // prepare the data in workOrder
                            // ---------------------------
                            if ( !workOrder_->data_.dimensions_equal(&dimReconSize) )
                            {
                                workOrder_->data_.create(&dimReconSize);
                                workOrder_->time_stamp_.create(&dimReconTimeStampSize);
                                Gadgetron::clear(workOrder_->time_stamp_);

                                workOrder_->physio_time_stamp_.create(&dimReconTimeStampSize);
                                Gadgetron::clear(workOrder_->physio_time_stamp_);
                            }

                            std::vector<size_t> indWorkOrder(5, 0);
                            for ( dim4=0; dim4<dimSize[4]; dim4++ )
                            {
                                this->ismrmrdDimIndex10D(ind, dims[4], dim4);
                                indWorkOrder[4] = dim4;

                                for ( dim3=0; dim3<dimSize[3]; dim3++ )
                                {
                                    this->ismrmrdDimIndex10D(ind, dims[3], dim3);
                                    indWorkOrder[3] = dim3;

                                    if ( dims[2] == DIM_Channel )
                                    {
                                        long long offset = data_->calculate_offset(ind);
                                        long long offsetWorkOrder = workOrder_->data_.calculate_offset(indWorkOrder);
                                        memcpy(workOrder_->data_.begin()+offsetWorkOrder, data_->begin()+offset, sizeof(T)*N3D);

                                        if ( time_stamp_ != NULL )
                                        {
                                            offset = time_stamp_->calculate_offset(ind);
                                            offsetWorkOrder = workOrder_->time_stamp_.calculate_offset(indWorkOrder);
                                            memcpy(workOrder_->time_stamp_.begin()+offsetWorkOrder, time_stamp_->begin()+offset, sizeof(real_value_type)*dimReconSize[1]);
                                            if ( physio_time_stamp_ != NULL )
                                            {
                                                memcpy(workOrder_->physio_time_stamp_.begin()+offsetWorkOrder, physio_time_stamp_->begin()+offset, sizeof(real_value_type)*dimReconSize[1]);
                                            }
                                        }
                                    }
                                    else
                                    {
                                        GWARN_STREAM("dims[2] != DIM_Channel, the time stamps will not be copied ... ");

                                        for ( dim2=0; dim2<dimSize[2]; dim2++ )
                                        {
                                            this->ismrmrdDimIndex10D(ind, dims[2], dim2);
                                            indWorkOrder[2] = dim2;

                                            long long offset = data_->calculate_offset(ind);
                                            long long offsetWorkOrder = workOrder_->data_.calculate_offset(indWorkOrder);
                                            memcpy(workOrder_->data_.begin()+offsetWorkOrder, data_->begin()+offset, sizeof(T)*N2D);
                                        }
                                    }
                                }
                            }

                            // ---------------------------
                            // prepare the ref in workOrder
                            // ---------------------------
                            if ( (ref_ != NULL) && (ref_->get_number_of_elements()>0) )
                            {
                                std::vector<size_t> indRef(11, 0);

                                if ( dim9 < dimSizeRef[9] )
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[9], dim9);
                                }
                                else
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[9], dimSizeRef[9]-1);
                                }

                                if ( dim8 < dimSizeRef[8] )
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[8], dim8);
                                }
                                else
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[8], dimSizeRef[8]-1);
                                }

                                if ( dim7 < dimSizeRef[7] )
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[7], dim7);
                                }
                                else
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[7], dimSizeRef[7]-1);
                                }

                                if ( dim6 < dimSizeRef[6] )
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[6], dim6);
                                }
                                else
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[6], dimSizeRef[6]-1);
                                }

                                if ( dim5 < dimSizeRef[5] )
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[5], dim5);
                                }
                                else
                                {
                                    this->ismrmrdDimIndex10D(indRef, dims[5], dimSizeRef[5]-1);
                                }

                                if ( !workOrder_->ref_.dimensions_equal(&dimReconSizeRef) )
                                {
                                    workOrder_->ref_.create(&dimReconSizeRef);
                                }

                                std::vector<size_t> indRefWorkOrder(11, 0);
                                for ( dim4=0; dim4<dimSize[4]; dim4++ )
                                {
                                    size_t dim4_ref = dim4;
                                    if ( dim4 < dimSizeRef[4] )
                                    {
                                        this->ismrmrdDimIndex10D(indRef, dims[4], dim4);
                                    }
                                    else
                                    {
                                        this->ismrmrdDimIndex10D(indRef, dims[4], dimSizeRef[4]-1);
                                        dim4_ref = dimSizeRef[4]-1;
                                    }
                                    indRefWorkOrder[4] = dim4_ref;

                                    for ( dim3=0; dim3<dimSize[3]; dim3++ )
                                    {
                                        size_t dim3_ref = dim3;
                                        if ( dim3 < dimSizeRef[3] )
                                        {
                                            this->ismrmrdDimIndex10D(indRef, dims[3], dim3);
                                        }
                                        else
                                        {
                                            this->ismrmrdDimIndex10D(indRef, dims[3], dimSizeRef[3]-1);
                                            dim3_ref = dimSizeRef[3]-1;
                                        }
                                        indRefWorkOrder[3] = dim3_ref;

                                        if ( dims[2] == DIM_Channel )
                                        {
                                            long long offset = ref_->calculate_offset(indRef);
                                            long long offsetWorkOrder = workOrder_->ref_.calculate_offset(indRefWorkOrder);
                                            memcpy(workOrder_->ref_.begin()+offsetWorkOrder, ref_->begin()+offset, sizeof(T)*N3DRef);
                                        }
                                        else
                                        {
                                            for ( dim2=0; dim2<dimSize[2]; dim2++ )
                                            {
                                                size_t dim2_ref = dim2;
                                                if ( dim2 < dimSizeRef[2] )
                                                {
                                                    this->ismrmrdDimIndex10D(indRef, dims[2], dim2);
                                                }
                                                else
                                                {
                                                    this->ismrmrdDimIndex10D(indRef, dims[2], dimSizeRef[2]-1);
                                                    dim2_ref = dimSizeRef[2]-1;
                                                }
                                                indRefWorkOrder[2] = dim2_ref;

                                                long long offset = ref_->calculate_offset(indRef);
                                                long long offsetWorkOrder = workOrder_->ref_.calculate_offset(indRefWorkOrder);
                                                memcpy(workOrder_->ref_.begin()+offsetWorkOrder, ref_->begin()+offset, sizeof(T)*N2DRef);
                                            }
                                        }
                                    }
                                }
                            }

                            // ---------------------------
                            // handle shared work order
                            // ---------------------------
                            if ( !shareAcrossWorkOrders && workOrder_->workFlow_BufferKernel_ && !workOrder_->workFlow_use_BufferedKernel_ )
                            {
                                GADGET_CHECK_RETURN_FALSE(workOrder_->reset());
                            }

                            if ( shareAcrossWorkOrders && !workOrder_->workFlow_use_BufferedKernel_ )
                            {
                                if ( dim5==0 )
                                {
                                    workOrder_->workFlow_use_BufferedKernel_ = false;
                                }
                                else
                                {
                                    workOrder_->workFlow_use_BufferedKernel_ = true;
                                }
                            }

                            if ( !debugFolder_.empty() )
                            {
                                gt_exporter_.exportArrayComplex(workOrder_->data_, debugFolder_+"workOrder_data");
                                gt_exporter_.exportArray(workOrder_->time_stamp_, debugFolder_+"workOrder_time_stamp");
                                gt_exporter_.exportArray(workOrder_->physio_time_stamp_, debugFolder_+"workOrder_physio_time_stamp");
                                gt_exporter_.exportArrayComplex(workOrder_->ref_, debugFolder_+"workOrder_ref");
                            }

                            // ---------------------------
                            // clean the input array
                            // ---------------------------
                            num_recon++;

                            if (num_recon == numOfRecon && E2>1)
                            {
                                this->data_->clear();
                                this->ref_->clear();
                            }

                            // ---------------------------
                            // perform the recon
                            // ---------------------------
                            GADGET_CHECK_RETURN_FALSE(worker_->performRecon(workOrder_));

                            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder_->complexIm_, debugFolder_+"workOrder_complexIm"); }

                            if ( workOrder_->complexIm_second_.get_number_of_elements()>0 )
                            {
                                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(workOrder_->complexIm_second_, debugFolder_+"workOrder_complexImSecond"); }
                            }

                            if ( shareAcrossWorkOrders )
                            {
                                workOrder_->workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;
                            }

                            // ---------------------------
                            // copy the recon complexIm
                            // ---------------------------
                            std::vector<size_t> indRes(ind);
                            indRes[0] = 0;
                            indRes[1] = 0;
                            indRes[2] = 0;
                            indRes[3] = 0;
                            indRes[4] = 0;
                            indRes[5] = dim5;
                            indRes[6] = dim6;
                            indRes[7] = dim7;
                            indRes[8] = dim8;
                            indRes[9] = dim9;

                            long long offset = res_.calculate_offset(indRes);
                            memcpy(res_.begin()+offset, workOrder_->complexIm_.begin(), workOrder_->complexIm_.get_number_of_bytes());

                            // ---------------------------
                            // copy the recon time stamp
                            // ---------------------------
                            if ( workOrder_->recon_time_stamp_.get_number_of_elements()>0 )
                            {
                                has_recon_time_stamp = true;

                                offset = res_time_stamp_.calculate_offset(indRes);
                                memcpy(res_time_stamp_.begin()+offset, workOrder_->recon_time_stamp_.begin(), workOrder_->recon_time_stamp_.get_number_of_bytes());
                            }

                            // ---------------------------
                            // copy the recon physio time stamp
                            // ---------------------------
                            if ( workOrder_->recon_physio_time_stamp_.get_number_of_elements()>0 )
                            {
                                has_recon_physio_time_stamp = true;

                                offset = res_physio_time_stamp_.calculate_offset(indRes);
                                memcpy(res_physio_time_stamp_.begin()+offset, workOrder_->recon_physio_time_stamp_.begin(), workOrder_->recon_physio_time_stamp_.get_number_of_bytes());
                            }

                            // ---------------------------
                            // copy the second set of recon complexIm
                            // ---------------------------
                            GADGET_CHECK_RETURN_FALSE(this->copyReconResultsSecond(dim5, dim6, dim7, dim8, dim9));

                            if ( workOrder_->complexIm_second_.get_number_of_elements()>0 )
                            {
                                has_second_res = true;
                            }

                            if ( workOrder_->recon_time_stamp_second_.get_number_of_elements()>0 )
                            {
                                has_recon_time_stamp_second = true;
                            }

                            if ( workOrder_->recon_physio_time_stamp_second_.get_number_of_elements()>0 )
                            {
                                has_recon_physio_time_stamp_second = true;
                            }

                            // ---------------------------
                            // copy the gfactor
                            // ---------------------------
                            GADGET_CHECK_RETURN_FALSE(this->copyGFactor(dim5, dim6, dim7, dim8, dim9, gfactor_needed));

                            // ---------------------------
                            // copy the wrap-round map
                            // ---------------------------
                            GADGET_CHECK_RETURN_FALSE(this->copyWrapAroundMap(dim5, dim6, dim7, dim8, dim9, wrap_around_map_needed));

                            // if not sharing across work order
                            if ( !shareAcrossWorkOrders && !workOrder_->workFlow_use_BufferedKernel_ && !workOrder_->workFlow_BufferKernel_ )
                            {
                                GADGET_CHECK_RETURN_FALSE(workOrder_->reset());
                            }
                        }

                        // in the outter dimensions, the work order is always reset
                        if ( !workOrder_->workFlow_use_BufferedKernel_ && !workOrder_->workFlow_BufferKernel_ )
                        {
                            GADGET_CHECK_RETURN_FALSE(workOrder_->reset());
                        }
                    }
                }
            }
        }

        if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res_, debugFolder_+"res_afterunwrapping"); }

        if ( has_second_res )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res_second_, debugFolder_+"res_second_afterunwrapping"); }
        }

        if ( has_recon_time_stamp )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_time_stamp_, debugFolder_+"res_time_stamp"); }
        }

        if ( has_recon_physio_time_stamp )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_physio_time_stamp_, debugFolder_+"res_physio_time_stamp"); }
        }

        if ( has_recon_time_stamp_second )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_time_stamp_second_, debugFolder_+"res_time_stamp_second"); }
        }

        if ( has_recon_physio_time_stamp_second )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_physio_time_stamp_second_, debugFolder_+"res_physio_time_stamp_second"); }
        }

        if ( gfactor_needed )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(gfactor_, debugFolder_+"gfactor_afterunwrapping"); }
        }

        if ( wrap_around_map_needed )
        {
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(wrap_around_map_, debugFolder_+"wrap_around_map_afterunwrapping"); }
        }

        // permute the res_ to the correct dimension order
        if (   ( (res_.get_number_of_elements()>dimResSize[0]*dimResSize[1]) && (dims[2]!=DIM_Channel) ) 
            || ( (res_.get_number_of_elements()>dimResSize[0]*dimResSize[1]*dimResSize[2])             ) )
        {
            std::vector<size_t> order;
            GADGET_CHECK_RETURN_FALSE(this->findISMRMRDPermuteOrder(dimsRes, dimsRes_, order));

            GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(res_, order));
            if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res_, debugFolder_+"res_afterPermute"); }

            if ( has_recon_time_stamp )
            {
                GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(res_time_stamp_, order));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_time_stamp_, debugFolder_+"res_time_stamp_afterPermute"); }
            }

            if ( has_recon_physio_time_stamp )
            {
                GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(res_physio_time_stamp_, order));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_physio_time_stamp_, debugFolder_+"res_physio_time_stamp_afterPermute"); }
            }

            if ( gfactor_needed )
            {
                GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(gfactor_, order));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(gfactor_, debugFolder_+"gfactor_afterPermute"); }
            }

            if ( wrap_around_map_needed )
            {
                GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(wrap_around_map_, order));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(wrap_around_map_, debugFolder_+"wrap_around_map_afterPermute"); }
            }
        }

        if ( has_second_res )
        {
            if (   ( (res_second_.get_number_of_elements()>dimResSize[0]*dimResSize[1]) && (dims[2]!=DIM_Channel) ) 
                || ( (res_second_.get_number_of_elements()>dimResSize[0]*dimResSize[1]*dimResSize[2])             ) )
            {
                std::vector<size_t> order;
                GADGET_CHECK_RETURN_FALSE(this->findISMRMRDPermuteOrder(dimsRes, dimsRes_, order));

                GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(res_second_, order));
                if ( !debugFolder_.empty() ) { gt_exporter_.exportArrayComplex(res_second_, debugFolder_+"res_second_afterPermute"); }

                if ( has_recon_time_stamp_second )
                {
                    GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(res_time_stamp_second_, order));
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_time_stamp_, debugFolder_+"res_time_stamp_second_afterPermute"); }
                }

                if ( has_recon_physio_time_stamp_second )
                {
                    GADGET_CHECK_RETURN_FALSE(this->permuteArrayOrder(res_physio_time_stamp_second_, order));
                    if ( !debugFolder_.empty() ) { gt_exporter_.exportArray(res_physio_time_stamp_second_, debugFolder_+"res_physio_time_stamp_second_afterPermute"); }
                }
            }
        }
        else
        {
            res_second_.clear();
        }

        if ( !has_recon_time_stamp )
        {
            res_time_stamp_.clear();
        }

        if ( !has_recon_physio_time_stamp )
        {
            res_physio_time_stamp_.clear();
        }

        if ( !has_recon_time_stamp_second )
        {
            res_time_stamp_second_.clear();
        }

        if ( !has_recon_physio_time_stamp_second )
        {
            res_physio_time_stamp_second_.clear();
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::configureWorkOrder(const std::vector<ISMRMRDDIM>& dims) ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
copyReconResultsSecond(size_t dim5, size_t dim6, size_t dim7, size_t dim8, size_t dim9)
{
    try
    {

        if ( workOrder_->complexIm_second_.get_number_of_elements()>0 )
        {
            std::vector<size_t> indRes(10);

            size_t RO = workOrder_->complexIm_second_.get_size(0);
            size_t E1 = workOrder_->complexIm_second_.get_size(1);
            size_t N = workOrder_->complexIm_second_.get_size(2);
            size_t S = workOrder_->complexIm_second_.get_size(3);

            std::vector<size_t> dims;

            bool hasTimeStamp = false;
            if ( workOrder_->recon_time_stamp_second_.get_number_of_elements()>0 )
            {
                hasTimeStamp = true;

                res_time_stamp_second_.get_dimensions(dims);
                if ( dims[3] != N ) dims[3] = N;
                if ( dims[4] != S ) dims[4] = S;

                res_time_stamp_second_.create(dims);
                Gadgetron::clear(res_time_stamp_second_);
            }

            bool hasPhysioTimeStamp = false;
            if ( workOrder_->recon_physio_time_stamp_second_.get_number_of_elements()>0 )
            {
                hasPhysioTimeStamp = true;

                res_physio_time_stamp_second_.get_dimensions(dims);
                if ( dims[3] != N ) dims[3] = N;
                if ( dims[4] != S ) dims[4] = S;

                res_physio_time_stamp_second_.create(dims);
                Gadgetron::clear(res_physio_time_stamp_second_);
            }

            res_second_.get_dimensions(dims);
            if ( dims[3] != N ) dims[3] = N;
            if ( dims[4] != S ) dims[4] = S;

            res_second_.create(dims);
            Gadgetron::clear(res_second_);

            size_t n, s;
            for ( s=0; s<S; s++ )
            {
                for ( n=0; n<N; n++ )
                {
                    indRes[0] = 0;
                    indRes[1] = 0;
                    indRes[2] = 0;
                    indRes[3] = n;
                    indRes[4] = s;
                    indRes[5] = dim5;
                    indRes[6] = dim6;
                    indRes[7] = dim7;
                    indRes[8] = dim8;
                    indRes[9] = dim9;

                    size_t offset = res_second_.calculate_offset(indRes);
                    memcpy(res_second_.begin()+offset, workOrder_->complexIm_second_.begin()+n*RO*E1+s*RO*E1*N, sizeof(T)*RO*E1);

                    if ( hasTimeStamp )
                    {
                        offset = res_time_stamp_second_.calculate_offset(indRes);
                        res_time_stamp_second_(offset) = workOrder_->recon_time_stamp_second_(0, 0, 0, n, s);
                    }

                    if ( hasPhysioTimeStamp )
                    {
                        offset = res_physio_time_stamp_second_.calculate_offset(indRes);
                        res_physio_time_stamp_second_(offset) = workOrder_->recon_physio_time_stamp_second_(0, 0, 0, n, s);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::copyReconResultsSecond() ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
copyGFactor(size_t dim5, size_t dim6, size_t dim7, size_t dim8, size_t dim9, bool gfactor_needed)
{
    try
    {
        if ( gfactor_needed && (workOrder_->gfactor_.get_size(0)==res_.get_size(0)) && (workOrder_->gfactor_.get_size(1) == res_.get_size(1)) )
        {
            size_t RO = gfactor_.get_size(0);
            size_t E1 = gfactor_.get_size(1);
            size_t N = gfactor_.get_size(3);
            size_t S = gfactor_.get_size(4);

            size_t gfactor_N = workOrder_->gfactor_.get_size(2);
            size_t gfactor_S = workOrder_->gfactor_.get_size(3);

            if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(workOrder_->gfactor_, debugFolder_ + "workOrder_gfactor_afterunwrapping"); }

            std::vector<size_t> indRes(10);
            indRes[0] = 0;
            indRes[1] = 0;
            indRes[2] = 0;
            indRes[3] = 0;
            indRes[4] = 0;
            indRes[5] = dim5;
            indRes[6] = dim6;
            indRes[7] = dim7;
            indRes[8] = dim8;
            indRes[9] = dim9;

            if ( (gfactor_N == N) && (gfactor_S == S) )
            {
                size_t offset = gfactor_.calculate_offset(indRes);
                memcpy(gfactor_.begin()+offset, workOrder_->gfactor_.begin(), workOrder_->gfactor_.get_number_of_bytes());
            }
            else
            {
                std::vector<size_t> indGfactor(9);
                indGfactor[0] = 0;
                indGfactor[1] = 0;
                indGfactor[2] = 0;
                indGfactor[3] = 0;
                indGfactor[4] = dim5;
                indGfactor[5] = dim6;
                indGfactor[6] = dim7;
                indGfactor[7] = dim8;
                indGfactor[8] = dim9;

                size_t n, s;
                for ( s=0; s<S; s++ )
                {
                    for ( n=0; n<N; n++ )
                    {
                        indRes[3] = n;
                        indRes[4] = s;
                        size_t offset = gfactor_.calculate_offset(indRes);

                        if ( n < gfactor_N )
                        {
                            indGfactor[2] = n;
                        }
                        else
                        {
                            indGfactor[2] = gfactor_N-1;
                        }

                        if ( s < gfactor_S )
                        {
                            indGfactor[3] = s;
                        }
                        else
                        {
                            indGfactor[3] = gfactor_S-1;
                        }

                        size_t offset2 = workOrder_->gfactor_.calculate_offset(indGfactor);

                        memcpy(gfactor_.begin()+offset, workOrder_->gfactor_.begin()+offset2, sizeof(T)*RO*E1);
                    }
                }
            }

            if (!debugFolder_.empty()) { gt_exporter_.exportArrayComplex(gfactor_, debugFolder_ + "gfactor_after_copyGFactor"); }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::copyGFactor() ... ");
        return false;
    }

    return true;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
copyWrapAroundMap(size_t dim5, size_t dim6, size_t dim7, size_t dim8, size_t dim9, bool wrap_around_map_needed)
{
    try
    {
        if ( wrap_around_map_needed && (workOrder_->wrap_around_map_.get_size(0)==res_.get_size(0)) && (workOrder_->wrap_around_map_.get_size(1) == res_.get_size(1)) )
        {
            size_t RO = wrap_around_map_.get_size(0);
            size_t E1 = wrap_around_map_.get_size(1);
            size_t N = wrap_around_map_.get_size(3);
            size_t S = wrap_around_map_.get_size(4);

            size_t wrap_around_map_CHA = workOrder_->wrap_around_map_.get_size(2);
            size_t wrap_around_map_N = workOrder_->wrap_around_map_.get_size(3);
            size_t wrap_around_map_S = workOrder_->wrap_around_map_.get_size(4);

            std::vector<size_t> indRes(10);
            size_t offset;

            indRes[0] = 0;
            indRes[1] = 0;
            indRes[2] = 0;
            indRes[3] = 0;
            indRes[4] = 0;
            indRes[5] = dim5;
            indRes[6] = dim6;
            indRes[7] = dim7;
            indRes[8] = dim8;
            indRes[9] = dim9;

            if ( (wrap_around_map_N == N) && (wrap_around_map_S == S) )
            {
                offset = wrap_around_map_.calculate_offset(indRes);
                memcpy(wrap_around_map_.begin()+offset, workOrder_->wrap_around_map_.begin(), workOrder_->wrap_around_map_.get_number_of_bytes());
            }
            else
            {
                std::vector<size_t> indWrapAroundMap(10);
                indWrapAroundMap[0] = 0;
                indWrapAroundMap[1] = 0;
                indWrapAroundMap[2] = 0;
                indWrapAroundMap[3] = 0;
                indWrapAroundMap[4] = 0;
                indWrapAroundMap[5] = dim5;
                indWrapAroundMap[6] = dim6;
                indWrapAroundMap[7] = dim7;
                indWrapAroundMap[8] = dim8;
                indWrapAroundMap[9] = dim9;

                size_t n, s;
                for ( s=0; s<S; s++ )
                {
                    for ( n=0; n<N; n++ )
                    {
                        indRes[3] = n;
                        indRes[4] = s;
                        offset = wrap_around_map_.calculate_offset(indRes);

                        if ( n < wrap_around_map_N )
                        {
                            indWrapAroundMap[3] = n;
                        }
                        else
                        {
                            indWrapAroundMap[3] = wrap_around_map_N-1;
                        }

                        if ( s < wrap_around_map_S )
                        {
                            indWrapAroundMap[4] = s;
                        }
                        else
                        {
                            indWrapAroundMap[4] = wrap_around_map_S-1;
                        }

                        size_t offset2 = workOrder_->wrap_around_map_.calculate_offset(indWrapAroundMap);

                        memcpy(wrap_around_map_.begin()+offset, workOrder_->wrap_around_map_.begin()+offset2, sizeof(T)*RO*E1*wrap_around_map_CHA);
                    }
                }
            }
        }
    }
    catch(...)
    {
        GERROR_STREAM("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::copyWrapAroundMap() ... ");
        return false;
    }

    return true;
}

}}
