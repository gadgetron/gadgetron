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

    gtPlusISMRMRDReconWorkFlowCartesian();
    virtual ~gtPlusISMRMRDReconWorkFlowCartesian();

    void printInfo(std::ostream& os);

    virtual bool preProcessing();

    virtual bool postProcessing();

    virtual bool configureWorkOrder(const std::vector<ISMRMRDDIM>& dims);

    // resize or cut the reconstruected images to the recon space
    // res_ [RO E1 CHA SLC E2 ...]
    virtual bool convertToReconSpace2D(hoNDArray<T>& input_, hoNDArray<T>& output_, bool isKSpace);
    // res_ [RO E1 E2 CHA ...]
    virtual bool convertToReconSpace3D(hoNDArray<T>& input_, hoNDArray<T>& output_, bool isKSpace);

    // predict the workOrder dimensions
    virtual bool predictDimensions() = 0;

    using BaseClass::data_;
    using BaseClass::ref_;
    using BaseClass::gfactor_;
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

    using BaseClass::gtPlus_util_;
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
    os << "data: image kspace data, 10D array [RO E1 CHA SLC E2 CON PHS REP SET SEG]" << endl;
    os << "ref: calibration data, 10D array [RO E1 CHA SLC E2 CON PHS REP SET SEG]" << endl;
    os << "----------------------------------------------------------" << endl;
}

template <typename T> 
bool gtPlusISMRMRDReconWorkFlowCartesian<T>::
preProcessing()
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(data_!=NULL);

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *data_, "incomingKSpace");

        // combine the segment dimension
        if ( SEG_.second > 1 )
        {
            GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(*data_, dataCurr_));
            *data_ = dataCurr_;
            SEG_.second = 1;

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *data_, "incomingKSpace_SEGCombined");
        }

        if ( (ref_ != NULL) && (ref_->get_number_of_elements()>0) ) { GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *ref_, "incomingRef"); }

        if ( ref_!=NULL && SEG_ref_.second>1 )
        {
            GADGET_CHECK_RETURN_FALSE(Gadgetron::sumOverLastDimension(*ref_, refCurr_));
            *ref_ = refCurr_;
            SEG_ref_.second = 1;
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *data_, "incomingRef_SEGCombined");
        }

        // if needed, remove the readout oversampling
        if ( overSamplingRatioRO_ > 1.0 )
        {
            GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(*data_));
            GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().cutpad2D(*data_, (size_t)(data_->get_size(0)/overSamplingRatioRO_), data_->get_size(1), dataCurr_));
            GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(dataCurr_));
            *data_ = dataCurr_;
            RO_.second = data_->get_size(0);

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *data_, "kspace_oversamplingRORemoved");

            if ( ref_ != NULL && ref_remove_oversampling_RO_ )
            {
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft1c(*ref_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().cutpad2D(*ref_, (size_t)(ref_->get_size(0)/overSamplingRatioRO_), ref_->get_size(1), refCurr_));
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft1c(refCurr_));
                *ref_ = refCurr_;
                RO_ref_.second = ref_->get_size(0);

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *ref_, "ref_oversamplingRORemoved");
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

            // GADGET_CHECK_PERFORM(!debugFolder_.empty(), prewhiteningMatrix.print(std::cout));

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *data_, "kspace_noiseprewhitenned");

            if ( ref_!=NULL && ref_apply_noisePreWhitening_ )
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().performNoisePrewhitening(*ref_, prewhiteningMatrix));
                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *ref_, "ref_noiseprewhitenned");
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

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *data_, "incomingKSpace_RO_setzeros");
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::preProcessing() ... ");
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
        if ( (GT_ABS(encodingFOV_RO_/2 - reconFOV_RO_)<0.1) && (GT_ABS(encodingFOV_E1_-reconFOV_E1_)<0.1) )
        {
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2DOnKSpace(input_, reconSizeRO_, reconSizeE1_, output_));
            }
            else
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().zpadResize2D(input_, reconSizeRO_, reconSizeE1_, output_));
            }
        }
        else if (encodingFOV_E1_>=reconFOV_E1_)
        {
            size_t encodingE1 = reconSizeE1_;
            if ( encodingFOV_E1_ > reconFOV_E1_ )
            {
                float spacingE1 = reconFOV_E1_/reconSizeE1_;
                encodingE1 = (size_t)(encodingFOV_E1_/spacingE1);
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

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }
            else if ( encodingE1 < E1 )
            {
                if ( isKSpace )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad2D(*pSrc, RO, encodingE1, *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(*pSrc, buffer2D));
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad2D(buffer2D, RO, encodingE1, *pDst));
                }

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(*pDst));

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
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad2D(*pSrc, reconSizeRO_, pSrc->get_size(1), *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(*pSrc, buffer2D));
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad2D(buffer2D, reconSizeRO_, pSrc->get_size(1), *pDst));
                }

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(*pDst));

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }

            // final cut
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad2D(*pSrc, reconSizeRO_, reconSizeE1_, *pDst));
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(*pDst));
            }
            else
            {
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(*pSrc, buffer2D));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad2D(buffer2D, reconSizeRO_, reconSizeE1_, *pDst));
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(*pDst));
            }

            if ( pDst != &output_ )
            {
                output_ = *pDst;
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::convertToReconSpace2D(const hoNDArray& input_, hoNDArray& output_, bool isKSpace) ... ");
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
        size_t RO = res_.get_size(0);
        size_t E1 = res_.get_size(1);
        size_t E2 = res_.get_size(2);

        output_ = input_;

        // if encoded FOV are the same as recon FOV
        if ( (GT_ABS(encodingFOV_RO_/2 - reconFOV_RO_)<0.1) && (GT_ABS(encodingFOV_E1_-reconFOV_E1_)<0.1) && (GT_ABS(encodingFOV_E2_-reconFOV_E2_)<0.1) )
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
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(*pSrc, RO, encodingE1, E2, *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pSrc, buffer3D));
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(buffer3D, RO, encodingE1, E2, *pDst));
                }

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst));

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
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(*pSrc, RO, pSrc->get_size(1), encodingE2, *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pSrc, buffer3D));
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(buffer3D, RO, pSrc->get_size(1), encodingE2, *pDst));
                }

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst));

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
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(*pSrc, reconSizeRO_, pSrc->get_size(1), pSrc->get_size(2), *pDst));
                }
                else
                {
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pSrc, buffer3D));
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(buffer3D, reconSizeRO_, pSrc->get_size(1), pSrc->get_size(2), *pDst));
                }

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst));

                isKSpace = false;
                pTmp = pSrc; pSrc = pDst; pDst = pTmp;
            }

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *pSrc, "res_beforeCut");

            // final cut on image
            if ( isKSpace )
            {
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pSrc, buffer3D));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(buffer3D, reconSizeRO_, reconSizeE1_, reconSizeE2_, *pDst));
            }
            else
            {
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtilComplex<T>().cutpad3D(*pSrc, reconSizeRO_, reconSizeE1_, reconSizeE2_, *pDst));
            }
            // GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(*pDst));

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, *pDst, "res_AfterCut");

            if ( pDst != &output_ )
            {
                output_ = *pDst;
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::convertToReconSpace3D(const hoNDArray& input_, hoNDArray& output_, bool isKSpace) ... ");
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
        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "complexIm_afterRecon");

        size_t RO = res_.get_size(0);
        size_t E1 = res_.get_size(1);
        size_t E2 = res_.get_size(4);

        bool has_gfactor = false;
        if ( (gfactor_.get_size(0)==RO) && (gfactor_.get_size(1)==E1) )
        {
            has_gfactor = true;
        }

        if ( E2_.second > 1 )
        {
            // dataCurr_ = res_;

            // need to permute the matrix order
            //size_t NDim = dataCurr_.get_number_of_dimensions();
            //std::vector<size_t> order(NDim, 1);

            //size_t ii;
            //for ( ii=0; ii<NDim; ii++ )
            //{
            //    order[ii] = ii;
            //}

            //order[0] = 0;
            //order[1] = 1;
            //order[2] = 4;
            //order[3] = 2;
            //order[4] = 3;

            GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.start("postProcessing - permute res array ... "));
            // boost::shared_ptr< hoNDArray<T> > data_permuted = Gadgetron::permute(const_cast<hoNDArray<T>*>(&dataCurr_), &order);
            GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteE2To3rdDimension(res_, dataCurr_));

            if ( has_gfactor )
            {
                GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteE2To3rdDimension(gfactor_, gfactorCurr_));
            }

            GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.stop());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, dataCurr_, "data_permuted");

            // dataCurr_ = *data_permuted;

            res_.reshape(dataCurr_.get_dimensions());

            bool inKSpace = false;

            if ( workOrder_->filterROE1E2_.get_size(0)==RO 
                    && workOrder_->filterROE1E2_.get_size(1)==E1 
                    && workOrder_->filterROE1E2_.get_size(2)==E2 )
            {
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(dataCurr_, res_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspace3DfilterROE1E2(res_, workOrder_->filterROE1E2_, dataCurr_));
                inKSpace = true;
            }
            else if ( (workOrder_->filterRO_.get_number_of_elements() == RO) 
                        && (workOrder_->filterE1_.get_number_of_elements() == E1) 
                        && (workOrder_->filterE2_.get_number_of_elements() == E2) )
            {
                GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.start("postProcessing - fft3c ... "));
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(dataCurr_, res_));
                GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.stop());

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "kspace_beforefiltered");

                GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.start("postProcessing - 3D kspace filter ... "));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspace3DfilterROE1E2(res_, workOrder_->filterRO_, workOrder_->filterE1_, workOrder_->filterE2_, dataCurr_));
                GADGET_CHECK_PERFORM(performTiming_, gt_timer1_.stop());

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, dataCurr_, "kspace_afterfiltered");
                inKSpace = true;
            }
            else
            {
                hoNDArray<T>* pSrc = &res_;
                hoNDArray<T>* pDst = &dataCurr_;

                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft3c(*pDst, *pSrc));

                bool filterPerformed = false;

                if ( workOrder_->filterRO_.get_number_of_elements() == RO )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterRO(*pSrc, workOrder_->filterRO_, *pDst));
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if ( workOrder_->filterE1_.get_number_of_elements() == E1 )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterE1(*pSrc, workOrder_->filterE1_, *pDst));
                    std::swap(pSrc, pDst);
                    filterPerformed = true;
                }

                if ( workOrder_->filterE2_.get_number_of_elements() == E2 )
                {
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspace3DfilterE2(*pSrc, workOrder_->filterE2_, *pDst));
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
                    dataCurr_ = res_;
                }

                inKSpace = true;
            }

            if ( inKSpace )
            {
                if ( !debugFolder_.empty() )
                {
                    hoNDArray<T> Im(dataCurr_);
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft3c(Im));
                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, Im, "complexIm_filtered");
                }
            }
            else
            {
                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "complexIm_filtered");
            }

            GADGET_CHECK_RETURN_FALSE(convertToReconSpace3D(dataCurr_, res_, inKSpace));

            GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteE2To5thDimension(res_, dataCurr_));

            //order[0] = 0;
            //order[1] = 1;
            //order[2] = 3;
            //order[3] = 4;
            //order[4] = 2;

            //data_permuted = Gadgetron::permute(const_cast<hoNDArray<T>*>(&res_), &order);
            //res_ = *data_permuted;

            res_.reshape(dataCurr_.get_dimensions());
            memcpy(res_.begin(), dataCurr_.begin(), res_.get_number_of_bytes());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "complexIm_zpadResize3D");

            if ( has_gfactor )
            {
                GADGET_CHECK_RETURN_FALSE(convertToReconSpace3D(gfactorCurr_, gfactor_, false));
                GADGET_CHECK_RETURN_FALSE(Gadgetron::permuteE2To5thDimension(gfactor_, gfactorCurr_));

                gfactor_.reshape(gfactorCurr_.get_dimensions());
                memcpy(gfactor_.begin(), gfactorCurr_.begin(), gfactor_.get_number_of_bytes());

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, gfactor_, "gfactor_zpadResize3D");
            }
        }
        else
        {
            dataCurr_ = res_;
            bool inKSpace = false;

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, dataCurr_, "complexIm_before_filtered");

            if ( workOrder_->filterROE1_.get_size(0)==RO && workOrder_->filterROE1_.get_size(1)==E1 )
            {
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterROE1(res_, workOrder_->filterROE1_, dataCurr_));
                inKSpace = true;
            }
            else if ( (workOrder_->filterRO_.get_number_of_elements() == RO) && (workOrder_->filterE1_.get_number_of_elements() == E1) )
            {
                GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res_));
                GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterROE1(res_, workOrder_->filterRO_, workOrder_->filterE1_, dataCurr_));
                inKSpace = true;
            }
            else
            {
                if ( (workOrder_->filterRO_.get_number_of_elements() == RO) && (workOrder_->filterE1_.get_number_of_elements() != E1) )
                {
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res_));
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterRO(res_, workOrder_->filterRO_, dataCurr_));
                    inKSpace = true;
                }

                if ( (workOrder_->filterRO_.get_number_of_elements() != RO) && (workOrder_->filterE1_.get_number_of_elements() == E1) )
                {
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->fft2c(dataCurr_, res_));
                    GADGET_CHECK_RETURN_FALSE(gtPlusISMRMRDReconUtil<T>().kspacefilterE1(res_, workOrder_->filterE1_, dataCurr_));
                    inKSpace = true;
                }
            }

            if ( inKSpace )
            {
                if ( !debugFolder_.empty() )
                {
                    hoNDArray<T> Im(dataCurr_);
                    GADGET_CHECK_RETURN_FALSE(Gadgetron::hoNDFFT<typename realType<T>::Type>::instance()->ifft2c(Im));
                    GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, Im, "complexIm_after_filtered");
                }
            }
            else
            {
                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "complexIm_after_filtered");
            }

            GADGET_CHECK_RETURN_FALSE(convertToReconSpace2D(dataCurr_, res_, inKSpace));

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "complexIm_zpadResize2D");

            if ( has_gfactor )
            {
                gfactorCurr_ = gfactor_;
                GADGET_CHECK_RETURN_FALSE(convertToReconSpace2D(gfactorCurr_, gfactor_, false));

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, gfactor_, "gfactor_zpadResize2D");
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::postProcessing() ... ");
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

        GADGET_CONDITION_MSG(!debugFolder_.empty(), "Recon dimensions : " << this->printISMRMRDDimensions(dims));
        GADGET_CONDITION_MSG(!debugFolder_.empty(), "Recon size       : " << this->printISMRMRDDimensionSize(dimSize));
        GADGET_CONDITION_MSG(!debugFolder_.empty(), "Recon ref size   : " << this->printISMRMRDDimensionSize(dimSizeRef));

        bool gfactor_needed = workOrder_->gfactor_needed_;

        // recon workOrder size
        std::vector<size_t> dimReconSize(5);
        dimReconSize[0] = dimSize[0];
        dimReconSize[1] = dimSize[1];
        dimReconSize[2] = dimSize[2];
        dimReconSize[3] = dimSize[3];
        dimReconSize[4] = dimSize[4];

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

        dimResSize[indChannelDim] = num_channels_res;
        res_.create(&dimResSize);

        std::vector<ISMRMRDDIM> dimsRes(dims);

        GADGET_CONDITION_MSG(!debugFolder_.empty(), "Recon res dimensions : " << this->printISMRMRDDimensions(dimsRes));
        GADGET_CONDITION_MSG(!debugFolder_.empty(), "Recon res size       : " << this->printISMRMRDDimensionSize(dimResSize));

        bool shareAcrossWorkOrders = (WorkOrderShareDim_!=DIM_NONE);

        if ( !debugFolder_.empty() )
        {
            gt_exporter_.exportArrayComplex(*data_, debugFolder_ + "data_");
            gt_exporter_.exportArrayComplex(*ref_, debugFolder_ + "ref_");
        }

        bool workFlow_use_BufferedKernel_ = workOrder_->workFlow_use_BufferedKernel_;

        // call up the recon
        size_t dim8, dim7, dim6, dim5, dim4, dim3, dim2;
        for ( dim8=0; dim8<dimSize[8]; dim8++ )
        {
            for ( dim7=0; dim7<dimSize[7]; dim7++ )
            {
                for ( dim6=0; dim6<dimSize[6]; dim6++ )
                {
                    for ( dim5=0; dim5<dimSize[5]; dim5++ )
                    {
                        std::vector<size_t> ind(10, 0);
                        this->ismrmrdDimIndex9D(ind, dims[8], dim8);
                        this->ismrmrdDimIndex9D(ind, dims[7], dim7);
                        this->ismrmrdDimIndex9D(ind, dims[6], dim6);
                        this->ismrmrdDimIndex9D(ind, dims[5], dim5);

                        if ( !workOrder_->data_.dimensions_equal(&dimReconSize) )
                        {
                            workOrder_->data_.create(&dimReconSize);
                        }

                        std::vector<size_t> indWorkOrder(5, 0);
                        for ( dim4=0; dim4<dimSize[4]; dim4++ )
                        {
                            this->ismrmrdDimIndex9D(ind, dims[4], dim4);
                            indWorkOrder[4] = dim4;

                            for ( dim3=0; dim3<dimSize[3]; dim3++ )
                            {
                                this->ismrmrdDimIndex9D(ind, dims[3], dim3);
                                indWorkOrder[3] = dim3;

                                if ( dims[2] == DIM_Channel )
                                {
                                    long long offset = data_->calculate_offset(ind);

                                    long long offsetWorkOrder = workOrder_->data_.calculate_offset(indWorkOrder);

                                    memcpy(workOrder_->data_.begin()+offsetWorkOrder, data_->begin()+offset, sizeof(T)*N3D);
                                }
                                else
                                {
                                    for ( dim2=0; dim2<dimSize[2]; dim2++ )
                                    {
                                        this->ismrmrdDimIndex9D(ind, dims[2], dim2);
                                        indWorkOrder[2] = dim2;

                                        long long offset = data_->calculate_offset(ind);

                                        long long offsetWorkOrder = workOrder_->data_.calculate_offset(indWorkOrder);

                                        memcpy(workOrder_->data_.begin()+offsetWorkOrder, data_->begin()+offset, sizeof(T)*N2D);
                                    }
                                }
                            }
                        }

                        if ( (ref_ != NULL) && (ref_->get_number_of_elements()>0) )
                        {
                            std::vector<size_t> indRef(10, 0);
                            if ( dim8 < dimSizeRef[8] )
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[8], dim8);
                            }
                            else
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[8], dimSizeRef[8]-1);
                            }

                            if ( dim7 < dimSizeRef[7] )
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[7], dim7);
                            }
                            else
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[7], dimSizeRef[7]-1);
                            }

                            if ( dim6 < dimSizeRef[6] )
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[6], dim6);
                            }
                            else
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[6], dimSizeRef[6]-1);
                            }

                            if ( dim5 < dimSizeRef[5] )
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[5], dim5);
                            }
                            else
                            {
                                this->ismrmrdDimIndex9D(indRef, dims[5], dimSizeRef[5]-1);
                            }

                            if ( !workOrder_->ref_.dimensions_equal(&dimReconSizeRef) )
                            {
                                workOrder_->ref_.create(&dimReconSizeRef);
                            }

                            std::vector<size_t> indRefWorkOrder(10, 0);
                            for ( dim4=0; dim4<dimSize[4]; dim4++ )
                            {
                                size_t dim4_ref = dim4;
                                if ( dim4 < dimSizeRef[4] )
                                {
                                    this->ismrmrdDimIndex9D(indRef, dims[4], dim4);
                                }
                                else
                                {
                                    this->ismrmrdDimIndex9D(indRef, dims[4], dimSizeRef[4]-1);
                                    dim4_ref = dimSizeRef[4]-1;
                                }
                                indRefWorkOrder[4] = dim4_ref;

                                for ( dim3=0; dim3<dimSize[3]; dim3++ )
                                {
                                    size_t dim3_ref = dim3;
                                    if ( dim3 < dimSizeRef[3] )
                                    {
                                        this->ismrmrdDimIndex9D(indRef, dims[3], dim3);
                                    }
                                    else
                                    {
                                        this->ismrmrdDimIndex9D(indRef, dims[3], dimSizeRef[3]-1);
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
                                                this->ismrmrdDimIndex9D(indRef, dims[2], dim2);
                                            }
                                            else
                                            {
                                                this->ismrmrdDimIndex9D(indRef, dims[2], dimSizeRef[2]-1);
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

                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder_->data_, "workOrder_data");
                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder_->ref_, "workOrder_ref");

                        // trigger the recon
                        GADGET_CHECK_RETURN_FALSE(worker_->performRecon(workOrder_));

                        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, workOrder_->complexIm_, "workOrder_complexIm");

                        if ( shareAcrossWorkOrders )
                        {
                            workOrder_->workFlow_use_BufferedKernel_ = workFlow_use_BufferedKernel_;
                        }

                        // copy the results
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

                        long long offset = res_.calculate_offset(indRes);
                        memcpy(res_.begin()+offset, workOrder_->complexIm_.begin(), workOrder_->complexIm_.get_number_of_bytes());

                        if ( gfactor_needed && (workOrder_->gfactor_.get_size(0)==res_.get_size(0)) && (workOrder_->gfactor_.get_size(1) == res_.get_size(1)) )
                        {
                            size_t RO = gfactor_.get_size(0);
                            size_t E1 = gfactor_.get_size(1);
                            size_t N = gfactor_.get_size(3);
                            size_t S = gfactor_.get_size(4);

                            size_t gfactor_N = workOrder_->gfactor_.get_size(2);
                            size_t gfactor_S = workOrder_->gfactor_.get_size(3);

                            if ( (gfactor_N == N) && (gfactor_S == S) )
                            {
                                offset = gfactor_.calculate_offset(indRes);
                                memcpy(gfactor_.begin()+offset, workOrder_->gfactor_.begin(), workOrder_->gfactor_.get_number_of_bytes());
                            }
                            else
                            {
                                std::vector<size_t> indGfactor(8);
                                indGfactor[0] = 0;
                                indGfactor[1] = 0;
                                indGfactor[2] = 0;
                                indGfactor[3] = 0;
                                indGfactor[4] = dim5;
                                indGfactor[5] = dim6;
                                indGfactor[6] = dim7;
                                indGfactor[7] = dim8;

                                size_t n, s;
                                for ( s=0; s<S; s++ )
                                {
                                    for ( n=0; n<N; n++ )
                                    {
                                        indRes[3] = n;
                                        indRes[4] = s;
                                        offset = gfactor_.calculate_offset(indRes);

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
                        }

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

        GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "res_afterunwrapping");

        if ( gfactor_needed )
        {
            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, gfactor_, "gfactor_afterunwrapping");
        }

        // permute the res_ to the correct dimension order
        if (   ( (res_.get_number_of_elements()>dimResSize[0]*dimResSize[1]) && (dims[2]!=DIM_Channel) ) 
            || ( (res_.get_number_of_elements()>dimResSize[0]*dimResSize[1]*dimResSize[2])             ) )
        {
            std::vector<size_t> order;
            GADGET_CHECK_RETURN_FALSE(this->findISMRMRDPermuteOrder(dimsRes, dimsRes_, order));

            boost::shared_ptr< hoNDArray<T> > res_permuted = Gadgetron::permute(&res_, &order);
            res_.reshape(res_permuted->get_dimensions());
            memcpy(res_.begin(), res_permuted->begin(), res_permuted->get_number_of_bytes());

            GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, res_, "res_afterPermute");

            if ( gfactor_needed )
            {
                boost::shared_ptr< hoNDArray<T> > gfactor_permuted = Gadgetron::permute(&gfactor_, &order);
                gfactor_.reshape(gfactor_permuted->get_dimensions());
                memcpy(gfactor_.begin(), gfactor_permuted->begin(), gfactor_permuted->get_number_of_bytes());

                GADGET_EXPORT_ARRAY_COMPLEX(debugFolder_, gt_exporter_, gfactor_, "gfactor_afterPermute");
            }
        }
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in gtPlusISMRMRDReconWorkFlowCartesian<T>::configureWorkOrder(const std::vector<ISMRMRDDIM>& dims) ... ");
        return false;
    }

    return true;
}

}}
