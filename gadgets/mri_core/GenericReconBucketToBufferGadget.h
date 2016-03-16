/** \file   GenericReconBucketToBufferGadget.h
    \brief  This is the class gadget for both 2DT and 3DT generic recon bucket to buffer conversion. 
            This conversion utlizes the encoding space matrix size and makes sure the alignment of kspace center.
    \author Hui Xue
*/

#pragma once

#include "BucketToBufferGadget.h"

namespace Gadgetron {

    class EXPORTGADGETSMRICORE GenericReconBucketToBufferGadget : public BucketToBufferGadget
    {
    public:
        GADGET_DECLARE(GenericReconBucketToBufferGadget);

        typedef BucketToBufferGadget BaseClass;

        GenericReconBucketToBufferGadget();
        virtual ~GenericReconBucketToBufferGadget();

    protected:

        virtual int process_config(ACE_Message_Block* mb);

        virtual void allocateDataArrays(IsmrmrdDataBuffered &  dataBuffer, ISMRMRD::AcquisitionHeader & acqhdr, ISMRMRD::Encoding encoding, IsmrmrdAcquisitionBucketStats & stats, bool forref);
        virtual void fillSamplingDescription(SamplingDescription & sampling, ISMRMRD::Encoding & encoding, IsmrmrdAcquisitionBucketStats & stats, ISMRMRD::AcquisitionHeader & acqhdr, bool forref);
        virtual void stuff(std::vector<IsmrmrdAcquisitionData>::iterator it, IsmrmrdDataBuffered & dataBuffer, ISMRMRD::Encoding encoding, bool forref);
    };
}
