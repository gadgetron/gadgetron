#pragma once
#include "Node.h"
#include "gadgetron_mricore_export.h"
#include "hoNDArray.h"

#include "mri_core_acquisition_bucket.h"
#include "mri_core_data.h"
#include <complex>
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

namespace Gadgetron {

    // TODO the ignore_segment_ flag is a hack for some EPI sequences
    // should be fixed on the converter side.

    // This gadget fills the IsmrmrdReconData structures with kspace readouts and sets up the sampling limits
    // For the cartesian sampling, the filled kspace ensures its center (N/2) is aligned with the specified center in
    // the encoding limits For the non-cartesian sampling, this "center alignment" constraint is not applied and kspace
    // lines are filled as their E1 and E2 indexes

    // Since the order of data can be changed from its acquried time order, there is no easy way to resort waveform data
    // Therefore, the waveform data was copied and passed with every buffer

    class BucketToBufferGadget : public Core::ChannelGadget<AcquisitionBucket> {
    public:
        BucketToBufferGadget(const Core::Context& context, const Core::GadgetProperties& props);
        enum class Dimension { average, contrast, phase, repetition, set, segment, slice, none };

    struct BufferKey {
        uint16_t average,slice,contrast,phase,repetition,set,segment;
        BufferKey(const BufferKey&) = default;
        BufferKey(const ISMRMRD::EncodingCounters& idx) : average{idx.average}, slice{idx.slice},contrast{idx.contrast}, phase{idx.phase},repetition{idx.repetition},set{idx.set},segment{idx.segment}{
            
        }
    };
    protected:
        NODE_PROPERTY(N_dimension, Dimension, "N-Dimensions", Dimension::none);
        NODE_PROPERTY(S_dimension, Dimension, "S-Dimensions", Dimension::none);

        NODE_PROPERTY(split_slices, bool, "Split slices", false);
        NODE_PROPERTY(ignore_segment, bool, "Ignore segment", false);
        NODE_PROPERTY(verbose, bool, "Whether to print more information", false);

        ISMRMRD::IsmrmrdHeader header;

        void process(Core::InputChannel<AcquisitionBucket>& in, Core::OutputChannel& out) override;
        BufferKey getKey(const ISMRMRD::EncodingCounters& idx) const;


        IsmrmrdDataBuffered makeDataBuffer(const ISMRMRD::AcquisitionHeader& acqhdr, ISMRMRD::Encoding encoding,
            const AcquisitionBucketStats& stats, bool forref) const;
        SamplingDescription createSamplingDescription(const ISMRMRD::Encoding& encoding,
            const AcquisitionBucketStats& stats, const ISMRMRD::AcquisitionHeader& acqhdr, bool forref) const ;
        void add_acquisition(IsmrmrdDataBuffered& dataBuffer, const Core::Acquisition& acq, ISMRMRD::Encoding encoding,
            const AcquisitionBucketStats& stats, bool forref);
        uint16_t getNE0(const ISMRMRD::AcquisitionHeader& acqhdr, const ISMRMRD::Encoding& encoding) const;
        uint16_t getNE1(const ISMRMRD::Encoding& encoding, const AcquisitionBucketStats& stats, bool forref) const;
        uint16_t getNE2(const ISMRMRD::Encoding& encoding, const AcquisitionBucketStats& stats, bool forref) const;
        uint16_t getNLOC(const ISMRMRD::Encoding& encoding, const AcquisitionBucketStats& stats) const;
    };

    void from_string(const std::string&, BucketToBufferGadget::Dimension&);
}
