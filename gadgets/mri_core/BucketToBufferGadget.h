#pragma once
#include "Node.h"
#include "hoNDArray.h"
#include "gadgetron_mricore_export.h"

#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>
#include <complex>
#include <map>
#include "mri_core_data.h"
#include "mri_core_acquisition_bucket.h"

namespace Gadgetron{

    // TODO the ignore_segment_ flag is a hack for some EPI sequences
    // should be fixed on the converter side.

    // This gadget fills the IsmrmrdReconData structures with kspace readouts and sets up the sampling limits
    // For the cartesian sampling, the filled kspace ensures its center (N/2) is aligned with the specified center in the encoding limits
    // For the non-cartesian sampling, this "center alignment" constraint is not applied and kspace lines are filled as their E1 and E2 indexes

    // Since the order of data can be changed from its acquried time order, there is no easy way to resort waveform data
    // Therefore, the waveform data was copied and passed with every buffer

  class BucketToBufferGadget :
public Core::ChannelGadget<AcquisitionBucket>
    {
    public:
        BucketToBufferGadget(const Core::Context& context, const Core::GadgetProperties& props);
        enum class Dimension {
            average,contrast,phase,repetition,set,segment,slice,none;
        };

    protected:
      NODE_PROPERTY(N_dimension,Dimension , "N-Dimensions", "" );
      NODE_PROPERTY(S_dimension,Dimension , "S-Dimensions", "" );

      NODE_PROPERTY(split_slices, bool, "Split slices", false);
      NODE_PROPERTY(ignore_segment, bool, "Ignore segment", false);
      NODE_PROPERTY(verbose, bool, "Whether to print more information", false);

      bool split_slices_;
      bool ignore_segment_;
      ISMRMRD::IsmrmrdHeader header;
      
      void process(Core::InputChannel<AcquisitionBucket>& in, Core::OutputChannel& out) override;
      size_t getKey(ISMRMRD::ISMRMRD_EncodingCounters idx) const;
      size_t getSlice(ISMRMRD::ISMRMRD_EncodingCounters idx) const ;
      size_t getN(ISMRMRD::ISMRMRD_EncodingCounters idx) const;
      size_t getS(ISMRMRD::ISMRMRD_EncodingCounters idx)const ;

      IsmrmrdReconBit & getRBit(std::map<size_t, IsmrmrdReconData > & recon_data_buffers, size_t key, uint16_t espace) const ;
      void allocateDataArrays(IsmrmrdDataBuffered &  dataBuffer, ISMRMRD::AcquisitionHeader & acqhdr, ISMRMRD::Encoding encoding, AcquisitionBucketStats& stats, bool forref);
      void fillSamplingDescription(SamplingDescription & sampling, ISMRMRD::Encoding & encoding,
          AcquisitionBucketStats& stats, ISMRMRD::AcquisitionHeader & acqhdr, bool forref);
      void stuff(std::vector<Core::Acquisition>::iterator it, IsmrmrdDataBuffered & dataBuffer, ISMRMRD::Encoding encoding, AcquisitionBucketStats& stats, bool forref);
    };

    void from_string(const std::string&, BucketToBufferGadget::Dimension& );
}
