//
// Created by dchansen on 12/19/18.
//

#include <GadgetMRIHeaders.h>
#include "GadgetIsmrmrdWriter.h"
#include "io/primitives.h"

void Gadgetron::GadgetIsmrmrdWaveformMessageWriter::serialize(std::ostream &stream,
                                                              const ISMRMRD::WaveformHeader& header,
                                                              const hoNDArray<uint32_t>& array) {

    using namespace Core;
    IO::write(stream,GADGET_MESSAGE_ISMRMRD_WAVEFORM);
    IO::write(stream,header);
    IO::write(stream,array);

}


void Gadgetron::GadgetIsmrmrdAcquisitionMessageWriter::serialize(std::ostream &stream,
                                                                 const ISMRMRD::AcquisitionHeader& header,
                                                                 const boost::optional<Gadgetron::hoNDArray<float>>& trajectory,
                                                                 const Gadgetron::hoNDArray<std::complex<float>>& data) {
    using namespace Core;

    IO::write(stream,GADGET_MESSAGE_ISMRMRD_ACQUISITION);
    IO::write(stream,header);
    if (trajectory)
        IO::write(stream,*trajectory);

    IO::write(stream,data);
}
