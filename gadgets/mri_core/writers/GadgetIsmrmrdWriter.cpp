//
// Created by dchansen on 12/19/18.
//

#include <GadgetMRIHeaders.h>
#include "GadgetIsmrmrdWriter.h"
#include "io/writers.h"

void Gadgetron::GadgetIsmrmrdWaveformMessageWriter::serialize(std::ostream &stream,
                                                              std::unique_ptr<ISMRMRD::WaveformHeader> header,
                                                              std::unique_ptr<hoNDArray<uint32_t>> array) {

    using namespace Core;
    IO::write(stream,GADGET_MESSAGE_ISMRMRD_WAVEFORM);
    IO::write(stream,*header);
    IO::write(stream,*array);

}


void Gadgetron::GadgetIsmrmrdAcquisitionMessageWriter::serialize(std::ostream &stream,
                                                                 std::unique_ptr<ISMRMRD::AcquisitionHeader> header,
                                                                 std::unique_ptr<boost::optional<Gadgetron::hoNDArray<float>>> trajectory,
                                                                 std::unique_ptr<Gadgetron::hoNDArray<std::complex<float>>> data) {
    using namespace Core;

    IO::write(stream,GADGET_MESSAGE_ISMRMRD_ACQUISITION);
    IO::write(stream,*header);
    if (*trajectory)
        IO::write(stream,**trajectory);

    IO::write(stream,*data);
}
