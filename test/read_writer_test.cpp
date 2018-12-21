
#include <gtest/gtest.h>
#include "Message.h"
#include "Channel.h"
#include "writers/GadgetIsmrmrdWriter.h"
#include "readers/GadgetIsmrmrdReader.h"
#include <sstream>

TEST(ReadWriteTest,AcquisitionTest){
    using namespace Gadgetron;

    auto stream = std::stringstream{};

    auto acquisition_header = std::make_unique<ISMRMRD::AcquisitionHeader>();
    acquisition_header->number_of_samples = 32;
    acquisition_header->active_channels = 1;
    acquisition_header->available_channels = 1;

    auto data= std::make_unique<hoNDArray<std::complex<float>>>(32);
    data->fill(42.0f);


    auto message = std::make_unique<Core::MessageTuple>(std::move(acquisition_header),std::move(data));


    auto reader = GadgetIsmrmrdAcquisitionMessageReader();
    auto writer = GadgetIsmrmrdAcquisitionMessageWriter();


    ASSERT_TRUE(writer.accepts(*message));









}
