
#include <gtest/gtest.h>
#include "Message.h"
#include "Channel.h"
#include "writers/GadgetIsmrmrdWriter.h"
#include "readers/GadgetIsmrmrdReader.h"
#include <sstream>
#include "hoNDArray_elemwise.h"
#include "io/primitives.h"

TEST(ReadWriteTest,AcquisitionTest){
    using namespace Gadgetron;

    auto stream = std::stringstream{};

    auto acquisition_header = ISMRMRD::AcquisitionHeader();
    acquisition_header.number_of_samples = 32;
    acquisition_header.active_channels = 1;
    acquisition_header.available_channels = 1;

    auto data= hoNDArray<std::complex<float>>(32,1);
    data.fill(42.0f);


    auto message = Core::Message(acquisition_header,data);


    auto reader = GadgetIsmrmrdAcquisitionMessageReader();
    auto writer = GadgetIsmrmrdAcquisitionMessageWriter();


    ASSERT_TRUE(writer.accepts(message));

    writer.write(stream,std::move(message));


   auto message_id = Core::IO::read<uint16_t>(stream);

   ASSERT_EQ(message_id,reader.slot());

    auto unpacked = Core::unpack<Core::Acquisition>(reader.read(stream));

    EXPECT_TRUE(bool(unpacked));

    auto data_unpacked = std::get<hoNDArray<std::complex<float>>>(*unpacked);
    ASSERT_EQ(data,data_unpacked);
}
