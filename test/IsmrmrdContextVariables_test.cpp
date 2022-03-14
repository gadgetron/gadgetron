#include "ismrmrd/xml.h"
#include <gtest/gtest.h>

#include <stdlib.h>

#include "IsmrmrdContextVariables.h"

using namespace Gadgetron;

const std::string ismrmrd_head_example =
    "<?xml version=\"1.0\"?>"
    "<ismrmrdHeader xmlns=\"http://www.ismrm.org/ISMRMRD\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema\" xsi:schemaLocation=\"http://www.ismrm.org/ISMRMRD ismrmrd.xsd\">"
    "  <subjectInformation>"
    "    <patientID>925046864</patientID>"
    "  </subjectInformation>"
    "  <studyInformation>"
    "    <studyID>1076436037</studyID>"
    "  </studyInformation>"
    "  <measurementInformation>"
    "    <measurementID>LLJGHH888986</measurementID>"
    "    <seriesDate>2012-08-13</seriesDate>"
    "    <seriesTime>09:10:12</seriesTime>"
    "    <patientPosition>HFS</patientPosition>"
    "    <initialSeriesNumber>1</initialSeriesNumber>"
    "    <protocolName>ExampleProt</protocolName>"
    "    <seriesDescription>MRIStudy1</seriesDescription>"
    "    <measurementDependency>"
    "      <dependencyType>Noise</dependencyType>"
    "      <measurementID>HHJJHJL000977889</measurementID>"
    "    </measurementDependency>"
    "    <measurementDependency>"
    "      <dependencyType>SurfaceCoilCorrection</dependencyType>"
    "      <measurementID>HHJJHJL000977810</measurementID>"
    "    </measurementDependency>"
    "  </measurementInformation>"
    "  <acquisitionSystemInformation>"
    "    <systemFieldStrength_T>1.494</systemFieldStrength_T>"
    "    <relativeReceiverNoiseBandwidth>0.79</relativeReceiverNoiseBandwidth>"
    "    <deviceID>20434</deviceID>"
    "  </acquisitionSystemInformation>"
    "  <experimentalConditions>"
    "    <H1resonanceFrequency_Hz>63642459</H1resonanceFrequency_Hz>"
    "  </experimentalConditions>"
    "  <encoding>"
    "    <encodedSpace>"
    "      <matrixSize>"
    "        <x>256</x>"
    "        <y>140</y>"
    "        <z>80</z>"
    "      </matrixSize>"
    "      <fieldOfView_mm>"
    "        <x>600</x>"
    "        <y>328.153125</y>"
    "        <z>160</z>"
    "      </fieldOfView_mm>"
    "    </encodedSpace>"
    "    <reconSpace>"
    "      <matrixSize>"
    "        <x>128</x>"
    "        <y>116</y>"
    "        <z>64</z>"
    "      </matrixSize>"
    "      <fieldOfView_mm>"
    "        <x>300</x>"
    "        <y>271.875</y>"
    "        <z>128</z>"
    "      </fieldOfView_mm>"
    "    </reconSpace>"
    "    <encodingLimits>"
    "      <kspace_encoding_step_1>"
    "        <minimum>0</minimum>"
    "        <maximum>83</maximum>"
    "        <center>28</center>"
    "      </kspace_encoding_step_1>"
    "    </encodingLimits>"
    "    <trajectory>cartesian</trajectory>"
    "  </encoding>"
    "</ismrmrdHeader>";

TEST(IsmrmrdContextVariablesTest, initialize_from_ismrmrd_header) {

  ISMRMRD::IsmrmrdHeader head;
  ISMRMRD::deserialize(ismrmrd_head_example.c_str(), head);

  IsmrmrdContextVariables ctx_vars(head);

  EXPECT_EQ(ctx_vars.get_device_id(), "20434");
  EXPECT_EQ(ctx_vars.get_subject_id(), "925046864");
  EXPECT_EQ(ctx_vars.get_session_id(), "1076436037");
  EXPECT_EQ(ctx_vars.get_measurement_id(), "LLJGHH888986");
}

TEST(IsmrmrdContextVariablesTest, initialize_from_structured_id) {
  std::string const structured_id("45387_925046864_1076436037_393");

  IsmrmrdContextVariables ctx_vars(structured_id);
  EXPECT_EQ(ctx_vars.get_device_id(), "45387");
  EXPECT_EQ(ctx_vars.get_subject_id(), "925046864");
  EXPECT_EQ(ctx_vars.get_session_id(), "1076436037");
  EXPECT_EQ(ctx_vars.get_measurement_id(), "393");
}

TEST(IsmrmrdContextVariablesTest, malformed_structured_ids_handled) {
  IsmrmrdContextVariables ctx_vars("45387_925046864_925046864_1076436037_393");
  IsmrmrdContextVariables ctx_vars("453879250468649250468641076436037393")
  IsmrmrdContextVariables ctx_vars("45387__925046864_1076436037_393")
}
