#include "ismrmrd/xml.h"
#include <gtest/gtest.h>

#include <stdlib.h>

#include "IsmrmrdContextVariables.h"

using namespace Gadgetron;

const std::string ismrmrd_header_with_variables =
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

const std::string ismrmrd_header_with_structured_measurement_id = R"(
  <?xml version="1.0"?>
  <ismrmrdHeader xmlns="http://www.ismrm.org/ISMRMRD" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xs="http://www.w3.org/2001/XMLSchema" xsi:schemaLocation="http://www.ismrm.org/ISMRMRD ismrmrd.xsd">
    <studyInformation>
      <studyTime>16:13:23</studyTime>
    </studyInformation>
    <measurementInformation>
      <measurementID>45387_925046864_1076436037_393</measurementID>
      <patientPosition>HFS</patientPosition>
      <protocolName>R1_PEFOV100_PERes100</protocolName>
      <measurementDependency>
        <dependencyType>SenMap</dependencyType>
        <measurementID>45387_925046864_1076436037_392</measurementID>
      </measurementDependency>
      <measurementDependency>
        <dependencyType>Noise</dependencyType>
        <measurementID>45387_925046864_1076436037_392</measurementID>
      </measurementDependency>
      <frameOfReferenceUID>1.3.12.2.1107.5.2.19.45387.1.20160311172423043.0.0.0</frameOfReferenceUID>
    </measurementInformation>
    <acquisitionSystemInformation>
      <systemVendor>SIEMENS</systemVendor>
      <systemModel>Skyra</systemModel>
      <systemFieldStrength_T>2.893620</systemFieldStrength_T>
      <relativeReceiverNoiseBandwidth>0.793000</relativeReceiverNoiseBandwidth>
      <receiverChannels>26</receiverChannels>
      <coilLabel>
        <coilNumber>58</coilNumber>
        <coilName>Body_18:1:B26</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>57</coilNumber>
        <coilName>Body_18:1:B25</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>26</coilNumber>
        <coilName>Spine_32:1:S21</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>25</coilNumber>
        <coilName>Spine_32:1:S22</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>17</coilNumber>
        <coilName>Body_18:1:B22</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>18</coilNumber>
        <coilName>Body_18:1:B21</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>49</coilNumber>
        <coilName>Spine_32:1:S34</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>50</coilNumber>
        <coilName>Spine_32:1:S33</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>54</coilNumber>
        <coilName>Body_18:1:B14</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>53</coilNumber>
        <coilName>Body_18:1:B13</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>30</coilNumber>
        <coilName>Body_18:1:B16</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>29</coilNumber>
        <coilName>Body_18:1:B15</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>34</coilNumber>
        <coilName>Body_18:1:B35</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>33</coilNumber>
        <coilName>Body_18:1:B36</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>22</coilNumber>
        <coilName>Body_18:1:B24</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>21</coilNumber>
        <coilName>Body_18:1:B23</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>42</coilNumber>
        <coilName>Body_18:1:B33</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>41</coilNumber>
        <coilName>Body_18:1:B34</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>2</coilNumber>
        <coilName>Spine_32:1:S23</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>1</coilNumber>
        <coilName>Spine_32:1:S24</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>61</coilNumber>
        <coilName>Body_18:1:B12</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>62</coilNumber>
        <coilName>Body_18:1:B11</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>45</coilNumber>
        <coilName>Spine_32:1:S32</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>46</coilNumber>
        <coilName>Spine_32:1:S31</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>37</coilNumber>
        <coilName>Body_18:1:B31</coilName>
      </coilLabel>
      <coilLabel>
        <coilNumber>38</coilNumber>
        <coilName>Body_18:1:B32</coilName>
      </coilLabel>
      <institutionName>NIH</institutionName>
    </acquisitionSystemInformation>
    <experimentalConditions>
      <H1resonanceFrequency_Hz>123159171</H1resonanceFrequency_Hz>
    </experimentalConditions>
    <encoding>
      <encodedSpace>
        <matrixSize>
          <x>384</x>
          <y>192</y>
          <z>1</z>
        </matrixSize>
        <fieldOfView_mm>
          <x>880.000000</x>
          <y>440.000000</y>
          <z>8.000000</z>
        </fieldOfView_mm>
      </encodedSpace>
      <reconSpace>
        <matrixSize>
          <x>192</x>
          <y>192</y>
          <z>1</z>
        </matrixSize>
        <fieldOfView_mm>
          <x>440.000000</x>
          <y>440.000000</y>
          <z>8.000000</z>
        </fieldOfView_mm>
      </reconSpace>
      <encodingLimits>
        <kspace_encoding_step_1>
          <minimum>0</minimum>
          <maximum>191</maximum>
          <center>96</center>
        </kspace_encoding_step_1>
        <kspace_encoding_step_2>
          <minimum>0</minimum>
          <maximum>0</maximum>
          <center>0</center>
        </kspace_encoding_step_2>
        <average>
          <minimum>0</minimum>
          <maximum>0</maximum>
          <center>0</center>
        </average>
        <slice>
          <minimum>0</minimum>
          <maximum>0</maximum>
          <center>0</center>
        </slice>
        <contrast>
          <minimum>0</minimum>
          <maximum>0</maximum>
          <center>0</center>
        </contrast>
        <phase>
          <minimum>0</minimum>
          <maximum>0</maximum>
          <center>0</center>
        </phase>
        <repetition>
          <minimum>0</minimum>
          <maximum>39</maximum>
          <center>0</center>
        </repetition>
        <set>
          <minimum>0</minimum>
          <maximum>0</maximum>
          <center>0</center>
        </set>
        <segment>
          <minimum>0</minimum>
          <maximum>0</maximum>
          <center>0</center>
        </segment>
      </encodingLimits>
      <trajectory>cartesian</trajectory>
      <parallelImaging>
        <accelerationFactor>
          <kspace_encoding_step_1>1</kspace_encoding_step_1>
          <kspace_encoding_step_2>1</kspace_encoding_step_2>
        </accelerationFactor>
        <calibrationMode>other</calibrationMode>
      </parallelImaging>
    </encoding>
    <sequenceParameters>
      <TR>570.239990</TR>
      <TE>1.320000</TE>
      <TI>300.000000</TI>
      <flipAngle_deg>12.000000</flipAngle_deg>
      <sequence_type>Flash</sequence_type>
      <echo_spacing>2.970000</echo_spacing>
    </sequenceParameters>
  </ismrmrdHeader>
)";

TEST(IsmrmrdContextVariablesTest, initialize_from_ismrmrd_header_with_all_elements) {

  ISMRMRD::IsmrmrdHeader head;
  ISMRMRD::deserialize(ismrmrd_header_with_variables.c_str(), head);

  IsmrmrdContextVariables ctx_vars(head);

  EXPECT_EQ(ctx_vars.get_device_id(), "20434");
  EXPECT_EQ(ctx_vars.get_subject_id(), "925046864");
  EXPECT_EQ(ctx_vars.get_session_id(), "1076436037");
  EXPECT_EQ(ctx_vars.get_measurement_id(), "LLJGHH888986");
}

TEST(IsmrmrdContextVariablesTest, initialize_from_ismrmrd_header_with_structured_measurement_id) {

  ISMRMRD::IsmrmrdHeader head;
  ISMRMRD::deserialize(ismrmrd_header_with_structured_measurement_id.c_str(), head);

  IsmrmrdContextVariables ctx_vars(head);
  
  EXPECT_EQ(ctx_vars.get_device_id(), "45387");
  EXPECT_EQ(ctx_vars.get_subject_id(), "925046864");
  EXPECT_EQ(ctx_vars.get_session_id(), "1076436037");
  EXPECT_EQ(ctx_vars.get_measurement_id(), "45387_925046864_1076436037_393");
}

TEST(IsmrmrdContextVariablesTest, malformed_structured_ids_handled) {
  auto measurements = std::vector<std::string> {
    "45387_925046864_925046864_1076436037_393",
    "453879250468649250468641076436037393",
    "45387__925046864_1076436037_393",
    };

  for (auto measurement_id : measurements) {
    ISMRMRD::IsmrmrdHeader head;
    ISMRMRD::deserialize(ismrmrd_header_with_structured_measurement_id.c_str(), head);
    head.measurementInformation.get().measurementID = measurement_id;

    IsmrmrdContextVariables ctx_vars(head);
    EXPECT_EQ(ctx_vars.get_subject_id(), "");
    EXPECT_EQ(ctx_vars.get_device_id(), "");
    EXPECT_EQ(ctx_vars.get_session_id(), "");
    EXPECT_EQ(ctx_vars.get_measurement_id(), measurement_id);
  }
}
