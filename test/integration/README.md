# Integration Test Conversion

This directory is now considered "legacy", but the test cases remain here until we have a proper "Siemens -> MRDv2" converter.

For now, if the MRDv2 model changes, or Siemens test data is added/updated, we can use the `convert_test_data.py` tool in this directory.

The converter tool reads each of the legacy "integration" test cases (`./cases/`) and test data (`./data`) and generates new end-to-end tests in `test/e2e/cases/`.

Once we have a proper "Siemens -> MRDv2" converter, we can skip the intermediate ISMRMRD v1 conversion and remove this directory.