#pragma once

#include <string>

#include <ismrmrd/xml.h>

namespace Gadgetron {

/**
 * Used to extract subject, device, session, and measurement IDs from an ISMRMRD
 * header, falling back to extracting values from a the measurement ID if
 * it is formatted as DEVICE_SUBJECT_SESSION_MEASUREMENT.
 */
class IsmrmrdContextVariables {
  public:
  IsmrmrdContextVariables(ISMRMRD::IsmrmrdHeader const& head);

  std::string get_subject_id() {
    return subject_id;
  }

  std::string get_device_id() {
    return device_id;
  }

  std::string get_session_id() {
    return session_id;
  }

  std::string get_measurement_id() {
    return measurement_id;
  }

  private:
  std::string subject_id = "";
  std::string device_id = "";
  std::string session_id = "";
  std::string measurement_id = "";
};

} // namespace Gadgetron
