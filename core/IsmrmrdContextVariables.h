#pragma once

#include <string>

#include <ismrmrd/xml.h>

namespace Gadgetron {

/**
 * Used to extract subject, device, session, and measurement IDs from an ISMRMRD
 * header and/or a structured measurement ID string. 
 * Structured measurement IDs are expected to be formatted as 
 * DEVICE_SUBJECT_SESSION_MEASUREMENT
 */
class IsmrmrdContextVariables {
  public:

  /**
   * Reads IDs from the given header. For any missing values, it attempts
   */
  IsmrmrdContextVariables(ISMRMRD::IsmrmrdHeader const& head);
  IsmrmrdContextVariables(std::string const& structured_measurement_id);

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
