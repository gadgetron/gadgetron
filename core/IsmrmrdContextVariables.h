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

    IsmrmrdContextVariables(std::string const& subject_id, std::string const& device_id, std::string const& session_id,
                            std::string const& measurement_id)
        : subject_id_(subject_id), device_id_(device_id), session_id_(session_id), measurement_id_(measurement_id){};

    std::string const& subject_id() const { return subject_id_; }

    std::string const& device_id() const { return device_id_; }

    std::string const& session_id() const { return session_id_; }

    std::string const& measurement_id() const { return measurement_id_; }

  private:
    std::string subject_id_ = "";
    std::string device_id_ = "";
    std::string session_id_ = "";
    std::string measurement_id_ = "";
};

} // namespace Gadgetron
