#include "IsmrmrdContextVariables.h"

#include <optional>
#include <regex>

#include "log.h"

namespace Gadgetron {

struct ExtractedContextVariables {
  std::string SubjectId;
  std::string DeviceId;
  std::string SessionId;
};

std::optional<ExtractedContextVariables> ExtractContextVariables(std::string const& structured_id) {
  std::regex reg("^([^_]*)_([^_]*)_([^_]*)_([^_]*)$");
  std::smatch match;

  if (std::regex_match(structured_id, match, reg)) {
    return ExtractedContextVariables{match[2], match[1], match[3]};
  } else {
    return std::nullopt;
  }
}

IsmrmrdContextVariables::IsmrmrdContextVariables(ISMRMRD::IsmrmrdHeader const& head) {
  if (head.acquisitionSystemInformation.is_present()) {
    auto acq_system = head.acquisitionSystemInformation.get();
    if (acq_system.deviceID.has_value()) {
      device_id = acq_system.deviceID.get();
    }
  }

  if (head.measurementInformation.is_present()) {
    auto meas_info = head.measurementInformation.get();
    if (meas_info.measurementID.is_present()) {
      measurement_id = meas_info.measurementID.get();
    }
  }

  if (head.subjectInformation.is_present()) {
    auto subject_info = head.subjectInformation.get();
    if (subject_info.patientID.is_present()) {
      subject_id = subject_info.patientID.get();
    }
  }

  if (head.studyInformation.is_present()) {
    auto study_info = head.studyInformation.get();
    if (study_info.studyID.is_present()) {
      session_id = study_info.studyID.get();
    }
  }

  if (measurement_id.empty()) {
    throw std::runtime_error("Empty measurement ID not allowed");
  }

  if (session_id.empty() || device_id.empty() || subject_id.empty()) {
    // We were unable to find all that we need, we can attempt to extract from measurement id.
    auto extracted = ExtractContextVariables(measurement_id);
    if (!extracted)  {
      GWARN_STREAM("WARNING: Attempted to extract context variables from measurement id, but failed.");
      return;
    }

    if (device_id.empty()) {
      device_id = extracted->DeviceId;
    }

    if (subject_id.empty()) {
      subject_id = extracted->SubjectId;
    }

    if (session_id.empty()) {
      session_id = extracted->SessionId;
    }
  }
}

IsmrmrdContextVariables::IsmrmrdContextVariables(std::string const& structured_measurement_id)
  : measurement_id(structured_measurement_id) {
  if (auto extracted = ExtractContextVariables(structured_measurement_id); extracted) {
    device_id = extracted->DeviceId;
    subject_id = extracted->SubjectId;
    session_id = extracted->SessionId;
  }
}

} // namespace Gadgetron
