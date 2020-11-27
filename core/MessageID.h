#pragma once

#include <cstdint>

namespace Gadgetron::Core {
    enum MessageID : uint16_t {
        FILENAME                                           = 1,
        CONFIG                                             = 2,
        HEADER                                             = 3,
        CLOSE                                              = 4,
        TEXT                                               = 5,
        QUERY                                              = 6,
        RESPONSE                                           = 7,
        ERROR                                              = 8,
		GADGET_MESSAGE_EXT_ID_MIN                          = 1000,
		GADGET_MESSAGE_ISMRMRD_ACQUISITION                 = 1008,
		GADGET_MESSAGE_DICOM_WITHNAME                      = 1018,
		GADGET_MESSAGE_DEPENDENCY_QUERY                    = 1019,
		GADGET_MESSAGE_ISMRMRD_IMAGE                       = 1022,
		GADGET_MESSAGE_RECONDATA                           = 1023,
		GADGET_MESSAGE_ISMRMRD_IMAGE_ARRAY                 = 1024,
		GADGET_MESSAGE_ISMRMRD_WAVEFORM                    = 1026,
		GADGET_MESSAGE_BUCKET                              = 1050,
		GADGET_MESSAGE_BUNDLE                              = 1051
    };
}
