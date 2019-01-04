
#include "Response.h"

Gadgetron::Core::Response::Response(uint16_t correlation_id, const std::string &response)
: correlation_id(correlation_id), response(response) {}
