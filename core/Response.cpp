
#include "Response.h"

Gadgetron::Core::Response::Response(uint64_t correlation_id, const std::string &response)
: correlation_id(correlation_id), response(response) {}
