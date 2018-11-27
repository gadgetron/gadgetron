#pragma once

#include "AST.h"
#include "Stream.h"

namespace Gadgetron::Core {


    Stream build_stream(const AST::Stream& stream_ast );
}