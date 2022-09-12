#pragma once

#include "writers/AcquisitionWriter.h"
#include "writers/WaveformWriter.h"

namespace Gadgetron {
    using GadgetIsmrmrdWaveformMessageWriter = Core::Writers::WaveformWriter;
    using GadgetIsmrmrdAcquisitionMessageWriter = Core::Writers::AcquisitionWriter;
}