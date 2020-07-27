#ifndef METRICSSINK_H
#define METRICSSINK_H

#include "metrics_export.h"

namespace Gadgetron
{
    class EXPORTGADGETRONMETRICS MetricsSink
    {
    public:
        virtual void ReconStart() = 0;
        virtual void ReconFinish() = 0;
        virtual ~MetricsSink() {}
    };
}

#endif // METRICSSINK_H