#ifndef NOOPMETRICSSINK_H
#define NOOPMETRICSSINK_H

#include "metrics_export.h"
#include "MetricsSink.h"

namespace Gadgetron
{
    class EXPORTGADGETRONMETRICS NoOpMetricsSink : public MetricsSink
    {
    public:
        virtual void ReconStart()
        {
            // Do nothing
        }

        virtual void ReconFinish()
        {
            // Do nothing
        }
    };
}

#endif // NOOPMETRICSSINK_H