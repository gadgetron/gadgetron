/** \file metrics.h
    \brief Singleton for Prometheus metrics
*/

#ifndef METRICS_H
#define METRICS_H

#include <memory>
#include "metrics_export.h"
#include "MetricsSink.h"

namespace Gadgetron
{
    class EXPORTGADGETRONMETRICS Metrics
    {
    public:
        static Metrics* instance();
        void ReconStart();
        void ReconFinish();

    protected:
        Metrics();
        static Metrics* instance_;
        std::shared_ptr<MetricsSink> metricsSink_;
    };
}


#endif // METRICS_H