#ifndef PROMETHEUSMETRICSSINK_H
#define PROMETHEUSMETRICSSINK_HPP

#include "metrics_export.h"
#include "MetricsSink.h"

#include <prometheus/counter.h>
#include <prometheus/exposer.h>
#include <prometheus/registry.h>

using namespace prometheus;

namespace Gadgetron
{
    class EXPORTGADGETRONMETRICS PrometheusMetricsSink : public MetricsSink
    {
    public:
        PrometheusMetricsSink();
        virtual void ReconStart();
        virtual void ReconFinish();

    protected:
        Exposer exposer_;
        std::shared_ptr<Registry> registry_;
        Counter* startCounter_;
        Counter* finishCounter_;
    };
}

#endif // PROMETHEUSMETRICSSINK_HPP