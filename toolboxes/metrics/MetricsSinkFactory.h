#ifndef METRICSSINKFACTORY_H
#define METRICSSINKFACTORY_H

#include <memory>
#include "MetricsSink.h"

#ifdef PROMETHEUS_METRICS
#include "PrometheusMetricsSink.h"
#else
#include "NoOpMetricsSink.h"
#endif

namespace Gadgetron
{
    std::shared_ptr<MetricsSink> CreateMetricsSink()
    {
#ifdef PROMETHEUS_METRICS
    return std::make_shared<PrometheusMetricsSink>();
#else
    return std::make_shared<NoOpMetricsSink>();
#endif
    }
}

#endif