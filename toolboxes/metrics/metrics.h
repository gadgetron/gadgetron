/** \file metrics.h
    \brief Singleton for Prometheus metrics
*/

#ifndef METRICS_H
#define METRICS_H

#include "metrics_export.h"

#ifdef PROMETHEUS_METRICS
#include <prometheus/counter.h>
#include <prometheus/exposer.h>
#include <prometheus/registry.h>
#endif

namespace Gadgetron
{
    class EXPORTGADGETRONMETRICS Metrics
    {
    public:
        static Metrics* instance();

    protected:
        Metrics();
        static Metrics* instance_;
#ifdef PROMETHEUS_METRICS
        prometheus::Exposer exposer_;
#endif
    };
}


#endif // METRICS_H