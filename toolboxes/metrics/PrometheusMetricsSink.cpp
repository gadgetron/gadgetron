#include "PrometheusMetricsSink.h"
#include <iostream>

using namespace prometheus;

namespace Gadgetron
{
    PrometheusMetricsSink::PrometheusMetricsSink()
        : exposer_{"127.0.0.1:8080"}
        , registry_(std::make_shared<Registry>())
        , startCounter_(0)
        , finishCounter_(0)
    {

        auto& startFamily = BuildCounter()
                                .Name("resconstructions_started_total")
                                .Help("Total number of reconstructions started")
                                .Register(*registry_);

        auto& start = startFamily.Add({});
        startCounter_ = &start;


        auto& finishFamily = BuildCounter()
                                .Name("resconstructions_finished_total")
                                .Help("Total number of reconstructions finished")
                                .Register(*registry_);

        auto& finish = finishFamily.Add({});
        finishCounter_ = &finish;

        exposer_.RegisterCollectable(registry_);

    }

    void PrometheusMetricsSink::ReconStart()
    {
        if (startCounter_ != 0)
        {
            startCounter_->Increment();
        }
    }

    void PrometheusMetricsSink::ReconFinish()
    {
        if (finishCounter_ != 0)
        {
            finishCounter_->Increment();
        }
    }
}