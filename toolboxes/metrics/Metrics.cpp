#include "Metrics.h"
#include "MetricsSinkFactory.h"


namespace Gadgetron
{
  Metrics* Metrics::instance()
  {
    if (!instance_) instance_ = new Metrics();
    return instance_;
  }

  Metrics* Metrics::instance_ = NULL;

  Metrics::Metrics()
    : metricsSink_(CreateMetricsSink())
  {
  }

  void Metrics::ReconStart()
  {
    this->metricsSink_->ReconStart();
  }

  void Metrics::ReconFinish()
  {
    this->metricsSink_->ReconFinish();
  }
}