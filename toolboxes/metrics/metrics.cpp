#include <metrics.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include <time.h>
#include <cstring>
#include <chrono>
#include <memory>
#include <thread>

using namespace prometheus;

namespace Gadgetron
{
  Metrics* Metrics::instance()
  {
    if (!instance_) instance_ = new Metrics();
    return instance_;
  }

  Metrics* Metrics::instance_ = NULL;

  Metrics::Metrics()
    : exposer_{"127.0.0.1:8080"}
  {
    auto registry = std::make_shared<prometheus::Registry>();

    // add a new counter family to the registry (families combine values with the
    // same name, but distinct label dimensions)
    auto& counter_family = BuildCounter()
                              .Name("time_running_seconds_total")
                              .Help("How many seconds is this server running?")
                              .Labels({{"label", "value"}})
                              .Register(*registry);

      // add a counter to the metric family
      auto& second_counter = counter_family.Add(
          {{"another_label", "value"}, {"yet_another_label", "value"}});

      // ask the exposer to scrape the registry on incoming scrapes
    exposer_.RegisterCollectable(registry);

  }
}