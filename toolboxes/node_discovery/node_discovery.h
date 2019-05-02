#ifndef GADGETRON_NODEDISCOVERY_H
#define GADGETRON_NODEDISCOVERY_H

#include "node_discovery_export.h"
#include "log.h"

#include <string>
#include <vector>

#define NODEDISCOVERY_ENVIRONMENT_VARIABLE "GADGETRON_NODEDISCOVERY"

namespace Gadgetron
{
  
  struct NodeDiscoveryEntry
  {
    std::string id;
    std::string host;
    int port;
    int restPort;
    int activeReconstructions;
    bool hasGpu;
  };

  
  EXPORTGADGETRONNODEDISCOVERY bool IsNodeDiscoveryConfigured();
  EXPORTGADGETRONNODEDISCOVERY bool DiscoverNodes(std::vector<NodeDiscoveryEntry>& nodes);
  EXPORTGADGETRONNODEDISCOVERY bool DiscoverNodes(std::vector<NodeDiscoveryEntry>& nodes, const char* command);
  
}

#endif //GADGETRON_NODEDISCOVERY_H
