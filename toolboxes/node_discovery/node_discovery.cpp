#include "node_discovery.h"

#include <cstdlib>
#include <sstream>
#include <iostream>

#include "boost/process.hpp"
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

namespace pt = boost::property_tree;
namespace bp = boost::process;

namespace Gadgetron
{
  
  bool IsNodeDiscoveryConfigured()
  {
    return std::getenv(NODEDISCOVERY_ENVIRONMENT_VARIABLE);
  }

  bool DiscoverNodes(std::vector<NodeDiscoveryEntry>& nodes)
  {
    if (IsNodeDiscoveryConfigured()) {
      return DiscoverNodes(nodes, std::getenv(NODEDISCOVERY_ENVIRONMENT_VARIABLE));
    } else {
      GERROR("Unable to do node discovery. Environment variable %s not set\n");
      nodes.clear();
      return false;
    }
  }

  bool DiscoverNodes(std::vector<NodeDiscoveryEntry>& nodes, const char* command)
  {

    //Setup and run command
    bp::ipstream is; //reading pipe-stream
    bp::child c(command, bp::std_out > is);

    //Wait for completion
    c.wait();

    //Parse output
    pt::ptree tree;
    pt::read_json(is, tree);

    //Version is a requirement
    nodes.clear();
    if (tree.get<std::string>("version") != "1.0.0")
    {
      GERROR("Version not found in json node list\n");
      return false;
    }

    //Find as many nodes as we can (could be empty)
    for (pt::ptree::value_type &node : tree.get_child("nodes"))
    {
	
      NodeDiscoveryEntry entry;
      entry.id = node.second.get<std::string>("id");
      entry.host = node.second.get<std::string>("host");
      entry.port = node.second.get<int>("port", 9002);
      entry.restPort = node.second.get<int>("restPort", 9080);
      entry.hasGpu = node.second.get<bool>("hasGpu", false);
      entry.activeReconstructions = node.second.get<int>("activeReconstructions", 0);
      
      nodes.push_back(entry);
    }

    return true;
  }
  
} //Namespace Gadgetron
