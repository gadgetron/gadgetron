#include "gtest/gtest.h"

#include "node_discovery.h"
#include <vector>
#include <fstream>

TEST (NodeDiscoveryTest, ReadNodeDiscoveryList) {

  //Create a file with node discovery data
  
  std::ofstream f("node_list.json");
  f << "{" << std::endl;
  f << "  \"version\": \"1.0.0\"," << std::endl;
  f << "  \"nodes\": [" << std::endl;
  f << "	{" << std::endl;
  f << "	    \"id\": \"1993a160-6227-46be-872f-8719ee703eec\"," << std::endl;
  f << "	    \"host\": \"10.0.0.1\"," << std::endl;
  f << "	    \"port\": 9002," << std::endl;
  f << "	    \"restPort\": 9080," << std::endl;
  f << "	    \"activeReconstructions\": 0," << std::endl;
  f << "	    \"hasGpu\": false" << std::endl;
  f << "	}," << std::endl;
  f << "	{" << std::endl;
  f << "	    \"id\": \"794d8d45-e7f0-4c96-880e-8a6f09a18ce2\"," << std::endl;
  f << "	    \"host\": \"10.0.0.2\"," << std::endl;
  f << "	    \"port\": 9002," << std::endl;
  f << "	    \"restPport\": 9080," << std::endl;
  f << "	    \"activeReconstructions\": 0," << std::endl;
  f << "	    \"hasGpu\": false" << std::endl;
  f << "	}" << std::endl;
  f << "    ]" << std::endl;
  f << "}" << std::endl;
  
  
  std::vector<Gadgetron::NodeDiscoveryEntry> nodes;

#ifdef _WIN32
  bool ret = DiscoverNodes(nodes, "type node_list.json");
#else
  bool ret = DiscoverNodes(nodes, "cat node_list.json");
#endif
  
  ASSERT_EQ(true, ret);
  ASSERT_EQ (2, nodes.size());
  ASSERT_EQ ("794d8d45-e7f0-4c96-880e-8a6f09a18ce2", nodes[1].id);
  ASSERT_EQ ("10.0.0.1", nodes[0].host);
}
