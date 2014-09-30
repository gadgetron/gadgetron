#include "gadgetron_config.h"

#include <iostream>

int main(int argc, char** argv)
{
  std::cout << "Gadgetron Version Info" << std::endl;
  std::cout << "  -- Version  : " << GADGETRON_VERSION_STRING << std::endl;
  std::cout << "  -- Git SHA1 : " << GADGETRON_GIT_SHA1_HASH << std::endl;
  return 0;
}
