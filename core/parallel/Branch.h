#pragma once

#include <map>
#include <memory>
#include <boost/dll.hpp>

#include "Channel.h"

namespace Gadgetron::Core::Parallel {

    class Branch {
    public:
        virtual ~Branch() = default;
        virtual void process(std::shared_ptr<Channel>, std::map<std::string, std::shared_ptr<Channel>>) = 0;
    };
}


#define GADGETRON_BRANCH_EXPORT(BranchClass)                                    \
std::unique_ptr<Gadgetron::Core::Parallel::Branch>                              \
branch_factory_##BranchClass() {                                                \
    return std::make_unique<BranchClass>();                                     \
}                                                                               \
                                                                                \
BOOST_DLL_ALIAS(                                                                \
        branch_factory_##BranchClass,                                           \
        branch_factory_export_##BranchClass                                     \
)                                                                               \
