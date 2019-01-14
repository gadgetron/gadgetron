#pragma once

#include <map>
#include <memory>
#include <boost/dll.hpp>

#include "Channel.h"

namespace Gadgetron::Core::Parallel {

    class Merge {
    public:
        virtual ~Merge() = default;
        virtual void process(std::map<std::string, std::shared_ptr<Channel>>, std::shared_ptr<Channel>) = 0;
    };
}


#define GADGETRON_MERGE_EXPORT(MergeClass)                                    \
std::unique_ptr<Gadgetron::Core::Parallel::Merge>                             \
merge_factory_##MergeClass() {                                                \
    return std::make_unique<MergeClass>();                                    \
}                                                                             \
                                                                              \
BOOST_DLL_ALIAS(                                                              \
        merge_factory_##MergeClass,                                           \
        merge_factory_export_##MergeClass                                     \
)                                                                             \
