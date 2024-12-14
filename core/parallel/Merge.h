#pragma once

#include <map>
#include <memory>
#include <boost/dll.hpp>

#include "Channel.h"
#include "Context.h"
#include "PropertyMixin.h"

namespace Gadgetron::Core::Parallel {

    class Merge : public PropertyMixin {
    public:
        explicit Merge(const GadgetProperties &props);

        virtual ~Merge() = default;
        virtual void process(std::map<std::string, GenericInputChannel>, OutputChannel) = 0;
    };
}

#define GADGETRON_MERGE_EXPORT(MergeClass)                                      \
std::unique_ptr<Gadgetron::Core::Parallel::Merge>                               \
merge_factory_##MergeClass(                                                     \
        const Gadgetron::Core::Context &context,                                \
        const std::string& name,                                                \
        const Gadgetron::Core::GadgetProperties &props                          \
) {                                                                             \
    return std::make_unique<MergeClass>(context, props);                        \
}                                                                               \
                                                                                \
BOOST_DLL_ALIAS(                                                                \
        merge_factory_##MergeClass,                                             \
        merge_factory_export_##MergeClass                                       \
)
