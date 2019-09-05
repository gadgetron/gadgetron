#pragma once

#include "../Channel.h"
#include "ChannelCreator.h"
#include "PropertyMixin.h"
#include <boost/dll/alias.hpp>

namespace Gadgetron::Core::Distributed {
    class Distributor : public PropertyMixin {

    public:
        virtual void process(GenericInputChannel input, ChannelCreator &creator, OutputChannel bypass) = 0;

        virtual ~Distributor() = default;

    private:
        explicit Distributor(const GadgetProperties &properties);

        template<class...> friend class TypedDistributor;
    };

    template<class... ARGS>
    class TypedDistributor : public Distributor {

    public:

        explicit TypedDistributor(const GadgetProperties &properties);

    void process(GenericInputChannel input, ChannelCreator &creator, OutputChannel bypass) final;

        virtual void process(InputChannel<ARGS...> &input, ChannelCreator &creator) = 0;

        ~TypedDistributor() override = default;
    };

    template<class... ARGS>
    void TypedDistributor<ARGS...>::process(GenericInputChannel input, ChannelCreator &creator,
                                            OutputChannel bypass) {
        InputChannel< ARGS...> typed_input(input, bypass);
        process(typed_input, creator);
    }

    template<class... ARGS>
    Gadgetron::Core::Distributed::TypedDistributor<ARGS...>::TypedDistributor(const GadgetProperties &properties)
            : Distributor(properties) {}
}


#define GADGETRON_DISTRIBUTOR_EXPORT(DistributorClass)                                                            \
std::unique_ptr<Gadgetron::Core::Distributed::Distributor> distributor_factory_##DistributorClass(                \
        const Gadgetron::Core::Context &context,                                                                  \
        const Gadgetron::Core::GadgetProperties &props                                                            \
) {                                                                                 \
    return std::make_unique<DistributorClass>(context, props);                      \
}                                                                                   \
                                                                                    \
BOOST_DLL_ALIAS(                                                                    \
        distributor_factory_##DistributorClass,                                     \
        distributor_factory_export_##DistributorClass                               \
)                                                                                   \


