#pragma once

#include <map>
#include <memory>
#include <boost/dll.hpp>

#include "Channel.h"
#include "Context.h"
#include "PropertyMixin.h"

namespace Gadgetron::Core::Parallel {

    class Branch : public PropertyMixin {
    public:
        virtual ~Branch() = default;
        virtual void process(
                std::shared_ptr<Channel> input,
                const std::map<std::string, std::shared_ptr<Channel>> &output,
                std::shared_ptr<Channel> bypass
        ) = 0;

    private:
        explicit Branch(const GadgetProperties &props);

        template<class...> friend class TypedBranch;
    };

    template<class ...ARGS>
    class TypedBranch : public Branch {
    public:
        explicit TypedBranch(const GadgetProperties &props);

        void process(
                std::shared_ptr<Channel> input,
                const std::map<std::string, std::shared_ptr<Channel>> &output,
                std::shared_ptr<Channel> bypass
        ) final;

        virtual void process(TypedInputChannel<ARGS...> &, std::map<std::string, std::shared_ptr<OutputChannel>>) = 0;
    };

}

#include "Branch.hpp"

#define GADGETRON_BRANCH_EXPORT(BranchClass)                                    \
std::unique_ptr<Gadgetron::Core::Parallel::Branch>                              \
branch_factory_##BranchClass(                                                   \
        const Gadgetron::Core::Context &context,                                \
        const Gadgetron::Core::GadgetProperties &props                          \
) {                                                                             \
    return std::make_unique<BranchClass>(context, props);                       \
}                                                                               \
                                                                                \
BOOST_DLL_ALIAS(                                                                \
        branch_factory_##BranchClass,                                           \
        branch_factory_export_##BranchClass                                     \
)                                                                               \
