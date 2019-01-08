#include "Branch.h"

#include "log.h"

namespace {


}

namespace Gadgetron::Core::Parallel {

    class FanoutBranch : public Branch {

    public:
        void process(std::shared_ptr<Channel> in, std::map<std::string, std::shared_ptr<Channel>> out) override {

            GINFO_STREAM("Fanout Branch process running.");

            InputChannel &input = *in;
            OutputChannel &output = *out.at("grappa-recon");

            for (auto message : input) {
                output.push_message(std::move(message));
            }

            GINFO_STREAM("Input stream closed. Closing output.");

            out.at("grappa-recon")->close();
            out.at("sense-spirit")->close();
        }
    };

    GADGETRON_BRANCH_EXPORT(FanoutBranch);
}

