
#include "Merge.h"

#include "log.h"

namespace {


}

namespace Gadgetron::Core::Parallel {

    class SimpleMerge : public Merge {

    public:
        void process(std::map<std::string, std::shared_ptr<Channel>> in, std::shared_ptr<Channel> out) override {

            GINFO_STREAM("Simple Merge process running.");

            std::vector<std::thread> threads;

            for (auto &pair : in) {
                threads.emplace_back(std::thread(
                        [](std::shared_ptr<InputChannel> in, std::shared_ptr<OutputChannel> out) {
                                for (auto message : *in) {
                                    out->push_message(std::move(message));
                                }
                        },
                        pair.second,
                        out
                ));
            }

            for (auto &thread : threads) { thread.join(); }

            GINFO_STREAM("Input streams closed. Closing output.");

            out->close();
        }
    };

    GADGETRON_MERGE_EXPORT(SimpleMerge);
}
