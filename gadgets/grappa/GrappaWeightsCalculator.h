#pragma once

#include "gadgetron_grappa_export.h"
#include "GrappaWeights.h"

#include <list>

namespace Gadgetron {


    class EXPORTGADGETSGRAPPA GrappaWeightsCalculator : public ACE_Task<ACE_MT_SYNCH> {
        using inherited = ACE_Task<ACE_MT_SYNCH>;

    public:
        GrappaWeightsCalculator()
                : inherited(), target_coils_(0) {}

        virtual ~GrappaWeightsCalculator() = default;

        virtual int init() {
            return 0;
        }

        int open(void * = nullptr) {
            return this->activate(THR_NEW_LWP | THR_JOINABLE, 1);
        }

        int close(unsigned long flags);

        int svc(void);

        int add_job(hoNDArray<std::complex<float>> *ref_data,
                    std::vector<std::pair<unsigned int, unsigned int> > sampled_region,
                    unsigned int acceleration_factor,
                    boost::shared_ptr<GrappaWeights<float>> destination,
                    std::vector<unsigned int> uncombined_channel_weights,
                    bool include_uncombined_channels_in_combined_weights = true);

        int add_uncombined_channel(unsigned int channel_id);
        int remove_uncombined_channel(unsigned int channel_id);

        int get_number_of_uncombined_channels() {
            return uncombined_channels_.size();
        }

        int get_number_of_target_coils() {
            return target_coils_;
        }

        void set_number_of_target_coils(int n) {
            target_coils_ = n;
        }

        bool get_use_gpu() {
            return use_gpu_;
        }

        void set_use_gpu(bool v) {
            use_gpu_ = v;
        }

    private:
        std::list<unsigned int> uncombined_channels_;
        int target_coils_;
        bool use_gpu_;

        hoNDArray<std::complex<float>> target_acs_;
        hoNDArray<std::complex<float>> complex_im_;
        hoNDArray<std::complex<float>> conv_ker_;
        hoNDArray<std::complex<float>> kIm_;
        hoNDArray<std::complex<float>> coil_map_;
        hoNDArray<std::complex<float>> unmixing_;
        hoNDArray<float> gFactor_;

        class WeightsDescription {
        public:
            boost::shared_ptr<GrappaWeights<float>> destination;

            std::vector<std::pair<unsigned int, unsigned int>> sampled_region;
            std::vector<unsigned int> uncombined_channel_weights;

            unsigned int acceleration_factor = 1;
            bool include_uncombined_channels_in_combined_weights = true;
        };

        void push_update(const WeightsDescription *description, const hoNDArray <std::complex<float>> *host_data) const;
    };
}
