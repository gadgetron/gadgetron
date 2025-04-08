#pragma once
#include "Gadget.h"
#include "hoNDArray.h"
#include "ismrmrd/meta.h"
#include "ismrmrd/ismrmrd.h"

namespace Gadgetron
{

    /**
    * This Gadget adds MRD meta attributes from Gadgetron attributes defined in mri_core_def.h.
    */

    template <typename T >
    class GenericMRDCommentGadget : public Core::ChannelGadget<Core::Image<T>>
    {
        public:
            using Core::ChannelGadget<Core::Image<T>>::ChannelGadget;

            GenericMRDCommentGadget(const Core::Context& context, const Core::GadgetProperties& props);
            ~GenericMRDCommentGadget() override = default;

            void process(Core::InputChannel<Core::Image<T>>& input, Core::OutputChannel& output) override;

        protected:

            // the MRD convention is defined at https://ismrmrd.readthedocs.io/en/latest/mrd_image_data.html
            std::map<std::string, std::string> dict_gt_2_mrd_;
    };

    class GenericMRDCommentShortGadget :public GenericMRDCommentGadget<short>
    {
    public:
        using GenericMRDCommentGadget<short>::GenericMRDCommentGadget;
        ~GenericMRDCommentShortGadget() override = default;
    };

    class GenericMRDCommentUShortGadget :public GenericMRDCommentGadget<unsigned short>
    {
    public:
        using GenericMRDCommentGadget<unsigned short>::GenericMRDCommentGadget;
        ~GenericMRDCommentUShortGadget() override = default;
    };

    class GenericMRDCommentIntGadget :public GenericMRDCommentGadget<int>
    {
    public:
        using GenericMRDCommentGadget<int>::GenericMRDCommentGadget;
        ~GenericMRDCommentIntGadget() override = default;
    };

    class GenericMRDCommentUIntGadget :public GenericMRDCommentGadget<unsigned int>
    {
    public:
        using GenericMRDCommentGadget<unsigned int>::GenericMRDCommentGadget;
        ~GenericMRDCommentUIntGadget() override = default;
    };

    class GenericMRDCommentFloatGadget :public GenericMRDCommentGadget<float>
    {
    public:
        using GenericMRDCommentGadget<float>::GenericMRDCommentGadget;
        ~GenericMRDCommentFloatGadget() override = default;
    };

    class GenericMRDCommentCxFloatGadget :public GenericMRDCommentGadget<std::complex<float>>
    {
    public:
        using GenericMRDCommentGadget<std::complex<float>>::GenericMRDCommentGadget;
        ~GenericMRDCommentCxFloatGadget() override = default;
    };

    class GenericMRDCommentCxDoubleGadget :public GenericMRDCommentGadget<std::complex<double>>
    {
    public:
        using GenericMRDCommentGadget<std::complex<double>>::GenericMRDCommentGadget;
        ~GenericMRDCommentCxDoubleGadget() override = default;
    };
}

