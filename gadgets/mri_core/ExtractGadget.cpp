/*
 * ExtractMagnitudeGadget.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: Michael S. Hansen
 */

#include "ExtractGadget.h"
#include <bitset>
#include <unordered_map>

#include <boost/math/constants/constants.hpp>


namespace Gadgetron {

    namespace {
        using IMTYPE = mrd::ImageType;

        const std::map<IMTYPE, size_t> series_offset{ { IMTYPE::kMagnitude, 0 },
            { IMTYPE::kReal, 1000 }, { IMTYPE::kImag, 2000 },
            { IMTYPE::kPhase, 3000 } };

        const std::map<int, IMTYPE> mask_to_imtype{ { 0, IMTYPE::kMagnitude },
            { 1, IMTYPE::kReal }, { 2, IMTYPE::kImag },
            { 3, IMTYPE::kPhase } };

        template <class FUNCTION>
        hoNDArray<float> extract(const hoNDArray<std::complex<float>>& data, FUNCTION&& extractor) {
            hoNDArray<float> output(data.dimensions());
            std::transform(data.begin(), data.end(), output.begin(), extractor);
            return output;
        }

        mrd::ImageData<float> extract_image(const mrd::ImageData<std::complex<float>>& data, IMTYPE imtype, float offset) {
            switch (imtype) {
            case IMTYPE::kMagnitude:
                return extract(data, [](auto& val) { return std::abs(val); });
            case IMTYPE::kPhase:
                return extract(data, [](auto& val) { return std::arg(val) + boost::math::constants::pi<float>(); });
            case IMTYPE::kReal:
                return extract(data, [&](auto& val) { return std::real(val) + offset; });
            case IMTYPE::kImag:
                return extract(data, [&](auto& val) { return std::imag(val) + offset; });
            default:
                throw std::runtime_error("Illegal image type encountered in extract_image");
            }
        }

    }

    void ExtractGadget::process(Core::InputChannel<mrd::Image<std::complex<float>>>& in, Core::OutputChannel& out) {
        size_t count = 0;
        for (auto image : in) {
            for (auto imtype : image_types) {
                mrd::Image<float> extracted;
                extracted.head = image.head;
                extracted.head.image_type = imtype;
                extracted.head.image_series_index = image.head.image_series_index.value_or(0) + series_offset.at(imtype);
                extracted.data = extract_image(image.data, imtype, real_imag_offset);
                extracted.meta = image.meta;

                out.push(std::move(extracted));
            }
        }
    }

    ExtractGadget::ExtractGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : Core::ChannelGadget<mrd::Image<std::complex<float>>>(context, props) {

        for (int i = 0; i < extract_mask.size(); i++) {
            if (extract_mask[i])
                image_types.insert(mask_to_imtype.at(i));
        }
        if (extract_magnitude)
            image_types.insert(mrd::ImageType::kMagnitude);
        if (extract_real)
            image_types.insert(mrd::ImageType::kReal);
        if (extract_imag)
            image_types.insert(mrd::ImageType::kImag);
        if (extract_phase)
            image_types.insert(mrd::ImageType::kPhase);

        if (image_types.empty())
            throw std::runtime_error("ExctractGadget: No valid extract functions specified");
    }

    GADGETRON_GADGET_EXPORT(ExtractGadget)

}
