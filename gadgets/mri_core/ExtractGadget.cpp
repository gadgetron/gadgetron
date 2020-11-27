/*
 * ExtractMagnitudeGadget.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: Michael S. Hansen
 */

#include "ExtractGadget.h"
#include <bitset>
#include <cpu/math/hoNDArray_math.h>
#include <unordered_map>

#include <boost/math/constants/constants.hpp>

#include "hoNDArray_fileio.h"

namespace Gadgetron {

    namespace {
        using IMTYPE = ISMRMRD::ISMRMRD_ImageTypes;

        const std::map<IMTYPE, size_t> series_offset{ { IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE, 0 },
            { IMTYPE::ISMRMRD_IMTYPE_REAL, 1000 }, { IMTYPE::ISMRMRD_IMTYPE_IMAG, 2000 },
            { IMTYPE::ISMRMRD_IMTYPE_PHASE, 3000 } };

        const std::map<int, IMTYPE> mask_to_imtype{ { 0, IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE },
            { 1, IMTYPE::ISMRMRD_IMTYPE_REAL }, { 2, IMTYPE::ISMRMRD_IMTYPE_IMAG },
            { 3, IMTYPE::ISMRMRD_IMTYPE_PHASE } };

        template <class FUNCTION>
        hoNDArray<float> extract(const hoNDArray<std::complex<float>>& data, FUNCTION&& extractor) {
            hoNDArray<float> output(data.dimensions());
            std::transform(data.begin(), data.end(), output.begin(), extractor);
            return output;
        }

        hoNDArray<float> extract_image(const hoNDArray<std::complex<float>>& data, IMTYPE imtype, float offset) {
            switch (imtype) {
            case ISMRMRD::ISMRMRD_IMTYPE_MAGNITUDE: return extract(data, [](auto& val) { return std::abs(val); });
            case ISMRMRD::ISMRMRD_IMTYPE_PHASE:
                return extract(data, [](auto& val) { return std::arg(val) + boost::math::constants::pi<float>(); });
            case ISMRMRD::ISMRMRD_IMTYPE_REAL: return extract(data, [&](auto& val) { return std::real(val) + offset; });
            case ISMRMRD::ISMRMRD_IMTYPE_IMAG: return extract(data, [&](auto& val) { return std::imag(val) + offset; });
            default: throw std::runtime_error("Illegal image tpye encounted in extract_image");
            }
        }

    }

    void ExtractGadget::process(Core::InputChannel<Core::Image<std::complex<float>>>& in, Core::OutputChannel& out) {
        for (auto image : in) {
            const auto& head = std::get<ISMRMRD::ImageHeader>(image);
            const auto& data = std::get<hoNDArray<std::complex<float>>>(image);
            const auto& meta = std::get<2>(image);

            for (auto imtype : image_types) {
                auto data_copy       = extract_image(data, imtype, real_imag_offset);
                auto head_copy       = head;
                head_copy.data_type  = ISMRMRD::ISMRMRD_FLOAT;
                head_copy.image_type = imtype;
                head_copy.image_series_index += series_offset.at(imtype);
                out.push(head_copy, std::move(data_copy), meta);
            }
        }
    }

    ExtractGadget::ExtractGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : Core::ChannelGadget<Core::Image<std::complex<float>>>(context, props) {

        for (int i = 0; i < extract_mask.size(); i++) {
            if (extract_mask[i])
                image_types.insert(mask_to_imtype.at(i));
        }
        if (extract_magnitude)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_MAGNITUDE);
        if (extract_real)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_REAL);
        if (extract_imag)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_IMAG);
        if (extract_phase)
            image_types.insert(IMTYPE::ISMRMRD_IMTYPE_PHASE);

        if (image_types.empty())
            throw std::runtime_error("ExctractGadget: No valid extract functions specified");
    }

    GADGETRON_GADGET_EXPORT(ExtractGadget)

}
