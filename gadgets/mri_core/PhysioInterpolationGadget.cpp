#include "PhysioInterpolationGadget.h"
#include "GadgetronTimer.h"
#include "Spline.h"
#include "mri_core_def.h"
#include "ismrmrd/meta.h"
#include "hoNDBSpline.h"
#include "ismrmrd/xml.h"

#include <numeric>
#include <queue>
#include <range/v3/view/transform.hpp>
#include <range/v3/range/conversion.hpp>

#ifdef USE_OMP

#include <omp.h>

#endif

namespace Gadgetron {

    namespace {

        struct Intervals {
            std::vector<float> intervals;
            std::vector<size_t> cycle_starts;
        };

        Intervals find_intervals(const std::vector<float> &time_stamps, bool first_beat_on_trigger) {

            float previous = -100.0;
            std::vector<float> intervals;
            std::vector<size_t> cycle_starts;
            for (size_t i = 0; i < time_stamps.size(); i++) {
                if ((time_stamps[i] < previous) || (first_beat_on_trigger && i == 0)) {
                    cycle_starts.push_back(i);
                } else if (i > 0) {
                    intervals.push_back(time_stamps[i] - time_stamps[i - 1]);
                }
                previous = time_stamps[i];
            }
            return {intervals, cycle_starts};
        };

        std::vector<float>
        calculate_cycle_lengths(const std::vector<float> &time_stamps, float median_interval,
                                const std::vector<size_t> &cycle_starts) {
            std::vector<float> cycle_lengths;
            float count = -1;
            for (size_t i = 0; i < cycle_starts.size(); i++) {
                float clength = time_stamps[cycle_starts[i] - 2] + median_interval - time_stamps[cycle_starts[i]];
                cycle_lengths.push_back(clength);
            }

            if (cycle_lengths.empty()) {
                float clength = time_stamps.back();
                cycle_lengths.push_back(clength);
            }

            std::sort(cycle_lengths.begin(), cycle_lengths.end());

            return cycle_lengths;
        }


        std::vector<float> calculate_relative_cycle_time(const std::vector<size_t> &cycle_starts,
                                                         const std::vector<float> &cycle_lengths,
                                                         const std::vector<float> &time_stamps) {//Calculate relative time stamps
            std::vector<float> relative_cycle_time;
            size_t current_cycle = 0;
            for (size_t i = 0; i < time_stamps.size(); i++) {
                if ((current_cycle < cycle_starts.size()) && (i >= cycle_starts[current_cycle]) &&
                    (current_cycle < cycle_starts.size())) {
                    current_cycle++;
                }
                relative_cycle_time.push_back(time_stamps[i] / cycle_lengths[current_cycle - 1] + current_cycle);
            }
            return relative_cycle_time;
        }

        auto spline_interpolate_series = [](const auto &buffer, const auto &relative_cycle_time,
                                            const auto &recon_cycle_time,
                                            std::vector<Core::Image<std::complex<float>>> &output) {
            GadgetronTimer interptime("Interpolation Time");


            //Let's interpolate the images
            size_t inelem = relative_cycle_time.size();
            size_t outelem = recon_cycle_time.size();
            size_t imageelem = std::get<1>(buffer.front()).size();

#ifdef USE_OMP
#pragma omp parallel for
#endif
            for (long long p = 0; p < (long long) imageelem; p++) {
                std::vector<std::complex<float> > data_in(inelem);
                //Get the input data for this pixel
                for (size_t i = 0; i < inelem; i++) data_in[i] = std::get<1>(buffer[i])[p];

                //Interpolate the data
                Spline<float, std::complex<float> > sp(relative_cycle_time, data_in);
                std::vector<std::complex<float> > data_out = sp[recon_cycle_time];

                //Copy it to the images
                for (size_t i = 0; i < outelem; i++) get<1>(output[i])[p] = data_out[i];
            }
        };

        auto bspline_interpolate_series = [](const auto &buffer, const auto &relative_cycle_time,
                                             const auto &recon_cycle_time,
                                             std::vector<Core::Image<std::complex<float>>> &output) {

            GadgetronTimer interptime("Interpolation Time using BSpline");
            size_t inelem = relative_cycle_time.size();
            size_t outelem = recon_cycle_time.size();
            size_t imageelem = std::get<1>(buffer.front()).size();
            size_t SplineDegree = 5;

            long long p;
#pragma omp parallel default(none) shared(SplineDegree, imageelem, inelem, outelem, buffer, relative_cycle_time, recon_cycle_time, output) private(p)
            {
                hoNDArray<std::complex<float> > data_in(inelem);
                hoNDArray<std::complex<float> > data_out(outelem);

                hoNDArray<std::complex<float> > coeff(inelem);

                hoNDBSpline<std::complex<float>, 1> interp;

                size_t i;

                size_t num = relative_cycle_time.size();

#pragma omp for
                for (p = 0; p < (long long) imageelem; p++) {
                    //Get the input data for this pixel
                    for (i = 0; i < inelem; i++) data_in[i] = std::get<1>(buffer[i])[p];

                    // compute the coefficient
                    interp.computeBSplineCoefficients(data_in, SplineDegree, coeff);

                    //Interpolate the data
                    for (i = 0; i < outelem; i++) {
                        float x = (num - 1) * (recon_cycle_time[i] - relative_cycle_time[0]) /
                                  (relative_cycle_time[num - 1] - relative_cycle_time[0]);
                        data_out(i) = interp.evaluateBSpline(coeff.begin(), inelem, SplineDegree, 0, x);
                    }

                    //Copy it to the images
                    for (i = 0; i < outelem; i++) std::get<1>(output[i])[p] = data_out[i];
                }
            }
        };

    }


    void PhysioInterpolationGadget::process(Core::InputChannel<Core::Image<std::complex<float>>> &in,
                                            Core::OutputChannel &out) {


        ISMRMRD::EncodingLimits e_limits = header.encoding[0].encodingLimits;
        auto slc_limit = e_limits.slice ? e_limits.slice->maximum + 1 : 1;


        auto buffers = std::map<int, std::vector<Core::Image<std::complex<float>>>>{};
        auto time_stamp_buffer = std::map<int, std::vector<float>>{};

        for (auto[hdr, data, meta]: in) {
            buffers[hdr.slice].emplace_back(hdr, data, meta);
            time_stamp_buffer[hdr.slice].push_back((float) (hdr.physiology_time_stamp[physiology_time_index]));
            out.push(Core::Image<std::complex<float>>{hdr, std::move(data), std::move(meta)});
        }

        for (auto &key_val: buffers) {
            auto slc = key_val.first;
            auto &buffer = key_val.second;
            GDEBUG("Processing slice: %d ... \n", slc);
            GDEBUG("Number of items on Q: %d\n", buffer.size());
            GDEBUG("Image with attribute flag : %d\n", bool(std::get<2>(buffer.front())));


            auto time_stamps = time_stamp_buffer[slc];

            auto[intervals, cycle_starts] = find_intervals(time_stamps, first_beat_on_trigger);

            if (intervals.empty()) continue;

            std::sort(intervals.begin(), intervals.end());

            float mean_interval = std::accumulate(intervals.begin(), intervals.end(), 0.0f) / intervals.size();
            float median_interval = intervals[(intervals.size() >> 1)];
            std::vector<float> cycle_lengths = calculate_cycle_lengths(time_stamps, median_interval, cycle_starts);

            float mean_cycle_length =
                    std::accumulate(cycle_lengths.begin(), cycle_lengths.end(), 0.0f) / cycle_lengths.size();
            float median_cycle_length = cycle_lengths[(cycle_lengths.size() >> 1)];
            //Make sure we have cycle lengths for all the cycles we have covered
            cycle_lengths.insert(cycle_lengths.begin(), median_cycle_length);
            cycle_lengths.push_back(median_cycle_length);


            GDEBUG("We have %d full cyles, first one starting at %d\n", cycle_starts.size() - 1, cycle_starts[0]);
            GDEBUG("Mean/Median frame width %f/%f\n", mean_interval, median_interval);
            GDEBUG("Mean/Median cycle_length %f/%f\n", mean_cycle_length, median_cycle_length);

            //Correct the first cycle assuming it is of median length:
            if (!first_beat_on_trigger) {
                float first_cycle_offset = (median_cycle_length - median_interval) + time_stamps[cycle_starts[0]] -
                                           time_stamps[cycle_starts[0] - 1];
                std::transform(time_stamps.begin(), time_stamps.end(), time_stamps.begin(),
                               [=](auto val) { return val + first_cycle_offset; });
            }
            std::vector<float> relative_cycle_time = calculate_relative_cycle_time(cycle_starts, cycle_lengths,
                                                                                   time_stamps);


            //Let's figure out which time points we would like to interpolate on:
            ///TODO: Deal with mode 1 and other future modes, we are only implementing mode 0 at the moment
            float phase_interval = 1.0f / static_cast<float>(phases);
            float max_time = std::floor(relative_cycle_time[relative_cycle_time.size() - 1]);
            std::vector<float> recon_cycle_time;
            for (float t = 1.0; t < (max_time - 0.001); t += phase_interval) {
                recon_cycle_time.push_back(t);
            }

            if (mode == PhysioInterpolationMode::complete) {
                recon_cycle_time.resize(phases);
            }

            //Now we can loop over each pixel and estimate the new frames, but first we have to have somewhere to put the data


            using namespace ranges;

            auto image_generator = [&](float cycle_time) -> Core::Image<std::complex<float>> {
                const auto&[ref_header, ref_data, ref_meta] = buffer.front();
                auto header = ref_header;
                auto data = hoNDArray<std::complex<float>>(ref_data.dimensions());
                auto meta = ref_meta;

                unsigned short current_cycle = static_cast<unsigned short>(std::floor(cycle_time + 0.0001));
                unsigned short current_phase = static_cast<unsigned short>(
                        (cycle_time + 0.0001 - current_cycle) / (1.0 / static_cast<float>(phases)) + 0.0001);

                header.physiology_time_stamp[physiology_time_index] = static_cast<unsigned>(std::floor(
                        (cycle_time + 0.0001 - current_cycle) * cycle_lengths[current_cycle]));
                header.phase = current_phase;
                header.image_index = current_phase + 1;
                header.image_series_index = current_cycle * 10 + header.slice;

                if (header.phase + 1 >= time_stamps.size()) header.phase = uint16_t(time_stamps.size() - 1);
                if (ref_meta) {
                    meta->set("PHS", long(header.phase));
                    meta->set(GADGETRON_IMAGENUMBER, long(header.image_index));
                    meta->append(GADGETRON_DATA_ROLE, "PhysionInterp");

                    double cycle_length_in_ms = time_stamp_resolution_ * cycle_lengths[current_cycle];
                    std::ostringstream ostr;
                    if (slc_limit > 1) {
                        ostr << "_SLC_" << header.slice << "_RR" << cycle_length_in_ms << "ms";
                    } else {
                        ostr << "_RR" << cycle_length_in_ms << "ms";
                    }

                    std::string imageComment = "PhysioInterp" + ostr.str();
                    meta->append(GADGETRON_IMAGECOMMENT, imageComment.c_str());

                    std::string seqDescription = "_PhysioInterp" + ostr.str();
                    meta->append(GADGETRON_SEQUENCEDESCRIPTION, seqDescription.c_str());

                    meta->append(GADGETRON_IMAGEPROCESSINGHISTORY, "Interp");
                }
                return {header, data, meta};
            };

            auto output = ranges::transform_view(recon_cycle_time, image_generator) | to<std::vector>;


            if ((interp_method == PhysioInterpolationMethod::Spline) || (mode != PhysioInterpolationMode::complete)) {
                spline_interpolate_series(buffer, relative_cycle_time, recon_cycle_time, output);

            } else {
                bspline_interpolate_series(buffer, relative_cycle_time,recon_cycle_time,output);

            }

            std::move(output.begin(),output.end(),out.begin());
        }
    }

    GADGETRON_GADGET_EXPORT(PhysioInterpolationGadget)
}
