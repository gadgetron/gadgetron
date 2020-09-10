//
// Created by david on 24/08/2020.
//

#include "FlagTriggerGadget.h"
#include <boost/algorithm/string.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include "io/from_string.h"
#include "ChannelAlgorithms.h"
#include "mri_core_data.h"
#include "mri_core_acquisition_bucket.h"



#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

using namespace Gadgetron;
static const std::unordered_map<std::string, Gadgetron::FlagTriggerGadget::TriggerFlags> flag_map =
    {

        {"first_in_encode_step1", FlagTriggerGadget::TriggerFlags::first_in_encode_step1},
        {"last_in_encode_step1", FlagTriggerGadget::TriggerFlags::last_in_encode_step1},
        {"first_in_encode_step2", FlagTriggerGadget::TriggerFlags::first_in_encode_step2},
        {"last_in_encode_step2", FlagTriggerGadget::TriggerFlags::last_in_encode_step2},
        {"first_in_average", FlagTriggerGadget::TriggerFlags::first_in_average},
        {"last_in_average", FlagTriggerGadget::TriggerFlags::last_in_average},
        {"first_in_slice", FlagTriggerGadget::TriggerFlags::first_in_slice},
        {"last_in_slice", FlagTriggerGadget::TriggerFlags::last_in_slice},
        {"first_in_contrast", FlagTriggerGadget::TriggerFlags::first_in_contrast},
        {"last_in_contrast", FlagTriggerGadget::TriggerFlags::last_in_contrast},
        {"first_in_phase", FlagTriggerGadget::TriggerFlags::first_in_phase},
        {"last_in_phase", FlagTriggerGadget::TriggerFlags::last_in_phase},
        {"first_in_repetition", FlagTriggerGadget::TriggerFlags::first_in_repetition},
        {"last_in_repetition", FlagTriggerGadget::TriggerFlags::last_in_repetition},
        {"first_in_set", FlagTriggerGadget::TriggerFlags::first_in_set},
        {"last_in_set", FlagTriggerGadget::TriggerFlags::last_in_set},
        {"first_in_segment", FlagTriggerGadget::TriggerFlags::first_in_segment},
        {"last_in_segment", FlagTriggerGadget::TriggerFlags::last_in_segment},
        {"is_noise_measurement", FlagTriggerGadget::TriggerFlags::is_noise_measurement},
        {"is_parallel_calibration", FlagTriggerGadget::TriggerFlags::is_parallel_calibration},
        {"is_parallel_calibration_and_imaging",
         FlagTriggerGadget::TriggerFlags::is_parallel_calibration_and_imaging},
        {"is_reverse", FlagTriggerGadget::TriggerFlags::is_reverse},
        {"is_navigation_data", FlagTriggerGadget::TriggerFlags::is_navigation_data},
        {"is_phasecorr_data", FlagTriggerGadget::TriggerFlags::is_phasecorr_data},
        {"last_in_measurement", FlagTriggerGadget::TriggerFlags::last_in_measurement},
        {"is_hpfeedback_data", FlagTriggerGadget::TriggerFlags::is_hpfeedback_data},
        {"is_dummyscan_data", FlagTriggerGadget::TriggerFlags::is_dummyscan_data},
        {"is_rtfeedback_data", FlagTriggerGadget::TriggerFlags::is_rtfeedback_data},
        {"is_surfacecoilcorrectionscan_data",
         FlagTriggerGadget::TriggerFlags::is_surfacecoilcorrectionscan_data},

        {"compression1", FlagTriggerGadget::TriggerFlags::compression1},
        {"compression2", FlagTriggerGadget::TriggerFlags::compression2},
        {"compression3", FlagTriggerGadget::TriggerFlags::compression3},
        {"compression4", FlagTriggerGadget::TriggerFlags::compression4},
        {"user1", FlagTriggerGadget::TriggerFlags::user1},
        {"user2", FlagTriggerGadget::TriggerFlags::user2},
        {"user3", FlagTriggerGadget::TriggerFlags::user3},
        {"user4", FlagTriggerGadget::TriggerFlags::user4},
        {"user5", FlagTriggerGadget::TriggerFlags::user5},
        {"user6", FlagTriggerGadget::TriggerFlags::user6},
        {"user7", FlagTriggerGadget::TriggerFlags::user7},
        {"user8", FlagTriggerGadget::TriggerFlags::user8}};

void Gadgetron::from_string(const std::string& str, FlagTriggerGadget::TriggerFlags& flag) {

    auto lower = boost::algorithm::to_lower_copy(str);
    using namespace std::string_literals;
    if (flag_map.count(lower) < 1)
        throw std::runtime_error("The string \""s + str + "\" is not a valid TriggerFlag.");

    flag = flag_map.at(lower);
}

void Gadgetron::FlagTriggerGadget::process(Core::InputChannel<Core::Acquisition>& in,
                                           Core::OutputChannel& out) {

    for (const auto& group : Core::Algorithm::buffer(in, [this](const Core::Acquisition& acq) {
             const auto& [header, data, traj] = acq;
             return (std::bitset<64>(header.flags) & flags).any();
         })) {
        auto bucket = AcquisitionBucket();
        for (auto acq : group){
            bucket.add_acquisition(std::move(acq));
        }
        out.push(std::move(bucket));
    }
}


namespace {
struct AcquisitionFlags_ : boost::spirit::qi::symbols<char,FlagTriggerGadget::TriggerFlags>{
    AcquisitionFlags_() {
        for (auto [key,value] : flag_map) this->add(key,value);
    }

};

template<typename Iterator>
std::bitset<64> parse_triggers(Iterator first, Iterator last) {
    using namespace boost::spirit::qi;
    namespace qi = boost::spirit::qi;
    using boost::spirit::ascii::space;
    AcquisitionFlags_ flag;
    std::vector<FlagTriggerGadget::TriggerFlags> v;
    bool r = phrase_parse(first, last,flag % '|',  space, v);
    if (!r || first != last) // fail if we did not get a full match
        throw std::runtime_error("Unable to parse trigger flags");
    using namespace ranges;
    return accumulate( v | views::transform([](auto triggerflag ){return static_cast<uint64_t>(triggerflag);}), uint64_t(0), std::bit_or());


}
}

Gadgetron::FlagTriggerGadget::FlagTriggerGadget(const Core::Context& context,
                                                const Core::GadgetProperties& props)
    : ChannelGadget(context, props) {
    using namespace ranges;
    flags = parse_triggers(trigger_flags.begin(),trigger_flags.end());
}