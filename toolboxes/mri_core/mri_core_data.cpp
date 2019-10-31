#include "mri_core_data.h"

#include <unordered_map>

using namespace Gadgetron;

static const std::unordered_map<std::string,Gadgetron::IsmrmrdCondition> condition_map = {
    { "kspace_encode_step_1",IsmrmrdCONDITION::KSPACE_ENCODE_STEP_1},
    { "kspace_encode_step_2",IsmrmrdCONDITION::KSPACE_ENCODE_STEP_2},
    {"average",IsmrmrdCONDITION::AVERAGE},
    {"slice",IsmrmrdCONDITION::SLICE},
    {"contrast",IsmrmrdCONDITION::CONTRAST},
    {"phase",IsmrmrdCONDITION::PHASE},
    {"repetition",IsmrmrdCONDITION::REPETITION},
    {"set",IsmrmrdCONDITION::SET},
    {"segment",IsmrmrdCONDITION::SEGMENT},
    {"user_0",IsmrmrdCONDITION::USER_0},
    {"user_1",IsmrmrdCONDITION::USER_1},
    {"user_2",IsmrmrdCONDITION::USER_2},
    {"user_3",IsmrmrdCONDITION::USER_3},
    {"user_4",IsmrmrdCONDITION::USER_4},
    {"user_5",IsmrmrdCONDITION::USER_5},
    {"user_6",IsmrmrdCONDITION::USER_6},
    {"user_7",IsmrmrdCONDITION::USER_7},
    {"n_acquisitions",IsmrmrdCONDITION::N_ACQUISITIONS},
    {"none",IsmrmrdCONDITION::NONE},
    {"",IsmrmrdCONDITION::NONE}

};

void Gadgetron::from_string(const std::string& stringval, IsmrmrdCondition& condition){
    std::string stringval_lower = stringval;
    for(auto& s : stringval_lower) s = std::to_lower(s);
    if (!condition_map.count(stringval_lower)) throw std::runtime_error("Could not convert string to IsmrmrdCondition");
    return condition_map[stringval_lower];


}
