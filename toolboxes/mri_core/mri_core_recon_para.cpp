/** \file   mri_core_recon_para.cpp
    \brief  Define the structures storing the reconstruction parameters for different algorithms
    \author Hui Xue
*/

#include "mri_core_recon_para.h"

namespace Gadgetron
{
    // --------------------------------------------------------------------

    std::string get_ismrmrd_dimension_name(IsmrmrdCONDITION v)
    {
        std::string name;

        switch (v)
        {
        case KSPACE_ENCODE_STEP_1:
            name = "E1";
            break;

        case KSPACE_ENCODE_STEP_2:
            name = "E2";
            break;

        case AVERAGE:
            name = "average";
            break;

        case SLICE:
            name = "slice";
            break;

        case CONTRAST:
            name = "constrast";
            break;

        case PHASE:
            name = "phase";
            break;

        case REPETITION:
            name = "repetition";
            break;

        case SET:
            name = "set";
            break;

        case SEGMENT:
            name = "segment";
            break;

        case USER_0:
            name = "user_0";
            break;

        case USER_1:
            name = "user_1";
            break;

        case USER_2:
            name = "user_2";
            break;

        case USER_3:
            name = "user_3";
            break;

        case USER_4:
            name = "user_4";
            break;

        case USER_5:
            name = "user_5";
            break;

        case USER_6:
            name = "user_6";
            break;

        case USER_7:
            name = "user_7";
            break;

        case NONE:
            name = "none";
            break;

        default:
            GERROR_STREAM("Unrecognized ISMRMRD dimension name : " << v);
        }

        return name;
    }

    IsmrmrdCONDITION get_ismrmrd_dimension(const std::string& name)
    {
        IsmrmrdCONDITION v;

        if (name == "E1")
        {
            v = KSPACE_ENCODE_STEP_1;
        }
        else if (name == "E2")
        {
            v = KSPACE_ENCODE_STEP_2;
        }
        else if (name == "average")
        {
            v = AVERAGE;
        }
        else if (name == "slice")
        {
            v = SLICE;
        }
        else if (name == "contrast")
        {
            v = CONTRAST;
        }
        else if (name == "phase")
        {
            v = PHASE;
        }
        else if (name == "repetition")
        {
            v = REPETITION;
        }
        else if (name == "set")
        {
            v = SET;
        }
        else if (name == "segment")
        {
            v = SEGMENT;
        }
        else if (name == "user_0")
        {
            v = USER_0;
        }
        else if (name == "user_1")
        {
            v = USER_1;
        }
        else if (name == "user_2")
        {
            v = USER_2;
        }
        else if (name == "user_3")
        {
            v = USER_3;
        }
        else if (name == "user_4")
        {
            v = USER_4;
        }
        else if (name == "user_5")
        {
            v = USER_5;
        }
        else if (name == "user_6")
        {
            v = USER_6;
        }
        else if (name == "user_7")
        {
            v = USER_7;
        }
        else if (name == "none")
        {
            v = NONE;
        }
        else
        {
            GERROR_STREAM("Unrecognized ISMRMRD dimension name : " << name);
        }

        return v;
    }

    // --------------------------------------------------------------------

    std::string get_ismrmrd_calib_mode_name(ismrmrdCALIBMODE v)
    {
        std::string name;

        switch (v)
        {
        case ISMRMRD_embedded:
            name = "embedded";
            break;

        case ISMRMRD_interleaved:
            name = "interleaved";
            break;

        case ISMRMRD_separate:
            name = "seperate";
            break;

        case ISMRMRD_external:
            name = "external";
            break;

        case ISMRMRD_other:
            name = "other";
            break;

        case ISMRMRD_noacceleration:
            name = "noacceleration";
            break;

        default:
            GERROR_STREAM("Unrecognized ISMRMRD calibration mode type : " << v);
        }

        return name;
    }

    ismrmrdCALIBMODE get_ismrmrd_calib_mode(const std::string& name)
    {
        ismrmrdCALIBMODE v;

        if (name == "embedded")
        {
            v = ISMRMRD_embedded;
        }
        else if (name == "interleaved")
        {
            v = ISMRMRD_interleaved;
        }
        else if (name == "seperate")
        {
            v = ISMRMRD_separate;
        }
        else if (name == "external")
        {
            v = ISMRMRD_external;
        }
        else if (name == "other")
        {
            v = ISMRMRD_other;
        }
        else if (name == "noacceleration")
        {
            v = ISMRMRD_noacceleration;
        }
        else
        {
            GERROR_STREAM("Unrecognized ISMRMRD calibration mode name : " << name);
        }

        return v;
    }

    // --------------------------------------------------------------------

}
