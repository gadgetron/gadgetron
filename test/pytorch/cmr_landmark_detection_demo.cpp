
#include <boost/filesystem.hpp>
#include "python_toolbox.h"
#include "ismrmrd/ismrmrd.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDArray_reductions.h"
#include "ImageIOAnalyze.h"

using namespace Gadgetron;

void prep_cml_ml(char* gt_home_str, std::string& gt_home, std::string& cmr_ml_home, std::string& model_dir)
{
    Gadgetron::initialize_python();

    GDEBUG_STREAM("=============================================================================================");

    if (gt_home_str != NULL) {
        size_t pos = std::string(gt_home_str).rfind("gadgetron");
        gt_home_str[pos - 1] = '\0';

        std::string gadgetron_python_path = std::string(gt_home_str) + "/share/gadgetron/python";
        Gadgetron::add_python_path(gadgetron_python_path);
        gt_home = gadgetron_python_path;
        GDEBUG_STREAM("Set up python path : " << gt_home);

        cmr_ml_home = gt_home + "/" + "cmr_ml";
        Gadgetron::add_python_path(cmr_ml_home);
        GDEBUG_STREAM("Set up python path : " << cmr_ml_home);

        Gadgetron::add_python_path(cmr_ml_home + "/utils");
        GDEBUG_STREAM("Set up python path : " << cmr_ml_home + "/utils");

        model_dir = cmr_ml_home + "/models";
        GDEBUG_STREAM("model_dir : " << model_dir);
    }

    GDEBUG_STREAM("=============================================================================================");
}

int main(int argc, char** argv)
{
    ImageIOAnalyze gt_io;

    std::string gt_home, cmr_ml_home, model_dir;
    char* gt_home_str = std::getenv("GADGETRON_HOME");
    prep_cml_ml(gt_home_str, gt_home, cmr_ml_home, model_dir);

    char* gt_ut_folder = std::getenv("GT_UNITTEST_DIRECTORY");

    std::string gt_ut_data_folder, gt_ut_res_folder;
    if (gt_ut_folder != NULL)
    {
        gt_ut_data_folder = std::string(gt_ut_folder);
        gt_ut_res_folder = std::string(gt_ut_folder) + "/result/";
        GDEBUG_STREAM("gt_ut_data_folder is " << gt_ut_data_folder);
        GDEBUG_STREAM("gt_ut_res_folder is " << gt_ut_res_folder);
    }

    // CH4
    if (gt_home_str != NULL)
    {
        GILLock lg;

        GDEBUG_STREAM(" --------------------------------------------------------------------------------------------------");
        GDEBUG_STREAM("Test for cine CH4 landmark detection, CMR ML");

        // load model
        PythonFunction<boost::python::object> load_model_cmr_landmark_detection("gadgetron_cmr_landmark_detection", "load_model_cmr_landmark_detection");
        boost::python::object model = load_model_cmr_landmark_detection(model_dir, "CMR_landmark_network_RO_352_E1_352_ch2_ch3_ch4_myo_pts_with_T1_LGE_LossMultiSoftProb_KLD_Dice_Pytorch_1.5.0_2020-06-17_20200617_111642.pts");
        bp::incref(model.ptr());

        hoNDArray<float> ED, ES;

        if (!gt_ut_data_folder.empty())
        {
            gt_io.import_array(ED, gt_ut_data_folder + "cmr_landmark_detection/RetroCine/CH4/20180104_1462193_ch4_ED");
            ED.print(std::cout);
            GDEBUG_STREAM("ED cine image = " << Gadgetron::nrm2(ED));

            gt_io.import_array(ES, gt_ut_data_folder + "cmr_landmark_detection/RetroCine/CH4/20180104_1462193_ch4_ES");
            ES.print(std::cout);
            GDEBUG_STREAM("ES cine image = " << Gadgetron::nrm2(ES));
        }
        else
        {
            ED.create(352, 352);
            ES.create(352, 352);
            Gadgetron::fill(ED, float(0));
            Gadgetron::fill(ES, float(0));
        }

        size_t RO = ED.get_size(0);
        size_t E1 = ES.get_size(1);

        hoNDArray<float> im;
        im.create(RO, E1, 2);
        memcpy(&im(0, 0, 0), ED.begin(), ED.get_number_of_bytes());
        memcpy(&im(0, 0, 1), ES.begin(), ES.get_number_of_bytes());

        boost::python::object* pModel = &model;

        PythonFunction< hoNDArray<float>, hoNDArray<float> > perform_cmr_landmark_detection("gadgetron_cmr_landmark_detection", "perform_cmr_landmark_detection");

        hoNDArray<float> pts, probs;
        std::tie(pts, probs) = perform_cmr_landmark_detection(im, *pModel);

        if (!gt_ut_res_folder.empty())
        {
            std::ostringstream ostr;
            ostr << "20180104_1462193_ch4_pts";
            GDEBUG_STREAM("save " << gt_ut_res_folder + ostr.str());
            gt_io.export_array(pts, gt_ut_res_folder + ostr.str());
        }

        if (!gt_ut_res_folder.empty())
        {
            std::ostringstream ostr;
            ostr << "20180104_1462193_ch4_probs";
            GDEBUG_STREAM("save " << gt_ut_res_folder + ostr.str());
            gt_io.export_array(probs, gt_ut_res_folder + ostr.str());
        }
    }

    return 0;
}
