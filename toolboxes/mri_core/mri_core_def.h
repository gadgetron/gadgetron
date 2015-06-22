/** \file   mri_core_def.h
    \brief  Define the symbols for mri_core toolbox
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"

namespace Gadgetron
{
    /// defination of image meta attributes
    /// user can set these attributes to record some properties of generated imaging results
    /// how to interpret these attributes depends on the client side
    #define GADGETRON_IMAGENUMBER                          "GADGETRON_ImageNumber"
    #define GADGETRON_IMAGECOMMENT                         "GADGETRON_ImageComment"
    #define GADGETRON_IMAGEPROCESSINGHISTORY               "GADGETRON_ImageProcessingHistory"
    #define GADGETRON_IMAGE_CATEGORY                       "GADGETRON_ImageCategory"
    #define GADGETRON_SEQUENCEDESCRIPTION                  "GADGETRON_SeqDescription"
    #define GADGETRON_IMAGE_WINDOWCENTER                   "GADGETRON_WindowCenter"
    #define GADGETRON_IMAGE_WINDOWWIDTH                    "GADGETRON_WindowWidth"
    #define GADGETRON_IMAGE_SCALE_RATIO                    "GADGETRON_ScaleRatio"
    #define GADGETRON_IMAGE_SCALE_OFFSET                   "GADGETRON_ScaleOffset"
    #define GADGETRON_IMAGE_COLORMAP                       "GADGETRON_ColorMap"
    #define GADGETRON_IMAGE_ECHOTIME                       "GADGETRON_TE"
    #define GADGETRON_IMAGE_INVERSIONTIME                  "GADGETRON_TI"

    /// role of image data
    #define GADGETRON_DATA_ROLE                            "GADGETRON_DataRole"
    #define GADGETRON_IMAGE_REGULAR                        "GADGETRON_Image"
    #define GADGETRON_IMAGE_RETRO                          "GADGETRON_ImageRetro"
    #define GADGETRON_IMAGE_MOCORECON                      "GADGETRON_ImageMoCo"
    #define GADGETRON_IMAGE_GFACTOR                        "GADGETRON_Gfactor"
    #define GADGETRON_IMAGE_SNR_MAP                        "GADGETRON_SNR_MAP"
    #define GADGETRON_IMAGE_STD_MAP                        "GADGETRON_STD_MAP"
    #define GADGETRON_IMAGE_WRAPAROUNDMAP                  "GADGETRON_WrapAround_MAP"
    #define GADGETRON_IMAGE_PHASE                          "GADGETRON_Phase"
    #define GADGETRON_IMAGE_INTENSITY_UNCHANGED            "GADGETRON_Image_Intensity_Unchanged"
    #define GADGETRON_IMAGE_AIF                            "GADGETRON_AIF"
    #define GADGETRON_IMAGE_AIF_LV_MASK                    "GADGETRON_AIFLVMASK"
    #define GADGETRON_IMAGE_AIF_Gd_CONCENTRATION           "GADGETRON_AIF_Gd_Concentration"
    #define GADGETRON_IMAGE_PERF_FLOW_MAP                  "GADGETRON_Perf_Flow_Map"
    #define GADGETRON_IMAGE_PERF_Gd_CONCENTRATION          "GADGETRON_Perf_Gd_Concentration"

    // other images than the regular reconstruction results
    #define GADGETRON_IMAGE_OTHER                          "GADGETRON_Image_Other"
    // other data roles
    #define GADGETRON_IMAGE_T2W                            "T2W"
    #define GADGETRON_IMAGE_PD                             "PD"
    #define GADGETRON_IMAGE_MAGIR                          "MAGIR"
    #define GADGETRON_IMAGE_PSIR                           "PSIR"

    #define GADGETRON_IMAGE_T1MAP                          "T1"
    #define GADGETRON_IMAGE_T1SDMAP                        "T1SD"
    #define GADGETRON_IMAGE_T2MAP                          "T2"
    #define GADGETRON_IMAGE_T2SDMAP                        "T2SD"
    #define GADGETRON_IMAGE_T2STARMAP                      "T2STAR"
    #define GADGETRON_IMAGE_T2STARMASKMAP                  "T2SMASKMAP"
    #define GADGETRON_IMAGE_T2STARSDMAP                    "T2STARSD"
    #define GADGETRON_IMAGE_T2STARAMAP                     "T2STARAMAP"
    #define GADGETRON_IMAGE_T2STARTRUNCMAP                 "T2STARTRUNCMAP"
    #define GADGETRON_IMAGE_R2STARMAP                      "R2STAR"

    #define GADGETRON_IMAGE_FAT                            "FAT"
    #define GADGETRON_IMAGE_WATER                          "WATER"
    #define GADGETRON_IMAGE_FREQMAP                        "FREQMAP"
    #define GADGETRON_IMAGE_B1MAP                          "B1MAP"
    #define GADGETRON_IMAGE_FLIPANGLEMAP                   "FLIPANGLEMAP"
    #define GADGETRON_IMAGE_FLOWMAP                        "FLOWMAP"

    //MSH: Interventional MRI (Interactive Real Time, IRT)
    #define GADGETRON_IMAGE_IRT_IMAGE                      "IRT_IMAGE"
    #define GADGETRON_IMAGE_IRT_DEVICE                     "IRT_DEVICE"
    #define GADGETRON_IMAGE_NUM_DEVICE_CHA                 "IRT_NUM_DEVICE_CHA"
    #define GADGETRON_IMAGE_CUR_DEVICE_CHA                 "IRT_CUR_DEVICE_CHA"

    /// data processing tag, used with ImageProcessingHistory
    #define GADGETRON_IMAGE_SURFACECOILCORRECTION           "NORM"
    #define GADGETRON_IMAGE_FILTER                          "FIL"
    #define GADGETRON_IMAGE_MOCO                            "MOCO"
    #define GADGETRON_IMAGE_AVE                             "AVE"

    /// figures created during reconstruction
    #define GADGETRON_IMAGE_RECON_FIGURE                   "FIG"
}
