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
    #define GADGETRON_IMAGENUMBER                          "GT_ImageNumber"
    #define GADGETRON_IMAGECOMMENT                         "GT_ImageComment"
    #define GADGETRON_IMAGEPROCESSINGHISTORY               "GT_ImageProcessingHistory"
    #define GADGETRON_IMAGE_CATEGORY                       "GT_ImageCategory"
    #define GADGETRON_SEQUENCEDESCRIPTION                  "GT_SeqDescription"
    #define GADGETRON_IMAGE_WINDOWCENTER                   "GT_WindowCenter"
    #define GADGETRON_IMAGE_WINDOWWIDTH                    "GT_WindowWidth"
    #define GADGETRON_IMAGE_SCALE_RATIO                    "GT_ScaleRatio"
    #define GADGETRON_IMAGE_SCALE_OFFSET                   "GT_ScaleOffset"
    #define GADGETRON_IMAGE_COLORMAP                       "GT_ColorMap"
    #define GADGETRON_IMAGE_ECHOTIME                       "GT_TE"
    #define GADGETRON_IMAGE_INVERSIONTIME                  "GT_TI"

    /// role of image data
    #define GADGETRON_DATA_ROLE                            "GT_DataRole"
    #define GADGETRON_IMAGE_REGULAR                        "GT_Image"
    #define GADGETRON_IMAGE_RETRO                          "GT_ImageRetro"
    #define GADGETRON_IMAGE_MOCORECON                      "GT_ImageMoCo"
    #define GADGETRON_IMAGE_GFACTOR                        "GT_Gfactor"
    #define GADGETRON_IMAGE_SNR_MAP                        "GT_SNR_MAP"
    #define GADGETRON_IMAGE_STD_MAP                        "GT_STD_MAP"
    #define GADGETRON_IMAGE_WRAPAROUNDMAP                  "GT_WrapAround_MAP"
    #define GADGETRON_IMAGE_PHASE                          "GT_Phase"
    #define GADGETRON_IMAGE_INTENSITY_UNCHANGED            "GT_Image_Intensity_Unchanged"
    #define GADGETRON_IMAGE_AIF                            "GT_AIF"
    // other images than the regular reconstruction results
    #define GADGETRON_IMAGE_OTHER                          "GT_Image_Other"
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

    #define GADGETRON_IMAGE_FAT                            "FAT"
    #define GADGETRON_IMAGE_WATER                          "WATER"
    #define GADGETRON_IMAGE_FREQMAP                        "FREQMAP"
    #define GADGETRON_IMAGE_B1MAP                          "B1MAP"
    #define GADGETRON_IMAGE_FLIPANGLEMAP                   "FLIPANGLEMAP"

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
}
