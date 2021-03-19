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
    #define GADGETRON_IMAGE_SATURATIONTIME                 "GADGETRON_TS"

    /// role of image data
    #define GADGETRON_DATA_ROLE                            "GADGETRON_DataRole"
    #define GADGETRON_IMAGE_REGULAR                        "Image"
    #define GADGETRON_IMAGE_RETRO                          "ImageRetro"
    #define GADGETRON_IMAGE_MOCORECON                      "ImageMoCo"
    #define GADGETRON_IMAGE_GFACTOR                        "Gfactor"
    #define GADGETRON_IMAGE_SNR_MAP                        "SNR_MAP"
    #define GADGETRON_IMAGE_STD_MAP                        "STD_MAP"
    #define GADGETRON_IMAGE_WRAPAROUNDMAP                  "WrapAround_MAP"
    #define GADGETRON_IMAGE_PHASE                          "Phase"
    #define GADGETRON_IMAGE_INTENSITY_UNCHANGED            "Image_Intensity_Unchanged"
    #define GADGETRON_IMAGE_AIF                            "AIF"
    #define GADGETRON_IMAGE_AIF_LV_MASK                    "AIFMASK"
    #define GADGETRON_IMAGE_AIF_Gd_CONCENTRATION           "Gd"
    #define GADGETRON_IMAGE_PERF_FLOW_MAP                  "Flow_Map"
    #define GADGETRON_IMAGE_PERF_MAP                       "Perf_Map"
    #define GADGETRON_IMAGE_PERF_MEANTRANSITTIME_MAP       "MTT_Map"
    #define GADGETRON_IMAGE_PERF_INTERSITITAL_VOLUME_MAP   "Interstitial_Volume_Map"
    #define GADGETRON_IMAGE_PERF_VASCULAR_VOLUME_MAP       "Vascular_Volume_Map"
    #define GADGETRON_IMAGE_PERF_Gd_Extraction_MAP         "Gd_Extraction_Map"
    #define GADGETRON_IMAGE_PERF_PERMEABILITY_SURFACE_AREA_MAP "PS_Map"
    #define GADGETRON_IMAGE_PERF_Gd_CONCENTRATION          "Gd"
    #define GADGETRON_IMAGE_PERF_AHA_SEGMENT_MODEL         "AHA"

    // other images than the regular reconstruction results
    #define GADGETRON_IMAGE_OTHER                          "Image_Other"
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
    #define GADGETRON_IMAGE_PERFUSIONMAP                   "PERFUSIONMAP"
    #define GADGETRON_IMAGE_MEANTRANSITTIMEMAP             "MTT"
    #define GADGETRON_IMAGE_INTERSTITIALVOLUMEMAP          "INTERSTITIALVOLUME"
    #define GADGETRON_IMAGE_VASCULARVOLUMEMAP              "VASCULARVOLUME"
    #define GADGETRON_IMAGE_GD_EXTRACTIONMAP               "EXTRACTIONMAP"
    #define GADGETRON_IMAGE_PERMEABILITY_SURFACE_AREAMAP   "PSMAP"
    #define GADGETRON_IMAGE_RADIAL_STRAINMAP               "RADIAL_STRAINMAP"
    #define GADGETRON_IMAGE_CIRCUM_STRAINMAP               "CIRCUM_STRAINMAP"
    #define GADGETRON_IMAGE_MAX_STRAINMAP                  "MAX_STRAINMAP"
    #define GADGETRON_IMAGE_ACTIVATIONTIME_STRAINMAP       "TA_STRAINMAP"

    /// image annotation tags
    #define GADGETRON_CMR_2D_ENDO_CONTOUR                   "ENDO"
    #define GADGETRON_CMR_2D_EPI_CONTOUR                    "EPI"

    /// 2D closed region of interest contour
    /// stored as 2*N_pts + 4 array
    /// first three values: ROI color in tuple
    /// the forth value: line thickness
    /// then [px1 py1 px2 py2 ...] for N_Pts points
    #define GADGETRON_2D_ROI                                "GT_ROI"

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
    #define GADGETRON_IMAGE_RECON_FIGURE                    "FIG"

    /// calculation comments
    #define GADGETRON_SUBTRACTION                           "SUB"
    #define GADGETRON_AI                                    "AI"

    /// control flags for image recon and other processing
    /// if set, skip the processing steps after the recon, e.g. partial fourier handling and kspace filter
    #define GADGETRON_SKIP_PROCESSING_AFTER_RECON           "Skip_processing_after_recon"
    #define GADGETRON_USE_DEDICATED_SCALING_FACTOR          "Use_dedicated_scaling_factor"

    /// if set, auto correct image orientation by MR convention
    #define GADGETRON_CORRECT_IMAGE_ORIENTATION             "Correct_image_orientation"

    /// instruct the client to keep image geometry
    #define GADGETRON_KEEP_IMAGE_GEOMETRY                   "Keep_image_geometry"

    /// instruct the client to send images to database without further processing
    #define GADGETRON_DIRECT_IMAGE_SEND                     "Gadgetron_ImageDirectSend"

    /// data flow tag
    /// if this flag is set to be 1 for a image, the image is immediately passed to the next gadget
    /// if this flag is 0, this image is a stored image by the accummulator
    /// whether to pass a stored image to the next gadget is determined by the processing gadget itself
    #define GADGETRON_PASS_IMMEDIATE                       "GT_PASSIMAGE_IMMEDIATE"
}
