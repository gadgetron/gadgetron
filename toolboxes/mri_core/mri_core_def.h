/** \file   mri_core_def.h
    \brief  Define the symbols for mri_core toolbox
    \author Hui Xue
*/

#pragma once

#include "ismrmrd/ismrmrd.h"

namespace Gadgetron
{
    // define the dimensions of ISMRMRD
    enum ISMRMRDDIM
    {
        DIM_ReadOut = 32,
        DIM_Encoding1,
        DIM_Channel,
        DIM_Slice,
        DIM_Encoding2,
        DIM_Contrast,
        DIM_Phase,
        DIM_Repetition,
        DIM_Set,
        DIM_Segment,
        DIM_Average,
        DIM_other1,
        DIM_other2,
        DIM_other3,
        DIM_NONE
    };

    // define the reconstruction algorithms
    enum ISMRMRDALGO
    {
        ISMRMRD_GRAPPA = 64,
        ISMRMRD_SENSE,
        ISMRMRD_SPIRIT,
        ISMRMRD_L1SPIRIT,
        ISMRMRD_SOFTSENSE,
        ISMRMRD_L1SOFTSENSE,
        ISMRMRD_2DTBINNING,
        ISMRMRD_2DTBINNING_FLOW,
        ISMRMRD_L1SPIRIT_SLEP,
        ISMRMRD_L1SPIRIT_SLEP_MOTION_COMP,
        ISMRMRD_NONE
    };

    // define the coil sensitivity map estimation algorithms
    enum ISMRMRDCOILMAPALGO
    {
        ISMRMRD_SOUHEIL = 96,
        ISMRMRD_SOUHEIL_ITER
    };

    // define the partial fourier/asymmetric echo handling algorithms
    enum ISMRMRDPFALGO
    {
        ISMRMRD_PF_HOMODYNE = 128,          // iterative homodyne
        ISMRMRD_PF_POCS,                    // POCS
        ISMRMRD_PF_FENGHUANG,               // convolution based method
        ISMRMRD_PF_ZEROFILLING_FILTER,      // zero-filling with partial fourier filter
        ISMRMRD_PF_ZEROFILLING,             // zero-filling without partial fourier filter
        ISMRMRD_PF_NONE
    };

    // define the kspace filter type
    enum ISMRMRDKSPACEFILTER
    {
        ISMRMRD_FILTER_GAUSSIAN = 160,
        ISMRMRD_FILTER_HANNING,
        ISMRMRD_FILTER_TUKEY,
        ISMRMRD_FILTER_TAPERED_HANNING,
        ISMRMRD_FILTER_NONE
    };

    // define the calibration mode of ISMRMRD
    enum ISMRMRDCALIBMODE
    {
        ISMRMRD_embedded = 256,
        ISMRMRD_interleaved,
        ISMRMRD_separate,
        ISMRMRD_external,
        ISMRMRD_other,
        ISMRMRD_noacceleration
    };

    // define the interpolation method
    enum ISMRMRDINTERP
    {
        ISMRMRD_INTERP_LINEAR = 512,
        ISMRMRD_INTERP_SPLINE,
        ISMRMRD_INTERP_BSPLINE
    };

    // define the interpolation method for retro-gating
    enum ISMRMRDINTERPRETROGATING
    {
        ISMRMRD_INTERP_RETRO_GATING_LINEAR = 600,
        ISMRMRD_INTERP_RETRO_GATING_CUBIC, 
        ISMRMRD_INTERP_RETRO_GATING_BSPLINE
    };

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

    /// data flow tag
    /// if this flag is set to be 1 for a image, the image is immediately passed to the next gadget
    /// if this flag is 0, this image is a stored image by the accummulator
    /// whether to pass a stored image to the next gadget is determined by the processing gadget itself
    #define GADGETRON_PASS_IMMEDIATE                       "GT_PASSIMAGE_IMMEDIATE"

    /// data processing tag, used with ImageProcessingHistory
    #define GADGETRON_IMAGE_SURFACECOILCORRECTION           "NORM"
    #define GADGETRON_IMAGE_FILTER                          "FIL"
    #define GADGETRON_IMAGE_MOCO                            "MOCO"
    #define GADGETRON_IMAGE_AVE                             "AVE"

    /// dimension string
    #define GADGETRON_RO                                    "RO"
    #define GADGETRON_E1                                    "E1"
    #define GADGETRON_CHA                                   "CHA"
    #define GADGETRON_SLC                                   "SLC"
    #define GADGETRON_E2                                    "E2"
    #define GADGETRON_CONTRAST                              "CON"
    #define GADGETRON_PHASE                                 "PHS"
    #define GADGETRON_REP                                   "REP"
    #define GADGETRON_SET                                   "SET"
    #define GADGETRON_SEGMENT                               "SEG"
    #define GADGETRON_AVERAGE                               "AVE"
    #define GADGETRON_OTHER1                                "OTH1"
    #define GADGETRON_OTHER2                                "OTH2"
    #define GADGETRON_OTHER3                                "OTH3"
    #define GADGETRON_NONE                                  "NONE"

    /// ISMRMRD Image fields
    #define ISMRMRD_IMAGE_version                       "ISMRMRD_IMAGE_version"
    #define ISMRMRD_IMAGE_flags                         "ISMRMRD_IMAGE_flags"
    #define ISMRMRD_IMAGE_measurement_uid               "ISMRMRD_IMAGE_measurement_uid"
    #define ISMRMRD_IMAGE_matrix_size                   "ISMRMRD_IMAGE_matrix_size"
    #define ISMRMRD_IMAGE_field_of_view                 "ISMRMRD_IMAGE_field_of_view"
    #define ISMRMRD_IMAGE_channels                      "ISMRMRD_IMAGE_channels"
    #define ISMRMRD_IMAGE_position                      "ISMRMRD_IMAGE_position"
    #define ISMRMRD_IMAGE_read_dir                      "ISMRMRD_IMAGE_read_dir"
    #define ISMRMRD_IMAGE_phase_dir                     "ISMRMRD_IMAGE_phase_dir"
    #define ISMRMRD_IMAGE_slice_dir                     "ISMRMRD_IMAGE_slice_dir"
    #define ISMRMRD_IMAGE_patient_table_position        "ISMRMRD_IMAGE_patient_table_position"
    #define ISMRMRD_IMAGE_average                       "ISMRMRD_IMAGE_average"
    #define ISMRMRD_IMAGE_slice                         "ISMRMRD_IMAGE_slice"
    #define ISMRMRD_IMAGE_contrast                      "ISMRMRD_IMAGE_contrast"
    #define ISMRMRD_IMAGE_phase                         "ISMRMRD_IMAGE_phase"
    #define ISMRMRD_IMAGE_repetition                    "ISMRMRD_IMAGE_repetition"
    #define ISMRMRD_IMAGE_set                           "ISMRMRD_IMAGE_set"
    #define ISMRMRD_IMAGE_acquisition_time_stamp        "ISMRMRD_IMAGE_acquisition_time_stamp"
    #define ISMRMRD_IMAGE_physiology_time_stamp         "ISMRMRD_IMAGE_physiology_time_stamp"
    #define ISMRMRD_IMAGE_image_data_type               "ISMRMRD_IMAGE_image_data_type"
    #define ISMRMRD_IMAGE_image_type                    "ISMRMRD_IMAGE_image_type"
    #define ISMRMRD_IMAGE_image_index                   "ISMRMRD_IMAGE_image_index"
    #define ISMRMRD_IMAGE_image_series_index            "ISMRMRD_IMAGE_image_series_index"
    #define ISMRMRD_IMAGE_user_int                      "ISMRMRD_IMAGE_user_int"
    #define ISMRMRD_IMAGE_user_float                    "ISMRMRD_IMAGE_user_float"
}
