/** \file       cmr_image_container_util.h
    \brief      Util functions for image and image containers
    \author     Hui Xue
*/

#pragma once

#define NOMINMAX

#include <cstdio>
#include <complex>

#include "cmr_export.h"

#include "hoNDImage_util.h"
#include "hoMRImage.h"
#include "hoNDImageContainer2D.h"

namespace Gadgetron
{
    /// resample image container
    template <typename T, unsigned int D>
    EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<T, D> >& input, const std::vector<size_t>& dim_out, hoNDImageContainer2D<hoMRImage<T, D> >& output);

    template <typename T, unsigned int D>
    EXPORTCMR void resample_image_container(const hoNDImageContainer2D<hoMRImage<T, D> >& input, double scale_ratio, hoNDImageContainer2D<hoMRImage<T, D> >& output);

    /// convert a 4D array to image container
    /// input: [RO E1 PHS SLC]
    /// output: image array with SLC rows and PHS columns
    template <typename T>
    EXPORTCMR void convert_4D_array_to_container(const hoNDArray<T>& input, hoNDImageContainer2D < hoMRImage<T, 2> > & output);

    template <typename T>
    EXPORTCMR void convert_container_to_4D_array(const hoNDImageContainer2D < hoMRImage<T, 2> > & input, hoNDArray<T>& output);

    /// resort image container and 4D array
    /// input : [RO E1 PHS SLC], resort along SLC
    template <typename T>
    EXPORTCMR void sort_4D_array(const hoNDArray<T>& input, const std::vector<size_t> basal_to_apex_slices, hoNDArray<T>& output);

    template <typename T>
    EXPORTCMR void sort_image_container(const hoNDImageContainer2D < hoMRImage<T, 2> > & input, const std::vector<size_t> basal_to_apex_slices, hoNDImageContainer2D < hoMRImage<T, 2> > & output);

    /// serialize meta container for all images in the container
    template <typename ImageType>
    EXPORTCMR void serialize_contrainer_meta_attributes(hoNDImageContainer2D<ImageType>& input, std::vector<std::vector<std::string> >& attribs);

    /// write attributes to hard drive
    EXPORTCMR void write_contrainer_meta_attributes(const std::vector<std::vector<std::string> >& attribs, const std::string& prefix);

}
