/**
\file       morphology.h
\brief      Morphological image processing related functions
\author     Hui Xue
*/

#ifndef IMAGE_MORPHOLOGY_H
#define IMAGE_MORPHOLOGY_H

#include "hoNDArray.h"
#include "ho2DArray.h"
#include "hoNDImage.h"

#include "hoNDArray_elemwise.h"
#include "hoNDArray_utils.h"

namespace Gadgetron
{
    /// perfrom 2D region growing
    /// input: a 2D array, with object pixels equal to object_value
    /// labelArray: 0 is background
    /// x, y: the seed indices
    /// label: the label of growed region in labelArray, >0
    /// is_8_connected: whether to label 8-connected objects; if false, 4-connected objects are labelled
    template <typename T>
    void region_growing_2d(const hoNDArray<T>& input, T object_value, hoNDArray<unsigned int>& label_array, size_t x, size_t y, unsigned int label, bool is_8_connected);

    /// perfrom connected component labelling
    /// input: a 2D array, with object pixels equal to object_value
    /// label: connected component label matrix, 0 is background
    template <typename T>
    void bwlabel_2d(const hoNDArray<T>& input, T object_value, hoNDArray<unsigned int>& label, bool is_8_connected);

    /// for the labelled array, find all regions and their areas
    void bwlabel_area_2d(const hoNDArray<unsigned int>& label_array, std::vector<unsigned int>& labels, std::vector<unsigned int>& areas);

    /// clean foreground and background using bwlabel
    template <typename T>
    void bwlabel_clean_fore_and_background(const hoNDArray<T>& input, T object_value, T bg_value, size_t obj_thres, size_t bg_size, bool is_8_connected, hoNDArray<T>& output);
}

#endif // IMAGE_MORPHOLOGY_H
