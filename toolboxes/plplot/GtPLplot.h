/** \file       GtPLplot.h
    \brief      PLplot functions to draw gadgetron data
    \author     Hui Xue
*/

#ifndef GT_PLPLOT_H
#define GT_PLPLOT_H

#include "hoNDArray.h"
#include "hoNDImage.h"

#include "PLplotExport.h"

#include "hoNDArray_elemwise.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_reductions.h"

#include "plConfig.h"
#include "plplot.h"

namespace Gadgetron
{
    // plot a 1D hoNDArray as a curve
    template <typename T> EXPORTGTPLPLOT
    bool plotCurve(const hoNDArray<T>& x, const hoNDArray<T>& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, size_t xsize, size_t ysize, bool trueColor, hoNDArray<float>& plotIm);
}

#endif // GT_PLPLOT_H
