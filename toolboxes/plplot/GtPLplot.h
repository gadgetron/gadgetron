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

namespace Gadgetron
{
    /// plot 1D hoNDArray as curves
    /// x: x coordinate of 1D curve
    /// y: values of 1D curve
    /// xlabel, ylabel: x and y label strings
    /// title: title string
    /// legend: legend string for every curve, if empty, no legends will be plotted
    /// xsize, ysize: the plotted image size
    /// plotIm: the plotted image, if trueColor = true, [xsize ysize 3] array
    /// if trueColor = false, the grey image will be generated [xsize ysize]
    /// if drawLine = false, the line will be not plotted, onl glyphs
    template <typename T> EXPORTGTPLPLOT
    bool plotCurves(const std::vector<hoNDArray<T> >& x, const std::vector<hoNDArray<T> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend, const std::vector<std::string>& symbols, 
        size_t xsize, size_t ysize, bool trueColor, bool drawLine, const std::vector<int>& lineStyple, const std::vector<int>& lineWidth, hoNDArray<float>& plotIm);

    /// user can supply x/y axis limits to plot
    /// xlim: xmin to xmax
    /// ylim: ymin to ymax
    template <typename T> EXPORTGTPLPLOT
        bool plotCurves(const std::vector<hoNDArray<T> >& x, const std::vector<hoNDArray<T> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend, const std::vector<std::string>& symbols, 
            size_t xsize, size_t ysize, T xlim[2], T ylim[2], bool trueColor, bool drawLine, const std::vector<int>& lineStyple, const std::vector<int>& lineWidth, hoNDArray<float>& plotIm);

    /// plot noise std
    template <typename T> EXPORTGTPLPLOT
    bool plotNoiseStandardDeviation(const hoNDArray< std::complex<T> >& m, const std::vector<std::string>& coilStrings, const std::string& xlabel, const std::string& ylabel, const std::string& title, size_t xsize, size_t ysize, bool trueColor, hoNDArray<float>& plotIm);
}

#endif // GT_PLPLOT_H
