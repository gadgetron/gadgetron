/** \file       GtPLplot.cpp
    \brief      Implement Gt PLplot functions.
    \author     Hui Xue
*/

#include "GtPLplot.h"
#include <iostream>
#include <stack>
#include <cmath>

namespace Gadgetron
{

template <typename T> EXPORTGTPLPLOT
bool plotCurve(const hoNDArray<T>& x, const hoNDArray<T>& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, size_t xsize, size_t ysize, bool trueColor, hoNDArray<float>& plotIm)
{
    try
    {
        size_t N = y.get_size(0);

        T maxV = Gadgetron::max(const_cast<hoNDArray<T>*>(&y));
        T minV = Gadgetron::min(const_cast<hoNDArray<T>*>(&y));

        T maxX = Gadgetron::max( const_cast<hoNDArray<T>*>(&x) );
        T minX = Gadgetron::min( const_cast<hoNDArray<T>*>(&x) );

        plsdev("mem");

        hoNDArray<unsigned char> im;
        im.create(3, xsize, ysize);
        Gadgetron::clear(im);

        plsmem(im.get_size(1), im.get_size(2), im.begin());

        plinit();
        plfont(2);

        pladv(0);
        plvpor(0.15, 0.85, 0.1, 0.9);

        T spaceX = 0.01*(maxX - minX);
        T spaceY = 0.05*(maxV - minV);

        plwind(minX - spaceX, maxX + spaceX, minV - spaceY, maxV + spaceY);

        plcol0(15);
        plbox("bcnst", 0.0, 0, "bcnstv", 0.0, 0);

        int mark[2], space[2];

        mark[0] = 4000;
        space[0] = 2500;
        plstyl(1, mark, space);

        hoNDArray<double> xd;
        xd.copyFrom(x);

        hoNDArray<double> yd;
        yd.copyFrom(y);

        plcol0(2);
        plline(N, xd.begin(), yd.begin());
        plstring(N, xd.begin(), yd.begin(), "#(728)");

        plcol0(15);
        plmtex("b", 3.2, 0.5, 0.5, xlabel.c_str());
        plmtex("t", 2.0, 0.5, 0.5, title.c_str());
        plmtex("l", 5.0, 0.5, 0.5, ylabel.c_str());

        plend();

        plotIm.copyFrom(im);

        if (trueColor)
        {
            std::vector<size_t> dim_order(3);
            dim_order[0] = 1;
            dim_order[1] = 2;
            dim_order[2] = 0;

            hoNDArray<float> plotImPermuted;
            plotImPermuted.copyFrom(plotIm);

            plotIm.create(xsize, ysize, 3);

            Gadgetron::permute(&plotImPermuted, &plotIm, &dim_order);
        }
        else
        {
            hoNDArray<float> plotIm2D;
            Gadgetron::sum_over_dimension(plotIm, plotIm2D, 0);

            plotIm2D.squeeze();

            std::vector<size_t> dim_order(2);
            dim_order[0] = 0;
            dim_order[1] = 1;

            plotIm.create(xsize, ysize);
            Gadgetron::permute(&plotIm2D, &plotIm, &dim_order);
        }
    }
    catch (...)
    {
        GERROR_STREAM("Errors happened in plotCurve(...) ... ");
        return false;
    }

    return true;
}

template EXPORTGTPLPLOT bool plotCurve(const hoNDArray<float>& x, const hoNDArray<float>& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, size_t xsize, size_t ysize, bool trueColor, hoNDArray<float>& plotIm);
template EXPORTGTPLPLOT bool plotCurve(const hoNDArray<double>& x, const hoNDArray<double>& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, size_t xsize, size_t ysize, bool trueColor, hoNDArray<float>& plotIm);

}
