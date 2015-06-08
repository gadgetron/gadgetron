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

// -----------------------------------------------

template <typename T> 
void findDataRange(const std::vector<hoNDArray<T> >& x, const std::vector<hoNDArray<T> >& y, T& minX, T& maxX, T& minY, T& maxY)
{
    size_t n;

    maxY = Gadgetron::max(const_cast<hoNDArray<T>*>(&y[0]));
    minY = Gadgetron::min(const_cast<hoNDArray<T>*>(&y[0]));

    maxX = Gadgetron::max(const_cast<hoNDArray<T>*>(&x[0]));
    minX = Gadgetron::min(const_cast<hoNDArray<T>*>(&x[0]));

    for (n = 1; n < x.size(); n++)
    {
        T v = Gadgetron::max(const_cast<hoNDArray<T>*>(&y[n]));
        if (v>maxY) maxY = v;

        v = Gadgetron::min(const_cast<hoNDArray<T>*>(&y[n]));
        if (v<minY) minY = v;

        v = Gadgetron::max(const_cast<hoNDArray<T>*>(&x[n]));
        if (v>maxX) maxX = v;

        v = Gadgetron::min(const_cast<hoNDArray<T>*>(&x[n]));
        if (v<minX) minX = v;
    }
}

//0    black(default background)
//1    red(default foreground)
//2    yellow
//3    green
//4    aquamarine
//5    pink
//6    wheat
//7    grey
//8    brown
//9    blue
//10   BlueViolet
//11   cyan
//12   turquoise
//13   magenta
//14   salmon
//15   white

void getPlotColor(size_t n, int& color)
{
    switch (n%15)
    {
        case 0:
            color = 1;
            break;

        case 1:
            color = 2;
            break;

        case 2:
            color = 3;
            break;

        case 3:
            color = 4;
            break;

        case 5:
            color = 6;
            break;

        case 6:
            color = 7;
            break;

        case 7:
            color = 8;
            break;

        case 8:
            color = 9;
            break;

        case 9:
            color = 10;
            break;

        case 10:
            color = 11;
            break;

        case 11:
            color = 12;
            break;

        case 12:
            color = 13;
            break;

        case 13:
            color = 14;
            break;

        case 14:
            color = 15;
            break;

        default:
            color = 15;
    }
}

void getPlotGlyph(size_t n, std::string& gly)
{
    switch (n%15)
    {
        case 0:
            gly = "#(718)"; // o
            break;

        case 1:
            gly = "#(728)"; // *
            break;

        case 2:
            gly = "#(725)"; // +
            break;

        case 3:
            gly = "#(727)"; // x
            break;

        case 5:
            gly = "#(729)"; // .
            break;

        case 6:
            gly = "#(852)"; // triangle (up)
            break;

        case 7:
            gly = "#(854)"; // triangle (down)
            break;

        case 8:
            gly = "#(853)"; // triangle (left)
            break;

        case 9:
            gly = "#(855)"; // triangle (right)
            break;

        case 10:
            gly = "#(857)"; // flag
            break;

        case 11:
            gly = "#(212)"; // :
            break;

        case 12:
            gly = "#(221)"; // <
            break;

        case 13:
            gly = "#(222)"; // >
            break;

        case 14:
            gly = "#(851)"; // square
            break;

        default:
            gly = "#(856)"; // star
    }
}

template <typename T> EXPORTGTPLPLOT
bool plotCurves(const std::vector<hoNDArray<T> >& x, const std::vector<hoNDArray<T> >& y, 
                const std::string& xlabel, const std::string& ylabel, 
                const std::string& title, const std::vector<std::string>& legend, 
                size_t xsize, size_t ysize, 
                bool trueColor, bool drawLine, 
                hoNDArray<float>& plotIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(x.size()>0);
        GADGET_CHECK_RETURN_FALSE(y.size()>0);
        GADGET_CHECK_RETURN_FALSE(x.size() == y.size());

        T minX, maxX, minY, maxY;
        findDataRange(x, y, minX, maxX, minY, maxY);

        plsdev("mem");

        hoNDArray<unsigned char> im;
        im.create(3, xsize, ysize);
        Gadgetron::clear(im);

        plsmem(im.get_size(1), im.get_size(2), im.begin());

        plinit();
        plfont(2);

        pladv(0);
        plvpor(0.15, 0.75, 0.1, 0.9);

        T spaceX = 0.01*(maxX - minX);
        T spaceY = 0.05*(maxY - minY);

        plwind(minX - spaceX, maxX + spaceX, minY - spaceY, maxY + spaceY);

        plcol0(15);
        plbox("bgcnst", 0.0, 0, "bgcnstv", 0.0, 0);

        int mark[2], space[2];

        mark[0] = 4000;
        space[0] = 2500;
        plstyl(1, mark, space);

        size_t num = x.size();

        size_t n;

        hoNDArray<double> xd, yd;

        // draw lines
        for (n = 0; n < num; n++)
        {
            size_t N = y[n].get_size(0);

            xd.copyFrom(x[n]);
            yd.copyFrom(y[n]);

            if (drawLine)
            {
                int c;
                getPlotColor(n, c);
                plcol0(c);
                plline(N, xd.begin(), yd.begin());
            }

            std::string gly;
            getPlotGlyph(n, gly);
            plstring(N, xd.begin(), yd.begin(), gly.c_str());
        }

        plcol0(15);
        plmtex("b", 3.2, 0.5, 0.5, xlabel.c_str());
        plmtex("t", 2.0, 0.5, 0.5, title.c_str());
        plmtex("l", 5.0, 0.5, 0.5, ylabel.c_str());

        // draw the legend
        if (legend.size() == x.size())
        {
            std::vector<PLINT> opt_array(num), text_colors(num), line_colors(num), line_styles(num), symbol_numbers(num), symbol_colors(num);
            std::vector<PLFLT> symbol_scales(num), line_widths(num), box_scales(num, 1);

            std::vector<std::string> glyphs(num);
            std::vector<const char*> symbols(num);
            PLFLT legend_width, legend_height;

            std::vector<const char*> legend_text(num);

            for (n = 0; n < num; n++)
            {
                int c;
                getPlotColor(n, c);
                getPlotGlyph(n, glyphs[n]);

                opt_array[n] = PL_LEGEND_SYMBOL;
                text_colors[n] = 15;
                line_colors[n] = c;
                line_styles[n] = 1;
                line_widths[n] = 1.;
                symbol_colors[n] = c;
                symbol_scales[n] = 1.;
                symbol_numbers[n] = 1;
                symbols[n] = glyphs[n].c_str();
                legend_text[n] = legend[n].c_str();
            }

            // plscol0a(15, 32, 32, 32, 0);
            pllegend(&legend_width, &legend_height,
                PL_LEGEND_BACKGROUND,
                PL_POSITION_OUTSIDE | PL_POSITION_RIGHT | PL_POSITION_TOP,
                0.02, 0.0, 0.15, 
                0, 15, 1, 0, 0,
                num, &opt_array[0],
                0.0, 0.5, 1.0, 1, &text_colors[0], (const char **)(&legend_text[0]),
                NULL, NULL, &box_scales[0], NULL,
                &line_colors[0], &line_styles[0], &line_widths[0],
                &symbol_colors[0], &symbol_scales[0], &symbol_numbers[0], (const char **)(&symbols[0]));
        }

        plend();


        // output plots
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
        GERROR_STREAM("Errors happened in plotCurves(...) ... ");
        return false;
    }

    return true;
}

template EXPORTGTPLPLOT bool plotCurves(const std::vector<hoNDArray<float> >& x, const std::vector<hoNDArray<float> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend,
    size_t xsize, size_t ysize, bool trueColor, bool drawLine, hoNDArray<float>& plotIm);

template EXPORTGTPLPLOT bool plotCurves(const std::vector<hoNDArray<double> >& x, const std::vector<hoNDArray<double> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend,
    size_t xsize, size_t ysize, bool trueColor, bool drawLine, hoNDArray<float>& plotIm);

}
