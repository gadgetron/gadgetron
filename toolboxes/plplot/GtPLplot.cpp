/** \file       GtPLplot.cpp
    \brief      Implement Gt PLplot functions.
    \author     Hui Xue
*/

#include "GtPLplot.h"
#include <iostream>
#include <stack>
#include <cmath>

#include "plConfig.h"
#include "plplot.h"

namespace Gadgetron
{

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

void outputPlotIm(const hoNDArray<unsigned char>& im, bool trueColor, hoNDArray<float>& plotIm)
{
    size_t xsize = im.get_size(1);
    size_t ysize = im.get_size(2);

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

        Gadgetron::permute(plotImPermuted, plotIm, dim_order);
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
        Gadgetron::permute(plotIm2D, plotIm, dim_order);
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
            color = 15;
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
            color = 1;
            break;

        default:
            color = 15;
    }
}

void getPlotGlyph(size_t n, std::string& gly)
{
    switch (n%26)
    {
        case 0:
            gly = "#(840)"; // circle;
            break;

        case 1:
            gly = "#(752)"; // *
            break;

        case 2:
            gly = "#(225)"; // +
            break;

        case 3:
            gly = "#(048)"; // x
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

        case 15:
            gly = "#(751)"; // rectangle
            break;

        case 16:
            gly = "#(828)"; 
            break;

        case 17:
            gly = "#(227)"; // x
            break;

        case 18:
            gly = "#(229)"; // .
            break;

        case 19:
            gly = "#(233)"; // #
            break;

        case 20:
            gly = "#(850)"; // filled circle
            break;

        case 21:
            gly = "#(844)"; // unfilled star
            break;

        case 22:
            gly = "#(843)"; // unfilled diamond
            break;

        case 23:
            gly = "#(842)"; // unfilled triangle
            break;

        case 24:
            gly = "#(841)"; // big square
            break;

        case 25:
            gly = "#(2367)"; // .
            break;

        default:
            gly = "#(856)"; // star
    }
}

template <typename T> EXPORTGTPLPLOT
bool plotCurves(const std::vector<hoNDArray<T> >& x, const std::vector<hoNDArray<T> >& y,
    const std::string& xlabel, const std::string& ylabel,
    const std::string& title, const std::vector<std::string>& legend,
    const std::vector<std::string>& symbols,
    size_t xsize, size_t ysize,
    bool trueColor, bool drawLine, 
    const std::vector<int>& lineStyple, 
    const std::vector<int>& lineWidth,
    hoNDArray<float>& plotIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(x.size()>0);
        GADGET_CHECK_RETURN_FALSE(y.size()>0);
        GADGET_CHECK_RETURN_FALSE(x.size() == y.size());

        T xlim[2], ylim[2];
        findDataRange(x, y, xlim[0], xlim[1], ylim[0], ylim[1]);

        return Gadgetron::plotCurves(x, y, xlabel, ylabel, title, legend, symbols, xsize, ysize, xlim, ylim, trueColor, drawLine, lineStyple, lineWidth, plotIm);
    }
    catch(...)
    {
        GERROR_STREAM("Errors happened in plotCurves(...) ... ");
        return false;
    }

    return true;
}

template EXPORTGTPLPLOT bool plotCurves(const std::vector<hoNDArray<float> >& x, const std::vector<hoNDArray<float> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend, const std::vector<std::string>& symbols, 
    size_t xsize, size_t ysize, bool trueColor, bool drawLine, const std::vector<int>& lineStyple, const std::vector<int>& lineWidth, hoNDArray<float>& plotIm);

template EXPORTGTPLPLOT bool plotCurves(const std::vector<hoNDArray<double> >& x, const std::vector<hoNDArray<double> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend, const std::vector<std::string>& symbols, 
    size_t xsize, size_t ysize, bool trueColor, bool drawLine, const std::vector<int>& lineStyple, const std::vector<int>& lineWidth, hoNDArray<float>& plotIm);

template <typename T> EXPORTGTPLPLOT
bool plotCurves(const std::vector<hoNDArray<T> >& x, const std::vector<hoNDArray<T> >& y, 
                const std::string& xlabel, const std::string& ylabel, 
                const std::string& title, const std::vector<std::string>& legend, 
                const std::vector<std::string>& symbols, 
                size_t xsize, size_t ysize, 
                T xlim[2], T ylim[2], 
                bool trueColor, bool drawLine, 
                const std::vector<int>& lineStyple, 
                const std::vector<int>& lineWidth,
                hoNDArray<float>& plotIm)
{
    try
    {
        GADGET_CHECK_RETURN_FALSE(x.size()>0);
        GADGET_CHECK_RETURN_FALSE(y.size()>0);
        GADGET_CHECK_RETURN_FALSE(x.size() == y.size());

        T minX = xlim[0];
        T maxX = xlim[1];
        T minY = ylim[0];
        T maxY = ylim[1];

        plsdev("mem");

        hoNDArray<unsigned char> im;
        im.create(3, xsize, ysize);
        Gadgetron::clear(im);

        plsmem(im.get_size(1), im.get_size(2), im.begin());

        plinit();
        plfont(2);

        pladv(0);

        if (legend.size() == x.size())
        {
            plvpor(0.11, 0.75, 0.1, 0.9);
        }
        else
        {
            plvpor(0.15, 0.85, 0.1, 0.9);
        }

        /*T spaceX = 0.01*(maxX - minX);
        T spaceY = 0.05*(maxY - minY);
        plwind(minX - spaceX, maxX + spaceX, minY - spaceY, maxY + spaceY);*/

        plvsta();
        plwind(minX, maxX, minY, maxY);
        plsmaj(0.0, 0.0);
        plsmin(0.0, 0.0);

        plcol0(15);
        plbox("bgcnst", 0.0, 0, "bgcnstv", 0.0, 0);

        // int mark[2], space[2];

        //mark[0] = 4000;
        //space[0] = 2500;
        //plstyl(1, mark, space);

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
                if (lineStyple.size() > n)
                {
                    pllsty(lineStyple[n]);
                    plwidth(lineWidth[n]);
                }
                else
                {
                    pllsty(1);
                    plwidth(2);
                }
                plline(N, xd.begin(), yd.begin());
            }

            std::string gly;
            if(symbols.size()>n)
            {
                gly = symbols[n];
                plstring(N, xd.begin(), yd.begin(), gly.c_str());
            }
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
            std::vector<const char*> syms;
            PLFLT legend_width, legend_height;

            std::vector<const char*> legend_text(num);

            for (n = 0; n < num; n++)
            {
                int c;
                getPlotColor(n, c);
                getPlotGlyph(n, glyphs[n]);

                opt_array[n] = PL_LEGEND_SYMBOL | PL_LEGEND_LINE;
                text_colors[n] = 15;
                line_colors[n] = c;

                if (lineStyple.size() > n)
                    line_styles[n] = lineStyple[n];
                else
                    line_styles[n] = 1;

                line_widths[n] = 0.4;
                symbol_colors[n] = c;
                symbol_scales[n] = 0.75;
                symbol_numbers[n] = 1;

                if (symbols.size() > n)
                {    
                    syms.push_back(symbols[n].c_str());
                }
                else
                {
                    syms.push_back("");
                }
                legend_text[n] = legend[n].c_str();
            }

            pllegend(&legend_width, 
                    &legend_height,
                    PL_LEGEND_BACKGROUND,
                    PL_POSITION_OUTSIDE | PL_POSITION_RIGHT | PL_POSITION_TOP,
                    0.02,                                       // x
                    0.0,                                        // y
                    0.05,                                       // plot_width
                    0,                                          // bg_color
                    15,                                         // bb_color
                    1,                                          // bb_style
                    0,                                          // nrow
                    0,                                          // ncolumn
                    num,                                        // nlegend
                    &opt_array[0], 
                    0.05,                                       // text_offset
                    0.35,                                       // text_scale
                    1.0,                                        // text_spacing
                    0.5,                                        // text_justification
                    &text_colors[0], 
                    (const char **)(&legend_text[0]), 
                    NULL,                                       // box_colors
                    NULL,                                       // box_patterns
                    &box_scales[0],                             // box_scales
                    NULL,                                       // box_line_widths
                    &line_colors[0], 
                    &line_styles[0], 
                    &line_widths[0],
                    &symbol_colors[0], 
                    &symbol_scales[0], 
                    &symbol_numbers[0], 
                    (const char **)(&syms[0])
                    );
        }

        plsmin(0.0, 1.0);
        plsmaj(0.0, 1.0);

        plend();

        outputPlotIm(im, trueColor, plotIm);
    }
    catch (...)
    {
        GERROR_STREAM("Errors happened in plotCurves(xlim, ylim) ... ");
        return false;
    }

    return true;
}

template EXPORTGTPLPLOT bool plotCurves(const std::vector<hoNDArray<float> >& x, const std::vector<hoNDArray<float> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend, const std::vector<std::string>& symbols, 
    size_t xsize, size_t ysize, float xlim[2], float ylim[2], bool trueColor, bool drawLine, const std::vector<int>& lineStyple, const std::vector<int>& lineWidth, hoNDArray<float>& plotIm);

template EXPORTGTPLPLOT bool plotCurves(const std::vector<hoNDArray<double> >& x, const std::vector<hoNDArray<double> >& y, const std::string& xlabel, const std::string& ylabel, const std::string& title, const std::vector<std::string>& legend, const std::vector<std::string>& symbols, 
    size_t xsize, size_t ysize, double xlim[2], double ylim[2], bool trueColor, bool drawLine, const std::vector<int>& lineStyple, const std::vector<int>& lineWidth, hoNDArray<float>& plotIm);

// ---------------------------------------------------

template <typename T> 
bool plotNoiseStandardDeviation(const hoNDArray< std::complex<T> >& m, const std::vector<std::string>& coilStrings,
                    const std::string& xlabel, const std::string& ylabel, const std::string& title,
                    size_t xsize, size_t ysize, bool trueColor,
                    hoNDArray<float>& plotIm)
{
    try
    {
        size_t CHA = m.get_size(0);
        GADGET_CHECK_RETURN_FALSE(coilStrings.size() == CHA);

        hoNDArray<double> xd, yd, yd2;

        xd.create(CHA);
        yd.create(CHA);

        size_t c;
        for (c = 0; c < CHA; c++)
        {
            xd(c) = c+1;
            yd(c) = std::sqrt( std::abs(m(c, c)) );
        }

        double maxY = Gadgetron::max(&yd);

        yd2 = yd;
        std::sort(yd2.begin(), yd2.end());
        double medY = yd2(CHA / 2);

        // increase dot line to be 1 sigma ~= 33%
        double medRange = 0.33;

        if (maxY < medY*(1 + medRange))
        {
            maxY = medY*(1 + medRange);
        }

        hoNDArray<unsigned char> im;
        im.create(3, xsize, ysize);
        Gadgetron::clear(im);

        plsdev("mem");

        plsmem(im.get_size(1), im.get_size(2), im.begin());

        plinit();
        plfont(2);
        pladv(0);
        plvpor(0.15, 0.75, 0.1, 0.8);

        plwind(0, CHA+1, 0, maxY*1.05);

        plcol0(15);
        plbox("bcnst", 0.0, 0, "bcnstv", 0.0, 0);

        std::string gly;
        getPlotGlyph(0, gly); // circle
        plstring(CHA, xd.begin(), yd.begin(), gly.c_str());

        // draw the median line
        pllsty(1);

        double px[2], py[2];

        px[0] = 0;
        px[1] = CHA+1;

        py[0] = medY;
        py[1] = medY;

        plline(2, px, py);

        pllsty(2);

        py[0] = medY*(1 - medRange);
        py[1] = medY*(1 - medRange);

        plline(2, px, py);

        py[0] = medY*(1 + medRange);
        py[1] = medY*(1 + medRange);

        plline(2, px, py);

        plmtex("b", 3.2, 0.5, 0.5, xlabel.c_str());
        plmtex("t", 2.0, 0.5, 0.5, title.c_str());
        plmtex("l", 5.0, 0.5, 0.5, ylabel.c_str());

        // draw the legend
        std::vector<PLINT> opt_array(CHA), text_colors(CHA), line_colors(CHA), line_styles(CHA), symbol_numbers(CHA), symbol_colors(CHA);
        std::vector<PLFLT> symbol_scales(CHA), line_widths(CHA), box_scales(CHA, 1);

        std::vector<const char*> symbols(CHA);
        PLFLT legend_width, legend_height;

        std::vector<const char*> legend_text(CHA);

        std::vector<std::string> legends(CHA);

        size_t n;
        for (n = 0; n < CHA; n++)
        {
            opt_array[n] = PL_LEGEND_SYMBOL;
            text_colors[n] = 15;
            line_colors[n] = 15;
            line_styles[n] = (n % 8 + 1);
            line_widths[n] = 0.2;
            symbol_colors[n] = 15;
            symbol_scales[n] = 0.75;
            symbol_numbers[n] = 1;
            symbols[n] = gly.c_str();

            std::ostringstream ostr;
            ostr << n+1 << ":" << coilStrings[n];

            legends[n] = ostr.str();

            legend_text[n] = legends[n].c_str();
        }

        pllegend(&legend_width,
            &legend_height,
            PL_LEGEND_BACKGROUND,
            PL_POSITION_OUTSIDE | PL_POSITION_RIGHT,
            0.02,                                       // x
            0.0,                                        // y
            0.05,                                       // plot_width
            0,                                          // bg_color
            15,                                         // bb_color
            1,                                          // bb_style
            0,                                          // nrow
            0,                                          // ncolumn
            CHA,                                        // nlegend
            &opt_array[0],
            0.05,                                       // text_offset
            0.5,                                        // text_scale
            1.0,                                        // text_spacing
            0.5,                                        // text_justification
            &text_colors[0],
            (const char **)(&legend_text[0]),
            NULL,                                       // box_colors
            NULL,                                       // box_patterns
            &box_scales[0],                             // box_scales
            NULL,                                       // box_line_widths
            &line_colors[0],
            &line_styles[0],
            &line_widths[0],
            &symbol_colors[0],
            &symbol_scales[0],
            &symbol_numbers[0],
            (const char **)(&symbols[0])
            );

        plend();

        outputPlotIm(im, trueColor, plotIm);
    }
    catch (...)
    {
        GERROR_STREAM("Errors happened in plotNoiseStandardDeviation(...) ... ");
        return false;
    }

    return true;
}

template EXPORTGTPLPLOT bool plotNoiseStandardDeviation(const hoNDArray< std::complex<float> >& m, const std::vector<std::string>& coilStrings,
    const std::string& xlabel, const std::string& ylabel, const std::string& title, size_t xsize, size_t ysize, bool trueColor, hoNDArray<float>& plotIm);

template EXPORTGTPLPLOT bool plotNoiseStandardDeviation(const hoNDArray< std::complex<double> >& m, const std::vector<std::string>& coilStrings,
    const std::string& xlabel, const std::string& ylabel, const std::string& title, size_t xsize, size_t ysize, bool trueColor, hoNDArray<float>& plotIm);

}
