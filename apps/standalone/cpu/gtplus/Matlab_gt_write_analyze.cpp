
#include <matrix.h>
#include <mat.h>
#ifdef MATLAB_DLL_EXPORT_SYM
    #define DLL_EXPORT_SYM extern "C" __declspec(dllexport)
#endif // MATLAB_DLL_EXPORT_SYM
#include <mex.h>

// Gadgetron includes
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDPoint.h"
#include "hoNDImage.h"
#include "ImageIOAnalyze.h"

#include "MatlabUtils.h"

#define MEXPRINTF(name) mexPrintf(#name);

#define NIn 3
#define NOut 0

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: Matlab_gt_write_analyze \n";
    outs << "Write out the Gadgetron produced Analyze image file" << endl;
    outs << "Support 2D/3D/4D/5D/6D images with short/float/double data types" << endl;
    outs << "3 Input paras:" << endl;
    outs << '\t' << "data       : image data" << endl;
    outs << '\t' << "header     : image header" << endl;
    outs << '\t' << "filename   : file name of the analyze image, no .hdr or .img extension needed" << endl;

    outs << "0 Output para" << endl;
    outs << "==============================================================================================" << endl;
    outs << std::ends; 

    std::string msg = outs.str();
    mexPrintf("%s\n", msg.c_str() );
}

template <typename T, unsigned int D> 
bool convert_from_matlab(const mxArray* aMx, const mxArray* aHeader, Gadgetron::hoNDImage<T, D>& data)
{
    try
    {
        Gadgetron::MatlabToHoNDImage(aMx, aHeader, data);
    }
    catch(...)
    {
        return false;
    }

    return true;
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    try
    {
        // ---------------------------------------------------------------
        // consistency check
        // ---------------------------------------------------------------    
        if (nrhs != NIn) 
        {
            mexWarnMsgTxt("3 input arguments are required ...");
            usage();
            return;
        }

        using namespace Gadgetron;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        ImageIOAnalyze gt_io;

        const mxArray* aMx = prhs[0];
        const mxArray* aHeader = prhs[1];

        std::string filename = Gadgetron::MatlabToStdString(prhs[2]);

        try
        {
            hoNDImage<float, 2> data;
            if ( convert_from_matlab(aMx, aHeader, data) )
            {
                gt_io.export_image(data, filename);
            }
            else
            {
                hoNDImage<double, 2> data;
                if ( convert_from_matlab(aMx, aHeader, data) )
                {
                    gt_io.export_image(data, filename);
                }
                else
                {
                    hoNDImage<short, 2> data;
                    if ( convert_from_matlab(aMx, aHeader, data) )
                    {
                        gt_io.export_image(data, filename);
                    }
                    else
                    {
                        throw("not 2D ... ");
                    }
                }
            }
        }
        catch(...)
        {
            try
            {
                hoNDImage<float, 3> data;
                if ( convert_from_matlab(aMx, aHeader, data) )
                {
                    gt_io.export_image(data, filename);
                }
                else
                {
                    hoNDImage<double, 3> data;
                    if ( convert_from_matlab(aMx, aHeader, data) )
                    {
                        gt_io.export_image(data, filename);
                    }
                    else
                    {
                        hoNDImage<short, 3> data;
                        if ( convert_from_matlab(aMx, aHeader, data) )
                        {
                            gt_io.export_image(data, filename);
                        }
                        else
                        {
                            throw("not 3D ... ");
                        }
                    }
                }
            }
            catch(...)
            {
                try
                {
                    hoNDImage<float, 4> data;
                    if ( convert_from_matlab(aMx, aHeader, data) )
                    {
                        gt_io.export_image(data, filename);
                    }
                    else
                    {
                        hoNDImage<double, 4> data;
                        if ( convert_from_matlab(aMx, aHeader, data) )
                        {
                            gt_io.export_image(data, filename);
                        }
                        else
                        {
                            hoNDImage<short, 4> data;
                            if ( convert_from_matlab(aMx, aHeader, data) )
                            {
                                gt_io.export_image(data, filename);
                            }
                            else
                            {
                                throw("not 4D ... ");
                            }
                        }
                    }
                }
                catch(...)
                {
                    try
                    {
                        hoNDImage<float, 5> data;
                        if ( convert_from_matlab(aMx, aHeader, data) )
                        {
                            gt_io.export_image(data, filename);
                        }
                        else
                        {
                            hoNDImage<double, 5> data;
                            if ( convert_from_matlab(aMx, aHeader, data) )
                            {
                                gt_io.export_image(data, filename);
                            }
                            else
                            {
                                hoNDImage<short, 5> data;
                                if ( convert_from_matlab(aMx, aHeader, data) )
                                {
                                    gt_io.export_image(data, filename);
                                }
                                else
                                {
                                    throw("not 5D ... ");
                                }
                            }
                        }
                    }
                    catch(...)
                    {
                        try
                        {
                            hoNDImage<float, 6> data;
                            if ( convert_from_matlab(aMx, aHeader, data) )
                            {
                                gt_io.export_image(data, filename);
                            }
                            else
                            {
                                hoNDImage<double, 6> data;
                                if ( convert_from_matlab(aMx, aHeader, data) )
                                {
                                    gt_io.export_image(data, filename);
                                }
                                else
                                {
                                    hoNDImage<short, 6> data;
                                    if ( convert_from_matlab(aMx, aHeader, data) )
                                    {
                                        gt_io.export_image(data, filename);
                                    }
                                    else
                                    {
                                        throw("not 6D ... ");
                                    }
                                }
                            }
                        }
                        catch(...)
                        {
                            mexWarnMsgTxt("Images must be 2D/3D/4D/5D/6D ...");
                            return;
                        }
                    }
                }
            }
        }
    }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab_gt_write_analyze(...) ...");
        return;
    }

    return;
}
