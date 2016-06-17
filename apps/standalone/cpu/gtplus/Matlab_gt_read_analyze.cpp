
#include <matrix.h>
#include <mat.h>
#ifdef MATLAB_DLL_EXPORT_SYM
    #define DLL_EXPORT_SYM extern "C" __declspec(dllexport)
#endif // MATLAB_DLL_EXPORT_SYM
#include <mex.h>

// Gadgetron includes
#include "ImageIOAnalyze.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDPoint.h"
#include "hoNDImage.h"

#include "MatlabUtils.h"

#define MEXPRINTF(name) mexPrintf(#name);

#define NIn 1
#define NOut 2

static void usage()
{
    using namespace std;
    std::stringstream outs;

    outs << "==============================================================================================" << endl;
    outs << "Usage: Matlab_gt_read_analyze \n";
    outs << "Read in the Gadgetron produced Analyze image file" << endl;
    outs << "Support 2D/3D/4D/5D/6D images with short/float/double data types" << endl;
    outs << "1 Input paras:" << endl;
    outs << '\t' << "filename   : file name of the analyze image, no .hdr or .img extension needed" << endl;

    outs << "2 Output para:" << endl;
    outs << '\t' << "data       : image data" << endl;
    outs << '\t' << "header     : image header" << endl;
    outs << "==============================================================================================" << endl;
    outs << std::ends; 

    std::string msg = outs.str();
    mexPrintf("%s\n", msg.c_str() );
}

template <typename T, unsigned int D>
bool import_image_call(Gadgetron::ImageIOAnalyze& gt_io, Gadgetron::hoNDImage<T, D>& data, std::string& filename)
{
    try
    {
        gt_io.import_image(data, filename);
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
            mexWarnMsgTxt("1 input arguments are required ...");
            usage();
            return;
        }

        if (nlhs < NOut )
        {
            mexWarnMsgTxt("2 output argument is required ...");
            usage();
            return;
        }

        using namespace Gadgetron;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // file name
        std::string filename = Gadgetron::MatlabToStdString(prhs[0]);

        mxArray* aMx = NULL;
        mxArray* aHeader = NULL;

        ImageIOAnalyze gt_io;

        try
        {
            hoNDImage<float, 2> data;
            if ( import_image_call(gt_io, data, filename) )
            {
                aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
            }
            else
            {
                hoNDImage<double, 2> data;
                if ( import_image_call(gt_io, data, filename) )
                {
                    aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                }
                else
                {
                    hoNDImage<short, 2> data;
                    if ( import_image_call(gt_io, data, filename) )
                    {
                        aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
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
                if ( import_image_call(gt_io, data, filename) )
                {
                    aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                }
                else
                {
                    hoNDImage<double, 3> data;
                    if ( import_image_call(gt_io, data, filename) )
                    {
                        aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                    }
                    else
                    {
                        hoNDImage<short, 3> data;
                        if ( import_image_call(gt_io, data, filename) )
                        {
                            aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
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
                    if ( import_image_call(gt_io, data, filename) )
                    {
                        aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                    }
                    else
                    {
                        hoNDImage<double, 4> data;
                        if ( import_image_call(gt_io, data, filename) )
                        {
                            aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                        }
                        else
                        {
                            hoNDImage<short, 4> data;
                            if ( import_image_call(gt_io, data, filename) )
                            {
                                aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
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
                        if ( import_image_call(gt_io, data, filename) )
                        {
                            aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                        }
                        else
                        {
                            hoNDImage<double, 5> data;
                            if ( import_image_call(gt_io, data, filename) )
                            {
                                aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                            }
                            else
                            {
                                hoNDImage<short, 5> data;
                                if ( import_image_call(gt_io, data, filename) )
                                {
                                    aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
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
                            if ( import_image_call(gt_io, data, filename) )
                            {
                                aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                            }
                            else
                            {
                                hoNDImage<double, 6> data;
                                if ( import_image_call(gt_io, data, filename) )
                                {
                                    aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
                                }
                                else
                                {
                                    hoNDImage<short, 6> data;
                                    if ( import_image_call(gt_io, data, filename) )
                                    {
                                        aMx = Gadgetron::hoNDImageToMatlab(&data, aHeader);
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

        // ---------------------------------------------------------------
        // output parameter
        // ---------------------------------------------------------------
        plhs[0] = aMx;
        plhs[1] = aHeader;
    }
    catch(...)
    {
        mexWarnMsgTxt("Exceptions happened in Matlab_gt_read_analyze(...) ...");
        return;
    }

    return;
}
