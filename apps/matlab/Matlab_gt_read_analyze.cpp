
#include <matrix.h>
#include <mat.h>

#ifdef MATLAB_DLL_EXPORT_SYM
    #define DLL_EXPORT_SYM extern "C" __declspec(dllexport)
#endif // MATLAB_DLL_EXPORT_SYM
#include <mex.h>

// Gadgetron includes
#include "gtMatlab.h"
#include "gtPlusISMRMRDReconUtil.h"
#include "hoNDArray.h"
#include "hoNDArray_fileio.h"
#include "hoNDPoint.h"
#include "hoNDImage.h"

#include "gtMatlabConverter.h"
#include "gtMatlabConverterComplex.h"

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
    printAuthorInfo(outs);
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
        using namespace Gadgetron::gtPlus;

        Gadgetron::gtMatlabConverter<float> converterFloat;
        Gadgetron::gtMatlabConverter<double> converterDouble;
        Gadgetron::gtMatlabConverter<short> converterShort;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        // file name
        std::string filename;
        converterFloat.Matlab2Str(prhs[0], filename);

        gtPlusIOAnalyze gt_io;

        mxArray* aMx = NULL;
        mxArray* aHeader = NULL;

        try
        {
            hoNDImage<float, 2> data;
            if ( gt_io.importImage(data, filename) )
            {
                converterFloat.hoNDImage2Matlab(data, aMx, aHeader);
            }
            else
            {
                hoNDImage<double, 2> data;
                if ( gt_io.importImage(data, filename) )
                {
                    converterDouble.hoNDImage2Matlab(data, aMx, aHeader);
                }
                else
                {
                    hoNDImage<short, 2> data;
                    if ( gt_io.importImage(data, filename) )
                    {
                        converterShort.hoNDImage2Matlab(data, aMx, aHeader);
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
                if ( gt_io.importImage(data, filename) )
                {
                    converterFloat.hoNDImage2Matlab(data, aMx, aHeader);
                }
                else
                {
                    hoNDImage<double, 3> data;
                    if ( gt_io.importImage(data, filename) )
                    {
                        converterDouble.hoNDImage2Matlab(data, aMx, aHeader);
                    }
                    else
                    {
                        hoNDImage<short, 3> data;
                        if ( gt_io.importImage(data, filename) )
                        {
                            converterShort.hoNDImage2Matlab(data, aMx, aHeader);
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
                    if ( gt_io.importImage(data, filename) )
                    {
                        converterFloat.hoNDImage2Matlab(data, aMx, aHeader);
                    }
                    else
                    {
                        hoNDImage<double, 4> data;
                        if ( gt_io.importImage(data, filename) )
                        {
                            converterDouble.hoNDImage2Matlab(data, aMx, aHeader);
                        }
                        else
                        {
                            hoNDImage<short, 4> data;
                            if ( gt_io.importImage(data, filename) )
                            {
                                converterShort.hoNDImage2Matlab(data, aMx, aHeader);
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
                        if ( gt_io.importImage(data, filename) )
                        {
                            converterFloat.hoNDImage2Matlab(data, aMx, aHeader);
                        }
                        else
                        {
                            hoNDImage<double, 5> data;
                            if ( gt_io.importImage(data, filename) )
                            {
                                converterDouble.hoNDImage2Matlab(data, aMx, aHeader);
                            }
                            else
                            {
                                hoNDImage<short, 5> data;
                                if ( gt_io.importImage(data, filename) )
                                {
                                    converterShort.hoNDImage2Matlab(data, aMx, aHeader);
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
                            if ( gt_io.importImage(data, filename) )
                            {
                                converterFloat.hoNDImage2Matlab(data, aMx, aHeader);
                            }
                            else
                            {
                                hoNDImage<double, 6> data;
                                if ( gt_io.importImage(data, filename) )
                                {
                                    converterDouble.hoNDImage2Matlab(data, aMx, aHeader);
                                }
                                else
                                {
                                    hoNDImage<short, 6> data;
                                    if ( gt_io.importImage(data, filename) )
                                    {
                                        converterShort.hoNDImage2Matlab(data, aMx, aHeader);
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
