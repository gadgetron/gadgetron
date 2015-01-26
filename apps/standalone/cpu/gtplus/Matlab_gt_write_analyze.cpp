
#include <matrix.h>
#include <mat.h>
#ifdef _WIN32
    #include <mexGT.h>
#else
    #include <mex.h>
#endif // _WIN32

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
    printAuthorInfo(outs);
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
        using namespace Gadgetron::gtPlus;

        Gadgetron::gtMatlabConverter<float> converterFloat;
        Gadgetron::gtMatlabConverter<double> converterDouble;
        Gadgetron::gtMatlabConverter<short> converterShort;

        // ---------------------------------------------------------------
        // input parameters
        // ---------------------------------------------------------------    
        gtPlusIOAnalyze gt_io;

        const mxArray* aMx = prhs[0];
        const mxArray* aHeader = prhs[1];

        std::string filename;
        converterFloat.Matlab2Str(prhs[2], filename);

        try
        {
            hoNDImage<float, 2> data;
            if ( converterFloat.Matlab2hoNDImage(aMx, aHeader, data) )
            {
                gt_io.exportImage(data, filename);
            }
            else
            {
                hoNDImage<double, 2> data;
                if ( converterDouble.Matlab2hoNDImage(aMx, aHeader, data) )
                {
                    gt_io.exportImage(data, filename);
                }
                else
                {
                    hoNDImage<short, 2> data;
                    if ( converterShort.Matlab2hoNDImage(aMx, aHeader, data) )
                    {
                        gt_io.exportImage(data, filename);
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
                if ( converterFloat.Matlab2hoNDImage(aMx, aHeader, data) )
                {
                    gt_io.exportImage(data, filename);
                }
                else
                {
                    hoNDImage<double, 3> data;
                    if ( converterDouble.Matlab2hoNDImage(aMx, aHeader, data) )
                    {
                        gt_io.exportImage(data, filename);
                    }
                    else
                    {
                        hoNDImage<short, 3> data;
                        if ( converterShort.Matlab2hoNDImage(aMx, aHeader, data) )
                        {
                            gt_io.exportImage(data, filename);
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
                    if ( converterFloat.Matlab2hoNDImage(aMx, aHeader, data) )
                    {
                        gt_io.exportImage(data, filename);
                    }
                    else
                    {
                        hoNDImage<double, 4> data;
                        if ( converterDouble.Matlab2hoNDImage(aMx, aHeader, data) )
                        {
                            gt_io.exportImage(data, filename);
                        }
                        else
                        {
                            hoNDImage<short, 4> data;
                            if ( converterShort.Matlab2hoNDImage(aMx, aHeader, data) )
                            {
                                gt_io.exportImage(data, filename);
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
                        if ( converterFloat.Matlab2hoNDImage(aMx, aHeader, data) )
                        {
                            gt_io.exportImage(data, filename);
                        }
                        else
                        {
                            hoNDImage<double, 5> data;
                            if ( converterDouble.Matlab2hoNDImage(aMx, aHeader, data) )
                            {
                                gt_io.exportImage(data, filename);
                            }
                            else
                            {
                                hoNDImage<short, 5> data;
                                if ( converterShort.Matlab2hoNDImage(aMx, aHeader, data) )
                                {
                                    gt_io.exportImage(data, filename);
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
                            if ( converterFloat.Matlab2hoNDImage(aMx, aHeader, data) )
                            {
                                gt_io.exportImage(data, filename);
                            }
                            else
                            {
                                hoNDImage<double, 6> data;
                                if ( converterDouble.Matlab2hoNDImage(aMx, aHeader, data) )
                                {
                                    gt_io.exportImage(data, filename);
                                }
                                else
                                {
                                    hoNDImage<short, 6> data;
                                    if ( converterShort.Matlab2hoNDImage(aMx, aHeader, data) )
                                    {
                                        gt_io.exportImage(data, filename);
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
