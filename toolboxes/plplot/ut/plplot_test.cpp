
#ifdef USE_OMP
#include "omp.h"
#endif // USE_OMP

#include "ismrmrd/ismrmrd.h"
#include "complext.h"

#include <gtest/gtest.h>

#include "hoNDArray_utils.h"

#include "ImageIOAnalyze.h"

#include "GadgetronTimer.h"

#include "GtPLplot.h"

#ifdef max
#undef max
#endif // max

using namespace Gadgetron;
using testing::Types;

template <typename T> class gt_plplot_Test : public ::testing::Test 
{
protected:
    virtual void SetUp()
    {
        GDEBUG_STREAM("=============================================================================================");
        gt_ut_folder_ = std::string(::getenv("GADGETRON_UNITTEST_DIRECTORY"));
        GDEBUG_STREAM("=============================================================================================");
        GDEBUG_STREAM("Unit Test for GtPlus");
        gt_ut_data_folder_ = gt_ut_folder_ + "/data/";
        gt_ut_res_folder_ = gt_ut_folder_ + "/result/";
        GDEBUG_STREAM("gt_ut_data_folder_ is " << gt_ut_data_folder_);
        GDEBUG_STREAM("gt_ut_res_folder_ is " << gt_ut_res_folder_);

        timer_.set_timing_in_destruction(false);

#ifdef WIN32
    #ifdef USE_OMP
        /// lock the threads
        #pragma omp parallel default(shared)
        {
            int tid = omp_get_thread_num();
            // std::cout << tid << std::endl;
            DWORD_PTR mask = (1 << tid);
            SetThreadAffinityMask( GetCurrentThread(), mask );
        }
    #endif // USE_OMP
#endif // WIN32
    }

    std::string gt_ut_folder_;
    std::string gt_ut_data_folder_;
    std::string gt_ut_res_folder_;

    ImageIOAnalyze gt_io_;
    GadgetronTimer timer_;
};

typedef Types<float, double> realImplementations;

typedef Types< std::complex<float> > cpfloatImplementations;

typedef Types<std::complex<float>, std::complex<double>, float_complext, double_complext> cplxImplementations;
typedef Types<std::complex<float>, std::complex<double> > stdCplxImplementations;
typedef Types<float_complext, double_complext> cplxtImplementations;

TYPED_TEST_SUITE(gt_plplot_Test, cpfloatImplementations);

TYPED_TEST(gt_plplot_Test, plplot_noise_covariance_test)
{
    typedef std::complex<float> T;

    ImageIOAnalyze gt_io;

    float v;

    std::string xlabel = "Channel";
    std::string ylabel = "Noise STD";
    std::string title = "Gadgetron SNR toolkit, Noise STD Plot";
    size_t xsize = 2048;
    size_t ysize = 2048;
    hoNDArray<float> plotIm;

    bool trueColor = true;

    hoNDArray< std::complex<float> > m;
    m.create(32, 32);
    Gadgetron::fill(m, std::complex<float>(1.0, 0) );

    std::vector<std::string> coilStrings(32);

    for (size_t n=0; n < m.get_size(0); n++)
    {
        std::ostringstream ostr;
        ostr << "Channel" << "_" << n;

        coilStrings[n] = ostr.str();
    }

    Gadgetron::plotNoiseStandardDeviation(m, coilStrings, xlabel, ylabel, title, xsize, ysize, trueColor, plotIm);

    gt_io.export_array(plotIm, this->gt_ut_res_folder_ + "plplot_trueColor_NoiseSTD");
}

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

TYPED_TEST(gt_plplot_Test, plplot_curves_test)
{
    typedef float T;

    ImageIOAnalyze gt_io;

    float v;

    std::string xlabel = "Heart Beat";
    std::string ylabel = "RR interval (ms)";
    std::string title = "Acqusition heart beat plot";
    size_t xsize = 2048;
    size_t ysize = 2048;
    hoNDArray<float> plotIm;

    bool trueColor = true;

    std::vector<hoNDArray<T> > xs(1);
    std::vector<hoNDArray<T> > ys(1);
    std::vector<std::string> legend(1);

    size_t num_hb = 60;

    xs[0].create(num_hb);
    ys[0].create(num_hb);

    for (size_t n = 0; n < num_hb; n++)
    {
        xs[0](n) = n + 1;
        ys[0](n) = 1000 + (T)std::rand()/RAND_MAX * 100.0;
    }

    legend.clear();

    std::ostringstream ostr;
    ostr << "Acquisition median RR interval " << 1000 << " ms ";

    T xlim[2], ylim[2];
    xlim[0] = 1;
    xlim[1] = num_hb;

    ylim[0] = 0;
    ylim[1] = 2000;

    std::vector<std::string> symbols(1);
    symbols[0] = "#(225)";

    Gadgetron::plotCurves(xs, ys, "Heart Beat", "RR interval (ms)", ostr.str(), legend, symbols, xsize, ysize, xlim, ylim, trueColor, false, plotIm);

    gt_io.export_array(plotIm, this->gt_ut_res_folder_ + "plplot_trueColor_HeartBeat");
}
