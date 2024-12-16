#ifndef gpuSpiralSensePrepGadget_H
#define gpuSpiralSensePrepGadget_H
#pragma once

#include "Gadget.h"
#include "cuCgSolver.h"
#include "cuNonCartesianSenseOperator.h"
#include "cuCgPreconditioner.h"
#include "cuNFFT.h"
#include "hoNDArray.h"
#include "vector_td.h"
#include "cuNFFT.h"

#include "../../toolboxes/mri/spiral/TrajectoryParameters.h"
#include <boost/optional.hpp>
#include <boost/shared_array.hpp>
#include <boost/shared_ptr.hpp>
#include <complex>

namespace Gadgetron {

    class gpuSpiralSensePrepGadget :
            public Gadget1<mrd::Acquisition> {

    public:
        gpuSpiralSensePrepGadget();

        virtual ~gpuSpiralSensePrepGadget();

    protected:
        GADGET_PROPERTY(deviceno, int, "GPU device number", 0);
        GADGET_PROPERTY(propagate_csm_from_set, int, "Which set to use for CSM", -1);
        GADGET_PROPERTY(buffer_using_solver, bool, "Use solver for buffer", false);
        GADGET_PROPERTY(use_multiframe_grouping, bool, "Use multiframe grouping", false);
        GADGET_PROPERTY(buffer_convolution_kernel_width, float, "Convolution kernel width for buffer", 5.5);
        GADGET_PROPERTY(buffer_convolution_oversampling_factor, float, "Oversampling used in buffer convolution", 1.25);
        GADGET_PROPERTY(reconstruction_os_factor_x, float, "Oversampling for reconstruction in x-direction", 1.0);
        GADGET_PROPERTY(reconstruction_os_factor_y, float, "Oversampling for reconstruction in y-direction", 1.0);

        virtual int process_config(const mrd::Header& header);

        virtual int process(GadgetContainerMessage<mrd::Acquisition> *m1);

    private:
        int samples_to_skip_start_;
        int samples_to_skip_end_;
        int samples_per_interleave_;
        int interleaves_;
        int slices_;
        int sets_;
        std::vector<long> image_counter_;
        Spiral::TrajectoryParameters trajectoryParameters;
        int device_number_;


        std::vector<long> interleaves_counter_singleframe_;
        std::vector<long> interleaves_counter_multiframe_;
        long acceleration_factor_;

        bool prepared_;
        bool use_multiframe_grouping_;
        bool buffer_using_solver_;

        int propagate_csm_from_set_;

        float kernel_width_;
        float oversampling_factor_;


        hoNDArray<floatd2> host_traj_;
        hoNDArray<float> host_weights_;

        std::vector<hoNDArray<float_complext>> host_data_buffer_;
        boost::shared_ptr<cuNDArray<float> > dcw_buffer_;

        std::vector<size_t> fov_vec_;
        std::vector<size_t> image_dimensions_recon_;
        uint64d2 image_dimensions_recon_os_;

        boost::shared_ptr<cuNFFT_plan<float, 2>> nfft_plan_;
        cuCgSolver<float_complext> cg_;
        boost::shared_ptr<cuNDArray<float_complext> > csm_;
        boost::shared_ptr<cuNonCartesianSenseOperator<float, 2> > E_;
        boost::shared_ptr<cuCgPreconditioner<float_complext> > D_;

        std::vector<std::vector<GadgetContainerMessage<mrd::Acquisition>*>> buffer_;
        std::vector<std::vector<mrd::ImageHeader>> image_headers_queue_;

        void prepare_nfft( const mrd::Acquisition& acq);

        mrd::ImageHeader make_image_header(const mrd::AcquisitionHeader &acq_header);

        void change_acceleration_factor(const mrd::Acquisition &acq);

        cuNDArray<float_complext> make_reg_image(const hoNDArray<float_complext> &buffer, size_t set, size_t num_coils);

        std::tuple<boost::shared_ptr<hoNDArray<float_complext>>,
                boost::shared_ptr<hoNDArray<floatd2>>,
                boost::shared_ptr<hoNDArray<float>>> get_data_from_queues(size_t set, size_t slice, size_t num_coils);

        void setup_buffers(const mrd::Acquisition &acq);

    };
}
#endif //gpuSpiralSensePrepGadget_H
