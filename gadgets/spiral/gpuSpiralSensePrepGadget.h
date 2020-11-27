#ifndef gpuSpiralSensePrepGadget_H
#define gpuSpiralSensePrepGadget_H
#pragma once

#include "gadgetron_spiral_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
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
#include <ismrmrd/ismrmrd.h>
#include <ismrmrd/xml.h>

namespace Gadgetron {

    class EXPORTGADGETS_SPIRAL gpuSpiralSensePrepGadget :
            public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float> > > {

    public:
        GADGET_DECLARE(gpuSpiralSensePrepGadget);

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

        virtual int process_config(ACE_Message_Block *mb);

        virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
                            GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);

        virtual GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *
        duplicate_profile(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *profile);

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

        std::vector<std::vector<std::pair<GadgetContainerMessage<ISMRMRD::AcquisitionHeader>*,GadgetContainerMessage<hoNDArray<std::complex<float>>>*>>> buffer_;
        std::vector<std::vector<ISMRMRD::ImageHeader>> image_headers_queue_;

        void prepare_nfft( const ISMRMRD::AcquisitionHeader& header);

        ISMRMRD::ImageHeader make_image_header(const ISMRMRD::AcquisitionHeader &acq_header);

        void change_acceleration_factor(const ISMRMRD::AcquisitionHeader &header);

        cuNDArray<float_complext> make_reg_image(const hoNDArray<float_complext> &buffer, size_t set, size_t num_coils);

        std::tuple<boost::shared_ptr<hoNDArray<float_complext>>,
                boost::shared_ptr<hoNDArray<floatd2>>,
                boost::shared_ptr<hoNDArray<float>>> get_data_from_queues(size_t set, size_t slice, size_t num_coils);

        void setup_buffers(const ISMRMRD::AcquisitionHeader &header);

    };
}
#endif //gpuSpiralSensePrepGadget_H
