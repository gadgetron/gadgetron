#include "FFTGadget.h"

namespace Gadgetron{

    FFTGadget::FFTGadget(const Core::Context& context, const Core::GadgetProperties& props) : ChannelGadget(context, props){
        image_counter_ = 0;
    }

    mrd::Image<std::complex<float>> CreateAndFFTImage(mrd::ReconBuffer dbuff, uint16_t n, uint16_t s, uint16_t loc, long long image_index){

        //7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
        uint16_t E0 = dbuff.data.get_size(0);
        uint16_t E1 = dbuff.data.get_size(1);
        uint16_t E2 = dbuff.data.get_size(2);
        uint16_t CHA = dbuff.data.get_size(3);
        uint16_t N = dbuff.data.get_size(4);
        uint16_t S = dbuff.data.get_size(5);
        uint16_t LOC = dbuff.data.get_size(6);

        //Each image will have size [E0,E1,E2,CHA]
        std::vector<size_t> img_dims(4);
        img_dims[0] = E0;
        img_dims[1] = E1;
        img_dims[2] = E2;
        img_dims[3] = CHA;

        //Set some information into the image header, and use the middle header for some info [E1, E2, N, S, LOC]
        mrd::AcquisitionHeader & acqhdr = dbuff.headers(dbuff.sampling.sampling_limits.kspace_encoding_step_1.center,
                                                             dbuff.sampling.sampling_limits.kspace_encoding_step_2.center,
                                                             n, s, loc);

        mrd::Image<std::complex<float>> image{};
        auto& imghdr = image.head;
        image.data.create(img_dims);

        /** TODO: This Gadget does NOT copy all Acquisition header fields to the Image header!
         *  e.g. measurement_uid and physiology_time_stamp are missing.
         *
         *  HOWEVER, the e2e test (epi_2d.yml) using this Gadget expects these fields to be empty.
         */
        imghdr.image_type = mrd::ImageType::kComplex;

        imghdr.field_of_view[0]   = dbuff.sampling.recon_fov.x;
        imghdr.field_of_view[1]   = dbuff.sampling.recon_fov.y;
        imghdr.field_of_view[2]   = dbuff.sampling.recon_fov.z;

        imghdr.average = acqhdr.idx.average;
        imghdr.slice = acqhdr.idx.slice;
        imghdr.contrast = acqhdr.idx.contrast;
        imghdr.phase = acqhdr.idx.phase;
        imghdr.repetition = acqhdr.idx.repetition;
        imghdr.set = acqhdr.idx.set;
        imghdr.acquisition_time_stamp = acqhdr.acquisition_time_stamp;

        imghdr.position = acqhdr.position;
        imghdr.col_dir = acqhdr.read_dir;
        imghdr.line_dir = acqhdr.phase_dir;
        imghdr.slice_dir = acqhdr.slice_dir;
        imghdr.patient_table_position = acqhdr.patient_table_position;
        imghdr.image_index = image_index;

        //Copy the 4D data block [E0,E1,E2,CHA] for this loc, n, and s into the output image
        memcpy(image.data.get_data_ptr(), &dbuff.data(0,0,0,0,n,s,loc), E0*E1*E2*CHA*sizeof(std::complex<float>));
        hoNDFFT<float>::instance()->ifft3c(image.data);

        return std::move(image);
    }

    void FFTGadget::process(Core::InputChannel<mrd::ReconData>& input, Core::OutputChannel& out) {
        for (auto reconData : input){
            //Iterate over all the recon bits
            for(auto& rbit : reconData.buffers)
            {
                //Grab a reference to the buffer containing the imaging data
                mrd::ReconBuffer & dbuff = rbit.data;

                uint16_t N = dbuff.data.get_size(4);
                uint16_t S = dbuff.data.get_size(5);
                uint16_t LOC = dbuff.data.get_size(6);

                //Loop over S and N and LOC
                for (uint16_t loc=0; loc < LOC; loc++) {
                    for (uint16_t s=0; s < S; s++) {
                        for (uint16_t n=0; n < N; n++) {
                            ++image_counter_;
                            out.push(CreateAndFFTImage(dbuff,n,s,loc,image_counter_));
                        }
                    }
                }
            }
        }
    }
    GADGETRON_GADGET_EXPORT(FFTGadget);
}
