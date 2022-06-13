#include "FFTGadget.h"

namespace Gadgetron{

    FFTGadget::FFTGadget(const Core::Context& context, const Core::GadgetProperties& props) : ChannelGadget(context, props){
        header = context.header;
        image_counter_ = 0;
    }

    Core::Image<std::complex<float>> CreateAndFFTImage(IsmrmrdDataBuffered dbuff, uint16_t n, uint16_t s, uint16_t loc, long long image_index){
        
        //7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
        uint16_t E0 = dbuff.data_.get_size(0);
        uint16_t E1 = dbuff.data_.get_size(1);
        uint16_t E2 = dbuff.data_.get_size(2);
        uint16_t CHA = dbuff.data_.get_size(3);
        uint16_t N = dbuff.data_.get_size(4);
        uint16_t S = dbuff.data_.get_size(5);
        uint16_t LOC = dbuff.data_.get_size(6);

        //Each image will have size [E0,E1,E2,CHA]
        std::vector<size_t> img_dims(4);
        img_dims[0] = E0;
        img_dims[1] = E1;
        img_dims[2] = E2;
        img_dims[3] = CHA;

        //Set some information into the image header, and use the middle header for some info [E1, E2, N, S, LOC]
        ISMRMRD::AcquisitionHeader & acqhdr = dbuff.headers_(dbuff.sampling_.sampling_limits_[1].center_, 
                                                             dbuff.sampling_.sampling_limits_[2].center_,
                                                             n, s, loc);

        auto imageHeader = ISMRMRD::ImageHeader();
        auto imageData = hoNDArray<std::complex<float>>(img_dims);

        imageHeader.matrix_size[0]     = img_dims[0];
        imageHeader.matrix_size[1]     = img_dims[1];
        imageHeader.matrix_size[2]     = img_dims[2];
        imageHeader.field_of_view[0]   = dbuff.sampling_.recon_FOV_[0];
        imageHeader.field_of_view[1]   = dbuff.sampling_.recon_FOV_[1];
        imageHeader.field_of_view[2]   = dbuff.sampling_.recon_FOV_[2];
        imageHeader.channels           = img_dims[3];
        
        imageHeader.average = acqhdr.idx.average;
        imageHeader.slice = acqhdr.idx.slice;
        imageHeader.contrast = acqhdr.idx.contrast;
        imageHeader.phase = acqhdr.idx.phase;
        imageHeader.repetition = acqhdr.idx.repetition;
        imageHeader.set = acqhdr.idx.set;
        imageHeader.acquisition_time_stamp = acqhdr.acquisition_time_stamp;

        memcpy(imageHeader.position, acqhdr.position, sizeof(float)*3);
        memcpy(imageHeader.read_dir, acqhdr.read_dir, sizeof(float)*3);
        memcpy(imageHeader.phase_dir, acqhdr.phase_dir, sizeof(float)*3);
        memcpy(imageHeader.slice_dir, acqhdr.slice_dir, sizeof(float)*3);
        memcpy(imageHeader.patient_table_position, acqhdr.patient_table_position, sizeof(float)*3);
        imageHeader.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
        imageHeader.image_index = image_index;

        //Copy the 4D data block [E0,E1,E2,CHA] for this loc, n, and s into the output image
        memcpy(imageData.get_data_ptr(), &dbuff.data_(0,0,0,0,n,s,loc), E0*E1*E2*CHA*sizeof(std::complex<float>));
        hoNDFFT<float>::instance()->ifft3c(imageData);

        return Core::Image<std::complex<float>>(std::move(imageHeader), std::move(imageData), std::nullopt));

    }

    Core::Image<std::complex<float>> PerformFFT(Core::Image<std::complex<float>> image){
        auto data = std::get<hoNDArray<std::complex<float>>>(image);
        hoNDFFT<float>::instance()->ifft3c(data);
        return image;
    }

    void FFTGadget::process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) {
        for (IsmrmrdReconData reconData : input){
            //Iterate over all the recon bits
            for(std::vector<IsmrmrdReconBit>::iterator it = reconData.rbit_.begin(); it != reconData.rbit_.end(); ++it)
            {
                //Grab a reference to the buffer containing the imaging data
                IsmrmrdDataBuffered & dbuff = it->data_;

                uint16_t N = dbuff.data_.get_size(4);
                uint16_t S = dbuff.data_.get_size(5);
                uint16_t LOC = dbuff.data_.get_size(6);

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
