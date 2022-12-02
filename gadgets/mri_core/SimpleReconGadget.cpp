#include "SimpleReconGadget.h"

namespace Gadgetron {

    SimpleReconGadget::SimpleReconGadget(const Core::Context& context, const Core::GadgetProperties& props)
        : ChannelGadget(context, props) {
        header = context.header;
        image_counter_ = 0;
    }

    void SimpleReconGadget::process(Core::InputChannel<IsmrmrdReconData>& input, Core::OutputChannel& out) {

        for (IsmrmrdReconData reconData : input) {
            // Iterate over all the recon bits
            for (std::vector<IsmrmrdReconBit>::iterator it = reconData.rbit_.begin(); it != reconData.rbit_.end(); ++it) {
                // Grab a reference to the buffer containing the imaging data
                // We are ignoring the reference data
                IsmrmrdDataBuffered& dbuff = it->data_;

                // Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
                uint16_t E0 = dbuff.data_.get_size(0);
                uint16_t E1 = dbuff.data_.get_size(1);
                uint16_t E2 = dbuff.data_.get_size(2);
                uint16_t CHA = dbuff.data_.get_size(3);
                uint16_t N = dbuff.data_.get_size(4);
                uint16_t S = dbuff.data_.get_size(5);
                uint16_t LOC = dbuff.data_.get_size(6);

                //Create an image array message
                IsmrmrdImageArray* cm1 = new IsmrmrdImageArray();

                //Grab references to the image array data and headers
                IsmrmrdImageArray & imarray = *cm1;

                // The image array data will be [E0,E1,E2,1,N,S,LOC] big
                // Will collapse across coils at the end
                std::vector<size_t> data_dims(7);
                data_dims[0] = E0;
                data_dims[1] = E1;
                data_dims[2] = E2;
                data_dims[3] = 1;
                data_dims[4] = N;
                data_dims[5] = S;
                data_dims[6] = LOC;
                imarray.data_.create(data_dims);

                // ImageHeaders will be [N, S, LOC]
                std::vector<size_t> header_dims(3);
                header_dims[0] = N;
                header_dims[1] = S;
                header_dims[2] = LOC;
                imarray.headers_.create(header_dims);

                // Loop over S and N and LOC
                for (uint16_t loc = 0; loc < LOC; loc++) {
                    for (uint16_t s = 0; s < S; s++) {
                        for (uint16_t n = 0; n < N; n++) {
                            // Set some information into the image header
                            // Use the middle acquisition header for some info
                            //[E1, E2, N, S, LOC]
                            ISMRMRD::AcquisitionHeader& acqhdr =
                                dbuff.headers_(dbuff.sampling_.sampling_limits_[1].center_,
                                            dbuff.sampling_.sampling_limits_[2].center_, n, s, loc);
                            imarray.headers_(n, s, loc).measurement_uid = acqhdr.measurement_uid;
                            imarray.headers_(n, s, loc).matrix_size[0] = E0;
                            imarray.headers_(n, s, loc).matrix_size[1] = E1;
                            imarray.headers_(n, s, loc).matrix_size[2] = E2;
                            imarray.headers_(n, s, loc).field_of_view[0] = dbuff.sampling_.recon_FOV_[0];
                            imarray.headers_(n, s, loc).field_of_view[1] = dbuff.sampling_.recon_FOV_[1];
                            imarray.headers_(n, s, loc).field_of_view[2] = dbuff.sampling_.recon_FOV_[2];
                            imarray.headers_(n, s, loc).channels = 1;
                            imarray.headers_(n, s, loc).average = acqhdr.idx.average;
                            imarray.headers_(n, s, loc).slice = acqhdr.idx.slice;
                            imarray.headers_(n, s, loc).contrast = acqhdr.idx.contrast;
                            imarray.headers_(n, s, loc).phase = acqhdr.idx.phase;
                            imarray.headers_(n, s, loc).repetition = acqhdr.idx.repetition;
                            imarray.headers_(n, s, loc).set = acqhdr.idx.set;
                            imarray.headers_(n, s, loc).acquisition_time_stamp = acqhdr.acquisition_time_stamp;
                            imarray.headers_(n, s, loc).position[0] = acqhdr.position[0];
                            imarray.headers_(n, s, loc).position[1] = acqhdr.position[1];
                            imarray.headers_(n, s, loc).position[2] = acqhdr.position[2];
                            imarray.headers_(n, s, loc).read_dir[0] = acqhdr.read_dir[0];
                            imarray.headers_(n, s, loc).read_dir[1] = acqhdr.read_dir[1];
                            imarray.headers_(n, s, loc).read_dir[2] = acqhdr.read_dir[2];
                            imarray.headers_(n, s, loc).phase_dir[0] = acqhdr.phase_dir[0];
                            imarray.headers_(n, s, loc).phase_dir[1] = acqhdr.phase_dir[1];
                            imarray.headers_(n, s, loc).phase_dir[2] = acqhdr.phase_dir[2];
                            imarray.headers_(n, s, loc).slice_dir[0] = acqhdr.slice_dir[0];
                            imarray.headers_(n, s, loc).slice_dir[1] = acqhdr.slice_dir[1];
                            imarray.headers_(n, s, loc).slice_dir[2] = acqhdr.slice_dir[2];
                            imarray.headers_(n, s, loc).patient_table_position[0] = acqhdr.patient_table_position[0];
                            imarray.headers_(n, s, loc).patient_table_position[1] = acqhdr.patient_table_position[1];
                            imarray.headers_(n, s, loc).patient_table_position[2] = acqhdr.patient_table_position[2];
                            imarray.headers_(n, s, loc).data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                            imarray.headers_(n, s, loc).image_index = ++image_counter_;

                            // Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
                            // Each chunk will be [E0,E1,E2,CHA] big
                            std::vector<size_t> chunk_dims(4);
                            chunk_dims[0] = E0;
                            chunk_dims[1] = E1;
                            chunk_dims[2] = E2;
                            chunk_dims[3] = CHA;
                            hoNDArray<std::complex<float>> chunk =
                                hoNDArray<std::complex<float>>(chunk_dims, &dbuff.data_(0, 0, 0, 0, n, s, loc));

                            // Do the FFTs in place
                            hoNDFFT<float>::instance()->ifft3c(chunk);

                            // Square root of the sum of squares
                            // Each image will be [E0,E1,E2,1] big
                            std::vector<size_t> img_dims(3);
                            img_dims[0] = E0;
                            img_dims[1] = E1;
                            img_dims[2] = E2;
                            hoNDArray<std::complex<float>> output =
                                hoNDArray<std::complex<float>>(img_dims, &imarray.data_(0, 0, 0, 0, n, s, loc));

                            // Compute d* d in place
                            multiplyConj(chunk, chunk, chunk);
                            // Add up
                            for (size_t c = 0; c < CHA; c++) {
                                output += hoNDArray<std::complex<float>>(img_dims, &chunk(0, 0, 0, c));
                            }
                            // Take the square root in place
                            sqrt_inplace(&output);
                        }
                    }
                }
                out.push(std::move(imarray));
            }
        }
    }
    GADGETRON_GADGET_EXPORT(SimpleReconGadget);

}
