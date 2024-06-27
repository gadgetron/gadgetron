#include "RemoveROOversamplingGadget.h"

namespace Gadgetron {

RemoveROOversamplingGadget::RemoveROOversamplingGadget(const Core::Context& context,
                                                       const Core::GadgetProperties& props)
    : Core::ChannelGadget<Core::Acquisition>(context, props) {
    auto h = (context.header);

    if (h.encoding.size() == 0) {
        GDEBUG("Number of encoding spaces: %d\n", h.encoding.size());
        GERROR("This Gadget needs an encoding description\n");
        return;
    }

    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;

    encodeNx_ = e_space.matrixSize.x;
    encodeFOV_ = e_space.fieldOfView_mm.x;
    reconNx_ = r_space.matrixSize.x;
    reconFOV_ = r_space.fieldOfView_mm.x;

#ifdef USE_OMP  // TODO: Should this be removed? Its from the old version
    omp_set_num_threads(1);
    GDEBUG_STREAM("RemoveROOversamplingGadget:omp_set_num_threads(1) ... ");
#endif // USE_OMP

    // If the encoding and recon matrix size and FOV are the same
    // then the data is not oversampled and we can safely pass
    // the data onto the next gadget
    if ((encodeNx_ == reconNx_) && (encodeFOV_ == reconFOV_)) {
        dowork_ = false;
    } else {
        dowork_ = true;
    }
}

void RemoveROOversamplingGadget::process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out) {
  for (auto [header, acq, traj] : in) {
    if(dowork_){
      hoNDArray<std::complex<float>>* temp = new hoNDArray<std::complex<float>>();
      if (!temp)
      {
        GERROR("Error creating new temp  array");
        return;
      }
      std::vector<size_t> data_out_dims = *acq.get_dimensions();
      if (!ifft_buf_.dimensions_equal(data_out_dims) )
      {
          ifft_buf_.create(data_out_dims);
          ifft_res_.create(data_out_dims);
      }
      float ratioFOV = encodeFOV_/reconFOV_;
      data_out_dims[0] = (size_t)(data_out_dims[0]/ratioFOV);
      if ( !fft_buf_.dimensions_equal(data_out_dims) )
      {
          fft_buf_.create(data_out_dims);
          fft_res_.create(data_out_dims);
      }
      try{ temp->create(data_out_dims);}
      catch (std::runtime_error &err)
      {
          GEXCEPTION(err,"Unable to create new data array for downsampled data\n");
          return;
      }
      size_t sRO = acq.get_size(0);
      size_t start = (size_t)((acq.get_size(0)-data_out_dims[0]) / 2);

      size_t dRO = temp->get_size(0);
      size_t numOfBytes = data_out_dims[0]*sizeof(std::complex<float>);

      int c;   
      int CHA = (int)(data_out_dims[1]);

      std::complex<float>* data_in, *data_out;
      hoNDFFT<float>::instance()->ifft1c(acq, ifft_res_, ifft_buf_);
      data_in  = ifft_res_.get_data_ptr();
      data_out = temp->get_data_ptr();

      for (c=0; c<CHA; c++)
      {
          memcpy( data_out+c*dRO, data_in+c*sRO+start, numOfBytes );
      }

      hoNDFFT<float>::instance()->fft1c(*temp, fft_res_, fft_buf_);

      memcpy(temp->begin(), fft_res_.begin(), fft_res_.get_number_of_bytes());

      header.number_of_samples = data_out_dims[0];
      header.center_sample = (uint16_t)(header.center_sample/ratioFOV);
      header.discard_pre = (uint16_t)(header.discard_pre / ratioFOV);
      header.discard_post = (uint16_t)(header.discard_post / ratioFOV);

      acq = *temp;
    }
    out.push(Core::Acquisition{header, std::move(acq), std::move(traj)});
  }
}
GADGETRON_GADGET_EXPORT(RemoveROOversamplingGadget)
} // namespace Gadgetron
