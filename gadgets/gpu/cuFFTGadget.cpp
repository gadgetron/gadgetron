
#include "cuFFTGadget.h"
#include "cuNDFFT.h"
#include "cuFFTPlan.h"
namespace Gadgetron{

  void cuFFTGadget::process(Core::InputChannel<mrd::Image<std::complex<float>>>& in, Core::OutputChannel& out) {

      for (auto [header, data, meta] : in){
          auto * tmp = (hoNDArray<complext<float>>*) &data;
          cuNDArray< complext<float> > cu_data(*tmp);

          cu_data.squeeze();
          cuNDFFT<float>::instance()->ifft(&cu_data);
          cu_data.to_host(tmp);


          out.push(header, std::move(data),std::move(meta));
      }

  }

GADGETRON_GADGET_EXPORT(cuFFTGadget);

  }

