
#include "cuFFTGadget.h"
#include "cuNDFFT.h"
#include "cuFFTCachedPlan.h"
namespace Gadgetron{

  void cuFFTGadget::process(Core::InputChannel<Core::Image<std::complex<float>>>& in, Core::OutputChannel& out) {

      auto plan = cuFFTCachedPlan<complext<float>>{};
      for (auto [header, data, meta] : in){
          auto * tmp = (hoNDArray<complext<float>>*) &data;
          cuNDArray< complext<float> > cu_data(*tmp);

          cu_data.squeeze();
          plan.ifft3c(cu_data);
          cu_data.to_host(tmp);


          out.push(header, std::move(data),std::move(meta));
      }

  }

GADGETRON_GADGET_EXPORT(cuFFTGadget);

  }

