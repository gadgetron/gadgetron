#include "CropAndCombineGadget.h"
#include "hoNDArray_math.h"

using namespace Gadgetron;
using namespace Gadgetron::Core;

namespace {
  template <class T> Image<T> cropAndCombine(Image<T>& image) {
      GDEBUG("CropAndCombineGadget is not well defined for real-valued images. Doing nothing.");
      return image;
  }

  template <class T> Image<std::complex<T>> cropAndCombine(Image<std::complex<T>>& image) {
      auto& header = std::get<ISMRMRD::ImageHeader>(image);
      auto& data = std::get<hoNDArray<std::complex<T>>>(image);

      hoNDArray<std::complex<T>> m3 = hoNDArray<std::complex<T>>();

      std::vector<size_t> new_dimensions(3);
      new_dimensions[0] = data.get_size(0) >> 1;
      new_dimensions[1] = data.get_size(1);
      new_dimensions[2] = data.get_size(2);

      try {
          m3.create(new_dimensions);
      } 
      catch (std::runtime_error& err) {
          GEXCEPTION(err, "CropAndCombineGadget, failed to allocate new array\n");
      }

    size_t dimx     = m3.get_size(0);
    size_t dimx_old = data.get_size(0);

    size_t dimy = m3.get_size(1);
    size_t dimz = m3.get_size(2);

    size_t channels = data.get_size(3);

    std::complex<T>* d1 = data.get_data_ptr();
    std::complex<T>* d2 = m3.get_data_ptr();

    size_t img_block_old = dimx_old*dimy*dimz;

    for (size_t z = 0; z < dimz; z++) {
      for (size_t y = 0; y < dimy; y++) {
        for (size_t x = 0; x < dimx; x++) {
          float mag = 0;
          float phase = 0;
          size_t offset_1 = z*dimy*dimx_old+y*dimx_old+x+((dimx_old-dimx)>>1);
          size_t offset_2 = z*dimy*dimx+y*dimx+x;
          for (size_t c = 0; c < channels; c++) {
            float mag_tmp = norm(d1[offset_1 + c*img_block_old]);
            phase += mag_tmp*arg(d1[offset_1 + c*img_block_old]);
            mag += mag_tmp;
          }
          d2[offset_2] = std::polar(std::sqrt(mag),phase);
        }
      }
    }

    //Modify header to match
    header.matrix_size[0] = header.matrix_size[0]>>1;
    header.channels = 1;

    header.field_of_view[0] = header.field_of_view[0]/2;

    //Now add the new array to the outgoing message
    return Image<std::complex<T>>(std::get<ISMRMRD::ImageHeader>(image), std::move(m3), std::get<optional<ISMRMRD::MetaContainer>>(image));
  }
} // namespace

namespace Gadgetron {
  AnyImage CropAndCombineGadget::process_function(AnyImage image) const {
      return visit([&](auto& image) -> AnyImage { return cropAndCombine(image); }, image);
  }

  GADGETRON_GADGET_EXPORT(CropAndCombineGadget);
} // namespace Gadgetron
