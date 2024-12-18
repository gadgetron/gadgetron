#include "DeviceChannelSplitterGadget.h"

 //This is needed for things such as data role, which should NOT be defined in gtPlus
#include "mri_core_def.h"

namespace Gadgetron{

template <typename T>
int DeviceChannelSplitterGadget<T>::process(GadgetContainerMessage<mrd::Image<T>>* m1)
{

  auto& head = m1->getObjectPtr()->head;
  auto& data = m1->getObjectPtr()->data;

  size_t array_channels = data.get_size(data.get_number_of_dimensions() - 1);
  size_t dim_x = data.get_size(0);
  size_t dim_y = data.get_size(1);
  size_t dim_z = data.get_size(2);
  size_t image_elements = dim_x * dim_y * dim_z;

  for (int i = 0; i < array_channels; i++) {
    auto im1 = new GadgetContainerMessage<mrd::Image<T>>();

    im1->getObjectPtr()->head = head;
    im1->getObjectPtr()->data.create(dim_x, dim_y, dim_z, 1);

    memcpy(im1->getObjectPtr()->data.get_data_ptr(), data.get_data_ptr() + i * image_elements, sizeof(T)*image_elements);

    auto& meta_out = im1->getObjectPtr()->meta;
    if (i == 0) {
      meta_out[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_IRT_IMAGE };
    } else {
      meta_out[GADGETRON_DATA_ROLE] = { GADGETRON_IMAGE_IRT_DEVICE };
      meta_out[GADGETRON_IMAGE_CUR_DEVICE_CHA] = { (long)i };

    }
    meta_out[GADGETRON_DATA_ROLE].push_back(GADGETRON_IMAGE_INTENSITY_UNCHANGED);

    if (array_channels > 1) {
      meta_out[GADGETRON_IMAGE_NUM_DEVICE_CHA] = { (long)(array_channels-1) };
    } else {
      meta_out[GADGETRON_IMAGE_NUM_DEVICE_CHA] = { (long)(-1) };
    }

    if (this->next()->putq(im1) == -1) {
      im1->release();
      GERROR("DeviceChannelSplitterGadget::process, passing data on to next gadget\n");
      return -1;
    }
  }

  // We are done with the original data
  m1->release();

  return GADGET_OK;
}

// Declare factories for the various template instances
GADGET_FACTORY_DECLARE(DeviceChannelSplitterGadgetFLOAT);
GADGET_FACTORY_DECLARE(DeviceChannelSplitterGadgetUSHORT);
GADGET_FACTORY_DECLARE(DeviceChannelSplitterGadgetCPLX);

} // namespace Gadgetron
