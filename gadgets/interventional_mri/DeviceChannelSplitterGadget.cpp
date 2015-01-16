#include "DeviceChannelSplitterGadget.h"
#include "ismrmrd/meta.h"

 //This is needed for things such as data role, which should NOT be defined in gtPlus
#include "mri_core_def.h"

namespace Gadgetron{

template <typename T>
int DeviceChannelSplitterGadget<T>
::process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
	  GadgetContainerMessage< hoNDArray< T > >* m2)
{
  
  //Some consistency checking
  unsigned int header_channels = m1->getObjectPtr()->channels;
  unsigned int array_channels = m2->getObjectPtr()->get_size(m2->getObjectPtr()->get_number_of_dimensions()-1);
  unsigned int dim_x = m2->getObjectPtr()->get_size(0);
  unsigned int dim_y = m2->getObjectPtr()->get_size(1);
  unsigned int dim_z = m2->getObjectPtr()->get_size(2);
  size_t image_elements = dim_x*dim_y*dim_z;

  if (header_channels != array_channels) {
    GDEBUG("Inconsistent number of header channels (%d) and array channels (%d)\n", header_channels, array_channels);
    m1->release();
    return GADGET_FAIL;
  }
  

  for (int i = 0; i < array_channels; i++) {
    

    GadgetContainerMessage<ISMRMRD::ImageHeader>* im1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
    *(im1->getObjectPtr()) = *(m1->getObjectPtr());
    im1->getObjectPtr()->channels = 1;
    
    /*
    GDEBUG("Image with matrix (cha=%d): %d, %d, %d and fov %f, %f, %f\n", 
		  i,
		  im1->getObjectPtr()->matrix_size[0], 
		  im1->getObjectPtr()->matrix_size[1], 
		  im1->getObjectPtr()->matrix_size[2],
		  im1->getObjectPtr()->field_of_view[0],
		  im1->getObjectPtr()->field_of_view[1],
		  im1->getObjectPtr()->field_of_view[2]);

    */

    GadgetContainerMessage< hoNDArray< T > >* im2 = new GadgetContainerMessage< hoNDArray< T > >();
    im2->getObjectPtr()->create(dim_x,dim_y,dim_z,1);
    memcpy(im2->getObjectPtr()->get_data_ptr(), m2->getObjectPtr()->get_data_ptr() + i*image_elements, sizeof(T)*image_elements);
    
    im1->cont(im2);
    
    Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>* im3 = new Gadgetron::GadgetContainerMessage<ISMRMRD::MetaContainer>();
    if (i == 0) {
      im3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_IRT_IMAGE);
    } else {
      im3->getObjectPtr()->set(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_IRT_DEVICE);
      im3->getObjectPtr()->set(GADGETRON_IMAGE_CUR_DEVICE_CHA, (long)i);

    }
    im3->getObjectPtr()->append(GADGETRON_DATA_ROLE, GADGETRON_IMAGE_INTENSITY_UNCHANGED);

    if (array_channels > 1) {
      im3->getObjectPtr()->set(GADGETRON_IMAGE_NUM_DEVICE_CHA, (long)(array_channels-1));
    } else {
      im3->getObjectPtr()->set(GADGETRON_IMAGE_NUM_DEVICE_CHA, (long)(-1));
    }

    im2->cont(im3);

    if (this->next()->putq(im1) == -1) {
      m1->release();
      GERROR("DeviceChannelSplitterGadget::process, passing data on to next gadget\n");
      return -1;
    }
  }

  //We are done with the original data
  m1->release();


  return GADGET_OK;
}

//Declare factories for the various template instances
GADGET_FACTORY_DECLARE(DeviceChannelSplitterGadgetFLOAT);
GADGET_FACTORY_DECLARE(DeviceChannelSplitterGadgetUSHORT);
GADGET_FACTORY_DECLARE(DeviceChannelSplitterGadgetCPLX);

}
