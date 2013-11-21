#include "../mri_core/GadgetIsmrmrdReadWrite.h"
#include "Gadgetron.h"
#include "GrappaGadget.h"
#include "GrappaUnmixingGadget.h"
#include "GadgetIsmrmrdReadWrite.h"

#include <ace/OS_NS_stdlib.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

namespace Gadgetron{

  GrappaGadget::GrappaGadget()
  : image_counter_(0)
  , image_series_(0)
  , first_call_(true)
  , target_coils_(0)
  {
  }

  GrappaGadget::~GrappaGadget()
  {
    for (unsigned int i = 0; i < buffers_.size(); i++) {
      if (buffers_[i]) delete buffers_[i];
      buffers_[i] = 0;


      if (image_data_[i]) {
	image_data_[i]->release();
	image_data_[i] = 0;
      }
    }
  }

  int GrappaGadget::close(unsigned long flags) {
    int ret = Gadget::close(flags);
    GADGET_DEBUG1("Shutting down GRAPPA Gadget\n");

    if (weights_calculator_.close(flags) < 0) {
      GADGET_DEBUG1("Failed to close down weights calculator\n");
      return GADGET_FAIL;
    }

    return ret;
  }

  int GrappaGadget::process_config(ACE_Message_Block* mb)
  {
    boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

    ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
    if (e_seq.size() != 1) {
      GADGET_DEBUG2("Number of encoding spaces: %d\n", e_seq.size());
      GADGET_DEBUG1("This Gadget only supports one encoding space\n");
      return GADGET_FAIL;
    }

    ISMRMRD::encodingSpaceType e_space = (*e_seq.begin()).encodedSpace();
    ISMRMRD::encodingSpaceType r_space = (*e_seq.begin()).reconSpace();
    ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

    unsigned int slices = e_limits.slice().present() ? e_limits.slice().get().maximum() + 1 : 1;
    dimensions_.push_back(e_space.matrixSize().x());
    dimensions_.push_back(e_space.matrixSize().y());
    dimensions_.push_back(e_space.matrixSize().z());
    dimensions_.push_back((cfg->acquisitionSystemInformation().present() && cfg->acquisitionSystemInformation().get().receiverChannels().present()) ?
			  cfg->acquisitionSystemInformation().get().receiverChannels().get() : 1);
    dimensions_.push_back(slices);

    fov_.push_back(r_space.fieldOfView_mm().x());
    fov_.push_back(r_space.fieldOfView_mm().y());
    fov_.push_back(r_space.fieldOfView_mm().z());

    line_offset_ = (dimensions_[1]>>1)-e_limits.kspace_encoding_step_1().get().center();

    return GADGET_OK;
  }


  int GrappaGadget::initial_setup()
  {

    GADGET_DEBUG2("Dimensions %d, %d, %d, %d, %d\n", dimensions_[0], dimensions_[1], dimensions_[2], dimensions_[3], dimensions_[4]);

    image_dimensions_.push_back(dimensions_[0] / 2); //TODO: fix this in general
    image_dimensions_.push_back(dimensions_[1]);
    image_dimensions_.push_back(dimensions_[2]);
    image_dimensions_.push_back(dimensions_[3]);


    weights_ = std::vector< boost::shared_ptr<GrappaWeights<float> > >(dimensions_[4]);

    buffers_ = std::vector<GrappaCalibrationBuffer* >(dimensions_[4],0);
    time_stamps_ = std::vector<ACE_UINT32>(dimensions_[4],0);

    //Let's figure out the number of target coils
    target_coils_ = this->get_int_value("target_coils");
    if ((target_coils_ <= 0) || (target_coils_ > dimensions_[3])) {
      target_coils_ = dimensions_[3];
    }

    GADGET_DEBUG2("Running GRAPPA recon with %d source channels and %d target channels\n", dimensions_[3], target_coils_);

    weights_calculator_.set_number_of_target_coils(target_coils_);

    //Let's figure out if we have channels that are supposed to be uncombined
    boost::shared_ptr<std::string> uncomb_str = this->get_string_value("uncombined_channels");
    std::vector<std::string> uncomb;
    boost::split(uncomb, *uncomb_str, boost::is_any_of(","));
    for (unsigned int i = 0; i < uncomb.size(); i++) {
      std::string ch = boost::algorithm::trim_copy(uncomb[i]);
      if (ch.size() > 0) {
	unsigned int channel_id = static_cast<unsigned int>(ACE_OS::atoi(ch.c_str()));
	weights_calculator_.add_uncombined_channel(channel_id);
      }
    }

    for (unsigned int i = 0; i < buffers_.size(); i++) {
      weights_[i] = boost::shared_ptr<GrappaWeights<float> >(new GrappaWeights<float>());

      //Let's set some default GRAPPA weights, so that we have something to work with the first couple of frames.
      /*
	std::vector<unsigned int> wdims = image_dimensions_;
	if (weights_calculator_.get_number_of_uncombined_channels()) {
	wdims.push_back(weights_calculator_.get_number_of_uncombined_channels()+1);
	}

	hoNDArray< std::complex<float> > tmp_w;
	if (!tmp_w.create(&wdims)) {
	GADGET_DEBUG1("Unable to create temporary array with dimensions\n");
	return GADGET_FAIL;
	}
	tmp_w.clear(std::complex<float>(1.0,0));
	weights_[i]->update(&tmp_w);
      */

      buffers_[i] = new GrappaCalibrationBuffer(image_dimensions_,
						weights_[i],
						&weights_calculator_);
    }


    if (weights_calculator_.open() < 0) {
      GADGET_DEBUG1("Failed to open GrappaWeightsCalculator\n");
      return GADGET_FAIL;
    }

    image_data_ = std::vector< GadgetContainerMessage< hoNDArray< std::complex<float> > >* >(dimensions_[4],0);
    for (unsigned int i = 0; i < image_data_.size(); i++) {
      if (create_image_buffer(i) != GADGET_OK) {
	GADGET_DEBUG1("Unable to create image buffers");
	return GADGET_FAIL;
      }
    }

    image_series_ = this->get_int_value("image_series");

    return GADGET_OK;
  }


  int GrappaGadget::
  process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader>* m1,
	  GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
  {

    if (first_call_) {
      if (m1->getObjectPtr()->active_channels != dimensions_[3]) {
	GADGET_DEBUG1("Detected coil number change. Maybe due to upstream channel reduction\n");
	dimensions_[3] = m1->getObjectPtr()->active_channels;
      }

      if (initial_setup() != GADGET_OK) {
	GADGET_DEBUG1("Initial Setup Failed\n");
	m1->release();
	return GADGET_FAIL;
      }
      first_call_ = false;
    }

    ISMRMRD::AcquisitionHeader* acq_head = m1->getObjectPtr();

    unsigned int samples =  acq_head->number_of_samples;
    unsigned int line = acq_head->idx.kspace_encode_step_1 + line_offset_;
    unsigned int partition = acq_head->idx.kspace_encode_step_2;
    unsigned int slice = acq_head->idx.slice;

    if (samples != image_dimensions_[0]) {
      GADGET_DEBUG1("GrappaGadget: wrong number of samples received\n");
      return GADGET_FAIL;
    }

    if (slice >= image_data_.size()) {
      GADGET_DEBUG1("Invalid slice number received\n");
      return GADGET_FAIL;
    }

    if (!image_data_[0]) {
      if (create_image_buffer(slice) != GADGET_OK) {
	GADGET_DEBUG1("Failed to allocate new slice buffer\n");
	return GADGET_FAIL;
      }
    }

    std::complex<float>* b = image_data_[slice]->getObjectPtr()->get_data_ptr();
    std::complex<float>* d = m2->getObjectPtr()->get_data_ptr();

    size_t offset= 0;
    //Copy the data for all the channels
    for (int c = 0; c < m1->getObjectPtr()->active_channels; c++) {
      offset =
	c*image_dimensions_[0]*image_dimensions_[1]*image_dimensions_[2] +
	partition*image_dimensions_[0]*image_dimensions_[1] +
	line*image_dimensions_[0];

      memcpy(b+offset,d+c*samples,sizeof(std::complex<float>)*samples);
    }


    bool is_last_scan_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_LAST_IN_SLICE).isSet(m1->getObjectPtr()->flags);

    bool is_first_scan_in_slice = ISMRMRD::FlagBit(ISMRMRD::ACQ_FIRST_IN_SLICE).isSet(m1->getObjectPtr()->flags);

    if (is_first_scan_in_slice) {
      time_stamps_[slice] = m1->getObjectPtr()->acquisition_time_stamp;
    }

    if (is_last_scan_in_slice) {

      GadgetContainerMessage<GrappaUnmixingJob>* cm0 =
	new GadgetContainerMessage<GrappaUnmixingJob>();

      GadgetContainerMessage<ISMRMRD::ImageHeader>* cm1 =
	new GadgetContainerMessage<ISMRMRD::ImageHeader>();


      /*
	GadgetContainerMessage< hoNDArray<std::complex<float> > >* cm2 =
	new GadgetContainerMessage< hoNDArray<std::complex<float> > >();

	std::vector<unsigned int> combined_dims(3,0);
	combined_dims[0] = image_dimensions_[0];
	combined_dims[1] = image_dimensions_[1];
	combined_dims[2] = image_dimensions_[2];

	if (weights_calculator_.get_number_of_uncombined_channels()) {
	combined_dims.push_back(weights_calculator_.get_number_of_uncombined_channels()+1);
	}

	if (!cm2->getObjectPtr()->create(&combined_dims)) {
	GADGET_DEBUG1("Unable to create combined image array\n");
	return GADGET_FAIL;
	}

	cm1->cont(cm2);
      */

      cm1->getObjectPtr()->matrix_size[0] = image_dimensions_[0];
      cm1->getObjectPtr()->matrix_size[1] = image_dimensions_[1];
      cm1->getObjectPtr()->matrix_size[2] = image_dimensions_[2];

      cm1->getObjectPtr()->field_of_view[0] = fov_[0];
      cm1->getObjectPtr()->field_of_view[1] = fov_[1];
      cm1->getObjectPtr()->field_of_view[2] = fov_[2];

      cm1->getObjectPtr()->channels       = 1+weights_calculator_.get_number_of_uncombined_channels();
      cm1->getObjectPtr()->slice              = m1->getObjectPtr()->idx.slice;
      cm1->getObjectPtr()->acquisition_time_stamp         = time_stamps_[slice];

      memcpy(cm1->getObjectPtr()->position,m1->getObjectPtr()->position,
	     sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->read_dir,m1->getObjectPtr()->read_dir,
	     sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->phase_dir,m1->getObjectPtr()->phase_dir,
	     sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->slice_dir,m1->getObjectPtr()->slice_dir,
	     sizeof(float)*3);

      memcpy(cm1->getObjectPtr()->patient_table_position,m1->getObjectPtr()->patient_table_position, sizeof(float)*3);

      cm1->getObjectPtr()->image_index = ++image_counter_;
      cm1->getObjectPtr()->image_series_index = image_series_;


      cm0->getObjectPtr()->weights_ = weights_[slice];
      cm0->cont(cm1);
      cm1->cont(image_data_[slice]);

      image_data_[slice] = 0;
      if (create_image_buffer(slice) != GADGET_OK) {
	GADGET_DEBUG1("Failed to create image buffer");
	return GADGET_FAIL;
      }

      if (this->next()->putq(cm0) < 0) {
	GADGET_DEBUG1("Failed to pass image on to next Gadget in chain\n");
	return GADGET_FAIL;
      }

      /*
	hoFFT<float>::instance()->ifft(image_data_[slice]->getObjectPtr(),0);
	hoFFT<float>::instance()->ifft(image_data_[slice]->getObjectPtr(),1);
	hoFFT<float>::instance()->ifft(image_data_[slice]->getObjectPtr(),2);

	//apply weights
	float scale_factor = (dimensions_[0] *dimensions_[1] *dimensions_[0] *dimensions_[1])/10;

	int appl_result = weights_[slice]->apply(image_data_[slice]->getObjectPtr(), cm2->getObjectPtr(), scale_factor);
	if (appl_result < 0) {
	GADGET_DEBUG2("Failed to apply GRAPPA weights: error code %d\n", appl_result);
	return GADGET_FAIL;
	}

	if (this->next()->putq(cm1) < 0) {
	GADGET_DEBUG1("Failed to pass image on to next Gadget in chain\n");
	return GADGET_FAIL;
	}
	image_data_[slice]->getObjectPtr()->clear(std::complex<float>(0.0f,0.0f));
      */
    }

    if (buffers_[slice]->add_data(m1->getObjectPtr(),m2->getObjectPtr()) < 0) {
      GADGET_DEBUG1("Failed to add incoming data to grappa calibration buffer\n");
      return GADGET_FAIL;
    }

    m1->release();
    return GADGET_OK;
  }


  int GrappaGadget::create_image_buffer(unsigned int slice)
  {
    if (slice >= image_data_.size()) {
      return GADGET_FAIL;
    }

    if (image_data_[slice] != 0) {
      image_data_[slice]->release();
      image_data_[slice] = 0;
    }

    image_data_[slice] = new GadgetContainerMessage< hoNDArray< std::complex<float> > >();
    try{ image_data_[slice]->getObjectPtr()->create(&image_dimensions_);}
    catch (std::runtime_error &err){
      GADGET_DEBUG_EXCEPTION(err,"Unable to create image buffers");
      return GADGET_FAIL;
    }

    std::fill(image_data_[slice]->getObjectPtr()->get_data_ptr(),
	      image_data_[slice]->getObjectPtr()->get_data_ptr()+image_data_[slice]->getObjectPtr()->get_number_of_elements(),
	      std::complex<float>(0.0f,0.0f));

    return GADGET_OK;

  }

  GADGET_FACTORY_DECLARE(GrappaGadget)
}
