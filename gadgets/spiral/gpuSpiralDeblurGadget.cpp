/*
Gadget for deblurring of real-time spiral images
Processes B0 field map and then calls the MFI operator to do deblurring based
on multi-frequency interpolation.
Input: Image data and B0 field map attached as reference
 - Image data in rbit_[0].data_.data_
 - Map data in rbit_[0].ref_->data_ (N dimension should be "set", 0th N-dim is TE0, 1st N-dim is TE1)
 Output: three image series
 - Deblurred image
 - Uncorrected image
 - B0 map
*/
#include "gpuSpiralDeblurGadget.h"
#include "vds.h"
#include "GPUTimer.h"
#include "GenericReconJob.h"
#include "b1_map.h"
#include "check_CUDA.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_fileio.h"
#include "cuNDArray_reductions.h"
#include "cuNDArray_utils.h"
#include "hoArmadillo.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_utils.h"
#include "sense_utilities.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "vector_td_utilities.h"
#include <algorithm>
#include <vector>

#include "MFIOperator.h"

#define ARMA_64BIT_WORD

#ifdef USE_ARMADILLO
	#include <armadillo>
#endif

#ifdef USE_OMP
    #include <omp.h>
#endif

using namespace std;
// using namespace Gadgetron;

namespace Gadgetron {

// Define desired precision
typedef float _real;
typedef complext<_real> _complext;
typedef reald<_real,2>::Type _reald2;
typedef cuNFFT_impl<_real,2> plan_type;

  gpuSpiralDeblurGadget::gpuSpiralDeblurGadget()
    : samples_to_skip_start_(0)
    , samples_to_skip_end_(0)
    , prepared_(false)
		, prepared_B0_(false)
  {
  }

  gpuSpiralDeblurGadget::~gpuSpiralDeblurGadget() {}

  int gpuSpiralDeblurGadget::process_config(const mrd::Header& header)
  {

    int number_of_devices = 0;
    if (cudaGetDeviceCount(&number_of_devices)!= cudaSuccess) {
      GDEBUG( "Error: unable to query number of CUDA devices.\n" );
      return GADGET_FAIL;
    }

    if (number_of_devices == 0) {
      GDEBUG( "Error: No available CUDA devices.\n" );
      return GADGET_FAIL;
    }

    device_number_ = deviceno.value();

    if (device_number_ >= number_of_devices) {
      GDEBUG("Adjusting device number from %d to %d\n", device_number_,  (device_number_%number_of_devices));
      device_number_ = (device_number_%number_of_devices);
    }

    if (cudaSetDevice(device_number_)!= cudaSuccess) {
      GDEBUG( "Error: unable to set CUDA device.\n" );
      return GADGET_FAIL;
    }

    cudaDeviceProp deviceProp;
    if( cudaGetDeviceProperties( &deviceProp, device_number_ ) != cudaSuccess) {
      GDEBUG( "Error: unable to query device properties.\n" );
      return GADGET_FAIL;
    }

	unsigned int warp_size = deviceProp.warpSize;

	auto& h = header;
	// Get the encoding space and trajectory description
    mrd::EncodingSpaceType e_space = h.encoding[0].encoded_space;
    mrd::EncodingSpaceType r_space = h.encoding[0].recon_space;
    mrd::EncodingLimitsType e_limits = h.encoding[0].encoding_limits;
    mrd::TrajectoryDescriptionType traj_desc;

    if (h.encoding[0].trajectory_description) {
      traj_desc = *h.encoding[0].trajectory_description;
    } else {
      GDEBUG("Trajectory description missing");
      return GADGET_FAIL;
    }

    if (traj_desc.identifier != "HargreavesVDS2000") {
      GDEBUG("Expected trajectory description identifier 'HargreavesVDS2000', not found.");
      return GADGET_FAIL;
    }

    long interleaves = -1;
    long fov_coefficients = -1;
    long sampling_time_ns = -1;
    double max_grad = -1.0;
    double max_slew = -1.0;
    double fov_coeff = -1.0;
    double kr_max = -1.0;


    for (std::vector<mrd::UserParameterLongType>::iterator i (traj_desc.user_parameter_long.begin()); i != traj_desc.user_parameter_long.end(); ++i) {
      if (i->name == "interleaves") {
        interleaves = i->value;
      } else if (i->name == "fov_coefficients") {
        fov_coefficients = i->value;
      } else if (i->name == "SamplingTime_ns") {
        sampling_time_ns = i->value;
      } else {
        GDEBUG("WARNING: unused trajectory parameter %s found\n", i->name.c_str());
      }
    }

    for (std::vector<mrd::UserParameterDoubleType>::iterator i (traj_desc.user_parameter_double.begin()); i != traj_desc.user_parameter_double.end(); ++i) {
      if (i->name == "MaxGradient_G_per_cm") {
	max_grad = i->value;
      } else if (i->name == "MaxSlewRate_G_per_cm_per_s") {
	max_slew = i->value;
      } else if (i->name == "FOVCoeff_1_cm") {
	fov_coeff = i->value;
      } else if (i->name == "krmax_per_cm") {
	kr_max= i->value;
      } else {
	GDEBUG("WARNING: unused trajectory parameter %s found\n", i->name.c_str());
      }
    }

	krmax_ = kr_max/2;
	sample_time = (1.0*sampling_time_ns) * 1e-9;
    gmax_ = max_grad;
    smax_ = max_slew;
    krmax_ = kr_max*1.1;
    fov_ = fov_coeff;

    // Determine reconstruction matrix sizes
    //
	fov_vec_.push_back(r_space.field_of_view_mm.x);
    fov_vec_.push_back(r_space.field_of_view_mm.y);
    fov_vec_.push_back(r_space.field_of_view_mm.z);

	kernel_width_ = buffer_convolution_kernel_width.value();
    oversampling_factor_ = buffer_convolution_oversampling_factor.value();

    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrix_size.x*reconstruction_os_factor_x.value()))+warp_size-1)/warp_size)*warp_size);
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrix_size.y*reconstruction_os_factor_y.value()))+warp_size-1)/warp_size)*warp_size);

    image_dimensions_recon_os_ = uint64d2
      (((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
       ((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);

	deblur_ = do_deblurring.value();

	return GADGET_OK;
  }



    int gpuSpiralDeblurGadget::process(Gadgetron::GadgetContainerMessage< mrd::ReconData >* m1)
    {
		//Image data in rbit_[0].data_.data_
		//Map data in rbit_[0].ref_->data_ (N dimension should be "set", 0th N-dim is TE0, 1st N-dim is TE1)

		mrd::ReconData* recon_bit_ = m1->getObjectPtr();

		// Allocate various counters if they are NULL
		if( !image_counter_.get() ){
			image_counter_ = boost::shared_array<long>(new long[1]);
			image_counter_[0] = 0;
		}

		//Get data from ReconDataBuffer if exists
		hoNDArray<std::complex<float>> host_data = recon_bit_->buffers[0].data.data;
		hoNDArray<std::complex<float>> B0_data;
		if(recon_bit_->buffers[0].ref){
					B0_data = recon_bit_->buffers[0].ref->data;
		}

		//Set up image containers, trajectory, and NFFT_plan - if not already prepared
		mrd::AcquisitionHeader& curr_header = recon_bit_->buffers[0].data.headers(0,0,0,0,0);
		if(!prepared_ && host_data.get_size(0) > 0){ //TODO: move to process_config?
		    Prepare_Plan(recon_bit_->buffers[0].data);
		}

		//Set up B0 map containers, trajectory, and NFFT_plan - if not already prepared
		if(!prepared_B0_ && recon_bit_->buffers[0].ref){
				Prepare_B0_Plan(*recon_bit_->buffers[0].ref);
		}

		//If there is reference data then we need to re-compute the B0 map
		if(recon_bit_->buffers[0].ref){
				Calc_B0Map(B0_data, &B0_map);
		}

		//If there is image data (i.e. reconbit isn't only map data), then we reconstruct and deblur
		if(host_data.get_number_of_elements() > 0){

			cuNDArray<complext<float>> gpu_data((hoNDArray<float_complext>&)host_data);
			nfft_plan_->compute( gpu_data, image, &gpu_weights, NFFT_comp_mode::BACKWARDS_NC2C );
			csm_ = estimate_b1_map<float,2>( image );
			csm_mult_MH<float,2>(&image, &reg_image, &csm_);
			host_image = *(reg_image.to_host());

			gpu_traj = host_traj;
			if(MFI.prepare(nfft_plan_, gpu_traj, host_data.get_dimensions(), image_dimensions_recon_,
                                image_dimensions_recon_os_, sample_time, .001)){
				output_image = MFI.MFI_apply(host_image, B0_map);
			}
			host_image = *(reg_image.to_host()); //calling MFI_apply corrupts host_image, Recall it from GPU;

			// Make sure the output images have the correct dimensions [X Y Z CHA]
			auto host_image_dimensions = host_image.dimensions();
			host_image_dimensions.push_back(1);
			host_image_dimensions.push_back(1);

			// Queue deblurred im
			GadgetContainerMessage<mrd::ImageComplexFloat> *deblurred = new GadgetContainerMessage<mrd::ImageComplexFloat>();
			deblurred->getObjectPtr()->head = get_image_header(curr_header, 1);
			deblurred->getObjectPtr()->data.create(host_image_dimensions);
			memcpy(deblurred->getObjectPtr()->data.get_data_ptr(), output_image.get_data_ptr(), output_image.get_number_of_elements()*sizeof(std::complex<float>));
			if (this->next()->putq(deblurred) < 0) {
			  GDEBUG("Failed to put job on queue.\n");
			  deblurred->release();
			  return GADGET_FAIL;
			}

			// Queue original im
			GadgetContainerMessage<mrd::ImageComplexFloat> *orig = new GadgetContainerMessage<mrd::ImageComplexFloat>();
			orig->getObjectPtr()->head = get_image_header(curr_header, 10);
			orig->getObjectPtr()->data.create(host_image_dimensions);
			memcpy(orig->getObjectPtr()->data.get_data_ptr(), host_image.get_data_ptr(), host_image.get_number_of_elements()*sizeof(std::complex<float>));
			if (this->next()->putq(orig) < 0) {
			  GDEBUG("Failed to put job on queue.\n");
			  orig->release();
			  return GADGET_FAIL;
			}

			// Queue B0 map
			hoNDArray< std::complex<float> >B0_image;
			B0_image.create(B0_map.dimensions());
			for(int i = 0; i<B0_map.get_number_of_elements(); i++){
				B0_image[i] = std::complex<float>(B0_map[i],0.0);
			}
			GadgetContainerMessage<mrd::ImageComplexFloat> *b0map = new GadgetContainerMessage<mrd::ImageComplexFloat>();
			b0map->getObjectPtr()->head = get_image_header(curr_header, 20);

			auto B0_dimensions = B0_image.dimensions();
			B0_dimensions.push_back(1);
			B0_dimensions.push_back(1);
			b0map->getObjectPtr()->data.create(B0_dimensions);
			memcpy(b0map->getObjectPtr()->data.get_data_ptr(), B0_image.get_data_ptr(), B0_image.get_number_of_elements()*sizeof(std::complex<float>));
			if (this->next()->putq(b0map) < 0) {
			  GDEBUG("Failed to put job on queue.\n");
			  b0map->release();
			  return GADGET_FAIL;
			}

			image_counter_[0]++; //increment image counter
		}

		m1->release();
		return GADGET_OK;
	}

	mrd::ImageHeader gpuSpiralDeblurGadget::get_image_header(mrd::AcquisitionHeader& curr_header, int series_index)
	{
		// Prepare an image header for this frame
		mrd::ImageHeader header{};

		header.field_of_view[0] = fov_vec_[0];
		header.field_of_view[1] = fov_vec_[1];
		header.field_of_view[2] = fov_vec_[2];

		header.slice = curr_header.idx.slice;
		header.set = curr_header.idx.set;

		header.acquisition_time_stamp = curr_header.acquisition_time_stamp;
		header.physiology_time_stamp = curr_header.physiology_time_stamp;

		header.position = curr_header.position;
		header.col_dir = curr_header.read_dir;
		header.line_dir = curr_header.phase_dir;
		header.slice_dir = curr_header.slice_dir;
		header.patient_table_position = curr_header.patient_table_position;

		header.image_index = image_counter_[0];
		header.image_series_index = series_index;

		return header;
	}

	void gpuSpiralDeblurGadget::Calc_B0Map(hoNDArray<std::complex<float>>& B0_data, hoNDArray<float>* B0_map){
		size_t R0 = B0_data.get_size(0);
		size_t E1 = B0_data.get_size(1);
		size_t E2 = B0_data.get_size(2);
		size_t CHA = B0_data.get_size(3);
		size_t N = B0_data.get_size(4);

		//Low pass filter map data since conjugate phase reconstruction assumes smoothly varying B0 map
		for(int i =0; i < B0_data.get_number_of_elements(); i++){
			B0_data[i] *= exp(-.5*pow((i%R0)/(R0/3),2.) );
		}
		//Recon B0 images
		hoNDArray<std::complex<float>>* B0_data_1 = new hoNDArray<std::complex<float>>(R0,E1,E2,CHA,B0_data.get_data_ptr());
		cuNDArray<complext<float>> gpu_B0_data((hoNDArray<float_complext>&)(*B0_data_1));
		nfft_plan_B0_->compute( gpu_B0_data, image, &gpu_weights_B0, NFFT_comp_mode::BACKWARDS_NC2C );
		csm_ = estimate_b1_map<float,2>( image );
		csm_mult_MH<float,2>(&image, &reg_image, &csm_);
		hoNDArray<complext<float>> B0_temp_0 = *reg_image.to_host();
		hoNDArray<std::complex<float>>* B0_data_2 = new hoNDArray<std::complex<float>>(R0,E1,E2,CHA,B0_data.get_data_ptr()+R0*E1*E2*CHA);//Start at index N = 1
		gpu_B0_data = *((hoNDArray<float_complext>*)B0_data_2);
		nfft_plan_B0_->compute( gpu_B0_data, image, &gpu_weights_B0, NFFT_comp_mode::BACKWARDS_NC2C );
		csm_mult_MH<float,2>(&image, &reg_image, &csm_);
		hoNDArray<complext<float>> B0_temp_1 = *reg_image.to_host();
		//Compute map
		B0_map->clear();
		B0_map->create(B0_temp_0.dimensions());
		B0_map->fill(0.0);
		auto map_ptr = B0_map->get_data_ptr();
		for (int i = 0; i < B0_temp_0.get_number_of_elements(); i++) {
				if(abs(B0_temp_1[i]) > 3){ //assumes SNR units and filters out noise < 3 std
					map_ptr[i] = Gadgetron::arg(B0_temp_0[i]*Gadgetron::conj(B0_temp_1[i]))/( 2*M_PI*.001 );//delTE = 1 ms
				}
		}

	}

	void gpuSpiralDeblurGadget::Prepare_Plan(mrd::ReconBuffer& data){
			size_t R0 = data.data.get_size(0);
			size_t E1 = data.data.get_size(1);
			size_t E2 = data.data.get_size(2);
			size_t CHA = data.data.get_size(3);

			//Setup image arrays
			image_dimensions_recon_.push_back(CHA);
			image.create(image_dimensions_recon_);
			image_dimensions_recon_.pop_back();
			reg_image.create(image_dimensions_recon_);

			host_traj.create(R0*E1);
			host_weights.create(R0*E1);

			//Trajectory should be attached, but if it isn't we call calc_vds
			mrd::AcquisitionHeader& curr_header = data.headers(0,0,0,0,0);
			size_t trajectory_dimensions = data.trajectory.get_size(1);
			if (trajectory_dimensions != 3) {
				//Setup calc_vds parameters
				int     nfov   = 1;
				int     ngmax  = 1e7;       /*  maximum number of gradient samples      */
				double  *xgrad;             /*  x-component of gradient.                */
				double  *ygrad;             /*  y-component of gradient.                */
				double  *x_trajectory;
				double  *y_trajectory;
				double  *weighting;
				int     ngrad;
				//Map trajecotry is different, parameters defined in user_floats
				calc_vds(smax_,gmax_,sample_time,sample_time,E1,&fov_,nfov,krmax_,ngmax,&xgrad,&ygrad,&ngrad);
				calc_traj(xgrad, ygrad, R0, E1, sample_time, krmax_, &x_trajectory, &y_trajectory, &weighting);

				for (int i = 0; i < (R0*E1); i++) {
					host_traj[i]   = floatd2(-x_trajectory[i]/2.,-y_trajectory[i]/2.);
					host_weights[i] = weighting[i];
				}

				delete [] xgrad;
				delete [] ygrad;
				delete [] x_trajectory;
				delete [] y_trajectory;
				delete [] weighting;
			}
			//If traj attached:
			else{
				for (int i = 0; i < (R0*E1); i++) {
					auto trajectory = data.trajectory.get_data_ptr();
					host_traj[i]   = floatd2(trajectory[i*3],trajectory[i*3+1]);
					host_weights[i] = trajectory[i*3+2];
				}
			}

		//upload to gpu
		gpu_traj = host_traj;
		gpu_weights = host_weights;
		//pre-process
		nfft_plan_ = NFFT<cuNDArray,_real,2>::make_plan( from_std_vector<size_t,2>(image_dimensions_recon_), image_dimensions_recon_os_, kernel_width_ );

		nfft_plan_->preprocess(gpu_traj, NFFT_prep_mode::ALL);
		prepared_ = true;
	}

	void gpuSpiralDeblurGadget::Prepare_B0_Plan(mrd::ReconBuffer& data){
		mrd::AcquisitionHeader& B0_header = data.headers(0,0,0,0,0);
		size_t R0 = data.data.get_size(0);
		size_t E1 = data.data.get_size(1);

		B0_traj.create(R0*E1);
		B0_weights.create(R0*E1);

		float krmaxB0_ = 2.*(B0_header.user_float[4])/10000.; //TODO: B0 acquisition info currently embedded into user parameters. In the future this should be in Encoding[1]. Requires amending XSL.
		size_t trajectory_dimensions = data.trajectory.get_size(1);
		if (trajectory_dimensions != 3) {
			//Setup calc_vds parameters
			const int     nfov   = 2;
			int     ngmax  =1e7;       /*  maximum number of gradient samples      */
			double  *xgrad;             /*  x-component of gradient.                */
			double  *ygrad;             /*  y-component of gradient.                */
			double  *x_trajectory;
			double  *y_trajectory;
			double  *weighting;
			int     ngrad;
			//Map trajecotry is different, parameters defined in user_floats
			GINFO("fov %f\n", B0_header.user_float[5]);
			GINFO("smax %f\n", 3*B0_header.user_float[3]/10.);
			GINFO("gmax %f\n", B0_header.user_float[1]/10.);
			GINFO("krmax %f\n", 2*((B0_header.user_float[4])/10000.));
			GINFO("sampling_time %f\n", sample_time);
			double fov2_[nfov] = {B0_header.user_float[5], -1/(1.1*2*((B0_header.user_float[4])/10000.))*B0_header.user_float[5]};
			calc_vds(3.*((B0_header.user_float[3])/10.),(B0_header.user_float[1])/10.,sample_time,sample_time,E1,&fov2_[0],nfov,2.*((B0_header.user_float[4])/10000.),ngmax,&xgrad,&ygrad,&ngrad);
			calc_traj(xgrad, ygrad, R0, E1, sample_time, krmaxB0_, &x_trajectory, &y_trajectory, &weighting);

			for (int i = 0; i < (R0*E1); i++) {
				B0_traj[i]   = floatd2(-x_trajectory[i]/(2.),-y_trajectory[i]/(2.));
				B0_weights[i] = weighting[i];
			}

			delete [] xgrad;
			delete [] ygrad;
			delete [] x_trajectory;
			delete [] y_trajectory;
			delete [] weighting;
		}
		else{
			for (int i = 0; i < (R0*E1); i++) {
				auto B0_trajectory = data.trajectory.get_data_ptr();
				B0_traj[i]   = floatd2(B0_trajectory[i*3],B0_trajectory[i*3+1]);
				B0_weights[i] = B0_trajectory[i*3+2];
			}
		}

		//upload to gpu
		gpu_traj = B0_traj;
		gpu_weights_B0 = B0_weights;

		//pre-process
		nfft_plan_B0_ = NFFT<cuNDArray,_real,2>::make_plan( from_std_vector<size_t,2>(image_dimensions_recon_), image_dimensions_recon_os_, kernel_width_ );

		nfft_plan_B0_->preprocess(gpu_traj, NFFT_prep_mode::NC2C);
		prepared_B0_= true;
	}

  GADGET_FACTORY_DECLARE(gpuSpiralDeblurGadget)
}
