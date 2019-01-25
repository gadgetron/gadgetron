#include "gpuSpiralDeblurGadget.h"
#include "GenericReconJob.h"
#include "cuNDArray_utils.h"
#include "hoNDArray_utils.h"
#include "cuNDArray_elemwise.h"
#include "cuNDArray_reductions.h"
#include "vector_td_utilities.h"
#include "hoNDArray_fileio.h"
#include "cuNDArray_fileio.h"
#include "vector_td.h"
#include "vector_td_operators.h"
#include "sense_utilities.h"
#include "check_CUDA.h"
#include "b1_map.h"
#include "hoArmadillo.h"
#include "GPUTimer.h"
#include "vds.h"
#include "ismrmrd/xml.h"
#include <algorithm>
#include <vector>

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

  int gpuSpiralDeblurGadget::process_config(ACE_Message_Block* mb)
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

    // Start parsing the ISMRMRD XML header
    //

    ISMRMRD::IsmrmrdHeader h;
    ISMRMRD::deserialize(mb->rd_ptr(),h);

	// Get the encoding space and trajectory description
    ISMRMRD::EncodingSpace e_space = h.encoding[0].encodedSpace;
    ISMRMRD::EncodingSpace r_space = h.encoding[0].reconSpace;
    ISMRMRD::EncodingLimits e_limits = h.encoding[0].encodingLimits;
    ISMRMRD::TrajectoryDescription traj_desc;

    if (h.encoding[0].trajectoryDescription) {
      traj_desc = *h.encoding[0].trajectoryDescription;
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
    
    
    for (std::vector<ISMRMRD::UserParameterLong>::iterator i (traj_desc.userParameterLong.begin()); i != traj_desc.userParameterLong.end(); ++i) {
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

    for (std::vector<ISMRMRD::UserParameterDouble>::iterator i (traj_desc.userParameterDouble.begin()); i != traj_desc.userParameterDouble.end(); ++i) {
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
	fov_vec_.push_back(r_space.fieldOfView_mm.x);
    fov_vec_.push_back(r_space.fieldOfView_mm.y);
    fov_vec_.push_back(r_space.fieldOfView_mm.z);

	kernel_width_ = buffer_convolution_kernel_width.value();
    oversampling_factor_ = buffer_convolution_oversampling_factor.value();
    
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize.x*reconstruction_os_factor_x.value()))+warp_size-1)/warp_size)*warp_size);  
    image_dimensions_recon_.push_back(((static_cast<unsigned int>(std::ceil(e_space.matrixSize.y*reconstruction_os_factor_y.value()))+warp_size-1)/warp_size)*warp_size);
      
    image_dimensions_recon_os_ = uint64d2
      (((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
       ((static_cast<unsigned int>(std::ceil(image_dimensions_recon_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);

	deblur_ = do_deblurring.value();

	return GADGET_OK;
  }



    int gpuSpiralDeblurGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {

		IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();		
		
		// Allocate various counters if they are NULL
		if( !image_counter_.get() ){
			image_counter_ = boost::shared_array<long>(new long[1]);
			image_counter_[0] = 0;
		}

		hoNDArray<std::complex<float>> host_data = recon_bit_->rbit_[0].data_.data_;
        hoNDArray<float> trajectory = *(recon_bit_->rbit_[0].data_.trajectory_);
        hoNDArray<float> B0_trajectory;
        if(recon_bit_->rbit_.size() > 1){
            B0_trajectory = *(recon_bit_->rbit_[1].data_.trajectory_);
        }
		//host_data = host_data2;
		ISMRMRD::AcquisitionHeader& curr_header = recon_bit_->rbit_[0].data_.headers_(0,0,0,0,0);
				std::cout << "next" << std::endl;
		if(!prepared_ && host_data.get_size(0) > 0){
		    size_t R0 = host_data.get_size(0);
		    size_t E1 = host_data.get_size(1);
		    size_t E2 = host_data.get_size(2);
			size_t CHA = host_data.get_size(3);
			//size_t N = host_data.get_size(4);*/
			//size_t S = host_data.get_size(5);
			//size_t SLC = host_data.get_size(6);			
	
			//Setup image arrays
			std::vector<size_t> image_dims;
			image_dims.push_back(image_dimensions_recon_[0]);
			image_dims.push_back(image_dimensions_recon_[1]);
			image_dims.push_back(CHA);
			image.create(&image_dims);
			host_image.create(&image_dims);	
			host_image.fill(0.0f);	

			// Setup output image array
			image_dims.pop_back();
			reg_image.create(&image_dims);
			host_image.fill(0.0f);
			B0_map.create(&image_dims);
			B0_map.fill(0.0f);
			output_image.create(&image_dims);

			host_traj.create(R0*E1);
			host_weights.create(R0*E1);

			if (curr_header.trajectory_dimensions != 3) {
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
			else{
				for (int i = 0; i < (R0*E1); i++) {
					host_traj[i]   = floatd2(trajectory[i*3],trajectory[i*3+1]);
					host_weights[i] = trajectory[i*3+2];
				}		
			}		

			//upload to gpu
			gpu_traj = host_traj;
			gpu_weights = host_weights;
			//pre-process
			nfft_plan_ = NFFT<cuNDArray,_real,2>::make_plan( from_std_vector<size_t,2>(image_dimensions_recon_), image_dimensions_recon_os_, kernel_width_ );
			nfft_plan_->preprocess(&gpu_traj, NFFT_prep_mode::ALL);
			prepared_ = true;
		}

		if(!prepared_B0_ && recon_bit_->rbit_.size() > 1){

			ISMRMRD::AcquisitionHeader& B0_header = recon_bit_->rbit_[1].data_.headers_(0,0,0,0,0);
			size_t R0 = recon_bit_->rbit_[1].data_.data_.get_size(0);
			size_t E1 = recon_bit_->rbit_[1].data_.data_.get_size(1);	

			B0_traj.create(R0*E1);
			B0_weights.create(R0*E1);

			float krmaxB0_ = 2.*(B0_header.user_float[4])/10000.; //Hack #mcr
			if (B0_header.trajectory_dimensions != 3) {
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
				std::cout << "fov " << B0_header.user_float[5] << std::endl;
				std::cout << "smax " << 3*B0_header.user_float[3]/10. << std::endl;
				std::cout << "gmax " << B0_header.user_float[1]/10. << std::endl;
				std::cout << "krmax " << 2*((B0_header.user_float[4])/10000.) << std::endl;
				std::cout << "sampling_time " << sample_time << std::endl;
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
					B0_traj[i]   = floatd2(B0_trajectory[i*3],B0_trajectory[i*3+1]);
					B0_weights[i] = B0_trajectory[i*3+2];
				}		
			}

			//upload to gpu
			gpu_traj = B0_traj;
			gpu_weights_B0 = B0_weights;
			//write_nd_array<floatd2>( &gpu_traj, "maptraj.real" );

			//pre-process
			nfft_plan_B0_ = NFFT<cuNDArray,_real,2>::make_plan( from_std_vector<size_t,2>(image_dimensions_recon_), image_dimensions_recon_os_, kernel_width_ );
			nfft_plan_B0_->preprocess(&gpu_traj, NFFT_prep_mode::NC2C);
			prepared_B0_= true;
		}

		if(recon_bit_->rbit_.size() > 1){
			size_t R0 = recon_bit_->rbit_[1].data_.data_.get_size(0);
			//#ifdef USE_OMP
			//#pragma omp parallel for 
			//#endif
			for(int i =0; i < recon_bit_->rbit_[1].data_.data_.get_number_of_elements(); i++){
				recon_bit_->rbit_[1].data_.data_[i] *= exp(-.5*pow((i%R0)/(R0/2),3.) );
				recon_bit_->rbit_[2].data_.data_[i] *= exp(-.5*pow((i%R0)/(R0/2),3.) );
			}
			cuNDArray<complext<float>> gpu_B0_data((hoNDArray<float_complext>*)&recon_bit_->rbit_[1].data_.data_);
			nfft_plan_B0_->compute( &gpu_B0_data, image, &gpu_weights_B0, NFFT_comp_mode::BACKWARDS_NC2C );
			csm_ = estimate_b1_map<float,2>( &image );		
			csm_mult_MH<float,2>(&image, &reg_image, &csm_);
			hoNDArray<complext<float>> B0_temp_0 = *reg_image.to_host();
			gpu_B0_data = *((hoNDArray<float_complext>*)&recon_bit_->rbit_[2].data_.data_);
			nfft_plan_B0_->compute( &gpu_B0_data, image, &gpu_weights_B0, NFFT_comp_mode::BACKWARDS_NC2C );	
			csm_mult_MH<float,2>(&image, &reg_image, &csm_);
			hoNDArray<complext<float>> B0_temp_1 = *reg_image.to_host();
			float B0_max;
			for (int i = 0; i < B0_temp_0.get_number_of_elements(); i++) {
			  if(abs(B0_temp_1[i]) > 5){ //assumes SNR units and filters out noise > 5 std
			    B0_map[i] = _real(arg(B0_temp_0[i]*conj(B0_temp_1[i]))/( 2*M_PI*.001 ));//delTE = 1 ms
			  }
		
			}
		}
		
		size_t R0 = host_data.get_size(0);
		size_t E1 = host_data.get_size(1);
		size_t CHA = host_data.get_size(3);
		
		if(R0 > 0){
			

		//Deblur using Multi-frequency Interpolation
		int fmax = (1/(2*.001))*1.2; 
		int L = std::ceil(2.5*fmax*R0*sample_time);
		if(L%2 == 0){ L++; }
		std::complex<float> om (0.0,2*M_PI);
		std::cout << "L = " << L << std::endl;

		//Compute MFI coefficients if not they do not already exist
		if( MFI_C.get_number_of_elements() == 0 ){
			std::cout << "doing something" << std::endl;
			//Setup some arma matrices
			arma::cx_fmat demod( R0 , L );
			MFI_C = hoNDArray<_complext>( fmax*2+1 , L );
			//arma::cx_fvec b( R0 );
			//arma::cx_fvec x( L );
			arma::cx_fmat b( R0, fmax*2+1);
			arma::cx_fmat x( L, fmax*2+1 );
			//Setup and solve least-squares probelm
			int j = 0;
			int i = 0;
			float f = 0.0;
			#ifdef USE_OMP
			#pragma omp parallel for private(f,i,j) shared(fmax, R0, L, om, demod) 
			#endif
			for(j = 0; j<L; j++){
				f = -fmax+j*fmax*2./(L-1);
				for(i = 0; i < R0; i++) {
					demod(i,j) = exp(om*(float)i*(float)sample_time*f);
				}
			}
			j = 0;
			i = 0;
			f = 0.0;
			#ifdef USE_OMP
			#pragma omp parallel for private(f,i,j) shared(fmax, R0, L, om, demod, b, x) 
			#endif
			for(j = 0; j<fmax*2+1; j++){
				f = -fmax+j;
				for(i = 0; i < R0; i++) {
					b(i,j) = exp(om*(float)i*(float)sample_time*f);
				}
				x.col(j) = arma::solve(demod, b.col(j));
				memcpy(MFI_C.get_data_ptr()+j*L, x.colptr(j), L*sizeof(std::complex<float>));		
			}
		}

		cuNDArray<complext<float>> gpu_data((hoNDArray<float_complext>*)&host_data);
		nfft_plan_->compute( &gpu_data, image, &gpu_weights, NFFT_comp_mode::BACKWARDS_NC2C );
		csm_ = estimate_b1_map<float,2>( &image );		
		csm_mult_MH<float,2>(&image, &reg_image, &csm_);
		host_image = *(reg_image.to_host());
		nfft_plan_->fft(image, NFFT_fft_mode::FORWARDS);
		auto gridded_data_0 = *(image.to_host());
		nfft_plan_->fft(image, NFFT_fft_mode::BACKWARDS);


		if( phase_mask.get_number_of_elements() == 0 ) {
			phase_mask.create(gridded_data_0.get_dimensions());
			phase_mask.fill(0.0f);
			float f_step = fmax/((L-1)/2.);
			std::complex<float> omega = std::complex<float>(sample_time,0.0);
			std::cout << omega << std::endl;
			//#ifdef USE_OMP
			//#pragma omp parallel for
			//#endif
			for(int r = 0; r < R0*E1*CHA; r++) {
				host_data[r] = omega*float(r%R0);
			}
			std::cout << sample_time*R0*2.*M_PI*f_step << std::endl;
			gpu_data = *((hoNDArray<float_complext>*)&host_data);
            hoNDArray<float> ones;
            ones.create(R0*E1);
            ones.fill(1.0);
            cuNDArray<float> cu_ones(&ones);
			nfft_plan_->compute( &gpu_data, image, &cu_ones, NFFT_comp_mode::BACKWARDS_NC2C );	
			//csm_mult_MH<float,2>(&image, &reg_image, &deref_csm);
			nfft_plan_->fft(image, NFFT_fft_mode::FORWARDS);
			//sum(&image,2);
			auto gridded_data_1 = *(image.to_host());
            auto offset = gridded_data_1[gridded_data_1.get_size(0)/2];
			for (int i = 0; i < gridded_data_1.get_number_of_elements(); i++) {
				phase_mask[i] = exp(_complext(0.0,2.*M_PI*f_step)*(abs(gridded_data_1[i]-offset)));///_complext(std::sqrt(gridded_data_1.get_size(0)),0.0));
			}
            std::cout << offset << std::endl;
			for(int x = 0; x < gridded_data_1.get_size(0); x++){
				float kx;
                if(x<gridded_data_1.get_size(0)/2){
                    kx = (x*2./gridded_data_1.get_size(0));
                }else{
                    kx = (-2+x*2.0/gridded_data_1.get_size(0));
                }
				for(int y = 0; y < gridded_data_1.get_size(1); y++){
				    float ky;
                    if(y<gridded_data_1.get_size(1)/2){
                        ky = (y*2./gridded_data_1.get_size(1));
                    }else{
                        ky = (-2.+y*2.0/gridded_data_1.get_size(1));
                    }
					phase_mask[x+y*gridded_data_1.get_size(0)] *= .5+std::atan(100.*(1.-std::sqrt(kx*kx+ky*ky)/.9))/M_PI;
				}
			}
            write_nd_array<_complext>( &gridded_data_1, "time_grid.cplx" );
		}		

		for (int j = 0; j < (L-1)/2; j++){
			#ifdef USE_OMP
			#pragma omp parallel for
			#endif
			for (int i = 0; i < gridded_data_0.get_number_of_elements(); i++) {
				gridded_data_0[i] *= conj(phase_mask[i]);
			}
		}
		
		output_image.fill(0.0f);	
		hoNDArray<_complext> temp_image(host_image.get_dimensions());
		int i;
		int j;
		for(j = 0; j<L; j++){
			std::cout << j << std::endl;
			//Update output image
			int mfc_index;
			if(j != 0){
				#ifdef USE_OMP
				#pragma omp parallel for
				#endif
				for (i = 0; i < gridded_data_0.get_number_of_elements(); i++) {
					gridded_data_0[i] *= phase_mask[i];
				}
			}
			image = gridded_data_0;
			nfft_plan_->fft(image, NFFT_fft_mode::BACKWARDS);
			csm_mult_MH<float,2>(&image, &reg_image, &csm_);
			temp_image = *(reg_image.to_host());
			#ifdef USE_OMP
			#pragma omp parallel for private(i,mfc_index)
			#endif
			for (i = 0; i < temp_image.get_number_of_elements(); i++) {
				mfc_index = int(B0_map[i]+fmax)*L+j;
				output_image[i] += (MFI_C[mfc_index]*temp_image[i]);
			}
		}

		std::cout << "finshes" << std::endl;
		// Prepare an image header for this frame
		GadgetContainerMessage<ISMRMRD::ImageHeader> *header = new GadgetContainerMessage<ISMRMRD::ImageHeader>();

		  {
		// Initialize header to all zeroes (there are a few fields we do not set yet)
		ISMRMRD::ImageHeader tmp;
		*(header->getObjectPtr()) = tmp;
		  }

		header->getObjectPtr()->version = curr_header.version;

		header->getObjectPtr()->matrix_size[0] = image_dimensions_recon_[0];
		header->getObjectPtr()->matrix_size[1] = image_dimensions_recon_[1];
		header->getObjectPtr()->matrix_size[2] = 1;

		header->getObjectPtr()->field_of_view[0] = fov_vec_[0];
		header->getObjectPtr()->field_of_view[1] = fov_vec_[1];
		header->getObjectPtr()->field_of_view[2] = fov_vec_[2];
		header->getObjectPtr()->channels = 1;//base_head->active_channels;

		header->getObjectPtr()->slice = curr_header.idx.slice;
		header->getObjectPtr()->set = curr_header.idx.set;

		header->getObjectPtr()->acquisition_time_stamp = curr_header.acquisition_time_stamp;
		memcpy(header->getObjectPtr()->physiology_time_stamp, curr_header.physiology_time_stamp, sizeof(uint32_t)*ISMRMRD::ISMRMRD_PHYS_STAMPS);

		memcpy(header->getObjectPtr()->position, curr_header.position, sizeof(float)*3);
		memcpy(header->getObjectPtr()->read_dir, curr_header.read_dir, sizeof(float)*3);
		memcpy(header->getObjectPtr()->phase_dir, curr_header.phase_dir, sizeof(float)*3);
		memcpy(header->getObjectPtr()->slice_dir, curr_header.slice_dir, sizeof(float)*3);
		memcpy(header->getObjectPtr()->patient_table_position, curr_header.patient_table_position, sizeof(float)*3);

		header->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
		header->getObjectPtr()->image_index = image_counter_[0]++; 
		header->getObjectPtr()->image_series_index = 1;

		GadgetContainerMessage< hoNDArray< std::complex<float> > >* cm2 = new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
		cm2->getObjectPtr()->create(host_image.get_dimensions());
		if(!deblur_){
			//memcpy(cm2->getObjectPtr()->get_data_ptr(), host_image.get_data_ptr(), host_image.get_number_of_elements()*sizeof(std::complex<float>));
			memcpy(cm2->getObjectPtr()->get_data_ptr(), host_image.get_data_ptr(), host_image.get_number_of_elements()*sizeof(std::complex<float>));
		}
		if(deblur_){
			memcpy(cm2->getObjectPtr()->get_data_ptr(), output_image.get_data_ptr(), output_image.get_number_of_elements()*sizeof(std::complex<float>));
		}
		header->cont(cm2);

		if (this->next()->putq(header) < 0) {
		  GDEBUG("Failed to put job on queue.\n");
		  header->release();
		  return GADGET_FAIL;
		}

		GadgetContainerMessage<ISMRMRD::ImageHeader> *header2 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		{
		ISMRMRD::ImageHeader tmp;
		*(header2->getObjectPtr()) = tmp;
		  }

		header2->getObjectPtr()->version = curr_header.version;

		header2->getObjectPtr()->matrix_size[0] = image_dimensions_recon_[0];
		header2->getObjectPtr()->matrix_size[1] = image_dimensions_recon_[1];
		header2->getObjectPtr()->matrix_size[2] = 1;

		header2->getObjectPtr()->field_of_view[0] = fov_vec_[0];
		header2->getObjectPtr()->field_of_view[1] = fov_vec_[1];
		header2->getObjectPtr()->field_of_view[2] = fov_vec_[2];
		header2->getObjectPtr()->channels = 1;//base_head->active_channels;

		header2->getObjectPtr()->slice = curr_header.idx.slice;
		header2->getObjectPtr()->set = curr_header.idx.set;

		header2->getObjectPtr()->acquisition_time_stamp = curr_header.acquisition_time_stamp;
		memcpy(header2->getObjectPtr()->physiology_time_stamp, curr_header.physiology_time_stamp, sizeof(uint32_t)*ISMRMRD::ISMRMRD_PHYS_STAMPS);

		memcpy(header2->getObjectPtr()->position, curr_header.position, sizeof(float)*3);
		memcpy(header2->getObjectPtr()->read_dir, curr_header.read_dir, sizeof(float)*3);
		memcpy(header2->getObjectPtr()->phase_dir, curr_header.phase_dir, sizeof(float)*3);
		memcpy(header2->getObjectPtr()->slice_dir, curr_header.slice_dir, sizeof(float)*3);
		memcpy(header2->getObjectPtr()->patient_table_position, curr_header.patient_table_position, sizeof(float)*3);
		header2->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
		//header->getObjectPtr()->image_index = image_counter_[0]++; 
		header2->getObjectPtr()->image_series_index = 10;
		
		GadgetContainerMessage<hoNDArray< std::complex<float> > >* cm3 = new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
		cm3->getObjectPtr()->create(host_image.get_dimensions());
		memcpy(cm3->getObjectPtr()->get_data_ptr(), host_image.get_data_ptr(), host_image.get_number_of_elements()*sizeof(std::complex<float>));
		header2->cont(cm3);

		if (this->next()->putq(header2) < 0) {
		  GDEBUG("Failed to put job on queue.\n");
		  header2->release();
		  return GADGET_FAIL;
		}

		}
		GadgetContainerMessage<ISMRMRD::ImageHeader> *header3 = new GadgetContainerMessage<ISMRMRD::ImageHeader>();
		{
		ISMRMRD::ImageHeader tmp;
		*(header3->getObjectPtr()) = tmp;
		  }

		header3->getObjectPtr()->version = curr_header.version;

		header3->getObjectPtr()->matrix_size[0] = image_dimensions_recon_[0];
		header3->getObjectPtr()->matrix_size[1] = image_dimensions_recon_[1];
		header3->getObjectPtr()->matrix_size[2] = 1;

		header3->getObjectPtr()->field_of_view[0] = fov_vec_[0];
		header3->getObjectPtr()->field_of_view[1] = fov_vec_[1];
		header3->getObjectPtr()->field_of_view[2] = fov_vec_[2];
		header3->getObjectPtr()->channels = 1;//base_head->active_channels;

		header3->getObjectPtr()->slice = curr_header.idx.slice;
		header3->getObjectPtr()->set = curr_header.idx.set;

		header3->getObjectPtr()->acquisition_time_stamp = curr_header.acquisition_time_stamp;
		memcpy(header3->getObjectPtr()->physiology_time_stamp, curr_header.physiology_time_stamp, sizeof(uint32_t)*ISMRMRD::ISMRMRD_PHYS_STAMPS);

		memcpy(header3->getObjectPtr()->position, curr_header.position, sizeof(float)*3);
		memcpy(header3->getObjectPtr()->read_dir, curr_header.read_dir, sizeof(float)*3);
		memcpy(header3->getObjectPtr()->phase_dir, curr_header.phase_dir, sizeof(float)*3);
		memcpy(header3->getObjectPtr()->slice_dir, curr_header.slice_dir, sizeof(float)*3);
		memcpy(header3->getObjectPtr()->patient_table_position, curr_header.patient_table_position, sizeof(float)*3);
		header3->getObjectPtr()->data_type = ISMRMRD::ISMRMRD_CXFLOAT;
		//header->getObjectPtr()->image_index = image_counter_[0]++; 
		header3->getObjectPtr()->image_series_index = 20;
		
		GadgetContainerMessage<hoNDArray< std::complex<float> > >* cm4 = new GadgetContainerMessage<hoNDArray< std::complex<float> > >();
		hoNDArray< std::complex<float> >B0_image;
		B0_image.create(B0_map.get_dimensions());
		for(int i = 0; i<B0_map.get_number_of_elements(); i++){
			B0_image[i] = std::complex<float>(B0_map[i],0.0);
		}
		cm4->getObjectPtr()->create(B0_image.get_dimensions());
		memcpy(cm4->getObjectPtr()->get_data_ptr(), B0_image.get_data_ptr(), B0_image.get_number_of_elements()*sizeof(std::complex<float>));
		header3->cont(cm4);

		if (this->next()->putq(header3) < 0) {
		  GDEBUG("Failed to put job on queue.\n");
		  header3->release();
		  return GADGET_FAIL;
		}

		m1->release();
		return GADGET_OK;
	}

  GADGET_FACTORY_DECLARE(gpuSpiralDeblurGadget)
}






























