#include "GriddingReconGadget.h"
#include "mri_core_grappa.h"
#include "cuNFFTOperator.h"
#include "cuNFFT.h"
#include "vector_td_utilities.h"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include "b1_map.h"
#include "cuCgSolver.h"
#include "cuNDArray_math.h"
#include "hoNDArray_math.h"
#include "hoNDArray_utils.h"
#include "cuNDArray_fileio.h"
#include "cudaDeviceManager.h"
#include <numeric>
#include <random>

namespace Gadgetron {

	GriddingReconGadget::GriddingReconGadget() : BaseClass()
	{
	}

	GriddingReconGadget::~GriddingReconGadget()
	{
	}

	int GriddingReconGadget::process_config(ACE_Message_Block* mb)
	{
		GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

		// -------------------------------------------------

		ISMRMRD::IsmrmrdHeader h;
		try
		{
			deserialize(mb->rd_ptr(), h);
		}
		catch (...)
		{
			GDEBUG("Error parsing ISMRMRD Header");
		}

		auto matrixsize = h.encoding.front().encodedSpace.matrixSize;


		kernel_width_ = kernel_width.value();
		oversampling_factor_ = gridding_oversampling_factor.value();
		
		image_dims_.push_back(matrixsize.x);
		image_dims_.push_back(matrixsize.y);
		
		//Figure out what the oversampled matrix size should be taking the warp size into consideration. 
		unsigned int warp_size = cudaDeviceManager::Instance()->warp_size();
		image_dims_os_ = uint64d2
			(((static_cast<size_t>(std::ceil(image_dims_[0]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size,
			 ((static_cast<size_t>(std::ceil(image_dims_[1]*oversampling_factor_))+warp_size-1)/warp_size)*warp_size);
		
		// In case the warp_size constraint kicked in
		oversampling_factor_ = float(image_dims_os_[0])/float(image_dims_[0]);
		
		return GADGET_OK;
	}

	int GriddingReconGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
	{
		if (perform_timing.value()) { gt_timer_local_.start("GriddingReconGadget::process"); }

		process_called_times_++;

		IsmrmrdReconData* recon_bit_ = m1->getObjectPtr();
		if (recon_bit_->rbit_.size() > num_encoding_spaces_)
		{
			GWARN_STREAM("Incoming recon_bit has more encoding spaces than the protocol : " << recon_bit_->rbit_.size() << " instead of " << num_encoding_spaces_);
		}

		// for every encoding space
		for (size_t e = 0; e < recon_bit_->rbit_.size(); e++)
		{

			GDEBUG_CONDITION_STREAM(verbose.value(), "Calling " << process_called_times_ << " , encoding space : " << e);
			GDEBUG_CONDITION_STREAM(verbose.value(), "======================================================================");

			IsmrmrdDataBuffered* buffer = &(recon_bit_->rbit_[e].data_);

			size_t RO = buffer->data_.get_size(0);
			size_t E1 = buffer->data_.get_size(1);
			size_t E2 = buffer->data_.get_size(2);
			size_t CHA = buffer->data_.get_size(3);
			size_t N = buffer->data_.get_size(4);
			size_t S = buffer->data_.get_size(5);
			size_t SLC = buffer->data_.get_size(6);

			if (E2 > 1) {
				GERROR("3D data is not supported in GriddingReconGadget\n");
				m1->release();
				return GADGET_FAIL;
			}
			
			if (buffer->trajectory_ == boost::none) {
				GERROR("Trajectories not found. Bailing out.\n");
				m1->release();
				return GADGET_FAIL;
			}


			std::vector<size_t> new_order = {0,1,2,4,5,6,3};

			boost::shared_ptr<cuNDArray<float>> dcw;
			boost::shared_ptr<cuNDArray<floatd2>> traj;

			auto & trajectory = *buffer->trajectory_;
			
			if (buffer->headers_[0].trajectory_dimensions == 3){
				auto traj_dcw = separate_traj_and_dcw(&trajectory);
				dcw = boost::make_shared<cuNDArray<float>>(std::get<1>(traj_dcw).get());
				traj = boost::make_shared<cuNDArray<floatd2>>(std::get<0>(traj_dcw).get());
			} else if (buffer->headers_[0].trajectory_dimensions == 2){
				auto old_traj_dims = *trajectory.get_dimensions();
				std::vector<size_t> traj_dims (old_traj_dims.begin()+1,old_traj_dims.end()); //Remove first element
				hoNDArray<floatd2> tmp_traj(traj_dims,(floatd2*)trajectory.get_data_ptr());
				traj = boost::make_shared<cuNDArray<floatd2>>(tmp_traj);
			} else {
				throw std::runtime_error("Unsupported number of trajectory dimensions");
			}

			auto permuted = permute((hoNDArray<float_complext>*)&buffer->data_,&new_order);
			cuNDArray<float_complext> data(*permuted);
			
			if (dcw){
				float scale_factor = float(prod(image_dims_os_))/asum(dcw.get());
				*dcw *= scale_factor;
			}

			//Gridding
			auto images = reconstruct(&data,traj.get(),dcw.get(),CHA);

			//Calculate coil sensitivity map
			auto csm = estimate_b1_map<float,2>(images.get());

                        //Coil combine
			*images *= *conj(csm.get());
			auto combined = sum(images.get(),images->get_number_of_dimensions()-1);
				
			auto host_img = combined->to_host();


			IsmrmrdImageArray imarray;

			auto elements = imarray.data_.get_number_of_elements();
			imarray.data_ = std::move(*boost::reinterpret_pointer_cast<decltype(imarray.data_)>(host_img));
//			memcpy(imarray.data_.get_data_ptr(), host_img->get_data_ptr(), host_img->get_number_of_bytes());

			this->compute_image_header(recon_bit_->rbit_[e], imarray, e);
			this->send_out_image_array(recon_bit_->rbit_[e], imarray, e, ((int)e + 1), GADGETRON_IMAGE_REGULAR);
			

			//Is this where we measure SNR?
			if (replicas.value() > 0 && snr_frame.value() == process_called_times_) {
						
				hoNDArray<std::complex<float> > rep_array(image_dims_[0], image_dims_[1], replicas.value());
				
				std::mt19937 engine;
				std::normal_distribution<float> distribution;
				for (size_t r = 0; r < replicas.value(); ++r) {

					if (r % 10 == 0) {
						GDEBUG("Running pseudo replics %d of %d\n", r, replicas.value());
					}
					hoNDArray<std::complex<float> > dtmp = buffer->data_;
					auto permuted_rep = permute((hoNDArray<float_complext>*)&dtmp,&new_order);
					auto dataptr = permuted_rep->get_data_ptr();

					for (size_t k =0; k <  permuted_rep->get_number_of_elements(); k++){
						dataptr[k] += std::complex<float>(distribution(engine),distribution(engine));
					}
					
					cuNDArray<float_complext> data_rep(*permuted_rep);
					
					images = reconstruct(&data_rep,traj.get(),dcw.get(),CHA);
					
					//Coil combine
					*images *= *conj(csm.get());
					auto combined = sum(images.get(),images->get_number_of_dimensions()-1);
					
					auto host_img = combined->to_host();
					
					auto elements = imarray.data_.get_number_of_elements();
					size_t offset = image_dims_[0]*image_dims_[1]*r;
					
					memcpy(rep_array.get_data_ptr()+offset, host_img->get_data_ptr(), sizeof(float)*2*elements);
				}
				

				hoNDArray<float> mag(rep_array.get_dimensions());
				hoNDArray<float> mean(image_dims_[0],image_dims_[1]);
				hoNDArray<float> std(image_dims_[0],image_dims_[1]);
				
				Gadgetron::abs(rep_array, mag);
				Gadgetron::sum_over_dimension(mag,mean,2);
				Gadgetron::scal(1.0f/replicas.value(), mean);
				
				mag -= mean;
				mag *= mag;
				
				Gadgetron::sum_over_dimension(mag,std,2);
				Gadgetron::scal(1.0f/(replicas.value()-1), std);
				Gadgetron::sqrt_inplace(&std);
				
				//SNR image
				mean /= std;
				imarray.data_ = *real_to_complex< std::complex<float> >(&mean);
				
				this->compute_image_header(recon_bit_->rbit_[e], imarray, e);
				this->send_out_image_array(recon_bit_->rbit_[e], imarray, e, image_series.value() + 100 * ((int)e + 3), GADGETRON_IMAGE_SNR_MAP);
			}
		}
		
		m1->release();

		if (perform_timing.value()) { gt_timer_local_.stop(); }

		return GADGET_OK;
	}

	boost::shared_ptr<cuNDArray<float_complext> > GriddingReconGadget::reconstruct(
		cuNDArray<float_complext>* data,
		cuNDArray<floatd2>* traj,
		cuNDArray<float>* dcw,
		size_t ncoils ) {

		//We have density compensation and iteration is set to false
		if (!iterate.value() && dcw) { 

			cuNFFT_plan<float,2> plan(from_std_vector<size_t,2>(image_dims_),image_dims_os_,kernel_width_);
			std::vector<size_t> recon_dims = image_dims_;
			recon_dims.push_back(ncoils);
			auto result = new cuNDArray<float_complext>(recon_dims);

			std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
			cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());
			
			plan.preprocess(&flat_traj,cuNFFT_plan<float,2>::NFFT_PREP_NC2C);
			plan.compute(data,result,dcw,cuNFFT_plan<float,2>::NFFT_BACKWARDS_NC2C);

			return boost::shared_ptr<cuNDArray<float_complext>>(result);
			
		} else { //No density compensation, we have to do iterative reconstruction.
			std::vector<size_t> recon_dims = image_dims_;
			recon_dims.push_back(ncoils);

			auto E = boost::make_shared<cuNFFTOperator<float,2>>();

			E->setup(from_std_vector<size_t,2>(image_dims_),image_dims_os_,kernel_width_);
			std::vector<size_t> flat_dims = {traj->get_number_of_elements()};
			cuNDArray<floatd2> flat_traj(flat_dims,traj->get_data_ptr());

			E->set_domain_dimensions(&recon_dims);
			cuCgSolver<float_complext> solver;
			solver.set_max_iterations(iteration_max.value());
			solver.set_encoding_operator(E);
			solver.set_tc_tolerance(iteration_tol.value());
			solver.set_output_mode(cuCgSolver<float_complext>::OUTPUT_SILENT);
			E->set_codomain_dimensions(data->get_dimensions().get());
			E->preprocess(&flat_traj);
			auto res = solver.solve(data);
			return res;
		}
	}


	std::tuple<boost::shared_ptr<hoNDArray<floatd2 > >, boost::shared_ptr<hoNDArray<float >>> GriddingReconGadget::separate_traj_and_dcw(
		hoNDArray<float >* traj_dcw) {
		std::vector<size_t> dims = *traj_dcw->get_dimensions();
		std::vector<size_t> reduced_dims(dims.begin()+1,dims.end()); //Copy vector, but leave out first dim
		auto  dcw = boost::make_shared<hoNDArray<float>>(reduced_dims);

		auto traj = boost::make_shared<hoNDArray<floatd2>>(reduced_dims);

		auto dcw_ptr = dcw->get_data_ptr();
		auto traj_ptr = traj->get_data_ptr();
		auto ptr = traj_dcw->get_data_ptr();
		for (size_t i = 0; i < traj_dcw->get_number_of_elements()/3; i++){
			traj_ptr[i][0] = ptr[i*3];
			traj_ptr[i][1] = ptr[i*3+1];
			dcw_ptr[i] = ptr[i*3+2];
		}

		return std::make_tuple(traj,dcw);
	}

	void GriddingReconGadget::compute_image_header(IsmrmrdReconBit &recon_bit, IsmrmrdImageArray &res, size_t e) {

        size_t RO = res.data_.get_size(0);
        size_t E1 = res.data_.get_size(1);
        size_t E2 = res.data_.get_size(2);
        size_t CHA = res.data_.get_size(3);
        size_t N = res.data_.get_size(4);
        size_t S = res.data_.get_size(5);
        size_t SLC = res.data_.get_size(6);


        res.headers_.create(N, S, SLC);
        res.meta_.resize(N*S*SLC);

        size_t n, s, slc;

        for (slc = 0; slc < SLC; slc++)
        {
            for (s = 0; s < S; s++)
            {
                for (n = 0; n < N; n++)
                {

                    const ISMRMRD::AcquisitionHeader& acq_header = recon_bit.data_.headers_(0,0,n,s,slc);
                    ISMRMRD::ImageHeader& im_header = res.headers_(n, s, slc);
                    ISMRMRD::MetaContainer& meta = res.meta_[n + s*N + slc*N*S];

                    im_header.version = acq_header.version;
                    im_header.data_type = ISMRMRD::ISMRMRD_CXFLOAT;
                    im_header.flags = acq_header.flags;
                    im_header.measurement_uid = acq_header.measurement_uid;

                    im_header.matrix_size[0] = (uint16_t)RO;
                    im_header.matrix_size[1] = (uint16_t)E1;
                    im_header.matrix_size[2] = (uint16_t)E2;

                    std::copy(recon_bit.data_.sampling_.recon_FOV_,std::end(recon_bit.data_.sampling_.recon_FOV_),im_header.field_of_view);

                    im_header.channels = (uint16_t)CHA;

                    std::copy(acq_header.position,std::end(acq_header.position),im_header.position);

					std::copy(acq_header.read_dir,std::end(acq_header.read_dir),im_header.read_dir);

					std::copy(acq_header.phase_dir,std::end(acq_header.phase_dir),im_header.phase_dir);
					std::copy(acq_header.slice_dir,std::end(acq_header.slice_dir),im_header.slice_dir);
					std::copy(acq_header.patient_table_position,std::end(acq_header.patient_table_position),im_header.patient_table_position);


                    im_header.average = acq_header.idx.average;
                    im_header.slice = acq_header.idx.slice;
                    im_header.contrast = acq_header.idx.contrast;
                    im_header.phase = acq_header.idx.phase;
                    im_header.repetition = acq_header.idx.repetition;
                    im_header.set = acq_header.idx.set;

                    im_header.acquisition_time_stamp = acq_header.acquisition_time_stamp;

					std::copy(acq_header.physiology_time_stamp,std::end(acq_header.physiology_time_stamp),im_header.physiology_time_stamp);

                    im_header.image_type = ISMRMRD::ISMRMRD_IMTYPE_COMPLEX;
                    im_header.image_index = (uint16_t)(n + s*N + slc*N*S);
                    im_header.image_series_index = 0;

                    std::copy(acq_header.user_float,std::end(acq_header.user_float),im_header.user_float);
					std::copy(acq_header.user_int,std::end(acq_header.user_int),im_header.user_int);

                    im_header.attribute_string_len = 0;

                    meta.set("encoding", (long)e);

                    meta.set("encoding_FOV"         , recon_bit.data_.sampling_.encoded_FOV_[0]);
                    meta.append("encoding_FOV"      , recon_bit.data_.sampling_.encoded_FOV_[1]);
                    meta.append("encoding_FOV"      , recon_bit.data_.sampling_.encoded_FOV_[2]);

                    meta.set("recon_FOV"            , recon_bit.data_.sampling_.recon_FOV_[0]);
                    meta.append("recon_FOV"         , recon_bit.data_.sampling_.recon_FOV_[1]);
                    meta.append("recon_FOV"         , recon_bit.data_.sampling_.recon_FOV_[2]);

                    meta.set("encoded_matrix"       , (long)recon_bit.data_.sampling_.encoded_matrix_[0]);
                    meta.append("encoded_matrix"    , (long)recon_bit.data_.sampling_.encoded_matrix_[1]);
                    meta.append("encoded_matrix"    , (long)recon_bit.data_.sampling_.encoded_matrix_[2]);

                    meta.set("recon_matrix"         , (long)recon_bit.data_.sampling_.recon_matrix_[0]);
                    meta.append("recon_matrix"      , (long)recon_bit.data_.sampling_.recon_matrix_[1]);
                    meta.append("recon_matrix"      , (long)recon_bit.data_.sampling_.recon_matrix_[2]);

                    meta.set("sampling_limits_RO"   , (long)recon_bit.data_.sampling_.sampling_limits_[0].min_);
                    meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].center_);
                    meta.append("sampling_limits_RO", (long)recon_bit.data_.sampling_.sampling_limits_[0].max_);

                    meta.set("sampling_limits_E1"   , (long)recon_bit.data_.sampling_.sampling_limits_[1].min_);
                    meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].center_);
                    meta.append("sampling_limits_E1", (long)recon_bit.data_.sampling_.sampling_limits_[1].max_);

                    meta.set("sampling_limits_E2"   , (long)recon_bit.data_.sampling_.sampling_limits_[2].min_);
                    meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].center_);
                    meta.append("sampling_limits_E2", (long)recon_bit.data_.sampling_.sampling_limits_[2].max_);

                    meta.set("PatientPosition", (double)res.headers_(n, s, slc).position[0]);
                    meta.append("PatientPosition", (double)res.headers_(n, s, slc).position[1]);
                    meta.append("PatientPosition", (double)res.headers_(n, s, slc).position[2]);

                    meta.set("read_dir", (double)res.headers_(n, s, slc).read_dir[0]);
                    meta.append("read_dir", (double)res.headers_(n, s, slc).read_dir[1]);
                    meta.append("read_dir", (double)res.headers_(n, s, slc).read_dir[2]);

                    meta.set("phase_dir", (double)res.headers_(n, s, slc).phase_dir[0]);
                    meta.append("phase_dir", (double)res.headers_(n, s, slc).phase_dir[1]);
                    meta.append("phase_dir", (double)res.headers_(n, s, slc).phase_dir[2]);

                    meta.set("slice_dir", (double)res.headers_(n, s, slc).slice_dir[0]);
                    meta.append("slice_dir", (double)res.headers_(n, s, slc).slice_dir[1]);
                    meta.append("slice_dir", (double)res.headers_(n, s, slc).slice_dir[2]);

                    meta.set("patient_table_position", (double)res.headers_(n, s, slc).patient_table_position[0]);
                    meta.append("patient_table_position", (double)res.headers_(n, s, slc).patient_table_position[1]);
                    meta.append("patient_table_position", (double)res.headers_(n, s, slc).patient_table_position[2]);

                    meta.set("acquisition_time_stamp", (long)res.headers_(n, s, slc).acquisition_time_stamp);

                    meta.set("physiology_time_stamp", (long)res.headers_(n, s, slc).physiology_time_stamp[0]);
                    meta.append("physiology_time_stamp", (long)res.headers_(n, s, slc).physiology_time_stamp[1]);
                    meta.append("physiology_time_stamp", (long)res.headers_(n, s, slc).physiology_time_stamp[2]);
                }
            }
        }
	}

	GADGET_FACTORY_DECLARE(GriddingReconGadget)
}
