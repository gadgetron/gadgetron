/******************************************************************************
 * Description: Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
 * Authors: 
 *   Mahamadou Diakite, PhD. [1]
 *   Nguyen Damien, PhD. [2]
 *   Francesco Santini, PhD. [2]
 * Institutions: 
 *   [1] National Institutes of Health (NIH)
 *   [2] University of Basel, Switzerland
 * Lang: C++
 * Date: 08/10/2018
 * Version: 1.5.0
 ******************************************************************************/

#include "bartgadget.h"
#include <boost/tokenizer.hpp>
#include <sstream>
#include <utility>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <memory>
#include <functional>
#include <mutex>

#include <cerrno>
#ifdef _WIN32
    #include <direct.h>
    #define getcwd _getcwd
    #define chdir _chdir
#else
    #include <unistd.h>
#endif

#include "bart_helpers.h"

enum bart_debug_levels { BART_DP_ERROR, BART_DP_WARN, BART_DP_INFO, BART_DP_DEBUG1, BART_DP_DEBUG2, BART_DP_DEBUG3, BART_DP_DEBUG4, BART_DP_TRACE, BART_DP_ALL };

// =============================================================================

std::vector<size_t> Gadgetron::read_BART_hdr(fs::path &filename)
{
     return read_BART_hdr<size_t>(filename);
}

std::pair<std::vector<size_t>, std::vector<std::complex<float>>>
Gadgetron::read_BART_files(fs::path &filename)
{
     return read_BART_files<size_t>(filename);
}

// =============================================================================

namespace Gadgetron {

     BartGadget::BartGadget() :
	  dp{}
     {}

     void BartGadget::replace_default_parameters(std::string & str)
     {
	  std::string::size_type pos = 0u;
	  while ((pos = str.find('$', pos)) != std::string::npos)
	  {
	       auto pos_end = str.find(' ', pos);
	       auto pos_diff = pos_end - pos;
	       std::string tmp = str.substr(pos, pos_diff);
	       tmp.erase(0, 1);
	       if (tmp == std::string("recon_matrix_x"))
		    str.replace(pos, pos_diff, std::to_string(dp.recon_matrix_x));
	       else if (tmp == std::string("recon_matrix_y"))
		    str.replace(pos, pos_diff, std::to_string(dp.recon_matrix_y));
	       else if (tmp == std::string("recon_matrix_z"))
		    str.replace(pos, pos_diff, std::to_string(dp.recon_matrix_z));
	       else if (tmp == std::string("FOV_x"))
		    str.replace(pos, pos_diff, std::to_string(dp.FOV_x));
	       else if (tmp == std::string("FOV_y"))
		    str.replace(pos, pos_diff, std::to_string(dp.FOV_y));
	       else if (tmp == std::string("FOV_z"))
		    str.replace(pos, pos_diff, std::to_string(dp.FOV_z));
	       else if (tmp == std::string("acc_factor_PE1"))
		    str.replace(pos, pos_diff, std::to_string(dp.acc_factor_PE1));
	       else if (tmp == std::string("acc_factor_PE2"))
		    str.replace(pos, pos_diff, std::to_string(dp.acc_factor_PE2));
	       else if (tmp == std::string("reference_lines_PE1"))
		    str.replace(pos, pos_diff, std::to_string(dp.reference_lines_PE1));
	       else if (tmp == std::string("reference_lines_PE2"))
		    str.replace(pos, pos_diff, std::to_string(dp.reference_lines_PE2));
	       else if (tmp == std::string("input_data"))
		    str.replace(pos, pos_diff, dp.input_data);
	       else if (tmp == std::string("reference_data"))
		    str.replace(pos, pos_diff, dp.reference_data);
	       else if (tmp == std::string("traj_data"))
		    str.replace(pos, pos_diff, dp.traj_data);
	       else {
		    GERROR( "Unknown default parameter, please see the complete list of available parameters...");
	       }
	       pos = pos_end;
	  }
     }
     
	
     int BartGadget::process_config(ACE_Message_Block * mb)
     {
	  GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

	  /** Let's get some information about the incoming data **/
	  ISMRMRD::IsmrmrdHeader h;
	  try {
	       deserialize(mb->rd_ptr(), h);
	  }
	  catch (...) {
	       GDEBUG("BartGadget::process_config: Failed to parse incoming ISMRMRD Header");
	  }

	  // ===================================================================
	  /* Data provided to the user might or might not be in-memory
	   *
	   * - if -DMEMONLY_CFL then we're always in-memory without extension
	   * - otherwise we add the *.mem extension if
	   *     + if the memory behaviour is BART_ALL_IN_MEM
	   *     + or if the memory behaviour is BART_MIX_DISK_MEM and the user requests it
	   */
	  const auto append_mem_ext_in(!memonly_cfl
				       && (memory_behaviour_ == BART_ALL_IN_MEM
					   || (memory_behaviour_ == BART_MIX_DISK_MEM
					       && BartStoreGadgetronInputInMemory.value())));
	  
	  dp.reference_data  = std::string("reference_data") + (append_mem_ext_in ? ".mem" : "");
	  dp.input_data = std::string("input_data")     + (append_mem_ext_in ? ".mem" : "");
	  dp.traj_data = std::string("traj_data")      + (append_mem_ext_in ? ".mem" : "");
	  
	  // ===================================================================
	  // Adjust BART's debug level
	  
	  GDEBUG_STREAM("BartGadget::process_config: setting BART debug level to "
			<< BartPrintfDebugLevel.value());
	  
	  auto& bart_debug_level(debug_level);
	  if (BartPrintfDebugLevel.value() == "DP_ERROR") {
	       bart_debug_level = BART_DP_ERROR;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_WARN") {
	       bart_debug_level = BART_DP_WARN;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_INFO") {
	       bart_debug_level = BART_DP_INFO;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_DEBUG1") {
	       bart_debug_level = BART_DP_DEBUG1;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_DEBUG2") {
	       bart_debug_level = BART_DP_DEBUG2;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_DEBUG3") {
	       bart_debug_level = BART_DP_DEBUG3;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_DEBUG4") {
	       bart_debug_level = BART_DP_DEBUG4;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_TRACE") {
	       bart_debug_level = BART_DP_TRACE;
	  }
	  else if (BartPrintfDebugLevel.value() == "DP_ALL") {
	       bart_debug_level = BART_DP_ALL;
	  }
	  else {
	       GWARN("Unknown debug level! defaulting to BART_DP_INFO");
	       bart_debug_level = BART_DP_INFO;
	  }

	  // ===================================================================
	  // Check status of bart commands script

	  if (AbsoluteBartCommandScript_path.value().empty())
	      AbsoluteBartCommandScript_path.value(context.paths.gadgetron_home.string() + "/share/gadgetron/bart");

	  command_script_ = AbsoluteBartCommandScript_path.value();
	  command_script_ /= BartCommandScript_name.value();

	  if (!fs::exists(command_script_)) {
	       GERROR("Can't find bart commands script: %s!\n", command_script_.c_str());
	       return GADGET_FAIL;
	  }

	  // ===================================================================

	  memory_behaviour_ = BART_MIX_DISK_MEM;
	  if (BartFileBehaviour.value() == "BART_ALL_IN_MEM") {
	       memory_behaviour_ = BART_ALL_IN_MEM;
	  }
	  else if (BartFileBehaviour.value() == "BART_MIX_DISK_MEM") {
	       memory_behaviour_ = BART_MIX_DISK_MEM;
	  }
	  else {
	       GERROR_STREAM("Invalid value specified for BartFileBehaviour: " << BartFileBehaviour.value());
	       return GADGET_FAIL;
	  }

	  // ===================================================================

	  for (const auto& enc: h.encoding) {
	       auto recon_space = enc.reconSpace;

	       GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Encoding matrix size: " << recon_space.matrixSize.x << " " << recon_space.matrixSize.y << " " << recon_space.matrixSize.z);
	       GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Encoding field_of_view : " << recon_space.fieldOfView_mm.x << " " << recon_space.fieldOfView_mm.y << " " << recon_space.fieldOfView_mm.z);
	       dp.recon_matrix_x = recon_space.matrixSize.x;
	       dp.recon_matrix_y = recon_space.matrixSize.y;
	       dp.recon_matrix_z = recon_space.matrixSize.z;
	       dp.FOV_x = static_cast<uint16_t>(recon_space.fieldOfView_mm.x);
	       dp.FOV_y = static_cast<uint16_t>(recon_space.fieldOfView_mm.y);
	       dp.FOV_z = static_cast<uint16_t>(recon_space.fieldOfView_mm.z);

	       if (!enc.parallelImaging)
	       {
		    GDEBUG_STREAM("BartGadget::process_config: Parallel Imaging not enable...");
	       }
	       else
	       {
		    auto p_imaging = *h.encoding.front().parallelImaging;
		    GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: acceleration Factor along PE1 is " << p_imaging.accelerationFactor.kspace_encoding_step_1);
		    GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: acceleration Factor along PE2 is " << p_imaging.accelerationFactor.kspace_encoding_step_2);
		    dp.acc_factor_PE1 = p_imaging.accelerationFactor.kspace_encoding_step_1;
		    dp.acc_factor_PE2 = p_imaging.accelerationFactor.kspace_encoding_step_2;

		    if (p_imaging.accelerationFactor.kspace_encoding_step_2 > 1) {
			 GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Limits of the size of the calibration region (PE1) " << h.userParameters->userParameterLong[0].name << " is " << h.userParameters->userParameterLong[0].value);
			 GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Limits of the size of the calibration region (PE2) " << h.userParameters->userParameterLong[1].name << " is " << h.userParameters->userParameterLong[1].value);
			 dp.reference_lines_PE1 = h.userParameters->userParameterLong[0].value;
			 dp.reference_lines_PE2 = h.userParameters->userParameterLong[1].value;
		    }
		    else if (p_imaging.accelerationFactor.kspace_encoding_step_1 > 1) {
			 GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Limits of the size of the calibration region (PE1) " << h.userParameters->userParameterLong[0].name << " is " << h.userParameters->userParameterLong[0].value);
			 dp.reference_lines_PE1 = h.userParameters->userParameterLong[0].value;
		    }

		    auto calib = *p_imaging.calibrationMode;
            auto separate = (calib == std::string("separate"));
            auto embedded = (calib == std::string("embedded"));
            auto external = (calib == std::string("external"));
            auto interleaved = (calib == std::string("interleaved"));
            auto other = (calib == std::string("other"));

		    if (p_imaging.accelerationFactor.kspace_encoding_step_1 > 1 || p_imaging.accelerationFactor.kspace_encoding_step_2 > 1)
		    {
			 if (interleaved) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode INTERLEAVE ");
			 }
			 else if (embedded) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode EMBEDDED");
			 }
			 else if (separate) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode SEPERATE");
			 }
			 else if (external) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode EXTERNAL");
			 }
			 else if (other) {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Calibration mode OTHER");
			 }
			 else {
			      GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Something went terribly wrong, this should never happen!");
			      return GADGET_FAIL;
			 }
		    }

	       }

	  }
	  return GADGET_OK;
     }

     int BartGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
     {
	  static std::mutex mtx;
	  std::lock_guard<std::mutex> guard(mtx);

	  constexpr auto memonly_cfl = 
#ifdef MEMONLY_CFL
	       true
#else
	       false
#endif /* MEMONLY_CFL */
	       ;
	  
	  GINFO_STREAM("Process start");

  
	  auto generated_files_folder(internal::generate_unique_folder(BartWorkingDirectory_path.value()));

	  bool create_folder(!memonly_cfl
			     || memory_behaviour_ != BART_ALL_IN_MEM);
	  if (create_folder) {
	       if (!(fs::exists(generated_files_folder)
		     && fs::is_directory(generated_files_folder))) {
		    if (fs::create_directories(generated_files_folder)){
			 GDEBUG_STREAM("Folder to store *.hdr & *.cfl files is " << generated_files_folder);
		    }
		    else {
			 GERROR("Folder to store *.hdr & *.cfl files doesn't exist...\n");
			 return GADGET_FAIL;
		    }
	       }
	  }
	  	  
	  internal::ScopeGuard cleanup_guard(generated_files_folder);
	  if (!create_folder) {
	       cleanup_guard.dismiss();	       
	  }

	  // Gadgetron data is always in-memory
	  const auto ref_filename_src  = std::string("meas_gadgetron_ref")   + (!memonly_cfl ? ".mem" : "");
	  const auto data_filename_src = std::string("meas_gadgetron_input") + (!memonly_cfl ? ".mem" : "");
	  const auto traj_filename_src = std::string("meas_gadgetron_traj")  + (!memonly_cfl ? ".mem" : "");

	  
	  /*** PROCESS EACH DATASET ***/
	  	  
	  auto it(0UL);
	  for(auto& recon_bit: m1->getObjectPtr()->rbit_) {

	       bool has_traj(recon_bit.data_.trajectory_);
	       
	       // Grab a reference to the buffer containing the reference data
	       auto& input_ref = (*recon_bit.ref_).data_;
	       // Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
	       std::vector<long> DIMS_ref{static_cast<long>(input_ref.get_size(0)),
					  static_cast<long>(input_ref.get_size(1)),
					  static_cast<long>(input_ref.get_size(2)),
					  static_cast<long>(input_ref.get_size(3)),
					  static_cast<long>(input_ref.get_size(4)),
					  static_cast<long>(input_ref.get_size(5)),
					  static_cast<long>(input_ref.get_size(6))};

	       // Grab a reference to the buffer containing the image data
	       auto& input = recon_bit.data_.data_;
	       // Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
	       std::vector<long> DIMS{static_cast<long>(input.get_size(0)),
				      static_cast<long>(input.get_size(1)),
				      static_cast<long>(input.get_size(2)),
				      static_cast<long>(input.get_size(3)),
				      static_cast<long>(input.get_size(4)),
				      static_cast<long>(input.get_size(5)),
				      static_cast<long>(input.get_size(6))};

	       // Grab a reference to the buffer containing the image trajectory data (if present)
	       std::unique_ptr<hoNDArray<std::complex<float>>> traj(nullptr);
	       if (has_traj) {
		    auto& traj_real = *recon_bit.data_.trajectory_;
		    traj = std::make_unique<hoNDArray<std::complex<float>>>(traj_real.get_dimensions());
		    // Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
		    std::vector<long> DIMS_traj{static_cast<long>(traj_real.get_size(0)),
						static_cast<long>(traj_real.get_size(1)),
						static_cast<long>(traj_real.get_size(2)),
						static_cast<long>(traj_real.get_size(3)),
						static_cast<long>(traj_real.get_size(4)),
						static_cast<long>(traj_real.get_size(5)),
						static_cast<long>(traj_real.get_size(6))};
		    std::transform(traj_real.begin(), traj_real.end(), 
				   traj->begin(), 
				   [] (float r) { return std::complex<float>(r, 0.); });

		    register_mem_cfl_non_managed(traj_filename_src.c_str(), DIMS_traj.size(), &DIMS_traj[0], &(*traj)[0]);
	       }

	       /* The reference data will be pointing to the image data if there is
		  no reference scan. Therefore, we won't write the reference data
		  into files if it's pointing to the raw data.*/
	       if (DIMS_ref != DIMS)
	       {
		    register_mem_cfl_non_managed(ref_filename_src.c_str(), DIMS_ref.size(), &DIMS_ref[0], &input_ref[0]);
	       }

	       register_mem_cfl_non_managed(data_filename_src.c_str(), DIMS.size(), &DIMS[0], &input[0]);

	       
	       if (DIMS_ref != DIMS)
	       {
		    std::ostringstream cmd;
		    cmd << "bart resize -c 0 " << DIMS[0] << " 1 " << DIMS[1] << " 2 " << DIMS[2]
			<< " " << ref_filename_src << " " << dp.reference_data;

		    GDEBUG_STREAM("Gadgetron data is loaded to " << ref_filename_src);
		    GDEBUG_STREAM("BART filename for reference data is " << dp.reference_data);
		    
		    if (!call_BART(cmd.str()))
		    {
			 return GADGET_FAIL;
		    }
	       }

	       std::ostringstream cmd2;
	       if (DIMS[4] != 1)
		    cmd2 << "bart reshape 1023 " << DIMS[0] << " " << DIMS[1] << " " << DIMS[2] << " " << DIMS[3] << " 1 1 1 " << DIMS[5] << " " << DIMS[6] << " " << DIMS[4];
	       else	
		    cmd2 << "bart copy";

	       cmd2 << " " << data_filename_src << " " << dp.input_data;
	       
	       GDEBUG_STREAM("Gadgetron data is loaded to " << data_filename_src);
	       GDEBUG_STREAM("BART filename for data is " << dp.input_data);
	       
	       if (!call_BART(cmd2.str()))
	       {
		    return GADGET_FAIL;
	       }

	       if (has_traj) {
		    std::ostringstream cmd3;
		    if (DIMS[4] != 1)
			 cmd3 << "bart reshape 1023 " << DIMS[0] << " " << DIMS[1] << " " << DIMS[2] << " " << DIMS[3] << " 1 1 1 " << DIMS[5] << " " << DIMS[6] << " " << DIMS[4];
		    else	
			 cmd3 << "bart copy";
		    
		    cmd3 << " " << traj_filename_src << " " << dp.traj_data;
	       
		    GDEBUG_STREAM("Gadgetron trajectory is loaded to " << traj_filename_src);
		    GDEBUG_STREAM("BART filename for trajectory is " << dp.traj_data);
		    
		    if (!call_BART(cmd3.str()))
		    {
			 return GADGET_FAIL;
		    }
	       }

	       /*** CALL BART COMMAND LINE from the scripting file ***/
	       GDEBUG("Starting processing user script\n");

	       std::string Commands_Line;
           std::ifstream inputFile(command_script_.string());
           if (inputFile.is_open())
	       {
		    std::string Line;
		    while (getline(inputFile, Line))
		    {
			 // crop comment
             Line = Line.substr(0, Line.find_first_of('#'));

			 internal::trim(Line);
			 if (Line.empty() || Line.compare(0, 4, "bart") != 0)
			      continue;
				
			 replace_default_parameters(Line);
			 if (!call_BART(Line))
			 {
			      return GADGET_FAIL;
			 }

			 Commands_Line = Line;
		    }
	       }
	       else
	       {
		    GERROR("Unable to open %s\n", command_script_.c_str());
		    return GADGET_FAIL;
	       }

	       // ==============================================================
	       
	       fs::path outputFile(internal::get_output_filename(Commands_Line));
	       GDEBUG_STREAM("Detected last output file: " << outputFile);

	       // Reshaped data is always in-memory
	       fs::path outputFileReshape(outputFile);

	       // Remove the extension from the user's output if needed
	       if (outputFileReshape.extension() == ".mem") {
		    outputFileReshape.replace_extension();
	       }
	       outputFileReshape += std::string("_reshape") + (!memonly_cfl ? ".mem" : "");
	       
	       // Reformat the data back to Gadgetron format
	       std::vector<long> header(16);
	       if (memonly_cfl || outputFile.extension() == ".mem") {
		    load_mem_cfl(outputFile.c_str(), header.size(), header.data()); 
	       }
	       else {
            auto tmp(generated_files_folder / outputFile);
            header = read_BART_hdr<long>(tmp);
	       }
	       
	       std::ostringstream cmd4;
	       cmd4 << "bart reshape 1023 " << header[0] << " " << header[1] << " " << header[2] << " " << header[3] << " " << header[9] * header[4] <<
		    " " << header[5] << " " << header[6] << " " << header[7] << " " << header[8] << " 1 " << outputFile.c_str() << " " << outputFileReshape.c_str();

	       if (!call_BART(cmd4.str()))
	       {
		    return GADGET_FAIL;
	       }
	  
	       /**** READ FROM BART FILES ***/
	       std::vector<long> DIMS_OUT(16);
	       auto data = reinterpret_cast<std::complex<float>*>(load_mem_cfl(outputFileReshape.c_str(), DIMS_OUT.size(), DIMS_OUT.data()));
	       
           if (data == nullptr) {
		    GERROR("Failed to retrieve data from in-memory CFL file!");
		    return GADGET_FAIL;
	       }


	       if (isBartFileBeingStored.value()) {
		    cleanup_guard.dismiss();
	       }

	       IsmrmrdImageArray imarray;

	       // Grab data from BART files
	       std::vector<size_t> BART_DATA_dims{
		    static_cast<size_t>(std::accumulate(DIMS_OUT.begin(), DIMS_OUT.end(), 1, std::multiplies<size_t>()))};
	       hoNDArray<std::complex<float>> DATA(BART_DATA_dims, data);

	       // The image array data will be [E0,E1,E2,1,N,S,LOC]
	       std::vector<size_t> data_dims(DIMS_OUT.begin(), DIMS_OUT.begin()+7);
	       DATA.reshape(data_dims);

	       // Extract the first image from each time frame (depending on the number of maps generated by the user)
	       std::vector<size_t> data_dims_Final{static_cast<size_t>(DIMS_OUT[0]),
						   static_cast<size_t>(DIMS_OUT[1]),
						   static_cast<size_t>(DIMS_OUT[2]),
						   static_cast<size_t>(DIMS_OUT[3]),
						   static_cast<size_t>(DIMS_OUT[4] / header[4]),
						   static_cast<size_t>(DIMS_OUT[5]),
						   static_cast<size_t>(DIMS_OUT[6])};
	       assert(header[4] > 0);
	       imarray.data_.create(data_dims_Final);

	       std::vector<std::complex<float> > DATA_Final;
	       DATA_Final.reserve(std::accumulate(data_dims_Final.begin(), data_dims_Final.end(), 1, std::multiplies<size_t>()));

	       //Each chunk will be [E0,E1,E2,CHA] big
	       std::vector<size_t> chunk_dims{data_dims_Final[0], data_dims_Final[1], data_dims_Final[2], data_dims_Final[3]};
	       const std::vector<size_t> Temp_one_1d(1, chunk_dims[0] * chunk_dims[1] * chunk_dims[2] * chunk_dims[3]);
	  
	       for (uint16_t loc = 0; loc < data_dims[6]; ++loc) {
		    for (uint16_t s = 0; s < data_dims[5]; ++s) {
			 for (uint16_t n = 0; n < data_dims[4]; n += header[4]) {
			      //Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
			      auto chunk = hoNDArray<std::complex<float> >(chunk_dims, &DATA(0, 0, 0, 0, n, s, loc));
			      chunk.reshape(Temp_one_1d);
			      DATA_Final.insert(DATA_Final.end(), chunk.begin(), chunk.end());
			 }
		    }
	       }

	       std::copy(DATA_Final.begin(), DATA_Final.end(), imarray.data_.begin());

	       compute_image_header(recon_bit, imarray, it);
	       send_out_image_array(imarray, it, image_series.value() + (static_cast<int>(it) + 1), GADGETRON_IMAGE_REGULAR);
	       ++it;
	  }

	  m1->release();
	  return GADGET_OK;
     }

     bool call_BART(const std::string &cmdline)
     {
	  GINFO_STREAM("Executing BART command: " << cmdline);
	  enum { MAX_ARGS = 256 };

	  int argc(0);
	  char* argv[MAX_ARGS];

	  auto cmdline_s = std::make_unique<char[]>(cmdline.size()+1);
	  strcpy(cmdline_s.get(), cmdline.c_str());

	  char *p2 = strtok(cmdline_s.get(), " ");
	  while (p2 && argc < MAX_ARGS-1)
	  {
	       argv[argc++] = p2;
           p2 = strtok(nullptr, " ");
	  }
	  argv[argc] = nullptr;
	  
	  // boost::char_separator<char> sep(" ");
	  // boost::tokenizer<boost::char_separator<char>> tokens(cmdline.begin(),
	  // 						       cmdline.end(),
	  // 						       sep);
	  // std::vector<std::unique_ptr<char[]>> tmp;
	  // auto k(0UL);
	  // for (auto tok: tokens) {
	  //      tmp.push_back(std::make_unique<char[]>(tok.size()+1));
	  //      strcpy(tmp.back().get(), tok.c_str());
	  //      argv[k++] = tmp.back().get();
	  // }
	  // argc = tmp.size();

	  char out_str[512] = {'\0'};
	  auto ret(bart_command(512, out_str, argc, argv));
	  if (ret == 0) {
	       if (strlen(out_str) > 0) {
		    GINFO(out_str);
	       }
	       return true;
	  }
     // else {
	       GERROR_STREAM("BART command failed with return code: " << ret);
	       return false;
     // }
     }

     GADGET_FACTORY_DECLARE(BartGadget)
}
