/****************************************************************************************************************************
* Description: Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
* Author: Mahamadou Diakite, PhD.
* Institution: National Institutes of Health (NIH)
* Lang: C++
* Date: 10/15/2017
* Version: 1.0.0
****************************************************************************************************************************/

#include "bartgadget.h"
#include <gadgetron_paths.h>
#include <boost/tokenizer.hpp>
#include <sstream>
#include <utility>
#include <numeric>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <random>
#include <functional>
#include <mutex>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>



namespace Gadgetron {

	BartGadget::BartGadget() :
		BaseClass(),
		dp{}
	{}

	std::vector<size_t> BartGadget::read_BART_hdr(const char *filename)
	{
		std::string filename_s = std::string(filename) + std::string(".hdr");
		std::ifstream infile(filename_s);
		if (!infile.is_open())
			GERROR("Failed to open file: %s\n", filename_s.c_str());

		std::vector<size_t> DIMS;
		if (infile.is_open())
		{
			std::vector<std::string> tokens;
			std::string line;

			while (std::getline(infile, line, '\n'))
			{
				tokens.push_back(line);
			}

			// Parse the dimensions
			const std::string s = tokens[1];
			std::stringstream ss(s);
			std::string items;
			while (getline(ss, items, ' ')) {
				DIMS.push_back(std::stoi(items, nullptr, 10));
			}
			infile.close();
		}

		return(DIMS);
	}


	std::pair< std::vector<size_t>, std::vector<std::complex<float> > >
		BartGadget::read_BART_files(const char *filename)
	{
		// Load the header file
		auto DIMS = read_BART_hdr(filename);

		// Load the cfl file
		std::string filename_s = std::string(filename) + std::string(".cfl");
		std::ifstream infile(filename_s, std::ifstream::binary);
		if (!infile.is_open())
			GERROR("Failed to open file: %s\n", filename_s.c_str());

		infile.seekg(0, infile.end);
		size_t size = static_cast<long>(infile.tellg());
		infile.seekg(0);
		std::vector<float> buffer(size);
		infile.read(reinterpret_cast<char*>(&buffer[0]), size);
		infile.close();
		// Reformat the data
		std::vector< std::complex<float> > Data;
		Data.reserve(buffer.size() / 2);

		for (size_t count = 0, buffer_size = buffer.size(); count != buffer_size; count += 2)
			Data.push_back(std::complex<float>(buffer[count], buffer[count + 1]));

		return (std::pair< std::vector<size_t>, std::vector< std::complex<float> > >(DIMS, Data));
	}

	std::string & BartGadget::getOutputFilename(const std::string & bartCommandLine)
	{
		static std::vector<std::string> outputFile;
		boost::char_separator<char> sep(" ");
		boost::tokenizer<boost::char_separator<char> > tokens(bartCommandLine, sep);
		for (auto itr = tokens.begin(), tokens_end = tokens.end(); itr != tokens_end; ++itr)
			outputFile.push_back(*itr);
		return (outputFile.back());
	}

	inline void BartGadget::cleanup(std::string &createdFiles)
	{
		boost::filesystem::remove_all(createdFiles);
	}

	inline void BartGadget::ltrim(std::string &str)
	{
		str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](int s) {return !std::isspace(s); }));
	}

	inline void BartGadget::rtrim(std::string &str)
	{
		str.erase(std::find_if(str.rbegin(), str.rend(), [](int s) {return !std::isspace(s);}).base(), str.end());
	}

	inline void BartGadget::trim(std::string &str)
	{
		ltrim(str);
		rtrim(str);
	}
	
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
		try
		{
			deserialize(mb->rd_ptr(), h);
		}
		catch (...)
		{
			GDEBUG("BartGadget::process_config: Failed to parse incoming ISMRMRD Header");
		}


		for (size_t ii = 0, h_encoding_size = h.encoding.size(); ii < h_encoding_size; ++ii)
		{
			ISMRMRD::EncodingSpace recon_space = h.encoding[ii].reconSpace;

			GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Encoding matrix size: " << recon_space.matrixSize.x << " " << recon_space.matrixSize.y << " " << recon_space.matrixSize.z);
			GDEBUG_CONDITION_STREAM(isVerboseON.value(), "BartGadget::process_config: Encoding field_of_view : " << recon_space.fieldOfView_mm.x << " " << recon_space.fieldOfView_mm.y << " " << recon_space.fieldOfView_mm.z);
			dp.recon_matrix_x = recon_space.matrixSize.x;
			dp.recon_matrix_y = recon_space.matrixSize.y;
			dp.recon_matrix_z = recon_space.matrixSize.z;
			dp.FOV_x = static_cast<uint16_t>(recon_space.fieldOfView_mm.x);
			dp.FOV_y = static_cast<uint16_t>(recon_space.fieldOfView_mm.y);
			dp.FOV_z = static_cast<uint16_t>(recon_space.fieldOfView_mm.z);

			if (!h.encoding[ii].parallelImaging)
			{
				GDEBUG_STREAM("BartGadget::process_config: Parallel Imaging not enable...");
			}
			else
			{
				ISMRMRD::ParallelImaging p_imaging = *h.encoding[0].parallelImaging;
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

				std::string calib = *p_imaging.calibrationMode;

				auto separate = (calib.compare("separate") == 0);
				auto embedded = (calib.compare("embedded") == 0);
				auto external = (calib.compare("external") == 0);
				auto interleaved = (calib.compare("interleaved") == 0);
				auto other = (calib.compare("other") == 0);

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

		// Check status of bart commands script
		std::string CommandScript = AbsoluteBartCommandScript_path.value() + "/" + BartCommandScript_name.value();
		if (!boost::filesystem::exists(CommandScript))
		{
			GERROR("Can't find bart commands script: %s!\n", CommandScript.c_str());
			return GADGET_FAIL;
		}

		// set the permission for the script
#ifdef _WIN32
		try
		{
			boost::filesystem::permissions(CommandScript, boost::filesystem::all_all);
		}
		catch (...)
		{
			GERROR("Error changing the permission of the command script.\n");
			return GADGET_FAIL;
		}
#else
		// in case an older version of boost is used in non-win system
		// the system call is used
		int res = chmod(CommandScript.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IWGRP | S_IXGRP | S_IROTH | S_IWOTH | S_IXOTH);
		if (res != 0)
		{
			GERROR("Error changing the permission of the command script.\n");
			return GADGET_FAIL;
		}
#endif // _WIN32

		// Check status of the folder containing the generated files (*.hdr & *.cfl)
		std::string generatedFilesFolder;
		static std::string outputFolderPath;
        if (BartWorkingDirectory_path.value().empty()){
            GERROR("Error: No BART working directory provided!");
            return GADGET_FAIL;
        }

		time_t rawtime;
		char buff[80];
		time(&rawtime);
		strftime(buff, sizeof(buff), "%H_%M_%S__", localtime(&rawtime));
		std::mt19937::result_type seed = static_cast<unsigned long>(time(0));
		auto dice_rand = std::bind(std::uniform_int_distribution<int>(1, 10000), std::mt19937(seed));
		std::string time_id(buff + std::to_string(dice_rand()));
                // Get the current process ID
 		std::string threadId = boost::lexical_cast<std::string>(boost::this_thread::get_id());
                unsigned long threadNumber = 0;
                sscanf(threadId.c_str(), "%lx", &threadNumber);
                outputFolderPath = BartWorkingDirectory_path.value() + "bart_" + time_id + "_" + std::to_string(threadNumber) + "/";
                generatedFilesFolder = std::string(outputFolderPath);
                
		generatedFilesFolder.pop_back();
	
		boost::filesystem::path dir(generatedFilesFolder);
		if (!boost::filesystem::exists(dir) || !boost::filesystem::is_directory(dir)){
			if (boost::filesystem::create_directories(dir)){
				GDEBUG("Folder to store *.hdr & *.cfl files is %s\n", generatedFilesFolder.c_str());
                        }
			else {
				GERROR("Folder to store *.hdr & *.cfl files doesn't exist...\n");
				return GADGET_FAIL;
			}
                }

		generatedFilesFolder += "/";

		/*USE WITH CAUTION*/
		if (boost::filesystem::exists(generatedFilesFolder) && isBartFolderBeingCachedToVM.value() && !isBartFileBeingStored.value())
		{
			std::ostringstream cmd;
			cmd << "mount -t tmpfs -o size" << AllocateMemorySizeInMegabytes.value() << "M, mode=0755 tmpfs " << generatedFilesFolder;
			if (system(cmd.str().c_str())) {
				cleanup(outputFolderPath);
				return GADGET_FAIL;
			}
		}


		std::vector<uint16_t> DIMS_ref, DIMS;

		/*** WRITE REFERENCE AND RAW DATA TO FILES ***/

		size_t encoding = 0;
		for (std::vector<IsmrmrdReconBit>::iterator it = m1->getObjectPtr()->rbit_.begin(), rbit_end =  m1->getObjectPtr()->rbit_.end(); it != rbit_end; ++it)
		{
			std::stringstream os;
			os << "_encoding_" << encoding;

			// Grab a reference to the buffer containing the reference data
			auto  & dbuff_ref = it->ref_;
			hoNDArray< std::complex<float> >& ref = (*dbuff_ref).data_;
			// Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
			uint16_t E0_ref = static_cast<uint16_t>(ref.get_size(0));
			uint16_t E1_ref = static_cast<uint16_t>(ref.get_size(1));
			uint16_t E2_ref = static_cast<uint16_t>(ref.get_size(2));
			uint16_t CHA_ref = static_cast<uint16_t>(ref.get_size(3));
			uint16_t N_ref = static_cast<uint16_t>(ref.get_size(4));
			uint16_t S_ref = static_cast<uint16_t>(ref.get_size(5));
			uint16_t LOC_ref = static_cast<uint16_t>(ref.get_size(6));
			DIMS_ref = { E0_ref, E1_ref, E2_ref, CHA_ref, N_ref, S_ref, LOC_ref };

			// Grab a reference to the buffer containing the image data
			IsmrmrdDataBuffered & dbuff = it->data_;
			// Data 7D, fixed order [E0, E1, E2, CHA, N, S, LOC]
			uint16_t E0 = static_cast<uint16_t>(dbuff.data_.get_size(0));
			uint16_t E1 = static_cast<uint16_t>(dbuff.data_.get_size(1));
			uint16_t E2 = static_cast<uint16_t>(dbuff.data_.get_size(2));
			uint16_t CHA = static_cast<uint16_t>(dbuff.data_.get_size(3));
			uint16_t N = static_cast<uint16_t>(dbuff.data_.get_size(4));
			uint16_t S = static_cast<uint16_t>(dbuff.data_.get_size(5));
			uint16_t LOC = static_cast<uint16_t>(dbuff.data_.get_size(6));
			// Set up data to be written into files
			DIMS = { E0, E1, E2, CHA, N, S, LOC };

			/* The reference data will be pointing to the image data if there is
				no reference scan. Therefore, we won't write the reference data
				into files if it's pointing to the raw data.*/
			if (DIMS_ref != DIMS) {
				std::vector<float> Temp_ref;
				Temp_ref.reserve(2 * E0_ref*E1_ref*E2_ref*CHA_ref*N_ref*S_ref*LOC_ref);

				for (uint16_t loc = 0; loc < LOC_ref; ++loc) {
					for (uint16_t s = 0; s < S_ref; ++s) {
						for (uint16_t n = 0; n < N_ref; ++n) {

							//Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
							//Each chunk will be [E0,E1,E2,CHA] big
							std::vector<size_t> chunk_dims(4);
							chunk_dims[0] = E0_ref;
							chunk_dims[1] = E1_ref;
							chunk_dims[2] = E2_ref;
							chunk_dims[3] = CHA_ref;
							hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &(*dbuff_ref).data_(0, 0, 0, 0, n, s, loc));
							std::vector<size_t> new_chunk_dims(1);
							new_chunk_dims[0] = E0_ref*E1_ref*E2_ref*CHA_ref;
							chunk.reshape(new_chunk_dims);
							// Fill BART container
							for (auto e : chunk) {
								Temp_ref.push_back(std::move(e.real()));
								Temp_ref.push_back(std::move(e.imag()));
							}
						}
					}
				}

				write_BART_Files(std::string(generatedFilesFolder + "meas_gadgetron_ref").c_str(), DIMS_ref, Temp_ref);
			}

			auto cm1 = std::make_unique<GadgetContainerMessage<IsmrmrdImageArray>>();
			std::vector<float> Temp;
			Temp.reserve(2 * E0*E1*E2*CHA*N*S*LOC);

			for (uint16_t loc = 0; loc < LOC; loc++) {
				for (uint16_t s = 0; s < S; s++) {
					for (uint16_t n = 0; n < N; n++) {

						//Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
						//Each chunk will be [E0,E1,E2,CHA] big
						std::vector<size_t> chunk_dims(4);
						chunk_dims[0] = E0;
						chunk_dims[1] = E1;
						chunk_dims[2] = E2;
						chunk_dims[3] = CHA;
						hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &dbuff.data_(0, 0, 0, 0, n, s, loc));

						std::vector<size_t> new_chunk_dims(1);
						new_chunk_dims[0] = E0*E1*E2*CHA;
						chunk.reshape(new_chunk_dims);
						// Fill BART container
						for (auto e : chunk) {
							Temp.push_back(std::move(e.real()));
							Temp.push_back(std::move(e.imag()));
						}
					}
				}
			}

			write_BART_Files(std::string(generatedFilesFolder + "meas_gadgetron").c_str(), DIMS, Temp);

			encoding++;
		}

		/* Before calling Bart let's do some bookkeeping */
		std::ostringstream cmd1, cmd2, cmd3;
		std::replace(generatedFilesFolder.begin(), generatedFilesFolder.end(), '\\', '/');

		if (DIMS_ref != DIMS)
		{
			cmd1 << "bart resize -c 0 " << DIMS[0] << " 1 " << DIMS[1] << " 2 " << DIMS[2] << " meas_gadgetron_ref reference_data";
			// Pass commands to Bart
			if (system(std::string("cd " + generatedFilesFolder + "&&" + cmd1.str()).c_str())) {
				cleanup(outputFolderPath);
				return GADGET_FAIL;
			}
		}

		if (DIMS[4] != 1)
			cmd2 << "bart reshape 1023 " << DIMS[0] << " " << DIMS[1] << " " << DIMS[2] << " " << DIMS[3] << " 1 1 1 " << DIMS[5] << " " << DIMS[6] << " " << DIMS[4] << " meas_gadgetron input_data";
		else	
			cmd2 << "bart scale 1.0 meas_gadgetron input_data";

		if (system(std::string("cd " + generatedFilesFolder + "&&" + cmd2.str()).c_str())) {
			cleanup(outputFolderPath);
			return GADGET_FAIL;
		}

		/*** CALL BART COMMAND LINE from the scripting file***/

		std::string Line, Commands_Line;
		std::fstream inputFile(CommandScript);
		if (inputFile.is_open())
		{
			while (getline(inputFile, Line))
			{
				// crop comment
				Line = Line.substr(0, Line.find_first_of("#"));

				trim(Line);
				if (Line.empty() || Line.compare(0, 4, "bart") != 0)
					continue;
				
				replace_default_parameters(Line);
				GDEBUG("%s\n", Line.c_str());
				
				auto ret = system(std::string("cd " + generatedFilesFolder + "&&" + Line).c_str());
				(void)ret;
				
				Commands_Line = Line;
			}
			inputFile.close();
		}
		else
		{
			GERROR("Unable to open %s\n", CommandScript.c_str());
			cleanup(outputFolderPath);
			return GADGET_FAIL;
		}

		std::string outputFile = getOutputFilename(Commands_Line);

		// Reformat the data back to gadgetron format
		auto header = read_BART_hdr(std::string(generatedFilesFolder + outputFile).c_str());
		cmd3 << "bart reshape 1023 " << header[0] << " " << header[1] << " " << header[2] << " " << header[3] << " " << header[9] * header[4] <<
			" " << header[5] << " " << header[6] << " " << header[7] << " " << header[8] << " 1 " << outputFile << " " << outputFile + std::string("_reshape");
		const auto cmd_s = cmd3.str();
		GDEBUG("%s\n", cmd_s.c_str());
		if (system(std::string("cd " + generatedFilesFolder + "&&" + cmd_s).c_str())) {
			cleanup(outputFolderPath);
			return GADGET_FAIL;
		}

		std::string outputfile_2 = getOutputFilename(cmd_s);
		/**** READ FROM BART FILES ***/
		std::pair< std::vector<size_t>, std::vector<std::complex<float> > > BART_DATA = read_BART_files(std::string(generatedFilesFolder + outputfile_2).c_str());

		if (!isBartFileBeingStored.value())
			cleanup(outputFolderPath);

		auto ims = std::make_unique<GadgetContainerMessage<IsmrmrdImageArray>>();
		IsmrmrdImageArray & imarray = *ims->getObjectPtr();

		// Grab data from BART files
		std::vector<size_t> BART_DATA_dims(1);
		BART_DATA_dims[0] = std::accumulate(BART_DATA.first.begin(), BART_DATA.first.end(), 1, std::multiplies<size_t>());
		hoNDArray<std::complex<float>> DATA = hoNDArray<std::complex<float>>(BART_DATA_dims, &BART_DATA.second[0]);

		// The image array data will be [E0,E1,E2,1,N,S,LOC]
		std::vector<size_t> data_dims(7);
		data_dims[0] = BART_DATA.first[0];
		data_dims[1] = BART_DATA.first[1];
		data_dims[2] = BART_DATA.first[2];
		data_dims[3] = BART_DATA.first[3];
		data_dims[4] = BART_DATA.first[4];
		data_dims[5] = BART_DATA.first[5];
		data_dims[6] = BART_DATA.first[6];

		DATA.reshape(data_dims);

		// Extract the first image from each time frame (depending on the number of maps generated by the user)
		std::vector<size_t> data_dims_Final(7);
		data_dims_Final[0] = BART_DATA.first[0];
		data_dims_Final[1] = BART_DATA.first[1];
		data_dims_Final[2] = BART_DATA.first[2];
		data_dims_Final[3] = BART_DATA.first[3];
		assert(header[4] > 0);
		data_dims_Final[4] = BART_DATA.first[4] / header[4];
		data_dims_Final[5] = BART_DATA.first[5];
		data_dims_Final[6] = BART_DATA.first[6];
		imarray.data_.create(&data_dims_Final);

		std::vector<std::complex<float> > DATA_Final;
		DATA_Final.reserve(std::accumulate(data_dims_Final.begin(), data_dims_Final.end(), 1, std::multiplies<size_t>()));

		for (uint16_t loc = 0; loc < data_dims[6]; ++loc) {
			for (uint16_t s = 0; s < data_dims[5]; ++s) {
				for (uint16_t n = 0; n < data_dims[4]; n += header[4]) {

					//Grab a wrapper around the relevant chunk of data [E0,E1,E2,CHA] for this loc, n, and s
					//Each chunk will be [E0,E1,E2,CHA] big
					std::vector<size_t> chunk_dims(4), Temp_one_1d(1);
					chunk_dims[0] = data_dims_Final[0];
					chunk_dims[1] = data_dims_Final[1];
					chunk_dims[2] = data_dims_Final[2];
					chunk_dims[3] = data_dims_Final[3];
					hoNDArray<std::complex<float> > chunk = hoNDArray<std::complex<float> >(chunk_dims, &DATA(0, 0, 0, 0, n, s, loc));

					Temp_one_1d[0] = chunk_dims[0] * chunk_dims[1] * chunk_dims[2] * chunk_dims[3];
					chunk.reshape(Temp_one_1d);
					DATA_Final.insert(DATA_Final.end(), chunk.begin(), chunk.end());
				}
			}
		}

		std::copy(DATA_Final.begin(), DATA_Final.end(), imarray.data_.begin());

		// Fill image header 
		for (size_t it = 0, rbit_size = m1->getObjectPtr()->rbit_.size(); it != rbit_size; ++it)
		{
			compute_image_header(m1->getObjectPtr()->rbit_[it], imarray, it);
			send_out_image_array(m1->getObjectPtr()->rbit_[it], imarray, it, image_series.value() + (static_cast<int>(it) + 1), GADGETRON_IMAGE_REGULAR);
		}

		m1->release();
		return GADGET_OK;
	}

	GADGET_FACTORY_DECLARE(BartGadget)
}
