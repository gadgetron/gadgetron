/****************************************************************************************************************************
* Description: Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
* Auhtor: Mahamadou Diakite, PhD.
* Institution: National Institutes of Health (NIH)
* Lang: C++
* Date: 10/19/2016
* Version: 0.0.1
****************************************************************************************************************************/

#include "bartgadget.h"
#include "gadgetron/hoNDFFT.h"
#include "gadgetron/hoNDArray_math.h"
#include <sstream>
#include <utility>
#include <numeric>
#include <cstdlib>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>


namespace Gadgetron {

	BartGadget::BartGadget() :
		image_counter_(0)
	{}

	BartGadget::~BartGadget()
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
		unsigned long count = 0;
		std::vector< std::complex<float> > Data;
		Data.reserve(buffer.size() / 2);

		for (size_t count = 0; count < buffer.size(); count += 2)
			Data.push_back(std::complex<float>(buffer[count], buffer[count + 1]));

		return (std::pair< std::vector<size_t>, std::vector< std::complex<float> > >(DIMS, Data));
	}

	std::string & BartGadget::getOutputFilename(const std::string & bartCommandLine)
	{
		static std::vector<std::string> outputFile;
		boost::char_separator<char> sep(" ");
		boost::tokenizer<boost::char_separator<char> > tokens(bartCommandLine, sep);
		for (auto itr = tokens.begin(); itr != tokens.end(); ++itr)
			outputFile.push_back(*itr);
		return (outputFile.back());
	}

	void BartGadget::cleanup(std::vector<std::string>& createdFiles)
	{
		if (createdFiles.empty())
			return;
		for (auto e : createdFiles) {
			std::ostringstream cmd;
			cmd << "rm -f " << e << ".hdr" << " " << e << ".cfl";
			system(cmd.str().c_str());
		}
	}

	int BartGadget::process(GadgetContainerMessage<IsmrmrdReconData>* m1)
	{
		// Check the status of bart commands script
		std::string CommandScript = BartCommandScript_path.value() + "/" + BartCommandScript_name.value();
		if (!boost::filesystem::exists(CommandScript))
		{
			GERROR("Can't find bart commands script: %s!\n", CommandScript.c_str());
			return GADGET_FAIL;
		}

		// Let's keep track of all files yet to be created
		std::vector<std::string> createdFiles;

		std::vector<uint16_t> DIMS_ref, DIMS;

		/*** WRITE REFERENCE AND RAW DATA TO FILES ***/

		for (std::vector<IsmrmrdReconBit>::iterator it = m1->getObjectPtr()->rbit_.begin(); it != m1->getObjectPtr()->rbit_.end(); ++it)
		{
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
								Temp_ref.push_back(e.real());
								Temp_ref.push_back(e.imag());
							}
						}
					}
				}

				if (!write_BART_Files("meas_gadgetron_ref", DIMS_ref, Temp_ref))
					createdFiles.push_back("meas_gadgetron_ref");
			}

			GadgetContainerMessage<IsmrmrdImageArray>* cm1 = new GadgetContainerMessage<IsmrmrdImageArray>();
			IsmrmrdImageArray & imarray = *cm1->getObjectPtr();

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
							Temp.push_back(e.real());
							Temp.push_back(e.imag());
						}
					}
				}
			}

			if (!write_BART_Files("meas_gadgetron", DIMS, Temp))
				createdFiles.push_back("meas_gadgetron");
		}

		/* Before calling Bart let's do some bookkeeping */
		std::ostringstream cmd1, cmd2, cmd3, cmd4, cmd5;
		if (DIMS_ref != DIMS)
		{
			cmd1 << "bart reshape 15 " << DIMS_ref[2] << " " << DIMS_ref[0] << " " << DIMS_ref[1] << " " << DIMS_ref[3] << " meas_gadgetron_ref meas_gadgetron_ref_reshape";
			cmd2 << "bart resize -c 1 " << DIMS[0] << " 2 " << DIMS[1] << " meas_gadgetron_ref_reshape reference_data";

			// Pass commands to Bart
			if (!system(cmd1.str().c_str()))
				createdFiles.push_back("meas_gadgetron_ref_reshape");

			if (!system(cmd2.str().c_str()))
				createdFiles.push_back("reference_data");
		}

		if (DIMS[4] != 1)
			cmd3 << "bart reshape 1023 " << DIMS[2] << " " << DIMS[0] << " " << DIMS[1] << " " << DIMS[3] << " 1 1 1 1 1 " << DIMS[4] << " meas_gadgetron input_data";
		else
			cmd3 << "bart reshape 15 " << DIMS[2] << " " << DIMS[0] << " " << DIMS[1] << " " << DIMS[3] << " meas_gadgetron input_data";

		if (!system(cmd3.str().c_str()));
			createdFiles.push_back("input_data");

		/*** CALL BART COMMAND LINE from the scripting file***/
		std::string Commands_Line;
		ifstream inputFile(CommandScript);
		if (inputFile.is_open())
		{
			while (getline(inputFile, Commands_Line))
			{
				if (Commands_Line.empty())
					continue;
				GDEBUG("%s\n", Commands_Line.c_str());
				if (!system(Commands_Line.c_str()))
					createdFiles.push_back(getOutputFilename(Commands_Line));
			}
			inputFile.close();
		}
		else
		{
			GERROR("Unable to open %s\n", CommandScript.c_str());
			cleanup(createdFiles);
			return GADGET_FAIL;
		}

		std::string outputFile = getOutputFilename(Commands_Line);
		// Reformat the data back to gadgetron format
		auto header = read_BART_hdr(outputFile.c_str());
		cmd4 << "bart reshape 1023 " << header[1] << " " << header[2] << " " << header[0] << " " << header[3] << " " << header[9] * header[4]
			<< " 1 1 1 1 1 " << outputFile << " " << outputFile + std::string("_reshape");

		const auto cmd_s = cmd4.str();
		GDEBUG("%s\n", cmd_s.c_str());
		if (!system(cmd_s.c_str()))
			createdFiles.push_back(outputFile + std::string("_reshape"));

		std::string outputfile_2 = getOutputFilename(cmd_s);
		/**** READ FROM BART FILES ***/
		std::pair< std::vector<size_t>, std::vector<std::complex<float> > > BART_DATA = read_BART_files(outputfile_2.c_str());

		cleanup(createdFiles);

		GadgetContainerMessage<IsmrmrdImageArray>* ims = new GadgetContainerMessage<IsmrmrdImageArray>();
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
		for (std::vector<IsmrmrdReconBit>::iterator it = m1->getObjectPtr()->rbit_.begin(); it != m1->getObjectPtr()->rbit_.end(); ++it)
		{
			IsmrmrdDataBuffered & dbuff = it->data_;
			// ImageHeaders will be [N, S, LOC]
			std::vector<size_t> header_dims(3);
			header_dims[0] = data_dims_Final[4];
			header_dims[1] = data_dims_Final[5];
			header_dims[2] = data_dims_Final[6];
			imarray.headers_.create(&header_dims);

			for (uint16_t loc = 0; loc < data_dims_Final[6]; ++loc) {
				for (uint16_t s = 0; s < data_dims_Final[5]; ++s) {
					for (uint16_t n = 0; n < data_dims_Final[4]; ++n) {

						ISMRMRD::AcquisitionHeader & acqhdr = dbuff.headers_(dbuff.sampling_.sampling_limits_[1].center_,
							dbuff.sampling_.sampling_limits_[2].center_,
							n, s, loc);
						imarray.headers_(n, s, loc).matrix_size[0] = data_dims_Final[0];
						imarray.headers_(n, s, loc).matrix_size[1] = data_dims_Final[1];
						imarray.headers_(n, s, loc).matrix_size[2] = data_dims_Final[2];
						imarray.headers_(n, s, loc).field_of_view[0] = dbuff.sampling_.recon_FOV_[0];
						imarray.headers_(n, s, loc).field_of_view[1] = dbuff.sampling_.recon_FOV_[1];
						imarray.headers_(n, s, loc).field_of_view[2] = dbuff.sampling_.recon_FOV_[2];
						imarray.headers_(n, s, loc).channels = 1;
						imarray.headers_(n, s, loc).average = acqhdr.idx.average;
						imarray.headers_(n, s, loc).slice = acqhdr.idx.slice;
						imarray.headers_(n, s, loc).contrast = acqhdr.idx.contrast;
						imarray.headers_(n, s, loc).phase = acqhdr.idx.phase;
						imarray.headers_(n, s, loc).repetition = acqhdr.idx.repetition;
						imarray.headers_(n, s, loc).set = acqhdr.idx.set;
						imarray.headers_(n, s, loc).acquisition_time_stamp = acqhdr.acquisition_time_stamp;
						imarray.headers_(n, s, loc).position[0] = acqhdr.position[0];
						imarray.headers_(n, s, loc).position[1] = acqhdr.position[1];
						imarray.headers_(n, s, loc).position[2] = acqhdr.position[2];
						imarray.headers_(n, s, loc).read_dir[0] = acqhdr.read_dir[0];
						imarray.headers_(n, s, loc).read_dir[1] = acqhdr.read_dir[1];
						imarray.headers_(n, s, loc).read_dir[2] = acqhdr.read_dir[2];
						imarray.headers_(n, s, loc).phase_dir[0] = acqhdr.phase_dir[0];
						imarray.headers_(n, s, loc).phase_dir[1] = acqhdr.phase_dir[1];
						imarray.headers_(n, s, loc).phase_dir[2] = acqhdr.phase_dir[2];
						imarray.headers_(n, s, loc).slice_dir[0] = acqhdr.slice_dir[0];
						imarray.headers_(n, s, loc).slice_dir[1] = acqhdr.slice_dir[1];
						imarray.headers_(n, s, loc).slice_dir[2] = acqhdr.slice_dir[2];
						imarray.headers_(n, s, loc).patient_table_position[0] = acqhdr.patient_table_position[0];
						imarray.headers_(n, s, loc).patient_table_position[1] = acqhdr.patient_table_position[1];
						imarray.headers_(n, s, loc).patient_table_position[2] = acqhdr.patient_table_position[2];
						imarray.headers_(n, s, loc).data_type = ISMRMRD::ISMRMRD_CXFLOAT;
						imarray.headers_(n, s, loc).image_index = ++image_counter_;

					}
				}
			}

			if (this->next()->putq(ims) < 0) {
				m1->release();
				return GADGET_FAIL;
			}
		}
		m1->release();
		return GADGET_OK;
	}

	GADGET_FACTORY_DECLARE(BartGadget)
}
