/****************************************************************************************************************************
* Description: Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
* Author: Mahamadou Diakite, PhD.
* Institution: National Institutes of Health (NIH)
* Lang: C++
* Date: 10/15/2017
* Version: 1.0.0
****************************************************************************************************************************/

#ifndef BART_GADGET_H
#define BART_GADGET_H

#include "GenericReconGadget.h"
#include "gadgetron_mricore_export.h"
#include "mri_core_data.h"
#include <gadgetron_paths.h>
#include <boost/filesystem.hpp>
#include <vector>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <string>


#if defined (WIN32)
#ifdef __BUILD_GADGETRON_bartgadget__
#define EXPORTGADGETS_bartgadget __declspec(dllexport)
#else
#define EXPORTGADGETS_bartgadget __declspec(dllimport)
#endif
#else
#define EXPORTGADGETS_bartgadget
#endif


namespace Gadgetron {
	
	// The user is free to add more parameters as the need arises.
	struct Default_parameters 
	{
		uint16_t recon_matrix_x;
		uint16_t recon_matrix_y;
		uint16_t recon_matrix_z;
		uint16_t FOV_x;
		uint16_t FOV_y;
		uint16_t FOV_z;
		uint16_t acc_factor_PE1;
		uint16_t acc_factor_PE2;
		uint16_t reference_lines_PE1;
		uint16_t reference_lines_PE2;
	};

	class EXPORTGADGETS_bartgadget BartGadget final : public GenericReconGadget
	{

	public:
		GADGET_DECLARE(BartGadget)
		
		using BaseClass = GenericReconGadget;

		BartGadget();
		~BartGadget() = default;

	protected:
		GADGET_PROPERTY(isVerboseON, bool, "Display some information about the incoming data", false);
        GADGET_PROPERTY(BartWorkingDirectory_path, std::string, "Absolute path to temporary file location", "/tmp/gadgetron/");
		GADGET_PROPERTY(AbsoluteBartCommandScript_path, std::string, "Absolute path to bart script(s)", get_gadgetron_home() + "/share/gadgetron/bart");
		GADGET_PROPERTY(BartCommandScript_name, std::string, "Script file containing bart command(s) to be loaded", "");
		GADGET_PROPERTY(isBartFileBeingStored, bool, "Store Bart file", false);

		GADGET_PROPERTY(image_series, int, "Set image series", 0);

		/*Caution: this option must be enable only if the user has root privilege and able to allocation virtual memory*/
		GADGET_PROPERTY(isBartFolderBeingCachedToVM, bool, "Mount bart directory to the virtual memory for better performance", false);
		GADGET_PROPERTY(AllocateMemorySizeInMegabytes, int, "Allocate memory to bart directory", 50);

		int process_config(ACE_Message_Block* mb);
		int process(GadgetContainerMessage<IsmrmrdReconData>* m1);		


	private:
		Default_parameters dp;
		
		// Write BART files
		template<typename T>
		void write_BART_hdr(const char* filename, std::vector<T> &DIMS);
		template<typename T, typename U>
		void write_BART_Files(const char* filename, std::vector<T> &DIMS, std::vector<U> &DATA);

		// Read BART files
		std::vector<size_t> read_BART_hdr(const char *filename);
		std::pair< std::vector<size_t>, std::vector<std::complex<float> > > read_BART_files(const char *filename);

		// Utility functions
		std::string &getOutputFilename(const std::string &bartCommandLine);
		void cleanup(std::string &createdFiles);
		void ltrim(std::string &str);
		void rtrim(std::string &str);
		void trim(std::string &str);
		void replace_default_parameters(std::string &str);
	};


	template<typename T>
	void BartGadget::write_BART_hdr(const char* filename, std::vector<T> &DIMS)
	{
		const size_t MAX_DIMS = 16;
		std::string filename_s = std::string(filename) + std::string(".hdr");
		std::vector<size_t> v(MAX_DIMS, 1);
		assert(DIMS.size() < MAX_DIMS);
		std::copy(DIMS.cbegin(), DIMS.cend(), v.begin());
		std::ofstream pFile;
		pFile.open(filename_s, std::ofstream::out);
		if (!pFile.is_open())
			GERROR("Failed to write into file: %s\n", filename);
		pFile << "# Dimensions\n";
		std::copy(v.cbegin(), v.cend(), std::ostream_iterator<size_t>(pFile, " "));
		pFile.close();
	}

	template<typename T, typename U>
	void BartGadget::write_BART_Files(const char* filename, std::vector<T> &DIMS, std::vector<U> &DATA)
	{
		write_BART_hdr(filename, DIMS);
		std::string filename_s = std::string(filename) + std::string(".cfl");
		std::ofstream pFile(filename_s, std::ofstream::out | std::ofstream::binary);
		if (!pFile.is_open())
			GERROR("Failed to write into file: %s\n", filename);
		pFile.write(reinterpret_cast<char*>(&DATA[0]), DATA.size() * sizeof(float));
		pFile.close();
	}
} // namespace Gadgetron
#endif //BART_GADGET_H
