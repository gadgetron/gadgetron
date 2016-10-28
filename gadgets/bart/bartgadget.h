/****************************************************************************************************************************
* Description: Gadget using the Berkeley Advanced Reconstruction Toolbox (BART)
* Auhtor: Mahamadou Diakite, PhD.
* Institution: National Institutes of Health (NIH)
* Lang: C++
* Date: 10/19/2016
* Version: 0.0.1
****************************************************************************************************************************/

#ifndef BART_GADGET_H
#define BART_GADGET_H

#include "Gadget.h"
#include "gadgetron_mricore_export.h"
#include "gadgetron/mri_core_data.h"

#include <vector>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <iostream>


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

	class EXPORTGADGETS_bartgadget BartGadget : public Gadget1<IsmrmrdReconData>

	{

	public:
		GADGET_DECLARE(BartGadget);

		BartGadget();
		virtual ~BartGadget();

	protected:
		GADGET_PROPERTY(BartCommandScript_path, std::string, "Path(s) to add to the to the application search path", "/home/myuser/local/share/gadgetron/bart");
		GADGET_PROPERTY(BartCommandScript_name, std::string, "Script file containing bart command(s) to be loaded", "");

		virtual int process(GadgetContainerMessage<IsmrmrdReconData>* m1);
		long long image_counter_;

	private:
		// Write BART files
		template<typename T>
		void write_BART_hdr(const char* filename, std::vector<T> &DIMS);
		template<typename T, typename U>
		int write_BART_Files(const char* filename, std::vector<T> &DIMS, std::vector<U> &DATA);

		// Read BART files
		std::vector<size_t> read_BART_hdr(const char *filename);
		std::pair< std::vector<size_t>, std::vector<std::complex<float> > > read_BART_files(const char *filename);

		// Utility functions
		std::string &getOutputFilename(const std::string &bartCommandLine);
		void cleanup(std::vector<std::string> &createdFiles);
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
	int BartGadget::write_BART_Files(const char* filename, std::vector<T> &DIMS, std::vector<U> &DATA)
	{
		write_BART_hdr(filename, DIMS);
		std::string filename_s = std::string(filename) + std::string(".cfl");
		std::ofstream pFile(filename_s, std::ofstream::out | std::ofstream::binary);
		if (!pFile.is_open())
			GERROR("Failed to write into file: %s\n", filename);
		pFile.write(reinterpret_cast<char*>(&DATA[0]), DATA.size()*sizeof(float));
		pFile.close();
		return GADGET_OK;
	}
}
#endif //BART_GADGET_H
