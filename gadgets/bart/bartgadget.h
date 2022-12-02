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

#ifndef BART_GADGET_H
#define BART_GADGET_H

#include "generic_recon_gadgets/GenericReconGadget.h"
#include "gadgetron_mricore_export.h"
#include "mri_core_data.h"

#include <vector>
#include <cassert>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <string>

#include <boost/filesystem.hpp>

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
     namespace fs = boost::filesystem;

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
	  std::string reference_data;
	  std::string input_data;
	  std::string traj_data;		
     };

     class EXPORTGADGETS_bartgadget BartGadget final : public GenericReconGadget
     {
     public:
	  enum bart_memory_behaviour {BART_ALL_IN_MEM, BART_ALL_ON_DISK, BART_MIX_DISK_MEM};
	  
	  GADGET_DECLARE(BartGadget);
		
	  using BaseClass = GenericReconGadget;

	  BartGadget();
	  ~BartGadget() = default;

     protected:
	  GADGET_PROPERTY(isVerboseON, bool, "Display some information about the incoming data", false);
	  // This property has no effect if BartGadget is compiled with -DMEMONLY_CFL or if the memory behaviour is either BART_ALL_IN_MEM or BART_ALL_ON_DISK
	  GADGET_PROPERTY(BartStoreGadgetronInputInMemory, bool, "Whether BartGadget should always store incoming data in-memory (might append *.mem extension)", true);
	  // The property controls how BART stores CFL files. Possible values are: BART_ALL_IN_MEM, BART_MIX_MEM_DISK, BART_ALL_ON_DISK
	  GADGET_PROPERTY(BartFileBehaviour, std::string, "Controls how BART stores files: either all in memory, or mixed disk/memory behaviour", "BART_MIX_DISK_MEM");
	  GADGET_PROPERTY(BartWorkingDirectory_path, std::string, "Absolute path to a temporary file location for generated BART files", "/tmp/gadgetron/");
	  GADGET_PROPERTY(AbsoluteBartCommandScript_path, std::string, "Absolute path to BART script(s)", "");
	  GADGET_PROPERTY(BartCommandScript_name, std::string, "Script file containing BART command(s) to be loaded", "");
	  GADGET_PROPERTY(isBartFileBeingStored, bool, "Keep generated BART files on the disk after the processing ends (has no effect if BART_ALL_IN_MEM is specified as behaviour)", false);
	  GADGET_PROPERTY(image_series, int, "Set image series", 0);
	  GADGET_PROPERTY(BartPrintfDebugLevel, std::string, "Debug level for BART outputs (DP_ERROR, DP_WARN, DP_INFO, DP_DEBUG1,etc.)", "DP_INFO");

	  int process_config(ACE_Message_Block* mb);
	  int process(GadgetContainerMessage<IsmrmrdReconData>* m1);		

     private:

	  static constexpr auto memonly_cfl = 
#ifdef MEMONLY_CFL
	       true
#else
	       false
#endif /* MEMONLY_CFL */
	       ;
	
	  Default_parameters dp;
	  bart_memory_behaviour memory_behaviour_;
	  fs::path command_script_;

	  void replace_default_parameters(std::string &str);
	  
     };

     bool call_BART(const std::string &cmdline);

     // Read BART files     
     template <typename int_t>
     std::vector<int_t> read_BART_hdr(fs::path &filename)
     {
	  std::vector<int_t> DIMS;
	  
	  std::ifstream infile((filename += ".hdr").c_str());
	  if (infile)
	  {
	       std::vector<std::string> tokens;
	       std::string line;

	       while (std::getline(infile, line, '\n'))
	       {
		    tokens.push_back(line);
	       }

	       // Parse the dimensions
	       std::stringstream ss(tokens[1]);
	       std::string items;
	       while (getline(ss, items, ' ')) {
		    DIMS.push_back(std::stoi(items));
	       }
	  }
	  else
	  {
	       GERROR("Failed to open file: %s\n", filename.c_str());
	  }

	  return DIMS;
     }

     template <typename int_t>
     std::pair<std::vector<int_t>, std::vector<std::complex<float>>>
     read_BART_files(fs::path &filename)
     {
	  // Load the header file
	  auto DIMS = read_BART_hdr<int_t>(filename);

	  // Load the cfl file
	  std::ifstream infile((filename += ".cfl").c_str(), std::ifstream::binary);
	  if (!infile)
	  {
	       GERROR("Failed to open file: %s\n", filename.c_str());
	       return std::make_pair(DIMS, std::vector<std::complex<float>>());
	  }
	  else
	  {
	       // First determine data size
	       infile.seekg(0, infile.end);
	       size_t size(infile.tellg());

	       // Now read data from file
	       infile.seekg(0);
	       std::vector<std::complex<float>> Data(size/sizeof(float)/2.);
	       infile.read(reinterpret_cast<char*>(Data.data()), size);

	       return std::make_pair(DIMS, Data);
	  }
     }
     
     // For convenience
     std::vector<size_t> read_BART_hdr(boost::filesystem::path &filename);
     std::pair< std::vector<size_t>, std::vector<std::complex<float> > > read_BART_files(fs::path &filename);

     // ------------------------------------------------------------------------
     // Write BART files
     template<typename int_t>
     void write_BART_hdr(fs::path filename, const std::vector<int_t>& DIMS)
     {
	  constexpr size_t MAX_DIMS = 16;
	  std::vector<size_t> v(MAX_DIMS, 1);
	  assert(DIMS.size() < MAX_DIMS);
	  std::copy(DIMS.cbegin(), DIMS.cend(), v.begin());
	  std::ofstream pFile((filename += ".hdr").c_str());
	  if (!pFile)
	       GERROR("Failed to write into file: %s\n", filename);
	  pFile << "# Dimensions\n";
	  std::copy(v.cbegin(), v.cend(), std::ostream_iterator<size_t>(pFile, " "));
     }

     template<typename int_t>
     void write_BART_Files(fs::path filename, const std::vector<int_t>& DIMS, const std::vector<std::complex<float>>& DATA)
     {
	  GDEBUG("Writing CFL file to %s\n", filename.c_str());
	  write_BART_hdr(filename, DIMS);
	  std::ofstream pFile((filename += ".cfl").c_str(),
			      std::ofstream::out | std::ofstream::binary);
	  if (!pFile)
	       GERROR("Failed to write into file: %s\n", filename);
	  pFile.write(reinterpret_cast<const char*>(&DATA[0]), DATA.size() * sizeof(std::complex<float>));
     }
     
     template<typename int_t>
     void write_BART_Files(fs::path filename, const std::vector<int_t>&DIMS, const hoNDArray<std::complex<float>>& DATA)
     {
	  GDEBUG("Writing CFL file to %s\n", filename.c_str());
	  write_BART_hdr(filename, DIMS);
	  std::ofstream pFile((filename += ".cfl").c_str(),
			      std::ofstream::out | std::ofstream::binary);
	  if (!pFile)
	       GERROR("Failed to write into file: %s\n", filename);
	  pFile.write(reinterpret_cast<const char*>(&DATA[0]), DATA.get_number_of_bytes());
     }
} // namespace Gadgetron
#endif //BART_GADGET_H
