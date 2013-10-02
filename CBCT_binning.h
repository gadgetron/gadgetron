/**
 * Temporal binning for CBCT
 **/

#pragma once

#include "CBCT_acquisition.h"

#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>
#include <stdexcept>
#include <boost/iterator/counting_iterator.hpp>

namespace Gadgetron {

  class CBCT_binning
  {

  public:
    
    CBCT_binning() {}
    CBCT_binning( std::vector< std::vector<unsigned int> > binning ) : binning_(binning) {}

    ~CBCT_binning() {}

    inline unsigned int get_number_of_bins()
    {
      return binning_.size();
    }

    inline unsigned int get_number_of_projections()
    {
      unsigned int acc = 0;
      for( unsigned int i=0; i<binning_.size(); i++ )
	acc += binning_[i].size();
      return acc;
    }

    inline unsigned int get_number_of_projections( unsigned int bin )
    {
      if( bin >= binning_.size() )
	throw std::runtime_error("CBCT_binning::get_number_of_projections(int) : bin is out of range");
      else
	return binning_[bin].size();
    }

    inline int get_maximum_projection_index()
    {
      int max_proj = -1;
      for( unsigned int i=0; i<binning_.size(); i++ )
	for( unsigned int j=0; j<binning_[i].size(); j++ )
	  if( int(binning_[i][j]) > max_proj ) 
	    max_proj = binning_[i][j];
      return max_proj;
    }
    
    inline void set_bins( std::vector< std::vector<unsigned int> > &bins ) {
      binning_ = bins;
    }

    inline std::vector< std::vector<unsigned int> > get_bins() {
      return binning_;
    }

    inline void set_bin( std::vector<unsigned int> &bin, unsigned int bin_number )
    {
      if( bin_number > binning_.size() )
	throw std::runtime_error("CBCT_binning::set_bin() : bin is out of range");
      else if( bin_number == binning_.size() )
	binning_.push_back(bin);
      else
	binning_[bin_number] = bin;
    }

    inline std::vector<unsigned int> get_bin( unsigned int bin )
    {
      if( bin >= binning_.size() )
	throw std::runtime_error("CBCT_binning::get_bin() : bin is out of range");
      else
	return binning_[bin];
    }

    inline void set_as_default_3d_bin( unsigned int num_projections )
    {
      binning_.push_back( std::vector<unsigned int>( boost::counting_iterator<unsigned int>(0),
						     boost::counting_iterator<unsigned int>(num_projections) ));
    }

    void load( std::string filename )
    {
      // Open file and make sure it is the expected version
      //

      hid_t file_id = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

      unsigned int dataformat_version;
      herr_t errCode;
      errCode=H5LTread_dataset (file_id, "/binning_dataformat_version", H5T_NATIVE_UINT, &dataformat_version);

      if(errCode < 0)
	throw std::runtime_error("Error reading /binning_dataformat_version");

      unsigned int needed = 1;
      if (!(dataformat_version == needed)) {
	std::cerr << "wrong format of binning file, found: "
		  << dataformat_version << ", needed: " << needed << std::endl;
	exit(EXIT_FAILURE);
      }

      // And get the bins
      //

      binning_.clear();

      unsigned int numBins;
      errCode=H5LTread_dataset (file_id, "/numBins", H5T_NATIVE_UINT, &numBins);
      if(errCode < 0)
	throw std::runtime_error("Error reading /numBins_dataformat_version");
      //std::cout << "Found " << numBins << " bins in file" << filename << std::endl;

      // Ok, so this really isn't very elegant.
      // A folder in the hdf5 file containing the data would be better...
      //

      for (unsigned int i=1; i<=numBins; i++) {
	std::stringstream path;
	path << "/bin_" << i;
	hsize_t dim;
	errCode=H5LTget_dataset_info(file_id,path.str().c_str(),&dim,NULL,NULL);
	if(errCode < 0)
	  throw std::runtime_error("Error reading bin info");
	binning_.push_back(std::vector<unsigned int>(dim,0.0f));
	errCode=H5LTread_dataset (file_id, path.str().c_str(), H5T_NATIVE_UINT, &binning_.back()[0]);
	if(errCode < 0)
	  throw std::runtime_error("Error reading bin data");
      }
    }

    void print( std::ostream &os = std::cout )
    {
      os << "---------- BINNING DATA ----------" << std::endl;
      os << "Number of bins: " << binning_.size() << std::endl;
      for (unsigned int b=0; b<binning_.size(); b++)
	os << "Number of projections in bin[" << b
	   << "]: " << binning_[b].size() << std::endl;
      os << "----------------------------------" << std::endl;
    }

  protected:
    std::vector< std::vector<unsigned int> > binning_;
  };
}

