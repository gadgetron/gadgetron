/**
 * Data format for CBCT acquisition (data and geometry)
 **/

#pragma once

#include "vector_td_io.h"
#include "hoNDArray.h"
#include "hoNDArray_utils.h"

#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <boost/shared_ptr.hpp>

namespace Gadgetron{

class CBCT_geometry
{
public:

	CBCT_geometry() {
		SAD_ = 0.0f;
		SDD_ = 0.0f;
		FOV_ = floatd2(0.0f);
	}

	~CBCT_geometry() {}

	inline void set_SAD( float SAD ) { SAD_ = SAD; }
	inline float get_SAD() { return SAD_; }

	inline void set_SDD( float SDD ) { SDD_ = SDD; }
	inline float get_SDD() { return SDD_; }

	inline void set_FOV( floatd2 v ) { FOV_ = v; }
	inline floatd2 get_FOV() { return FOV_; }

	inline void set_angles( std::vector<float> &angles ) { angles_ = angles; }
	inline std::vector<float>& get_angles() { return angles_; }

	inline void set_offsets( std::vector<floatd2> &offsets ) { offsets_ = offsets; }
	inline std::vector<floatd2>& get_offsets() { return offsets_; }

	// Basic output support
	//

	void print( std::ostream& os )
	{
		os << "------------ GEOMETRY ------------" << std::endl;
		if (angles_.size() == 0)
			os << "Angles: " << "EMPTY" << std::endl;
		else {
			os << "Angles: ";
			os << "Angles: " << angles_.front() << " ... " << angles_.back()
	  										 << ", number of angles: " << angles_.size() << std::endl;
		}

		if (offsets_.size() == 0)
			os << "Offsets: " << "EMPTY" << std::endl;
		else {
			os << "Offsets: contains " << offsets_.size() << " elements" << std::endl;
		}

		os << "SDD: " << SDD_ << "mm" << std::endl;
		os << "SAD: " << SAD_ << "mm" << std::endl;
		os << "FOV: " << FOV_ << "mm" << std::endl;
		os << "----------------------------------" << std::endl;
	}

	void save( hid_t file_id )
	{
		{
			unsigned int dataformat_version=2;
			hsize_t dims[1] = {1};
			H5LTmake_dataset(file_id, "/geometry_dataformat_version", 1, dims, H5T_NATIVE_UINT, &dataformat_version);
		}
		{
			hsize_t dims[1] = {1};
			H5LTmake_dataset(file_id, "/SAD", 1, dims, H5T_NATIVE_FLOAT, &SAD_);
		}
		{
			hsize_t dims[1] = {1};
			H5LTmake_dataset(file_id, "/SDD", 1, dims, H5T_NATIVE_FLOAT, &SDD_);
		}
		{
			hsize_t dims[1] = {2};
			H5LTmake_dataset(file_id, "/FOV", 1, dims, H5T_NATIVE_FLOAT, &FOV_);
		}
		{
			hsize_t dims[1] = {angles_.size()};
			H5LTmake_dataset(file_id, "/angles", 1, dims, H5T_NATIVE_FLOAT, &angles_[0]);
		}
		{
			std::vector<float> offsetx, offsety;
			for( unsigned int i=0; i<offsets_.size(); i++ ){
				floatd2 offset = offsets_[i];
				offsetx.push_back(offset[0]);
				offsety.push_back(offset[1]);
			}
			hsize_t dims[1] = {offsets_.size()};
			H5LTmake_dataset(file_id, "/offsetx", 1, dims, H5T_NATIVE_FLOAT, &offsetx[0]);
			H5LTmake_dataset(file_id, "/offsety", 1, dims, H5T_NATIVE_FLOAT, &offsety[0]);
		}
	}

protected:

	float SAD_;
	float SDD_;
	floatd2 FOV_;
	std::vector<float> angles_;
	std::vector<floatd2> offsets_;
};

class CBCT_acquisition {

public:

	CBCT_acquisition() {}

	CBCT_acquisition( boost::shared_ptr< hoNDArray<float> > projections,
			boost::shared_ptr<CBCT_geometry> geometry )
	{
		geometry_ = geometry;
		projections_ = projections;
	}

	virtual ~CBCT_acquisition() {}

	inline void set_geometry( boost::shared_ptr<CBCT_geometry> geometry ) {
		geometry_ = geometry;
	}

	inline boost::shared_ptr<CBCT_geometry> get_geometry() {
		return geometry_; }

	inline void set_projections( boost::shared_ptr< hoNDArray<float> > projections ) {
		projections_ = projections;
	}

	inline boost::shared_ptr< hoNDArray<float> > get_projections() {
		return projections_;
	}

	void downsample( unsigned int num_downsamples )
	{
		for (int k = 0; k < num_downsamples; k++)
			projections_ = boost::make_shared<hoNDArray<float>>(Gadgetron::downsample<float,2>(projections_.get()));
	}

	void load( std::string filename )
	{
		// Open hdf5 file
		//

		hid_t file_id = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

		// Load geometry.
		// This loader is for version 2 of the format
		//

		unsigned int geom_dataformat_version;
		herr_t errCode;
		errCode = H5LTread_dataset (file_id, "/geometry_dataformat_version", H5T_NATIVE_UINT, &geom_dataformat_version);

		if (errCode < 0){
			throw std::runtime_error("Error reading /geometry_dataformat_version from file.");
		}
		unsigned int needed = 2;
		if (!(geom_dataformat_version == needed)) {
			std::stringstream ss;
			ss << "wrong geometry data version inside hdf5 file, found: "
					<< geom_dataformat_version << ", needed: " << needed << std::endl;
			throw std::runtime_error(ss.str());
		}

		// Allocate new geometry
		//

		geometry_ = boost::shared_ptr<CBCT_geometry>(new CBCT_geometry());

		// Get angles array
		//

		hsize_t dim;
		errCode = H5LTget_dataset_info(file_id,"/angles",&dim,NULL,NULL);
		if (errCode < 0) 	throw std::runtime_error("Error getting /angles dataset info from file.");

		std::vector<float> angles (dim,0.0f);
		geometry_->set_angles(angles);
		errCode=H5LTread_dataset (file_id, "/angles", H5T_NATIVE_FLOAT, &geometry_->get_angles()[0]);
		if (errCode < 0) 	throw std::runtime_error("Error reading /angles from file.");

		// Get offsets array
		//

		errCode=H5LTget_dataset_info(file_id,"/offsetx",&dim,NULL,NULL);
		if (errCode < 0) 	throw std::runtime_error("Error getting /offsetx dataset info from file.");
		std::vector<float> offsets_x = std::vector<float>(dim,0.0f);
		errCode=H5LTread_dataset (file_id, "/offsetx", H5T_NATIVE_FLOAT, &offsets_x[0]);
		if (errCode < 0) 	throw std::runtime_error("Error reading /offsetx from file.");

		errCode=H5LTget_dataset_info(file_id,"/offsety",&dim,NULL,NULL);
		if (errCode < 0) 	throw std::runtime_error("Error getting /offsety dataset info from file.");
		std::vector<float> offsets_y = std::vector<float>(dim,0.0f);
		errCode=H5LTread_dataset (file_id, "/offsety", H5T_NATIVE_FLOAT, &offsets_y[0]);
		if (errCode < 0) 	throw std::runtime_error("Error reading /offsety from file.");

		if( offsets_x.size() != offsets_y.size() ){
			throw std::runtime_error("CBCT_geometry::load : x/y offset arrays has different lengths");
		}

		geometry_->get_offsets().clear();
		for( unsigned int i=0; i<offsets_x.size(); i++ ){
			geometry_->get_offsets().push_back(floatd2(offsets_x[i], offsets_y[i]));
		}

		// Test data format of the projections
		//

		unsigned int proj_dataformat_version;
		errCode=H5LTread_dataset (file_id, "/projection_dataformat_version", H5T_NATIVE_UINT, &proj_dataformat_version);
		if (errCode < 0) 	throw std::runtime_error("Error reading /projection_dataformat_version from file.");

		needed = 1;
		if (!(proj_dataformat_version == needed)) {
			std::stringstream ss;
			ss << "wrong projection data format version inside hdf5 file, found: "
					<< proj_dataformat_version << ", needed: " << needed;
			throw std::runtime_error(ss.str());
		}

		hsize_t vec_dim[3];
		errCode=H5LTget_dataset_info(file_id,"/projections",vec_dim,NULL,NULL);
		if (errCode < 0) 	throw std::runtime_error("Error getting /projections dataset info from file.");
		std::vector<size_t> dims;
		dims.push_back(vec_dim[2]);
		dims.push_back(vec_dim[1]);
		dims.push_back(vec_dim[0]);

		projections_ = boost::shared_ptr<hoNDArray<float> >(new hoNDArray<float>(&dims));
		errCode=H5LTread_dataset (file_id,"/projections", H5T_NATIVE_FLOAT, projections_->get_data_ptr());
		if (errCode < 0) 	throw std::runtime_error("Error reading /projections from file.");

		// Get SAD / SDD / FOV
		//

		float SAD, SDD;
		floatd2 FOV;

		errCode=H5LTread_dataset (file_id, "/SAD", H5T_NATIVE_FLOAT, &SAD);
		if (errCode < 0) 	throw std::runtime_error("Error reading /SAD from file.");
		errCode=H5LTread_dataset (file_id, "/SDD", H5T_NATIVE_FLOAT, &SDD);
		if (errCode < 0) 	throw std::runtime_error("Error reading /SDD from file.");
		errCode=H5LTread_dataset (file_id, "/FOV", H5T_NATIVE_FLOAT, &FOV);
		if (errCode < 0){
			floatd2 spacing;
			errCode=H5LTread_dataset (file_id, "/spacing", H5T_NATIVE_FLOAT, &spacing);
			FOV[0] = spacing[0]*dims[0];
			FOV[1] = spacing[1]*dims[1];
			if (errCode < 0) throw std::runtime_error("Error reading /FOV from file.");
		}

		geometry_->set_SAD(SAD);
		geometry_->set_SDD(SDD);
		geometry_->set_FOV(FOV);
		H5Fclose (file_id);
	}

	void save( std::string filename )
	{
		hid_t file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

		unsigned int dataformat_version=1;
		hsize_t dims[1] = {1};
		H5LTmake_dataset(file_id,"/projection_dataformat_version", 1, dims, H5T_NATIVE_UINT, &dataformat_version);

		boost::shared_ptr<std::vector<size_t > > pdims = projections_->get_dimensions();
		hsize_t *dims2 = new hsize_t[pdims->size()];
		for (int i = 0; i < pdims->size(); i++)
			dims2[i] = pdims->at(pdims->size()-i-1);
		H5LTmake_dataset(file_id,"/projections", pdims->size(), dims2, H5T_NATIVE_FLOAT, projections_->get_data_ptr());
		delete[] dims2;

		geometry_->save(file_id);

		H5Fclose (file_id);
	}

protected:

	boost::shared_ptr<CBCT_geometry> geometry_;
	boost::shared_ptr< hoNDArray<float> > projections_;
};
}
