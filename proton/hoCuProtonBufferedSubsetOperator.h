#pragma once
#include "subsetOperator.h"
#include "hoCuNDArray.h"
#include "hoCuOperatorPathBackprojection.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <numeric>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include "GadgetronException.h"

#include <boost/thread.hpp>
/**
 * This class is inconsistent, messy and requires polish. It's also suspected to be unreliable
 */
namespace Gadgetron{
template<class REAL> class hoCuProtonBufferedSubsetOperator : public subsetOperator<hoCuNDArray<REAL> >{

	typedef hoCuProtonBufferedSubsetOperator<REAL> Type;
public:
	hoCuProtonBufferedSubsetOperator() : subsetOperator<hoCuNDArray<REAL> >(1), loaded_subset(-1), loading_subset(-1){
		this->subset_dimensions = std::vector<std::vector<unsigned int> > (1, std::vector<unsigned int>(1) );
	}
	hoCuProtonBufferedSubsetOperator(int subsets) : subsetOperator<hoCuNDArray<REAL> >(subsets), loaded_subset(-1), loading_subset(-1){
		this->subset_dimensions = std::vector<std::vector<unsigned int> > (subsets, std::vector<unsigned int>(1) );
	}

	virtual ~hoCuProtonBufferedSubsetOperator(){};
	//Not functional nor elegant... please rewrite when less tired
	boost::shared_ptr<hoCuNDArray<REAL> > load_data(std::string filename,vector_td<REAL,3> _physical_dims,vector_td<REAL,3> _origin, REAL background=0){
		physical_dims = _physical_dims;
		origin=_origin;
		data_file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		groupnames = group_paths("/",data_file);
		calculate_subset_dimensions(data_file,groupnames);
		size_t elements = 0;
		for (int i = 0; i < this->number_of_subsets; i++)
			elements += std::accumulate(this->subset_dimensions[i].begin(),this->subset_dimensions[i].end(),1,std::multiplies<unsigned int>());

		boost::shared_ptr<hoCuNDArray<REAL> > projections = loadProjections(data_file,groupnames, elements);

		return projections;


	}

	void loadSubset(int subset,hoCuNDArray<REAL>* projections){
		if (subset != loaded_subset){

			if(subset == loading_subset){
				loader.join();
			}else
				loadSplinesGroup(data_file,groupnames,subset);
			splines_lock.lock();
			loaded_splines = loading_splines;
			splines_lock.unlock();
			loaded_subset = subset;
			loader = boost::thread(&Type::loadSplinesGroup,this,data_file,groupnames,(subset+1)%this->number_of_subsets);
			loading_subset = (subset+1)%this->number_of_subsets;
			op = hoCuOperatorPathBackprojection<REAL>();
			boost::shared_ptr<hoCuNDArray<REAL> > fake_projections(new hoCuNDArray<REAL>(projections->get_dimensions()));
			op.setup(loaded_splines,physical_dims,fake_projections,origin,REAL(0));
		}

	}


	virtual void mult_M(hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, int subset, bool accumulate){
		std::stringstream ss;
		ss << "Subset " << subset << " out of bounds";
		if (subset >= this->number_of_subsets) BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
		loadSubset(subset,out);
		op.mult_M(in,out,accumulate);

	}
	virtual void mult_MH(hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, int subset, bool accumulate){
			std::stringstream ss;
			ss << "Subset " << subset << " out of bounds";
			if (subset >= this->number_of_subsets ) BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
			loadSubset(subset,in);
			op.mult_MH(in,out,accumulate);
	}
	virtual void mult_MH_M(hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, int subset, bool accumulate){
				BOOST_THROW_EXCEPTION(runtime_error("Rather embarrassingly mult_MH_M doesn't work with buffering right now."));
	}

	using subsetOperator<hoCuNDArray<REAL> >::mult_M;
	using subsetOperator<hoCuNDArray<REAL> >::mult_MH;
	using subsetOperator<hoCuNDArray<REAL> >::mult_MH_M;
  virtual boost::shared_ptr< linearOperator<hoCuNDArray<REAL> > > clone() {
  	BOOST_THROW_EXCEPTION(runtime_error("Due to dependency on boost thread."));
  	return boost::shared_ptr< linearOperator<hoCuNDArray<REAL> > >();
	 }

protected:

	struct Spline{
		float x,y,z,x2,y2,z2;
		float dirx,diry,dirz,dirx2,diry2,dirz2;
	};


	boost::shared_ptr<hoCuNDArray<REAL> > loadProjections(hid_t file_id,std::vector<std::string>& groupnames, size_t elements){

				std::vector<unsigned int> projection_dim(1,elements);

				boost::shared_ptr<hoCuNDArray<REAL> > projections(new hoCuNDArray<REAL>(&projection_dim));
				hid_t strtype;                     /* Datatype ID */
				herr_t status;
				const size_t float_size[1] = {sizeof(float) };
				const size_t float_offset[1] = {0};
				std::string projections_name = "projections";

				std::vector< boost::shared_ptr<hoCuNDArray<REAL> > > tmp_proj = this->projection_subsets(projections.get());
				std::vector<REAL*> ptrs;
				for (int i = 0; i< tmp_proj.size(); i++) ptrs.push_back(tmp_proj[i]->get_data_ptr());

				REAL* ptr = projections->get_data_ptr();
				for (int i = 0; i < groupnames.size(); i++){
					hsize_t nfields,nrecords;
					herr_t err = H5TBget_table_info (file_id, (groupnames[i]+projections_name).c_str(), &nfields, &nrecords );
					if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Illegal hdf5 dataset provided"));
					size_t extra = nrecords%this->number_of_subsets;
					hsize_t offset = 0;
					for (int subset =0; subset < this->number_of_subsets; subset++){
						hsize_t batchSize = nrecords/this->number_of_subsets;
						if (subset < extra)  batchSize += 1;
						err = H5TBread_records (file_id, (groupnames[i]+projections_name).c_str(), offset, batchSize, sizeof(float),  float_offset, float_size,  ptrs[subset] );
						if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Unable to read splines from hdf5 file"));
						offset += batchSize;
						ptrs[subset] += batchSize;
					}

				}
				return projections;
		}

	void loadSplinesGroup(hid_t file_id, std::vector<std::string>& groupnames, int subset){


		const size_t dst_sizes[12] = { sizeof(float) ,sizeof(float), sizeof(float),
				 	 	 	 	 	 	 	 	 sizeof(float) ,sizeof(float), sizeof(float),
				 	 	 	 	 	 	 	 	 sizeof(float) ,sizeof(float), sizeof(float),
				 	 	 	 	 	 	 	 	 sizeof(float) ,sizeof(float), sizeof(float)};
		const size_t dst_offset[12] = { HOFFSET( Spline, x ),HOFFSET( Spline, y),HOFFSET( Spline, z ),
				HOFFSET( Spline, x2 ),HOFFSET( Spline, y2),HOFFSET( Spline, z2 ),
	  		HOFFSET( Spline, dirx ),HOFFSET( Spline, diry ),HOFFSET( Spline, dirz ),
	  		HOFFSET( Spline, dirx2 ),HOFFSET( Spline, diry2 ),HOFFSET( Spline, dirz2 )};

		std::vector<unsigned int> spline_dims;
		spline_dims.push_back(4);
		spline_dims.push_back(std::accumulate(this->subset_dimensions[subset].begin(),this->subset_dimensions[subset].end(),1,std::multiplies<unsigned int>()));

		splines_lock.lock();
		loading_splines = boost::shared_ptr<hoCuNDArray<vector_td<REAL,3> > >(new hoCuNDArray<vector_td<REAL,3> >(&spline_dims));
		std::string splines_name = "splines";

		vector_td<REAL,3>* splinePtr = loading_splines->get_data_ptr();

		hid_t strtype;                     /* Datatype ID */
		herr_t status;


		for (int i = 0; i < groupnames.size(); i++){
			hsize_t nfields,nrecords;
			herr_t err = H5TBget_table_info (file_id, (groupnames[i]+splines_name).c_str(), &nfields, &nrecords );
			if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Illegal hdf5 dataset provided"));
			int extra = nrecords%this->number_of_subsets;

			hsize_t batchSize = nrecords/this->number_of_subsets;
			hsize_t offset = batchSize*subset+std::min(subset,extra);
			if (subset < extra ) batchSize += 1;
			err = H5TBread_records (file_id, (groupnames[i]+splines_name).c_str(), offset, batchSize, sizeof(Spline),  dst_offset, dst_sizes,  splinePtr );
			if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Unable to read splines from hdf5 file"));
			splinePtr += batchSize*sizeof(Spline)/sizeof(vector_td<REAL,3>);
			}
		splines_lock.unlock();
	}
	void calculate_subset_dimensions(hid_t file_id, std::vector<std::string>& groupnames){

		std::string projection_name = "projections";
		for (int i = 0; i < groupnames.size(); i++){
			hsize_t nfields,nrecords;
			herr_t err = H5TBget_table_info (file_id, (groupnames[i]+projection_name).c_str(), &nfields, &nrecords );
			if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Illegal hdf5 dataset provided"));
			size_t extra = nrecords%this->number_of_subsets;
			for (int subset =0; subset < this->number_of_subsets; subset++){
				this->subset_dimensions[subset][0] += nrecords/this->number_of_subsets;
				if (subset < extra)  this->subset_dimensions[subset][0] += 1;
			}

		}
	}

	std::vector<std::string> group_paths(std::string path,hid_t file_id){

		char node[2048];
		hsize_t nobj,len;
		herr_t err;
		hid_t group_id = H5Gopen1(file_id,path.c_str());

		err = H5Gget_num_objs(group_id, &nobj);

		std::vector<std::string> result;
		for(hsize_t i =0; i < nobj; i++){
			len = H5Gget_objname_by_idx(group_id, i,
						node, sizeof(node) );
			std::string nodestr = std::string(path).append(node).append("/");
			int otype =  H5Gget_objtype_by_idx(group_id, i );
			switch(otype){
			case H5G_GROUP:
				//cout << nodestr << " is a GROUP" << endl;
				result.push_back(nodestr);
				break;
			}

		}
		H5Gclose(group_id);
		return result;

	}

	hoCuOperatorPathBackprojection<REAL> op;

	boost::thread loader;

	int loaded_subset;
	int loading_subset;

	boost::mutex splines_lock;
	boost::shared_ptr<hoCuNDArray<vector_td<REAL,3> > > loading_splines;

	boost::shared_ptr<hoCuNDArray<vector_td<REAL,3> > > loaded_splines;
	vector_td<REAL,3> physical_dims;
	vector_td<REAL,3> origin;

	hid_t data_file;
	std::vector<std::string> groupnames;
};
}
