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
/**
 * This class is inconsistent, messy and requires polish. It's also suspected to be unreliable
 */
namespace Gadgetron{
template<class REAL> class hoCuProtonSubsetOperator : public subsetOperator<hoCuNDArray<REAL> >{

public:
	hoCuProtonSubsetOperator() : subsetOperator<hoCuNDArray<REAL> >(1){
		this->subset_dimensions = std::vector<std::vector<unsigned int> > (1, std::vector<unsigned int>(1) );
	}
	hoCuProtonSubsetOperator(int subsets) : subsetOperator<hoCuNDArray<REAL> >(subsets){
		this->subset_dimensions = std::vector<std::vector<unsigned int> > (subsets, std::vector<unsigned int>(1) );
	}

	virtual ~hoCuProtonSubsetOperator(){};
	//Not functional nor elegant... please rewrite when less tired
	boost::shared_ptr<hoCuNDArray<REAL> > load_data(std::string filename,vector_td<REAL,3> physical_dims,vector_td<REAL,3> origin, REAL background=0){
		operators.clear();
		spline_arrays.clear();

		hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		std::vector<std::string> groupnames = group_paths("/",file_id);
		calculate_subset_dimensions(file_id,groupnames);
		size_t elements = 0;
		for (int i = 0; i < this->number_of_subsets; i++)
			elements += std::accumulate(this->subset_dimensions[i].begin(),this->subset_dimensions[i].end(),1,std::multiplies<unsigned int>());

		loadSplines(file_id,groupnames,elements);

		boost::shared_ptr<hoCuNDArray<REAL> > projections = loadProjections(file_id,groupnames, elements);

		std::vector<boost::shared_ptr<hoCuNDArray<REAL> > > sub_projections = this->projection_subsets(projections.get());
		for (int i = 0; i < this->number_of_subsets; i++){
			operators.push_back(hoCuOperatorPathBackprojection<REAL>());
			operators.back().setup(spline_arrays[i],physical_dims,sub_projections[i],origin,background);
			std::cout << "Projections: " << i << " " << nrm2(sub_projections[i].get()) << " " << *sub_projections[i]->begin() << std::endl;
		}
		vector_td<REAL,3>* test = splines->get_data_ptr();
		std::cout << "First spline " << test[0] <<" " << test[1] << " " << test[2] << " " << test[3] << std::endl;


		return projections;


	}


	virtual void mult_M(hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, int subset, bool accumulate=false){
		std::stringstream ss;
		ss << "Subset " << subset << " out of bounds";
		if (subset >= operators.size() ) BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
		operators[subset].mult_M(in,out,accumulate);
	}
	virtual void mult_MH(hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, int subset, bool accumulate=false){
			std::stringstream ss;
			ss << "Subset " << subset << " out of bounds";
			if (subset >= operators.size() ) BOOST_THROW_EXCEPTION(runtime_error(ss.str()));
			operators[subset].mult_MH(in,out,accumulate);
	}
	virtual void mult_MH_M(hoCuNDArray<REAL>* in, hoCuNDArray<REAL>* out, int subset, bool accumulate=false){
				if (subset >= operators.size() ) BOOST_THROW_EXCEPTION(runtime_error("Subset out of bounds"));
				operators[subset].mult_MH_M(in,out,accumulate);
	}
  virtual boost::shared_ptr< linearOperator<hoCuNDArray<REAL> > > clone() {
       return linearOperator<hoCuNDArray<REAL> >::clone(this);
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
					err = H5TBread_records (file_id, (groupnames[i]+projections_name).c_str(), offset, batchSize, sizeof(float),  float_offset, float_size,  ptr );
					if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Unable to read splines from hdf5 file"));
					offset += batchSize;
					ptr += batchSize;
				}

			}
			return projections;
	}

	void loadSplines(hid_t file_id, std::vector<std::string>& groupnames, size_t elements){

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
		spline_dims.push_back(elements);

		splines = boost::shared_ptr<hoCuNDArray<vector_td<REAL,3> > >(new hoCuNDArray<vector_td<REAL,3> >(&spline_dims));
		std::string splines_name = "splines";

		vector_td<REAL,3>* splinePtr = splines->get_data_ptr();
		for (int subset = 0; subset < this->number_of_subsets; subset++){
			std::vector<unsigned int> subspline_dim =this->subset_dimensions[subset];
			subspline_dim.insert(subspline_dim.begin(),4);

			boost::shared_ptr<hoCuNDArray<vector_td<REAL,3> > > subSpline(new hoCuNDArray<vector_td<REAL,3> >(&subspline_dim,splinePtr));
			spline_arrays.push_back(subSpline);
			splinePtr += std::accumulate(subspline_dim.begin(),subspline_dim.end(),1,std::multiplies<unsigned int>());
		}

		hid_t strtype;                     /* Datatype ID */
		herr_t status;

		splinePtr = splines->get_data_ptr();
		for (int i = 0; i < groupnames.size(); i++){
			hsize_t nfields,nrecords;
			herr_t err = H5TBget_table_info (file_id, (groupnames[i]+splines_name).c_str(), &nfields, &nrecords );
			if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Illegal hdf5 dataset provided"));
			size_t extra = nrecords%this->number_of_subsets;
			hsize_t offset = 0;
			for (int subset =0; subset < this->number_of_subsets; subset++){
				hsize_t batchSize = nrecords/this->number_of_subsets;
				if (subset < extra)  batchSize += 1;
				err = H5TBread_records (file_id, (groupnames[i]+splines_name).c_str(), offset, batchSize, sizeof(Spline),  dst_offset, dst_sizes,  splinePtr );
				if (err < 0) BOOST_THROW_EXCEPTION(runtime_error("Unable to read splines from hdf5 file"));
				offset += batchSize;
				splinePtr += batchSize*sizeof(Spline)/sizeof(vector_td<REAL,3>);
			}

		}
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
		hid_t group_id = H5Gopen(file_id,path.c_str(),H5P_DEFAULT);

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

	std::vector< hoCuOperatorPathBackprojection<REAL> > operators;

	boost::shared_ptr<hoCuNDArray<vector_td<REAL,3> > > splines;
	std::vector< boost::shared_ptr<hoCuNDArray<vector_td<REAL,3> > > > spline_arrays;
};
}
