/**
 * Projection space data set
 **/

#pragma once

#include "PS_Geometry.h"
#include <vector>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <boost/iterator/counting_iterator.hpp>
namespace Gadgetron {


class PS_BinningData {
 protected:
    std::vector< std::vector<unsigned int> > binningData;
    bool clean;

 public:

    PS_BinningData() {

       clean = true;
    }

    PS_BinningData(std::vector< std::vector<unsigned int> >& bd) : binningData(bd) {
        clean = false;
    }

    ~PS_BinningData() {
    }

    void loadData(std::string inFile) {
    	hid_t file_id = H5Fopen (inFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        // test header field
    	unsigned int dataformat_version;
			H5LTread_dataset (file_id, "/binning_dataformat_version", H5T_NATIVE_UINT, &dataformat_version);

			unsigned int needed = 1;
			if (!(dataformat_version == needed)) {
					std::cerr << "wrong binning data format version inside hdf5 file, found: "
										<< dataformat_version << ", needed: " << needed << std::endl;
					exit(EXIT_FAILURE);
			}
			binningData=std::vector< std::vector<unsigned int> >();
			unsigned int numBins;
			H5LTread_dataset (file_id, "/binning_dataformat_version", H5T_NATIVE_UINT, &numBins);
			//Ok, so this really isn't very elegant. A folder in the hdf5 file containing the data would be better
			for (unsigned int i=1; i<=numBins; i++) {
					std::stringstream path;
					path << "/bin_" << i;
					hsize_t dim;
					H5LTget_dataset_info(file_id,path.str().c_str(),&dim,NULL,NULL);
					binningData.push_back(std::vector<unsigned int>(dim,0.0f));
					H5LTread_dataset (file_id, path.str().c_str(), H5T_NATIVE_UINT, &binningData.back()[0]);
			}
    }

    void generateData(PS_Geometry* ps_g) {
        unsigned int numProjections = ps_g->getAnglesArray().size();
        binningData.push_back(std::vector<unsigned int>(boost::counting_iterator<unsigned int>(1),
        		boost::counting_iterator<unsigned int>(numProjections+1)));
    } 
/*
    void loadDataDir(std::string inDir) {
        std::string num_bins_filename = inDir + "/numBins.txt";
        std::vector<unsigned int> numBinsVec = file2vector<unsigned int>(num_bins_filename);
        assert(numBinsVec.size() == 1);

        for (unsigned int p=1; p<=numBinsVec[0]; p++) {
            char number[6];
            sprintf(number, "%05i", p);
            std::string current_phase_filename = inDir + "/bin_"+ number + ".txt";

            std::vector<unsigned int> current_phase_vector = 
                file2vector<unsigned int>(current_phase_filename);
            Array1D<unsigned int>* array = new Array1D<unsigned int>(current_phase_vector);
            binningData->push_back(array);
        }
    }

    void saveData(std::string outFile) {
        if (!File::Exists(outFile)) {
            createEmptyFile(outFile);
        }
        appendVariable(outFile, "/binning_dataformat_version", 1);
            
        unsigned int numBins = binningData->size();
        appendVariable(outFile, "/numBins", numBins);
        for (unsigned int i=1; i<=numBins; i++) {
            std::stringstream path;
            path << "bin_" << i;
            Array1D<unsigned int>* array = (*binningData)[i-1];
            appendTypedArray1D(outFile, path.str(), *array );
            std::cout << "path: " << path.str() << std::endl;
        }
    }

 */
    unsigned int getNumberOfProjections() {
        unsigned int s = 0;
        for (unsigned int i= 0; i != binningData.size(); ++i)
            s += binningData[i].size();
        return s;
    }

    std::vector< std::vector<unsigned int> >& getBinningData() {
        return binningData;
    }

    void setBinningData(std::vector< std::vector<unsigned int> > b) {
       binningData = b;
    }

    void print(std::ostream& os) {
        os << "---------- BINNING DATA ----------" << std::endl;
        os << "Number of bins: " << binningData.size() << std::endl;
        for (unsigned int b=0; b<binningData.size(); b++)
            os << "Number of projections in bin[" << b 
               << "]: " << binningData[b].size() << std::endl;
        os << "----------------------------------" << std::endl;
    }
};
}
