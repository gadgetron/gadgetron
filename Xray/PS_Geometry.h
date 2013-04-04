/**
 * Projection space data set
 **/
#pragma once

#include "hdf5.h"
#include "hdf5_hl.h"
#include "vector_td_io.h"
namespace Gadgetron {

class PS_Geometry {
 protected:
    float SAD;
    float SDD;
    std::vector<float> anglesArray;
    floatd2 spacingArray;
    floatd3 SAGxArray;
    floatd3 SAGyArray;

 public:
    void setSAD( float SAD ) { this->SAD = SAD; }
    void setSDD( float SDD ) { this->SDD = SDD; }

    std::vector<float>& getAnglesArray() { return anglesArray; }
    float getSAD() { return SAD; }
    float getSDD() { return SDD; }

    floatd2 getSpacingArray() { return spacingArray; }
    floatd3 getSAGxArray() { return SAGxArray; }
    floatd3 getSAGyArray() { return SAGyArray; }

    PS_Geometry() {
        spacingArray = floatd2(0);
        SAGxArray = floatd3(0);

        SAGyArray = floatd3(0);
        SAD = 1000;
        SDD = 1500;

    }

    ~PS_Geometry() {

    }

    void loadData(std::string inFile) {
    	hid_t file_id = H5Fopen (inFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    	// test header field

			unsigned int dataformat_version;
			H5LTread_dataset (file_id, "/geometry_dataformat_version", H5T_NATIVE_UINT, &dataformat_version);


			unsigned int needed = 1;
			if (!(dataformat_version == needed)) {
					std::stringstream ss;
					ss << "wrong geometry data format version inside hdf5 file, found: "
						 << dataformat_version << ", needed: " << needed << std::endl;
					throw std::runtime_error(ss.str());
			}
			hsize_t dim;
			H5LTget_dataset_info(file_id,"/angles",&dim,NULL,NULL);
			anglesArray = std::vector<float>(dim,0.0f);
			H5LTread_dataset (file_id, "/angles", H5T_NATIVE_FLOAT, &anglesArray[0]);
			H5LTread_dataset (file_id, "/SAD", H5T_NATIVE_FLOAT, &SAD);
			H5LTread_dataset (file_id, "/SDD", H5T_NATIVE_FLOAT, &SDD);
			H5LTread_dataset (file_id, "/spacing", H5T_NATIVE_FLOAT, &spacingArray[0]);
			H5LTread_dataset (file_id, "/SAGx", H5T_NATIVE_FLOAT, &SAGxArray[0]);
			H5LTread_dataset (file_id, "/SAGy", H5T_NATIVE_FLOAT, &SAGyArray[0]);



    }
/*
    void saveData(std::string outFile) {
        if (!File::Exists(outFile)) {
            createEmptyFile(outFile);
        }
        appendVariable(outFile, "/geometry_dataformat_version", 1);
        if(anglesArray == NULL) {
            std::stringstream ss;
            ss << "could not save angles because of NULL pointer, while writing file: " << outFile;
            throw std::runtime_error(ss.str());
        }
        appendTypedArray1D(outFile,  "/angles", *anglesArray );

        appendVariable(outFile, "/SAD", SAD);
        appendVariable(outFile, "/SDD", SDD);

        appendTypedArray1D(outFile, "/spacing", *spacingArray );
        appendTypedArray1D(outFile, "/SAGx", *SAGxArray );
        appendTypedArray1D(outFile, "/SAGy", *SAGyArray );
    }
    */

    void setSAGx(float vx, float vy, float vz) {
        SAGxArray = floatd3(vx,vy,vz);
    }

    void setSAGy(float vx, float vy, float vz) {
    		SAGyArray = floatd3(vx,vy,vz);

    }

    void setSpacing(float vx, float vy) {
    		spacingArray = floatd2(vx,vy);
    }

    void setAngles(std::vector<float> v) {
    		anglesArray = v;
    }

    void setAngles(float* a, unsigned int size) {
    		anglesArray = std::vector<float>(a,a+size);
    }

    void print(std::ostream& os) {
        os << "------------ GEOMETRY ------------" << std::endl;
        if (anglesArray.size() == 0)
            os << "angles: " << "EMPTY" << std::endl;
        else {
        	os << "Angles: ";
        		for (int i = 0; i < anglesArray.size(); i++) os << anglesArray[i] << std::endl;
            /*os << "angles: " << anglesArray[0] << " ... " << anglesArray.back()
               << ", number of angles: " << anglesArray.size() << std::endl;*/
        }
        os << "SDD: " << SDD << "mm" << std::endl;
        os << "SAD: " << SAD << "mm" << std::endl;
        os << "Spacing: " << spacingArray << "mm" << std::endl;
        os << "SAGx: " << SAGxArray << std::endl;
        os << "SAGy: " << SAGyArray << std::endl;
        os << "----------------------------------" << std::endl;
    }
};

}
