/**
 * Projection space data set
 **/


#include "PS_Geometry.h"
#include "hoNDArray_fileio.h"
#include "hdf5.h"
#include "hdf5_hl.h"
//template <class T>
namespace Gadgetron{
class PS_Dataset {
    /*
    PSData(hoNDArray<T>* data, float* angles, float2 spacing, float2 origin, float SSD, float SAD)
     : data(data), angles(angles), spacing(spacing, origin(origin), SSD(SSD), SAD(SAD) {}
    */

 protected: 
    bool cleanup_projections;
    bool cleanup_geometry;

    PS_Geometry* geometry;
    boost::shared_ptr<hoNDArray<float> > projections;

 public:
    /*
    void setSAD( float SAD ) { geometry->setSAD(SAD); }
    void setSDD( float SDD ) { geometry->setSDD(SDD); }

    float getSAD() { return geometry->getSAD(); }
    float getSDD() { return geometry->getSDD(); }

    Array1D<float>* getSpacingArray() { return geometry->getSpacingArray(); }
    Array1D<float>* getSAGxArray() { return geometry->getSAGxArray(); }
    Array1D<float>* getSAGyArray() { return geometry->getSAGyArray(); }

    Array1D<float>* getAnglesArray() { return geometry->getAnglesArray(); }
    */
    boost::shared_ptr<hoNDArray<float> >  getProjections() { return projections; }

    PS_Dataset() {
        cleanup_projections = true;
        cleanup_geometry = true;
        geometry = NULL;
    }

    PS_Dataset(boost::shared_ptr<hoNDArray<float> > p) {
        cleanup_projections = false;
        cleanup_geometry = true;
        geometry = NULL;
        projections = p;
    }

    PS_Dataset(boost::shared_ptr<hoNDArray<float> > p, PS_Geometry* geo) {
        cleanup_projections = false;
        cleanup_geometry = false;
        geometry = geo;
        projections = p;
    }

    PS_Dataset(PS_Geometry* geo) {
        cleanup_projections = true;
        cleanup_geometry = false;
        geometry = geo;
    }

    virtual ~PS_Dataset() {
        if (cleanup_geometry)
            if (geometry != 0)
                delete geometry;
    }

    PS_Geometry* getGeometry() {
        return geometry;
    }

    void loadData(std::string inFile, bool loadGeometry = true, bool loadBinning = false) {
        bool isRAWfile = (strstr(inFile.c_str(),".real") != NULL);
        if (isRAWfile) {
            if (loadGeometry) {
                std::stringstream ss;
                ss << "file: " << inFile << " is real format and cannot contain geometry information";
                throw std::runtime_error(ss.str());
            } else if (loadBinning) {
                std::stringstream ss;
                ss << "file: " << inFile 
                   << " is real format and cannot contain binning information";
                throw std::runtime_error(ss.str());
            }
            
            projections =
                read_nd_array<float>(inFile.c_str());
        } else {
            // test header field
        	hid_t file_id = H5Fopen (inFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        	unsigned int dataformat_version;
					H5LTread_dataset (file_id, "/projection_dataformat_version", H5T_NATIVE_UINT, &dataformat_version);
					unsigned int needed = 1;
					if (!(dataformat_version == needed)) {
							std::stringstream ss;
							ss << "wrong projection data format version inside hdf5 file, found: "
								 << dataformat_version << ", needed: " << needed;
							throw std::runtime_error(ss.str());
					}

					hsize_t dim[3];
					H5LTget_dataset_info(file_id,"/projections",dim,NULL,NULL);
					std::vector<unsigned int> dims(dim,dim+3);
					projections = boost::shared_ptr<hoNDArray<float> >(new hoNDArray<float>(&dims));
					H5LTread_dataset (file_id,"/projections", H5T_NATIVE_FLOAT, projections->get_data_ptr());
					cleanup_projections = true;
					cleanup_geometry = true;

					if (loadGeometry)
							geometry->loadData(inFile);
					//if (loadBinning)
					//binning->loadData(inFile);
        }
    }
/*
    void saveData(std::string outFile) {
        createEmptyFile(outFile);
        appendVariable(outFile, "/projection_dataformat_version", 1);
        appendTypedArray3D(outFile, "/projections", *projections );        
        if (geometry != NULL)
            geometry->saveData(outFile);
    }
    */
    /*
    void setSAGx(float vx, float vy, float vz) {
        geometry->setSAGx(vx, vy, vz);
    }

    void setSAGy(float vx, float vy, float vz) {
        geometry->setSAGy(vx, vy, vz);
    }

    void setSpacing(float vx, float vy) {
        geometry->setSpacing(vx, vy);
    }

    void setAngles(float* a) {
        geometry->setAngles(a, projections->getNumberOfProjections());
    }
*/
}; // PS_Dataset
}
