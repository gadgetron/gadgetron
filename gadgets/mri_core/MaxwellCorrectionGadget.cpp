#include "MaxwellCorrectionGadget.h"
#include "GadgetronTimer.h"
#include "Spline.h"
#include "ismrmrd/xml.h"

#include <numeric>
#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

#ifdef M_PI
#undef M_PI
#endif // M_PI
#define M_PI 3.14159265358979323846

    MaxwellCorrectionGadget::MaxwellCorrectionGadget()
        : maxwell_coefficients_present_(false)
        , maxwell_coefficients_(4,0)
    {
    }

    MaxwellCorrectionGadget::~MaxwellCorrectionGadget() {}

    int MaxwellCorrectionGadget::process_config(ACE_Message_Block* mb)
    {

        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(),h);

        if (h.userParameters)
        {
            for (std::vector<ISMRMRD::UserParameterDouble>::const_iterator i (h.userParameters->userParameterDouble.begin()); 
                i != h.userParameters->userParameterDouble.end(); i++)
            {
                    if (i->name == "MaxwellCoefficient_0") {
                        maxwell_coefficients_[0] = i->value;
                    } else if (i->name == "MaxwellCoefficient_1") {
                        maxwell_coefficients_[1] = i->value;
                    } else if (i->name == "MaxwellCoefficient_2") {
                        maxwell_coefficients_[2] = i->value;
                    } else if (i->name == "MaxwellCoefficient_3") {
                        maxwell_coefficients_[3] = i->value;
                    } else {
                        GDEBUG("WARNING: unused user parameter parameter %s found\n", i->name.c_str());
                    }
            }
        } else {
            GDEBUG("MaxwellCorrection coefficients are supposed to be in the UserParameters. No user parameter section found\n");
            return GADGET_OK;
        }

        patient_position_ = h.measurementInformation->patientPosition;


        for (size_t par_idx = 0; par_idx < h.userParameters->userParameterLong.size(); par_idx++){

            std::string userParameterLong_tmp = h.userParameters->userParameterLong[par_idx].name;
            if ( userParameterLong_tmp.compare("Flow_Dir") == 0)
            {
            	this->FlowDirection_ = h.userParameters->userParameterLong[par_idx].value;
            }

        }

        maxwell_coefficients_present_ = true;

        GDEBUG("Maxwell Coefficients: %f, %f, %f, %f\n", maxwell_coefficients_[0], maxwell_coefficients_[1], maxwell_coefficients_[2], maxwell_coefficients_[3]);

        return GADGET_OK;
    }

    void MaxwellCorrectionGadget::patient_to_physical_coordinate(std::vector<float> &norm_vec, std::string patient_position){

    	std::vector<float> tmp_nor_vec = norm_vec;

     if (patient_position == "HFS"){
    	 norm_vec[0] = tmp_nor_vec[0];
    	 norm_vec[1] = -1*tmp_nor_vec[1];
    	 norm_vec[2] = -1*tmp_nor_vec[2];
     }else if (patient_position == "FFS"){
    	 norm_vec[0] = -1*tmp_nor_vec[0];
    	 norm_vec[1] = -1*tmp_nor_vec[1];
    	 norm_vec[2] = tmp_nor_vec[2];
     }else if (patient_position ==  "HFP"){
    	 norm_vec[0] = -1*tmp_nor_vec[0];
    	 norm_vec[1] = tmp_nor_vec[1];
    	 norm_vec[2] = -1*tmp_nor_vec[2];
     }else if (patient_position ==  "FFP"){
    	 norm_vec[0] = tmp_nor_vec[0];
    	 norm_vec[1] = tmp_nor_vec[1];
    	 norm_vec[2] = tmp_nor_vec[2];
     }else if (patient_position ==  "HFDL"){
    	 norm_vec[0] = -1*tmp_nor_vec[1];
    	 norm_vec[1] = -1*tmp_nor_vec[0];
    	 norm_vec[2] = -1*tmp_nor_vec[2];
     }else if (patient_position ==  "FFDL"){
    	 norm_vec[0] = tmp_nor_vec[1];
    	 norm_vec[1] = -1*tmp_nor_vec[0];
    	 norm_vec[2] = tmp_nor_vec[2];
     }else if(patient_position ==  "HFDR"){
    	 norm_vec[0] = tmp_nor_vec[1];
    	 norm_vec[1] = tmp_nor_vec[0];
    	 norm_vec[2] = -1*tmp_nor_vec[2];
     }else if (patient_position ==  "FFDR"){
    	 norm_vec[0] = -1*tmp_nor_vec[1];
    	 norm_vec[1] = tmp_nor_vec[0];
    	 norm_vec[2] = tmp_nor_vec[2];
     }


    }


    void MaxwellCorrectionGadget::
	find_flow_dir(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1){

        float flow_encoding_dir[3];

        if (this->FlowDirection_ == 1){
        	flow_encoding_dir[0] = m1->getObjectPtr()->phase_dir[0];
        	flow_encoding_dir[1] = m1->getObjectPtr()->phase_dir[1];
        	flow_encoding_dir[2] = m1->getObjectPtr()->phase_dir[2];
        }else if(this->FlowDirection_ == 2){
        	flow_encoding_dir[0] = m1->getObjectPtr()->read_dir[0];
        	flow_encoding_dir[1] = m1->getObjectPtr()->read_dir[1];
        	flow_encoding_dir[2] = m1->getObjectPtr()->read_dir[2];
        }else{
        	flow_encoding_dir[0] = m1->getObjectPtr()->read_dir[1]*m1->getObjectPtr()->phase_dir[2] - m1->getObjectPtr()->read_dir[2]*m1->getObjectPtr()->phase_dir[1];
        	flow_encoding_dir[1] = m1->getObjectPtr()->read_dir[0]*m1->getObjectPtr()->phase_dir[2] - m1->getObjectPtr()->read_dir[2]*m1->getObjectPtr()->phase_dir[0];
        	flow_encoding_dir[2] = m1->getObjectPtr()->read_dir[0]*m1->getObjectPtr()->phase_dir[1] - m1->getObjectPtr()->read_dir[1]*m1->getObjectPtr()->phase_dir[0];
        }

        float abs_flow_dir[3];
        abs_flow_dir[0] = std::abs(flow_encoding_dir[0]);
        abs_flow_dir[1] = std::abs(flow_encoding_dir[1]);
        abs_flow_dir[2] = std::abs(flow_encoding_dir[2]);

        int max_element_index = std::max_element(abs_flow_dir, abs_flow_dir+3) - abs_flow_dir;
        this->FlipPhaseDirection_ = std::signbit(flow_encoding_dir[max_element_index]);

    }


    int MaxwellCorrectionGadget::
        process(GadgetContainerMessage<ISMRMRD::ImageHeader>* m1,
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {
        if (maxwell_coefficients_present_) {
            //GDEBUG("Got coefficients\n");

            int Nx = m2->getObjectPtr()->get_size(0);
            int Ny = m2->getObjectPtr()->get_size(1);
            int Nz = m2->getObjectPtr()->get_size(2);

            float dx = m1->getObjectPtr()->field_of_view[0] / Nx;
            float dy = m1->getObjectPtr()->field_of_view[1] / Ny;
            float dz = m1->getObjectPtr()->field_of_view[2] / Nz;

            this->RO_dir_Physical_.resize(3);
            this->PE_dir_Physical_.resize(3);
            this->SLC_dir_Physical_.resize(3);
            this->SLC_position_Physical_.resize(3);
            for (int idx = 0; idx < 3; idx++){
            	this->RO_dir_Physical_[idx] = -1*m1->getObjectPtr()->read_dir[idx];
            	this->PE_dir_Physical_[idx] = -1*m1->getObjectPtr()->phase_dir[idx];
            	this->SLC_dir_Physical_[idx] = -1*m1->getObjectPtr()->slice_dir[idx];
            	this->SLC_position_Physical_[idx] = m1->getObjectPtr()->position[idx];
            }
            patient_to_physical_coordinate(this->RO_dir_Physical_, this->patient_position_);
            patient_to_physical_coordinate(this->PE_dir_Physical_, this->patient_position_);
            patient_to_physical_coordinate(this->SLC_dir_Physical_, this->patient_position_);
            patient_to_physical_coordinate(this->SLC_position_Physical_, this->patient_position_);

            //this->find_flow_dir(m1);


            /*
            GDEBUG("Nx = %d, Ny = %d, Nz = %d\n", Nx, Ny, Nz);
            GDEBUG("dx = %f, dy = %f, dz = %f\n", dx, dy, dz);
            GDEBUG("img_pos_x = %f, img_pos_y = %f, img_pos_z = %f\n", m1->getObjectPtr()->position[0], m1->getObjectPtr()->position[1], m1->getObjectPtr()->position[2]);
            */

            std::vector<float> dR(3,0);
            std::vector<float> dP(3,0);
            std::vector<float> dS(3,0);
            std::vector<float> p(3,0);

            for (int z = 0; z < Nz; z++) {
                for (int y = 0; y < Ny; y++) {
                    for (int x = 0; x < Nx; x++) {

                        dR[0] = (x-Nx/2+0.5) * dx * this->RO_dir_Physical_[0];
                        dR[1] = (x-Nx/2+0.5) * dx * this->RO_dir_Physical_[1];
                        dR[2] = (x-Nx/2+0.5) * dx * this->RO_dir_Physical_[2];

                        dP[0] = (y-Ny/2+0.5) * dy * this->PE_dir_Physical_[0];
                        dP[1] = (y-Ny/2+0.5) * dy * this->PE_dir_Physical_[1];
                        dP[2] = (y-Ny/2+0.5) * dy * this->PE_dir_Physical_[2];

                        if (Nz > 1) {
                            dS[0] = (z-Nz/2+0.5) * dz * this->SLC_dir_Physical_[0];
                            dS[1] = (z-Nz/2+0.5) * dz * this->SLC_dir_Physical_[1];
                            dS[2] = (z-Nz/2+0.5) * dz * this->SLC_dir_Physical_[2];
                        }

                        p[0] = this->SLC_position_Physical_[0] + dP[0] + dR[0] + dS[0];
                        p[1] = this->SLC_position_Physical_[1] + dP[1] + dR[1] + dS[1];
                        p[2] = this->SLC_position_Physical_[2] + dP[2] + dR[2] + dS[2];

                        //Convert to centimeters
                        p[0] = p[0]/1000.0;
                        p[1] = p[1]/1000.0;
                        p[2] = p[2]/1000.0;

                        float delta_phi = maxwell_coefficients_[0]*p[2]*p[2] +
                            maxwell_coefficients_[1]*(p[0]*p[0] + p[1]*p[1]) + 
                            maxwell_coefficients_[2]*p[0]*p[2] + 
                            maxwell_coefficients_[3]*p[1]*p[2];

                        //if (this->FlipPhaseDirection_) delta_phi *= -1.0;

                        long index = z*Ny*Nx+y*Nx+x;
                        std::complex<float>* data_ptr = m2->getObjectPtr()->get_data_ptr();

                        std::complex<float> correction = std::polar(1.0f,static_cast<float>(2*M_PI*delta_phi));

                        data_ptr[index] *= correction;
                        //data_ptr[index] = correction;
                    }
                }
            }

        }

        if (this->next()->putq(m1) < 0) {
            GDEBUG("Unable to put data on next Gadgets Q\n");
            return GADGET_FAIL;
        }
        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(MaxwellCorrectionGadget)
}
