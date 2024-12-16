#include "MaxwellCorrectionGadget.h"

namespace Gadgetron {

#ifdef M_PI
#undef M_PI
#endif // M_PI
#define M_PI 3.14159265358979323846

MaxwellCorrectionGadget::MaxwellCorrectionGadget(const Core::Context& context, const Core::GadgetProperties& props)
    : Core::ChannelGadget<mrd::Image<std::complex<float>>>(context, props) {
    maxwell_coefficients_present_ = false;
    maxwell_coefficients_ = std::vector<double>(4, 0);
    auto h = (context.header);

    if (h.user_parameters) {
        for (std::vector<mrd::UserParameterDoubleType>::const_iterator i(h.user_parameters->user_parameter_double.begin());
             i != h.user_parameters->user_parameter_double.end(); i++) {
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
        GDEBUG("MaxwellCorrection coefficients are supposed to be in the UserParameters. No user parameter section "
               "found\n");
        return;
    }

    patient_position_ = h.measurement_information->patient_position;

    for (size_t par_idx = 0; par_idx < h.user_parameters->user_parameter_long.size(); par_idx++) {

        std::string userParameterLong_tmp = h.user_parameters->user_parameter_long[par_idx].name;
        if (userParameterLong_tmp.compare("Flow_Dir") == 0) {
            this->FlowDirection_ = h.user_parameters->user_parameter_long[par_idx].value;
        }
    }

    maxwell_coefficients_present_ = true;

    GDEBUG("Maxwell Coefficients: %f, %f, %f, %f\n", maxwell_coefficients_[0], maxwell_coefficients_[1],
           maxwell_coefficients_[2], maxwell_coefficients_[3]);
}

void MaxwellCorrectionGadget::patient_to_physical_coordinate(std::vector<float>& norm_vec,
                                                             mrd::PatientPosition patient_position) {

    std::vector<float> tmp_nor_vec = norm_vec;

    if (patient_position == mrd::PatientPosition::kHFS) {
        norm_vec[0] = tmp_nor_vec[0];
        norm_vec[1] = -1 * tmp_nor_vec[1];
        norm_vec[2] = -1 * tmp_nor_vec[2];
    } else if (patient_position == mrd::PatientPosition::kFFS) {
        norm_vec[0] = -1 * tmp_nor_vec[0];
        norm_vec[1] = -1 * tmp_nor_vec[1];
        norm_vec[2] = tmp_nor_vec[2];
    } else if (patient_position == mrd::PatientPosition::kHFP) {
        norm_vec[0] = -1 * tmp_nor_vec[0];
        norm_vec[1] = tmp_nor_vec[1];
        norm_vec[2] = -1 * tmp_nor_vec[2];
    } else if (patient_position == mrd::PatientPosition::kFFP) {
        norm_vec[0] = tmp_nor_vec[0];
        norm_vec[1] = tmp_nor_vec[1];
        norm_vec[2] = tmp_nor_vec[2];
    } else if (patient_position == mrd::PatientPosition::kHFDL) {
        norm_vec[0] = -1 * tmp_nor_vec[1];
        norm_vec[1] = -1 * tmp_nor_vec[0];
        norm_vec[2] = -1 * tmp_nor_vec[2];
    } else if (patient_position == mrd::PatientPosition::kFFDL) {
        norm_vec[0] = tmp_nor_vec[1];
        norm_vec[1] = -1 * tmp_nor_vec[0];
        norm_vec[2] = tmp_nor_vec[2];
    } else if (patient_position == mrd::PatientPosition::kHFDR) {
        norm_vec[0] = tmp_nor_vec[1];
        norm_vec[1] = tmp_nor_vec[0];
        norm_vec[2] = -1 * tmp_nor_vec[2];
    } else if (patient_position == mrd::PatientPosition::kFFDR) {
        norm_vec[0] = -1 * tmp_nor_vec[1];
        norm_vec[1] = tmp_nor_vec[0];
        norm_vec[2] = tmp_nor_vec[2];
    }
}

void MaxwellCorrectionGadget::find_flow_dir(mrd::Image<std::complex<float>> image) {
    float flow_encoding_dir[3];

    if (this->FlowDirection_ == 1) {
        flow_encoding_dir[0] = image.head.line_dir[0];
        flow_encoding_dir[1] = image.head.line_dir[1];
        flow_encoding_dir[2] = image.head.line_dir[2];
    } else if (this->FlowDirection_ == 2) {
        flow_encoding_dir[0] = image.head.col_dir[0];
        flow_encoding_dir[1] = image.head.col_dir[1];
        flow_encoding_dir[2] = image.head.col_dir[2];
    } else {
        flow_encoding_dir[0] = image.head.col_dir[1] * image.head.line_dir[2] - image.head.col_dir[2] * image.head.line_dir[1];
        flow_encoding_dir[1] = image.head.col_dir[0] * image.head.line_dir[2] - image.head.col_dir[2] * image.head.line_dir[0];
        flow_encoding_dir[2] = image.head.col_dir[0] * image.head.line_dir[1] - image.head.col_dir[1] * image.head.line_dir[0];
    }

    float abs_flow_dir[3];
    abs_flow_dir[0] = std::abs(flow_encoding_dir[0]);
    abs_flow_dir[1] = std::abs(flow_encoding_dir[1]);
    abs_flow_dir[2] = std::abs(flow_encoding_dir[2]);

    int max_element_index = std::max_element(abs_flow_dir, abs_flow_dir + 3) - abs_flow_dir;
    this->FlipPhaseDirection_ = std::signbit(flow_encoding_dir[max_element_index]);
}

void MaxwellCorrectionGadget::process(Core::InputChannel<mrd::Image<std::complex<float>>>& in,
                                      Core::OutputChannel& out) {
    for (auto image : in) {
        if (maxwell_coefficients_present_) {

            int Nx = image.data.get_size(0);
            int Ny = image.data.get_size(1);
            int Nz = image.data.get_size(2);

            float dx = image.head.field_of_view[0] / Nx;
            float dy = image.head.field_of_view[1] / Ny;
            float dz = image.head.field_of_view[2] / Nz;

            this->RO_dir_Physical_.resize(3);
            this->PE_dir_Physical_.resize(3);
            this->SLC_dir_Physical_.resize(3);
            this->SLC_position_Physical_.resize(3);
            for (int idx = 0; idx < 3; idx++) {
                this->RO_dir_Physical_[idx] = -1 * image.head.col_dir[idx];
                this->PE_dir_Physical_[idx] = -1 * image.head.line_dir[idx];
                this->SLC_dir_Physical_[idx] = -1 * image.head.slice_dir[idx];
                this->SLC_position_Physical_[idx] = image.head.position[idx];
            }
            patient_to_physical_coordinate(this->RO_dir_Physical_, this->patient_position_);
            patient_to_physical_coordinate(this->PE_dir_Physical_, this->patient_position_);
            patient_to_physical_coordinate(this->SLC_dir_Physical_, this->patient_position_);
            patient_to_physical_coordinate(this->SLC_position_Physical_, this->patient_position_);

            std::vector<float> dR(3, 0);
            std::vector<float> dP(3, 0);
            std::vector<float> dS(3, 0);
            std::vector<float> p(3, 0);

            for (int z = 0; z < Nz; z++) {
                for (int y = 0; y < Ny; y++) {
                    for (int x = 0; x < Nx; x++) {

                        dR[0] = (x - Nx / 2 + 0.5) * dx * this->RO_dir_Physical_[0];
                        dR[1] = (x - Nx / 2 + 0.5) * dx * this->RO_dir_Physical_[1];
                        dR[2] = (x - Nx / 2 + 0.5) * dx * this->RO_dir_Physical_[2];

                        dP[0] = (y - Ny / 2 + 0.5) * dy * this->PE_dir_Physical_[0];
                        dP[1] = (y - Ny / 2 + 0.5) * dy * this->PE_dir_Physical_[1];
                        dP[2] = (y - Ny / 2 + 0.5) * dy * this->PE_dir_Physical_[2];

                        if (Nz > 1) {
                            dS[0] = (z - Nz / 2 + 0.5) * dz * this->SLC_dir_Physical_[0];
                            dS[1] = (z - Nz / 2 + 0.5) * dz * this->SLC_dir_Physical_[1];
                            dS[2] = (z - Nz / 2 + 0.5) * dz * this->SLC_dir_Physical_[2];
                        }

                        p[0] = this->SLC_position_Physical_[0] + dP[0] + dR[0] + dS[0];
                        p[1] = this->SLC_position_Physical_[1] + dP[1] + dR[1] + dS[1];
                        p[2] = this->SLC_position_Physical_[2] + dP[2] + dR[2] + dS[2];

                        // Convert to centimeters
                        p[0] = p[0] / 1000.0;
                        p[1] = p[1] / 1000.0;
                        p[2] = p[2] / 1000.0;

                        float delta_phi = maxwell_coefficients_[0] * p[2] * p[2] +
                                          maxwell_coefficients_[1] * (p[0] * p[0] + p[1] * p[1]) +
                                          maxwell_coefficients_[2] * p[0] * p[2] +
                                          maxwell_coefficients_[3] * p[1] * p[2];

                        long index = z * Ny * Nx + y * Nx + x;
                        std::complex<float>* data_ptr = image.data.data();
                        std::complex<float> correction = std::polar(1.0f, static_cast<float>(2 * M_PI * delta_phi));
                        data_ptr[index] *= correction;
                        // data_ptr[index] = correction;
                    }
                }
            }
        }
        out.push(std::move(image));
    }
}

GADGETRON_GADGET_EXPORT(MaxwellCorrectionGadget)

} // namespace Gadgetron