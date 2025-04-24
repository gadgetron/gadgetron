
#include "GenericMRDCommentGadget.h"
#include "mri_core_def.h"
#include "mri_core_utility.h"

namespace Gadgetron {


template <typename T > 
GenericMRDCommentGadget<T>::GenericMRDCommentGadget(const Core::Context& context, const Core::GadgetProperties& props) : Core::ChannelGadget<Core::Image<T>>(context, props)
{
    // initialize the dict
    dict_gt_2_mrd_[GADGETRON_DATA_ROLE] = "DataRole";
    dict_gt_2_mrd_[GADGETRON_SEQUENCEDESCRIPTION] = "SequenceDescriptionAdditional";
    dict_gt_2_mrd_[GADGETRON_IMAGECOMMENT] = "ImageComments";
    dict_gt_2_mrd_[GADGETRON_IMAGE_SCALE_OFFSET] = "RescaleIntercept";
    dict_gt_2_mrd_[GADGETRON_IMAGE_SCALE_RATIO] = "RescaleSlope";
    dict_gt_2_mrd_[GADGETRON_IMAGE_WINDOWCENTER] = "WindowCenter";
    dict_gt_2_mrd_[GADGETRON_IMAGE_WINDOWWIDTH] = "WindowWidth";
    dict_gt_2_mrd_[GADGETRON_IMAGE_COLORMAP] = "LUTFileName";
    dict_gt_2_mrd_[GADGETRON_IMAGE_ECHOTIME] = "EchoTime";
    dict_gt_2_mrd_[GADGETRON_IMAGE_INVERSIONTIME] = "InversionTime";
    dict_gt_2_mrd_[GADGETRON_IMAGE_SATURATIONTIME] = "SaturationTime";
    dict_gt_2_mrd_[GADGETRON_2D_ROI] = "ROI";
    dict_gt_2_mrd_[GADGETRON_IMAGEPROCESSINGHISTORY] = "ImageTypeValue4";
    dict_gt_2_mrd_[GADGETRON_DIRECT_IMAGE_SEND] = "DirectSend";
}

template <typename T > 
void GenericMRDCommentGadget<T>::process(Core::InputChannel<Core::Image<T>>& in, Core::OutputChannel& out)
{
    std::vector<std::string> simple_str_fields;
    // simple_str_fields.push_back(GADGETRON_DATA_ROLE); // handle this seperately
    simple_str_fields.push_back(GADGETRON_SEQUENCEDESCRIPTION);
    simple_str_fields.push_back(GADGETRON_IMAGECOMMENT);
    simple_str_fields.push_back(GADGETRON_IMAGE_COLORMAP);
    simple_str_fields.push_back(GADGETRON_IMAGEPROCESSINGHISTORY);

    std::vector<std::string> simple_long_fields;
    simple_long_fields.push_back(GADGETRON_IMAGE_WINDOWCENTER);
    simple_long_fields.push_back(GADGETRON_IMAGE_WINDOWWIDTH);
    simple_long_fields.push_back(GADGETRON_DIRECT_IMAGE_SEND);
    //simple_long_fields.push_back(GADGETRON_IMAGE_SCALE_OFFSET);
    //simple_long_fields.push_back(GADGETRON_IMAGE_SCALE_RATIO);

    std::vector<std::string> simple_double_fields;
    simple_double_fields.push_back(GADGETRON_IMAGE_ECHOTIME);
    simple_double_fields.push_back(GADGETRON_IMAGE_INVERSIONTIME);
    simple_double_fields.push_back(GADGETRON_IMAGE_SATURATIONTIME);

    for (auto [header,data,meta] : in)
    {
        if (meta)
        {
            for (auto id : simple_str_fields)
            {
                if (meta->length(id.c_str()) > 0)
                {
                    std::vector<std::string> v;
                    Gadgetron::get_ismrmrd_meta_values(*meta, id, v);
                    Gadgetron::set_ismrmrd_meta_values(*meta, dict_gt_2_mrd_[id], v);
                }
            }

            for (auto id : simple_long_fields)
            {
                if (meta->length(id.c_str()) > 0)
                {
                    std::vector<long> v;
                    Gadgetron::get_ismrmrd_meta_values(*meta, id, v);
                    Gadgetron::set_ismrmrd_meta_values(*meta, dict_gt_2_mrd_[id], v);
                }
            }

            for (auto id : simple_double_fields)
            {
                if (meta->length(id.c_str()) > 0)
                {
                    std::vector<double> v;
                    Gadgetron::get_ismrmrd_meta_values(*meta, id, v);
                    Gadgetron::set_ismrmrd_meta_values(*meta, dict_gt_2_mrd_[id], v);
                }
            }

            // handle ROI here
            typedef std::map<std::string, std::vector<ISMRMRD::MetaValue> > map_t;

            std::vector<std::string> rois;

            map_t::iterator iter;
            for (iter=meta->begin(); iter!=meta->end(); iter++)
            {
                std::string key = (*iter).first;
                if(key.find(GADGETRON_2D_ROI)!=std::string::npos)
                {
                    rois.push_back(key);
                }
            }

            for (auto id : rois)
            {
                std::vector<double> v;
                Gadgetron::get_ismrmrd_meta_values(*meta, id, v);

                std::string new_id(id);
                new_id.replace(0, std::string(GADGETRON_2D_ROI).length(), dict_gt_2_mrd_[GADGETRON_2D_ROI]);

                auto it = v.begin() + 4; // line style
                v.insert(it, 0);

                auto it2 = v.begin() + 5; // visibility
                v.insert(it2, 1);

                Gadgetron::set_ismrmrd_meta_values(*meta, new_id, v);
            }

            // handle data role
            std::vector<std::string> v;
            Gadgetron::get_ismrmrd_meta_values(*meta, GADGETRON_DATA_ROLE, v);

            bool set_data_role = true;
            for (auto id : v)
            {
                if (id.compare(GADGETRON_IMAGE_REGULAR)==0)
                {
                    set_data_role = false;
                    break;
                }

                if (id.compare(GADGETRON_IMAGE_RETRO)==0)
                {
                    set_data_role = false;
                    break;
                }

                if (id.compare(GADGETRON_IMAGE_MOCORECON)==0)
                {
                    set_data_role = false;
                    break;
                }
            }

            if (set_data_role)
            {
                Gadgetron::set_ismrmrd_meta_values(*meta, dict_gt_2_mrd_[GADGETRON_DATA_ROLE], v);
            }
        }

        out.push(std::move(header), std::move(data), std::move(meta));
    }
}

GADGETRON_GADGET_EXPORT(GenericMRDCommentShortGadget)
GADGETRON_GADGET_EXPORT(GenericMRDCommentUShortGadget)
GADGETRON_GADGET_EXPORT(GenericMRDCommentIntGadget)
GADGETRON_GADGET_EXPORT(GenericMRDCommentUIntGadget)
GADGETRON_GADGET_EXPORT(GenericMRDCommentFloatGadget)
GADGETRON_GADGET_EXPORT(GenericMRDCommentCxFloatGadget)
GADGETRON_GADGET_EXPORT(GenericMRDCommentCxDoubleGadget)

} // namespace Gadgetron