#include "GrappaWeightsCalculator.h"
#include "GadgetContainerMessage.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "hoNDArray_fileio.h"
#include "hoNDArray_reductions.h"
#include "GadgetronTimer.h"

#ifdef USE_CUDA
    #include "GPUTimer.h"
    #include "cuNDFFT.h"
    #include "b1_map.h"
    #include "htgrappa.h"
    #include <cuComplex.h>
#endif // USE_CUDA

#include "complext.h"

#include "hoNDFFT.h"
#include "hoNDArray_elemwise.h"
#include "mri_core_grappa.h"
#include "mri_core_coil_map_estimation.h"

namespace Gadgetron{

template <class T> class EXPORTGADGETSGRAPPA GrappaWeightsDescription
{

public:
    std::vector< std::pair<unsigned int, unsigned int> > sampled_region;
    unsigned int acceleration_factor;
    boost::shared_ptr<GrappaWeights<T> > destination;
    std::vector<unsigned int> uncombined_channel_weights;
    bool include_uncombined_channels_in_combined_weights;
};

template <class T> int GrappaWeightsCalculator<T>::svc(void)  {
    ACE_Message_Block *mb;

    while (this->getq(mb) >= 0) {
        if (mb->msg_type() == ACE_Message_Block::MB_HANGUP) {
            GDEBUG("Hanging up in weights calculator\n");
            if (this->putq(mb) == -1) {
              GERROR("GrappaWeightsCalculator::svc, putq");
              return -1;
            }
            break;
        }

        GadgetContainerMessage< GrappaWeightsDescription<T> >* mb1
        = AsContainerMessage< GrappaWeightsDescription<T> >(mb);

        if (!mb1) {
            mb->release();
            return -2;
        }

        GadgetContainerMessage< hoNDArray< std::complex<T> > >* mb2
        = AsContainerMessage< hoNDArray< std::complex<T> > >(mb1->cont());

        if (!mb2) {
            mb->release();
            return -3;
        }

        hoNDArray<float_complext>* host_data =
                reinterpret_cast< hoNDArray<float_complext>* >(mb2->getObjectPtr());

        size_t ks = 5;
        size_t power = 3;

#ifndef USE_CUDA
        use_gpu_ = false;
#endif // USE_CUDA

        if (use_gpu_)
        {
#ifdef USE_CUDA
            // Copy the image data to the device
            cuNDArray<float_complext> device_data(host_data);
            device_data.squeeze();

            std::vector<size_t> ftdims(2,0); ftdims[1] = 1;

            //Go to image space
             cuNDFFT<float>::instance()->ifft( &device_data, &ftdims);

            size_t RO = device_data.get_size(0);
            size_t E1 = device_data.get_size(1);
            size_t CHA = device_data.get_size(2);

            boost::shared_ptr< cuNDArray<float_complext> > csm;
            {
                //GPUTimer timer("GRAPPA CSM");
                csm = boost::make_shared<cuNDArray<float_complext>>(estimate_b1_map<float,2>( device_data, target_coils_ ));

                // estimate_b1_map_2D_NIH_Souheil( &device_data, &csm, ks, power, D, DH_D, V1, U1 );

                //GDEBUG("Coils in csm: %d\n", csm->get_size(2));
            }
            //Go back to kspace
            cuNDFFT<float>::instance()->fft(&device_data, &ftdims);

            cuNDArray<complext<float> > unmixing_dev;
            boost::shared_ptr< std::vector<size_t> > data_dimensions = device_data.get_dimensions();

            if (uncombined_channels_.size() > 0) {
                data_dimensions->push_back(uncombined_channels_.size()+1);
            }

            try{unmixing_dev.create(data_dimensions.get());}
            catch (std::runtime_error &err){
                GEXCEPTION(err,"Unable to allocate device memory for unmixing coeffcients\n");
                return GADGET_FAIL;
            }

            {
                //GPUTimer unmix_timer("GRAPPA Unmixing");
                //GadgetronTimer timer("GRAPPA unmixing", true);
                std::vector<unsigned int> kernel_size;

                //TODO: Add parameters for kernel size
                kernel_size.push_back(5);
                kernel_size.push_back(4);
                if ( htgrappa_calculate_grappa_unmixing(reinterpret_cast< cuNDArray<complext<float> >* >(&device_data),
                        csm.get(),
                        (unsigned int)(mb1->getObjectPtr()->acceleration_factor),
                        &kernel_size,
                        &unmixing_dev,
                        &(mb1->getObjectPtr()->sampled_region),
                        &uncombined_channels_) < 0) {
                    GDEBUG("GRAPPA unmixing coefficients calculation failed\n");
                    return GADGET_FAIL;
                }
            }

            if (mb1->getObjectPtr()->destination) {
                boost::shared_ptr< hoNDArray<complext<float> > > unmixing_host = unmixing_dev.to_host();

                //TODO: This reshaping needs to take uncombined channels into account
                boost::shared_ptr< std::vector<size_t> > tmp_dims = mb2->getObjectPtr()->get_dimensions();
                if (uncombined_channels_.size()) tmp_dims->push_back((size_t)(uncombined_channels_.size() + 1));

                try {
                    unmixing_host->reshape(tmp_dims.get());
                }
                catch (std::runtime_error &err){
                    GEXCEPTION(err, "Reshaping of GRAPPA weights failed \n");

                }

                if (mb1->getObjectPtr()->destination->update(reinterpret_cast<hoNDArray<std::complex<float> >* >(unmixing_host.get())) < 0) {
                    GDEBUG("Update of GRAPPA weights failed\n");
                    return GADGET_FAIL;
                }
            }
            else {
                GDEBUG("Undefined GRAPPA weights destination\n");
                return GADGET_FAIL;
            }
#endif // USE_CUDA
        }
        else
        {
            host_data->squeeze();

            size_t RO = host_data->get_size(0);
            size_t E1 = host_data->get_size(1);
            size_t CHA = host_data->get_size(2);

            std::vector<size_t> data_dimensions;
            host_data->get_dimensions(data_dimensions);

            if (uncombined_channels_.size() > 0) {
                data_dimensions.push_back(uncombined_channels_.size() + 1);
            }

            try{ unmixing_.create(data_dimensions); }
            catch (std::runtime_error &err){
                GEXCEPTION(err, "Unable to allocate host memory for unmixing coeffcients\n");
                return GADGET_FAIL;
            }

            // compute the unmixing coefficients
            size_t numUnCombined = uncombined_channels_.size();

            double thres = 0.0005;
            size_t kRO = 5;
            size_t kNE1 = 4;

            if (numUnCombined==0)
            {
                hoNDArray< std::complex<float> > acs(RO, E1, target_coils_, reinterpret_cast< std::complex<float>* >(host_data->begin()));
                hoNDArray< std::complex<float> > target_acs(RO, E1, target_coils_, acs.begin());

                // estimate coil map
                if (!complex_im_.dimensions_equal(&target_acs))
                {
                    complex_im_.create(RO, E1, target_coils_);
                }

                hoNDFFT<float>::instance()->ifft2c(target_acs, complex_im_);
                Gadgetron::coil_map_2d_Inati(complex_im_, coil_map_, ks, power);

                // compute unmixing coefficients
                if (mb1->getObjectPtr()->acceleration_factor == 1)
                {
                    Gadgetron::conjugate(coil_map_, coil_map_);
                    Gadgetron::clear(unmixing_);
                    memcpy(unmixing_.begin(), coil_map_.begin(), coil_map_.get_number_of_bytes());
                }
                else
                {
                    size_t startRO = mb1->getObjectPtr()->sampled_region[0].first;
                    size_t endRO = mb1->getObjectPtr()->sampled_region[0].second;

                    size_t startE1 = mb1->getObjectPtr()->sampled_region[1].first;
                    size_t endE1 = mb1->getObjectPtr()->sampled_region[1].second;

                    Gadgetron::grappa2d_calib_convolution_kernel(acs, target_acs,
                        (size_t)(mb1->getObjectPtr()->acceleration_factor),
                        thres, kRO, kNE1, startRO, endRO, startE1, endE1, conv_ker_);

                    Gadgetron::grappa2d_image_domain_kernel(conv_ker_, RO, E1, kIm_);

                    Gadgetron::clear(unmixing_);

                    Gadgetron::grappa2d_unmixing_coeff(kIm_, coil_map_, (size_t)(mb1->getObjectPtr()->acceleration_factor), unmixing_, gFactor_);

                    // GDEBUG_STREAM("cpu triggered - unmixing_ : " << Gadgetron::norm2(unmixing_));
                }
            }
            else
            {
                hoNDArray< std::complex<float> > acs(RO, E1, CHA, reinterpret_cast< std::complex<float>* >(host_data->begin()));

                // handle the case that all channels are reconed
                size_t target_coils_with_uncombined = target_coils_ + numUnCombined;
                if (target_coils_with_uncombined > CHA) target_coils_with_uncombined = CHA;

                std::list<unsigned int>::iterator it;

                // compute unmixing coefficients
                if (mb1->getObjectPtr()->acceleration_factor == 1)
                {
                    // if no acceleration, the input data is used for coil map estimation
                    // no need to differentiate combined and uncombined channels
                    if (!complex_im_.dimensions_equal(&acs))
                    {
                        complex_im_.create(RO, E1, acs.get_size(2));
                    }

                    hoNDFFT<float>::instance()->ifft2c(acs, complex_im_);

                    Gadgetron::coil_map_2d_Inati(complex_im_, coil_map_, ks, power);

                    Gadgetron::conjugate(coil_map_, coil_map_);

                    Gadgetron::clear(unmixing_);

                    // copy back to unmixing
                    memcpy(unmixing_.begin(), coil_map_.begin(), sizeof(std::complex<float>)*RO*E1*CHA);

                    // set uncombined channels
                    size_t ind = 1;
                    for (it = uncombined_channels_.begin(); it != uncombined_channels_.end(); it++)
                    {
                        std::complex<float>* pUnmixing = unmixing_.begin() + ind*RO*E1*CHA + (*it)*RO*E1;
                        for (size_t p = 0; p<RO*E1; p++)
                        {
                            pUnmixing[p] = 1;
                        }

                        ind++;
                    }
                }
                else
                {
                    // first, assemble the target_acs
                    // the combined channel comes first and then all uncombined channels
                    std::vector<size_t> dimTarget(3);
                    dimTarget[0] = RO;
                    dimTarget[1] = E1;
                    dimTarget[2] = target_coils_with_uncombined;

                    if (!target_acs_.dimensions_equal(&dimTarget))
                    {
                        target_acs_.create(RO, E1, target_coils_with_uncombined);
                    }

                    // copy first combined channels and all uncombined channels to target_acs_
                    size_t sCha, ind(0), ind_uncombined(0);

                    // record from which src channel, a target channel is selected
                    std::vector<size_t> srcChaLoc(target_coils_with_uncombined);
                    for (sCha = 0; sCha<CHA; sCha++)
                    {
                        bool uncombined = false;
                        for (it = uncombined_channels_.begin(); it != uncombined_channels_.end(); it++)
                        {
                            if (sCha == *it)
                            {
                                uncombined = true;
                                break;
                            }
                        }

                        if (!uncombined)
                        {
                            if ( ind < (target_coils_with_uncombined - numUnCombined) )
                            {
                                memcpy(target_acs_.begin() + ind * RO*E1, acs.begin() + sCha * RO*E1, sizeof(std::complex<float>)*RO*E1);
                                srcChaLoc[ind] = sCha;
                                ind++;
                            }
                        }
                        else
                        {
                            memcpy(target_acs_.begin() + (target_coils_with_uncombined - numUnCombined + ind_uncombined) * RO*E1, acs.begin() + sCha * RO*E1, sizeof(std::complex<float>)*RO*E1);
                            srcChaLoc[target_coils_with_uncombined - numUnCombined + ind_uncombined] = sCha;
                            ind_uncombined++;
                        }
                    }

                    if (!complex_im_.dimensions_equal(&target_acs_))
                    {
                        complex_im_.create(RO, E1, target_acs_.get_size(2));
                    }

                    hoNDFFT<float>::instance()->ifft2c(target_acs_, complex_im_);

                    Gadgetron::coil_map_2d_Inati(complex_im_, coil_map_, ks, power);

                    Gadgetron::grappa2d_calib_convolution_kernel(acs, target_acs_,
                        (size_t)(mb1->getObjectPtr()->acceleration_factor),
                        thres, kRO, kNE1, conv_ker_);

                    Gadgetron::grappa2d_image_domain_kernel(conv_ker_, RO, E1, kIm_);

                    // kIm_ stored the unwrapping coefficients as [RO E1 CHA target_coils_with_uncombined]
                    // for the target_coils_with_uncombined dimension, combined channels come first and then uncombined channels

                    Gadgetron::clear(unmixing_);

                    hoNDArray< std::complex<float> > unmixing_all_channels(RO, E1, CHA, unmixing_.begin());
                    Gadgetron::grappa2d_unmixing_coeff(kIm_, coil_map_, (size_t)(mb1->getObjectPtr()->acceleration_factor), unmixing_all_channels, gFactor_);

                    // set unmixing coefficients for uncombined channels
                    ind = 1;
                    for (it = uncombined_channels_.begin(); it != uncombined_channels_.end(); it++)
                    {
                        memcpy(unmixing_.begin() + ind*RO*E1*CHA, kIm_.begin() + (target_coils_with_uncombined - numUnCombined + ind - 1)*RO*E1*CHA, sizeof(std::complex<float>)*RO*E1*CHA);
                        ind++;
                    }
                }
            }

            // pass the unmixing coefficients
            if (mb1->getObjectPtr()->destination)
            {
                boost::shared_ptr< hoNDArray< std::complex<float> > > unmixing_host(new hoNDArray< std::complex<float> >());
                boost::shared_ptr< std::vector<size_t> > tmp_dims = mb2->getObjectPtr()->get_dimensions();
                if (uncombined_channels_.size()) tmp_dims->push_back((size_t)(uncombined_channels_.size() + 1));

                try {
                    unmixing_host->create(tmp_dims.get());
                    Gadgetron::clear(*unmixing_host);
                }
                catch (std::runtime_error &err){
                    GEXCEPTION(err, "Reshaping of GRAPPA weights failed \n");

                }

                memcpy(unmixing_host->begin(), unmixing_.begin(), unmixing_.get_number_of_bytes());

                // GDEBUG_STREAM("cpu triggered ... : " << Gadgetron::norm2(*unmixing_host));

                if (mb1->getObjectPtr()->destination->update(unmixing_host.get()) < 0) {
                    GDEBUG("Update of GRAPPA weights failed\n");
                    return GADGET_FAIL;
                }
            }
            else {
                GDEBUG("Undefined GRAPPA weights destination\n");
                return GADGET_FAIL;
            }
        }

        mb->release();
    }

    return 0;
}

template <class T> int GrappaWeightsCalculator<T>::close(unsigned long flags) {
    int rval = 0;
    if (flags == 1) {
        ACE_Message_Block *hangup = new ACE_Message_Block();
        hangup->msg_type( ACE_Message_Block::MB_HANGUP );
        if (this->putq(hangup) == -1) {
            hangup->release();
            GERROR("GrappaWeightsCalculator::close, putq");
            return -1;
        }
        //GDEBUG("Waiting for weights calculator to finish\n");
        rval = this->wait();
        //GDEBUG("Weights calculator to finished\n");
    }
    return rval;
}


template <class T> int GrappaWeightsCalculator<T>::
add_job( hoNDArray< std::complex<T> >* ref_data,
        std::vector< std::pair<unsigned int, unsigned int> > sampled_region,
        unsigned int acceleration_factor,
        boost::shared_ptr< GrappaWeights<T> > destination,
        std::vector<unsigned int> uncombined_channel_weights,
        bool include_uncombined_channels_in_combined_weights)
        {

    GadgetContainerMessage< GrappaWeightsDescription<T> >* mb1 =
            new GadgetContainerMessage< GrappaWeightsDescription<T> >();

    if (!mb1) {
        return -1;
    }

    /*
  for (unsigned int i = 0; i < sampled_region.size(); i++) {
      GDEBUG("Sampled region %d: [%d, %d]\n", i, sampled_region[i].first, sampled_region[i].second);
  }
     */

    mb1->getObjectPtr()->sampled_region = sampled_region;
    mb1->getObjectPtr()->acceleration_factor = acceleration_factor;
    mb1->getObjectPtr()->destination = destination;
    mb1->getObjectPtr()->uncombined_channel_weights = uncombined_channel_weights;
    mb1->getObjectPtr()->include_uncombined_channels_in_combined_weights =
            include_uncombined_channels_in_combined_weights;


    GadgetContainerMessage< hoNDArray< std::complex<T> > >* mb2 =
            new GadgetContainerMessage< hoNDArray< std::complex<T> > >();

    if (!mb2) {
        mb1->release();
        return -2;
    }

    mb1->cont(mb2);

    try{mb2->getObjectPtr()->create(ref_data->get_dimensions().get());}
    catch (std::runtime_error &err ){
        mb1->release();
        return -3;
    }

    memcpy(mb2->getObjectPtr()->get_data_ptr(), ref_data->get_data_ptr(),
            ref_data->get_number_of_elements()*sizeof(T)*2);

    this->putq(mb1);

    return 0;
        }

template <class T> int GrappaWeightsCalculator<T>::add_uncombined_channel(unsigned int channel_id)
        {
    remove_uncombined_channel(channel_id);
    uncombined_channels_.push_back(channel_id);
    return 0;
        }

template <class T> int GrappaWeightsCalculator<T>::remove_uncombined_channel(unsigned int channel_id)
        {
    uncombined_channels_.remove(channel_id);
    return 0;
        }



template class EXPORTGADGETSGRAPPA GrappaWeightsDescription<float>;
template class EXPORTGADGETSGRAPPA GrappaWeightsCalculator<float>;
//template class EXPORTGADGETSGRAPPA GrappaWeightsCalculator<double>; //TOFO
//template class EXPORTGADGETSGRAPPA GrappaWeightsDescription<double>;

}
