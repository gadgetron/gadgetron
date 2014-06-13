#include "PhysioInterpolationGadget.h"
#include "Gadgetron.h"
#include "GadgetIsmrmrdReadWrite.h"
#include "GadgetronTimer.h"
#include "Spline.h"
#include "GtPlusDefinition.h"
#include "hoNDMetaAttributes.h"

#include <numeric>
#ifdef USE_OMP
    #include <omp.h>
#endif 

namespace Gadgetron{

    PhysioInterpolationGadget::PhysioInterpolationGadget() 
        : phys_time_index_(0)
        , phases_to_reconstruct_(30)
        , image_with_attrib_(false)
        , first_beat_on_trigger_(false)
    {
        set_parameter(std::string("physiology_time_index").c_str(), "0");
        set_parameter(std::string("mode").c_str(), "0");
        set_parameter(std::string("phases").c_str(), "30");
    }

    PhysioInterpolationGadget::~PhysioInterpolationGadget() {}

    int PhysioInterpolationGadget::process_config(ACE_Message_Block* mb)
    {
        phys_time_index_ = get_int_value("physiology_time_index");
        phases_to_reconstruct_ = get_int_value("phases");
        mode_ = get_int_value("mode");
        first_beat_on_trigger_ = get_bool_value("first_beat_on_trigger");

        boost::shared_ptr<ISMRMRD::ismrmrdHeader> cfg = Gadgetron::parseIsmrmrdXMLHeader(std::string(mb->rd_ptr()));

        ISMRMRD::ismrmrdHeader::encoding_sequence e_seq = cfg->encoding();
        ISMRMRD::encodingLimitsType e_limits = (*e_seq.begin()).encodingLimits();

        if (e_limits.slice().present())
        {
            slc_limit_ = e_limits.slice().get().maximum() + 1;
        }
        else
        {
            slc_limit_ = 1;
        }

        buffer_.resize(slc_limit_);

        size_t slc;
        for ( slc=0; slc<slc_limit_; slc++ )
        {
            buffer_[slc] = boost::shared_ptr< ACE_Message_Queue<ACE_MT_SYNCH> >(new ACE_Message_Queue<ACE_MT_SYNCH>(ACE_Message_Queue_Base::DEFAULT_HWM * 10, ACE_Message_Queue_Base::DEFAULT_LWM * 10) );
        }

        time_stamps_.resize(slc_limit_);

        return GADGET_OK;
    }

    int PhysioInterpolationGadget::close(unsigned long flags)
    {
        int ret = Gadget::close(flags);

        if ( flags != 0 )
        {
            GADGET_DEBUG1("PhysioInterpolationGadget::close...\n");

            size_t slc;
            for ( slc=0; slc<slc_limit_; slc ++ )
            {
                GADGET_DEBUG2("Processing slice: %d ... \n", slc);
                GADGET_DEBUG2("Number of items on Q: %d\n", buffer_[slc]->message_count());
                GADGET_DEBUG2("Image with attribute flag : %d\n", image_with_attrib_);

                if (time_stamps_[slc].size() != buffer_[slc]->message_count())
                {
                    GADGET_DEBUG1("Inconsistent number of messages and time stamps\n");
                    buffer_[slc]->flush();
                    return GADGET_FAIL;
                }

                float previous = -100.0;
                float sum_int  = 0.0; 
                std::vector<float> intervals;
                float int_count = 0.0;
                std::vector<size_t> cycle_starts;
                for (size_t i = 0; i < time_stamps_[slc].size(); i++)
                {
                    if ( (time_stamps_[slc][i] < previous) || (first_beat_on_trigger_ && i==0) )
                    {
                        cycle_starts.push_back(i);
                    }
                    else if (i > 0 )
                    {
                        sum_int += time_stamps_[slc][i]-time_stamps_[slc][i-1];
                        intervals.push_back(time_stamps_[slc][i]-time_stamps_[slc][i-1]);
                        int_count += 1.0;
                    }
                    previous = time_stamps_[slc][i];
                }

                std::sort(intervals.begin(),intervals.end());

                float mean_interval = sum_int/int_count;
                float median_interval = intervals[(intervals.size()>>1)];

                float average_cycle_length = 0.0f;
                std::vector<float> cycle_lengths;
                float count = 0;
                for (size_t i = 1; i < cycle_starts.size(); i++)
                {
                    float clength = time_stamps_[slc][cycle_starts[i]-1] + median_interval - time_stamps_[slc][cycle_starts[i]];
                    cycle_lengths.push_back(clength);
                }

                std::sort(cycle_lengths.begin(),cycle_lengths.end());
                float mean_cycle_length = std::accumulate(cycle_lengths.begin(), cycle_lengths.end(), 0.0f)/cycle_lengths.size();
                float median_cycle_length = cycle_lengths[(cycle_lengths.size()>>1)];

                GADGET_DEBUG2("We have %d full cyles, first one starting at %d\n", cycle_starts.size()-1, cycle_starts[0]);
                GADGET_DEBUG2("Mean/Median frame width %f/%f\n", mean_interval,median_interval);
                GADGET_DEBUG2("Mean/Median cycle_length %f/%f\n", mean_cycle_length,median_cycle_length);

                //Correct the first cycle assuming it is of median length:
                if ( !first_beat_on_trigger_ )
                {
                    float first_cycle_offset = (median_cycle_length-median_interval)+time_stamps_[slc][cycle_starts[0]]-time_stamps_[slc][cycle_starts[0]-1];
                    for (size_t i = 0; i < cycle_starts[0]; i++)
                    {
                        time_stamps_[slc][i] += first_cycle_offset;
                    }
                }

                //Calculate relative time stamps
                size_t current_cycle = 0;
                std::vector<float> relative_cycle_time;

                //Make sure we have cycle lengths for all the cycles we have covered
                cycle_lengths.insert(cycle_lengths.begin(),median_cycle_length);
                cycle_lengths.push_back(median_cycle_length);

                for (size_t i = 0; i < time_stamps_[slc].size(); i++)
                {
                    if ((current_cycle<cycle_starts.size()) && (i >= cycle_starts[current_cycle]) && (current_cycle < cycle_starts.size()))
                    {
                        current_cycle++;
                    }
                    relative_cycle_time.push_back(time_stamps_[slc][i]/cycle_lengths[current_cycle] + current_cycle);
                }

                //Make a temporary list of all the data pointers from the Q
                std::vector< ISMRMRD::ImageHeader* > hptrs;
                std::vector< hoNDArray< std::complex<float> > * > aptrs;
                std::vector< GtImageAttribType* > attribptrs;

                ACE_Message_Queue<ACE_MT_SYNCH>::ITERATOR it( *buffer_[slc] );
                for (ACE_Message_Block* entry = 0;
                    it.next (entry) != 0;
                    it.advance ()) 
                {
                    GadgetContainerMessage< ISMRMRD::ImageHeader >* tmpm1 =
                        AsContainerMessage< ISMRMRD::ImageHeader >(entry);

                    GadgetContainerMessage< hoNDArray< std::complex<float> > > * tmpm2 = 
                        AsContainerMessage< hoNDArray< std::complex<float> >  >(entry->cont());

                    if (!tmpm1 || !tmpm2 )
                    {
                        GADGET_DEBUG1("Failed to cast data on Q, bailing out\n");
                        buffer_[slc]->flush();
                        return GADGET_FAIL;
                    }

                    hptrs.push_back(tmpm1->getObjectPtr());
                    aptrs.push_back(tmpm2->getObjectPtr());

                    if ( image_with_attrib_ )
                    {
                        GadgetContainerMessage< GtImageAttribType > * tmpm3 = 
                            AsContainerMessage< GtImageAttribType >(entry->cont()->cont());

                        if ( !tmpm3 )
                        {
                            GADGET_DEBUG1("Failed to cast data on Q, bailing out\n");
                            buffer_[slc]->flush();
                            return GADGET_FAIL;
                        }

                        attribptrs.push_back(tmpm3->getObjectPtr());
                    }
                }

                //Let's figure out which time points we would like to interpolate on:
                ///TODO: Deal with mode 1 and other future modes, we are only implementing mode 0 at the moment
                float phase_interval = 1.0f/static_cast<float>(phases_to_reconstruct_);
                float max_time = floor(relative_cycle_time[relative_cycle_time.size()-1]);
                std::vector<float> recon_cycle_time;
                for (float t=1.0;t<(max_time-0.001);t+=phase_interval)
                {
                    recon_cycle_time.push_back(t);
                }

                if ( mode_ == 1 )
                {
                    std::vector<float> recon_cycle_time_first_beat(phases_to_reconstruct_);
                    memcpy(&recon_cycle_time_first_beat[0], &recon_cycle_time[0], sizeof(float)*phases_to_reconstruct_);
                    recon_cycle_time = recon_cycle_time_first_beat;
                }

                //Now we can loop over each pixel and estimate the new frames, but first we have to have somewhere to put the data
                std::vector< GadgetContainerMessage< ISMRMRD::ImageHeader >* > out_heads;
                std::vector< GadgetContainerMessage< hoNDArray< std::complex<float> > > * > out_data;
                std::vector< GadgetContainerMessage< GtImageAttribType> * > out_attrib;

                for (size_t i = 0; i < recon_cycle_time.size(); i++)
                {
                    GadgetContainerMessage<ISMRMRD::ImageHeader>* tmpm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>;
                    GadgetContainerMessage< hoNDArray< std::complex<float> > >* tmpm2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >;

                    tmpm1->cont(tmpm2);

                    (*tmpm1->getObjectPtr()) = (*hptrs[0]);
                    tmpm2->getObjectPtr()->create(aptrs[0]->get_dimensions());

                    out_heads.push_back(tmpm1);
                    out_data.push_back(tmpm2);

                    unsigned short current_cycle = static_cast<unsigned short>(floor(recon_cycle_time[i] + 0.0001));
                    unsigned short current_phase = static_cast<unsigned short>((recon_cycle_time[i]+0.0001-current_cycle)/(1.0/static_cast<float>(phases_to_reconstruct_)) + 0.0001);

                    tmpm1->getObjectPtr()->physiology_time_stamp[phys_time_index_] = static_cast<unsigned>(floor((recon_cycle_time[i]+0.0001-current_cycle)*cycle_lengths[current_cycle])); 
                    tmpm1->getObjectPtr()->phase = current_phase;
                    tmpm1->getObjectPtr()->image_index = current_phase+1 + (uint16_t)slc*phases_to_reconstruct_;
                    tmpm1->getObjectPtr()->image_series_index = current_cycle*10;

                    // make sure the phase is within the acquisition limit
                    if ( tmpm1->getObjectPtr()->phase+1 >= time_stamps_[slc].size() )
                    {
                        tmpm1->getObjectPtr()->phase = (uint16_t)(time_stamps_[slc].size()-1);
                    }

                    if ( image_with_attrib_ )
                    {
                        GadgetContainerMessage< GtImageAttribType >* tmpm3 = new GadgetContainerMessage< GtImageAttribType >;

                        tmpm2->cont(tmpm3);
                        (*tmpm3->getObjectPtr()) = (*attribptrs[0]);
                        out_attrib.push_back(tmpm3);

                        tmpm3->getObjectPtr()->attribute1_.set(GTPLUS_PHASE,      0, tmpm1->getObjectPtr()->phase);
                        tmpm3->getObjectPtr()->attribute1_.set(GTPLUS_IMAGENUMBER, 0, tmpm1->getObjectPtr()->image_index);

                        tmpm3->getObjectPtr()->attribute4_.set(GTPLUS_DATA_ROLE, "PhysioInterp");
                        tmpm3->getObjectPtr()->attribute4_.set(GTPLUS_IMAGECOMMENT, "PhysioInterp");
                        tmpm3->getObjectPtr()->attribute4_.set(GTPLUS_SEQUENCEDESCRIPTION, "_PhysioInterp");

                        tmpm3->getObjectPtr()->attribute4_.set(GTPLUS_IMAGEPROCESSINGHISTORY, "Interp");
                    }
                }

                //Let's interpolate the images
                size_t inelem = relative_cycle_time.size();
                size_t outelem = recon_cycle_time.size();
                size_t imageelem = aptrs[0]->get_number_of_elements();

                {
                    GadgetronTimer interptime("Interpolation Time");

        #ifdef USE_OMP
        #pragma omp parallel for
        #endif
                    for (long long p = 0; p < (long long)imageelem; p++)
                    {
                        std::vector< std::complex<float> > data_in(inelem);

                        //Get the input data for this pixel
                        for (size_t i = 0; i < inelem; i++) data_in[i] = aptrs[i]->get_data_ptr()[p];

                        //Interpolate the data
                        Spline<float, std::complex<float> > sp(relative_cycle_time, data_in);
                        std::vector<std::complex<float> > data_out = sp[recon_cycle_time];

                        //Copy it to the images
                        for (size_t i = 0; i < outelem; i++) out_data[i]->getObjectPtr()->get_data_ptr()[p] = data_out[i];
                    }
                }

                //Send out the images
                for (size_t i = 0; i < out_heads.size(); i++)
                {
                    if (this->next()->putq(out_heads[i]) < 0)
                    {
                        GADGET_DEBUG1("Unable to put data on next Gadgets Q\n");
                        return GADGET_FAIL;
                    }
                }

                //We can get rid of the buffered data now
                buffer_[slc]->flush();
            }
        }

        return ret;
    }

    int PhysioInterpolationGadget::
        process(GadgetContainerMessage< ISMRMRD::ImageHeader >* m1, GadgetContainerMessage< hoNDArray< std::complex<float> > >* m2)
    {
        GadgetContainerMessage<ISMRMRD::ImageHeader>* header = new GadgetContainerMessage<ISMRMRD::ImageHeader>;
        GadgetContainerMessage< hoNDArray< std::complex<float> > >* img = new GadgetContainerMessage< hoNDArray< std::complex<float> > >;

        (*header->getObjectPtr()) = (*m1->getObjectPtr());
        (*img->getObjectPtr()) = (*m2->getObjectPtr());
        header->cont(img);

        GadgetContainerMessage<GtImageAttribType>* m3 = 0;
        if (m2)
        {
            m3 = AsContainerMessage<GtImageAttribType>(m2->cont());
        }

        if ( m3 )
        {
            image_with_attrib_ = true;
        }
        else
        {
            image_with_attrib_ = false;
        }

        if ( image_with_attrib_ )
        {
            GadgetContainerMessage< GtImageAttribType >* attrib = new GadgetContainerMessage< GtImageAttribType >;
            (*attrib->getObjectPtr()) = *m3->getObjectPtr();
            img->cont(attrib);
        }

        uint16_t slc = header->getObjectPtr()->slice;

        if (buffer_[slc]->enqueue_tail(header) < 0)
        {
            GADGET_DEBUG1("Failed to add image to buffer\n");
            header->release();
            return GADGET_FAIL;
        }

        time_stamps_[slc].push_back( (float)(m1->getObjectPtr()->physiology_time_stamp[phys_time_index_]) );

        if (this->next()->putq(m1) < 0)
        {
            GADGET_DEBUG1("Unable to put data on next Gadgets Q\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(PhysioInterpolationGadget)
}
