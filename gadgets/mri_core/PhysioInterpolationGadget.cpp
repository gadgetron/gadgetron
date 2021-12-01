#include "PhysioInterpolationGadget.h"
#include "GadgetronTimer.h"
#include "Spline.h"
#include "mri_core_def.h"
#include "ismrmrd/meta.h"
#include "hoNDBSpline.h"
#include "ismrmrd/xml.h"

#include <numeric>
#include <queue>

#ifdef USE_OMP
#include <omp.h>
#endif 

namespace Gadgetron{

    namespace {

        struct Intervals {
            std::vector<float> intervals;
            std::vector<size_t> cycle_starts;
        };

         Intervals find_intervals(const std::vector<float>& time_stamps, bool first_beat_on_trigger) {

            float previous = -100.0;
            std::vector<float> intervals;
            std::vector<size_t> cycle_starts;
            for (size_t i = 0; i < time_stamps.size(); i++) {
                if ((time_stamps[i] < previous) || (first_beat_on_trigger && i == 0)) {
                    cycle_starts.push_back(i);
                } else if (i > 0) {
                    intervals.push_back(time_stamps[i] - time_stamps[i - 1]);
                }
                previous = time_stamps[i];
            }
            return {intervals,cycle_starts};
        };

    }


    void PhysioInterpolationGadget::process(Core::InputChannel<Core::Image<std::complex<float>>>& in,
                                            Core::OutputChannel& out) {


        ISMRMRD::EncodingLimits e_limits = header.encoding[0].encodingLimits;
        auto slc_limit = e_limits.slice ? e_limits.slice->maximum+1 : 1;


        auto buffers = std::map<int, std::queue<Core::Image<std::complex<float>>>>{};
        auto time_stamp_buffer = std::map<int, std::vector<float>>{};

        for (auto [hdr,data,meta] : in ){
            buffers[hdr.slice].emplace(hdr,data,meta);
            time_stamp_buffer[hdr.slice].push_back( (float)(hdr.physiology_time_stamp[physiology_time_index]) );
            out.push(Core::Image<std::complex<float>>{hdr,std::move(data),std::move(meta)});
        }

            for (auto& [slc,buffer] : buffers )
            {
                GDEBUG("Processing slice: %d ... \n", slc);
                GDEBUG("Number of items on Q: %d\n", buffer.size());
                GDEBUG("Image with attribute flag : %d\n", bool(std::get<2>(buffer.front())));



                const auto& time_stamp = time_stamp_buffer[slc];

                auto [intervals,cycle_stars] = find_intervals(time_stamp, first_beat_on_trigger);

                if ( intervals.empty() ) continue;

                std::sort(intervals.begin(),intervals.end());

                float mean_interval = std::accumulate(intervals.begin(), intervals.end(),0.0f)/intervals.size();
                float median_interval = intervals[(intervals.size()>>1)];

                std::vector<float> cycle_lengths;
                float count = 0;
                for (size_t i = 1; i < cycle_starts.size(); i++)
                {
                    float clength = time_stamps_[slc][cycle_starts[i]-1] + median_interval - time_stamps_[slc][cycle_starts[i]];
                    cycle_lengths.push_back(clength);
                }

                if ( cycle_lengths.empty() )
                {
                    size_t phs = time_stamps_[slc].size();
                    float clength = time_stamps_[slc][phs-1];
                    cycle_lengths.push_back(clength);
                }

                std::sort(cycle_lengths.begin(),cycle_lengths.end());
                float mean_cycle_length = std::accumulate(cycle_lengths.begin(), cycle_lengths.end(), 0.0f)/cycle_lengths.size();
                float median_cycle_length = cycle_lengths[(cycle_lengths.size()>>1)];

                GDEBUG("We have %d full cyles, first one starting at %d\n", cycle_starts.size()-1, cycle_starts[0]);
                GDEBUG("Mean/Median frame width %f/%f\n", mean_interval,median_interval);
                GDEBUG("Mean/Median cycle_length %f/%f\n", mean_cycle_length,median_cycle_length);

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
                    relative_cycle_time.push_back(time_stamps_[slc][i]/cycle_lengths[current_cycle-1] + current_cycle);
                }

                //Make a temporary list of all the data pointers from the Q
                std::vector< ISMRMRD::ImageHeader* > hptrs;
                std::vector< hoNDArray< std::complex<float> > * > aptrs;
                std::vector< ISMRMRD::MetaContainer* > attribptrs;

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
                        GDEBUG("Failed to cast data on Q, bailing out\n");
                        buffer_[slc]->flush();
                        return GADGET_FAIL;
                    }

                    hptrs.push_back(tmpm1->getObjectPtr());
                    aptrs.push_back(tmpm2->getObjectPtr());

                    if ( image_with_attrib_ )
                    {
                        GadgetContainerMessage< ISMRMRD::MetaContainer > * tmpm3 = 
                            AsContainerMessage< ISMRMRD::MetaContainer >(entry->cont()->cont());

                        if ( !tmpm3 )
                        {
                            GDEBUG("Failed to cast data on Q, bailing out\n");
                            buffer_[slc]->flush();
                            return GADGET_FAIL;
                        }

                        attribptrs.push_back(tmpm3->getObjectPtr());
                    }
                }

                //Let's figure out which time points we would like to interpolate on:
                ///TODO: Deal with mode 1 and other future modes, we are only implementing mode 0 at the moment
                float phase_interval = 1.0f/static_cast<float>(phases_to_reconstruct_);
                float max_time = std::floor(relative_cycle_time[relative_cycle_time.size()-1]);
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
                std::vector< GadgetContainerMessage< ISMRMRD::MetaContainer> * > out_attrib;

                for (size_t i = 0; i < recon_cycle_time.size(); i++)
                {
                    GadgetContainerMessage<ISMRMRD::ImageHeader>* tmpm1 = new GadgetContainerMessage<ISMRMRD::ImageHeader>;
                    GadgetContainerMessage< hoNDArray< std::complex<float> > >* tmpm2 = new GadgetContainerMessage< hoNDArray< std::complex<float> > >;

                    tmpm1->cont(tmpm2);

                    (*tmpm1->getObjectPtr()) = (*hptrs[0]);
                    tmpm2->getObjectPtr()->create(aptrs[0]->dimensions());

                    out_heads.push_back(tmpm1);
                    out_data.push_back(tmpm2);

                    unsigned short current_cycle = static_cast<unsigned short>(std::floor(recon_cycle_time[i] + 0.0001));
                    unsigned short current_phase = static_cast<unsigned short>((recon_cycle_time[i]+0.0001-current_cycle)/(1.0/static_cast<float>(phases_to_reconstruct_)) + 0.0001);

                    tmpm1->getObjectPtr()->physiology_time_stamp[phys_time_index_] = static_cast<unsigned>(std::floor((recon_cycle_time[i]+0.0001-current_cycle)*cycle_lengths[current_cycle]));
                    tmpm1->getObjectPtr()->phase = current_phase;
                    tmpm1->getObjectPtr()->image_index = current_phase+1;
                    tmpm1->getObjectPtr()->image_series_index = current_cycle * 10 + tmpm1->getObjectPtr()->slice;

                    // make sure the phase is within the acquisition limit
                    if ( tmpm1->getObjectPtr()->phase+1 >= time_stamps_[slc].size() )
                    {
                        tmpm1->getObjectPtr()->phase = (uint16_t)(time_stamps_[slc].size()-1);
                    }

                    if ( image_with_attrib_ )
                    {
                        GadgetContainerMessage< ISMRMRD::MetaContainer >* tmpm3 = new GadgetContainerMessage< ISMRMRD::MetaContainer >;

                        tmpm2->cont(tmpm3);
                        (*tmpm3->getObjectPtr()) = (*attribptrs[0]);
                        out_attrib.push_back(tmpm3);

                        tmpm3->getObjectPtr()->set("PHS",      (long)tmpm1->getObjectPtr()->phase);
                        tmpm3->getObjectPtr()->set(GADGETRON_IMAGENUMBER, (long)tmpm1->getObjectPtr()->image_index);

                        tmpm3->getObjectPtr()->append(GADGETRON_DATA_ROLE, "PhysioInterp");

                        double cycle_length_in_ms = time_stamp_resolution_.value()*cycle_lengths[current_cycle];

                        std::ostringstream ostr;
                        if (slc_limit_ > 1)
                        {
                            ostr << "_SLC_" << tmpm1->getObjectPtr()->slice << "_RR" << cycle_length_in_ms << "ms";
                        }
                        else
                        {
                            ostr << "_RR" << cycle_length_in_ms << "ms";
                        }

                        std::string imageComment = "PhysioInterp" + ostr.str();
                        tmpm3->getObjectPtr()->append(GADGETRON_IMAGECOMMENT, imageComment.c_str());

                        std::string seqDescription = "_PhysioInterp" + ostr.str();
                        tmpm3->getObjectPtr()->append(GADGETRON_SEQUENCEDESCRIPTION, seqDescription.c_str());

                        tmpm3->getObjectPtr()->append(GADGETRON_IMAGEPROCESSINGHISTORY, "Interp");
                    }
                }

                //Let's interpolate the images
                size_t inelem = relative_cycle_time.size();
                size_t outelem = recon_cycle_time.size();
                size_t imageelem = aptrs[0]->get_number_of_elements();

                if ( (interp_method_ == "Spline") || (mode_ != 1) )
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
                else
                {
                    GadgetronTimer interptime("Interpolation Time using BSpline");

                    size_t SplineDegree = 5;

                    long long p;
#pragma omp parallel default(none) shared(SplineDegree, imageelem, inelem, outelem, aptrs, relative_cycle_time, recon_cycle_time, out_data) private(p)
                    {
                        hoNDArray< std::complex<float> > data_in(inelem);
                        hoNDArray< std::complex<float> > data_out(outelem);

                        hoNDArray< std::complex<float> > coeff(inelem);

                        hoNDBSpline< std::complex<float>, 1 > interp;

                        size_t i;

                        size_t num = relative_cycle_time.size();

#pragma omp for
                        for (p = 0; p < (long long)imageelem; p++)
                        {
                            //Get the input data for this pixel
                            for (i = 0; i < inelem; i++) data_in(i) = aptrs[i]->get_data_ptr()[p];

                            // compute the coefficient
                            interp.computeBSplineCoefficients(data_in, SplineDegree, coeff);

                            //Interpolate the data
                            for (i = 0; i < outelem; i++)
                            {
                                float x = (num-1)*(recon_cycle_time[i]-relative_cycle_time[0])/(relative_cycle_time[num-1] - relative_cycle_time[0]);
                                data_out(i) = interp.evaluateBSpline(coeff.begin(), inelem, SplineDegree, 0, x);
                            }

                            //Copy it to the images
                            for (i = 0; i < outelem; i++) out_data[i]->getObjectPtr()->get_data_ptr()[p] = data_out[i];
                        }
                    }
                }

                //Send out the images
                for (size_t i = 0; i < out_heads.size(); i++)
                {
                    if (this->next()->putq(out_heads[i]) < 0)
                    {
                        GDEBUG("Unable to put data on next Gadgets Q\n");
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

        GadgetContainerMessage<ISMRMRD::MetaContainer>* m3 = 0;
        if (m2)
        {
            m3 = AsContainerMessage<ISMRMRD::MetaContainer>(m2->cont());
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
            GadgetContainerMessage< ISMRMRD::MetaContainer >* attrib = new GadgetContainerMessage< ISMRMRD::MetaContainer >;
            (*attrib->getObjectPtr()) = *m3->getObjectPtr();
            img->cont(attrib);
        }

        uint16_t slc = header->getObjectPtr()->slice;

        if (buffer_[slc]->enqueue_tail(header) < 0)
        {
            GDEBUG("Failed to add image to buffer\n");
            header->release();
            return GADGET_FAIL;
        }

        time_stamps_[slc].push_back( (float)(m1->getObjectPtr()->physiology_time_stamp[phys_time_index_]) );

        if (this->next()->putq(m1) < 0)
        {
            GDEBUG("Unable to put data on next Gadgets Q\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    GADGET_FACTORY_DECLARE(PhysioInterpolationGadget)
}
