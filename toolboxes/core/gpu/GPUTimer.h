/** file GPUTimer.h
    Utility to measure Cuda performance. 
*/

#ifndef __GPUTIMER_H
#define __GPUTIMER_H

#pragma once

#include <iostream>
#include <string>
#include <cuda_runtime_api.h>

namespace Gadgetron{

    class GPUTimer
    {
    public:
        GPUTimer() : name_("GPUTimer"), timing_in_destruction_(true)
        {
            start();
        }

        GPUTimer(bool timing) : name_("GPUTimer"), timing_in_destruction_(timing)
        {
            if ( timing_in_destruction_ )
            {
                start();
            }
        }

        GPUTimer(const char* name) : name_(name), timing_in_destruction_(true)
        {
            start();
        }

        virtual ~GPUTimer() 
        {
            if ( timing_in_destruction_ )
            {
                stop();
            }
        }

        virtual void start()
        {
            cudaEventCreate(&start_event_);
            cudaEventCreate(&stop_event_);
            cudaEventRecord( start_event_, 0 );
        }

        virtual void stop()
        {
            float time;
            cudaEventRecord( stop_event_, 0 );
            cudaEventSynchronize( stop_event_ );
            cudaEventElapsedTime( &time, start_event_, stop_event_ );
            cudaEventDestroy( start_event_ );
            cudaEventDestroy( stop_event_ );

            GDEBUG_STREAM(name_ << ": " << time << " ms" << std::endl);
        }

        void set_timing_in_destruction(bool timing) { timing_in_destruction_ = timing; }

        cudaEvent_t start_event_;
        cudaEvent_t stop_event_;

        std::string name_;
        bool timing_in_destruction_;
    };
}
#endif //__GPUTIMER_H
