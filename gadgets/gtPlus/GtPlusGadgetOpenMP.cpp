
#include "GtPlusGadgetOpenMP.h"

namespace Gadgetron
{

#ifdef USE_OMP

bool prepOpenMP()
{
    try
    {
        GADGET_MSG_DEPRECATED("--> OpenMP info <--");
        GADGET_MSG_DEPRECATED("--------------------------------------------------------");

        int numOpenMPProcs = omp_get_num_procs();
        GADGET_MSG_DEPRECATED("GtPlusRecon, numOpenMPProcs : " << numOpenMPProcs);

        #ifndef WIN32
            #ifndef GCC_OLD_FLAG
                int maxOpenMPLevels = omp_get_max_active_levels();
                GADGET_MSG_DEPRECATED("GtPlusRecon, maxOpenMPLevels : " << maxOpenMPLevels);
            #endif // GCC_OLD_FLAG
        #endif // WIN32

        int maxOpenMPThreads = omp_get_max_threads();
        GADGET_MSG_DEPRECATED("GtPlusRecon, maxOpenMPThreads : " << maxOpenMPThreads);

        if ( numOpenMPProcs != maxOpenMPThreads )
        {
            GADGET_MSG_DEPRECATED("GtPlusRecon, numOpenMPProcs != maxOpenMPThreads , hyperthreading must be disabled ... ");
            omp_set_num_threads(numOpenMPProcs);
        }

        // omp_set_nested(1);
        int allowOpenMPNested = omp_get_nested();
        GADGET_MSG_DEPRECATED("GtPlusRecon, allowOpenMPNested : " << allowOpenMPNested);

        #ifdef WIN32
            GADGET_MSG_DEPRECATED("----------------------------------");
            GADGET_MSG_DEPRECATED("GtPlus, set thread affinity ... ");

            /// lock the threads
            #pragma omp parallel default(shared)
            {
                int tid = omp_get_thread_num();
                DWORD_PTR mask = (1 << tid);
                GADGET_MSG_DEPRECATED("thread id : " << tid << " - mask : " << mask);
                SetThreadAffinityMask( GetCurrentThread(), mask );
            }
        #endif // WIN32

        GADGET_MSG_DEPRECATED("--------------------------------------------------------");
    }
    catch(...)
    {
        GADGET_ERROR_MSG("Errors in GtPlus prepOpenMP() ... ");
        return false;
    }

    return true;
}

#else

bool prepOpenMP()
{
    return true;
}

#endif // USE_OMP

}
