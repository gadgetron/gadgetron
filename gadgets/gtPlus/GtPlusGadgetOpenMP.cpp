
#include "GtPlusGadgetOpenMP.h"

namespace Gadgetron
{

#ifdef USE_OMP

bool prepOpenMP()
{
    try
    {
        GDEBUG_STREAM("--> OpenMP info <--");
        GDEBUG_STREAM("--------------------------------------------------------");

        int numOpenMPProcs = omp_get_num_procs();
        GDEBUG_STREAM("GtPlusRecon, numOpenMPProcs : " << numOpenMPProcs);

        #ifndef WIN32
            #ifndef GCC_OLD_FLAG
                int maxOpenMPLevels = omp_get_max_active_levels();
                GDEBUG_STREAM("GtPlusRecon, maxOpenMPLevels : " << maxOpenMPLevels);
            #endif // GCC_OLD_FLAG
        #endif // WIN32

        int maxOpenMPThreads = omp_get_max_threads();
        GDEBUG_STREAM("GtPlusRecon, maxOpenMPThreads : " << maxOpenMPThreads);

        if ( numOpenMPProcs != maxOpenMPThreads )
        {
            GDEBUG_STREAM("GtPlusRecon, numOpenMPProcs != maxOpenMPThreads , hyperthreading must be disabled ... ");
            omp_set_num_threads(numOpenMPProcs);
        }

        // omp_set_nested(1);
        int allowOpenMPNested = omp_get_nested();
        GDEBUG_STREAM("GtPlusRecon, allowOpenMPNested : " << allowOpenMPNested);

        #ifdef WIN32
            GDEBUG_STREAM("----------------------------------");
            GDEBUG_STREAM("GtPlus, set thread affinity ... ");

            /// lock the threads
            #pragma omp parallel default(shared)
            {
                int tid = omp_get_thread_num();
                DWORD_PTR mask = (1 << tid);
                GDEBUG_STREAM("thread id : " << tid << " - mask : " << mask);
                SetThreadAffinityMask( GetCurrentThread(), mask );
            }
        #endif // WIN32

        GDEBUG_STREAM("--------------------------------------------------------");
    }
    catch(...)
    {
        GERROR_STREAM("Errors in GtPlus prepOpenMP() ... ");
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
