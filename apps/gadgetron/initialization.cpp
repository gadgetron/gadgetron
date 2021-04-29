
#include "initialization.h"

#ifdef FORCE_LIMIT_OPENBLAS_NUM_THREADS
#include <cblas.h>
#endif
#include <locale>
namespace Gadgetron::Server {

    void configure_blas_libraries() {

#ifdef FORCE_LIMIT_OPENBLAS_NUM_THREADS
        /*
         * Certain distributions (Ubuntu) doesn't ship OpenMP OpenBLAS, so we have to make do with
         * a version compiled for pthreads. It'll work - maybe even without degrading performance
         * a lot, but only if we limit the number of threads to one.
         */
        openblas_set_num_threads(1);
#endif

    }

    void set_locale() {
        try {
           std::locale::global(std::locale(""));
        } catch (...) {
            std::locale::global(std::locale::classic());
        }

    }
}
