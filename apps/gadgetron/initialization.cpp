
#include "log.h"
#include "initialization.h"

#include <cstdlib>
#include <string>

#ifdef __clang__
    #define unary_function  __unary_function
#endif

#include <boost/algorithm/string.hpp>

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


    void check_environment_variables() {

        auto get_policy = []() -> std::string {
            auto raw = std::getenv("OMP_WAIT_POLICY");
            return boost::algorithm::to_lower_copy(raw ? std::string(raw) : std::string());
        };

        if ("passive" != get_policy()) {
            GWARN_STREAM("Environment variable 'OMP_WAIT_POLICY' not set to 'PASSIVE'.");
            GWARN_STREAM("Gadgetron may experience serious performance issues under heavy load " <<
                         "(multiple simultaneous reconstructions, etc.)")
        }
    }

    void set_locale() {
        try {
           std::locale::global(std::locale(""));
        } catch (...) {
            std::locale::global(std::locale::classic());
        }
    }
}
