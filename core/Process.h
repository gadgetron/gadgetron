#pragma once

#ifdef __clang__
   #define  unary_function  __unary_function
#endif

#include <mutex>
#include <boost/process.hpp>

namespace Gadgetron::Process {

    namespace detail{
        std::lock_guard<std::mutex> obtain_process_lock();
    }

    template<class ... ARGS>
    boost::process::child child(ARGS&&... args){
        auto guard = detail::obtain_process_lock();
        return boost::process::child(std::forward<ARGS>(args)...);
    }

    template<class... ARGS>
    auto system(ARGS&&... args){
        auto guard = detail::obtain_process_lock();
        return boost::process::system(std::forward<ARGS>(args)...);
    }

}
