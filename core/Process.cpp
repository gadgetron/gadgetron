//
// Created by dch on 6/23/21.
//

#include "Process.h"

static std::mutex process_mutex;

std::lock_guard<std::mutex> Gadgetron::Process::detail::obtain_process_lock() {
    return std::lock_guard(process_mutex);
}
