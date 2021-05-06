#include "log.h"
#include "paths.h"

#if defined _WIN32 || _WIN64
#include <libloaderapi.h>
#else
extern "C" {
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
}
#endif

namespace {

    namespace fs = boost::filesystem;

#if defined __APPLE__
    std::string get_executable_path() {
        char path[PATH_MAX];
        char resolved[PATH_MAX];
        uint32_t size = sizeof(path);

        if ((_NSGetExecutablePath(path, &size) == 0) && (realpath(path, resolved) != NULL)) {
            return std::string(resolved);
        } else {
            throw std::runtime_error("Could not determine location of Gadgetron binary.");
        }
    }
#elif defined _WIN32 || _WIN64
    #define MAX_GADGETRON_HOME_LENGTH 1024

    std::string get_executable_path() {
        // Full path to the executable (including the executable file)
        char buffer[MAX_GADGETRON_HOME_LENGTH];

        // When passing NULL to GetModuleHandle, it returns handle of exe itself
        HMODULE hModule = GetModuleHandle(NULL);

        if (hModule == NULL) {
            throw std::runtime_error("Could not determine location of Gadgetron binary.");
        }

        GetModuleFileName(hModule, buffer, sizeof(buffer));

        return std::string(buffer);
    }

    boost::filesystem::path get_data_directory(){
         auto appdata = std::getenv("APPDATA");

            return boost::filesystem::path(appdata) / "gadgetron";

    }
#else
    std::string get_executable_path(size_t buffer_size = 1024) {

        auto buffer = std::make_unique<char[]>(buffer_size);

        ssize_t len = readlink("/proc/self/exe", buffer.get(), buffer_size);

        if (len < 0) {
            throw std::runtime_error("Failed to read /proc/self/exe - cannot determine Gadgetron binary path.");
        }

        if (len == buffer_size) {
            // Allocated buffer was probably too small. Try again with a bigger buffer.
            return get_executable_path(buffer_size * 2);
        }

        return std::string(buffer.get(), size_t(len));
    }

    boost::filesystem::path get_data_directory(){
        auto home_dir = std::getenv("HOME");
        if (!home_dir){
            home_dir = getpwuid(getuid())->pw_dir;
        }
        return boost::filesystem::path(home_dir) / ".gadgetron";
    }
#endif
}


namespace Gadgetron::Server {

#ifdef _WIN32
    const boost::filesystem::path default_working_folder() {
        return "c:/temp/gadgetron/";
    }
#else
    const boost::filesystem::path default_working_folder() {
        return "/tmp/gadgetron/";
    }
#endif

    const boost::filesystem::path default_gadgetron_home() {

        const char *home = std::getenv("GADGETRON_HOME");

        if (home != nullptr) {
            return fs::path(home);
        }

        fs::path executable_path = get_executable_path();

        GDEBUG_STREAM("Executable path: " << executable_path);

        fs::path gadgetron_home = executable_path
                .parent_path()
                .parent_path();

        GDEBUG_STREAM("Default Gadgetron home: " << gadgetron_home);

        return gadgetron_home;
    }

    const boost::filesystem::path default_database_folder() {
        return default_storage_folder() / "database";
    }

    const boost::filesystem::path default_storage_folder() {
        return get_data_directory() / "storage";
    }
}


