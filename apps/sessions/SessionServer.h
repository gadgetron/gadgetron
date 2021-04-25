#pragma once
#include <string>
#include <thread>
#include <filesystem>
namespace Gadgetron::Sessions {

    std::thread start_session_server( unsigned int port, const std::filesystem::path& database_folder, const std::filesystem::path& blob_folder );
}