#ifndef BART_HELPERS_H_INCLUDED
#define BART_HELPERS_H_INCLUDED

#include "bartgadget.h"
#include "bart/bart_embed_api.h"

namespace internal {
     namespace fs = Gadgetron::fs;
     
     class ScopeGuard
     {
     public:
	  ScopeGuard(fs::path p) : p_(std::move(p))
	       {
		    char buf[1024] = { '\0' };
		    auto* ptr = getcwd(buf, 1024);
		    cwd_ = std::string(ptr);
		    auto r = chdir(p_.c_str());
	       }
	  ~ScopeGuard()
	       {
		    if (is_active_) {
			 fs::remove_all(p_);
		    }
		    auto r = chdir(cwd_.c_str());
		    deallocate_all_mem_cfl();
	       }

	  void dismiss() { is_active_ = false; }
     private:
	  bool is_active_;
	  const fs::path p_;
	  fs::path cwd_;
     };
     
     fs::path generate_unique_folder(const fs::path& working_directory);

     void ltrim(std::string &str);
     void rtrim(std::string &str);
     void trim(std::string &str);
	
     std::string get_output_filename(const std::string& command_line);
} // namespace internal

#endif //BART_HELPERS_H_INCLUDED
