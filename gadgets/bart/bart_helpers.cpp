#include "bart_helpers.h"

#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/tokenizer.hpp>

#include <random>

// =============================================================================

namespace fs = Gadgetron::fs;

fs::path internal::generate_unique_folder(const fs::path& working_directory)
{
     typedef std::chrono::system_clock clock_t;
	  
     char buff[80];
     auto now = clock_t::to_time_t(clock_t::now());
     std::strftime(buff, sizeof(buff), "%H_%M_%S__", std::localtime(&now));
     std::random_device rd;
     auto time_id(buff + std::to_string(std::uniform_int_distribution<>(1, 10000)(rd)));
     // Get the current process ID
     auto threadId = boost::lexical_cast<std::string>(boost::this_thread::get_id());
     auto threadNumber(0UL);
     sscanf(threadId.c_str(), "%lx", &threadNumber);
     return  working_directory / ("bart_"
				  + time_id
				  + "_"
				  + std::to_string(threadNumber));
}

// ========================================================================

void internal::ltrim(std::string &str)
{
     str.erase(str.begin(), std::find_if(str.begin(), str.end(),
					 [](int s) {return !std::isspace(s); }));
}

void internal::rtrim(std::string &str)
{
     str.erase(std::find_if(str.rbegin(), str.rend(),
			    [](int s) {return !std::isspace(s);}).base(),
	       str.end());
}

void internal::trim(std::string &str)
{
     ltrim(str);
     rtrim(str);
}

std::string internal::get_output_filename(const std::string& command_line)
{
     boost::char_separator<char> sep(" ");
     std::vector<std::string> tokens;
     boost::tokenizer<boost::char_separator<char>> tokenizer(command_line, sep);
     for (auto& tok : tokenizer) {
	  tokens.push_back(tok);
     }

     // Find begin of argument list
     auto num_args(tokens.size()-1);
     auto it(tokens.cbegin());
     if (*it == "bart") {
	  ++it;
	  --num_args;
     }
     // it now points to the name of the BART command

     if (*it == "bitmask"
	 || *it == "estdelay"
	 || *it == "estdims"
	 || *it == "estshift"
	 || *it == "estvar"
	 || *it == "nrmse"
	 || *it == "sdot"
	 || *it == "show"
	 || *it == "version") {
	  return "";
     }
     else {
	  const auto end(tokens.cend());
	  const auto arg1(*(it+1));
	  const auto arg1_not_option(arg1[0] != '-');

	  /* 
	   * Essentially, the algorithm considers only the BART functions where
	   * the output argument might not be the last one.
	   *
	   * For these commands, if a sufficient number of argument is provided
	   * there is no ambiguity possible since we can always find the last
	   * "flag-like" looking argument on the command line and then infer
	   * the position of the output from there.
	   * If there might be an ambiguity, we look at the first argument of 
	   * the BART command and depending on whether it is a flag or not, 
	   * we can infer the position of the output argument on the command
	   * line.
	   */
	  if (*it == "ecalib") {
	       if (num_args > 3) {
		    auto it2 = --tokens.cend();
		    for (; it2 != it && (*it2)[0] != '-'; --it2);

		    const char& first((*it2)[0]);
		    const char& second((*it2)[1]);
		    if (first == '-'
			&& (it2->size() > 2 // handle cases like -t12 (ie. -t 12)
			    || second == 'S'
			    || second == 'W'
			    || second == 'I'
			    || second == '1'
			    || second == 'P'
			    || second == 'a')) {
			 std::cout << "ecalib: flag\n";
			 return *(it2 + 2);
		    }
		    else {
			 return *(it2 + 3);
		    }
	       }
	       else if (num_args == 3 && arg1_not_option) {
		    return *(end - 2);			 
	       }

	       // Minimum arguments to ecalib
	       // 2       kspace sens
	       // 3 -I    kspace sens
	       // 3       kspace sens ev-maps
	       // 4 -t 11 kspace sens
	       // 4 -P    kspace sens ev-maps
	       // 4 -t11  kspace sens ev-maps
	  }
	  else if (*it == "ecaltwo") {
	       if (num_args > 6) {
		    auto it2 = --tokens.cend();
		    for (; it2 != it && (*it2)[0] != '-'; --it2);
		    
		    const char& first((*it2)[0]);
		    const char& second((*it2)[1]);
		    if (first == '-'
			&& (it2->size() > 2 // handle cases like -t12 (ie. -t 12)
			    || second == 'S')) {
			 return *(it2 + 5);
		    }
		    else {
			 return *(it2 + 6);
		    }
	       }
	       else if (num_args == 6 && arg1[0] != '-') {
		    return *(end - 2);		 
	       }
	       // Minimum arguments to ecaltwo
	       // 6      x y z <input> <sensitivities> [<ev_maps>]
	       // 6 -S   x y z <input> <sensitivities>
	       // 7 -m 3 x y z <input> <sensitivities>
	       // 7 -S   x y z <input> <sensitivities> [<ev_maps>]
	  }
	  else if (*it == "nlinv") {
	       if (num_args > 3) {
		    auto it2 = --tokens.end();
		    for (; it2 != it && (*it2)[0] != '-'; --it2);

		    const char& first((*it2)[0]);
		    const char& second((*it2)[1]);
		    if (first == '-'
			&& (it2->size() > 2 // handle cases like -t12 (ie. -t 12)
			    || second == 'c'
			    || second == 'N'
			    || second == 'U'
			    || second == 'g'
			    || second == 'S')) {
			 return *(it2 + 2);
		    }
		    else {
			 return *(it2 + 3);
		    }
	       }
	       else if (num_args == 3 && arg1[0] != '-') {
		    return *(end - 2);
	       }
	       // Minimum arguments to nlinv
	       // 2       kspace output
	       // 3       kspace output sens
	       // 3  -S   kspace output
	       // 4  -i d kspace output
	       // 4  -S   kspace output sens
	  }
	  else if (*it == "whiten") {
	       if (num_args > 5) {
		    auto it2 = --tokens.end();
		    for (; it2 != it && (*it2)[0] != '-'; --it2);

		    const char& first((*it2)[0]);
		    const char& second((*it2)[1]);
		    if (first == '-'
			&& (it2->size() > 2 // handle cases like -t12 (ie. -t 12)
			    || second == 'n')) {
			 return *(it2 + 3);
		    }
		    else {
			 return *(it2 + 4);
		    }
	       }
	       else if(num_args == 5 && arg1_not_option) {
		    return *(end - 3);
	       }
	       else if (num_args == 4 && arg1_not_option) {
		    return *(end - 2);
	       }
	       // Minimum arguments to whiten
	       // 3         input ndata output
	       // 4         input ndata output optmat
	       // 4 -n      input ndata output
	       // 5 -o asd  input ndata output
	       // 5 -n      input ndata output optmat
	       // 5         input ndata output optmat covar_out
	  }
     }
     return tokens.back();
}
