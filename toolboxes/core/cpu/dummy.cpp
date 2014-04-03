//
// THIS IS A TEMPORARY FILE
//
// This file is to be removed from the repository once Hui merges his branch into mem_ops/development
// Currently this is the only .cpp in this folder, and one needs to be present to satisfy cmake and generate the .lib
//

#ifdef WIN32
#include <stdio.h>

namespace Gadgetron {

	void __declspec(dllexport) __this_is_a_temp_dummy_to_force_lib_file__(void){
		printf("\n\nINSIDE DUMMY\n\n");
	}
}

#endif