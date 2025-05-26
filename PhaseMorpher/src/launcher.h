#pragma once
#include "modules/MainIterator.h"
namespace pf {

	inline void pm_run(int argc, char* argv[]) {

		main_iterator::init_modules(argc, argv);

		main_iterator::run();
		
	}

}