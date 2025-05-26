#pragma once
#include "../../../base/sysfiles.h"
#include "../../model_modules/Model_Params.h"
#include "../MicrostructureInit_Params.h"
#include "../../input_modules/ioFiles_Params.h"
namespace pf {
	namespace geometry_structure {
		// - init geometry structure and batch method
		void init();
		// - 
		void definiteNucleation();

	}
}