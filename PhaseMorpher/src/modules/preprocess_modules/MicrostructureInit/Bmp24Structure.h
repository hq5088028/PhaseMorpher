#pragma once
#include "../../../base/sysfiles.h"
#include "../MicrostructureInit_Params.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../../tools/BMP24Reader.h"
#include "../../input_modules/inputfiles/selectFile.h"
#include "../../input_modules/ioFiles_Params.h"
namespace pf {
	namespace bmp24_structure {
		void generate_structure_from_BMP_pic();
		void init();
	}
}