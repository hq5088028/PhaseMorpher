#pragma once
#include "TemperatureFieldPoints.h"
#include "../Model_Params.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
namespace pf {
	enum TemperatureFieldEquation {
		TFE_Const = 0, TFE_Standard
	};
	namespace temp_parameters {
		// Declare the models that need to be activated
		inline bool is_temp_field_on = false;
		inline TemperatureFieldEquation temp_equation = TemperatureFieldEquation::TFE_Const;
		// - temperature field data
		inline Mesh_Boundry<TemperatureFieldPoint> temperature_field;
		// - 
		inline REAL TEMP_MAX_VARIATION = 0.0;




	}
}