#pragma once
#include "../Model_Params.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
namespace pf {
	enum TemperatureFieldEquation {
		TFE_Const = 0, TFE_Standard
	};
	struct TemperatureFieldPoint
	{
		REAL temp;
		REAL lap_temp;
		Vector3 grad_temp;
		~TemperatureFieldPoint() {

		}
		TemperatureFieldPoint& operator=(const TemperatureFieldPoint& n) {
			temp = n.temp;
			lap_temp = n.lap_temp;
			grad_temp = n.grad_temp;
			return *this;
		}
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