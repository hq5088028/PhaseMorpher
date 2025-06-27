#pragma once
#include "ConcentrationFieldPoints.h"
#include "../Model_Params.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
namespace pf {
	enum ConcentrationFieldEquation {
		CFE_Const = 0, CFE_TotalConcentration, CFE_GrandPotential
	};
	namespace con_parameters {
		// Declare the models that need to be activated
		inline bool is_con_field_on = false;
		inline ConcentrationFieldEquation con_equation = ConcentrationFieldEquation::CFE_Const;
		// - concentration field data
		inline Mesh_Boundry<ConcentrationFieldPoint> concentration_field;
		// - 
		inline int con_number = 0;
		// - 
		inline REAL CON_MAX_VARIATION = 0.0;
	}
}