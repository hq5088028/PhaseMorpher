#pragma once
#include "../Model_Params.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
namespace pf {
	enum ConcentrationFieldEquation {
		CFE_Const = 0, CFE_TotalConcentration, CFE_GrandPotential
	};
	namespace model_parameters {
		// Declare the models that need to be activated
		inline bool is_con_field_on = false;
		inline ConcentrationFieldEquation con_equation = ConcentrationFieldEquation::CFE_Const;
		// - 
		inline int con_property_number = 0;
		inline int con_number = 0;
		// - 
		inline vector<int> con_property;
		// - 
		inline REAL CON_MAX_VARIATION = 0.0;
		// -
	}
}