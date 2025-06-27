#pragma once
#include "../Model_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
#include "PhaseFieldPoints.h"
#include "GrainsOrientations.h"
namespace pf {
	enum PhaseFieldEquation {
		PFE_Const = 0, PFE_AllenCahn_Pairwise
	};
	namespace phi_parameters {
		// Declare the models that need to be activated
		inline bool is_phi_field_on = false;
		inline PhaseFieldEquation phi_equation = PhaseFieldEquation::PFE_Const;
		// - phase field data
		inline Mesh_Boundry<PhaseFieldPoint> phase_field;
		// - 
		inline int phi_property_number = 0;
		inline int phi_number = 0;
		// - 
		inline vector<int> phi_property;
		// - 
		inline REAL PHI_MAX_VARIATION = 0.0;
		// -
	}

}