#pragma once
#include "../ConcentrationField_Params.h"
namespace pf {
	namespace grand_potential_functions {
		vector<double> pre_calculation_grand_potential_functional(double dt);

		double solve_grand_potential_functional(double dt);

	}
}