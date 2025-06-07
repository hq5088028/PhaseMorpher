#pragma once
#include "ConcentrationField_Params.h"
namespace pf {

	inline void init_concentration_field_equations() {
		int control_equation = 0;
		WriteDebugFile("# Solver.ControlEquation.ConcentrationField = 0 - Const , 1 - TotalConcentration, 2 - GrandPotential, ... \n");
		if (InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.ConcentrationField", control_equation, true)) {
			con_parameters::is_con_field_on = true;
			con_parameters::con_equation = ConcentrationFieldEquation(control_equation);
			switch (con_parameters::con_equation)
			{
			case ConcentrationFieldEquation::CFE_TotalConcentration: {
				break;
			}
			case ConcentrationFieldEquation::CFE_GrandPotential: {
				break;
			}
			default:
				break;
			}
		}

	}

}