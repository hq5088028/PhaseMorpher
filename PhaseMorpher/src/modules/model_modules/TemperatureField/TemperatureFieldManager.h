#pragma once
#include "TemperatureField_Params.h"
namespace pf {

	inline void init_temperature_field_equations() {
		int control_equation = 0;
		WriteDebugFile("# Solver.ControlEquation.TemperatureField = 0 - Const , 1 - Standard, ... \n");
		if (InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.TemperatureField", control_equation, true)) {
			model_parameters::temp_equation = TemperatureFieldEquation(control_equation);
			model_parameters::is_temp_field_on = true;
			switch (model_parameters::temp_equation)
			{
			case TemperatureFieldEquation::TFE_Standard:
				break;
			default:
				break;
			}
		}
	}

}