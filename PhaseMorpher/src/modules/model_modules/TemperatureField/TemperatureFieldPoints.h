#pragma once
#include "../../../base/MACRO_DEF.h"
#include "../../data_struct/VectorMatrix.h"
namespace pf {

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

}