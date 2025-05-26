#pragma once
#include "../data_struct/Mesh_0.h"
#include "../data_struct/Point_0.h"
#include "../data_struct/RotationMatrix.h"
#include "../MainIterator_Params.h"
#include "../../tools/mathTools.h"
#include "../input_modules/ioFiles_Params.h"
namespace pf {
	enum Dimension { One_Dimension, Two_Dimension, Three_Dimension };
	namespace mesh_parameters {
		// main mesh size
		inline int MESH_NX = 1;
		inline int MESH_NY = 1;
		inline int MESH_NZ = 1;
		inline Dimension dimention = Dimension::One_Dimension;
		// main mesh boundary condition
		inline BoundaryCondition x_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition y_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition z_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition x_up = BoundaryCondition::PERIODIC;
		inline BoundaryCondition y_up = BoundaryCondition::PERIODIC;
		inline BoundaryCondition z_up = BoundaryCondition::PERIODIC;
		// grid size
		inline REAL delt_r = 1.0;
	}
	namespace time_parameters {
		// real time
		inline REAL Real_Time = 0.0;
		// time interval for each simulation step
		inline REAL delt_t = 1.0;
	}
	namespace materials_system {
		// phases in this simulation
		inline std::vector<std::string> PHASES;
		inline int phi_property(std::string phi_name) {
			for (int index = 0; index < PHASES.size(); index++)
				if (phi_name.compare(PHASES[index]) == 0)
					return index;
			std::string error = "ERROR, phase name: " + phi_name + " has not been defined !";
			WriteError(error);
			SYS_PROGRAM_STOP;
		}
		// components in this simulation
		inline std::vector<std::string> COMPONENTS;
		inline int con_property(std::string con_name) {
			for (int index = 0; index < COMPONENTS.size(); index++)
				if (con_name.compare(COMPONENTS[index]) == 0)
					return index;
			std::string error = "ERROR, component name: " + con_name + " has not been defined !";
			WriteError(error);
			SYS_PROGRAM_STOP;
		}

	}
	namespace simulation_mesh {
		// - phase field data
		inline Mesh_Boundry<PhaseFieldPoint> phase_field;
		// - concentration field data
		inline Mesh_Boundry<ConcentrationFieldPoint> concentration_field;
		// - temperature field data
		inline Mesh_Boundry<TemperatureFieldPoint> temperature_field;

	}
}