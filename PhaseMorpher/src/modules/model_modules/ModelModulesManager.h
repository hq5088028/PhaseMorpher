/*
This file is a part of the microstructure intelligent design software project.

Created:     Qi Huang 2025.03

Modified:    Qi Huang 2025.03;

Copyright (c) 2019-2025 Science center for phase diagram, phase transition, material intelligent design and manufacture, Central South University, China

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free
	Software Foundation, either version 3 of the License, or (at your option) any later version.

	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
	FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

	You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#pragma once
#include "../Module.h"
#include "Model_Params.h"
#include "../MainIterator_Params.h"
#include "PhaseField/PhaseFieldManager.h"
#include "ConcentrationField/ConcentrationFieldManager.h"
#include "TemperatureField/TemperatureFieldManager.h"
namespace pf {
	inline void init_basic_time_mesh_parameters();
	inline void init_simulation_system();
	inline void init_model_modules() {
		// parallel
		InputFileReader::get_instance()->read_int_value("Solver.Loop.OpenMP_Thread", main_iterator::OpenMP_Thread_Counts, true);
		// - mesh and time parameters
		init_basic_time_mesh_parameters();
		// - init simulation system
		init_simulation_system();
		// - init all physical fields
		init_phase_field_equations();
		init_concentration_field_equations();
		init_temperature_field_equations();
	}

	inline void init_basic_time_mesh_parameters() {
		// - mesh and time parameters
		InputFileReader::get_instance()->read_int_value("Solver.Loop.begin_step", main_iterator::ITE_Begin_Step, true);
		InputFileReader::get_instance()->read_int_value("Solver.Loop.end_step", main_iterator::ITE_End_Step, true);
		InputFileReader::get_instance()->read_REAL_value("Solver.Loop.RealTime.init", time_parameters::Real_Time, true);
		InputFileReader::get_instance()->read_REAL_value("Solver.Loop.dt", time_parameters::delt_t, true);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.Nx", mesh_parameters::MESH_NX, true);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.Ny", mesh_parameters::MESH_NY, true);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.Nz", mesh_parameters::MESH_NZ, true);
		// dimension rules
		if (mesh_parameters::MESH_NX < 1 || mesh_parameters::MESH_NY < 1 || mesh_parameters::MESH_NZ < 1) {
			string error_report = "> ERROR, one edge length (Nx or Ny or Nz) of domain defined error, domain does not exist ! Please set Nx > 0 and Ny > 0 and Nz > 0 !\n";
			WriteError(error_report);
			std::exit(0);
		}
		else if (mesh_parameters::MESH_NX == 1 && (mesh_parameters::MESH_NY > 1 || mesh_parameters::MESH_NZ > 1)) {
			string error_report = "> ERROR, edge length (Nx, Ny, Nz) of domain should be set (Nx > 1, Ny = 1, Nz = 1) in 1 dimension or (Nx > 1, Ny > 1, Nz = 1) in 2 dimension !\n";
			WriteError(error_report);
			std::exit(0);
		}
		else if (mesh_parameters::MESH_NX > 1 && mesh_parameters::MESH_NY == 1 && mesh_parameters::MESH_NZ > 1) {
			string error_report = "> ERROR, edge length (Nx, Ny, Nz) of domain should be set as (Nx > 1, Ny > 1, Nz = 1) in 2 dimension !\n";
			WriteError(error_report);
			std::exit(0);
		}
		if (mesh_parameters::MESH_NY == 1 && mesh_parameters::MESH_NZ == 1)
			mesh_parameters::dimention = Dimension::One_Dimension;
		else if (mesh_parameters::MESH_NZ == 1)
			mesh_parameters::dimention = Dimension::Two_Dimension;
		else
			mesh_parameters::dimention = Dimension::Three_Dimension;
		InputFileReader::get_instance()->read_REAL_value("Solver.Mesh.dr", mesh_parameters::delt_r, true);
		WriteDebugFile("# Solver.Mesh.BoundaryCondition : 0 - FIXED , 1 - PERIODIC , 2 - ADIABATIC\n");
		int boundary_condition = 0;
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_up", boundary_condition, true);
		mesh_parameters::x_up = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.x_down", boundary_condition, true);
		mesh_parameters::x_down = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_up", boundary_condition, true);
		mesh_parameters::y_up = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.y_down", boundary_condition, true);
		mesh_parameters::y_down = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_up", boundary_condition, true);
		mesh_parameters::z_up = BoundaryCondition(boundary_condition);
		InputFileReader::get_instance()->read_int_value("Solver.Mesh.BoundaryCondition.z_dowm", boundary_condition, true);
		mesh_parameters::z_down = BoundaryCondition(boundary_condition);
	}
	inline void init_simulation_system() {
		// - phases
		WriteDebugFile("# Solver.SimulationSystem.PhaseNames = (name_0, name_1, ... ) \n");
		string phase_key = "Solver.SimulationSystem.PhaseNames", phase_input = "()";
		InputFileReader::get_instance()->read_string_value(phase_key, phase_input, true);
		vector<input_value> phase_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, phase_key, phase_input, true);
		for (int index = 0; index < phase_value.size(); index++) {
			bool is_exist = false;
			for (int index2 = 0; index2 < materials_system::PHASES.size(); index2++)
				if (materials_system::PHASES[index2].compare(phase_value[index].string_value) == 0)
					is_exist = true;
			if (is_exist) {
				WriteError("ERROR, Solver.SimulationSystem.PhaseNames = (), define repeated phase name ! ");
				SYS_PROGRAM_STOP;
			}
			materials_system::PHASES.push_back(phase_value[index].string_value);
		}
		// - components
		WriteDebugFile("# Solver.SimulationSystem.ComponentNames = (name_0, name_1, ... ) \n");
		string comp_key = "Solver.SimulationSystem.ComponentNames", comp_input = "()";
		InputFileReader::get_instance()->read_string_value(comp_key, comp_input, true);
		vector<input_value> comp_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_STRING, comp_key, comp_input, true);
		for (int index = 0; index < comp_value.size(); index++) {
			bool is_exist = false;
			for (int index2 = 0; index2 < materials_system::COMPONENTS.size(); index2++)
				if (materials_system::COMPONENTS[index2].compare(comp_value[index].string_value) == 0)
					is_exist = true;
			if (is_exist) {
				WriteError("ERROR, Solver.SimulationSystem.ComponentNames = (), define repeated component name ! ");
				SYS_PROGRAM_STOP;
			}
			materials_system::COMPONENTS.push_back(comp_value[index].string_value);
		}
	}
}
