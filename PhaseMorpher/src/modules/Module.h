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
#include <vector>
namespace pf {
	struct Solver_Module
	{
		// executes in the preprocess
		void(*exec_pre_i)();
		void(*exec_pre_ii)();
		void(*exec_pre_iii)();
		// executes in the model scheme
		void(*exec_i)();
		void(*exec_ii)();
		void(*exec_iii)();
		// executes in the postprocess
		void(*exec_pos_i)();
		void(*exec_pos_ii)();
		void(*exec_pos_iii)();
		// delete module
		void(*deinit)();
	};

	// module list
	inline std::vector<Solver_Module> module_list;

	// to creat a new module
	inline void load_a_new_module(void(*exec_pre_i)(), void(*exec_pre_ii)(), void(*exec_pre_iii)(),
		void(*exec_i)(), void(*exec_ii)(), void(*exec_iii)(), 
		void(*exec_pos_i)(), void(*exec_pos_ii)(), void(*exec_pos_iii)(), void(*deinit)()) {
		Solver_Module _module;
		_module.exec_pre_i = exec_pre_i;
		_module.exec_pre_ii = exec_pre_ii;
		_module.exec_pre_iii = exec_pre_iii;
		_module.exec_i = exec_i;
		_module.exec_ii = exec_ii;
		_module.exec_iii = exec_iii;
		_module.exec_pos_i = exec_pos_i;
		_module.exec_pos_ii = exec_pos_ii;
		_module.exec_pos_iii = exec_pos_iii;
		_module.deinit = deinit;
		module_list.push_back(_module);
	}

	inline void default_module_function() {
		return;
	}

}