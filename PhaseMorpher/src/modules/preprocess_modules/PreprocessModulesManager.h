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
#include "../../base/sysfiles.h"
#include "MicrostructureInit.h"
namespace pf {
	inline void init_preprocess_modules() {
		microstructure_init::init();
	}
}
