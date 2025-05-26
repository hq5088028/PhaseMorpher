#pragma once
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Module.h"
#include "../MainIterator_Params.h"
#include "../model_modules/Model_Params.h"
#include "../model_modules/PhaseField/PhaseField_Params.h"
namespace pf {
	namespace data_statistics_functions {

		inline PhaseFieldPoint statistical_phi() {
			PhaseFieldPoint phi_info;
			phi_info.init(model_parameters::phi_number);
			int SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (int x = simulation_mesh::phase_field.X_BEGIN(); x <= simulation_mesh::phase_field.X_END(); x++)
				for (int y = simulation_mesh::phase_field.Y_BEGIN(); y <= simulation_mesh::phase_field.Y_END(); y++)
					for (int z = simulation_mesh::phase_field.Z_BEGIN(); z <= simulation_mesh::phase_field.Z_END(); z++) {
						PhaseFieldPoint& point = simulation_mesh::phase_field(x, y, z);
						for (int index = 0; index < model_parameters::phi_number; index++)
							phi_info.phi[index] += point.phi[index] / SIZE;
					}
			return phi_info;
		}



	}

	inline void init() {

	}

}