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
			phi_info.init(phi_parameters::phi_number);
			int SIZE = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ;
			for (int x = phi_parameters::phase_field.X_BEGIN(); x <= phi_parameters::phase_field.X_END(); x++)
				for (int y = phi_parameters::phase_field.Y_BEGIN(); y <= phi_parameters::phase_field.Y_END(); y++)
					for (int z = phi_parameters::phase_field.Z_BEGIN(); z <= phi_parameters::phase_field.Z_END(); z++) {
						PhaseFieldPoint& point = phi_parameters::phase_field(x, y, z);
						for (int index = 0; index < phi_parameters::phi_number; index++)
							phi_info.phi[index] += point.phi[index] / SIZE;
					}
			return phi_info;
		}



	}

	inline void init() {

	}

}