#pragma once
#include "../../data_struct/RotationMatrix.h"
#include <vector>
namespace pf {
	class GrainsOrientations
	{
	public:
		GrainsOrientations() = default;
		~GrainsOrientations() {
			orientations.clear();
		}
		GrainsOrientations& operator=(const GrainsOrientations& n) {
			rotation_gauge = n.rotation_gauge;
			orientations = n.orientations;
			return *this;
		}
		void init(RotationGauge _rotation_gauge, int phi_number) {
			rotation_gauge = _rotation_gauge;
			orientations.resize(phi_number);
			for (int index = 0; index < phi_number; index++) {
				orientations[index][0] = REAL(0.0);
				orientations[index][1] = REAL(0.0);
				orientations[index][2] = REAL(0.0);
			}
		}
		void set_phi_orientation(int phi_index, Vector3 radian) {
			orientations[phi_index] = radian;
		}
		Vector3 get_phi_orientation(int phi_index) {
			return orientations[phi_index];
		}
		Matrix3x3 get_phi_rotationMatrix(int phi_index) {
			return RotationMatrix::rotationMatrix(orientations[phi_index], rotation_gauge);
		}
		RotationGauge rotation_gauge;
		std::vector<Vector3> orientations;
	};
}