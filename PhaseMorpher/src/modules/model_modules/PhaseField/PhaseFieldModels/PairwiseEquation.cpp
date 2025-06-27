#include "PairwiseEquation.h"

namespace pf {
	namespace pair_wise_functions {
		using namespace pf::pair_wise_parameters;
		// - 
		Flag currentFlag(int x, int y, int z, PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			if (point.phi[phi_index] >= Phi_Num_Cut_Off && point.phi[phi_index] <= (1.0 - Phi_Num_Cut_Off))
				return pf_INTERFACE;
			if (point.phi[phi_index] < Phi_Num_Cut_Off) {
				if (diff_method == DifferenceMethod::FIVE_POINT) {
					if (phase_field->at(x + 1, y, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x - 1, y, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y + 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y - 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y, z + 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y, z - 1).phi[phi_index] >= Phi_Num_Cut_Off)
						return pf_NEAR_INTERFACE;
				}
				else {
					if (phase_field->at(x + 1, y, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x - 1, y, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y + 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y - 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y, z + 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y, z - 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x + 1, y + 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x + 1, y - 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x + 1, y, z + 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x + 1, y, z - 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x - 1, y + 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x - 1, y - 1, z).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x - 1, y, z + 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x - 1, y, z - 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y + 1, z + 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y + 1, z - 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y - 1, z + 1).phi[phi_index] >= Phi_Num_Cut_Off
						|| phase_field->at(x, y - 1, z - 1).phi[phi_index] >= Phi_Num_Cut_Off)
						return pf_NEAR_INTERFACE;
				}
			}
			else if (point.phi[phi_index] > Phi_Num_Cut_Off_R) {
				if (diff_method == DifferenceMethod::FIVE_POINT) {
					if (phase_field->at(x + 1, y, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x - 1, y, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y + 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y - 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y, z + 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y, z - 1).phi[phi_index] <= Phi_Num_Cut_Off_R)
						return pf_NEAR_INTERFACE;
				}
				else {
					if (phase_field->at(x + 1, y, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x - 1, y, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y + 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y - 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y, z + 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y, z - 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x + 1, y + 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x + 1, y - 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x + 1, y, z + 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x + 1, y, z - 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x - 1, y + 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x - 1, y - 1, z).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x - 1, y, z + 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x - 1, y, z - 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y + 1, z + 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y + 1, z - 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y - 1, z + 1).phi[phi_index] <= Phi_Num_Cut_Off_R
						|| phase_field->at(x, y - 1, z - 1).phi[phi_index] <= Phi_Num_Cut_Off_R)
						return pf_NEAR_INTERFACE;
				}
			}
			return pf_BULK;
		}
		// - 
		Vector3 normals(REAL alpha_phi, REAL beta_phi, Vector3 alpha_grad, Vector3 beta_grad) {
			Vector3 normal;
			normal[0] = alpha_grad[0] * beta_phi - alpha_phi * beta_grad[0];
			normal[1] = alpha_grad[1] * beta_phi - alpha_phi * beta_grad[1];
			normal[2] = alpha_grad[2] * beta_phi - alpha_phi * beta_grad[2];
			REAL length = sqrt(normal * normal);
			if (length < SYS_EPSILON) {
				normal[0] = 0.0;
				normal[1] = 0.0;
				normal[2] = 0.0;
			}
			else {
				normal[0] = normal[0] / length;
				normal[1] = normal[1] / length;
				normal[2] = normal[2] / length;
			};
			return normal;
		};
		// - mobility
		REAL Lij_pair_wise_const_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			int check = 0;
			for (auto block = const_block_Lij.begin(); block < const_block_Lij.end(); block++) {
				check = 0;
				for (auto phiIndex = block->begin(); phiIndex < block->end(); phiIndex++)
					if (*phiIndex == phi_alpha_index || *phiIndex == phi_beta_index)
						check++;
				if (check == 2)
					return const_Lij;
			}
			return 0.0;
		};
		REAL Lij_pair_wise_anisotropic_const_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			int check = 0; REAL Lij = 0.0;
			for (auto block = const_block_Lij.begin(); block < const_block_Lij.end(); block++) {
				check = 0;
				for (auto phiIndex = block->begin(); phiIndex < block->end(); phiIndex++)
					if (*phiIndex == phi_alpha_index || *phiIndex == phi_beta_index)
						check++;
				if (check == 2)
					Lij = const_Lij;
			}
			// - anisotropic
			int phi_alpha_property = phi_property[phi_alpha_index], phi_beta_property = phi_property[phi_beta_index],
				anisotropic_model = intMob_anisotropic_model_matrix[phi_alpha_property][phi_beta_property];
			if (anisotropic_model == Int_Mobility_Anisotropic::IMA_ISO) {
				return Lij;
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_CUBIC) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index) 
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] 
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 + intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_SUN) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- intMobAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- intMobAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_YANG) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + intMobAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
			}
			return Lij;
		};
		REAL Lij_pair_wise_const_cross_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			for (auto phiIndex = const_block_Lij[0].begin(); phiIndex < const_block_Lij[0].end(); phiIndex++)
				if (*phiIndex == phi_alpha_index) {
					for (auto phiIndex = const_block_Lij[1].begin(); phiIndex < const_block_Lij[1].end(); phiIndex++)
						if (*phiIndex == phi_beta_index)
							return const_Lij;
				}
				else if (*phiIndex == phi_beta_index) {
					for (auto phiIndex = const_block_Lij[1].begin(); phiIndex < const_block_Lij[1].end(); phiIndex++)
						if (*phiIndex == phi_alpha_index)
							return const_Lij;
				}
			return 0.0;
		};
		REAL Lij_pair_wise_anisotropic_const_cross_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			int check = 0; REAL Lij = 0.0;
			for (auto phiIndex = const_block_Lij[0].begin(); phiIndex < const_block_Lij[0].end(); phiIndex++)
				if (*phiIndex == phi_alpha_index) {
					for (auto phiIndex = const_block_Lij[1].begin(); phiIndex < const_block_Lij[1].end(); phiIndex++)
						if (*phiIndex == phi_beta_index)
							Lij = const_Lij;
				}
				else if (*phiIndex == phi_beta_index) {
					for (auto phiIndex = const_block_Lij[1].begin(); phiIndex < const_block_Lij[1].end(); phiIndex++)
						if (*phiIndex == phi_alpha_index)
							Lij = const_Lij;
				}
			// - anisotropic
			int phi_alpha_property = phi_property[phi_alpha_index], phi_beta_property = phi_property[phi_beta_index],
				anisotropic_model = intMob_anisotropic_model_matrix[phi_alpha_property][phi_beta_property];
			if (anisotropic_model == Int_Mobility_Anisotropic::IMA_ISO) {
				return Lij;
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_CUBIC) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property]
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 + intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_SUN) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- intMobAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- intMobAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_YANG) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + intMobAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
			}
			return Lij;
		};
		REAL Lij_pair_wise_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			return const_Lij;
		};
		REAL Lij_pair_wise_anisotropic_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			REAL Lij = const_Lij;
			// - anisotropic
			int phi_alpha_property = phi_property[phi_alpha_index], phi_beta_property = phi_property[phi_beta_index],
				anisotropic_model = intMob_anisotropic_model_matrix[phi_alpha_property][phi_beta_property];
			if (anisotropic_model == Int_Mobility_Anisotropic::IMA_ISO) {
				return Lij;
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_CUBIC) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property]
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 + intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_SUN) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- intMobAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- intMobAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_YANG) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + intMobAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
			}
			return Lij;
		};
		REAL Lij_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			return matirx_Lij[phi_alpha_index][phi_beta_index];
		};
		REAL Lij_anisotropic_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			// to be defined
			REAL Lij = matirx_Lij[phi_alpha_index][phi_beta_index];
			// - anisotropic
			int phi_alpha_property = phi_property[phi_alpha_index], phi_beta_property = phi_property[phi_beta_index],
				anisotropic_model = intMob_anisotropic_model_matrix[phi_alpha_property][phi_beta_property];
			if (anisotropic_model == Int_Mobility_Anisotropic::IMA_ISO) {
				return Lij;
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_CUBIC) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property]
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 + intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_SUN) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- intMobAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- intMobAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_YANG) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + intMobAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
			}
			return Lij;
		};
		REAL Lij_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			for (int index = 0; index < block_Lij.size(); index++)
				if ((block_Lij[index][0] <= phi_alpha_index && phi_alpha_index <= block_Lij[index][1])
					&& (block_Lij[index][0] <= phi_beta_index && phi_beta_index <= block_Lij[index][1]))
					return block_value_Lij[index];
			return 0.0;
		};
		REAL Lij_anisotropic_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			// to be defined
			REAL Lij = 0.0;
			for (int index = 0; index < block_Lij.size(); index++)
				if ((block_Lij[index][0] <= phi_alpha_index && phi_alpha_index <= block_Lij[index][1])
					&& (block_Lij[index][0] <= phi_beta_index && phi_beta_index <= block_Lij[index][1]))
					Lij = block_value_Lij[index];
			// - anisotropic
			int phi_alpha_property = phi_property[phi_alpha_index], phi_beta_property = phi_property[phi_beta_index],
				anisotropic_model = intMob_anisotropic_model_matrix[phi_alpha_property][phi_beta_property];
			if (anisotropic_model == Int_Mobility_Anisotropic::IMA_ISO) {
				return Lij;
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_CUBIC) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property]
					* (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 + intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_SUN) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0)
						- 30.0 * norm[2] * norm[2] + 3.0)
					- intMobAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0)
						- 315.0 * std::pow(norm[2], 4.0) + 105.0 * norm[2] * norm[2] - 5.0)
					- intMobAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1]
						+ 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0)
						- std::pow(norm[1], 6.0)));
			}
			else if (anisotropic_model == Int_Mobility_Anisotropic::IMA_HEX_YANG) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				Lij = Lij * REAL(1.0 - intMobAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3 * norm[2] * norm[2] - 1.0), 2.0)
					- intMobAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + intMobAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
			}
			return Lij;
		};
		// - interface energy
		REAL xi_ab_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			return const_xi_ab;
		};
		REAL xi_ab_anisotropic_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			REAL loc_xi_ab = const_xi_ab;
			int phi_alpha_property = phi_property[phi_alpha_index], phi_beta_property = phi_property[phi_beta_index];
			int anisotropic_model = intEn_anisotropic_model_matrix[phi_alpha_property][phi_beta_property];
			if (anisotropic_model == Int_Energy_Anisotropic::IEA_ISO) {
				return loc_xi_ab;
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_CUBIC) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_HEX_BOETTGER) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 - intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] +
					15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_HEX_SUN) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					+ intEnAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0) - 30.0 * norm[2] * norm[2] + 3.0)
					+ intEnAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0) - 315.0 * std::pow(norm[2], 4.0)
						+ 105.0 * norm[2] * norm[2] - 5.0)
					+ intEnAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0) - std::pow(norm[1], 6.0)));
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_HEX_YANG) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3.0 * norm[2] * norm[2] - 1.0), 2.0)
					+ intEnAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + intEnAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
			}
			return loc_xi_ab;
		};
		REAL xi_abc_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index, int phi_gamma_index) {
			return const_xi_abc;
		};
		REAL xi_ab_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			return matirx_xi_ab[phi_property[phi_alpha_index]][phi_property[phi_beta_index]];
		};
		REAL xi_ab_anisotropic_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			REAL loc_xi_ab = matirx_xi_ab[phi_property[phi_alpha_index]][phi_property[phi_beta_index]];
			// - anisotropic
			int phi_alpha_property = phi_property[phi_alpha_index], phi_beta_property = phi_property[phi_beta_index];
			int anisotropic_model = intEn_anisotropic_model_matrix[phi_alpha_property][phi_beta_property];
			// - anisotropic
			if (anisotropic_model == Int_Energy_Anisotropic::IEA_ISO) {
				return loc_xi_ab;
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_CUBIC) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * (1.5 - 2.5 * (std::pow(norm[0], 4.0) + std::pow(norm[1], 4.0) + std::pow(norm[2], 4.0))));
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_HEX_BOETTGER) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 - intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * (std::pow(norm[0], 6.0) - std::pow(norm[1], 6.0) -
					15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] +
					15.0 * std::pow(norm[1], 4.0) * norm[0] * norm[0] +
					(5.0 * std::pow(norm[2], 4.0) - 5.0 * std::pow(norm[2], 2.0) + std::pow(norm[2], 6.0))));
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_HEX_SUN) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * sqrt(5.0 / 16.0 / PI) * (3.0 * norm[2] * norm[2] - 1.0)
					+ intEnAniso2_matrix[phi_alpha_property][phi_beta_property] * 3.0 / 16.0 / sqrt(PI) * (35.0 * std::pow(norm[2], 4.0) - 30.0 * norm[2] * norm[2] + 3.0)
					+ intEnAniso3_matrix[phi_alpha_property][phi_beta_property] * sqrt(13.0 / PI) / 32.0 * (231.0 * std::pow(norm[2], 6.0) - 315.0 * std::pow(norm[2], 4.0)
						+ 105.0 * norm[2] * norm[2] - 5.0)
					+ intEnAniso4_matrix[phi_alpha_property][phi_beta_property] * sqrt(6006.0 / PI) / 64.0 * (std::pow(norm[0], 6.0)
						- 15.0 * std::pow(norm[0], 4.0) * norm[1] * norm[1] + 15.0 * norm[0] * norm[0] * std::pow(norm[1], 4.0) - std::pow(norm[1], 6.0)));
			}
			else if (anisotropic_model == Int_Energy_Anisotropic::IEA_HEX_YANG) {
				Vector3 norm = grains_orientation.get_phi_rotationMatrix(phi_beta_index) * (grains_orientation.get_phi_rotationMatrix(phi_alpha_index)
					* normals(point.phi[phi_alpha_index], point.phi[phi_beta_index], point.grad_phi[phi_alpha_index], point.grad_phi[phi_beta_index]));
				loc_xi_ab = loc_xi_ab * REAL(1.0 + intEnAniso1_matrix[phi_alpha_property][phi_beta_property] * std::pow((3.0 * norm[2] * norm[2] - 1.0), 2.0)
					+ intEnAniso2_matrix[phi_alpha_property][phi_beta_property] * std::pow((norm[0] * norm[0] * norm[0] - 3.0 * norm[0] * norm[1] * norm[1]), 2.0)
					* std::pow((9.0 * norm[2] * norm[2] - 1.0 + intEnAniso3_matrix[phi_alpha_property][phi_beta_property]), 2.0));
			}
			return loc_xi_ab;
		};
		REAL xi_abc_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index, int phi_gamma_index) {
			return matirx_xi_abc[phi_property[phi_alpha_index]][phi_property[phi_beta_index]][phi_property[phi_gamma_index]];
		};
		// - crack resistence
		REAL Gc_acc(PhaseFieldPoint& point, PairwisePoint& point_ex) {
			REAL gc = 0.0, sum_h = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_index = point_ex.active_index[index];
				if (phi_index != PAIRWISE_ACC_STOP && phi_property[phi_index] != crack_phase_property) {
					REAL h = interpolation_func(point.phi[phi_index]);
					sum_h += h;
					//< Gc of each phase
					gc += h * phase_Gc[phi_property[phi_index]];
				}
			}
			if (sum_h < SYS_EPSILON)
				return 0.0;
			else
				return gc / sum_h;
		}
		REAL dGc_dphi_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			REAL dh_dphi = dinterpolation_func(point.phi[phi_index]), sum_h = 0.0, gc_phase = phase_Gc[phi_property[phi_index]];
			for (int index = 0; index < phi_acc_number; index++) {
				if (point_ex.active_index[index] != PAIRWISE_ACC_STOP &&
					phi_property[point_ex.active_index[index]] != crack_phase_property)
					sum_h += interpolation_func(point.phi[point_ex.active_index[index]]);
			}
			if (sum_h < SYS_EPSILON)
				return 0.0;
			return dh_dphi * gc_phase / sum_h - interpolation_func(point.phi[phi_index]) * gc_phase / sum_h / sum_h * dh_dphi;
		}
		// - interface energy model
		REAL dfint_dphi_grad_S1996_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			REAL grad = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_bIndex = point_ex.active_index[index];
				if (phi_bIndex != PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
					grad += 2 * interface_width * _xi_ab(point, phi_index, phi_bIndex) * (point.grad_phi[phi_bIndex] * point.grad_phi[phi_bIndex] * point.phi[phi_index]
						- point.grad_phi[phi_index] * point.grad_phi[phi_bIndex] * point.phi[phi_bIndex] + point.phi[phi_index] * point.phi[phi_bIndex] * point.lap_phi[phi_bIndex] 
						- point.phi[phi_bIndex] * point.phi[phi_bIndex] * point.lap_phi[phi_index]);
				}
			}
			return grad;
		};
		REAL dfint_dphi_grad_S1999_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			REAL grad = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_bIndex = point_ex.active_index[index];
				if (phi_bIndex != PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
					grad += interface_width * _xi_ab(point, phi_index, phi_bIndex) * point.lap_phi[phi_bIndex];
				}
			}
			return grad;
		};
		REAL dfint_dphi_pot_Nobstacle_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			REAL pot = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_bIndex = point_ex.active_index[index];
				if (phi_bIndex != PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
					pot += 16 / interface_width / PI / PI * _xi_ab(point, phi_index, phi_bIndex) * point.phi[phi_bIndex];
					for (int index2 = index + 1; index2 < phi_acc_number; index2++) {
						int phi_gIndex = point_ex.active_index[index2];
						if (phi_gIndex != PAIRWISE_ACC_STOP && phi_gIndex != phi_index)
							pot += _xi_abc(point, phi_index, phi_bIndex, phi_gIndex) * point.phi[phi_bIndex] * point.phi[phi_gIndex] / interface_width;
					}
				}
			}
			return pot;
		};
		REAL dfint_dphi_pot_Nwell_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			REAL pot = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_bIndex = point_ex.active_index[index];
				if (phi_bIndex != PAIRWISE_ACC_STOP && phi_bIndex != phi_index) {
					pot += 18 / interface_width * _xi_ab(point, phi_index, phi_bIndex) * point.phi[phi_index] * point.phi[phi_bIndex] * point.phi[phi_bIndex];
					for (int index2 = index + 1; index2 < phi_acc_number; index2++) {
						int phi_gIndex = point_ex.active_index[index2];
						if (phi_gIndex != PAIRWISE_ACC_STOP && phi_gIndex != phi_index)
							pot += 2 / interface_width * _xi_abc(point, phi_index, phi_bIndex, phi_gIndex) * point.phi[phi_index]
							* point.phi[phi_bIndex] * point.phi[phi_gIndex] * point.phi[phi_bIndex] * point.phi[phi_gIndex];
					}
				}
			}
			return pot;
		};
		REAL dfint_dphi_solid_D2016_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			// defined by users
			REAL crack = 0.0, abs_crackphase_grad = 0.0, crackphase_fraction = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_bIndex = point_ex.active_index[index];
				if (phi_bIndex != PAIRWISE_ACC_STOP && phi_property[phi_bIndex] == crack_phase_property) {
					abs_crackphase_grad = abs(point.grad_phi[phi_bIndex] * point.grad_phi[phi_bIndex]);
					crackphase_fraction = point.phi[phi_bIndex];
				}
			}
			crack = dGc_dphi_acc(point, point_ex, phi_index) * (crack_int_width * abs_crackphase_grad + crackphase_fraction * crack_K / crack_int_width);
			// copy from interfaceEnergy.h,		Steinbach1996
			REAL grad = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_bIndex = point_ex.active_index[index];
				if (phi_bIndex != PAIRWISE_ACC_STOP && phi_bIndex != phi_index && 
					point_ex.flag[phi_bIndex] && phi_property[phi_bIndex] != crack_phase_property) {
					grad += 2 * interface_width * _xi_ab(point, phi_index, phi_bIndex) * (point.grad_phi[phi_bIndex] * point.grad_phi[phi_bIndex] * point.phi[phi_index]
						- point.grad_phi[phi_index] * point.grad_phi[phi_bIndex] * point.phi[phi_bIndex] + point.phi[phi_index] * point.phi[phi_bIndex] * point.lap_phi[phi_bIndex]
						- point.phi[phi_bIndex] * point.phi[phi_bIndex] * point.lap_phi[phi_index]);
				}
			}
			// copy from interfaceEnergy.h		Nestler_Obstacle
			REAL pot = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_bIndex = point_ex.active_index[index];
				if (phi_bIndex != PAIRWISE_ACC_STOP && (phi_bIndex != phi_index || point_ex.flag[phi_bIndex] 
					|| phi_property[phi_bIndex] != crack_phase_property)) {
					pot += 16 / interface_width / PI / PI * _xi_ab(point, phi_index, phi_bIndex) * point.phi[phi_bIndex];
					for (int index2 = 0; index2 < phi_acc_number; index2++) {
						int phi_gIndex = point_ex.active_index[index2];
						if (phi_gIndex != PAIRWISE_ACC_STOP && phi_gIndex != phi_index
							&& point_ex.flag[phi_gIndex] && phi_property[phi_gIndex] != crack_phase_property)
							pot += _xi_abc(point, phi_index, phi_bIndex, phi_gIndex) * point.phi[phi_bIndex] * point.phi[phi_gIndex] / interface_width;
					}
				}
			}
			return crack + grad + pot;
		};
		REAL dfint_dphi_crack_D2016_acc(int x, int y, int z, PhaseFieldPoint& point, int phi_index) {
			REAL gc = phase_field_Gc(x, y, z);
			Vector3 delt_gc = { (phase_field_Gc(x + 1, y, z) - phase_field_Gc(x - 1, y, z)) / 2 / *delt_r,
			(phase_field_Gc(x, y + 1, z) - phase_field_Gc(x, y - 1, z)) / 2 / *delt_r,
			(phase_field_Gc(x, y, z + 1) - phase_field_Gc(x, y, z - 1)) / 2 / *delt_r };
			return gc / crack_int_width * crack_K - 2 * crack_int_width * (point.grad_phi[phi_index] * delt_gc + gc * point.lap_phi[phi_index]);
		};
		void pairwise_normalize_daniel_crack_acc(PhaseFieldPoint& point, PairwisePoint& point_ex) {
			bool is_damaged = false;
			for (int index = 0; index < phi_acc_number; index++) {
				int phi_index = point_ex.active_index[index];
				if (phi_index != PAIRWISE_ACC_STOP && phi_property[phi_index] == crack_phase_property
					&& point.phi[phi_index] > crack_num_cut_off)
					is_damaged = true;
			}
			if (is_damaged) {
				for (int index = 0; index < phi_acc_number; index++) {
					int phi_index = point_ex.active_index[index];
					if (phi_index != PAIRWISE_ACC_STOP) {
						if (phi_property[phi_index] == crack_phase_property) {
							point.phi[phi_index] = 1.0;
							point_ex.int_increment[phi_index] = 0.0;
							point_ex.bulk_increment[phi_index] = 0.0;
						}
						else {
							point.phi[phi_index] = 0.0;
							point_ex.int_increment[phi_index] = 0.0;
							point_ex.bulk_increment[phi_index] = 0.0;
						}
					}
				}
			}
		}
		REAL dfint_dphi_pairwise_acc(int x, int y, int z, PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1996_acc(point, point_ex, phi_index) + dfint_dphi_pot_Nobstacle_acc(point, point_ex, phi_index);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1996 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1996_acc(point, point_ex, phi_index) + dfint_dphi_pot_Nwell_acc(point, point_ex, phi_index);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Obstacle) {
				return dfint_dphi_grad_S1999_acc(point, point_ex, phi_index) + dfint_dphi_pot_Nobstacle_acc(point, point_ex, phi_index);
			}
			else if (interface_gradient == Int_Gradient::Steinbach_1999 && interface_potential == Int_Potential::Nestler_Well) {
				return dfint_dphi_grad_S1999_acc(point, point_ex, phi_index) + dfint_dphi_pot_Nwell_acc(point, point_ex, phi_index);
			}
			else if (interface_gradient == Int_Gradient::DanielCrack_G2016 || interface_potential == Int_Potential::DanielCrack_P2016) {
				if (phi_property[phi_index] == crack_phase_property)
					return dfint_dphi_crack_D2016_acc(x, y, z, point, phi_index);
				else
					return dfint_dphi_solid_D2016_acc(point, point_ex, phi_index);
			}
			else {
				return 0.0;
			}
		};
		REAL dfint_dphi_S2009(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index) {
			return _xi_ab(point, phi_alpha_index, phi_beta_index) * (PI * PI / 2 / interface_width * (point.phi[phi_alpha_index] - point.phi[phi_beta_index])
				+ interface_width * (point.phi[phi_beta_index] * point.lap_phi[phi_alpha_index] - point.phi[phi_alpha_index] * point.lap_phi[phi_beta_index]));
		};
		REAL dfbulk_dphi_pairwise_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index) {
			REAL dfbulk_sum = 0.0;
			for (auto dfbulk = delt_Fbulk_delt_phi.begin(); dfbulk < delt_Fbulk_delt_phi.end(); dfbulk++)
				dfbulk_sum += (*dfbulk)(point, point_ex, phi_index);
			return dfbulk_sum;
		}
		REAL source_pairwise_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_alpha_index, int phi_beta_index) {
			REAL source_sum = 0.0;
			for (auto source = source_alpha_beta.begin(); source < source_alpha_beta.end(); source++)
				source_sum += (*source)(point, point_ex, phi_alpha_index, phi_beta_index);
			return source_sum;
		}

		void pairwise_normalize_normal(pf::PhaseFieldPoint& point, PairwisePoint& point_ex) {
			REAL scale = 1.0, increment = 0.0;
			for (int index = 0; index < phi_number; index++) {
				if (point_ex.flag[index]) {
					increment = *delt_t * (point_ex.int_increment[index] + point_ex.bulk_increment[index]);
					if (isTwoNumEquality(increment, 0.0))
						continue;
					else if ((point.phi[index] + increment) > 1) {
						REAL p_scale = (1 - point.phi[index]) / increment;
						if (p_scale < scale)
							scale = p_scale;
					}
					else if ((point.phi[index] + increment) < 0) {
						REAL p_scale = abs(point.phi[index] / increment);
						if (p_scale < scale)
							scale = p_scale;
					}
				}
			}
			for (int index = 0; index < phi_number; index++) {
				point_ex.int_increment[index] *= scale;
				point_ex.bulk_increment[index] *= scale;
			}
		}

		inline void pairwise_normalize_acc(pf::PhaseFieldPoint& point, PairwisePoint& point_ex) {
			REAL scale = 1.0, increment = 0.0;
			for (int index = 0; index < phi_acc_number; index++) {
				int acc_index = point_ex.active_index[index];
				if (acc_index != PAIRWISE_ACC_STOP) {
					increment = *delt_t * (point_ex.int_increment[acc_index] + point_ex.bulk_increment[acc_index]);
					if (isTwoNumEquality(increment, 0.0))
						continue;
					else if ((point.phi[acc_index] + increment) > 1) {
						REAL p_scale = (1 - point.phi[acc_index]) / increment;
						if (p_scale < scale)
							scale = p_scale;
					}
					else if ((point.phi[acc_index] + increment) < 0) {
						REAL p_scale = abs(point.phi[acc_index] / increment);
						if (p_scale < scale)
							scale = p_scale;
					}
				}
			}
			for (int index = 0; index < phi_acc_number; index++) {
				int acc_index = point_ex.active_index[index];
				if (acc_index != PAIRWISE_ACC_STOP) {
					point_ex.int_increment[acc_index] *= scale;
					point_ex.bulk_increment[acc_index] *= scale;
				}
			}
		}

		//=========================================================================================================================================

		void pre_calculation_phi_pair_wise() {
#pragma omp parallel for
			for (int x = phase_field->X_BEGIN(); x <= phase_field->X_END(); x++)
				for (int y = phase_field->Y_BEGIN(); y <= phase_field->Y_END(); y++)
					for (int z = phase_field->Z_BEGIN(); z <= phase_field->Z_END(); z++) {
						PhaseFieldPoint& point = phase_field->at(x, y, z);
						PairwisePoint& point_ex = phase_field_pairwise(x, y, z);
						int actIndex = 0;
						for (int index = 0; index < phi_acc_number; index++)
							point_ex.active_index[index] = PAIRWISE_ACC_STOP;
						for (int index = 0; index < phi_acc_number; index++) {
							if (point_ex.flag[index] || point.phi[index] > Phi_Num_Cut_Off) {
								point_ex.active_index[actIndex] = index; // save the active phi index here , which edded by  PAIRWISE_ACC_STOP
								actIndex++;
							}
							if (point_ex.flag[index]) {
								point_ex.bulk_increment[index] = 0.0;
								point_ex.int_increment[index] = 0.0;
								point.grad_phi[index][0] = (phase_field->at(x + 1, y, z).phi[index] - phase_field->at(x - 1, y, z).phi[index]) / 2 / *delt_r;
								point.grad_phi[index][1] = (phase_field->at(x, y + 1, z).phi[index] - phase_field->at(x, y - 1, z).phi[index]) / 2 / *delt_r;
								point.grad_phi[index][2] = (phase_field->at(x, y, z + 1).phi[index] - phase_field->at(x, y, z - 1).phi[index]) / 2 / *delt_r;
								if (diff_method == DifferenceMethod::FIVE_POINT) {
									point.lap_phi[index] =
										(phase_field->at(x + 1, y, z).phi[index] + phase_field->at(x - 1, y, z).phi[index]
											+ phase_field->at(x, y + 1, z).phi[index] + phase_field->at(x, y - 1, z).phi[index]
											+ phase_field->at(x, y, z + 1).phi[index] + phase_field->at(x, y, z - 1).phi[index]
											- 6 * point.phi[index]) / *delt_r / *delt_r;
								}
								else if (diff_method == DifferenceMethod::NINE_POINT) {
									point.lap_phi[index] =
										((phase_field->at(x + 1, y, z).phi[index] + phase_field->at(x - 1, y, z).phi[index]
											+ phase_field->at(x, y + 1, z).phi[index] + phase_field->at(x, y - 1, z).phi[index]
											+ phase_field->at(x, y, z + 1).phi[index] + phase_field->at(x, y, z - 1).phi[index]) * 4
											+ phase_field->at(x - 1, y - 1, z).phi[index] + phase_field->at(x - 1, y + 1, z).phi[index]
											+ phase_field->at(x + 1, y - 1, z).phi[index] + phase_field->at(x + 1, y + 1, z).phi[index]
											+ phase_field->at(x - 1, y, z - 1).phi[index] + phase_field->at(x - 1, y, z + 1).phi[index]
											+ phase_field->at(x + 1, y, z - 1).phi[index] + phase_field->at(x + 1, y, z + 1).phi[index]
											+ phase_field->at(x, y - 1, z - 1).phi[index] + phase_field->at(x, y - 1, z + 1).phi[index]
											+ phase_field->at(x, y + 1, z - 1).phi[index] + phase_field->at(x, y + 1, z + 1).phi[index]
											- 36 * point.phi[index]) / 6 / *delt_r / *delt_r;
								}
							}
							else {
								point.grad_phi[index][0] = 0.0;
								point.grad_phi[index][1] = 0.0;
								point.grad_phi[index][2] = 0.0;
								point.lap_phi[index] = 0.0;
							}
						}
						// - for the crack propagation
						if (interface_gradient == Int_Gradient::DanielCrack_G2016 || interface_potential == Int_Potential::DanielCrack_P2016)
							phase_field_Gc(x, y, z) = Gc_acc(point, point_ex);
					}
#pragma omp parallel for
			for (int x = phase_field->X_BEGIN(); x <= phase_field->X_END(); x++)
				for (int y = phase_field->Y_BEGIN(); y <= phase_field->Y_END(); y++)
					for (int z = phase_field->Z_BEGIN(); z <= phase_field->Z_END(); z++) {
						PhaseFieldPoint& point = phase_field->at(x, y, z);
						PairwisePoint& point_ex = phase_field_pairwise(x, y, z);
						// interface energy
						if (interface_gradient == Int_Gradient::Steinbach_G2009 || interface_potential == Int_Potential::Steinbach_P2009) {
							for (int aIndex = 0; aIndex < phi_acc_number; aIndex++) {
								int alpha_index = point_ex.active_index[aIndex];
								if (alpha_index == PAIRWISE_ACC_STOP)
									continue;
								for (int bIndex = aIndex + 1; bIndex < phi_acc_number; bIndex++) {
									int beta_index = point_ex.active_index[bIndex];
									if (beta_index == PAIRWISE_ACC_STOP)
										continue;
									REAL int_incre_b_a = mobility(point, alpha_index, beta_index)
										/ interface_width * dfint_dphi_S2009(point, alpha_index, beta_index);
									if (is_phi_normalized) {
										if ((int_incre_b_a > SYS_EPSILON && (point.phi[alpha_index] > SYS_EPSILON_R || point.phi[beta_index] < SYS_EPSILON))
											|| (int_incre_b_a < -SYS_EPSILON && (point.phi[alpha_index] < SYS_EPSILON || point.phi[beta_index] > SYS_EPSILON_R)))
											int_incre_b_a = 0.0;
									}
									point_ex.int_increment[alpha_index] += int_incre_b_a;
									point_ex.int_increment[beta_index] -= int_incre_b_a;
#ifdef _DEBUG
									if (_isnan(int_incre_b_a)) {
										cout << "DEBUG: interface energy (pair-wise functions) error !" << endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
						}
						else {
							for (int aIndex = 0; aIndex < phi_acc_number; aIndex++) {
								int alpha_index = point_ex.active_index[aIndex];
								if (alpha_index == PAIRWISE_ACC_STOP)
									continue;
								for (int bIndex = aIndex + 1; bIndex < phi_acc_number; bIndex++) {
									int beta_index = point_ex.active_index[bIndex];
									if (beta_index == PAIRWISE_ACC_STOP)
										continue;
									REAL int_incre_b_a = mobility(point, alpha_index, beta_index) / interface_width
										* (dfint_dphi_pairwise_acc(x, y, z, point, point_ex, beta_index) - dfint_dphi_pairwise_acc(x, y, z, point, point_ex, alpha_index));
									if (is_phi_normalized) {
										if ((int_incre_b_a > SYS_EPSILON && (point.phi[alpha_index] > SYS_EPSILON_R || point.phi[beta_index] < SYS_EPSILON))
											|| (int_incre_b_a < -SYS_EPSILON && (point.phi[alpha_index] < SYS_EPSILON || point.phi[beta_index] > SYS_EPSILON_R)))
											int_incre_b_a = 0.0;
									}
									point_ex.int_increment[alpha_index] += int_incre_b_a;
									point_ex.int_increment[beta_index] -= int_incre_b_a;
#ifdef _DEBUG
									if (_isnan(int_incre_b_a)) {
										cout << "DEBUG: interface energy (pair-wise functions) error !" << endl;
										SYS_PROGRAM_STOP;
									}
#endif
								}
							}
						}
						// driving force and source term
						for (int aIndex = 0; aIndex < phi_acc_number; aIndex++) {
							int alpha_index = point_ex.active_index[aIndex];
							if (alpha_index == PAIRWISE_ACC_STOP)
								continue;
							for (int bIndex = aIndex + 1; bIndex < phi_acc_number; bIndex++) {
								int beta_index = point_ex.active_index[bIndex];
								if (beta_index == PAIRWISE_ACC_STOP)
									continue;
								if (point.phi[alpha_index] > PhiCon_Num_Cut_Off && point.phi[beta_index] > PhiCon_Num_Cut_Off) {
									REAL bulk_increment_b_a = mobility(point, alpha_index, beta_index)
										* (dfbulk_dphi_pairwise_acc(point, point_ex, beta_index) - dfbulk_dphi_pairwise_acc(point, point_ex, alpha_index))
										+ source_pairwise_acc(point, point_ex, alpha_index, beta_index);
									if (is_phi_normalized) {
										if ((bulk_increment_b_a > SYS_EPSILON && (point.phi[alpha_index] > SYS_EPSILON_R || point.phi[beta_index] < SYS_EPSILON))
											|| (bulk_increment_b_a < -SYS_EPSILON && (point.phi[alpha_index] < SYS_EPSILON || point.phi[beta_index] > SYS_EPSILON_R)))
											bulk_increment_b_a = 0.0;
									}
									point_ex.bulk_increment[alpha_index] += bulk_increment_b_a;
									point_ex.bulk_increment[beta_index] -= bulk_increment_b_a;
								}
							}
						}
						// numerical treatment to nomalize the phi
						if (is_phi_normalized) {
							pairwise_normalize_acc(point, point_ex);
							if (interface_gradient == Int_Gradient::DanielCrack_G2016 || interface_potential == Int_Potential::DanielCrack_P2016)
								pairwise_normalize_daniel_crack_acc(point, point_ex);
						}
					}
		}

		REAL solve_phi_pair_wise() {
			REAL MAX_PHI_INCREMENT = 0.0;
#pragma omp parallel for
			for (int x = phase_field->X_BEGIN(); x <= phase_field->X_END(); x++)
				for (int y = phase_field->Y_BEGIN(); y <= phase_field->Y_END(); y++)
					for (int z = phase_field->Z_BEGIN(); z <= phase_field->Z_END(); z++) {
						PhaseFieldPoint& point = phase_field->at(x, y, z);
						PairwisePoint& point_ex = phase_field_pairwise(x, y, z);
						bool phi_change = false;
						for (int index = 0; index < phi_acc_number; index++) {
							int phi_index = point_ex.active_index[index];
							if (phi_index != PAIRWISE_ACC_STOP && point_ex.flag[phi_index]) {
								phi_change = true;
								point.old_phi[phi_index] = point.phi[phi_index];
								point.phi[phi_index] += *delt_t * (point_ex.int_increment[phi_index] + point_ex.bulk_increment[phi_index]);
#ifdef _OPENMP
#pragma omp critical
#endif
								{
									if (abs(point.old_phi[phi_index] - point.phi[phi_index]) > MAX_PHI_INCREMENT)
										MAX_PHI_INCREMENT = abs(point.old_phi[phi_index] - point.phi[phi_index]);
								}
#ifdef _DEBUG
								if (_isnan(point.phi[phi_index])) {
									string error_report = "ERROR : Phase fraction is NaN in x = " + to_string(x) + " , y = " + to_string(y) 
										+ " , z = " + to_string(z) + " position, the phase index is " + to_string(phi_index) + "\n";
									cout << error_report << endl;
									SYS_PROGRAM_STOP;
								}
#endif
							}
						}
						if (phi_change) {
							// normalize the phi
							if (is_phi_normalized)
								point.normalize_phi();
							for (int index = 0; index < phi_number; index++) {
								// change _flag
								if (point.phi[index] >= Phi_Num_Cut_Off && point.phi[index] <= Phi_Num_Cut_Off_R) {
									if (point_ex.flag[index] != pf_INTERFACE) {
										point_ex.flag[index] = pf_INTERFACE;
										if (phase_field_pairwise(x + 1, y, z).flag[index] == pf_BULK)
											phase_field_pairwise(x + 1, y, z).flag[index] = pf_NEAR_INTERFACE;
										if (phase_field_pairwise(x - 1, y, z).flag[index] == pf_BULK)
											phase_field_pairwise(x - 1, y, z).flag[index] = pf_NEAR_INTERFACE;
										if (phase_field_pairwise(x, y + 1, z).flag[index] == pf_BULK)
											phase_field_pairwise(x, y + 1, z).flag[index] = pf_NEAR_INTERFACE;
										if (phase_field_pairwise(x, y - 1, z).flag[index] == pf_BULK)
											phase_field_pairwise(x, y - 1, z).flag[index] = pf_NEAR_INTERFACE;
										if (phase_field_pairwise(x, y, z + 1).flag[index] == pf_BULK)
											phase_field_pairwise(x, y, z + 1).flag[index] = pf_NEAR_INTERFACE;
										if (phase_field_pairwise(x, y, z - 1).flag[index] == pf_BULK)
											phase_field_pairwise(x, y, z - 1).flag[index] = pf_NEAR_INTERFACE;
										if (diff_method == DifferenceMethod::NINE_POINT) {
											if (phase_field_pairwise(x + 1, y + 1, z).flag[index] == pf_BULK)
												phase_field_pairwise(x + 1, y + 1, z).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x + 1, y - 1, z).flag[index] == pf_BULK)
												phase_field_pairwise(x + 1, y - 1, z).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x + 1, y, z + 1).flag[index] == pf_BULK)
												phase_field_pairwise(x + 1, y, z + 1).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x + 1, y, z - 1).flag[index] == pf_BULK)
												phase_field_pairwise(x + 1, y, z - 1).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x - 1, y + 1, z).flag[index] == pf_BULK)
												phase_field_pairwise(x - 1, y + 1, z).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x - 1, y - 1, z).flag[index] == pf_BULK)
												phase_field_pairwise(x - 1, y - 1, z).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x - 1, y, z + 1).flag[index] == pf_BULK)
												phase_field_pairwise(x - 1, y, z + 1).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x - 1, y, z - 1).flag[index] == pf_BULK)
												phase_field_pairwise(x - 1, y, z - 1).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x, y + 1, z + 1).flag[index] == pf_BULK)
												phase_field_pairwise(x, y + 1, z + 1).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x, y + 1, z - 1).flag[index] == pf_BULK)
												phase_field_pairwise(x, y + 1, z - 1).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x, y - 1, z + 1).flag[index] == pf_BULK)
												phase_field_pairwise(x, y - 1, z + 1).flag[index] = pf_NEAR_INTERFACE;
											if (phase_field_pairwise(x, y - 1, z - 1).flag[index] == pf_BULK)
												phase_field_pairwise(x, y - 1, z - 1).flag[index] = pf_NEAR_INTERFACE;
										}
									}
								}
								else if (point.phi[index] < Phi_Num_Cut_Off) {
									if (diff_method == DifferenceMethod::FIVE_POINT) {
										if (phase_field->at(x + 1, y, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x - 1, y, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y + 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y - 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y, z + 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y, z - 1).phi[index] >= Phi_Num_Cut_Off)
											point_ex.flag[index] = pf_NEAR_INTERFACE;
										else
											point_ex.flag[index] = pf_BULK;
									}
									else {
										if (phase_field->at(x + 1, y, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x - 1, y, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y + 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y - 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y, z + 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y, z - 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x + 1, y + 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x + 1, y - 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x + 1, y, z + 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x + 1, y, z - 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x - 1, y + 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x - 1, y - 1, z).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x - 1, y, z + 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x - 1, y, z - 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y + 1, z + 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y + 1, z - 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y - 1, z + 1).phi[index] >= Phi_Num_Cut_Off
											|| phase_field->at(x, y - 1, z - 1).phi[index] >= Phi_Num_Cut_Off)
											point_ex.flag[index] = pf_NEAR_INTERFACE;
										else
											point_ex.flag[index] = pf_BULK;
									}
								}
								else if (point.phi[index] > Phi_Num_Cut_Off_R) {
									if (diff_method == DifferenceMethod::FIVE_POINT) {
										if (phase_field->at(x + 1, y, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x - 1, y, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y + 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y - 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y, z + 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y, z - 1).phi[index] <= Phi_Num_Cut_Off_R)
											point_ex.flag[index] = pf_NEAR_INTERFACE;
										else
											point_ex.flag[index] = pf_BULK;
									}
									else {
										if (phase_field->at(x + 1, y, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x - 1, y, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y + 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y - 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y, z + 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y, z - 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x + 1, y + 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x + 1, y - 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x + 1, y, z + 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x + 1, y, z - 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x - 1, y + 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x - 1, y - 1, z).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x - 1, y, z + 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x - 1, y, z - 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y + 1, z + 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y + 1, z - 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y - 1, z + 1).phi[index] <= Phi_Num_Cut_Off_R
											|| phase_field->at(x, y - 1, z - 1).phi[index] <= Phi_Num_Cut_Off_R)
											point_ex.flag[index] = pf_NEAR_INTERFACE;
										else
											point_ex.flag[index] = pf_BULK;
									}
								}
							}
						}
					}
			return MAX_PHI_INCREMENT;
		}

		void boundary_condition_crack() {
			if (x_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(0, y, z)          = phase_field->at(MESH_NX, y, z);
						phase_field_pairwise(0, y, z) = phase_field_pairwise(MESH_NX, y, z);
						phase_field_Gc(0, y, z)       = phase_field_Gc(MESH_NX, y, z);
					}
			}
			else if (x_down == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(0, y, z) = phase_field->at(1, y, z);
						phase_field_pairwise(0, y, z) = phase_field_pairwise(1, y, z);
						phase_field_Gc(0, y, z) = phase_field_Gc(1, y, z);
					}
			}
			else {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(0, y, z);
					}
			}
			if (y_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, 0, z) = phase_field->at(x, MESH_NY, z);
						phase_field_pairwise(x, 0, z) = phase_field_pairwise(x, MESH_NY, z);
						phase_field_Gc(x, 0, z) = phase_field_Gc(x, MESH_NY, z);
					}
			}
			else if (y_down == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, 0, z) = phase_field->at(x, 1, z);
						phase_field_pairwise(x, 0, z) = phase_field_pairwise(x, 1, z);
						phase_field_Gc(x, 0, z) = phase_field_Gc(x, 1, z);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nz(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(x, 0, z);
					}
			}
			if (z_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, 0) = phase_field->at(x, y, MESH_NZ);
						phase_field_pairwise(x, y, 0) = phase_field_pairwise(x, y, MESH_NZ);
						phase_field_Gc(x, y, 0) = phase_field_Gc(x, y, MESH_NZ);
					}
			}
			else if (z_down == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, 0) = phase_field->at(x, y, 1);
						phase_field_pairwise(x, y, 0) = phase_field_pairwise(x, y, 1);
						phase_field_Gc(x, y, 0) = phase_field_Gc(x, y, 1);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						BoundaryCondition_PhiFix(x, y, 0);
					}
			}
			if (x_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(MESH_NX + 1, y, z) = phase_field->at(1, y, z);
						phase_field_pairwise(MESH_NX + 1, y, z) = phase_field_pairwise(1, y, z);
						phase_field_Gc(MESH_NX + 1, y, z) = phase_field_Gc(1, y, z);
					}
			}
			else if (x_up == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(MESH_NX + 1, y, z) = phase_field->at(MESH_NX, y, z);
						phase_field_pairwise(MESH_NX + 1, y, z) = phase_field_pairwise(MESH_NX, y, z);
						phase_field_Gc(MESH_NX + 1, y, z) = phase_field_Gc(MESH_NX, y, z);
					}
			}
			else {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(MESH_NX + 1, y, z);
					}
			}
			if (y_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, MESH_NY + 1, z) = phase_field->at(x, 1, z);
						phase_field_pairwise(x, MESH_NY + 1, z) = phase_field_pairwise(x, 1, z);
						phase_field_Gc(x, MESH_NY + 1, z) = phase_field_Gc(x, 1, z);
					}
			}
			else if (y_up == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, MESH_NY + 1, z) = phase_field->at(x, MESH_NY, z);
						phase_field_pairwise(x, MESH_NY + 1, z) = phase_field_pairwise(x, MESH_NY, z);
						phase_field_Gc(x, MESH_NY + 1, z) = phase_field_Gc(x, MESH_NY, z);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(x, MESH_NY + 1, z);
					}
			}
			if (z_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, MESH_NZ + 1) = phase_field->at(x, y, 1);
						phase_field_pairwise(x, y, MESH_NZ + 1) = phase_field_pairwise(x, y, 1);
						phase_field_Gc(x, y, MESH_NZ + 1) = phase_field_Gc(x, y, 1);
					}
			}
			else if (z_up == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, MESH_NZ + 1) = phase_field->at(x, y, MESH_NZ);
						phase_field_pairwise(x, y, MESH_NZ + 1) = phase_field_pairwise(x, y, MESH_NZ);
						phase_field_Gc(x, y, MESH_NZ + 1) = phase_field_Gc(x, y, MESH_NZ);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						BoundaryCondition_PhiFix(x, y, MESH_NZ + 1);
					}
			}
		}

		void boundary_condition() {
			if (x_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(0, y, z) = phase_field->at(MESH_NX, y, z);
						phase_field_pairwise(0, y, z) = phase_field_pairwise(MESH_NX, y, z);
					}
			}
			else if (x_down == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(0, y, z) = phase_field->at(1, y, z);
						phase_field_pairwise(0, y, z) = phase_field_pairwise(1, y, z);
					}
			}
			else {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(0, y, z);
					}
			}
			if (y_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, 0, z) = phase_field->at(x, MESH_NY, z);
						phase_field_pairwise(x, 0, z) = phase_field_pairwise(x, MESH_NY, z);
					}
			}
			else if (y_down == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, 0, z) = phase_field->at(x, 1, z);
						phase_field_pairwise(x, 0, z) = phase_field_pairwise(x, 1, z);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nz(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(x, 0, z);
					}
			}
			if (z_down == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, 0) = phase_field->at(x, y, MESH_NZ);
						phase_field_pairwise(x, y, 0) = phase_field_pairwise(x, y, MESH_NZ);
					}
			}
			else if (z_down == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, 0) = phase_field->at(x, y, 1);
						phase_field_pairwise(x, y, 0) = phase_field_pairwise(x, y, 1);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						BoundaryCondition_PhiFix(x, y, 0);
					}
			}
			if (x_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(MESH_NX + 1, y, z) = phase_field->at(1, y, z);
						phase_field_pairwise(MESH_NX + 1, y, z) = phase_field_pairwise(1, y, z);
					}
			}
			else if (x_up == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(MESH_NX + 1, y, z) = phase_field->at(MESH_NX, y, z);
						phase_field_pairwise(MESH_NX + 1, y, z) = phase_field_pairwise(MESH_NX, y, z);
					}
			}
			else {
#pragma omp parallel for
				for (int y = 0; y < phase_field->Ny(); y++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(MESH_NX + 1, y, z);
					}
			}
			if (y_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, MESH_NY + 1, z) = phase_field->at(x, 1, z);
						phase_field_pairwise(x, MESH_NY + 1, z) = phase_field_pairwise(x, 1, z);
					}
			}
			else if (y_up == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						phase_field->at(x, MESH_NY + 1, z) = phase_field->at(x, MESH_NY, z);
						phase_field_pairwise(x, MESH_NY + 1, z) = phase_field_pairwise(x, MESH_NY, z);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int z = 0; z < phase_field->Nz(); z++) {
						BoundaryCondition_PhiFix(x, MESH_NY + 1, z);
					}
			}
			if (z_up == BoundaryCondition::PERIODIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, MESH_NZ + 1) = phase_field->at(x, y, 1);
						phase_field_pairwise(x, y, MESH_NZ + 1) = phase_field_pairwise(x, y, 1);
					}
			}
			else if (z_up == BoundaryCondition::ADIABATIC) {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						phase_field->at(x, y, MESH_NZ + 1) = phase_field->at(x, y, MESH_NZ);
						phase_field_pairwise(x, y, MESH_NZ + 1) = phase_field_pairwise(x, y, MESH_NZ);
					}
			}
			else {
#pragma omp parallel for
				for (int x = 0; x < phase_field->Nx(); x++)
					for (int y = 0; y < phase_field->Ny(); y++) {
						BoundaryCondition_PhiFix(x, y, MESH_NZ + 1);
					}
			}
		}
	};
}