#include "GrandPotentialEquation.h"
namespace pf {
	namespace grand_potential_functions {

		vector<double> pre_calculation_grand_potential_functional(double dt) {
			vector<double> MAX_VARIATION; MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0); MAX_VARIATION.push_back(0.0);
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						node.customValues[ExternalFields::CON_Smooth_Phi] = node.cal_phases_fraction_by_index(phase_indexes);
						for (auto comp = node.x.begin(); comp < node.x.end(); comp++)
							comp->value = 0.0;
						for (auto phase = node.begin(); phase < node.end(); phase++) {
							for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++)
								comp->value = 0.0;
							phase_x(node, *phase);
							for (auto comp = phase->x.begin(); comp < phase->x.end(); comp++) {
								node.x[comp->index].value += phase->phi * comp->value;
								phase->potential.add_con(comp->index, node.potential[comp->index].value);
							}
						}
						for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
							p->gradient.set_to_zero();
							p->increment = 0.0;
							p->laplacian = 0.0;
							for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++)
								node.kinetics_coeff.set(p->index, p2->index, M_ij(node, p->index, p2->index));
						}
						init_grand_potential_on_moving_interface(node, ConEquationDomain::CEDomain_Standard, threshold);
					}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
							if (diff_method == DifferenceMethod::FIVE_POINT) {
								pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
								pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
								pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
								pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
								pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
								pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
								if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upx = &node;
								if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downx = &node;
								if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upy = &node;
								if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downy = &node;
								if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upz = &node;
								if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downz = &node;
								for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
									p->gradient[0] = (node_upx->potential[p->index].value -
										node_downx->potential[p->index].value) / 2.0 / phaseMesh->dr;
									p->gradient[1] = (node_upy->potential[p->index].value -
										node_downy->potential[p->index].value) / 2.0 / phaseMesh->dr;
									p->gradient[2] = (node_upz->potential[p->index].value -
										node_downz->potential[p->index].value) / 2.0 / phaseMesh->dr;
									p->laplacian = (node_downx->potential[p->index].value
										+ node_upx->potential[p->index].value
										+ node_downy->potential[p->index].value
										+ node_upy->potential[p->index].value
										+ node_downz->potential[p->index].value
										+ node_upz->potential[p->index].value - 6 * p->value) / phaseMesh->dr / phaseMesh->dr;
								}
								for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
									for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
										node.kinetics_coeff(p->index, p2->index).gradient[grad_x] = (node_upx->kinetics_coeff(p->index, p2->index).value -
											node_downx->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
										node.kinetics_coeff(p->index, p2->index).gradient[grad_y] = (node_upy->kinetics_coeff(p->index, p2->index).value -
											node_downy->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
										node.kinetics_coeff(p->index, p2->index).gradient[grad_z] = (node_upz->kinetics_coeff(p->index, p2->index).value -
											node_downz->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									}
								}
							}
							else if (diff_method == DifferenceMethod::NINE_POINT) {
								pf::PhaseNode* node_upx = &node.get_neighbor_node(Direction::x_up);
								pf::PhaseNode* node_downx = &node.get_neighbor_node(Direction::x_down);
								pf::PhaseNode* node_upy = &node.get_neighbor_node(Direction::y_up);
								pf::PhaseNode* node_downy = &node.get_neighbor_node(Direction::y_down);
								pf::PhaseNode* node_upz = &node.get_neighbor_node(Direction::z_up);
								pf::PhaseNode* node_downz = &node.get_neighbor_node(Direction::z_down);
								if (node_upx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upx = &node;
								if (node_downx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downx = &node;
								if (node_upy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upy = &node;
								if (node_downy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downy = &node;
								if (node_upz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upz = &node;
								if (node_downz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downz = &node;
								pf::PhaseNode* node_upxupy = &node.get_long_range_node(1, 1, 0);
								pf::PhaseNode* node_downxdowny = &node.get_long_range_node(-1, -1, 0);
								pf::PhaseNode* node_upydownx = &node.get_long_range_node(-1, 1, 0);
								pf::PhaseNode* node_downyupx = &node.get_long_range_node(1, -1, 0);
								pf::PhaseNode* node_upxupz = &node.get_long_range_node(1, 0, 1);
								pf::PhaseNode* node_downxdownz = &node.get_long_range_node(-1, 0, -1);
								pf::PhaseNode* node_upzdownx = &node.get_long_range_node(-1, 0, 1);
								pf::PhaseNode* node_downzupx = &node.get_long_range_node(1, 0, -1);
								pf::PhaseNode* node_upzupy = &node.get_long_range_node(0, 1, 1);
								pf::PhaseNode* node_downzdowny = &node.get_long_range_node(0, -1, -1);
								pf::PhaseNode* node_upydownz = &node.get_long_range_node(0, 1, -1);
								pf::PhaseNode* node_downyupz = &node.get_long_range_node(0, -1, 1);
								if (node_upxupy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upxupy = &node;
								if (node_downxdowny->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downxdowny = &node;
								if (node_upydownx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upydownx = &node;
								if (node_downyupx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downyupx = &node;
								if (node_upxupz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upxupz = &node;
								if (node_downxdownz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downxdownz = &node;
								if (node_upzdownx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upzdownx = &node;
								if (node_downzupx->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downzupx = &node;
								if (node_upzupy->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upzupy = &node;
								if (node_downzdowny->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downzdowny = &node;
								if (node_upydownz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_upydownz = &node;
								if (node_downyupz->customValues[ExternalFields::CON_Smooth_Phi] < threshold)
									node_downyupz = &node;

								for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
									p->gradient[0] = (node_upx->potential[p->index].value -
										node_downx->potential[p->index].value) / 2.0 / phaseMesh->dr;
									p->gradient[1] = (node_upy->potential[p->index].value -
										node_downy->potential[p->index].value) / 2.0 / phaseMesh->dr;
									p->gradient[2] = (node_upz->potential[p->index].value -
										node_downz->potential[p->index].value) / 2.0 / phaseMesh->dr;
									p->laplacian = (node_downx->potential[p->index].value
										+ node_upx->potential[p->index].value
										+ node_downy->potential[p->index].value
										+ node_upy->potential[p->index].value
										+ node_downz->potential[p->index].value
										+ node_upz->potential[p->index].value - 6 * p->value) / phaseMesh->dr / phaseMesh->dr;
								}
								for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
									for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
										node.kinetics_coeff(p->index, p2->index).gradient[grad_x] = (node_upx->kinetics_coeff(p->index, p2->index).value -
											node_downx->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
										node.kinetics_coeff(p->index, p2->index).gradient[grad_y] = (node_upy->kinetics_coeff(p->index, p2->index).value -
											node_downy->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
										node.kinetics_coeff(p->index, p2->index).gradient[grad_z] = (node_upz->kinetics_coeff(p->index, p2->index).value -
											node_downz->kinetics_coeff(p->index, p2->index).value) / 2.0 / phaseMesh->dr;
									}
								}
							}
						}
					}
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
							Vector3 s_phi_grad = node.cal_customValues_gradient(ExternalFields::CON_Smooth_Phi, phaseMesh->dr);
							double s_phi_grad_norm = s_phi_grad.abs();
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								double DIFF_FLUX = 0.0, REAC_FLUX = 0.0, PHASETRANS_FLUX = 0.0;
								// diffusion term
								for (auto p2 = node.potential.begin(); p2 < node.potential.end(); p2++) {
									Vector3 D_KineticCoef = node.kinetics_coeff.get_gradientVec3(p->index, p2->index);
									double kinetic_coef = node.kinetics_coeff(p->index, p2->index).value;

									DIFF_FLUX += kinetic_coef * (s_phi_grad * p2->gradient)
										+ (D_KineticCoef * p2->gradient + kinetic_coef * p2->laplacian) * node.customValues[ExternalFields::CON_Smooth_Phi];

								}
#ifdef _DEBUG
								if (_isnan(DIFF_FLUX)) {
									cout << "DEBUG: DIFF_FLUX error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif
								// reaction term
								if (s_phi_grad_norm > SYS_EPSILON)
									REAC_FLUX = s_phi_grad_norm * int_flux(node, p->index);
								REAC_FLUX += Source(node, p->index) * node.customValues[ExternalFields::CON_Smooth_Phi];

								// phase transtion term
								for (auto phase = node.begin(); phase < node.end(); phase++)
									PHASETRANS_FLUX -= (phase->phi - phase->old_phi) / dt * phase->x[p->index].value;
#ifdef _DEBUG
								if (_isnan(PHASETRANS_FLUX)) {
									cout << "DEBUG: PHASETRANS_FLUX error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif
								// summary
								p->increment = DIFF_FLUX + REAC_FLUX + PHASETRANS_FLUX;
#ifdef _DEBUG
								if (_isnan(p->increment)) {
									cout << "DEBUG: node.potential[x].increment error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif
								double pre_factor_con_i = 0.0;
								for (auto phase = node.begin(); phase < node.end(); phase++)
									pre_factor_con_i += phase->phi * dphase_x_du(node, *phase, p->index);
#ifdef _DEBUG
								if (Is_Equality(pre_factor_con_i, 0.0)) {
									cout << "DEBUG: pre_factor_con_i error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif

								p->increment /= pre_factor_con_i;

#ifdef _OPENMP
#pragma omp critical
#endif
								{
									if (abs(DIFF_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX]) // MAX_diffusionFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_DIFFUSION_FLUX] = abs(DIFF_FLUX * dt);
									if (abs(REAC_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX]) // MAX_reactionFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_REACTION_FLUX] = abs(REAC_FLUX * dt);
									if (abs(PHASETRANS_FLUX * dt) > MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX]) // MAX_phaseTransFlux
										MAX_VARIATION[CH_RETURN::CH_MAX_PHI_TRANS_FLUX] = abs(PHASETRANS_FLUX * dt);
								}

							}

						}
					}
			return MAX_VARIATION;
		}

		double solve_grand_potential_functional(double dt) {
			double MAX_VARIATION = 0.0;
#pragma omp parallel for
			for (int x = 0; x < phaseMesh->limit_x; x++)
				for (int y = 0; y < phaseMesh->limit_y; y++)
					for (int z = 0; z < phaseMesh->limit_z; z++) {
						PhaseNode& node = (*phaseMesh)(x, y, z);
						if (node.customValues[ExternalFields::CON_Smooth_Phi] > threshold) {
							for (auto p = node.potential.begin(); p < node.potential.end(); p++) {
								double old_p = p->value;
								p->value += p->increment * dt;
								Boundary_Condition_TotalX(node, p->index);
#ifdef _OPENMP
#pragma omp critical
#endif
								{
									if (abs(old_p - p->value) > MAX_VARIATION) // MAX_x_increment
										MAX_VARIATION = abs(old_p - p->value);
								}
#ifdef _DEBUG
								if (_isnan(p->value)) {
									cout << "DEBUG: node.potential[x].value error !" << endl;
									SYS_PROGRAM_STOP;
								}
#endif
							}
						}
						node.customValues[ExternalFields::CON_Smooth_Old_Phi] = node.customValues[ExternalFields::CON_Smooth_Phi];
					}
			return MAX_VARIATION;
		}

	}
}