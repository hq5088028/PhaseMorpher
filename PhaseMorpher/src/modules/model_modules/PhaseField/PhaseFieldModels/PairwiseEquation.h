#pragma once
#include <omp.h>
#include "../../../data_struct/Mesh_0.h"
#include "../../../../tools/mathTools.h"
#include "../PhaseFieldPoints.h"
#include "../GrainsOrientations.h"
namespace pf {
	using namespace std;
	// - 
	namespace pair_wise_parameters {
		// - 
		inline REAL* PHI_MAX_VARIATION;
		inline REAL* delt_t;
		inline REAL* delt_r;
		inline Mesh_Boundry<PhaseFieldPoint>* phase_field = nullptr;
		inline int (*phi_name_property)(std::string phi_name);
		// main mesh boundary condition
		inline int MESH_NX = 1;
		inline int MESH_NY = 1;
		inline int MESH_NZ = 1;
		inline BoundaryCondition x_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition y_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition z_down = BoundaryCondition::PERIODIC;
		inline BoundaryCondition x_up = BoundaryCondition::PERIODIC;
		inline BoundaryCondition y_up = BoundaryCondition::PERIODIC;
		inline BoundaryCondition z_up = BoundaryCondition::PERIODIC;
		// - 
		inline int phi_property_number = 0;
		inline int phi_number = 0;
		// - 
		inline vector<int> phi_property;
		// Statement
		enum Int_Gradient { Steinbach_1996, Steinbach_1999, Steinbach_G2009, DanielCrack_G2016 };
		enum Int_Potential { Nestler_Well, Nestler_Obstacle, Steinbach_P2009, DanielCrack_P2016 };
		enum DifferenceMethod { FIVE_POINT, NINE_POINT };
		enum Flag { pf_BULK, pf_NEAR_INTERFACE, pf_INTERFACE }; // where is the points
		// pair-wise field data,
		inline Mesh_Boundry<PairwisePoint> phase_field_pairwise;
		inline Mesh_Boundry<REAL> phase_field_Gc;
		// - accelerated
		inline int phi_acc_number = 9;
		// - normalization
		inline bool is_phi_normalized = false;
		// - 
		inline DifferenceMethod diff_method = DifferenceMethod::FIVE_POINT;
		inline GrainsOrientations grains_orientation;
		// - driving force functions
		inline vector<REAL(*)(pf::PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index)> delt_Fbulk_delt_phi;
		// - source functions
		inline vector<REAL(*)(pf::PhaseFieldPoint& point, PairwisePoint& point_ex, int alpha_index, int beta_index)> source_alpha_beta;
		// - boundary condition funtions for fix value
		inline REAL(*BoundaryCondition_PhiFix)(int x, int y, int z);
		// - parameters
		// - anisotropy
		// - mobility
		inline REAL const_Lij = REAL(0.0);
		inline vector<vector<int>> const_block_Lij; // <- phi index
		inline vector<vector<REAL>> matirx_Lij; // <- phi index
		inline vector<vector<int>> block_Lij; // <- phi index
		inline vector<REAL> block_value_Lij; // <- phi index
		// - interface mobility anisotropy
		enum Int_Mobility_Anisotropic { IMA_ISO, IMA_CUBIC, IMA_HEX_BOETTGER, IMA_HEX_SUN, IMA_HEX_YANG };
		inline vector<vector<int>>  intMob_anisotropic_model_matrix; // <- phi property
		inline vector<vector<REAL>> intMobAniso1_matrix; // <- phi property
		inline vector<vector<REAL>> intMobAniso2_matrix; // <- phi property
		inline vector<vector<REAL>> intMobAniso3_matrix; // <- phi property
		inline vector<vector<REAL>> intMobAniso4_matrix; // <- phi property
		// - interface energy models
		inline int interface_gradient = Int_Gradient::Steinbach_G2009;
		inline int interface_potential = Int_Potential::Steinbach_P2009;
		// interface energy
		inline REAL interface_width = REAL(4.0);
		inline REAL const_xi_ab = REAL(0.0);
		inline REAL const_xi_abc = REAL(0.0);
		inline vector<vector<REAL>> matirx_xi_ab; // <- phi property
		inline vector<vector<vector<REAL>>> matirx_xi_abc; // <- phi property
		// interface energy anisotropy
		enum Int_Energy_Anisotropic { IEA_ISO, IEA_CUBIC, IEA_HEX_BOETTGER, IEA_HEX_SUN, IEA_HEX_YANG };
		inline vector<vector<int>>  intEn_anisotropic_model_matrix; // <- phi property
		inline vector<vector<REAL>> intEnAniso1_matrix; // <- phi property
		inline vector<vector<REAL>> intEnAniso2_matrix; // <- phi property
		inline vector<vector<REAL>> intEnAniso3_matrix; // <- phi property
		inline vector<vector<REAL>> intEnAniso4_matrix; // <- phi property
		// multigrains crack
		const REAL crack_K = REAL(9.0 / 64.0);
		inline int crack_phase_property = -1;
		inline vector<REAL> phase_Gc;   // <- phi property
		inline REAL crack_dr = REAL(1.0);
		inline REAL crack_int_width = REAL(1.0);
		inline REAL crack_num_cut_off = REAL(0.99);
	}
	namespace pair_wise_functions {
		// -
		pair_wise_parameters::Flag currentFlag(int x, int y, int z, PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		// - 
		Vector3 normals(REAL alpha_phi, REAL beta_phi, Vector3 alpha_grad, Vector3 beta_grad);
		// - mobility
		REAL Lij_pair_wise_const_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_pair_wise_anisotropic_const_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_pair_wise_const_cross_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_pair_wise_anisotropic_const_cross_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_pair_wise_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_pair_wise_anisotropic_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_anisotropic_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL Lij_anisotropic_block(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		// - interface energy
		REAL xi_ab_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL xi_ab_anisotropic_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL xi_abc_const(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index, int phi_gamma_index);
		REAL xi_ab_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL xi_ab_anisotropic_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL xi_abc_matrix(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index, int phi_gamma_index);
		// - 
		inline REAL(*mobility)(pf::PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		inline REAL(*_xi_ab)(pf::PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		inline REAL(*_xi_abc)(pf::PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index, int phi_gamma_index);
		// - crack resistence
		REAL Gc_acc(PhaseFieldPoint& point, PairwisePoint& point_ex);
		REAL dGc_dphi_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		// - interface energy model
		REAL dfint_dphi_grad_S1996_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		REAL dfint_dphi_grad_S1999_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		REAL dfint_dphi_pot_Nobstacle_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		REAL dfint_dphi_pot_Nwell_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		REAL dfint_dphi_solid_D2016_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		REAL dfint_dphi_crack_D2016_acc(int x, int y, int z, PhaseFieldPoint& point, int phi_index);
		void pairwise_normalize_daniel_crack_acc(PhaseFieldPoint& point, PairwisePoint& point_ex);
		REAL dfint_dphi_pairwise_acc(int x, int y, int z, PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		void pairwise_normalize_normal(pf::PhaseFieldPoint& point, pf::PairwisePoint& point_ex);
		REAL dfint_dphi_S2009(PhaseFieldPoint& point, int phi_alpha_index, int phi_beta_index);
		REAL dfbulk_dphi_pairwise_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_index);
		REAL source_pairwise_acc(PhaseFieldPoint& point, PairwisePoint& point_ex, int phi_alpha_index, int phi_beta_index);

		void pre_calculation_phi_pair_wise();

		REAL solve_phi_pair_wise();

		void boundary_condition();

		void boundary_condition_crack();

		// prepare pair-wise field variables
		inline void exec_pre_iii() {
#pragma omp parallel for
			for (int x = pair_wise_parameters::phase_field->X_BEGIN(); x <= pair_wise_parameters::phase_field->X_END(); x++)
				for (int y = pair_wise_parameters::phase_field->Y_BEGIN(); y <= pair_wise_parameters::phase_field->Y_END(); y++)
					for (int z = pair_wise_parameters::phase_field->Z_BEGIN(); z <= pair_wise_parameters::phase_field->Z_END(); z++) {
						PhaseFieldPoint& point = pair_wise_parameters::phase_field->at(x, y, z);
						PairwisePoint& point_ex = pair_wise_parameters::phase_field_pairwise(x, y, z);
						if (pair_wise_parameters::is_phi_normalized) {
							REAL sum_phi = 0.0;
							for (int index = 0; index < pair_wise_parameters::phi_number; index++)
								sum_phi += point.phi[index];
							for (int index = 0; index < pair_wise_parameters::phi_number; index++)
								point.phi[index] /= sum_phi;
						}
						for (int index = 0; index < pair_wise_parameters::phi_number; index++) {
							point.old_phi[index] = point.phi[index];
							point_ex.flag[index] = currentFlag(x, y, z, point, point_ex, index);
						}
					}
		}

		// main calculation
		inline void exec_i() {
			pre_calculation_phi_pair_wise();
			*pair_wise_parameters::PHI_MAX_VARIATION = solve_phi_pair_wise();
			if (pair_wise_parameters::interface_gradient == pair_wise_parameters::Int_Gradient::DanielCrack_G2016
				|| pair_wise_parameters::interface_potential == pair_wise_parameters::Int_Potential::DanielCrack_P2016)
				boundary_condition_crack();
			else
				boundary_condition();
		}

		// release field and parameters
		inline void deinit() {
			for (int x = 0; x < pair_wise_parameters::phase_field_pairwise.Nx(); x++)
				for (int y = 0; y < pair_wise_parameters::phase_field_pairwise.Ny(); y++)
					for (int z = 0; z < pair_wise_parameters::phase_field_pairwise.Nz(); z++) {
						pair_wise_parameters::phase_field_pairwise(x, y, z).clear();
					}
			pair_wise_parameters::phase_field_pairwise.clear();
			pair_wise_parameters::phase_field_Gc.clear();
		}

		//// init field and parameters
		//inline void init_pairwise_equation(Mesh_Boundry<PhaseFieldPoint>& _phase_field, REAL& _PHI_MAX_VARIATION, REAL& _delt_t, REAL& _delt_r,
		//	int (*_phi_name_property)(std::string phi_name), int _phi_number, int _phi_property_number, vector<int> _phi_property, 
		//	int _MESH_NX, int _MESH_NY, int _MESH_NZ, BoundaryCondition _x_down, BoundaryCondition _y_down, BoundaryCondition _z_down, 
		//	BoundaryCondition _x_up, BoundaryCondition _y_up, BoundaryCondition _z_up) {
		//	// - 
		//	pair_wise_parameters::MESH_NX = _MESH_NX;
		//	pair_wise_parameters::MESH_NY = _MESH_NY;
		//	pair_wise_parameters::MESH_NZ = _MESH_NZ;
		//	pair_wise_parameters::x_down = _x_down;
		//	pair_wise_parameters::y_down = _y_down;
		//	pair_wise_parameters::z_down = _z_down;
		//	pair_wise_parameters::x_up = _x_up;
		//	pair_wise_parameters::y_up = _y_up;
		//	pair_wise_parameters::z_up = _z_up;
		//	// - 
		//	pair_wise_parameters::phase_field = &_phase_field;
		//	pair_wise_parameters::PHI_MAX_VARIATION = &_PHI_MAX_VARIATION;
		//	pair_wise_parameters::delt_t = &_delt_t;
		//	pair_wise_parameters::delt_r = &_delt_r;
		//	pair_wise_parameters::phi_name_property = _phi_name_property;   // func ptr
		//	pair_wise_parameters::phi_number = _phi_number;
		//	pair_wise_parameters::phi_property_number = _phi_property_number;
		//	pair_wise_parameters::phi_property = _phi_property;
		//	// - 
		//	InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.Acceleration.ContainerSize", pair_wise_parameters::phi_acc_number, true);
		//	pair_wise_parameters::phase_field_pairwise.init(pair_wise_parameters::MESH_NX, pair_wise_parameters::MESH_NY, pair_wise_parameters::MESH_NZ);
		//	pair_wise_parameters::phase_field_Gc.init(pair_wise_parameters::MESH_NX, pair_wise_parameters::MESH_NY, pair_wise_parameters::MESH_NZ);
		//	for (int x = 0; x < pair_wise_parameters::phase_field_pairwise.Nx(); x++)
		//		for (int y = 0; y < pair_wise_parameters::phase_field_pairwise.Ny(); y++)
		//			for (int z = 0; z < pair_wise_parameters::phase_field_pairwise.Nz(); z++) {
		//				pair_wise_parameters::phase_field_pairwise(x, y, z).init(pair_wise_parameters::phi_acc_number, pair_wise_parameters::phi_number);
		//			}
		//	InputFileReader::get_instance()->read_bool_value("Solver.ControlEquation.Pairwise.NormalizePhi", pair_wise_parameters::is_phi_normalized, true);

		//	int difference_method = pair_wise_parameters::DifferenceMethod::FIVE_POINT;
		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.DifferenceMethod : 0 - FIVE_POINT , 1 - NINE_POINT \n");
		//	InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.DifferenceMethod", difference_method, true);
		//	pair_wise_parameters::diff_method = pair_wise_parameters::DifferenceMethod(difference_method);

		//	WriteDebugFile("# .Pairwise.GrainsOrientations.rotation_gauge = 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n");
		//	WriteDebugFile("#                                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n");
		//	int rotation_gauge = RotationGauge::RG_ZXZ;
		//	InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.GrainsOrientations.rotation_gauge", rotation_gauge, true);
		//	pair_wise_parameters::grains_orientation.init(RotationGauge(rotation_gauge), pair_wise_parameters::phi_number);

		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.GrainsOrientations = {[(phi_index_0, phi_index_2, ... ),(rotation_angle_1, rotation_angle_2, rotation_angle_3)],  ... } \n");
		//	string grains_key = "Solver.ControlEquation.Pairwise.GrainsOrientations", grains_input = "{}";
		//	InputFileReader::get_instance()->read_string_value(grains_key, grains_input, true);
		//	vector<InputValueType> grains_structure; grains_structure.push_back(InputValueType::IVType_INT); grains_structure.push_back(InputValueType::IVType_REAL);
		//	vector<vector<vector<input_value>>> grains_value;
		//	grains_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(grains_structure, grains_key, grains_input, true);
		//	for (int i = 0; i < grains_value.size(); i++) {
		//		for (int j = 0; j < grains_value[i][0].size(); j++)
		//			pair_wise_parameters::grains_orientation.set_phi_orientation(grains_value[i][0][j].int_value,
		//				Vector3(AngleToRadians(grains_value[i][1][0].REAL_value), AngleToRadians(grains_value[i][1][1].REAL_value), AngleToRadians(grains_value[i][1][2].REAL_value)));
		//	}

		//	// - interface mobility
		//	mobility = Lij_pair_wise_const;
		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.Lij.const  = Lij_value \n");
		//	WriteDebugFile("#                                .matrix = [(phi_i, phi_j, Lij_value), ... ] \n");
		//	WriteDebugFile("#                                .block = [(phi_begin, phi_end, Lij_value), ... ] \n");
		//	int intMob_type = 0;
		//	string matrix_key = "Solver.ControlEquation.Pairwise.Lij.matrix", matrix_input = "[()]", block_key = "Solver.ControlEquation.Pairwise.Lij.block", block_input = "[()]";
		//	if (InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.Lij.const", pair_wise_parameters::const_Lij, true)) {
		//		mobility = Lij_pair_wise_const;
		//		intMob_type = 0;
		//		WriteDebugFile("# Solver.ControlEquation.Pairwise.Lij.Const.block       = [(phi_index1, phi_index2, ... ), ... ] \n");
		//		WriteDebugFile("#                                    .cross_block = [(phi_alpha1, phi_alpha2, ... ), (phi_beta1, phi_beta2, ... )] \n");
		//		string const_block_key = "Solver.ControlEquation.Pairwise.Lij.Const.block", const_cross_block_key = "Solver.ControlEquation.Pairwise.Lij.Const.cross_block", const_block_input = "[()]";
		//		if (InputFileReader::get_instance()->read_string_value(const_block_key, const_block_input, true)) {
		//			mobility = Lij_pair_wise_const_block;
		//			intMob_type = 1;
		//			vector<vector<input_value>> const_block_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
		//			for (int index = 0; index < const_block_value.size(); index++) {
		//				vector<int> const_block;
		//				for (int phi_index = 0; phi_index < const_block_value[index].size(); phi_index++)
		//					const_block.push_back(const_block_value[index][phi_index].int_value);
		//				pair_wise_parameters::const_block_Lij.push_back(const_block);
		//			}
		//		}
		//		else if (InputFileReader::get_instance()->read_string_value(const_cross_block_key, const_block_input, true)) {
		//			mobility = Lij_pair_wise_const_cross_block;
		//			intMob_type = 4;
		//			vector<vector<input_value>> const_block_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
		//			for (int index = 0; index < 2; index++) {
		//				vector<int> const_block;
		//				for (int phi_index = 0; phi_index < const_block_value[index].size(); phi_index++)
		//					const_block.push_back(const_block_value[index][phi_index].int_value);
		//				pair_wise_parameters::const_block_Lij.push_back(const_block);
		//			}
		//		}
		//	}
		//	else if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_input, true)) {
		//		mobility = Lij_matrix;
		//		pair_wise_parameters::matirx_Lij.resize(pair_wise_parameters::phi_number);
		//		for (int index = 0; index < pair_wise_parameters::phi_number; index++)
		//			pair_wise_parameters::matirx_Lij[index].resize(pair_wise_parameters::phi_number);
		//		intMob_type = 2;
		//		vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_REAL);
		//		vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key, matrix_input, true);
		//		for (int index = 0; index < matrix_value.size(); index++)
		//			pair_wise_parameters::matirx_Lij[matrix_value[index][0].int_value][matrix_value[index][1].int_value] = matrix_value[index][2].REAL_value;
		//	}
		//	else if (InputFileReader::get_instance()->read_string_value(block_key, block_input, true)) {
		//		mobility = Lij_block;
		//		intMob_type = 3;
		//		vector<InputValueType> block_structure; block_structure.push_back(InputValueType::IVType_INT); 
		//		block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_REAL);
		//		vector<vector<input_value>> block_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(block_structure, block_key, block_input, true);
		//		for (int index = 0; index < block_value.size(); index++) {
		//			pair_wise_parameters::block_Lij.push_back({ block_value[index][0].int_value, block_value[index][1].int_value });
		//			pair_wise_parameters::block_value_Lij.push_back(block_value[index][2].REAL_value);
		//		}
		//	}
		//	// - anisotropic
		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.Lij.Anisotropic.matrix = [(phi_a_name, phi_b_name, AnisotropicModelIndex ), ...]\n");
		//	WriteDebugFile("#                                     AnisotropicModelIndex : 0 - IEA_ISO, 1 - IEA_CUBIC, 2 - IEA_HEX_BOETTGER, 3 - IEA_HEX_SUN, 4 - IEA_HEX_YANG \n");
		//	string matrix_aniso_key = "Solver.ControlEquation.Pairwise.Anisotropic.matrix", matrix_aniso_input = "[()]";
		//	if (InputFileReader::get_instance()->read_string_value(matrix_aniso_key, matrix_aniso_input, true)) {
		//		if (intMob_type == 0)
		//			mobility = Lij_pair_wise_anisotropic_const;
		//		else if (intMob_type == 1)
		//			mobility = Lij_pair_wise_anisotropic_const_block;
		//		else if (intMob_type == 2)
		//			mobility = Lij_anisotropic_matrix;
		//		else if (intMob_type == 3)
		//			mobility = Lij_anisotropic_block;
		//		else if (intMob_type == 4)
		//			mobility = Lij_pair_wise_anisotropic_const_cross_block;
		//		pair_wise_parameters::intMob_anisotropic_model_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intMobAniso1_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intMobAniso2_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intMobAniso3_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intMobAniso4_matrix.resize(pair_wise_parameters::phi_property_number);
		//		for (int index = 0; index < pair_wise_parameters::phi_property_number; index++) {
		//			pair_wise_parameters::intMob_anisotropic_model_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intMobAniso1_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intMobAniso2_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intMobAniso3_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intMobAniso4_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//		}
		//		vector<InputValueType> matrix_aniso_structure; matrix_aniso_structure.push_back(InputValueType::IVType_STRING); matrix_aniso_structure.push_back(InputValueType::IVType_STRING); matrix_aniso_structure.push_back(InputValueType::IVType_INT);
		//		vector<vector<input_value>> matrix_aniso_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_aniso_structure, matrix_aniso_key, matrix_aniso_input, true);
		//		for (int index = 0; index < matrix_aniso_value.size(); index++) {
		//			int alpha_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][0].string_value),
		//				beta_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][1].string_value);
		//			pair_wise_parameters::intMob_anisotropic_model_matrix[alpha_property][beta_property] = matrix_aniso_value[index][2].int_value;
		//			pair_wise_parameters::intMob_anisotropic_model_matrix[beta_property][alpha_property] = matrix_aniso_value[index][2].int_value;
		//			if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_CUBIC) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//			else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//			else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_HEX_SUN) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
		//				string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
		//				string matrix_aniso_par4 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_4";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
		//				pair_wise_parameters::intMobAniso2_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso2_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
		//				pair_wise_parameters::intMobAniso3_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso3_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par4, aniso_par, true);
		//				pair_wise_parameters::intMobAniso4_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso4_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//			else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_HEX_YANG) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
		//				string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
		//				pair_wise_parameters::intMobAniso2_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso2_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
		//				pair_wise_parameters::intMobAniso3_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intMobAniso3_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//		}
		//	}

		//	// - interface energy
		//	_xi_ab = xi_ab_const;
		//	_xi_abc = xi_abc_const;
		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.InterfaceEnergy.int_gradient : 0 - Steinbach_1996 , 1 - Steinbach_1999 , 2 - Steinbach_G2009 , 3 - DanielCrack_G2016\n");
		//	InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.InterfaceEnergy.int_gradient", pair_wise_parameters::interface_gradient, true);
		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.InterfaceEnergy.int_potential : 0 - Nestler_Well , 1 - Nestler_Obstacle , 2 - Steinbach_P2009 , 3 - DanielCrack_P2016\n");
		//	InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.InterfaceEnergy.int_potential", pair_wise_parameters::interface_potential, true);
		//	InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.interface_width", pair_wise_parameters::interface_width, true);

		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.xi_ab.const  = xi_ab \n");
		//	WriteDebugFile("#                                      .matrix = [(phi_a_name, phi_b_name, xi_ab_value), ...] \n");
		//	
		//	int int_energy_type = 0;
		//	string matrix_key1 = "Solver.ControlEquation.Pairwise.xi_ab.matrix", matrix_input1 = "[()]";
		//	if (InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.xi_ab.const", pair_wise_parameters::const_xi_ab, true)) {
		//		_xi_ab = xi_ab_const;
		//		int_energy_type = 0;
		//	}
		//	else if (InputFileReader::get_instance()->read_string_value(matrix_key1, matrix_input1, true)) {
		//		_xi_ab = xi_ab_matrix;
		//		pair_wise_parameters::matirx_xi_ab.resize(pair_wise_parameters::phi_property_number);
		//		for (int index = 0; index < pair_wise_parameters::phi_property_number; index++)
		//			pair_wise_parameters::matirx_xi_ab[index].resize(pair_wise_parameters::phi_property_number);
		//		int_energy_type = 1;
		//		vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_REAL);
		//		vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key1, matrix_input1, true);
		//		for (int index = 0; index < matrix_value.size(); index++) {
		//			int alpha_property = pair_wise_parameters::phi_name_property(matrix_value[index][0].string_value),
		//				beta_property = pair_wise_parameters::phi_name_property(matrix_value[index][1].string_value);
		//			pair_wise_parameters::matirx_xi_ab[alpha_property][beta_property] = matrix_value[index][2].REAL_value;
		//		}
		//	}

		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.xi_ab.Anisotropic.matrix = [(phi_a_name, phi_b_name, AnisotropicModelIndex ), ...]\n");
		//	WriteDebugFile("#                                       AnisotropicModelIndex : 0 - IEA_ISO, 1 - IEA_CUBIC, 2 - IEA_HEX_BOETTGER, 3 - IEA_HEX_SUN, 4 - IEA_HEX_YANG \n");
		//	string matrix_aniso_key2 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic.matrix", matrix_aniso_input2 = "[()]";
		//	if (InputFileReader::get_instance()->read_string_value(matrix_aniso_key2, matrix_aniso_input2, true)) {
		//		if (int_energy_type == 0)
		//			_xi_ab = xi_ab_anisotropic_const;
		//		else if (int_energy_type == 1)
		//			_xi_ab = xi_ab_anisotropic_matrix;
		//		pair_wise_parameters::intEn_anisotropic_model_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intEnAniso1_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intEnAniso2_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intEnAniso3_matrix.resize(pair_wise_parameters::phi_property_number);
		//		pair_wise_parameters::intEnAniso4_matrix.resize(pair_wise_parameters::phi_property_number);
		//		for (int index = 0; index < pair_wise_parameters::phi_property_number; index++) {
		//			pair_wise_parameters::intEn_anisotropic_model_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intEnAniso1_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intEnAniso2_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intEnAniso3_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//			pair_wise_parameters::intEnAniso4_matrix[index].resize(pair_wise_parameters::phi_property_number);
		//		}
		//		vector<InputValueType> matrix_aniso_structure; matrix_aniso_structure.push_back(InputValueType::IVType_STRING); 
		//		matrix_aniso_structure.push_back(InputValueType::IVType_STRING); matrix_aniso_structure.push_back(InputValueType::IVType_INT);
		//		vector<vector<input_value>> matrix_aniso_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_aniso_structure, matrix_aniso_key, matrix_aniso_input, true);
		//		for (int index = 0; index < matrix_aniso_value.size(); index++) {
		//			int alpha_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][0].string_value),
		//				beta_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][1].string_value);
		//			pair_wise_parameters::intEn_anisotropic_model_matrix[alpha_property][beta_property] = matrix_aniso_value[index][2].int_value;
		//			pair_wise_parameters::intEn_anisotropic_model_matrix[beta_property][alpha_property] = matrix_aniso_value[index][2].int_value;
		//			if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_CUBIC) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//			else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_HEX_BOETTGER) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//			else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_HEX_SUN) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
		//				string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
		//				string matrix_aniso_par4 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_4";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
		//				pair_wise_parameters::intEnAniso2_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso2_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
		//				pair_wise_parameters::intEnAniso3_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso3_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par4, aniso_par, true);
		//				pair_wise_parameters::intEnAniso4_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso4_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//			else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_HEX_YANG) {
		//				string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
		//				string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
		//				string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
		//				REAL aniso_par = 0.0;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
		//				pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
		//				pair_wise_parameters::intEnAniso2_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso2_matrix[beta_property][alpha_property] = aniso_par;
		//				InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
		//				pair_wise_parameters::intEnAniso3_matrix[alpha_property][beta_property] = aniso_par;
		//				pair_wise_parameters::intEnAniso3_matrix[beta_property][alpha_property] = aniso_par;
		//			}
		//		}
		//	}

		//	WriteDebugFile("# Solver.ControlEquation.Pairwise.xi_abc.const  = xi_ab \n");
		//	WriteDebugFile("#                                       .matrix = [(phi_a_name, phi_b_name, phi_c_name, xi_abc_value), ...] \n");
		//	
		//	string matrix_key2 = "Solver.ControlEquation.Pairwise.xi_abc.matrix", matrix_input2 = "[()]";
		//	if (InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.xi_abc.const", pair_wise_parameters::const_xi_abc, true)) {
		//		_xi_abc = xi_abc_const;
		//	}
		//	else if (InputFileReader::get_instance()->read_string_value(matrix_key2, matrix_input2, true)) {
		//		_xi_abc = xi_abc_matrix;
		//		pair_wise_parameters::matirx_xi_abc.resize(pair_wise_parameters::phi_property_number);
		//		for (int index1 = 0; index1 < pair_wise_parameters::phi_property_number; index1++) {
		//			pair_wise_parameters::matirx_xi_abc[index1].resize(pair_wise_parameters::phi_property_number);
		//			for (int index2 = 0; index2 < pair_wise_parameters::phi_property_number; index2++)
		//				pair_wise_parameters::matirx_xi_abc[index1][index2].resize(pair_wise_parameters::phi_property_number);
		//		}
		//		vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_STRING)
		//			; matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_REAL);
		//		vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key2, matrix_input2, true);
		//		for (int index = 0; index < matrix_value.size(); index++) {
		//			int alpha_property = pair_wise_parameters::phi_name_property(matrix_value[index][0].string_value),
		//				beta_property = pair_wise_parameters::phi_name_property(matrix_value[index][1].string_value),
		//				gamma_property = pair_wise_parameters::phi_name_property(matrix_value[index][2].string_value);
		//			pair_wise_parameters::matirx_xi_abc[alpha_property][beta_property][gamma_property] = matrix_value[index][3].REAL_value;
		//			pair_wise_parameters::matirx_xi_abc[alpha_property][gamma_property][beta_property] = matrix_value[index][3].REAL_value;
		//			pair_wise_parameters::matirx_xi_abc[beta_property][alpha_property][gamma_property] = matrix_value[index][3].REAL_value;
		//			pair_wise_parameters::matirx_xi_abc[beta_property][gamma_property][alpha_property] = matrix_value[index][3].REAL_value;
		//			pair_wise_parameters::matirx_xi_abc[gamma_property][alpha_property][beta_property] = matrix_value[index][3].REAL_value;
		//			pair_wise_parameters::matirx_xi_abc[gamma_property][beta_property][alpha_property] = matrix_value[index][3].REAL_value;
		//		}
		//	}

		//	// - multiple crack
		//	if (pair_wise_parameters::interface_gradient == pair_wise_parameters::Int_Gradient::DanielCrack_G2016 || pair_wise_parameters::interface_potential == pair_wise_parameters::Int_Potential::DanielCrack_P2016) {

		//		pair_wise_parameters::is_phi_normalized = true;
		//		pair_wise_parameters::phase_Gc.resize(pair_wise_parameters::phi_property_number);

		//		string crack_name = "";
		//		if (InputFileReader::get_instance()->read_string_value("Solver.ControlEquation.Pairwise.crack_phi_name", crack_name, true))
		//			pair_wise_parameters::crack_phase_property = pair_wise_parameters::phi_name_property(crack_name);

		//		InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.crack_int_width", pair_wise_parameters::crack_int_width, true);

		//		InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.crack_phi_cut_off", pair_wise_parameters::crack_num_cut_off, true);

		//		WriteDebugFile("# Solver.ControlEquation.Pairwise.Gc = [(phi_name, Gc_value), ...] \n");
		//		string gc_key = "Solver.ControlEquation.Pairwise.Gc", gc_input = "[()]";
		//		if (InputFileReader::get_instance()->read_string_value(gc_key, gc_input, true)) {
		//			vector<vector<input_value>> gc_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_STRING,InputValueType::IVType_REAL }, gc_key, gc_input, true);
		//			for (int index = 0; index < gc_value.size(); index++)
		//				pair_wise_parameters::phase_Gc[pair_wise_parameters::phi_name_property(gc_value[index][0].string_value)] = gc_value[index][1].REAL_value;
		//		}
		//	}
		//	// - load this module
		//	load_a_new_module(default_module_function, default_module_function, exec_pre_iii,
		//		exec_i, default_module_function, default_module_function,
		//		default_module_function, default_module_function, default_module_function,
		//		deinit);
		//}

	};

}