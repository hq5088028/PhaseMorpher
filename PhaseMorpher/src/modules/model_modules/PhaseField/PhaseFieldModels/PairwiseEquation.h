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

	};

}