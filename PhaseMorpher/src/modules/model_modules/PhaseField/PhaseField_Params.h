#pragma once
#include "../Model_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
namespace pf {
	enum PhaseFieldEquation {
		PFE_Const = 0, PFE_AllenCahn_Pairwise
	};
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
		vector<Vector3> orientations;
	};
	struct PhaseFieldPoint
	{
		std::vector<REAL> phi;
		std::vector<REAL> old_phi;
		std::vector<REAL> lap_phi;
		std::vector<Vector3> grad_phi;
		void init(int phi_number) {
			phi.resize(phi_number);
			old_phi.resize(phi_number);
			lap_phi.resize(phi_number);
			grad_phi.resize(phi_number);
			for (int index = 0; index < phi_number; index++) {
				phi[index] = 0.0;
				old_phi[index] = 0.0;
				lap_phi[index] = 0.0;
				grad_phi[index].set_to_zero();
			}
		}
		void normalize_phi() {
			REAL sum = 0.0;
			for (int index = 0; index < phi.size(); index++) {
				REAL& PHI = phi[index];
				if (PHI > 1.0)
					PHI = 1.0;
				else if (PHI < 0.0)
					PHI = 0.0;
				sum += PHI;
			}
			for (int index = 0; index < phi.size(); index++)
				phi[index] /= sum;
		}
		void clear() {
			phi.clear();
			std::vector<REAL>().swap(phi);
			old_phi.clear();
			std::vector<REAL>().swap(old_phi);
			lap_phi.clear();
			std::vector<REAL>().swap(lap_phi);
			grad_phi.clear();
			std::vector<Vector3>().swap(grad_phi);
		}
		~PhaseFieldPoint() {
			clear();
		}
		PhaseFieldPoint& operator=(const PhaseFieldPoint& n) {
			phi = n.phi;
			old_phi = n.old_phi;
			lap_phi = n.lap_phi;
			grad_phi = n.grad_phi;
			return *this;
		}
	};
	namespace phi_parameters {
		// Declare the models that need to be activated
		inline bool is_phi_field_on = false;
		inline PhaseFieldEquation phi_equation = PhaseFieldEquation::PFE_Const;
		// - phase field data
		inline Mesh_Boundry<PhaseFieldPoint> phase_field;
		// - 
		inline int phi_property_number = 0;
		inline int phi_number = 0;
		// - 
		inline vector<int> phi_property;
		// - 
		inline REAL PHI_MAX_VARIATION = 0.0;
		// -
		namespace pair_wise_equation {
			// Statement
			enum Int_Gradient { Steinbach_1996, Steinbach_1999, Steinbach_G2009, DanielCrack_G2016 };
			enum Int_Potential { Nestler_Well, Nestler_Obstacle, Steinbach_P2009, DanielCrack_P2016 };
			enum DifferenceMethod { FIVE_POINT, NINE_POINT };
			enum Flag { pf_BULK, pf_NEAR_INTERFACE, pf_INTERFACE }; // where is the points
			const int PAIRWISE_ACC_STOP = -1;
			// pair-wise field data,
			struct PairwisePoint
			{
				std::vector<int>  active_index;    // <- phi_acc_number
				std::vector<int>  flag;            // <- phi_number
				std::vector<REAL> int_increment;   // <- phi_number
				std::vector<REAL> bulk_increment;  // <- phi_number
				void init(int phi_acc_number, int phi_number) {
					active_index.resize(phi_acc_number);
					flag.resize(phi_number);
					int_increment.resize(phi_number);
					bulk_increment.resize(phi_number);
					for (int index = 0; index < phi_acc_number; index++)
						active_index[index] = PAIRWISE_ACC_STOP;
					for (int index = 0; index < phi_number; index++) {
						flag[index] = 0;
						int_increment[index] = 0;
						bulk_increment[index] = 0;
					}
				}
				void clear() {
					active_index.clear();
					std::vector<int>().swap(active_index);
					flag.clear();
					std::vector<int>().swap(flag);
					int_increment.clear();
					std::vector<REAL>().swap(int_increment);
					bulk_increment.clear();
					std::vector<REAL>().swap(bulk_increment);
				}
			};
			inline Mesh_Boundry<PairwisePoint> phase_field_pairwise;
			inline Mesh_Boundry<REAL> phase_field_Gc;
			// - accelerated
			inline int phi_acc_number = 9;
			// - noemalization
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
			inline REAL const_xi_ab =     REAL(0.0);
			inline REAL const_xi_abc =    REAL(0.0);
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
	}

}