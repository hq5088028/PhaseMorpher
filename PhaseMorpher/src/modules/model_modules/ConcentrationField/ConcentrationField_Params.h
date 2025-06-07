#pragma once
#include "../Model_Params.h"
#include "../../input_modules/ioFiles_Params.h"
#include "../../input_modules/inputfiles/InputFileReader.h"
namespace pf {
	enum ConcentrationFieldEquation {
		CFE_Const = 0, CFE_TotalConcentration, CFE_GrandPotential
	};
	struct ConcentrationFieldPoint
	{
		std::vector<REAL> con;
		std::vector<REAL> lap_con;
		std::vector<Vector3> grad_con;
		void set_con_number(int Number) {
			con.resize(Number);
			lap_con.resize(Number);
			grad_con.resize(Number);
			for (int index = 0; index < Number; index++) {
				con[index] = 0.0;
				lap_con[index] = 0.0;
				grad_con[index].set_to_zero();
			}
		}
		void normalize_con() {
			for (int index = 0; index < con.size(); index++) {
				REAL& CON = con[index];
				if (CON >= 1.0)
					CON = SYS_EPSILON_R;
				else if (CON <= 0.0)
					CON = SYS_EPSILON;
			}
		}
		void clear() {
			con.clear();
			std::vector<REAL>().swap(con);
			lap_con.clear();
			std::vector<REAL>().swap(lap_con);
			grad_con.clear();
			std::vector<Vector3>().swap(grad_con);
		}
		~ConcentrationFieldPoint() {
			clear();
		}
		ConcentrationFieldPoint& operator=(const ConcentrationFieldPoint& n) {
			con = n.con;
			lap_con = n.lap_con;
			grad_con = n.grad_con;
			return *this;
		}
	};
	namespace con_parameters {
		// Declare the models that need to be activated
		inline bool is_con_field_on = false;
		inline ConcentrationFieldEquation con_equation = ConcentrationFieldEquation::CFE_Const;
		// - concentration field data
		inline Mesh_Boundry<ConcentrationFieldPoint> concentration_field;
		// - 
		inline int con_number = 0;
		// - 
		inline REAL CON_MAX_VARIATION = 0.0;
		// -
		namespace grand_potential_equation {
			// - statement
			enum DifferenceMethod { FIVE_POINT, NINE_POINT };
			// - grand potnetial equation data
			struct GrandPotentialPoint {
				// - effect region
				std::vector<REAL> effect_phi;  // <- to describe the effect regions
				std::vector<REAL> effect_phi_old;  // <- to describe the effect region change
				// - grand potential extra
				std::vector<REAL> grand_potential_ex;  // <- grand potential contributed by other energy
				// - concentration

			};
			inline Mesh_Boundry<GrandPotentialPoint> grand_potential_field;
			// - parameters
			inline DifferenceMethod diff_method = FIVE_POINT;
			inline int effect_region_size = 0;
			inline vector<int> effect_solvent;
			inline vector<vector<int>> effect_phi_indexes;
			inline vector<REAL> effect_threshold;
			inline vector<vector<REAL>> effect_grand_potential_range;
			// - diffusion potential functions
			inline vector<REAL(*)(pf::ConcentrationFieldPoint& point, GrandPotentialPoint& point_ex, int phi_index, int con_index)> delt_Fchem_delt_con;
			inline vector<REAL(*)(pf::ConcentrationFieldPoint& point, GrandPotentialPoint& point_ex, int phi_index, int con_index)> delt_Fother_delt_con;
			// - source functions
			inline vector<REAL(*)(pf::PhaseFieldPoint& point, PairwisePoint& point_ex, int alpha_index, int beta_index)> source_alpha_beta;
			// - boundary condition funtions for fix value
			inline REAL(*BoundaryCondition_ConFix)(int x, int y, int z);


		}
	}
}