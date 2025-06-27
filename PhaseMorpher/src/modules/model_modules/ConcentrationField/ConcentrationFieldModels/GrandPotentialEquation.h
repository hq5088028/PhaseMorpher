#pragma once
#include "../../../data_struct/Mesh_0.h"
#include "../ConcentrationFieldPoints.h"
#include "../../PhaseField/PhaseFieldPoints.h"
namespace pf {
	using namespace std;
	// -
	namespace grand_potential_parameters {
		// - statement
		enum DifferenceMethod { FIVE_POINT, NINE_POINT };
		// - grand potnetial equation data
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
		inline vector<REAL(*)(pf::PhaseFieldPoint& point, GrandPotentialPoint& point_ex, int alpha_index, int beta_index)> source_alpha_beta;
		// - boundary condition funtions for fix value
		inline REAL(*BoundaryCondition_ConFix)(int x, int y, int z);


	}

	namespace grand_potential_functions {

		std::vector<double> pre_calculation_grand_potential_functional(double dt);

		double solve_grand_potential_functional(double dt);

	}
}