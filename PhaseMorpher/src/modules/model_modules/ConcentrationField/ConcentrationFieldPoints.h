#pragma once
#include <vector>
#include "../../data_struct/VectorMatrix.h"
#include "../../../base/MACRO_DEF.h"
namespace pf {

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

	struct GrandPotentialPoint {
		// - effect region
		std::vector<REAL> effect_phi;  // <- to describe the effect regions
		std::vector<REAL> effect_phi_old;  // <- to describe the effect region change
		// - grand potential extra
		std::vector<REAL> grand_potential_ex;  // <- grand potential contributed by other energy
		// - concentration

	};

}