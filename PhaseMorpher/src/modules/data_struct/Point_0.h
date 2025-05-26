#pragma once
#include "../../base/sysfiles.h"
#include "VectorMatrix.h"
namespace pf {
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
	struct TemperatureFieldPoint
	{
		REAL temp;
		REAL lap_temp;
		Vector3 grad_temp;
		~TemperatureFieldPoint() {
			
		}
		TemperatureFieldPoint& operator=(const TemperatureFieldPoint& n) {
			temp = n.temp;
			lap_temp = n.lap_temp;
			grad_temp = n.grad_temp;
			return *this;
		}
	};
}