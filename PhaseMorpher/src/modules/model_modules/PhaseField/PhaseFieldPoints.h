#pragma once
#include "../../data_struct/VectorMatrix.h"
#include <vector>
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

	const int PAIRWISE_ACC_STOP = -1;

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


}