#pragma once
#include "../../base/sysfiles.h"
namespace pf {
	enum BoundaryCondition { FIXED, PERIODIC, ADIABATIC };

	template <typename T>
	class Mesh {
	public:
		// - init
		Mesh(int NX, int NY, int NZ)
			: nx(NX), ny(NY), nz(NZ), data(NX* NY* NZ) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
		}
		Mesh() = default;
		void init(int NX, int NY, int NZ) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			nx = NX;
			ny = NY;
			nz = NZ;
			data.resize(NX * NY * NZ);
		}
		~Mesh() {
			data.clear();
		}
		void clear() {
			nx = 0;
			ny = 0;
			nz = 0;
			data.clear();
			std::vector<T>().swap(data);
		}

		int Nx() const { return nx; }
		int Ny() const { return ny; }
		int Nz() const { return nz; }

		// - find element directly
		T& operator()(int x, int y, int z) {
			return data[MESH_INDEX(x, y, z, nx, ny)];
		}

	protected:
		int nx;
		int ny;
		int nz;
		std::vector<T> data;
	};

	template <typename T>
	class Mesh_Boundry : public Mesh<T> {
	public:
		// - init
		Mesh_Boundry(int NX, int NY, int NZ) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			Mesh<T>::nx = NX + 2;
			Mesh<T>::ny = NY + 2;
			Mesh<T>::nz = NZ + 2;
			Mesh<T>::data.resize((NX + 2) * (NY + 2) * (NZ + 2));
		}
		Mesh_Boundry() = default;

		void init(int NX, int NY, int NZ) {
			if (NX <= 0 || NY <= 0 || NZ <= 0) {
				std::cout << "ERROR ! size of mesh cannot be set as 0 !" << std::endl;
				SYS_PROGRAM_STOP;
			}
			Mesh<T>::nx = NX + 2;
			Mesh<T>::ny = NY + 2;
			Mesh<T>::nz = NZ + 2;
			Mesh<T>::data.resize((NX + 2) * (NY + 2) * (NZ + 2));
		}
		~Mesh_Boundry() {
			Mesh<T>::data.clear();
		}

		int X_BEGIN() { return 1; };
		int Y_BEGIN() { return 1; };
		int Z_BEGIN() { return 1; };
		int X_END() { return Mesh<T>::nx - 2; };
		int Y_END() { return Mesh<T>::ny - 2; };
		int Z_END() { return Mesh<T>::nz - 2; };

	};

}