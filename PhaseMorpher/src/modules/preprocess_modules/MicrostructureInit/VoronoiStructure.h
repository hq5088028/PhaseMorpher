#pragma once
#include "../MicrostructureInit_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
namespace pf {
	namespace voronoi_structure {
		// -
		void init();
		// -
		void generate_voronoi_structure();
		void generate_voronoi_structure_new_method();
		// -
		void generate_voronoi_structure_in_phis();
		// - 
		inline void write_scalar_voronoi_points(ofstream& fout) {
			fout << "<DataArray type = \"Float64\" Name = \"" << "voronoi_points" <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (int j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (int i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						bool is_point_here = false;
						for (int index = 0; index < microstructure_init::voronoi_points.size(); index++) {
							int pos_x = REAL_to_int(microstructure_init::voronoi_points[index].x),
								pos_y = REAL_to_int(microstructure_init::voronoi_points[index].y),
								pos_z = REAL_to_int(microstructure_init::voronoi_points[index].z);
							if (pos_x == i && pos_y == j && pos_z == k)
								is_point_here = true;
						}
						if (is_point_here)
							fout << 1 << endl;
						else
							fout << 0 << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
}