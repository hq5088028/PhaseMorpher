#pragma once
#include "../../base/sysfiles.h"
#include "../model_modules/Model_Params.h"
#include "../input_modules/ioFiles_Params.h"
#include "../model_modules/PhaseField/PhaseField_Params.h"
#include "../MainIterator_Params.h"
#include "../Module.h"
namespace pf {

	namespace write_vts {
		inline int output_frequence = -1;
		inline bool is_show_with_boundary = false;
		inline int x_begin = 0;
		inline int y_begin = 0;
		inline int z_begin = 0;
		inline int x_end = 0;
		inline int y_end = 0;
		inline int z_end = 0;
		inline std::vector<void(*)(std::ofstream& fout)> write_vts_scalar_list;
		inline std::vector<void(*)(std::ofstream& fout)> write_vts_vector_list;
		namespace default_functions {
			inline void open_vts_scalar_file(std::ofstream& fout, std::string tail) {
				std::string fname;
				fname = input_output_files_parameters::WorkingFolder_Path + dirSeparator + "scalar_variables_" + tail + ".vts";
				fout.open(fname);
				if (!fout) {
					std::cout << "Failed to write the vtk file..." << std::endl;
					fout.close();
					return;
				}
				fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << std::endl;
				fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
				fout << "<StructuredGrid WholeExtent=\""
					<< x_begin << " " << x_end << " "
					<< y_begin << " " << y_end << " "
					<< z_begin << " " << z_end << "\"> " << std::endl;
				fout << "<PointData Scalars= \"ScalarData\">" << std::endl;
			}
			inline void open_vts_vec3_file(std::ofstream& fout, std::string tail) {
				std::string fname;
				fname = input_output_files_parameters::WorkingFolder_Path + dirSeparator + "vec3_variables_" + tail + ".vts";
				fout.open(fname);
				if (!fout) {
					std::cout << "Failed to write the vtk file..." << std::endl;
					fout.close();
					return;
				}
				fout << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << std::endl;
				fout << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
				fout << "<StructuredGrid WholeExtent=\""
					<< 0 << " " << x_end - x_begin + 1 << " "
					<< 0 << " " << y_end - y_begin + 1 << " "
					<< 0 << " " << z_end - z_begin + 1 << "\"> " << std::endl;
				fout << "<PointData  Vectors= \"VectorData\">" << std::endl;
			}
			inline void close_vts_file(std::ofstream& fout) {
				fout << "</PointData>" << std::endl;
				fout << "<Points>" << std::endl;
				fout << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
				for (int k = z_begin; k <= z_end; ++k)
					for (int j = y_begin; j <= y_end; ++j)
						for (int i = x_begin; i <= x_end; ++i) {
							fout << i * mesh_parameters::delt_r << " " << j * mesh_parameters::delt_r << " " << k * mesh_parameters::delt_r << "\n";
						}
				fout << "</DataArray>" << std::endl;
				fout << "</Points>" << std::endl;
				fout << "</StructuredGrid>" << std::endl;
				fout << "</VTKFile>" << std::endl;
				fout.close();
			}
		}
		inline void write_vts_pre_iii() {
			ofstream fout;
			// - 
			default_functions::open_vts_scalar_file(fout, "step0");
			for (auto writer = write_vts_scalar_list.begin(); writer < write_vts_scalar_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
			// - 
			default_functions::open_vts_vec3_file(fout, "step0");
			for (auto writer = write_vts_vector_list.begin(); writer < write_vts_vector_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
		}

		inline void write_vts_pos_iii() {
			if (main_iterator::Current_ITE_step % output_frequence != 0)
				return;
			ofstream fout;
			// - 
			default_functions::open_vts_scalar_file(fout, "step" + to_string(main_iterator::Current_ITE_step));
			for (auto writer = write_vts_scalar_list.begin(); writer < write_vts_scalar_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
			// - 
			default_functions::open_vts_vec3_file(fout, "step" + to_string(main_iterator::Current_ITE_step));
			for (auto writer = write_vts_vector_list.begin(); writer < write_vts_vector_list.end(); writer++)
				(*writer)(fout);
			default_functions::close_vts_file(fout);
		}

		inline void init() {
			InputFileReader::get_instance()->read_int_value("Solver.Output.VTS.frequence", output_frequence, true);
			if (output_frequence == 0) {
				load_a_new_module(default_module_function, default_module_function, write_vts_pre_iii,
					default_module_function, default_module_function, default_module_function,
					default_module_function, default_module_function, default_module_function, default_module_function);
			}
			else if (output_frequence > 0) {
				load_a_new_module(default_module_function, default_module_function, write_vts_pre_iii,
					default_module_function, default_module_function, default_module_function,
					default_module_function, default_module_function, write_vts_pos_iii, default_module_function);
			}
			InputFileReader::get_instance()->read_bool_value("Solver.Output.VTS.show_with_boundary", is_show_with_boundary, true);
			if (is_show_with_boundary) {
				x_begin = 0;
				y_begin = 0;
				z_begin = 0;
				x_end = phi_parameters::phase_field.Nx() - 1;
				y_end = phi_parameters::phase_field.Ny() - 1;
				z_end = phi_parameters::phase_field.Nz() - 1;
			}
			else {
				x_begin = phi_parameters::phase_field.X_BEGIN();
				y_begin = phi_parameters::phase_field.Y_BEGIN();
				z_begin = phi_parameters::phase_field.Z_BEGIN();
				x_end = phi_parameters::phase_field.X_END();
				y_end = phi_parameters::phase_field.Y_END();
				z_end = phi_parameters::phase_field.Z_END();
			}
		}
	}

	inline void load_vts_scalar_func(void(*buff)(std::ofstream& fout)) {
		write_vts::write_vts_scalar_list.push_back(buff);
	}
	inline void load_vts_vector_func(void(*buff)(std::ofstream& fout)) {
		write_vts::write_vts_vector_list.push_back(buff);
	}
}