#pragma once
#include "../postprocess_modules/WriteMeshData.h"
#include "../../tools/mathTools.h"
#include "../input_modules/ioFiles_Params.h"
#include "../input_modules/inputfiles/selectFile.h"
#include "../model_modules/Model_Params.h"
#include "../model_modules/PhaseField/PhaseField_Params.h"
#include "../Module.h"
#include "MicrostructureInit_Params.h"
#include "MicrostructureInit/GeometryStructure.h"
#include "MicrostructureInit/Bmp24Structure.h"
#include "MicrostructureInit/PorousStructure.h"
#include "MicrostructureInit/VoronoiStructure.h"
namespace pf {
	namespace microstructure_init {
		namespace functions {
			inline void init_mesh_with_datafile(pf::write_mesh_data::Data_MeshInfo& report, std::string dat_path) {
				if (write_mesh_data::read_dataFile(datafile_path, datafile_report)) {
					string _report = "> Error : datafile_path can't be opened ! \n";
					WriteDebugFile(_report);
					WriteLog(_report);
					std::exit(0);
				}
				else {
					// check datafile with this input file
					vector<string> bool_type; bool_type.push_back("MISMATCH"); bool_type.push_back("MATCH");
					stringstream _report;
					_report << "> " << std::endl;
					_report << "> | init microstructure with datafile :                    (check in this simulation)" << std::endl;
					bool all_valid = true;
					if (datafile_report.is_phi_mesh) {
						bool line_valid = true;
						if (datafile_report.PNx == simulation_mesh::phase_field.Nx()
							&& datafile_report.PNy == simulation_mesh::phase_field.Ny()
							&& datafile_report.PNz == simulation_mesh::phase_field.Nz()
							&& datafile_report.phi_number == model_parameters::phi_number) {
							line_valid = true;
						}
						else {
							line_valid = false;
							all_valid = false;
						}
						_report << "> | phase-field size               : Nx - " << datafile_report.PNx - 2 
							    << ", Ny - " << datafile_report.PNy - 2 << ", Nz - " << datafile_report.PNz - 2 
							    << ", phi number - " << datafile_report.phi_number << " (" << bool_type[line_valid] << ")" << std::endl;
					}
					_report << ">" << std::endl;
					WriteLog(_report.str());
					if (all_valid == false) {
						string error_report = "> Wainning : mesh data structure from datafile and inputfile mismatch ! \n";
						WriteDebugFile(error_report);
						WriteLog(error_report);
					}
				}
			}
		}
		// -
		inline void init_microstructure_pre_i() {
			using namespace simulation_mesh;
			if (is_datafile_init) {
				WriteLog("> Open datafile : " + datafile_path + "\n");
				functions::init_mesh_with_datafile(datafile_report, datafile_path);
			}
			else {
				if (model_parameters::is_phi_field_on) {
#pragma omp parallel for
					for (int x = 0; x < phase_field.Nx(); x++)
						for (int y = 0; y < phase_field.Ny(); y++)
							for (int z = 0; z < phase_field.Nz(); z++)
								phase_field(x, y, z).phi[matrix_phi_index] = matrix_phi_value;
				}
				if (model_parameters::is_con_field_on) {
#pragma omp parallel for
					for (int x = 0; x < concentration_field.Nx(); x++)
						for (int y = 0; y < concentration_field.Ny(); y++)
							for (int z = 0; z < concentration_field.Nz(); z++)
								concentration_field(x, y, z).con = matrix_con;
				}
				if (model_parameters::is_temp_field_on) {
#pragma omp parallel for
					for (int x = 0; x < temperature_field.Nx(); x++)
						for (int y = 0; y < temperature_field.Ny(); y++)
							for (int z = 0; z < temperature_field.Nz(); z++)
								temperature_field(x, y, z).temp = matrix_temperature;
				}
				// - normal init structure
				geometry_structure::definiteNucleation();
				if (is_read_bmp24file) {
					bmp24_structure::generate_structure_from_BMP_pic();
					geometry_structure::definiteNucleation();
				}
				// - porous
				if (is_porous) {
					if (porous_phis_indexs.size() != 0) {
						porous_structure::quartet_structure_generation_in_phis();
					}
					else {
						porous_structure::quartet_structure_generation();
					}
					geometry_structure::definiteNucleation();
				}
				// - voronoi
				if (is_voronoi) {
					if (voronoi_phis_indexs.size() != 0) {
						voronoi_structure::generate_voronoi_structure_in_phis();
					}
					else {
						voronoi_structure::generate_voronoi_structure();
						// voronoi_structure::generate_voronoi_structure_new_method();
					}
					geometry_structure::definiteNucleation();
				}
				// - others
			}
		}
		inline void dinit() {
			nucleation_box.geometry_box.clear();
			nucleation_box.point_set_box.clear();
		}
		// -
		inline void init() {
			InputFileReader::get_instance()->read_bool_value("Preprocess.Microstructure.is_datafile_init", is_datafile_init, true);
			if (is_datafile_init) {
				WriteDebugFile("# Preprocess.Microstructure.datafile_path : relative path from infile folder.\n");
				if (InputFileReader::get_instance()->read_string_value("Preprocess.Microstructure.datafile_path", datafile_path, true))
					is_read_datafile_by_path = true;
				datafile_path = input_output_files_parameters::InFile_Path + dirSeparator + datafile_path;
#ifdef _WIN32
				if (!is_read_datafile_by_path) {
					WriteLog("> Please select a datafile (in .dat format) to initialize the simulation mesh ...");
					selectDataFile(datafile_path);
				}
#else
				if (!is_read_datafile_by_path) {
					WriteDebugFile("> Error : Preprocess.Microstructure.datafile_path should be defined ! \n");
					std::exit(0);
				}
#endif
			}
			else {
				// - init matrix 
				if (model_parameters::is_phi_field_on) {
					string matrix_key = "Preprocess.Microstructure.matrix_phi", matrix_string = "()";
					WriteDebugFile("# .matrix_phi = ( phi_index, phi_value ) \n");
					if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_string, true)) {
						vector<InputValueType> matrix_structure;
						matrix_structure.push_back(InputValueType::IVType_INT);
						matrix_structure.push_back(InputValueType::IVType_REAL);
						vector<input_value> matrix_value = 
							InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(matrix_structure, matrix_key, matrix_string, true);
						matrix_phi_index = matrix_value[0].int_value;
						matrix_phi_value = matrix_value[1].REAL_value;
						check_phi_index(matrix_phi_index);
					}
				}
				if (model_parameters::is_con_field_on) {
					string matrix_key = "Preprocess.Microstructure.matrix_con", matrix_string = "()";
					WriteDebugFile("# .matrix_con = ( comp_0_value, comp_1_value, ... ) \n");
					matrix_con.resize(model_parameters::con_number);
					if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_string, true)) {
						vector<input_value> matrix_value = 
							InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, matrix_key, matrix_string, true);
						for (int index = 0; index < matrix_value.size(); index++)
							matrix_con[index] = matrix_value[index].REAL_value;
						check_con_size(int(matrix_con.size()));
					}
				}
				if (model_parameters::is_temp_field_on) {
					string matrix_key = "Preprocess.Microstructure.matrix_temperature";
					InputFileReader::get_instance()->read_REAL_value(matrix_key, matrix_temperature, true);
				}
			}
			// - init geometry structure
			geometry_structure::init();
			// - init bmp24 structure
			bmp24_structure::init();
			// - init porous structure
			porous_structure::init();
			// - init voronoi structure
			voronoi_structure::init();
			// - others

			load_a_new_module(init_microstructure_pre_i, default_module_function, default_module_function,
				default_module_function, default_module_function, default_module_function,
				default_module_function, default_module_function, default_module_function, default_module_function);
		}
	}
}