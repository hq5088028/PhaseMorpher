#include "Bmp24Structure.h"
namespace pf {
	namespace bmp24_structure {
		void generate_structure_from_BMP_pic() {
#ifdef _WIN32
			using namespace simulation_mesh;
			if (mesh_parameters::MESH_NZ > 1) {
				string error_report = "> Warnning : generate structure error : Nz > 1, cant init structure by a BMP picture !\n";
				WriteDebugFile(error_report);
				return;
			}
			for (int layer_index = 0; layer_index < microstructure_init::bmp24_layer; layer_index++) {
				BMP24reader bmpReader;
				bmpReader.safe(microstructure_init::bmp24file_path);
				bmpReader.read(microstructure_init::bmp24file_path);
				PointSet set;
				for (int y = 0; y < mesh_parameters::MESH_NY; y++)
					for (int x = 0; x < mesh_parameters::MESH_NX; x++) {
						int xx = x * bmpReader.bmp_width / mesh_parameters::MESH_NX,
							yy = y * bmpReader.bmp_height / mesh_parameters::MESH_NY;
						REAL graypercent = bmpReader.getGrayPercentage(xx, yy);
						if (graypercent > microstructure_init::bmp24_threshold[layer_index][0] && graypercent < microstructure_init::bmp24_threshold[layer_index][1])
							set.add_point(x + 1, y + 1, 1, microstructure_init::bmp24_phi_value[layer_index]);
					}
				set.generate_step = 0;
				if (model_parameters::is_phi_field_on) {
					set.phaseIndex = microstructure_init::bmp24_phi_index[layer_index];
					set.is_normalized = microstructure_init::bmp24_phi_normalized[layer_index];
				}
				if (model_parameters::is_con_field_on)
					set.con = microstructure_init::bmp24_con[layer_index];
				if (model_parameters::is_temp_field_on)
					set.temperature = microstructure_init::bmp24_temperature[layer_index];
				microstructure_init::nucleation_box.point_set_box.push_back(set);
			}
#endif
		}

		void init() {
#ifdef _WIN32
			if (InputFileReader::get_instance()->read_string_value("Preprocess.Microstructure.bmp24.file_path", microstructure_init::bmp24file_path, true)) {
				microstructure_init::is_read_bmp24file = true;
				if (microstructure_init::bmp24file_path == "") {
					selectAllFile(microstructure_init::bmp24file_path);
				}
				else {
					microstructure_init::bmp24file_path = input_output_files_parameters::InFileFolder_Path + dirSeparator + microstructure_init::bmp24file_path;
				}
				InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.bmp24.layer_number", microstructure_init::bmp24_layer, true);
				for (int bmp_layer = 0; bmp_layer < microstructure_init::bmp24_layer; bmp_layer++) {
					// - gray_threshold
					WriteDebugFile("# .gray_threshold = (range_left, range_right) \n");
					string bmp24_threshold_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".gray_threshold", bmp24_threshold_input = "(0,0)";
					InputFileReader::get_instance()->read_string_value(bmp24_threshold_key, bmp24_threshold_input, true);
					vector<input_value> bmp24_threshold_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, bmp24_threshold_key, bmp24_threshold_input, true);
					microstructure_init::bmp24_threshold.push_back({ bmp24_threshold_value[0].REAL_value, bmp24_threshold_value[1].REAL_value });
					// - phi
					if (model_parameters::is_phi_field_on) {
						WriteDebugFile("# .phi = ( phi_index, phi_value, is_normalized ) \n");
						string bmp24_phi_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".phi", bmp24_phi_input = "(0,0,false)";
						InputFileReader::get_instance()->read_string_value(bmp24_phi_key, bmp24_phi_input, true);
						vector<InputValueType> bmp24_phi_structure; bmp24_phi_structure.push_back(InputValueType::IVType_INT);
						bmp24_phi_structure.push_back(InputValueType::IVType_REAL); bmp24_phi_structure.push_back(InputValueType::IVType_BOOL);
						vector<input_value> bmp24_phi_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(bmp24_phi_structure, bmp24_phi_key, bmp24_phi_input, true);
						microstructure_init::bmp24_phi_index.push_back(bmp24_phi_value[0].int_value);
						check_phi_index(bmp24_phi_value[0].int_value);
						microstructure_init::bmp24_phi_value.push_back(bmp24_phi_value[1].REAL_value);
						microstructure_init::bmp24_phi_normalized.push_back(bmp24_phi_value[2].bool_value);
					}
					// - con
					if (model_parameters::is_con_field_on) {
						vector<REAL> con;
						WriteDebugFile("# .con = ( comp_0_value, comp_1_value, ... ) \n");
						string bmp24_x_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".con", bmp24_x_input = "()";
						if (InputFileReader::get_instance()->read_string_value(bmp24_x_key, bmp24_x_input, true)) {
							vector<input_value> bmp24_x_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, bmp24_x_key, bmp24_x_input, true);
							for (int x_index = 0; x_index < bmp24_x_value.size(); x_index++)
								con.push_back(bmp24_x_value[x_index].REAL_value);
							check_con_size(int(con.size()));
							microstructure_init::bmp24_con.push_back(con);
						}
					}
					// - temperature
					if (model_parameters::is_temp_field_on) {
						string bmp24_temp_key = "Preprocess.Microstructure.bmp24_layer_" + to_string(bmp_layer) + ".temperature"; REAL bmp_temp = 0.0;
						InputFileReader::get_instance()->read_REAL_value(bmp24_temp_key, bmp_temp, true);
						microstructure_init::bmp24_temperature.push_back(bmp_temp);
					}
				}
			}
#endif
		}
	}
}