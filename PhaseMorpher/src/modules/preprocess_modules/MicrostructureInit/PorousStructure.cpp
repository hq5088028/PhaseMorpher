#include "PorousStructure.h"
namespace pf {
	namespace porous_structure {
		void init() {
			using namespace microstructure_init;
			WriteDebugFile("# Preprocess.Microstructure.porous = (first_phi_index, second_phi_index, porosity, noise_level) \n");
			string porous_key = "Preprocess.Microstructure.porous", porous_input = "()";
			if (InputFileReader::get_instance()->read_string_value(porous_key, porous_input, true)) {
				is_porous = true;
				vector<InputValueType> porous_structure = { InputValueType::IVType_INT ,
					InputValueType::IVType_INT ,InputValueType::IVType_REAL, InputValueType::IVType_REAL };
				vector<input_value> porous_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(porous_structure, porous_key, porous_input, true);
				porous_first_phi_index = porous_value[0].int_value;
				porous_second_phi_index = porous_value[1].int_value;
				check_phi_index(porous_first_phi_index);
				check_phi_index(porous_second_phi_index);
				porosity = porous_value[2].REAL_value;
				porous_init_noise = porous_value[3].REAL_value;
				if (model_parameters::is_con_field_on) {
					WriteDebugFile("# .con = [(first_comp_0_value, first_comp_0_value, ... ), (second_comp_0_value, second_comp_0_value, ... )] \n");
					string porous_x_key = "Preprocess.Microstructure.Porous.con", porous_x_input = "[()]";
					if (InputFileReader::get_instance()->read_string_value(porous_x_key, porous_x_input, true)) {
						vector<vector<input_value>> porous_x_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, porous_x_key, porous_x_input, true);
						for (int x_index = 0; x_index < porous_x_value[0].size(); x_index++)
							porous_first_con.push_back(porous_x_value[0][x_index].REAL_value);
						for (int x_index = 0; x_index < porous_x_value[1].size(); x_index++)
							porous_second_con.push_back(porous_x_value[1][x_index].REAL_value);
						check_con_size(int(porous_first_con.size()));
						check_con_size(int(porous_second_con.size()));
					}
				}
				if (model_parameters::is_temp_field_on) {
					WriteDebugFile("# .temperature = (first_temperature, second_temperature) \n");
					string porous_temp_key = "Preprocess.Microstructure.Porous.temperature", porous_temp_input = "(0,0)";
					InputFileReader::get_instance()->read_string_value(porous_temp_key, porous_temp_input, true);
					vector<input_value> porous_temp_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, porous_temp_key, porous_temp_input, true);
					porous_first_temperature = porous_temp_value[0].REAL_value;
					porous_second_temperature = porous_temp_value[1].REAL_value;
				}

				string porous_norm_key = "Preprocess.Microstructure.Porous.is_normalized";
				InputFileReader::get_instance()->read_bool_value(porous_norm_key, is_porous_normalized, true);

				if (InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.Porous.rand_seed", porous_rand_seed, true))
					is_porous_rand = false;

				WriteDebugFile("# .in_phi_indexs = ( phi_index_1, phi_index_2, ... ) \n");
				string porous_in_phis_key = "Preprocess.Microstructure.Porous.in_phi_indexs", porous_in_phis_input = "()";
				InputFileReader::get_instance()->read_string_value(porous_in_phis_key, porous_in_phis_input, true);
				vector<input_value> porous_in_phis_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, porous_in_phis_key, porous_in_phis_input, true);
				for (int in_phi_index = 0; in_phi_index < porous_in_phis_value.size(); in_phi_index++) {
					porous_phis_indexs.push_back(porous_in_phis_value[in_phi_index].int_value);
					check_phi_index(porous_in_phis_value[in_phi_index].int_value);
				}
			}
		}

		void quartet_structure_grow_2D(vector<vector<int>>& Solid, std::mt19937& gen) {
			using namespace microstructure_init;
			vector<vector<int>> buff_solid;
			const int Nx = int(Solid.size()), Ny = int(Solid[0].size());
			buff_solid.resize(Nx);
#pragma omp parallel for
			for (int i = 0; i < Nx; i++) {
				buff_solid[i].resize(Ny);
				for (int j = 0; j < Ny; j++) {
					buff_solid[i][j] = Solid[i][j];
				}
			}
			/// \brief Grow in eight directions.
			/// Grow directions
			///*****    6    2    5   *****
			///*****                  *****
			///*****    3    C    1   *****
			///*****                  *****
			///*****    7    4    8   *****
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
#pragma omp parallel for
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					if (i == 0) {
						if (j == 0) {
							if (Solid[i][j] == 1)
							{
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
							}
						}
						else if (j == Ny - 1) {
							if (Solid[i][j] == 1)
							{
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
						else {
							if (Solid[i][j] == 1)
							{
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
					}
					else if (i == Nx - 1) {
						if (j == 0) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
							}
						}
						else if (j == Ny - 1) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
							}
						}
						else {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
							}
						}
					}
					else {
						if (j == 0) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
							}
						}
						else if (j == Ny - 1) {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
						else {
							if (Solid[i][j] == 1) {
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i + 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i - 1][j] = 1;
								if (real_dist(gen) < porous_TwoD_d1) buff_solid[i][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j + 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i - 1][j - 1] = 1;
								if (real_dist(gen) < porous_TwoD_d5) buff_solid[i + 1][j - 1] = 1;
							}
						}
					}
				}
			}
#pragma omp parallel for
			for (int i = 0; i < Nx; i++) {
				for (int j = 0; j < Ny; j++) {
					Solid[i][j] = buff_solid[i][j];
				}
			}
		}

		void quartet_structure_grow_3D(vector<vector<vector<int>>>& arrgrid, vector<vector<int>>& soild, int& Tnumsoild, std::mt19937& gen) {
			using namespace microstructure_init;
			const int NX = int(arrgrid.size()), NY = int(arrgrid[0].size()), NZ = int(arrgrid[0][0].size());
			int numsoild = Tnumsoild;
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			for (int index_soild = 0; index_soild < Tnumsoild; index_soild++) {
				int index_i = soild[index_soild][0];
				int index_j = soild[index_soild][1];
				int index_k = soild[index_soild][2];
				//1向右方向生长
				if (index_j < NY - 1) {
					int i = index_i;
					int j = index_j + 1;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//2向后方向生长
				if (index_i > 0) {
					int i = index_i - 1;
					int j = index_j;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//3向左方向生长
				if (index_j > 0) {
					int i = index_i;
					int j = index_j - 1;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//4向前方向生长
				if (index_i < NX - 1) {
					int i = index_i + 1;
					int j = index_j;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//5向上方向生长
				if (index_k < NZ - 1) {
					int i = index_i;
					int j = index_j;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//6向下方向生长		
				if (index_k > 0) {
					int i = index_i;
					int j = index_j;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d1) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//7向水平方向右前生长
				if (index_i < NX - 1 && index_j < NY - 1) {
					int i = index_i + 1;
					int j = index_j + 1;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//8向水平方向左前生长
				if (index_i < NX - 1 && index_j > 0) {
					int i = index_i + 1;
					int j = index_j - 1;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//9向水平方向右后生长
				if (index_i > 0 && index_j < NY - 1) {
					int i = index_i - 1;
					int j = index_j + 1;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//10向水平方向左后生长
				if (index_i > 0 && index_j > 0) {
					int i = index_i - 1;
					int j = index_j - 1;
					int k = index_k;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//11向右上方向生长
				if (index_j < NY - 1 && index_k < NZ - 1) {
					int i = index_i;
					int j = index_j + 1;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//12向右下方向生长
				if (index_j < NY - 1 && index_k >0) {
					int i = index_i;
					int j = index_j + 1;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//13向左上方向生长
				if (index_j > 0 && index_k < NZ - 1) {
					int i = index_i;
					int j = index_j - 1;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//14向左下方向生长
				if (index_j > 0 && index_k > 0) {
					int i = index_i;
					int j = index_j - 1;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//15向前上方向生长
				if (index_i < NX - 1 && index_k < NZ - 1) {
					int i = index_i + 1;
					int j = index_j;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//16向前下方向生长
				if (index_i < NX - 1 && index_k >0) {
					int i = index_i + 1;
					int j = index_j;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//17向后上方向生长
				if (index_i > 0 && index_k < NZ - 1) {
					int i = index_i - 1;
					int j = index_j;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//18向后下方向生长
				if (index_i > 0 && index_k > 0) {
					int i = index_i - 1;
					int j = index_j;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d7) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//19向右前上对角线方向生长
				if (index_i < NX - 1 && index_j < NY - 1 && index_k < NZ - 1) {
					int i = index_i + 1;
					int j = index_j + 1;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//20向右后上对角线方向生长
				if (index_i > 0 && index_j < NY - 1 && index_k < NZ - 1) {
					int i = index_i - 1;
					int j = index_j + 1;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//21向左后上对角线方向生长
				if (index_i > 0 && index_j > 0 && index_k < NZ - 1) {
					int i = index_i - 1;
					int j = index_j - 1;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//22向左前上对角线方向生长
				if (index_i < NX - 1 && index_j>0 && index_k < NZ - 1) {
					int i = index_i + 1;
					int j = index_j - 1;
					int k = index_k + 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//23向右前下对角线方向生长
				if (index_i < NX - 1 && index_j < NY - 1 && index_k>0) {
					int i = index_i + 1;
					int j = index_j + 1;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//24向右后下对角线方向生长
				if (index_i > 0 && index_j < NY - 1 && index_k>0) {
					int i = index_i - 1;
					int j = index_j + 1;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//25向左后下对角线方向生长
				if (index_i > 0 && index_j > 0 && index_k > 0) {
					int i = index_i - 1;
					int j = index_j - 1;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
				//26向左前下对角线方向生长
				if (index_i < NX - 1 && index_j>0 && index_k > 0) {
					int i = index_i + 1;
					int j = index_j - 1;
					int k = index_k - 1;
					if (arrgrid[i][j][k] == 0 && real_dist(gen) < porous_ThreeD_d19) {
						numsoild = numsoild + 1;
						arrgrid[i][j][k] = 1;
						vector<int> new_solid = { i, j, k };
						soild.push_back(new_solid);
					}
				}
			}
			Tnumsoild = numsoild;
		}

		void quartet_structure_generation() {
			using namespace microstructure_init;
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!is_porous_rand) {
				gen.seed(porous_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			WriteLog("> Quartet structure generation:");
			if (mesh_parameters::dimention == Dimension::Two_Dimension) {
				vector<vector<int>> Solid;
				Solid.resize(mesh_parameters::MESH_NX);
				/// No solid in the beginning
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					Solid[i].resize(mesh_parameters::MESH_NY);
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						Solid[i][j] = 0;
					}
				}
				/// Produce growth core
				double core_p = 0.0;
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						if (real_dist(gen) < porous_init_noise) {
							Solid[i][j] = 1;
							core_p += 1.0;
						}
					}
				}
				core_p /= mesh_parameters::MESH_NX * mesh_parameters::MESH_NY;
				string out_ptr = "  percentage of core in 2D = " + to_string(core_p) + "\n";
				WriteLog(out_ptr);
				///  Produce process
				int step = 0;
				do {
					step++;
					quartet_structure_grow_2D(Solid, gen);
					/// Calculate porosity
					core_p = 0.0;
					for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
						for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
							core_p += Solid[i][j];
						}
					}
					core_p = core_p / (mesh_parameters::MESH_NX * mesh_parameters::MESH_NY);
					out_ptr = "  grow step: " + to_string(step) + ", percentage of particle/pore = " + to_string(core_p) + "\n";
					WriteLog(out_ptr);
				} while (core_p < porosity);
				out_ptr = "  End of growth ! \n";
				WriteLog(out_ptr);
				PointSet set_phi_1, set_phi_2;
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						if (Solid[i][j]) {
							set_phi_1.add_point(REAL(i + 1), REAL(j + 1), REAL(1), REAL(1.0));
						}
						else {
							set_phi_2.add_point(REAL(i + 1), REAL(j + 1), REAL(1), REAL(1.0));
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = microstructure_init::porous_first_phi_index;
				set_phi_1.temperature = microstructure_init::porous_first_temperature;
				set_phi_1.con = microstructure_init::porous_first_con;
				set_phi_1.is_normalized = microstructure_init::is_porous_normalized;
				
				nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = microstructure_init::porous_second_phi_index;
				set_phi_2.temperature = microstructure_init::porous_second_temperature;
				set_phi_2.con = microstructure_init::porous_second_con;
				set_phi_2.is_normalized = microstructure_init::is_porous_normalized;
				nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else if (mesh_parameters::dimention == Dimension::Three_Dimension) {
				vector<vector<vector<int>>> arrgrid; arrgrid.resize(mesh_parameters::MESH_NX);
				vector<vector<int>> soild;
				int numsoild = 0, Nxyz = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ, numtotal_need = int(porosity * Nxyz);
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					arrgrid[i].resize(mesh_parameters::MESH_NY);
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						arrgrid[i][j].resize(mesh_parameters::MESH_NZ);
						for (int k = 0; k < mesh_parameters::MESH_NZ; k++) {
							arrgrid[i][j][k] = 0;
							if (real_dist(gen) < porous_init_noise) {
								arrgrid[i][j][k] = 1;
								vector<int> buff = { i, j, k };
								soild.push_back(buff);
								numsoild = numsoild + 1;
							}
						}
					}
				}
				double core_p = double(numsoild) / double(Nxyz);
				string out_ptr = "  percentage of core in 3D = " + to_string(core_p) + "\n";
				WriteLog(out_ptr);
				int Tnumsoild = numsoild, step = 0;
				while (Tnumsoild < numtotal_need) {
					step++;
					quartet_structure_grow_3D(arrgrid, soild, Tnumsoild, gen);
					core_p = double(Tnumsoild) / double(Nxyz);
					out_ptr = "  grow step: " + to_string(step) + ", percentage of particle/pore = " + to_string(core_p) + "\n";
					WriteLog(out_ptr);
				}
				out_ptr = "  End of growth ! \n";
				WriteLog(out_ptr);
				PointSet set_phi_1, set_phi_2;
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						for (int k = 0; k < mesh_parameters::MESH_NZ; k++) {
							if (arrgrid[i][j][k]) {
								set_phi_1.add_point(REAL(i + 1), REAL(j + 1), REAL(k + 1), REAL(1.0));
							}
							else {
								set_phi_2.add_point(REAL(i + 1), REAL(j + 1), REAL(k + 1), REAL(1.0));
							}
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = microstructure_init::porous_first_phi_index;
				set_phi_1.temperature = microstructure_init::porous_first_temperature;
				set_phi_1.con = microstructure_init::porous_first_con;
				set_phi_1.is_normalized = microstructure_init::is_porous_normalized;
				nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = microstructure_init::porous_second_phi_index;
				set_phi_2.temperature = microstructure_init::porous_second_temperature;
				set_phi_2.con = microstructure_init::porous_second_con;
				set_phi_2.is_normalized = microstructure_init::is_porous_normalized;
				nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else {
				WriteLog("> Quartet structure generation Failed, can't generate at one dimension ! \n");
			}
		}

		void quartet_structure_generation_in_phis() {
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!microstructure_init::is_porous_rand) {
				gen.seed(microstructure_init::porous_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			WriteLog("> Quartet structure generation:");
			if (mesh_parameters::dimention == Dimension::Two_Dimension) {
				vector<vector<int>> Solid;
				Solid.resize(mesh_parameters::MESH_NX);
				/// No solid in the beginning
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					Solid[i].resize(mesh_parameters::MESH_NY);
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						Solid[i][j] = 0;
					}
				}
				/// Produce growth core
				REAL core_p = 0.0;
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						if (real_dist(gen) < microstructure_init::porous_init_noise) {
							Solid[i][j] = 1;
							core_p += 1.0;
						}
					}
				}
				core_p /= mesh_parameters::MESH_NX * mesh_parameters::MESH_NY;
				string out_ptr = "  percentage of core in 2D = " + to_string(core_p) + "\n";
				WriteLog(out_ptr);
				///  Produce process
				int step = 0;
				do {
					step++;
					quartet_structure_grow_2D(Solid, gen);
					/// Calculate porosity
					core_p = 0.0;
					for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
						for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
							core_p += Solid[i][j];
						}
					}
					core_p = core_p / (mesh_parameters::MESH_NX * mesh_parameters::MESH_NY);
					out_ptr = "  grow step: " + to_string(step) + ", percentage of particle/pore = " + to_string(core_p) + "\n";
					WriteLog(out_ptr);
				} while (core_p < microstructure_init::porosity);
				out_ptr = "  End of growth ! \n";
				WriteLog(out_ptr);
				PointSet set_phi_1, set_phi_2;
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						PhaseFieldPoint& point = simulation_mesh::phase_field(i + 1, j + 1, 0 + 1);
						REAL sum_phi = 0.0;
						for (int index = 0; index < model_parameters::phi_number; index++)
							for (int index2 = 0; index2 < microstructure_init::porous_phis_indexs.size(); index2++)
								if (index == microstructure_init::porous_phis_indexs[index2] && point.phi[index] > SYS_EPSILON) {
									sum_phi += point.phi[index];
									point.phi[index] = 0.0;
								}
						if (sum_phi > SYS_EPSILON) {
							if (Solid[i][j]) {
								set_phi_1.add_point(i + 1, j + 1, 0 + 1, sum_phi);
							}
							else {
								set_phi_2.add_point(i + 1, j + 1, 0 + 1, sum_phi);
							}
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = microstructure_init::porous_first_phi_index;
				set_phi_1.temperature = microstructure_init::porous_first_temperature;
				set_phi_1.con = microstructure_init::porous_first_con;
				set_phi_1.is_normalized = microstructure_init::is_porous_normalized;
				microstructure_init::nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = microstructure_init::porous_second_phi_index;
				set_phi_2.temperature = microstructure_init::porous_second_temperature;
				set_phi_2.con = microstructure_init::porous_second_con;
				set_phi_2.is_normalized = microstructure_init::is_porous_normalized;
				microstructure_init::nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else if (mesh_parameters::dimention == Dimension::Three_Dimension) {
				vector<vector<vector<int>>> arrgrid; arrgrid.resize(mesh_parameters::MESH_NX);
				vector<vector<int>> soild;
				int numsoild = 0, Nxyz = mesh_parameters::MESH_NX * mesh_parameters::MESH_NY * mesh_parameters::MESH_NZ, numtotal_need = int(microstructure_init::porosity * Nxyz);
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					arrgrid[i].resize(mesh_parameters::MESH_NY);
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						arrgrid[i][j].resize(mesh_parameters::MESH_NZ);
						for (int k = 0; k < mesh_parameters::MESH_NZ; k++) {
							arrgrid[i][j][k] = 0;
							if (real_dist(gen) < microstructure_init::porous_init_noise) {
								arrgrid[i][j][k] = 1;
								vector<int> buff = { i, j, k };
								soild.push_back(buff);
								numsoild = numsoild + 1;
							}
						}
					}
				}
				REAL core_p = REAL(numsoild) / REAL(Nxyz);
				string out_ptr = "  percentage of core in 3D = " + to_string(core_p) + "\n";
				WriteLog(out_ptr);
				int Tnumsoild = numsoild, step = 0;
				while (Tnumsoild < numtotal_need) {
					step++;
					quartet_structure_grow_3D(arrgrid, soild, Tnumsoild, gen);
					core_p = REAL(Tnumsoild) / REAL(Nxyz);
					out_ptr = "  grow step: " + to_string(step) + ", percentage of particle/pore = " + to_string(core_p) + "\n";
					WriteLog(out_ptr);
				}
				out_ptr = "  End of growth ! \n";
				WriteLog(out_ptr);
				PointSet set_phi_1, set_phi_2;
				for (int i = 0; i < mesh_parameters::MESH_NX; i++) {
					for (int j = 0; j < mesh_parameters::MESH_NY; j++) {
						for (int k = 0; k < mesh_parameters::MESH_NZ; k++) {
							PhaseFieldPoint& point = simulation_mesh::phase_field(i + 1, j + 1, k + 1);
							REAL sum_phi = 0.0;
							for (int index = 0; index < model_parameters::phi_number; index++)
								for (int index2 = 0; index2 < microstructure_init::porous_phis_indexs.size(); index2++)
									if (index == microstructure_init::porous_phis_indexs[index2] && point.phi[index] > SYS_EPSILON) {
										sum_phi += point.phi[index];
										point.phi[index] = 0.0;
									}
							if (sum_phi > SYS_EPSILON) {
								if (arrgrid[i][j][k]) {
									set_phi_1.add_point(i + 1, j + 1, k + 1, sum_phi);
								}
								else {
									set_phi_2.add_point(i + 1, j + 1, k + 1, sum_phi);
								}
							}
						}
					}
				}
				set_phi_1.generate_step = 0;
				set_phi_1.phaseIndex = microstructure_init::porous_first_phi_index;
				set_phi_1.temperature = microstructure_init::porous_first_temperature;
				set_phi_1.con = microstructure_init::porous_first_con;
				set_phi_1.is_normalized = microstructure_init::is_porous_normalized;
				microstructure_init::nucleation_box.point_set_box.push_back(set_phi_1);

				set_phi_2.generate_step = 0;
				set_phi_2.phaseIndex = microstructure_init::porous_second_phi_index;
				set_phi_2.temperature = microstructure_init::porous_second_temperature;
				set_phi_2.con = microstructure_init::porous_second_con;
				set_phi_2.is_normalized = microstructure_init::is_porous_normalized;
				microstructure_init::nucleation_box.point_set_box.push_back(set_phi_2);
			}
			else {
				WriteLog("> Quartet structure generation Failed, can't generate at one dimension ! \n");
			}
		}

	}
}