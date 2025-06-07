#include "VoronoiStructure.h"
namespace pf {
	namespace voronoi_structure {

		inline REAL voronoi_points_distance(Point point) {
			using namespace microstructure_init;
			REAL points_distance = -1;
			if (voronoi_type == VoronoiType::VDT_CONST) {
				points_distance = voronoi_const_pointsDistance;
			}
			else if (voronoi_type == VoronoiType::VDT_REF_DOT) {
				REAL pos = voronoi_reference_dot.to_length(point.x, point.y, point.z);
				if (pos > voronoi_reference_dot_distance) {
					points_distance = voronoi_reference_dot_max_pointsDistance;
				}
				else {
					points_distance = pos / voronoi_reference_dot_distance * voronoi_reference_dot_max_pointsDistance + (1 - pos / voronoi_reference_dot_distance) * voronoi_reference_dot_min_pointsDistance;
				}
			}
			else if (voronoi_type == VoronoiType::VDT_REF_SURFACE) {
				REAL pos = voronoi_reference_surface.distance(point.x, point.y, point.z);
				if (pos > voronoi_reference_surface_distance) {
					points_distance = voronoi_reference_surface_max_pointsDistance;
				}
				else {
					points_distance = pos / voronoi_reference_surface_distance * voronoi_reference_surface_max_pointsDistance + (1 - pos / voronoi_reference_surface_distance) * voronoi_reference_surface_min_pointsDistance;
				}
			}
			else if (voronoi_type == VoronoiType::VDT_DOTS_MATRIX) {
				vector<REAL> dots_lenghts;
				points_distance = 0;
				REAL total_lenght = 0;
				for (auto dot = voronoi_matrix_dots.begin(); dot < voronoi_matrix_dots.end(); dot++) {
					dots_lenghts.push_back(dot->to_length(point.x, point.y, point.z));
					total_lenght += dot->to_length(point.x, point.y, point.z);
				}
				total_lenght = total_lenght * REAL(voronoi_matrix_dots.size() - 1);
				if (total_lenght < SYS_EPSILON) {
					WriteDebugFile("ERROR: VDT_DOTS_MATRIX method need more than one dot ! \n");
					exit(0);
				}
				for (auto index = 0; index < voronoi_matrix_dots.size(); index++) {
					REAL potential = 0;
					for (auto index2 = 0; index2 < voronoi_matrix_dots.size(); index2++)
						if (index != index2)
							potential += dots_lenghts[index2];
					potential /= total_lenght;
					points_distance += potential * voronoi_matrix_dots_pointsDistance[index];
				}
			}
			return points_distance;
		}

		void init() {
			using namespace microstructure_init;
			WriteDebugFile("# .voronoi = [(phi_index_begin, phi_index_end), ( box_origin_point ),( box_end_point )] \n");
			string voronoi_property_key = "Preprocess.Microstructure.voronoi", voronoi_property_input = "[()]";
			if (InputFileReader::get_instance()->read_string_value(voronoi_property_key, voronoi_property_input, true)) {
				is_voronoi = true;
				vector<InputValueType> voronoi_property_structure = vector<InputValueType>{ InputValueType::IVType_INT, InputValueType::IVType_REAL, InputValueType::IVType_REAL };
				vector<vector<input_value>> voronoi_value = InputFileReader::get_instance()->trans_matrix_2d_array_const_to_input_value(voronoi_property_structure, voronoi_property_key, voronoi_property_input, true);
				voronoi_phi_index_range[0] = voronoi_value[0][0].int_value;
				voronoi_phi_index_range[1] = voronoi_value[0][1].int_value;
				check_phi_index(voronoi_phi_index_range[0]);
				check_phi_index(voronoi_phi_index_range[1]);
				voronoi_box_position[0] = voronoi_value[1][0].REAL_value;
				voronoi_box_position[1] = voronoi_value[1][1].REAL_value;
				voronoi_box_position[2] = voronoi_value[1][2].REAL_value;
				voronoi_box_size[0] = voronoi_value[2][0].REAL_value;
				voronoi_box_size[1] = voronoi_value[2][1].REAL_value;
				voronoi_box_size[2] = voronoi_value[2][2].REAL_value;

				InputFileReader::get_instance()->read_bool_value("Preprocess.Microstructure.Voronoi.is_box_periodic", is_voronoi_mirror_generation, true);

				WriteDebugFile("# .con = (comp_0_value, comp_1_value, ... )] \n");
				string voronoi_x_key = "Preprocess.Microstructure.Voronoi.con", voronoi_x_input = "()";
				InputFileReader::get_instance()->read_string_value(voronoi_x_key, voronoi_x_input, true);
				vector<input_value> voronoi_x_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, voronoi_x_key, voronoi_x_input, true);
				for (int x_index = 0; x_index < voronoi_x_value.size(); x_index++)
					voronoi_con.push_back(voronoi_x_value[x_index].REAL_value);
				check_con_size(int(voronoi_con.size()));

				string voronoi_temperature_key = "Preprocess.Microstructure.Voronoi.temperature";
				InputFileReader::get_instance()->read_REAL_value(voronoi_temperature_key, voronoi_temperature, true);
				if (InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.Voronoi.rand_seed", voronoi_rand_seed, true))
					is_voronoi_rand = false;

				WriteDebugFile("# Preprocess.Microstructure.Voronoi.const_distance = 0.0 \n");
				WriteDebugFile("#                                  .ref_dot = [(REF_DOT), (distance, min_points_distance, max_points_distance)] \n");
				WriteDebugFile("#                                  .ref_surface = [(REF_SURF_POINT_1), (REF_SURF_POINT_2), (REF_SURF_POINT_3), (distance, min_points_distance, max_points_distance)] \n");
				WriteDebugFile("#                                  .dots_matrix = [(DOT1, points_distances1), (DOT2, points_distances2), ... ] \n");
				string voronoi_const_distance_key = "Preprocess.Microstructure.Voronoi.const_distance";
				string voronoi_ref_dot_key = "Preprocess.Microstructure.Voronoi.ref_dot", voronoi_ref_dot_input = "[()]";
				string voronoi_ref_surface_key = "Preprocess.Microstructure.Voronoi.ref_surface", voronoi_ref_surface_input = "[()]";
				string voronoi_dots_matrix_key = "Preprocess.Microstructure.Voronoi.dots_matrix", voronoi_dots_matrix_input = "[()]";
				if (InputFileReader::get_instance()->read_REAL_value(voronoi_const_distance_key, voronoi_const_pointsDistance, true)) {
					voronoi_type = VoronoiType::VDT_CONST;
				}
				else if (InputFileReader::get_instance()->read_string_value(voronoi_ref_dot_key, voronoi_ref_dot_input, true)) {
					voronoi_type = VoronoiType::VDT_REF_DOT;
					vector<vector<input_value>> voronoi_ref_dot_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, voronoi_ref_dot_key, voronoi_ref_dot_input, true);
					voronoi_reference_dot.set(voronoi_ref_dot_value[0][0].REAL_value, voronoi_ref_dot_value[0][1].REAL_value, voronoi_ref_dot_value[0][2].REAL_value);
					voronoi_reference_dot_distance = voronoi_ref_dot_value[1][0].REAL_value;
					voronoi_reference_dot_min_pointsDistance = voronoi_ref_dot_value[1][1].REAL_value;
					voronoi_reference_dot_max_pointsDistance = voronoi_ref_dot_value[1][2].REAL_value;
				}
				else if (InputFileReader::get_instance()->read_string_value(voronoi_ref_surface_key, voronoi_ref_surface_input, true)) {
					voronoi_type = VoronoiType::VDT_REF_SURFACE;
					vector<vector<input_value>> voronoi_ref_surface_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, voronoi_ref_surface_key, voronoi_ref_surface_input, true);
					voronoi_reference_surface.init(voronoi_ref_surface_value[0][0].REAL_value, voronoi_ref_surface_value[0][1].REAL_value, voronoi_ref_surface_value[0][2].REAL_value,
						voronoi_ref_surface_value[1][0].REAL_value, voronoi_ref_surface_value[1][1].REAL_value, voronoi_ref_surface_value[1][2].REAL_value,
						voronoi_ref_surface_value[2][0].REAL_value, voronoi_ref_surface_value[2][1].REAL_value, voronoi_ref_surface_value[2][2].REAL_value);
					voronoi_reference_surface_distance = voronoi_ref_surface_value[3][0].REAL_value;
					voronoi_reference_surface_min_pointsDistance = voronoi_ref_surface_value[3][1].REAL_value;
					voronoi_reference_surface_max_pointsDistance = voronoi_ref_surface_value[3][2].REAL_value;
				}
				else if (InputFileReader::get_instance()->read_string_value(voronoi_dots_matrix_key, voronoi_dots_matrix_input, true)) {
					voronoi_type = VoronoiType::VDT_DOTS_MATRIX;
					vector<vector<input_value>> voronoi_dots_matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, voronoi_dots_matrix_key, voronoi_dots_matrix_input, true);
					if (voronoi_dots_matrix_value.size() < 2) {
						WriteDebugFile(" ERROR: Preprocess.Microstructure.Voronoi.dots_matrix, the number of DOTS should be larger than 1. \n");
						exit(0);
					}
					for (int index = 0; index < voronoi_dots_matrix_value.size(); index++) {
						voronoi_matrix_dots.push_back(Point(voronoi_dots_matrix_value[index][0].REAL_value, voronoi_dots_matrix_value[index][1].REAL_value, voronoi_dots_matrix_value[index][2].REAL_value));
						voronoi_matrix_dots_pointsDistance.push_back(voronoi_dots_matrix_value[index][3].REAL_value);
					}
				}
				WriteDebugFile("# .in_phi_indexs = ( phi_index_1, phi_index_2, ... ) \n");
				string voronoi_in_phis_key = "Preprocess.Microstructure.Voronoi.in_phi_indexs", voronoi_in_phis_input = "()";
				InputFileReader::get_instance()->read_string_value(voronoi_in_phis_key, voronoi_in_phis_input, true);
				vector<input_value> voronoi_in_phis_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, voronoi_in_phis_key, voronoi_in_phis_input, true);
				for (int in_phi_index = 0; in_phi_index < voronoi_in_phis_value.size(); in_phi_index++)
					voronoi_phis_indexs.push_back(voronoi_in_phis_value[in_phi_index].int_value);

				bool is_debug = false;
				InputFileReader::get_instance()->read_bool_value("Preprocess.Microstructure.Voronoi.VTS.debug", is_debug, true);
				if (is_debug)
					pf::load_vts_scalar_func(write_scalar_voronoi_points);
			}
		}

		void generate_voronoi_structure() {
			Dimension dimention = Dimension::Three_Dimension;
			if (microstructure_init::voronoi_box_size[0] == 0 && microstructure_init::voronoi_box_size[1] == 0 && microstructure_init::voronoi_box_size[2] == 0)
				return;
			else if (microstructure_init::voronoi_box_size[0] == 0 || microstructure_init::voronoi_box_size[1] == 0 || microstructure_init::voronoi_box_size[2] == 0)
				dimention = Dimension::Two_Dimension;
			int grain_number = microstructure_init::voronoi_phi_index_range[1] - microstructure_init::voronoi_phi_index_range[0] + 1;
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!microstructure_init::is_voronoi_rand) {
				gen.seed(microstructure_init::voronoi_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			int grain_index = 0;
			while (grain_index < grain_number) {
				bool is_point_add = true;
				REAL rand_x = REAL(real_dist(gen)), rand_y = REAL(real_dist(gen)), rand_z = REAL(real_dist(gen));
				Point p(rand_x * microstructure_init::voronoi_box_size[0] + microstructure_init::voronoi_box_position[0],
					rand_y * microstructure_init::voronoi_box_size[1] + microstructure_init::voronoi_box_position[1],
					rand_z * microstructure_init::voronoi_box_size[2] + microstructure_init::voronoi_box_position[2]);
				for (auto ip = microstructure_init::voronoi_points.begin(); ip < microstructure_init::voronoi_points.end(); ip++) {
					REAL D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
					REAL distance = voronoi_points_distance(Point((p.x + ip->x) / 2, (p.y + ip->y) / 2, (p.z + ip->z) / 2));
					if (distance < 0)
						distance = 0;
					if (D2 < (distance * distance)) {
						is_point_add = false;
					}
				}
				if (is_point_add) {
					microstructure_init::voronoi_points.push_back(p);
					grain_index++;
					string str_report = "> Voronoi: generate point index : " + to_string(grain_index) + ", at : (x, y, z) (" + to_string(p.x) + ", " + to_string(p.y) + ", " + to_string(p.z) + ")\n";
					WriteLog(str_report);
				}
			}
			// > periodic boundary condition
			vector<vector<Point>> mirror_points; // = 27 * points.size()
			int region_number = 0;
			if (dimention == Dimension::Three_Dimension)
				region_number = 27;
			else
				region_number = 9;
			if (not microstructure_init::is_voronoi_mirror_generation)
				region_number = 1;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				mirror_points[0][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, 0);
				if (microstructure_init::is_voronoi_mirror_generation) {
					using namespace microstructure_init;
					mirror_points[1][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, 0);
					mirror_points[2][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, 0);
					mirror_points[3][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), 0);
					mirror_points[4][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), 0);
					mirror_points[5][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[6][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[7][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					mirror_points[8][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					if (dimention == Dimension::Three_Dimension) {
						mirror_points[9][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, int(voronoi_box_size[2]));
						mirror_points[10][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, int(-voronoi_box_size[2]));
						mirror_points[11][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[12][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[13][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[14][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[15][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[16][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[17][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[18][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[19][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[20][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[23][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[24][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[26][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[21][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[22][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[25][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
					}
				}
			}
			for (int region = 0; region < region_number; region++)
				for (int grain = 0; grain < grain_number; grain++) {
					Polyhedron poly(mirror_points[region][grain]);
					vector<point_in_region_index> record_points;
					record_points.push_back(point_in_region_index(region, grain));
					for (unsigned int region_index = 0; region_index < mirror_points.size(); region_index++)
						for (unsigned int grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
							// Avoid inclusion points
							bool is_point_contained = false;
							for (auto re = record_points.begin(); re < record_points.end(); re++)
								if (re->region == region_index && re->grain_index == grain_index)
									is_point_contained = true;
							if (is_point_contained)
								continue;
							// prepare
							Vector3 norm = get_vector(mirror_points[region_index][grain_index], poly.point_inside_polyhedron);
							// Vector3 norm = get_vector(poly.point_inside_polyhedron, mirror_points[region_index][grain_index]);
							Point mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
							// Judged ipsilateral to poly center
							if (poly.check_point_inside_polyhedron(mid_point) == false)
								continue;
							// < add point
							poly.add_surf(norm, mid_point);
							// < Eliminate meaningless points in poly
							for (auto re = record_points.begin(); re < record_points.end();) {
								Point check = (microstructure_init::voronoi_points[grain] + mirror_points[re->region][re->grain_index]) / 2;
								if (poly.check_point_inside_polyhedron(check) == false) {
									re = record_points.erase(re);
								}
								else {
									++re;
								}
							}
						}
					string str_report = "> Voronoi: One polyhedron in region : " + to_string(region) + ", grain : " 
						+ to_string(grain) + " has been generated ! \n";
					WriteLog(str_report);
					GeometricRegion geo;
					geo.geometryProperty = Geometry::Geo_Polyhedron;
					geo.generate_step = 0;
					geo.polyhedron = poly;
					geo.phaseIndex = microstructure_init::voronoi_phi_index_range[0] + grain;
					geo.temperature = microstructure_init::voronoi_temperature;
					geo.con = microstructure_init::voronoi_con;
					geo.phi = 1;
					geo.isNormalized = true;
					microstructure_init::nucleation_box.geometry_box.push_back(geo);
				}

		}

		void generate_voronoi_structure_new_method() {
			Dimension dimention = Dimension::Three_Dimension;
			if (microstructure_init::voronoi_box_size[0] == 0 && microstructure_init::voronoi_box_size[1] == 0 && microstructure_init::voronoi_box_size[2] == 0)
				return;
			else if (microstructure_init::voronoi_box_size[0] == 0 || microstructure_init::voronoi_box_size[1] == 0 || microstructure_init::voronoi_box_size[2] == 0)
				dimention = Dimension::Two_Dimension;
			int grain_number = microstructure_init::voronoi_phi_index_range[1] - microstructure_init::voronoi_phi_index_range[0] + 1;
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!microstructure_init::is_voronoi_rand) {
				gen.seed(microstructure_init::voronoi_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			int grain_index = 0;
			while (grain_index < grain_number) {
				bool is_point_add = true;
				REAL rand_x = REAL(real_dist(gen)), rand_y = REAL(real_dist(gen)), rand_z = REAL(real_dist(gen));
				Point p(rand_x * microstructure_init::voronoi_box_size[0] + microstructure_init::voronoi_box_position[0],
					rand_y * microstructure_init::voronoi_box_size[1] + microstructure_init::voronoi_box_position[1],
					rand_z * microstructure_init::voronoi_box_size[2] + microstructure_init::voronoi_box_position[2]);
				for (auto ip = microstructure_init::voronoi_points.begin(); ip < microstructure_init::voronoi_points.end(); ip++) {
					REAL D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
					REAL distance = voronoi_points_distance(Point((p.x + ip->x) / 2, (p.y + ip->y) / 2, (p.z + ip->z) / 2));
					if (distance < 0)
						distance = 0;
					if (D2 < (distance * distance)) {
						is_point_add = false;
					}
				}
				if (is_point_add) {
					microstructure_init::voronoi_points.push_back(p);
					grain_index++;
					string str_report = "> Voronoi: generate point index : " + to_string(grain_index) + ", at : (x, y, z) (" + to_string(p.x) + ", " + to_string(p.y) + ", " + to_string(p.z) + ")\n";
					WriteLog(str_report);
				}
			}
			// > periodic boundary condition
			vector<vector<Point>> mirror_points; // = 27 * points.size()
			int region_number = 0;
			if (dimention == Dimension::Three_Dimension)
				region_number = 27;
			else
				region_number = 9;
			if (not microstructure_init::is_voronoi_mirror_generation)
				region_number = 1;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				mirror_points[0][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, 0);
				if (microstructure_init::is_voronoi_mirror_generation) {
					using namespace microstructure_init;
					mirror_points[1][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, 0);
					mirror_points[2][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, 0);
					mirror_points[3][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), 0);
					mirror_points[4][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), 0);
					mirror_points[5][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[6][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[7][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					mirror_points[8][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					if (dimention == Dimension::Three_Dimension) {
						mirror_points[9][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, int(voronoi_box_size[2]));
						mirror_points[10][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, int(-voronoi_box_size[2]));
						mirror_points[11][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[12][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[13][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[14][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[15][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[16][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[17][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[18][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[19][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[20][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[23][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[24][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[26][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[21][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[22][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[25][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
					}
				}
			}
			for (int region = 0; region < region_number; region++)
				for (int grain = 0; grain < grain_number; grain++) {
					Polyhedron poly(mirror_points[region][grain]);
					vector<point_in_region_index> record_points;
					record_points.push_back(point_in_region_index(region, grain));
					bool find_new_point = false;
					do
					{
						find_new_point = false;
						REAL point_distance2 = REAL_MAX();
						int buff_region = 0, buff_grain = 0;
						for (unsigned int region_index = 0; region_index < mirror_points.size(); region_index++)
							for (unsigned int grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
								// Avoid inclusion points
								bool is_point_contained = false;
								for (auto re = record_points.begin(); re < record_points.end(); re++)
									if (re->region == region_index && re->grain_index == grain_index)
										is_point_contained = true;
								if (is_point_contained)
									continue;
								// check mid point
								Point mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
								// Judged ipsilateral to poly center
								if (poly.check_point_inside_polyhedron(mid_point) == false)
									continue;
								find_new_point = true;
								REAL distance2 = (poly.point_inside_polyhedron.x - mirror_points[region_index][grain_index].x) * (poly.point_inside_polyhedron.x - mirror_points[region_index][grain_index].x)
									+ (poly.point_inside_polyhedron.y - mirror_points[region_index][grain_index].y) * (poly.point_inside_polyhedron.y - mirror_points[region_index][grain_index].y)
									+ (poly.point_inside_polyhedron.z - mirror_points[region_index][grain_index].z) * (poly.point_inside_polyhedron.z - mirror_points[region_index][grain_index].z);
								if (distance2 < point_distance2) {
									point_distance2 = distance2;
									buff_region = region_index;
									buff_grain = grain_index;
								}
							}
						if (find_new_point) {
							Point mid_point = (poly.point_inside_polyhedron + mirror_points[buff_region][buff_grain]) / 2;
							Vector3 norm = get_vector(mirror_points[buff_region][buff_grain], poly.point_inside_polyhedron);
							poly.add_surf(norm, mid_point);
							record_points.push_back(point_in_region_index(buff_region, buff_grain));
						}
					} while (find_new_point);
					string str_report = "> Voronoi: One polyhedron in region : " + to_string(region) + ", grain : "
						+ to_string(grain) + " has been generated ! \n";
					WriteLog(str_report);
					GeometricRegion geo;
					geo.geometryProperty = Geometry::Geo_Polyhedron;
					geo.generate_step = 0;
					geo.polyhedron = poly;
					geo.phaseIndex = microstructure_init::voronoi_phi_index_range[0] + grain;
					geo.temperature = microstructure_init::voronoi_temperature;
					geo.con = microstructure_init::voronoi_con;
					geo.phi = 1;
					geo.isNormalized = true;
					microstructure_init::nucleation_box.geometry_box.push_back(geo);
				}

		}

		void generate_voronoi_structure_in_phis() {
			Dimension dimention = Dimension::Three_Dimension;
			if (microstructure_init::voronoi_box_size[0] == 0 && microstructure_init::voronoi_box_size[1] == 0 && microstructure_init::voronoi_box_size[2] == 0)
				return;
			else if (microstructure_init::voronoi_box_size[0] == 0 || microstructure_init::voronoi_box_size[1] == 0 || microstructure_init::voronoi_box_size[2] == 0)
				dimention = Dimension::Two_Dimension;
			int grain_number = microstructure_init::voronoi_phi_index_range[1] - microstructure_init::voronoi_phi_index_range[0] + 1;
			// > generate points
			std::random_device rd; // 高质量随机数种子生成器
			// std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
			std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
			if (!microstructure_init::is_voronoi_rand) {
				gen.seed(microstructure_init::voronoi_rand_seed);
			}
			std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
			// std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
			// std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
			// REAL rand = real_dist(gen);  // random
			int grain = 0;
			while (grain < grain_number) {
				bool is_point_add = true;
				REAL rand_x = REAL(real_dist(gen)), rand_y = REAL(real_dist(gen)), rand_z = REAL(real_dist(gen));
				Point p(rand_x * microstructure_init::voronoi_box_size[0] + microstructure_init::voronoi_box_position[0],
					rand_y * microstructure_init::voronoi_box_size[1] + microstructure_init::voronoi_box_position[1],
					rand_z * microstructure_init::voronoi_box_size[2] + microstructure_init::voronoi_box_position[2]);
				PhaseFieldPoint& point = phi_parameters::phase_field(int(p.x), int(p.y), int(p.z));
				REAL sum_phi = 0;
				for (int index = 0; index < phi_parameters::phi_number; index++)
					for (int index2 = 0; index2 < microstructure_init::voronoi_phis_indexs.size(); index2++)
						if (index == microstructure_init::voronoi_phis_indexs[index2])
							sum_phi += point.phi[index];
				for (auto ip = microstructure_init::voronoi_points.begin(); ip < microstructure_init::voronoi_points.end(); ip++) {
					REAL D2 = (p.x - ip->x) * (p.x - ip->x) + (p.y - ip->y) * (p.y - ip->y) + (p.z - ip->z) * (p.z - ip->z);
					REAL distance = voronoi_points_distance(Point((p.x + ip->x) / 2, (p.y + ip->y) / 2, (p.z + ip->z) / 2));
					if (distance < 0)
						distance = 0;
					if (D2 < (distance * distance)) {
						is_point_add = false;
					}
				}
				if (is_point_add && sum_phi > SYS_EPSILON) {
					grain++;
					microstructure_init::voronoi_points.push_back(p);
					string str_report = "> Voronoi: generate point index : " + to_string(grain) + ",  at : (x, y, z) (" + to_string(p.x) + ", " + to_string(p.y) + ", " + to_string(p.z) + ")\n";
					WriteLog(str_report);
				}
			}
			// > periodic boundary condition
			vector<vector<Point>> mirror_points; // = 27 * points.size()
			int region_number = 0;
			if (dimention == Dimension::Three_Dimension)
				region_number = 27;
			else
				region_number = 9;
			if (not microstructure_init::is_voronoi_mirror_generation)
				region_number = 1;
			mirror_points.resize(region_number);
			for (auto region = mirror_points.begin(); region < mirror_points.end(); region++)
				region->resize(grain_number);
#pragma omp parallel for
			for (int grain = 0; grain < grain_number; grain++) {
				using namespace microstructure_init;
				mirror_points[0][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, 0);
				if (microstructure_init::is_voronoi_mirror_generation) {
					mirror_points[1][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, 0);
					mirror_points[2][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, 0);
					mirror_points[3][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), 0);
					mirror_points[4][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), 0);
					mirror_points[5][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[6][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), 0);
					mirror_points[7][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					mirror_points[8][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), 0);
					if (dimention == Dimension::Three_Dimension) {
						mirror_points[9][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, int(voronoi_box_size[2]));
						mirror_points[10][grain] = microstructure_init::voronoi_points[grain] + Point(0, 0, int(-voronoi_box_size[2]));
						mirror_points[11][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[12][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, int(voronoi_box_size[2]));
						mirror_points[13][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[14][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), 0, int(-voronoi_box_size[2]));
						mirror_points[15][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[16][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[17][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[18][grain] = microstructure_init::voronoi_points[grain] + Point(0, int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[19][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[20][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[23][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[24][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[26][grain] = microstructure_init::voronoi_points[grain] + Point(int(-voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[21][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(voronoi_box_size[2]));
						mirror_points[22][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(voronoi_box_size[1]), int(-voronoi_box_size[2]));
						mirror_points[25][grain] = microstructure_init::voronoi_points[grain] + Point(int(voronoi_box_size[0]), int(-voronoi_box_size[1]), int(-voronoi_box_size[2]));
					}
				}
			}
			for (int region = 0; region < region_number; region++)
				for (int grain = 0; grain < grain_number; grain++) {
					Polyhedron poly(mirror_points[region][grain]);
					vector<point_in_region_index> record_points;
					point_in_region_index rp(region, grain);
					record_points.push_back(rp);
					for (unsigned int region_index = 0; region_index < mirror_points.size(); region_index++)
						for (unsigned int grain_index = 0; grain_index < mirror_points[region_index].size(); grain_index++) {
							// Avoid inclusion points
							bool is_point_contained = false;
							for (auto re = record_points.begin(); re < record_points.end(); re++)
								if (re->region == region_index && re->grain_index == grain_index)
									is_point_contained = true;
							if (is_point_contained)
								continue;
							// prepare
							// Vector3 norm = get_vector(mirror_points[region_index][grain_index], poly.point_inside_polyhedron);
							Vector3 norm = get_vector(poly.point_inside_polyhedron, mirror_points[region_index][grain_index]);
							Point mid_point = (poly.point_inside_polyhedron + mirror_points[region_index][grain_index]) / 2;
							// Judged ipsilateral to poly center
							if (poly.check_point_inside_polyhedron(mid_point) == false)
								continue;
							// < add point
							poly.add_surf(norm, mid_point);
							// < Eliminate meaningless points in poly
							for (auto re = record_points.begin(); re < record_points.end();) {
								Point check = (microstructure_init::voronoi_points[grain] + mirror_points[re->region][re->grain_index]) / 2;
								if (poly.check_point_inside_polyhedron(check) == false) {
									re = record_points.erase(re);
								}
								else {
									++re;
								}
							}
						}
					PointSet set;
					REAL sum_sum_phi = 0;
					for (int z = 0; z < phi_parameters::phase_field.Nz(); z++)
						for (int y = 0; y < phi_parameters::phase_field.Ny(); y++)
							for (int x = 0; x < phi_parameters::phase_field.Nx(); x++) {
								PhaseFieldPoint& point = phi_parameters::phase_field(x, y, z);
								if (poly.check_point_inside_polyhedron(pf::Point(x, y, z)) == true) {
									REAL sum_phi = 0;
									for (int index = 0; index < phi_parameters::phi_number; index++)
										for (int index2 = 0; index2 < microstructure_init::voronoi_phis_indexs.size(); index2++)
											if (index == microstructure_init::voronoi_phis_indexs[index2] && point.phi[index] > SYS_EPSILON) {
												sum_phi += point.phi[index];
												point.phi[index] = 0;
											}
									if (sum_phi > SYS_EPSILON) {
										set.add_point(x, y, z, sum_phi);
										sum_sum_phi += sum_phi;
									}
								}
							}
					string str_report = "> Voronoi: One polyhedron in region : " + to_string(region) + ", grain : " + to_string(grain) + " could been generated ! \n";
					WriteLog(str_report);
					set.generate_step = 0;
					set.phaseIndex = microstructure_init::voronoi_phi_index_range[0] + grain;
					set.temperature = microstructure_init::voronoi_temperature;
					set.con = microstructure_init::voronoi_con;
					set.is_normalized = false;
					if (sum_sum_phi > SYS_EPSILON) {
						microstructure_init::nucleation_box.point_set_box.push_back(set);
					}
				}
		}

	}
}