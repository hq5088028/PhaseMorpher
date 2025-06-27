#include "GeometryStructure.h"
namespace pf {
	namespace geometry_structure {

		std::vector<std::array<std::array<REAL, 3>, 3>> center_length_to_center_surfaces(const std::array<REAL, 3> center, const std::array<REAL, 3> lengths) {
			std::array<bool, 3> infinity_direction{ false, false, false };
			std::array<REAL, 3> half_lengths{};
			// cal half length
			for (int index{ 0 }; index < 3; index++) {
				if (lengths.at(index) < 0) {
					half_lengths.at(index) = lengths.at(index);
					std::cout << "encounter infinity length at direction " << index << "." << std::endl;
					infinity_direction.at(index) = true;
				}
				else if (lengths.at(index) == 0) {
					std::invalid_argument("length should not be zero.");
				}
				else {
					half_lengths.at(index) = lengths.at(index) / 2;
				}
			}
			std::vector<std::array<std::array<REAL, 3>, 3>> surfaces{};
			std::array<std::array<REAL, 3>, 3> single_surface{};

			for (int index{ 0 }; index < 3; index++) {
				std::array<REAL, 3> single_point{ center };
				if (not infinity_direction.at(index)) {
					single_point.at(index) = center.at(index) - half_lengths.at(index); // create point_1 with coordinate 1
					single_point.at((index + 1) % 3) = center.at((index + 1) % 3);      // create point_1 with coordinate 2
					single_point.at((index + 2) % 3) = center.at((index + 2) % 3);      // create point_1 with coordinate 3
					single_surface.at(index) = single_point;                            // point_1 save
					// overwriiting point
					single_point.at((index + 1) % 3) += 1;             // point_2 via point_1 shift along direction 2
					single_surface.at((index + 1) % 3) = single_point; // point_2 save
					// overwriiting point
					single_point.at((index + 2) % 3) += 1;             // point_3 via point_2 shift along direction 3
					single_surface.at((index + 2) % 3) = single_point; // point_3 save
					// surface + finished
					surfaces.push_back(single_surface);

					single_point.at(index) += lengths.at(index); // 1_point via point_3 shift along direction -1
					single_surface.at(index) = single_point;     // 1_point save
					// overwriiting point
					single_point.at((index + 1) % 3) -= 1;             // 2_point via 1_point shift along direction 2
					single_surface.at((index + 1) % 3) = single_point; // 2_point save
					// overwriiting point
					single_point.at((index + 2) % 3) -= 1;             // 3_point shift along direction 3
					single_surface.at((index + 2) % 3) = single_point; // 3_point save
					// surface - finished
					surfaces.push_back(single_surface);
				}
			}
			return surfaces;
		}

		void init() {
			// geometry structure
			int nucleation_layer = 0;
			vector<int> geometry_layer;
			InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.geometry_layer_number", nucleation_layer, true);
			for (int layer_index = 0; layer_index < nucleation_layer; layer_index++) {
				GeometricRegion geo;
				if (layer_index == 0) {
					WriteDebugFile("# .property = (phi_index, geometry_type, rotation_gauge, reverse_region) \n");
					WriteDebugFile("#              geometry_type  : 0 - None, 1 - Ellipsoid, 2 - Polyhedron, 3 - Geo_SegmentedCylinder, 4 - RectangularCuboid \n");
					WriteDebugFile("#              rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n");
					WriteDebugFile("#                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n");
					WriteDebugFile("#  The origin of the simulation mesh is ( x = 1 , y = 1 , z = 1 ). \n");
				}
				string layer_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".property";
				string layer_input_property = "(0,0,1,false)";
				if (InputFileReader::get_instance()->read_string_value(layer_key, layer_input_property, true)) {
					vector<InputValueType> layer_input_property_structure; layer_input_property_structure.push_back(InputValueType::IVType_INT); layer_input_property_structure.push_back(InputValueType::IVType_INT); 
					layer_input_property_structure.push_back(InputValueType::IVType_INT); layer_input_property_structure.push_back(InputValueType::IVType_BOOL);
					vector<input_value> layer_input_property_value = InputFileReader::get_instance()->trans_matrix_1d_array_to_input_value(layer_input_property_structure, layer_key, layer_input_property, true);
					check_phi_index(layer_input_property_value[0].int_value);
					geo.init(pf::Geometry(layer_input_property_value[1].int_value), 0, layer_input_property_value[0].int_value, layer_input_property_value[3].bool_value);
					pf::RotationGauge rotation_gauge = pf::RotationGauge(layer_input_property_value[2].int_value);
					if (geo.geometryProperty == pf::Geometry::Geo_Ellipsoid) {
						WriteDebugFile("# .ellipsoid = [(core_x,core_y,core_z),(radius_x,radius_y,radius_z),(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n");
						string layer_ellipsoid_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".ellipsoid";
						string layer_input_ellipsoid = "[(0,0,0),(0,0,0),(0,0,0)]";
						InputFileReader::get_instance()->read_string_value(layer_ellipsoid_key, layer_input_ellipsoid, true);
						vector<vector<input_value>> layer_input_ellipsoid_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, layer_ellipsoid_key, layer_input_ellipsoid, true);
						geo.ellipSolid.set_core(layer_input_ellipsoid_value[0][0].REAL_value, layer_input_ellipsoid_value[0][1].REAL_value, layer_input_ellipsoid_value[0][2].REAL_value);
						geo.ellipSolid.set_radius(layer_input_ellipsoid_value[1][0].REAL_value, layer_input_ellipsoid_value[1][1].REAL_value, layer_input_ellipsoid_value[1][2].REAL_value);
						REAL radian[] = { AngleToRadians(layer_input_ellipsoid_value[2][0].REAL_value), AngleToRadians(layer_input_ellipsoid_value[2][1].REAL_value), AngleToRadians(layer_input_ellipsoid_value[2][2].REAL_value) };
						geo.ellipSolid.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
					}
					if (geo.geometryProperty == pf::Geometry::Geo_Polyhedron) {
						WriteDebugFile("# .polyhedron = {[inside_point],[surf_point,surf_point,surf_point], .... ,[(rotation_angle_1,rotation_angle_2,rotation_angle_3)]} \n");
						WriteDebugFile("#                surf_point = (position_x,position_y,position_z) \n");
						string layer_polyhedron_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".polyhedron";
						string layer_input_polyhedron = "{[(-1,0,0)],[(0,0,0),(0,1,0),(0,0,1)],[(0,0,0)]}";
						InputFileReader::get_instance()->read_string_value(layer_polyhedron_key, layer_input_polyhedron, true);
						vector<vector<vector<input_value>>> layer_input_polyhedron_value = InputFileReader::get_instance()->trans_matrix_3d_const_const_const_to_input_value(InputValueType::IVType_REAL, layer_polyhedron_key, layer_input_polyhedron, true);
						geo.polyhedron.set_a_point_inside_polyhedron(pf::Point(layer_input_polyhedron_value[0][0][0].REAL_value, layer_input_polyhedron_value[0][0][1].REAL_value, layer_input_polyhedron_value[0][0][2].REAL_value));
						for (int surf_index = 1; surf_index < layer_input_polyhedron_value.size() - 1; surf_index++) {
							geo.polyhedron.add_surf(pf::Point(layer_input_polyhedron_value[surf_index][0][0].REAL_value, layer_input_polyhedron_value[surf_index][0][1].REAL_value, layer_input_polyhedron_value[surf_index][0][2].REAL_value),
								pf::Point(layer_input_polyhedron_value[surf_index][1][0].REAL_value, layer_input_polyhedron_value[surf_index][1][1].REAL_value, layer_input_polyhedron_value[surf_index][1][2].REAL_value),
								pf::Point(layer_input_polyhedron_value[surf_index][2][0].REAL_value, layer_input_polyhedron_value[surf_index][2][1].REAL_value, layer_input_polyhedron_value[surf_index][2][2].REAL_value));
						}
						int rIndex = int(layer_input_polyhedron_value.size() - 1);
						REAL radian[] = { AngleToRadians(layer_input_polyhedron_value[rIndex][0][0].REAL_value), AngleToRadians(layer_input_polyhedron_value[rIndex][0][1].REAL_value), AngleToRadians(layer_input_polyhedron_value[rIndex][0][2].REAL_value) };
						geo.polyhedron.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
					}
					if (geo.geometryProperty == pf::Geometry::Geo_SegmentedCylinder) {
						WriteDebugFile("# .segmented_cylinder = [(radius),(central_axis_point_x,central_axis_point_y,central_axis_point_z), ... at least two point ... ,(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n");
						string layer_cylindricity_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".segmented_cylinder";
						string layer_input_cylindricity = "[(0),(0,0,0),(0,0,0),(0,0,0)]";
						InputFileReader::get_instance()->read_string_value(layer_cylindricity_key, layer_input_cylindricity, true);
						vector<vector<input_value>> layer_input_cylindricity_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, layer_cylindricity_key, layer_input_cylindricity, true);
						int vec_size = int(layer_input_cylindricity_value.size());
						if (vec_size < 4) {
							WriteDebugFile("> ERROR .segmented_cylinder : the central axis line need at least two point ! \n");
							exit(0);
						}
						geo.cylinder.set_radius(layer_input_cylindricity_value[0][0].REAL_value);
						for (int index = 1; index < vec_size - 1; index++) {
							geo.cylinder.add_point(Point(layer_input_cylindricity_value[index][0].REAL_value, layer_input_cylindricity_value[index][1].REAL_value, layer_input_cylindricity_value[index][2].REAL_value));
						}
						REAL radian[] = { AngleToRadians(layer_input_cylindricity_value[vec_size - 1][0].REAL_value),
							AngleToRadians(layer_input_cylindricity_value[vec_size - 1][1].REAL_value),
							AngleToRadians(layer_input_cylindricity_value[vec_size - 1][2].REAL_value) };
						geo.cylinder.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
					}
					if (geo.geometryProperty == pf::Geometry::Geo_RectangularCuboid) {
						geo.geometryProperty = pf::Geometry::Geo_Polyhedron;
						WriteDebugFile("# .rectangular_cuboid = [(central_point),(x_length,y_length,z_length),(rotation_angle_1,rotation_angle_2,rotation_angle_3)] \n");
						string layer_polyhedron_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".rectangular_cuboid";
						string layer_input_polyhedron = "[(-1,-1,-1),(1,1,1),(0,0,0)]";
						InputFileReader::get_instance()->read_string_value(layer_polyhedron_key, layer_input_polyhedron, true);
						vector<vector<input_value>> layer_input_polyhedron_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_REAL, layer_polyhedron_key, layer_input_polyhedron, true);
						pf::Point central_point(layer_input_polyhedron_value[0][0].REAL_value, layer_input_polyhedron_value[0][1].REAL_value, layer_input_polyhedron_value[0][2].REAL_value);
						geo.polyhedron.set_a_point_inside_polyhedron(central_point);

						REAL x_length = layer_input_polyhedron_value[1][0].REAL_value + SYS_EPSILON,
							y_length = layer_input_polyhedron_value[1][1].REAL_value + SYS_EPSILON,
							z_length = layer_input_polyhedron_value[1][2].REAL_value + SYS_EPSILON;
						// - 1
						pf::Point point_on_surf_1(central_point.x - x_length / 2, central_point.y, central_point.z);
						Vector3 norm_1 = get_vector(central_point, point_on_surf_1);
						geo.polyhedron.add_surf(norm_1, point_on_surf_1);
						// - 2
						pf::Point point_on_surf_2(central_point.x + x_length / 2, central_point.y, central_point.z);
						Vector3 norm_2 = get_vector(central_point, point_on_surf_2);
						geo.polyhedron.add_surf(norm_2, point_on_surf_2);
						// - 3
						pf::Point point_on_surf_3(central_point.x, central_point.y - y_length / 2, central_point.z);
						Vector3 norm_3 = get_vector(central_point, point_on_surf_3);
						geo.polyhedron.add_surf(norm_3, point_on_surf_3);
						// - 4
						pf::Point point_on_surf_4(central_point.x, central_point.y + y_length / 2, central_point.z);
						Vector3 norm_4 = get_vector(central_point, point_on_surf_4);
						geo.polyhedron.add_surf(norm_4, point_on_surf_4);
						// - 5
						pf::Point point_on_surf_5(central_point.x, central_point.y, central_point.z - z_length / 2);
						Vector3 norm_5 = get_vector(central_point, point_on_surf_5);
						geo.polyhedron.add_surf(norm_5, point_on_surf_5);
						// - 6
						pf::Point point_on_surf_6(central_point.x, central_point.y, central_point.z + z_length / 2);
						Vector3 norm_6 = get_vector(central_point, point_on_surf_6);
						geo.polyhedron.add_surf(norm_6, point_on_surf_6);
						// -
						REAL radian[] = { AngleToRadians(layer_input_polyhedron_value[2][0].REAL_value), AngleToRadians(layer_input_polyhedron_value[2][1].REAL_value), AngleToRadians(layer_input_polyhedron_value[2][2].REAL_value) };
						geo.polyhedron.set_rotation_radian_and_rotation_gauge(radian, rotation_gauge);
					}
					if (phi_parameters::is_phi_field_on) {
						string layer_phi_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".phi";
						InputFileReader::get_instance()->read_REAL_value(layer_phi_key, geo.phi, true);
						string layer_norm_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".is_normalized";
						InputFileReader::get_instance()->read_bool_value(layer_norm_key, geo.isNormalized, true);
					}
					if (con_parameters::is_con_field_on) {
						WriteDebugFile("# .con = ( comp_0_value, comp_1_value, ... ) \n");
						string layer_x_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".con", layer_x_input = "()";
						geo.con.resize(0);
						if (InputFileReader::get_instance()->read_string_value(layer_x_key, layer_x_input, true)) {
							vector<input_value> layer_input_x_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_REAL, layer_x_key, layer_x_input, true);
							for (int x_index = 0; x_index < layer_input_x_value.size(); x_index++)
								geo.con.push_back(layer_input_x_value[x_index].REAL_value);
							check_con_size(int(geo.con.size()));
						}
					}
					if (temp_parameters::is_temp_field_on) {
						string layer_temp_key = "Preprocess.Microstructure.geometry_layer_" + to_string(layer_index) + ".temp";
						InputFileReader::get_instance()->read_REAL_value(layer_temp_key, geo.temperature, true);
					}
					geometry_layer.push_back(int(microstructure_init::nucleation_box.geometry_box.size()));
					microstructure_init::nucleation_box.geometry_box.push_back(geo);
				}
			}
			// batch geometry
			if (nucleation_layer > 0) {
				int homo_rotation_gauge = 1;
				WriteDebugFile("# .rotation_gauge : 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n");
				WriteDebugFile("#                   6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n");
				InputFileReader::get_instance()->read_int_value("Preprocess.Microstructure.GeometryLayerBatch.rotation_gauge", homo_rotation_gauge, true);
				WriteDebugFile("# .transformation = {[(trans_layer_index), (move_x, move_y, move_z), (rotation_angle_1,rotation_angle_2,rotation_angle_3), (rotation_angle_1,rotation_angle_2,rotation_angle_3,rotation_center_x,rotation_center_y,rotation_center_z)], ... } \n");
				string layer_transform_key = "Preprocess.Microstructure.GeometryLayerBatch.transformation", layer_transform_input = "{[(0),(0,0,0),(0,0,0),(0,0,0,0,0,0)]}";
				if (InputFileReader::get_instance()->read_string_value(layer_transform_key, layer_transform_input, true)) {
					vector<InputValueType> layer_transform_structure = { InputValueType::IVType_INT, InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL };
					vector<vector<vector<input_value>>> layer_transform_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(layer_transform_structure, layer_transform_key, layer_transform_input, true);
					for (int index = 0; index < layer_transform_value.size(); index++) {
						int geo_index = geometry_layer[layer_transform_value[index][0][0].int_value];
						REAL move_x = layer_transform_value[index][1][0].REAL_value, move_y = layer_transform_value[index][1][1].REAL_value, move_z = layer_transform_value[index][1][2].REAL_value,
							rotate_1 = AngleToRadians(layer_transform_value[index][2][0].REAL_value), rotate_2 = AngleToRadians(layer_transform_value[index][2][1].REAL_value), rotate_3 = AngleToRadians(layer_transform_value[index][2][2].REAL_value),
							Rotate_1 = AngleToRadians(layer_transform_value[index][3][0].REAL_value), Rotate_2 = AngleToRadians(layer_transform_value[index][3][1].REAL_value), Rotate_3 = AngleToRadians(layer_transform_value[index][3][2].REAL_value),
							Rotate_x = layer_transform_value[index][3][3].REAL_value, Rotate_y = layer_transform_value[index][3][4].REAL_value, Rotate_z = layer_transform_value[index][3][5].REAL_value;
						GeometricRegion& geo = microstructure_init::nucleation_box.geometry_box[geo_index];
						if (geo.geometryProperty == pf::Geometry::Geo_Ellipsoid) {
							geo.ellipSolid.move(move_x, move_y, move_z);
							geo.ellipSolid.add_rotation_radian(rotate_1, rotate_2, rotate_3);
							Vector3 pv(geo.ellipSolid.core.x - Rotate_x, geo.ellipSolid.core.y - Rotate_y, geo.ellipSolid.core.z - Rotate_z);
							Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
							geo.ellipSolid.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
						}
						else if (geo.geometryProperty == pf::Geometry::Geo_Polyhedron) {
							geo.polyhedron.move(move_x, move_y, move_z);
							geo.polyhedron.add_rotation_radian(rotate_1, rotate_2, rotate_3);
							Vector3 pv(geo.polyhedron.point_inside_polyhedron.x - Rotate_x, geo.polyhedron.point_inside_polyhedron.y - Rotate_y, geo.polyhedron.point_inside_polyhedron.z - Rotate_z);
							Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
							geo.polyhedron.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
						}
						else if (geo.geometryProperty == pf::Geometry::Geo_SegmentedCylinder) {
							geo.cylinder.move(move_x, move_y, move_z);
							geo.cylinder.add_rotation_radian(rotate_1, rotate_2, rotate_3);
							Vector3 pv(geo.cylinder.geometric_center.x - Rotate_x, geo.cylinder.geometric_center.y - Rotate_y, geo.cylinder.geometric_center.z - Rotate_z);
							Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
							geo.cylinder.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
						}
					}
				}
				WriteDebugFile("# .copy = {[(copy_layer_index), (paste_phi_index), (move_x, move_y, move_z), (rotation_angle_1,rotation_angle_2,rotation_angle_3), (rotation_angle_1,rotation_angle_2,rotation_angle_3,rotation_center_x,rotation_center_y,rotation_center_z)], ... } \n");
				string layer_copy_key = "Preprocess.Microstructure.GeometryLayerBatch.copy", layer_copy_input = "{[(0),(0),(0,0,0),(0,0,0),(0,0,0,0,0,0)]}";
				if (InputFileReader::get_instance()->read_string_value(layer_copy_key, layer_copy_input, true)) {
					vector<InputValueType> layer_copy_structure = { InputValueType::IVType_INT, InputValueType::IVType_INT, InputValueType::IVType_REAL, InputValueType::IVType_REAL, InputValueType::IVType_REAL };
					vector<vector<vector<input_value>>> layer_copy_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(layer_copy_structure, layer_copy_key, layer_copy_input, true);
					for (int index = 0; index < layer_copy_value.size(); index++) {
						int geo_index = geometry_layer[layer_copy_value[index][0][0].int_value], paste_phi_index = layer_copy_value[index][1][0].int_value;
						REAL move_x = layer_copy_value[index][2][0].REAL_value, move_y = layer_copy_value[index][2][1].REAL_value, move_z = layer_copy_value[index][2][2].REAL_value,
							rotate_1 = AngleToRadians(layer_copy_value[index][3][0].REAL_value), rotate_2 = AngleToRadians(layer_copy_value[index][3][1].REAL_value), rotate_3 = AngleToRadians(layer_copy_value[index][3][2].REAL_value),
							Rotate_1 = AngleToRadians(layer_copy_value[index][4][0].REAL_value), Rotate_2 = AngleToRadians(layer_copy_value[index][4][1].REAL_value), Rotate_3 = AngleToRadians(layer_copy_value[index][4][2].REAL_value),
							Rotate_x = layer_copy_value[index][4][3].REAL_value, Rotate_y = layer_copy_value[index][4][4].REAL_value, Rotate_z = layer_copy_value[index][4][5].REAL_value;
						GeometricRegion geo = microstructure_init::nucleation_box.geometry_box[geo_index];
						geo.phaseIndex = paste_phi_index;
						check_phi_index(paste_phi_index);
						if (geo.geometryProperty == pf::Geometry::Geo_Ellipsoid) {
							geo.ellipSolid.move(move_x, move_y, move_z);
							geo.ellipSolid.add_rotation_radian(rotate_1, rotate_2, rotate_3);
							Vector3 pv(geo.ellipSolid.core.x - Rotate_x, geo.ellipSolid.core.y - Rotate_y, geo.ellipSolid.core.z - Rotate_z);
							Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
							geo.ellipSolid.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
						}
						else if (geo.geometryProperty == pf::Geometry::Geo_Polyhedron) {
							geo.polyhedron.move(move_x, move_y, move_z);
							geo.polyhedron.add_rotation_radian(rotate_1, rotate_2, rotate_3);
							Vector3 pv(geo.polyhedron.point_inside_polyhedron.x - Rotate_x, geo.polyhedron.point_inside_polyhedron.y - Rotate_y, geo.polyhedron.point_inside_polyhedron.z - Rotate_z);
							Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
							geo.polyhedron.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
						}
						else if (geo.geometryProperty == pf::Geometry::Geo_SegmentedCylinder) {
							geo.cylinder.move(move_x, move_y, move_z);
							geo.cylinder.add_rotation_radian(rotate_1, rotate_2, rotate_3);
							Vector3 pv(geo.cylinder.geometric_center.x - Rotate_x, geo.cylinder.geometric_center.y - Rotate_y, geo.cylinder.geometric_center.z - Rotate_z);
							Vector3 pv2 = pv.get_rotated_vec3(RotationMatrix::rotationMatrix(Vector3(Rotate_1, Rotate_2, Rotate_3), pf::RotationGauge(homo_rotation_gauge)));
							geo.cylinder.move(pv2[0] - pv[0], pv2[1] - pv[1], pv2[2] - pv[2]);
						}
						microstructure_init::nucleation_box.geometry_box.push_back(geo);
					}
				}
			}
		}

		void definiteNucleation() {
			using namespace microstructure_init;
			bool _creatNewStorage = true;
			for (auto geo = nucleation_box.geometry_box.begin(); geo < nucleation_box.geometry_box.end();) {
				if (geo->generate_step == main_iterator::Current_ITE_step) {
					if (geo->isNormalized) {
						if (geo->phi > 1.0)
							geo->phi = 1.0;
						else if (geo->phi < 0.0)
							geo->phi = 0.0;
					}
					if (geo->geometryProperty == Geometry::Geo_Ellipsoid) {
						if (phi_parameters::is_phi_field_on) {
#pragma omp parallel for
							for (int z = 0; z < phi_parameters::phase_field.Nz(); z++)
								for (int y = 0; y < phi_parameters::phase_field.Ny(); y++)
									for (int x = 0; x < phi_parameters::phase_field.Nx(); x++) {
										PhaseFieldPoint& point = phi_parameters::phase_field(x, y, z);
										// -
										Vector3 p(x - geo->ellipSolid.core.x, y - geo->ellipSolid.core.y, z - geo->ellipSolid.core.z);
										p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->ellipSolid.radian_x, geo->ellipSolid.radian_y, geo->ellipSolid.radian_z),
											geo->ellipSolid.rotationGauge));
										p[0] += geo->ellipSolid.core.x;
										p[1] += geo->ellipSolid.core.y;
										p[2] += geo->ellipSolid.core.z;
										Point po(p[0], p[1], p[2]);
										bool check0 = geo->ellipSolid.check_point_inside_ellipsoid(po);
										if (!geo->isReverseRegion && check0) {
											if (geo->isNormalized) {
												REAL sum_phis = 0.0;
												for (int index = 0; index < phi_parameters::phi_number; index++)
													if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
														sum_phis += point.phi[index];
												if (sum_phis > SYS_EPSILON) {
													for (int index = 0; index < phi_parameters::phi_number; index++) {
														if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															point.phi[index] *= (1 - geo->phi) / sum_phis;
													}
												}
											}
											point.phi[geo->phaseIndex] = geo->phi;
										}
										else if (geo->isReverseRegion && !check0) {
											if (geo->isNormalized) {
												REAL sum_phis = 0.0;
												for (int index = 0; index < phi_parameters::phi_number; index++)
													if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
														sum_phis += point.phi[index];
												if (sum_phis > SYS_EPSILON) {
													for (int index = 0; index < phi_parameters::phi_number; index++) {
														if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															point.phi[index] *= (1 - geo->phi) / sum_phis;
													}
												}
											}
											point.phi[geo->phaseIndex] = geo->phi;
										}
									}
						}
						if (con_parameters::is_con_field_on) {
#pragma omp parallel for
							for (int z = 0; z < con_parameters::concentration_field.Nz(); z++)
								for (int y = 0; y < con_parameters::concentration_field.Ny(); y++)
									for (int x = 0; x < con_parameters::concentration_field.Nx(); x++) {
										ConcentrationFieldPoint& point = con_parameters::concentration_field(x, y, z);
										// -
										Vector3 p(x - geo->ellipSolid.core.x, y - geo->ellipSolid.core.y, z - geo->ellipSolid.core.z);
										p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->ellipSolid.radian_x, geo->ellipSolid.radian_y, geo->ellipSolid.radian_z),
											geo->ellipSolid.rotationGauge));
										p[0] += geo->ellipSolid.core.x;
										p[1] += geo->ellipSolid.core.y;
										p[2] += geo->ellipSolid.core.z;
										Point po(p[0], p[1], p[2]);
										bool check0 = geo->ellipSolid.check_point_inside_ellipsoid(po);
										if (!geo->isReverseRegion && check0) {
											point.con = geo->con;
										}
										else if (geo->isReverseRegion && !check0) {
											point.con = geo->con;
										}
									}
						}
						if (temp_parameters::is_temp_field_on) {
#pragma omp parallel for
							for (int z = 0; z < temp_parameters::temperature_field.Nz(); z++)
								for (int y = 0; y < temp_parameters::temperature_field.Ny(); y++)
									for (int x = 0; x < temp_parameters::temperature_field.Nx(); x++) {
										TemperatureFieldPoint& point = temp_parameters::temperature_field(x, y, z);
										// -
										Vector3 p(x - geo->ellipSolid.core.x, y - geo->ellipSolid.core.y, z - geo->ellipSolid.core.z);
										p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->ellipSolid.radian_x, geo->ellipSolid.radian_y, geo->ellipSolid.radian_z),
											geo->ellipSolid.rotationGauge));
										p[0] += geo->ellipSolid.core.x;
										p[1] += geo->ellipSolid.core.y;
										p[2] += geo->ellipSolid.core.z;
										Point po(p[0], p[1], p[2]);
										bool check0 = geo->ellipSolid.check_point_inside_ellipsoid(po);
										if (!geo->isReverseRegion && check0) {
											point.temp = geo->temperature;
										}
										else if (geo->isReverseRegion && !check0) {
											point.temp = geo->temperature;
										}
									}
						}
						stringstream report;
						report << "> A new Ellipsoid for grain : " << to_string(geo->phaseIndex) << " phase : " << materials_system::PHASES[phi_parameters::phi_property[geo->phaseIndex]] << " has been initialized at position ( "
							<< geo->ellipSolid.core.x << ", " << geo->ellipSolid.core.y << ", " << geo->ellipSolid.core.z << " )." << std::endl;
						WriteLog(report.str());
					}
					else if (geo->geometryProperty == Geometry::Geo_Polyhedron) {
						if (phi_parameters::is_phi_field_on) {
#pragma omp parallel for
							for (int z = 0; z < phi_parameters::phase_field.Nz(); z++)
								for (int y = 0; y < phi_parameters::phase_field.Ny(); y++)
									for (int x = 0; x < phi_parameters::phase_field.Nx(); x++) {
										PhaseFieldPoint& point = phi_parameters::phase_field(x, y, z);
										// -
										Vector3 pv(x - geo->polyhedron.point_inside_polyhedron.x, y - geo->polyhedron.point_inside_polyhedron.y, z - geo->polyhedron.point_inside_polyhedron.z);
										pv.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->polyhedron.radian_x, geo->polyhedron.radian_y, geo->polyhedron.radian_z),
											geo->polyhedron.rotationGauge));
										pv[0] += geo->polyhedron.point_inside_polyhedron.x;
										pv[1] += geo->polyhedron.point_inside_polyhedron.y;
										pv[2] += geo->polyhedron.point_inside_polyhedron.z;
										Point p(pv[0], pv[1], pv[2]);
										//Point p(x, y, z);
										bool check0 = geo->polyhedron.check_point_inside_polyhedron(p);
										if (!geo->isReverseRegion && check0) {
											if (geo->isNormalized) {
												REAL sum_phis = 0.0;
												for (int index = 0; index < phi_parameters::phi_number; index++)
													if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
														sum_phis += point.phi[index];
												if (sum_phis > SYS_EPSILON) {
													for (int index = 0; index < phi_parameters::phi_number; index++) {
														if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															point.phi[index] *= (1 - geo->phi) / sum_phis;
													}
												}
											}
											point.phi[geo->phaseIndex] = geo->phi;
										}
										else if (geo->isReverseRegion && !check0) {
											if (geo->isNormalized) {
												REAL sum_phis = 0.0;
												for (int index = 0; index < phi_parameters::phi_number; index++)
													if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
														sum_phis += point.phi[index];
												if (sum_phis > SYS_EPSILON) {
													for (int index = 0; index < phi_parameters::phi_number; index++) {
														if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															point.phi[index] *= (1 - geo->phi) / sum_phis;
													}
												}
											}
											point.phi[geo->phaseIndex] = geo->phi;
										}
									};
						}
						if (con_parameters::is_con_field_on) {
#pragma omp parallel for
							for (int z = 0; z < con_parameters::concentration_field.Nz(); z++)
								for (int y = 0; y < con_parameters::concentration_field.Ny(); y++)
									for (int x = 0; x < con_parameters::concentration_field.Nx(); x++) {
										ConcentrationFieldPoint& point = con_parameters::concentration_field(x, y, z);
										Vector3 pv(x - geo->polyhedron.point_inside_polyhedron.x, y - geo->polyhedron.point_inside_polyhedron.y, z - geo->polyhedron.point_inside_polyhedron.z);
										pv.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->polyhedron.radian_x, geo->polyhedron.radian_y, geo->polyhedron.radian_z),
											geo->polyhedron.rotationGauge));
										pv[0] += geo->polyhedron.point_inside_polyhedron.x;
										pv[1] += geo->polyhedron.point_inside_polyhedron.y;
										pv[2] += geo->polyhedron.point_inside_polyhedron.z;
										Point p(pv[0], pv[1], pv[2]);
										//Point p(x, y, z);
										bool check0 = geo->polyhedron.check_point_inside_polyhedron(p);
										if (!geo->isReverseRegion && check0) {
											point.con = geo->con;
										}
										else if (geo->isReverseRegion && !check0) {
											point.con = geo->con;
										}
									};
						}
						if (temp_parameters::is_temp_field_on) {
#pragma omp parallel for
							for (int z = 0; z < temp_parameters::temperature_field.Nz(); z++)
								for (int y = 0; y < temp_parameters::temperature_field.Ny(); y++)
									for (int x = 0; x < temp_parameters::temperature_field.Nx(); x++) {
										TemperatureFieldPoint& point = temp_parameters::temperature_field(x, y, z);
										Vector3 pv(x - geo->polyhedron.point_inside_polyhedron.x, y - geo->polyhedron.point_inside_polyhedron.y, z - geo->polyhedron.point_inside_polyhedron.z);
										pv.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->polyhedron.radian_x, geo->polyhedron.radian_y, geo->polyhedron.radian_z),
											geo->polyhedron.rotationGauge));
										pv[0] += geo->polyhedron.point_inside_polyhedron.x;
										pv[1] += geo->polyhedron.point_inside_polyhedron.y;
										pv[2] += geo->polyhedron.point_inside_polyhedron.z;
										Point p(pv[0], pv[1], pv[2]);
										//Point p(x, y, z);
										/*p.do_boundary(phaseMesh._bc_x_up, phaseMesh._bc_y_up, phaseMesh._bc_z_up, phaseMesh._bc_x_down, phaseMesh._bc_y_down, phaseMesh._bc_z_down,
											phaseMesh.limit_x, phaseMesh.limit_y, phaseMesh.limit_z);*/
										bool check0 = geo->polyhedron.check_point_inside_polyhedron(p);
										if (!geo->isReverseRegion && check0) {
											point.temp = geo->temperature;
										}
										else if (geo->isReverseRegion && !check0) {
											point.temp = geo->temperature;
										}
									};
						}
						stringstream report;
						report << "> A new Polygon for grain : " << to_string(geo->phaseIndex) << " phase: " << materials_system::PHASES[phi_parameters::phi_property[geo->phaseIndex]] << " has been initialized at position ( "
							<< geo->polyhedron.point_inside_polyhedron.x << ", " << geo->polyhedron.point_inside_polyhedron.y << ", " << geo->polyhedron.point_inside_polyhedron.z << " )." << std::endl;
						WriteLog(report.str());
					}
					else if (geo->geometryProperty == Geometry::Geo_SegmentedCylinder) {
						if (phi_parameters::is_phi_field_on) {
#pragma omp parallel for
							for (int z = 0; z < phi_parameters::phase_field.Nz(); z++)
								for (int y = 0; y < phi_parameters::phase_field.Ny(); y++)
									for (int x = 0; x < phi_parameters::phase_field.Nx(); x++) {
										PhaseFieldPoint& point = phi_parameters::phase_field(x, y, z);
										// Vector3 p(x, y, z);
										Vector3 p(x - geo->cylinder.geometric_center.x, y - geo->cylinder.geometric_center.y, z - geo->cylinder.geometric_center.z);
										p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->cylinder.radian_x, geo->cylinder.radian_y, geo->cylinder.radian_z),
											geo->cylinder.rotationGauge));
										p[0] += geo->cylinder.geometric_center.x;
										p[1] += geo->cylinder.geometric_center.y;
										p[2] += geo->cylinder.geometric_center.z;
										Point po(p[0], p[1], p[2]);
										bool check0 = geo->cylinder.check_point_inside_segmented_cylinder(po);
										if (!geo->isReverseRegion && check0) {
											if (geo->isNormalized) {
												REAL sum_phis = 0.0;
												for (int index = 0; index < phi_parameters::phi_number; index++)
													if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
														sum_phis += point.phi[index];
												if (sum_phis > SYS_EPSILON) {
													for (int index = 0; index < phi_parameters::phi_number; index++) {
														if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															point.phi[index] *= (1 - geo->phi) / sum_phis;
													}
												}
											}
											point.phi[geo->phaseIndex] = geo->phi;
										}
										else if (geo->isReverseRegion && !check0) {
											if (geo->isNormalized) {
												REAL sum_phis = 0.0;
												for (int index = 0; index < phi_parameters::phi_number; index++)
													if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
														sum_phis += point.phi[index];
												if (sum_phis > SYS_EPSILON) {
													for (int index = 0; index < phi_parameters::phi_number; index++) {
														if (point.phi[index] > SYS_EPSILON && index != geo->phaseIndex)
															point.phi[index] *= (1 - geo->phi) / sum_phis;
													}
												}
											}
											point.phi[geo->phaseIndex] = geo->phi;
										}
									};
						}
						if (con_parameters::is_con_field_on) {
#pragma omp parallel for
							for (int z = 0; z < con_parameters::concentration_field.Nz(); z++)
								for (int y = 0; y < con_parameters::concentration_field.Ny(); y++)
									for (int x = 0; x < con_parameters::concentration_field.Nx(); x++) {
										ConcentrationFieldPoint& point = con_parameters::concentration_field(x, y, z);
										// Vector3 p(x, y, z);
										Vector3 p(x - geo->cylinder.geometric_center.x, y - geo->cylinder.geometric_center.y, z - geo->cylinder.geometric_center.z);
										p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->cylinder.radian_x, geo->cylinder.radian_y, geo->cylinder.radian_z),
											geo->cylinder.rotationGauge));
										p[0] += geo->cylinder.geometric_center.x;
										p[1] += geo->cylinder.geometric_center.y;
										p[2] += geo->cylinder.geometric_center.z;
										Point po(p[0], p[1], p[2]);
										bool check0 = geo->cylinder.check_point_inside_segmented_cylinder(po);
										if (!geo->isReverseRegion && check0) {
											point.con = geo->con;
										}
										else if (geo->isReverseRegion && !check0) {
											point.con = geo->con;
										}
									};
						}
						if (temp_parameters::is_temp_field_on) {
#pragma omp parallel for
							for (int z = 0; z < temp_parameters::temperature_field.Nz(); z++)
								for (int y = 0; y < temp_parameters::temperature_field.Ny(); y++)
									for (int x = 0; x < temp_parameters::temperature_field.Nx(); x++) {
										TemperatureFieldPoint& point = temp_parameters::temperature_field(x, y, z);
										// Vector3 p(x, y, z);
										Vector3 p(x - geo->cylinder.geometric_center.x, y - geo->cylinder.geometric_center.y, z - geo->cylinder.geometric_center.z);
										p.do_rotate(RotationMatrix::rotationMatrix(Vector3(geo->cylinder.radian_x, geo->cylinder.radian_y, geo->cylinder.radian_z),
											geo->cylinder.rotationGauge));
										p[0] += geo->cylinder.geometric_center.x;
										p[1] += geo->cylinder.geometric_center.y;
										p[2] += geo->cylinder.geometric_center.z;
										Point po(p[0], p[1], p[2]);
										bool check0 = geo->cylinder.check_point_inside_segmented_cylinder(po);
										if (!geo->isReverseRegion && check0) {
											point.temp = geo->temperature;
										}
										else if (geo->isReverseRegion && !check0) {
											point.temp = geo->temperature;
										}
									}
						}
						stringstream report;
						report << "> A new Segmented Cylinder for grain : " << to_string(geo->phaseIndex) << " phase: " << materials_system::PHASES[phi_parameters::phi_property[geo->phaseIndex]] << " has been initialized at position ( "
							<< geo->cylinder.geometric_center.x << ", " << geo->cylinder.geometric_center.y << ", " << geo->cylinder.geometric_center.z << " )." << std::endl;
						WriteLog(report.str());
					}
					geo = nucleation_box.geometry_box.erase(geo);
				}
				else {
					geo++;
				}
			}
			for (auto point_set = nucleation_box.point_set_box.begin(); point_set < nucleation_box.point_set_box.end();) {
				if (point_set->generate_step == main_iterator::Current_ITE_step) {
					if (point_set->is_normalized) {
#pragma omp parallel for
						for (int index = 0; index < point_set->points_phi.size(); index++) {
							if (point_set->points_phi[index] > 1.0)
								point_set->points_phi[index] = 1.0;
							else if (point_set->points_phi[index] < 0.0)
								point_set->points_phi[index] = 0.0;
						}
					}
					if (phi_parameters::is_phi_field_on) {
#pragma omp parallel for
						for (int point_index = 0; point_index < point_set->points.size(); point_index++) {
							PhaseFieldPoint& point = phi_parameters::phase_field(REAL_to_int(point_set->points[point_index].x),
								REAL_to_int(point_set->points[point_index].y), REAL_to_int(point_set->points[point_index].z));
							if (point_set->is_normalized) {
								REAL sum_phis = 0.0;
								for (int index = 0; index < phi_parameters::phi_number; index++)
									if (point.phi[index] > SYS_EPSILON && index != point_set->phaseIndex)
										sum_phis += point.phi[index];
								if (sum_phis > SYS_EPSILON) {
									for (int index = 0; index < phi_parameters::phi_number; index++) {
										if (point.phi[index] > SYS_EPSILON && index != point_set->phaseIndex)
											point.phi[index] *= (1 - point_set->points_phi[point_index]) / sum_phis;
									}
								}
							}
							point.phi[point_set->phaseIndex] = point_set->points_phi[point_index];
						}
					}
					if (con_parameters::is_con_field_on) {
#pragma omp parallel for
						for (int point_index = 0; point_index < point_set->points.size(); point_index++) {
							ConcentrationFieldPoint& point = con_parameters::concentration_field(REAL_to_int(point_set->points[point_index].x),
								REAL_to_int(point_set->points[point_index].y), REAL_to_int(point_set->points[point_index].z));
							point.con = point_set->con;
						}
					}
					if (temp_parameters::is_temp_field_on) {
#pragma omp parallel for
						for (int point_index = 0; point_index < point_set->points.size(); point_index++) {
							TemperatureFieldPoint& point = temp_parameters::temperature_field(REAL_to_int(point_set->points[point_index].x),
								REAL_to_int(point_set->points[point_index].y), REAL_to_int(point_set->points[point_index].z));
							point.temp = point_set->temperature;
						}
					}
					stringstream report;
					report << "> A new PointSet for grain : " << to_string(point_set->phaseIndex) << " phase : "
						<< materials_system::PHASES[phi_parameters::phi_property[point_set->phaseIndex]] << " has been initialized;" << std::endl;
					WriteLog(report.str());
					point_set = nucleation_box.point_set_box.erase(point_set);
				}
				else {
					point_set++;
				}
			}
		}

	}
}