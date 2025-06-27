#pragma once
#include "PhaseField_Params.h"
#include "../Model_Params.h"
#include "../../postprocess_modules/WriteVTS.h"
// - models
#include "PhaseFieldModels/PairwiseEquation.h"
namespace pf {
	inline void write_scalar_grains(ofstream& fout);
	inline void write_scalar_phi_index(ofstream& fout);
	inline void write_scalar_phi_property(ofstream& fout);
	inline void write_scalar_phi_name(ofstream& fout);
	inline void write_scalar_phi_all(ofstream& fout);
	inline void deinit_phase_field_equations() {
		for (int x = 0; x < phi_parameters::phase_field.Nx(); x++)
			for (int y = 0; y < phi_parameters::phase_field.Ny(); y++)
				for (int z = 0; z < phi_parameters::phase_field.Nz(); z++) {
					phi_parameters::phase_field(x, y, z).clear();
				}
		phi_parameters::phase_field.clear();
	}
	inline void init_pair_wise_equation_parameters();
	inline void init_phase_field_equations() {
		// - 
		int control_equation = 0;
		WriteDebugFile("# Solver.ControlEquation.PhaseField = 0 - Const, 1 - AllenCahn_Pairwise, ... \n");
		if (InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.PhaseField", control_equation, true)) {
			// - 
			phi_parameters::is_phi_field_on = true;
			// - 
			phi_parameters::phi_equation = PhaseFieldEquation(control_equation);
			// - 
			phi_parameters::phi_property_number = int(materials_system::PHASES.size());
			InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.PhaseField.PhaseNumber", phi_parameters::phi_number, true);
			phi_parameters::phi_property.resize(phi_parameters::phi_number);
			phi_parameters::phase_field.init(mesh_parameters::MESH_NX, mesh_parameters::MESH_NY, mesh_parameters::MESH_NZ);
			for (int x = 0; x < phi_parameters::phase_field.Nx(); x++)
				for (int y = 0; y < phi_parameters::phase_field.Ny(); y++)
					for (int z = 0; z < phi_parameters::phase_field.Nz(); z++) {
						phi_parameters::phase_field(x, y, z).init(phi_parameters::phi_number);
					}
			// - 
			load_a_new_module(default_module_function, default_module_function, default_module_function,
				default_module_function, default_module_function, default_module_function,
				default_module_function, default_module_function, default_module_function,
				deinit_phase_field_equations);
			// - 
			for (int index = 0; index < phi_parameters::phi_number; index++)
				phi_parameters::phi_property[index] = -1;
			if (phi_parameters::phi_property_number != 0)
				WriteDebugFile("# Solver.ControlEquation.PhaseField.Property.PhiName = ( Phi_index_1, Phi_index_2, ... ) \n");
			string phi_property_key = "Solver.ControlEquation.PhaseField.Property.", phi_input = "()";
			vector<int> not_define;
			for (int index = 0; index < phi_parameters::phi_property_number; index++) {
				InputFileReader::get_instance()->read_string_value(phi_property_key + materials_system::PHASES[index], phi_input, true);
				vector<input_value> phi_property_value = InputFileReader::get_instance()->trans_matrix_1d_const_to_input_value(InputValueType::IVType_INT, phi_property_key + materials_system::PHASES[index], phi_input, true);
				for (int index2 = 0; index2 < phi_property_value.size(); index2++)
					if (phi_parameters::phi_property[phi_property_value[index2].int_value] == -1)
						phi_parameters::phi_property[phi_property_value[index2].int_value] = index;
					else {
						WriteError("> ERROR, phi index = " + to_string(phi_property_value[index2].int_value) + " has been defined multiple times !\n");
						SYS_PROGRAM_STOP;
					}
			}
			for (int index = 0; index < phi_parameters::phi_number; index++)
				if (phi_parameters::phi_property[index] == -1)
					not_define.push_back(index);
			if (not_define.size() != 0) {
				string error = "> ERROR, phi index = ";
				for (int index = 0; index < not_define.size(); index++)
					error += to_string(not_define[index]) + " , ";
				error += " has not been defined !\n";
				WriteError(error);
				SYS_PROGRAM_STOP;
			}
			switch (phi_parameters::phi_equation)
			{
			case PhaseFieldEquation::PFE_AllenCahn_Pairwise: {
				init_pair_wise_equation_parameters();
				break;
			}
			default:
				break;
			}
			bool output_vts = false;
			InputFileReader::get_instance()->read_bool_value("Solver.ControlEquation.PhaseField.VTS.grains", output_vts, true);
			if (output_vts)
				load_vts_scalar_func(write_scalar_grains);
			output_vts = false;
			InputFileReader::get_instance()->read_bool_value("Solver.ControlEquation.PhaseField.VTS.index", output_vts, true);
			if (output_vts)
				load_vts_scalar_func(write_scalar_phi_index);
			output_vts = false;
			InputFileReader::get_instance()->read_bool_value("Solver.ControlEquation.PhaseField.VTS.property", output_vts, true);
			if (output_vts)
				load_vts_scalar_func(write_scalar_phi_property);
			output_vts = false;
			InputFileReader::get_instance()->read_bool_value("Solver.ControlEquation.PhaseField.VTS.name", output_vts, true);
			if (output_vts)
				load_vts_scalar_func(write_scalar_phi_name);
			output_vts = false;
			InputFileReader::get_instance()->read_bool_value("Solver.ControlEquation.PhaseField.VTS.all_phi", output_vts, true);
			if (output_vts)
				load_vts_scalar_func(write_scalar_phi_all);
		}
	}
	inline void write_scalar_grains(ofstream& fout) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "grains" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = write_vts::z_begin; k <= write_vts::z_end; ++k)
			for (int j = write_vts::y_begin; j <= write_vts::y_end; ++j)
				for (int i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
					PhaseFieldPoint& point = phi_parameters::phase_field(i, j, k);
					REAL fix = 0.0;
					for (int index = 0; index < phi_parameters::phi_number; index++)
						fix += point.phi[index] * point.phi[index];
					if (std::isnan(fix))
						fout << NaN() << endl;
					else
						fout << 1.0 - fix << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void write_scalar_phi_index(ofstream& fout) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "phi_index" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = write_vts::z_begin; k <= write_vts::z_end; ++k)
			for (int j = write_vts::y_begin; j <= write_vts::y_end; ++j)
				for (int i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
					PhaseFieldPoint& point = phi_parameters::phase_field(i, j, k);
					REAL fix = 0.0;
					for (int index = 0; index < phi_parameters::phi_number; index++)
						fix += point.phi[index] * index;
					fout << fix << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void write_scalar_phi_property(ofstream& fout) {
		fout << "<DataArray type = \"Float64\" Name = \"" << "phi_property" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = write_vts::z_begin; k <= write_vts::z_end; ++k)
			for (int j = write_vts::y_begin; j <= write_vts::y_end; ++j)
				for (int i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
					PhaseFieldPoint& point = phi_parameters::phase_field(i, j, k);
					REAL fix = 0.0;
					for (int index = 0; index < phi_parameters::phi_number; index++)
						fix += point.phi[index] * phi_parameters::phi_property[index];
					fout << fix << endl;
				}
		fout << "</DataArray>" << endl;
	}
	inline void write_scalar_phi_name(ofstream& fout) {
		for (int pindex = 0; pindex < materials_system::PHASES.size(); pindex++) {
			string phi_name = "phi_" + materials_system::PHASES[pindex];
			fout << "<DataArray type = \"Float64\" Name = \"" << phi_name <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (int j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (int i = write_vts::x_begin; i <= write_vts::x_end; ++i) {
						PhaseFieldPoint& point = phi_parameters::phase_field(i, j, k);
						REAL fix = 0.0;
						for (int index = 0; index < phi_parameters::phi_number; index++)
							if (pindex == phi_parameters::phi_property[index])
								fix += point.phi[index];
						fout << fix << endl;
					}
			fout << "</DataArray>" << endl;
		}
	}
	inline void write_scalar_phi_all(ofstream& fout) {
		for (int pindex = 0; pindex < phi_parameters::phi_number; pindex++) {
			string phi_name = "phi_" + to_string(pindex);
			fout << "<DataArray type = \"Float64\" Name = \"" << phi_name <<
				"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = write_vts::z_begin; k <= write_vts::z_end; ++k)
				for (int j = write_vts::y_begin; j <= write_vts::y_end; ++j)
					for (int i = write_vts::x_begin; i <= write_vts::x_end; ++i)
						fout << phi_parameters::phase_field(i, j, k).phi[pindex] << endl;
			fout << "</DataArray>" << endl;
		}
	}
	inline void init_pair_wise_equation_parameters() {
		// - 
		pair_wise_parameters::MESH_NX = mesh_parameters::MESH_NX;
		pair_wise_parameters::MESH_NY = mesh_parameters::MESH_NY;
		pair_wise_parameters::MESH_NZ = mesh_parameters::MESH_NZ;
		pair_wise_parameters::x_down = mesh_parameters::x_down;
		pair_wise_parameters::y_down = mesh_parameters::y_down;
		pair_wise_parameters::z_down = mesh_parameters::z_down;
		pair_wise_parameters::x_up = mesh_parameters::x_up;
		pair_wise_parameters::y_up = mesh_parameters::y_up;
		pair_wise_parameters::z_up = mesh_parameters::z_up;
		// - 
		pair_wise_parameters::phase_field = &phi_parameters::phase_field;
		pair_wise_parameters::PHI_MAX_VARIATION = &phi_parameters::PHI_MAX_VARIATION;
		pair_wise_parameters::delt_t = &time_parameters::delt_t;
		pair_wise_parameters::delt_r = &mesh_parameters::delt_r;
		pair_wise_parameters::phi_name_property = materials_system::phi_property;   // func ptr
		pair_wise_parameters::phi_number = phi_parameters::phi_number;
		pair_wise_parameters::phi_property_number = phi_parameters::phi_property_number;
		pair_wise_parameters::phi_property = phi_parameters::phi_property;
		// - 
		InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.Acceleration.ContainerSize", pair_wise_parameters::phi_acc_number, true);
		pair_wise_parameters::phase_field_pairwise.init(pair_wise_parameters::MESH_NX, pair_wise_parameters::MESH_NY, pair_wise_parameters::MESH_NZ);
		pair_wise_parameters::phase_field_Gc.init(pair_wise_parameters::MESH_NX, pair_wise_parameters::MESH_NY, pair_wise_parameters::MESH_NZ);
		for (int x = 0; x < pair_wise_parameters::phase_field_pairwise.Nx(); x++)
			for (int y = 0; y < pair_wise_parameters::phase_field_pairwise.Ny(); y++)
				for (int z = 0; z < pair_wise_parameters::phase_field_pairwise.Nz(); z++) {
					pair_wise_parameters::phase_field_pairwise(x, y, z).init(pair_wise_parameters::phi_acc_number, pair_wise_parameters::phi_number);
				}
		InputFileReader::get_instance()->read_bool_value("Solver.ControlEquation.Pairwise.NormalizePhi", pair_wise_parameters::is_phi_normalized, true);

		int difference_method = pair_wise_parameters::DifferenceMethod::FIVE_POINT;
		WriteDebugFile("# Solver.ControlEquation.Pairwise.DifferenceMethod : 0 - FIVE_POINT , 1 - NINE_POINT \n");
		InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.DifferenceMethod", difference_method, true);
		pair_wise_parameters::diff_method = pair_wise_parameters::DifferenceMethod(difference_method);

		WriteDebugFile("# .Pairwise.GrainsOrientations.rotation_gauge = 0 - XYX, 1 - XZX, 2 - YXY, 3 - YZY, 4  - ZXZ, 5  - ZYZ \n");
		WriteDebugFile("#                                               6 - XYZ, 7 - XZY, 8 - YXZ, 9 - YZX, 10 - ZXY, 11 - ZYX \n");
		int rotation_gauge = RotationGauge::RG_ZXZ;
		InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.GrainsOrientations.rotation_gauge", rotation_gauge, true);
		pair_wise_parameters::grains_orientation.init(RotationGauge(rotation_gauge), pair_wise_parameters::phi_number);

		WriteDebugFile("# Solver.ControlEquation.Pairwise.GrainsOrientations = {[(phi_index_0, phi_index_2, ... ),(rotation_angle_1, rotation_angle_2, rotation_angle_3)],  ... } \n");
		string grains_key = "Solver.ControlEquation.Pairwise.GrainsOrientations", grains_input = "{}";
		InputFileReader::get_instance()->read_string_value(grains_key, grains_input, true);
		vector<InputValueType> grains_structure; grains_structure.push_back(InputValueType::IVType_INT); grains_structure.push_back(InputValueType::IVType_REAL);
		vector<vector<vector<input_value>>> grains_value;
		grains_value = InputFileReader::get_instance()->trans_matrix_3d_const_array_const_to_input_value(grains_structure, grains_key, grains_input, true);
		for (int i = 0; i < grains_value.size(); i++) {
			for (int j = 0; j < grains_value[i][0].size(); j++)
				pair_wise_parameters::grains_orientation.set_phi_orientation(grains_value[i][0][j].int_value,
					Vector3(AngleToRadians(grains_value[i][1][0].REAL_value), AngleToRadians(grains_value[i][1][1].REAL_value), AngleToRadians(grains_value[i][1][2].REAL_value)));
		}

		// - interface mobility
		pair_wise_functions::mobility = pair_wise_functions::Lij_pair_wise_const;
		WriteDebugFile("# Solver.ControlEquation.Pairwise.Lij.const  = Lij_value \n");
		WriteDebugFile("#                                .matrix = [(phi_i, phi_j, Lij_value), ... ] \n");
		WriteDebugFile("#                                .block = [(phi_begin, phi_end, Lij_value), ... ] \n");
		int intMob_type = 0;
		string matrix_key = "Solver.ControlEquation.Pairwise.Lij.matrix", matrix_input = "[()]", block_key = "Solver.ControlEquation.Pairwise.Lij.block", block_input = "[()]";
		if (InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.Lij.const", pair_wise_parameters::const_Lij, true)) {
			pair_wise_functions::mobility = pair_wise_functions::Lij_pair_wise_const;
			intMob_type = 0;
			WriteDebugFile("# Solver.ControlEquation.Pairwise.Lij.Const.block       = [(phi_index1, phi_index2, ... ), ... ] \n");
			WriteDebugFile("#                                    .cross_block = [(phi_alpha1, phi_alpha2, ... ), (phi_beta1, phi_beta2, ... )] \n");
			string const_block_key = "Solver.ControlEquation.Pairwise.Lij.Const.block", const_cross_block_key = "Solver.ControlEquation.Pairwise.Lij.Const.cross_block", const_block_input = "[()]";
			if (InputFileReader::get_instance()->read_string_value(const_block_key, const_block_input, true)) {
				pair_wise_functions::mobility = pair_wise_functions::Lij_pair_wise_const_block;
				intMob_type = 1;
				vector<vector<input_value>> const_block_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
				for (int index = 0; index < const_block_value.size(); index++) {
					vector<int> const_block;
					for (int phi_index = 0; phi_index < const_block_value[index].size(); phi_index++)
						const_block.push_back(const_block_value[index][phi_index].int_value);
					pair_wise_parameters::const_block_Lij.push_back(const_block);
				}
			}
			else if (InputFileReader::get_instance()->read_string_value(const_cross_block_key, const_block_input, true)) {
				pair_wise_functions::mobility = pair_wise_functions::Lij_pair_wise_const_cross_block;
				intMob_type = 4;
				vector<vector<input_value>> const_block_value = InputFileReader::get_instance()->trans_matrix_2d_const_const_to_input_value(InputValueType::IVType_INT, const_block_key, const_block_input, true);
				for (int index = 0; index < 2; index++) {
					vector<int> const_block;
					for (int phi_index = 0; phi_index < const_block_value[index].size(); phi_index++)
						const_block.push_back(const_block_value[index][phi_index].int_value);
					pair_wise_parameters::const_block_Lij.push_back(const_block);
				}
			}
		}
		else if (InputFileReader::get_instance()->read_string_value(matrix_key, matrix_input, true)) {
			pair_wise_functions::mobility = pair_wise_functions::Lij_matrix;
			pair_wise_parameters::matirx_Lij.resize(pair_wise_parameters::phi_number);
			for (int index = 0; index < pair_wise_parameters::phi_number; index++)
				pair_wise_parameters::matirx_Lij[index].resize(pair_wise_parameters::phi_number);
			intMob_type = 2;
			vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_INT); matrix_structure.push_back(InputValueType::IVType_REAL);
			vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key, matrix_input, true);
			for (int index = 0; index < matrix_value.size(); index++)
				pair_wise_parameters::matirx_Lij[matrix_value[index][0].int_value][matrix_value[index][1].int_value] = matrix_value[index][2].REAL_value;
		}
		else if (InputFileReader::get_instance()->read_string_value(block_key, block_input, true)) {
			pair_wise_functions::mobility = pair_wise_functions::Lij_block;
			intMob_type = 3;
			vector<InputValueType> block_structure; block_structure.push_back(InputValueType::IVType_INT);
			block_structure.push_back(InputValueType::IVType_INT); block_structure.push_back(InputValueType::IVType_REAL);
			vector<vector<input_value>> block_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(block_structure, block_key, block_input, true);
			for (int index = 0; index < block_value.size(); index++) {
				pair_wise_parameters::block_Lij.push_back({ block_value[index][0].int_value, block_value[index][1].int_value });
				pair_wise_parameters::block_value_Lij.push_back(block_value[index][2].REAL_value);
			}
		}
		// - anisotropic
		WriteDebugFile("# Solver.ControlEquation.Pairwise.Lij.Anisotropic.matrix = [(phi_a_name, phi_b_name, AnisotropicModelIndex ), ...]\n");
		WriteDebugFile("#                                     AnisotropicModelIndex : 0 - IEA_ISO, 1 - IEA_CUBIC, 2 - IEA_HEX_BOETTGER, 3 - IEA_HEX_SUN, 4 - IEA_HEX_YANG \n");
		string matrix_aniso_key = "Solver.ControlEquation.Pairwise.Anisotropic.matrix", matrix_aniso_input = "[()]";
		if (InputFileReader::get_instance()->read_string_value(matrix_aniso_key, matrix_aniso_input, true)) {
			if (intMob_type == 0)
				pair_wise_functions::mobility = pair_wise_functions::Lij_pair_wise_anisotropic_const;
			else if (intMob_type == 1)
				pair_wise_functions::mobility = pair_wise_functions::Lij_pair_wise_anisotropic_const_block;
			else if (intMob_type == 2)
				pair_wise_functions::mobility = pair_wise_functions::Lij_anisotropic_matrix;
			else if (intMob_type == 3)
				pair_wise_functions::mobility = pair_wise_functions::Lij_anisotropic_block;
			else if (intMob_type == 4)
				pair_wise_functions::mobility = pair_wise_functions::Lij_pair_wise_anisotropic_const_cross_block;
			pair_wise_parameters::intMob_anisotropic_model_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intMobAniso1_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intMobAniso2_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intMobAniso3_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intMobAniso4_matrix.resize(pair_wise_parameters::phi_property_number);
			for (int index = 0; index < pair_wise_parameters::phi_property_number; index++) {
				pair_wise_parameters::intMob_anisotropic_model_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intMobAniso1_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intMobAniso2_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intMobAniso3_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intMobAniso4_matrix[index].resize(pair_wise_parameters::phi_property_number);
			}
			vector<InputValueType> matrix_aniso_structure; matrix_aniso_structure.push_back(InputValueType::IVType_STRING); matrix_aniso_structure.push_back(InputValueType::IVType_STRING); matrix_aniso_structure.push_back(InputValueType::IVType_INT);
			vector<vector<input_value>> matrix_aniso_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_aniso_structure, matrix_aniso_key, matrix_aniso_input, true);
			for (int index = 0; index < matrix_aniso_value.size(); index++) {
				int alpha_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][0].string_value),
					beta_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][1].string_value);
				pair_wise_parameters::intMob_anisotropic_model_matrix[alpha_property][beta_property] = matrix_aniso_value[index][2].int_value;
				pair_wise_parameters::intMob_anisotropic_model_matrix[beta_property][alpha_property] = matrix_aniso_value[index][2].int_value;
				if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_CUBIC) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
				}
				else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_HEX_BOETTGER) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
				}
				else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_HEX_SUN) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
					string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
					string matrix_aniso_par4 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_4";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
					pair_wise_parameters::intMobAniso2_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso2_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
					pair_wise_parameters::intMobAniso3_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso3_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par4, aniso_par, true);
					pair_wise_parameters::intMobAniso4_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso4_matrix[beta_property][alpha_property] = aniso_par;
				}
				else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Mobility_Anisotropic::IMA_HEX_YANG) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
					string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.Lij.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intMobAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso1_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
					pair_wise_parameters::intMobAniso2_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso2_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
					pair_wise_parameters::intMobAniso3_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intMobAniso3_matrix[beta_property][alpha_property] = aniso_par;
				}
			}
		}

		// - interface energy
		pair_wise_functions::_xi_ab = pair_wise_functions::xi_ab_const;
		pair_wise_functions::_xi_abc = pair_wise_functions::xi_abc_const;
		WriteDebugFile("# Solver.ControlEquation.Pairwise.InterfaceEnergy.int_gradient : 0 - Steinbach_1996 , 1 - Steinbach_1999 , 2 - Steinbach_G2009 , 3 - DanielCrack_G2016\n");
		InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.InterfaceEnergy.int_gradient", pair_wise_parameters::interface_gradient, true);
		WriteDebugFile("# Solver.ControlEquation.Pairwise.InterfaceEnergy.int_potential : 0 - Nestler_Well , 1 - Nestler_Obstacle , 2 - Steinbach_P2009 , 3 - DanielCrack_P2016\n");
		InputFileReader::get_instance()->read_int_value("Solver.ControlEquation.Pairwise.InterfaceEnergy.int_potential", pair_wise_parameters::interface_potential, true);
		InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.interface_width", pair_wise_parameters::interface_width, true);

		WriteDebugFile("# Solver.ControlEquation.Pairwise.xi_ab.const  = xi_ab \n");
		WriteDebugFile("#                                      .matrix = [(phi_a_name, phi_b_name, xi_ab_value), ...] \n");

		int int_energy_type = 0;
		string matrix_key1 = "Solver.ControlEquation.Pairwise.xi_ab.matrix", matrix_input1 = "[()]";
		if (InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.xi_ab.const", pair_wise_parameters::const_xi_ab, true)) {
			pair_wise_functions::_xi_ab = pair_wise_functions::xi_ab_const;
			int_energy_type = 0;
		}
		else if (InputFileReader::get_instance()->read_string_value(matrix_key1, matrix_input1, true)) {
			pair_wise_functions::_xi_ab = pair_wise_functions::xi_ab_matrix;
			pair_wise_parameters::matirx_xi_ab.resize(pair_wise_parameters::phi_property_number);
			for (int index = 0; index < pair_wise_parameters::phi_property_number; index++)
				pair_wise_parameters::matirx_xi_ab[index].resize(pair_wise_parameters::phi_property_number);
			int_energy_type = 1;
			vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_REAL);
			vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key1, matrix_input1, true);
			for (int index = 0; index < matrix_value.size(); index++) {
				int alpha_property = pair_wise_parameters::phi_name_property(matrix_value[index][0].string_value),
					beta_property = pair_wise_parameters::phi_name_property(matrix_value[index][1].string_value);
				pair_wise_parameters::matirx_xi_ab[alpha_property][beta_property] = matrix_value[index][2].REAL_value;
			}
		}

		WriteDebugFile("# Solver.ControlEquation.Pairwise.xi_ab.Anisotropic.matrix = [(phi_a_name, phi_b_name, AnisotropicModelIndex ), ...]\n");
		WriteDebugFile("#                                       AnisotropicModelIndex : 0 - IEA_ISO, 1 - IEA_CUBIC, 2 - IEA_HEX_BOETTGER, 3 - IEA_HEX_SUN, 4 - IEA_HEX_YANG \n");
		string matrix_aniso_key2 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic.matrix", matrix_aniso_input2 = "[()]";
		if (InputFileReader::get_instance()->read_string_value(matrix_aniso_key2, matrix_aniso_input2, true)) {
			if (int_energy_type == 0)
				pair_wise_functions::_xi_ab = pair_wise_functions::xi_ab_anisotropic_const;
			else if (int_energy_type == 1)
				pair_wise_functions::_xi_ab = pair_wise_functions::xi_ab_anisotropic_matrix;
			pair_wise_parameters::intEn_anisotropic_model_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intEnAniso1_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intEnAniso2_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intEnAniso3_matrix.resize(pair_wise_parameters::phi_property_number);
			pair_wise_parameters::intEnAniso4_matrix.resize(pair_wise_parameters::phi_property_number);
			for (int index = 0; index < pair_wise_parameters::phi_property_number; index++) {
				pair_wise_parameters::intEn_anisotropic_model_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intEnAniso1_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intEnAniso2_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intEnAniso3_matrix[index].resize(pair_wise_parameters::phi_property_number);
				pair_wise_parameters::intEnAniso4_matrix[index].resize(pair_wise_parameters::phi_property_number);
			}
			vector<InputValueType> matrix_aniso_structure; matrix_aniso_structure.push_back(InputValueType::IVType_STRING);
			matrix_aniso_structure.push_back(InputValueType::IVType_STRING); matrix_aniso_structure.push_back(InputValueType::IVType_INT);
			vector<vector<input_value>> matrix_aniso_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_aniso_structure, matrix_aniso_key, matrix_aniso_input, true);
			for (int index = 0; index < matrix_aniso_value.size(); index++) {
				int alpha_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][0].string_value),
					beta_property = pair_wise_parameters::phi_name_property(matrix_aniso_value[index][1].string_value);
				pair_wise_parameters::intEn_anisotropic_model_matrix[alpha_property][beta_property] = matrix_aniso_value[index][2].int_value;
				pair_wise_parameters::intEn_anisotropic_model_matrix[beta_property][alpha_property] = matrix_aniso_value[index][2].int_value;
				if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_CUBIC) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
				}
				else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_HEX_BOETTGER) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
				}
				else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_HEX_SUN) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
					string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
					string matrix_aniso_par4 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_4";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
					pair_wise_parameters::intEnAniso2_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso2_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
					pair_wise_parameters::intEnAniso3_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso3_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par4, aniso_par, true);
					pair_wise_parameters::intEnAniso4_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso4_matrix[beta_property][alpha_property] = aniso_par;
				}
				else if (matrix_aniso_value[index][2].int_value == pair_wise_parameters::Int_Energy_Anisotropic::IEA_HEX_YANG) {
					string matrix_aniso_par1 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_1";
					string matrix_aniso_par2 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_2";
					string matrix_aniso_par3 = "Solver.ControlEquation.Pairwise.xi_ab.Anisotropic." + matrix_aniso_value[index][0].string_value + "|" + matrix_aniso_value[index][1].string_value + ".parameter_3";
					REAL aniso_par = 0.0;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par1, aniso_par, true);
					pair_wise_parameters::intEnAniso1_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso1_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par2, aniso_par, true);
					pair_wise_parameters::intEnAniso2_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso2_matrix[beta_property][alpha_property] = aniso_par;
					InputFileReader::get_instance()->read_REAL_value(matrix_aniso_par3, aniso_par, true);
					pair_wise_parameters::intEnAniso3_matrix[alpha_property][beta_property] = aniso_par;
					pair_wise_parameters::intEnAniso3_matrix[beta_property][alpha_property] = aniso_par;
				}
			}
		}

		WriteDebugFile("# Solver.ControlEquation.Pairwise.xi_abc.const  = xi_ab \n");
		WriteDebugFile("#                                       .matrix = [(phi_a_name, phi_b_name, phi_c_name, xi_abc_value), ...] \n");

		string matrix_key2 = "Solver.ControlEquation.Pairwise.xi_abc.matrix", matrix_input2 = "[()]";
		if (InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.xi_abc.const", pair_wise_parameters::const_xi_abc, true)) {
			pair_wise_functions::_xi_abc = pair_wise_functions::xi_abc_const;
		}
		else if (InputFileReader::get_instance()->read_string_value(matrix_key2, matrix_input2, true)) {
			pair_wise_functions::_xi_abc = pair_wise_functions::xi_abc_matrix;
			pair_wise_parameters::matirx_xi_abc.resize(pair_wise_parameters::phi_property_number);
			for (int index1 = 0; index1 < pair_wise_parameters::phi_property_number; index1++) {
				pair_wise_parameters::matirx_xi_abc[index1].resize(pair_wise_parameters::phi_property_number);
				for (int index2 = 0; index2 < pair_wise_parameters::phi_property_number; index2++)
					pair_wise_parameters::matirx_xi_abc[index1][index2].resize(pair_wise_parameters::phi_property_number);
			}
			vector<InputValueType> matrix_structure; matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_STRING)
				; matrix_structure.push_back(InputValueType::IVType_STRING); matrix_structure.push_back(InputValueType::IVType_REAL);
			vector<vector<input_value>> matrix_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value(matrix_structure, matrix_key2, matrix_input2, true);
			for (int index = 0; index < matrix_value.size(); index++) {
				int alpha_property = pair_wise_parameters::phi_name_property(matrix_value[index][0].string_value),
					beta_property = pair_wise_parameters::phi_name_property(matrix_value[index][1].string_value),
					gamma_property = pair_wise_parameters::phi_name_property(matrix_value[index][2].string_value);
				pair_wise_parameters::matirx_xi_abc[alpha_property][beta_property][gamma_property] = matrix_value[index][3].REAL_value;
				pair_wise_parameters::matirx_xi_abc[alpha_property][gamma_property][beta_property] = matrix_value[index][3].REAL_value;
				pair_wise_parameters::matirx_xi_abc[beta_property][alpha_property][gamma_property] = matrix_value[index][3].REAL_value;
				pair_wise_parameters::matirx_xi_abc[beta_property][gamma_property][alpha_property] = matrix_value[index][3].REAL_value;
				pair_wise_parameters::matirx_xi_abc[gamma_property][alpha_property][beta_property] = matrix_value[index][3].REAL_value;
				pair_wise_parameters::matirx_xi_abc[gamma_property][beta_property][alpha_property] = matrix_value[index][3].REAL_value;
			}
		}

		// - multiple crack
		if (pair_wise_parameters::interface_gradient == pair_wise_parameters::Int_Gradient::DanielCrack_G2016 || pair_wise_parameters::interface_potential == pair_wise_parameters::Int_Potential::DanielCrack_P2016) {

			pair_wise_parameters::is_phi_normalized = true;
			pair_wise_parameters::phase_Gc.resize(pair_wise_parameters::phi_property_number);

			string crack_name = "";
			if (InputFileReader::get_instance()->read_string_value("Solver.ControlEquation.Pairwise.crack_phi_name", crack_name, true))
				pair_wise_parameters::crack_phase_property = pair_wise_parameters::phi_name_property(crack_name);

			InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.crack_int_width", pair_wise_parameters::crack_int_width, true);

			InputFileReader::get_instance()->read_REAL_value("Solver.ControlEquation.Pairwise.crack_phi_cut_off", pair_wise_parameters::crack_num_cut_off, true);

			WriteDebugFile("# Solver.ControlEquation.Pairwise.Gc = [(phi_name, Gc_value), ...] \n");
			string gc_key = "Solver.ControlEquation.Pairwise.Gc", gc_input = "[()]";
			if (InputFileReader::get_instance()->read_string_value(gc_key, gc_input, true)) {
				vector<vector<input_value>> gc_value = InputFileReader::get_instance()->trans_matrix_2d_const_array_to_input_value({ InputValueType::IVType_STRING,InputValueType::IVType_REAL }, gc_key, gc_input, true);
				for (int index = 0; index < gc_value.size(); index++)
					pair_wise_parameters::phase_Gc[pair_wise_parameters::phi_name_property(gc_value[index][0].string_value)] = gc_value[index][1].REAL_value;
			}
		}
		// - load this module
		load_a_new_module(default_module_function, default_module_function, pair_wise_functions::exec_pre_iii,
			pair_wise_functions::exec_i, default_module_function, default_module_function,
			default_module_function, default_module_function, default_module_function,
			pair_wise_functions::deinit);
		// - 
	}
}