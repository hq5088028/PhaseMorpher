#pragma once
#include "PhaseFieldModels/PairwiseEquation.h"
#include "../../Module.h"
#include "../../postprocess_modules/WriteVTS.h"
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
				// - 
				pair_wise_functions::init_pairwise_equation();
				// - 
				load_a_new_module(default_module_function, default_module_function, pair_wise_functions::exec_pre_iii,
					pair_wise_functions::exec_i, default_module_function, default_module_function,
					default_module_function, default_module_function, default_module_function,
					pair_wise_functions::deinit);
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
}