#pragma once
#include "../../base/sysfiles.h"
#include "../Module.h"
#include "../model_modules/Model_Params.h"
#include "../input_modules/ioFiles_Params.h"
#include "../model_modules/PhaseField/PhaseField_Params.h"
#include "../model_modules/ConcentrationField/ConcentrationField_Params.h"
#include "../model_modules/TemperatureField/TemperatureField_Params.h"
namespace pf {
	namespace write_mesh_data {
		const std::string mainName = "MeshData";
		const std::string format = ".dat";
		inline int output_frequence = -1;
		inline bool is_phi_mesh = true;
		inline bool is_con_mesh = true;
		inline bool is_temp_mesh = true;
		inline bool is_fluid_mesh = false;
		struct Data_MeshInfo {
			bool is_phi_mesh;
			int phi_number;
			int PNx;
			int PNy;
			int PNz;
			bool is_con_mesh;
			int con_number;
			int CNx;
			int CNy;
			int CNz;
			bool is_temp_mesh;
			int TNx;
			int TNy;
			int TNz;
			bool is_fluid_mesh;
			int FNx;
			int FNy;
			int FNz;
			void operator=(const Data_MeshInfo& n) {
				PNx = n.PNx;
				PNy = n.PNy;
				PNz = n.PNz;
				CNx = n.CNx;
				CNy = n.CNy;
				CNz = n.CNz;
				TNx = n.TNx;
				TNy = n.TNy;
				TNz = n.TNz;
				FNx = n.FNx;
				FNy = n.FNy;
				FNz = n.FNz;
				is_phi_mesh = n.is_phi_mesh;
				phi_number = n.phi_number;
				is_con_mesh = n.is_con_mesh;
				con_number = n.con_number;
				is_temp_mesh = n.is_temp_mesh;
				is_fluid_mesh = n.is_fluid_mesh;
			}
		};
		inline bool write_dataFile(std::string mark) {
			std::string fname;
			if (input_output_files_parameters::WorkingFolder_Path == "")
				fname = mainName + mark + format;
			else
				fname = input_output_files_parameters::WorkingFolder_Path + dirSeparator + mainName + mark + format;
			std::ofstream fout(fname, std::ios::binary);
			if (!fout) {
				std::cout << "Failed to write the data file!" << std::endl;
				fout.close();
				return false;
			}
			{ ///< defined in sequence(same with read)
				Data_MeshInfo mesh_info;
				mesh_info.is_phi_mesh = is_phi_mesh;
				mesh_info.phi_number = phi_parameters::phi_number;
				mesh_info.PNx = phi_parameters::phase_field.Nx();
				mesh_info.PNy = phi_parameters::phase_field.Ny();
				mesh_info.PNz = phi_parameters::phase_field.Nz();
				mesh_info.is_con_mesh = is_con_mesh;
				mesh_info.con_number = con_parameters::con_number;
				mesh_info.CNx = con_parameters::concentration_field.Nx();
				mesh_info.CNy = con_parameters::concentration_field.Ny();
				mesh_info.CNz = con_parameters::concentration_field.Nz();
				mesh_info.is_temp_mesh = is_temp_mesh;
				mesh_info.TNx = temp_parameters::temperature_field.Nx();
				mesh_info.TNy = temp_parameters::temperature_field.Ny();
				mesh_info.TNz = temp_parameters::temperature_field.Nz();
				mesh_info.is_fluid_mesh = is_fluid_mesh;
				mesh_info.FNx = 0;
				mesh_info.FNy = 0;
				mesh_info.FNz = 0;
				fout.write((const char*)&mesh_info, sizeof(Data_MeshInfo));
				///< storage for mesh
				if (mesh_info.is_phi_mesh) {
					for (int x = 0; x < mesh_info.PNx; x++)
						for (int y = 0; y < mesh_info.PNy; y++)
							for (int z = 0; z < mesh_info.PNz; z++) {
								PhaseFieldPoint& point = phi_parameters::phase_field(x, y, z);
								for (int index = 0; index < mesh_info.phi_number; index++)
									fout.write((const char*)&point.phi[index], sizeof(REAL));
							}
				}
				if (mesh_info.is_con_mesh) {
					for (int x = 0; x < mesh_info.CNx; x++)
						for (int y = 0; y < mesh_info.CNy; y++)
							for (int z = 0; z < mesh_info.CNz; z++) {
								ConcentrationFieldPoint& point = con_parameters::concentration_field(x, y, z);
								for (int index = 0; index < mesh_info.con_number; index++)
									fout.write((const char*)&point.con[index], sizeof(REAL));
							}
				}
				if (mesh_info.is_temp_mesh) {
					for (int x = 0; x < mesh_info.TNx; x++)
						for (int y = 0; y < mesh_info.TNy; y++)
							for (int z = 0; z < mesh_info.TNz; z++)
								fout.write((const char*)&temp_parameters::temperature_field(x, y, z).temp, sizeof(REAL));
				}
				if (mesh_info.is_fluid_mesh) {
					// - 
				}
			}
			fout.close();
			return true;
		}
		inline bool read_dataFile(string name, Data_MeshInfo& mesh_info) {
			std::fstream fin(name, std::ios::binary | ios::in);
			if (!fin) {
				std::cout << "Failed to read the aim file!" << std::endl;
				fin.close();
				return false;
			}
			fin.read((char*)&mesh_info, sizeof(Data_MeshInfo));
			///< read for mesh
			REAL data_buff = 0;
			if (mesh_info.is_phi_mesh) {
				int Nx = phi_parameters::phase_field.Nx(), Ny = phi_parameters::phase_field.Ny(),
					Nz = phi_parameters::phase_field.Nz(),
					phi_number = phi_parameters::phi_number;
				for (int x = 0; x < mesh_info.PNx; x++)
					for (int y = 0; y < mesh_info.PNy; y++)
						for (int z = 0; z < mesh_info.PNz; z++) {
							if (x < Nx && y < Ny && z < Nz) {
								PhaseFieldPoint& point = phi_parameters::phase_field(x, y, z);
								for (int index = 0; index < mesh_info.phi_number; index++) {
									fin.read((char*)&data_buff, sizeof(REAL));
									if (index < phi_number)
										point.phi[index] = data_buff;
								}
							}
							else {
								for (int index = 0; index < mesh_info.phi_number; index++)
									fin.read((char*)&data_buff, sizeof(REAL));
							}
						}
			}
			if (mesh_info.is_con_mesh) {
				int Nx = con_parameters::concentration_field.Nx(), Ny = con_parameters::concentration_field.Ny(),
					Nz = con_parameters::concentration_field.Nz(),
					con_number = con_parameters::con_number;
				for (int x = 0; x < mesh_info.CNx; x++)
					for (int y = 0; y < mesh_info.CNy; y++)
						for (int z = 0; z < mesh_info.CNz; z++) {
							if (x < Nx && y < Ny && z < Nz) {
								ConcentrationFieldPoint& point = con_parameters::concentration_field(x, y, z);
								for (int index = 0; index < mesh_info.con_number; index++) {
									fin.read((char*)&data_buff, sizeof(REAL));
									if (index < con_number)
										point.con[index] = data_buff;
								}
							}
							else {
								for (int index = 0; index < mesh_info.con_number; index++)
									fin.read((char*)&data_buff, sizeof(REAL));
							}
						}
			}
			if (mesh_info.is_temp_mesh) {
				int Nx = temp_parameters::temperature_field.Nx(), Ny = temp_parameters::temperature_field.Ny(),
					Nz = temp_parameters::temperature_field.Nz();
				for (int x = 0; x < mesh_info.TNx; x++)
					for (int y = 0; y < mesh_info.TNy; y++)
						for (int z = 0; z < mesh_info.TNz; z++) {
							if (x < Nx && y < Ny && z < Nz) {
								fin.read((char*)&temp_parameters::temperature_field(x, y, z).temp, sizeof(REAL));
							}
							else {
								fin.read((char*)&data_buff, sizeof(REAL));
							}
						}
			}
			if (mesh_info.is_fluid_mesh) {
				// - 
			}
			fin.close();
			return true;
		}
		inline void write_data_pre_iii() {
			write_dataFile("_init");
		}
		inline void write_data_pos_iii() {
			if (main_iterator::Current_ITE_step % output_frequence != 0)
				return;
			write_dataFile("_step" + to_string(main_iterator::Current_ITE_step));
		}
		inline void init() {
			InputFileReader::get_instance()->read_int_value("Solver.Output.MeshData.frequence", output_frequence, true);
			if (output_frequence == 0) {
				load_a_new_module(default_module_function, default_module_function, write_data_pre_iii,
					default_module_function, default_module_function, default_module_function,
					default_module_function, default_module_function, default_module_function, default_module_function);
			}
			else if (output_frequence > 0) {
				load_a_new_module(default_module_function, default_module_function, write_data_pre_iii,
					default_module_function, default_module_function, default_module_function,
					default_module_function, default_module_function, write_data_pos_iii, default_module_function);
			}
		}
	}
}