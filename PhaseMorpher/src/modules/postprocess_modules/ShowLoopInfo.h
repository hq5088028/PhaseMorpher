#pragma once
#include "../input_modules/inputfiles/InputFileReader.h"
#include "../Module.h"
#include "../MainIterator_Params.h"
#include "../model_modules/Model_Params.h"
#include "../model_modules/PhaseField/PhaseField_Params.h"
#include "../model_modules/ConcentrationField/ConcentrationField_Params.h"
#include "../model_modules/TemperatureField/TemperatureField_Params.h"
#include "../../tools/timer.h"
#include "../postprocess_modules/DataStatistics.h"
namespace pf {
	namespace show_loop_information {
		inline int screen_loop_step   = INT_MAX;
		inline int screen_output_step = INT_MAX;
		inline void exec_pre_iii() {
			pf::PhaseFieldPoint phi_info = data_statistics_functions::statistical_phi();
			stringstream log;
			log << "> PHI CON TEMP information :" << endl;
			for (int index = 0; index < phi_parameters::phi_number; index++)
				log << ">  Phase_" << index << "(" << materials_system::PHASES[phi_parameters::phi_property[index]] << ")" << " : "
				<< setprecision(5) << phi_info.phi[index] * 100 << " %" << endl;
			WriteLog(log.str());
		}

		inline void exec_pos_i() {
			if (main_iterator::Current_ITE_step % screen_output_step == 0) {
				pf::PhaseFieldPoint phi_info = data_statistics_functions::statistical_phi();
				stringstream log;
				log << "#====================================================================================================" << endl;
				log << timer::return_cunrrent_time_by_string();
				log.setf(ios::fixed);
				log << "# Simulation step " << main_iterator::Current_ITE_step << " has been finished!" << endl;
				log << "# This " << screen_output_step << " steps used " << setprecision(3) << timer::interval_end(main_iterator::t_interval_begin) << "(secs.), " << endl;
				log << "# Total " << main_iterator::Current_ITE_step << " steps used " << setprecision(3) << timer::total_duration_sec(main_iterator::t_total_begin) << "(secs.)." << endl;
				log << "# Real simulation time is " << setprecision(3) << time_parameters::Real_Time << " (secs.)" << endl;
				log << "#----------------------------------------------------------------------------------------------------" << endl;
				for (int index = 0; index < phi_parameters::phi_number; index++)
					log << "#  Phase_" << index << "(" << materials_system::PHASES[phi_parameters::phi_property[index]] << ")" << " : "
					<< setprecision(5) << phi_info.phi[index] * 100 << " %" << endl;
				WriteLog(log.str());
				timer::interval_begin(main_iterator::t_interval_begin);
			}
		}

		inline void exec_pos_ii() {
			if (main_iterator::Current_ITE_step % screen_loop_step == 0 || main_iterator::Current_ITE_step % screen_output_step == 0) {
				stringstream log;
				log << "#------------------------------------------ PCT Field -----------------------------------------------" << endl;
				log << "# CURRENT STEP = " << main_iterator::Current_ITE_step << ", REAL TIME = " << time_parameters::Real_Time << endl;
				log << "# MAX PHI INCREMENT  = " << setprecision(5) << phi_parameters::PHI_MAX_VARIATION << endl;
				log << "# MAX CON INCREMENT  = " << setprecision(5) << con_parameters::CON_MAX_VARIATION << endl;
				log << "# MAX TEMP INCREMENT = " << setprecision(5) << temp_parameters::TEMP_MAX_VARIATION << endl;
				log << "#----------------------------------------------------------------------------------------------------" << endl;
				if (main_iterator::Current_ITE_step % screen_output_step != 0)
					log << endl << endl;
				WriteLog(log.str());
			}
		}

		// info end
		inline void exec_pos_iii() {
			if (main_iterator::Current_ITE_step % screen_output_step == 0) {
				stringstream log;
				log << "#====================================================================================================" << endl;
				log << endl << endl;
				WriteLog(log.str());
			}
		}

		inline void init() {
			InputFileReader::get_instance()->read_int_value("Solver.Output.LOG.loop_info_step", screen_loop_step, true);
			InputFileReader::get_instance()->read_int_value("Solver.Output.LOG.screen_output_step", screen_output_step, true);
			if (screen_loop_step < 1)
				screen_loop_step = INT_MAX;
			if (screen_output_step < 1)
				screen_output_step = INT_MAX;
			load_a_new_module(default_module_function, default_module_function, exec_pre_iii,
				default_module_function, default_module_function, default_module_function,
				exec_pos_i, exec_pos_ii, exec_pos_iii, default_module_function);
		}
	}
}