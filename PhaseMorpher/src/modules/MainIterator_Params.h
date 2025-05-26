#pragma once
namespace pf {

	namespace main_iterator {
		inline int ITE_Begin_Step = 0;
		inline int ITE_End_Step = 0;
		inline int OpenMP_Thread_Counts = 1;
		inline int Current_ITE_step = 0;

		inline double t_interval_modules_init{};
		inline double t_interval_modules_pre_exec{};
		inline double t_interval_modules_exec{};
		inline double t_interval_modules_pos_exec{};
		inline double t_interval_modules_deinit{};

		inline double t_total_begin{};
		inline double t_interval_begin{};
	}

}