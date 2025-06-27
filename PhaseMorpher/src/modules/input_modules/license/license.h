#pragma once
#include "../../../tools/timer.h"
#include "../../../tools/ioModify.h"
#include "computer_info.h"
#include "ProgramPath.h"
#include <sstream>
namespace pf {
	using namespace std;
	const string license_path = program_path.string() + "PhaseMorpher.license";
	namespace license {
		void init_license();
		void deinit_license();
		bool check_license(bool _debug = true);
		bool generate_a_license_for_this_computer();
		bool generate_a_license_for_other_computer(string cpu_id);
		void print_license_info();
		//bool get_base_board_info(string& base_board_code);
		//bool get_cpu_info(string& cpu_id);
		bool is_license();
		inline void print_current_time_on_screen() {
			timer::print_cunrrent_time_on_screen();
		}
		inline void clean_cmd() {
#ifdef _WIN32
			system("cls");
#else
			system("clear");
#endif
		}
	}
}