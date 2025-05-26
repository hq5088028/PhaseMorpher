#pragma once
#include<iostream>
#include<fstream>

namespace pf {
	inline void print_info() {
		std::cout << "" << std::endl;
		std::cout << " @@@@@@@@@@@@    @@@@@@        @@@@@@ " << std::endl;
		std::cout << " @@@@@@@@@@@@    @@@@@@        @@@@@@ " << std::endl;
		std::cout << " @@@@@@@@@@@    @@@@@@@@      @@@@@@@@ " << std::endl;
		std::cout << " @@@@@@@@@@@    @@@@@@@@      @@@@@@@@ " << std::endl;
		std::cout << " @@@@@@@@@@    @@@@@@@@@@    @@@@@@@@@@  " << std::endl;
		std::cout << " @@@@@@        @@@@@@@@@@    @@@@@@@@@@@ " << std::endl;
		std::cout << " @@@@@@       @@@@@@@@@@@@@@@@@@@@@@@@@@ " << std::endl;
		std::cout << " @@@@@@       @@@@@@@@@@@@@@@@@@@@@@@@@@ " << std::endl;
		std::cout << " @@@@@@      @@@@@@@@@@@@@@@@@@@@@@@@@@@ " << std::endl;
		std::cout << " @@@@@@      @@@@@@@@@@@@@@@@@@@@@@@@@@@ " << std::endl;
		std::cout << "" << std::endl;
	};
	inline void print_help_doc() {
		std::cout << "************************************ HELP DOC ************************************" << std::endl;
		std::cout << "-H      this help" << std::endl;
		std::cout << "-I      Infomation" << std::endl;
		std::cout << "-D      Print current working directory" << std::endl;
		std::cout << "-L      Use the last simu's infomation" << std::endl;
		std::cout << "-C      Choose a input file" << std::endl;
		std::cout << "\"path\"  input the path to start simu; if path contains space, please quote it" << std::endl;
		std::cout << "-S      Sequentially simus (default) (under progress)" << std::endl;
		std::cout << "-P      Parallelly simus, follows with a number (under progress)" << std::endl;
		std::cout << "-Q      Exit the program" << std::endl;
		std::cout << "**********************************************************************************" << std::endl;
	};
}
