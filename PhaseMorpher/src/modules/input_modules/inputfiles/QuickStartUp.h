/*

Created:     Qi Huang 2025.3.21

Modified:    Qi Huang 2025.3.21;

Email:       

*/

#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <filesystem>
#if !defined(_WIN32)
#include <readline/readline.h>
#include <readline/history.h>
#endif //_WIN32
#include "../../../tools/StringModify.h"
#include "../license/ProgramPath.h"
namespace pf {
	inline bool Quick_StartUp(std::string& infile_path) {
		std::ifstream quick_simu_info(program_path / "start.in");
		if (quick_simu_info) {
			infile_path = "";
			getline(quick_simu_info, infile_path);
			quick_simu_info.close();
			std::filesystem::path p_last_path(infile_path);
			if (std::filesystem::exists(infile_path)) {
				std::filesystem::remove(program_path / "start.in");
				return true;
			}
		}
		return false;
	}
}