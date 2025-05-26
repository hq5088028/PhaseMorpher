/*

Created:     Qi Huang 2025.03

Modified:    Qi Huang 2025.03;

*/
#pragma once
#include "../Module.h"
#include "inputfiles/UserStartUp.h"
#include "inputfiles/QuickStartUp.h"
#include "inputfiles/infile_selector.h"
#include "ioFiles_Params.h"
#include "inputfiles/InputFileReader.h"
#include "license/license.h"
namespace pf {
	static void init_input_modules() {
		// license init
		license::init_license();
		// debug input file
		bool input_file_debug = false;
		InputFileReader::get_instance()->read_bool_value("InputFile.debug", input_file_debug, true);
		stringstream out; 
		out << "============================== Parameters Definition Format =============================" << endl;
		out << "                             Parameter_Name = Parameter_Value" << endl;
		out << "======================================= M a c r o =======================================" << endl;
		out << "1.tube               " << ", $TUBE[1,10,2]$ = 1,3,5,7,9" << endl;
		out << "2.rand               " << ", $RAND_INT[1,10]$ = 1 - 10" << endl;
		out << "                     " << ", $RAND_REAL[1,10]$ = 1.000000 - 10.000000" << endl;
		out << "========================= Define Custom Variables and Functions =========================" << endl;
		out << "Define.Var = A,0.1" << endl;
		out << "Define.Func = ABC@{[A*pow(A,2)]}@" << endl;
		out << "# default functions      : \"pow(val, ord)\", \"sqrt(val)\", \"abs(val)\", \"exp(val)\", \"ln(val)\"," << endl;
		out << "#                          \"log(base_val, val)\", \"sin(val)\", \"cos(val)\", \"tan(val)\"," << endl;
		out << "#                          \"asin(val)\", \"acos(val)\", \"atan(val)\", \"cos(val)\", " << endl;
		out << "#                          \"tan(val)\", \"asin(val)\", \"acos(val)\", \"atan(val)\"" << endl;
		out << "=========================================================================================" << endl;
		add_string_to_file(out.str(), input_output_files_parameters::DebugFile_Path);
		if (input_file_debug) {
			InputFileReader::get_instance()->debug_infile_and_valid_words();
			InputFileReader::get_instance()->debug_custom_variavle_and_funcs();
		}
		out.str("");
		out << "=========================================================================================" << endl;
		out << "================================= Parameters Definition =================================" << endl;
		out << "=========================================================================================" << endl;
		add_string_to_file(out.str(), input_output_files_parameters::DebugFile_Path);
		WriteLog("> MODULE INIT : Input File (.pm) Ready !\n");
	}
}