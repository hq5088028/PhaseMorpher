#pragma once
#include "../../../tools/ioModify.h"
#include "../license/license.h"
namespace pf {
	inline char input_single_char();
	inline void license_info();
	inline void permissions_info();
	inline void cpu_id_info();
	inline void about_info();
	inline void pm_selection();
	inline void pm_info();

	inline char input_single_char() {
#ifdef _READLINE_H_ 
		char* c_str{ readline("") };
		while (strlen(c_str) != 1) {
			pf::printf_color_on_control("only one character accepted, please re-enter."); std::cout << std::endl;
			free(c_str);
			c_str = readline("");
		}
		char single_char{ c_str[0] };
		free(c_str);
#else
		std::string str{};
		std::getline(std::cin, str);
		while (str.length() != 1) {
			pf::printf_color_on_control("only one character accepted, please re-enter."); std::cout << std::endl;
			std::getline(std::cin, str);
		}
		char single_char{ str.at(0) };
#endif
		return single_char;
	}
	inline void license_info() {
		license::check_license(false);
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 1 PhaseMorpher license information :");
		if (license::is_license())
			pf::printf_color_on_control("     activated", 32);
		else
			pf::printf_color_on_control("     inactivated", 31);
		std::cout << std::endl;
		std::cout << std::endl;
		license::print_license_info();
		std::cout << std::endl;
		pf::printf_color_on_control("> press \"G/g\" to generate a license for this computer/server;");
		std::cout << std::endl;
		pf::printf_color_on_control("> press \"H/h\" to generate a license for other computer/server;");
		std::cout << std::endl;
		pf::printf_color_on_control("> press enter to confirm;");
		std::cout << std::endl;
		pf::printf_color_on_control("> press any other keys to go back to previous selection.");
		std::cout << std::endl;
		char option{ input_single_char() }; //get enter
		switch (option) {
		case 'g':
		case 'G': {
			license::generate_a_license_for_this_computer();
			break;
		}
		case 'h':
		case 'H': {
			pf::printf_color_on_control("> write the cpu ID here:");
			std::string base_board_code;
			std::cin >> base_board_code;
			char enter = getchar(); //get enter
			license::generate_a_license_for_other_computer(base_board_code);
			break;
		}
		default: {
			pm_selection();
			break;
		}
		}
	}
	inline void permissions_info() {
		license::check_license();
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 2 PhaseMorpher user permissions :");
		std::cout << "\n" << std::endl;
		if (license::is_license()) {
			pf::printf_color_on_control("> You get ");
			pf::printf_color_on_control("full access", 32);
			pf::printf_color_on_control(" to PhaseMorpher.");
		}
		else {
			pf::printf_color_on_control("> You cannot use the ");
			pf::printf_color_on_control("parallel", 31);
			pf::printf_color_on_control(" function in PhaseMorpher.");
		}
		std::cout << "\n" << std::endl;
		pf::printf_color_on_control("> press any other keys to go back to previous selection.");
		std::cout << std::endl;
		char option{ input_single_char() }; //get enter
		pm_selection();
	}
	inline void cpu_id_info() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 3 PhaseMorpher cpu ID :");
		std::cout << "\n" << std::endl;
		std::string cpu_id = get_cpuId();
		pf::printf_color_on_control("> cpu ID of this computer is:     ");
		pf::printf_color_on_control(cpu_id, 37, 42);
		std::cout << "\n" << std::endl;
		pf::printf_color_on_control("> press any other keys to go back to previous selection.");
		std::cout << std::endl;

		char option{ input_single_char() }; //get enter
		pm_selection();
	}
	inline void about_info() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------", 34);
		std::cout << std::endl;
		pf::printf_color_on_control("> 4 PhaseMorpher about :"); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> Developers:");
		pf::printf_color_on_control("[1] Qi Huang; [2] Zihang Wang;", 34); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> Emails:");
		pf::printf_color_on_control("[1] qihuang0908@163.com; [2] w.zihang@qq.com;", 34); std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("    へ　　　　　／|"); std::cout << std::endl;
		pf::printf_color_on_control("　　/＼7　　　 ∠＿/"); std::cout << std::endl;
		pf::printf_color_on_control("　 /　│　　 ／　／"); std::cout << std::endl;
		pf::printf_color_on_control("　│　Z ＿,＜　／　　 /`ヽ"); std::cout << std::endl;
		pf::printf_color_on_control("　│　　　　　ヽ　　 /　　〉"); std::cout << std::endl;
		pf::printf_color_on_control("　 Y　　　　　`　 /　　/"); std::cout << std::endl;
		pf::printf_color_on_control("　/●　 　●　　 |〈　　/"); std::cout << std::endl;
		pf::printf_color_on_control("　() へ　　（） |　＼〈"); std::cout << std::endl;
		pf::printf_color_on_control("　　>   _　 ィ　 │ ／／"); std::cout << std::endl;
		pf::printf_color_on_control("　 / へ　　 /　/＜| ＼＼       Pikachu says: let's go ! it's time to start a simulation !"); std::cout << std::endl;
		pf::printf_color_on_control("　 ヽ_/　　(_／　 │／／"); std::cout << std::endl;
		pf::printf_color_on_control("　　7　　　　　　　|／"); std::cout << std::endl;
		pf::printf_color_on_control("　　＞―r￣￣` ―＿) "); std::cout << std::endl;
		std::cout << std::endl;
		pf::printf_color_on_control("> press any other keys to go back to previous selection.");
		std::cout << std::endl;
		char option{ input_single_char() }; //get enter
		pm_selection();
	}
	inline void pm_selection() {
		pf::printf_color_on_control("------------------------------------------------------------------------------------------------------\n", 34);
		license::check_license(false);
		pf::printf_color_on_control("> Select the one you want (press corresponding number and press enter):");
		std::cout << "\n" << std::endl;

		// 1 license information
		pf::printf_color_on_control("> 1 PhaseMorpher: license information :");
		license::check_license(false);
		if (license::is_license())
			pf::printf_color_on_control("     activated", 32);
		else
			pf::printf_color_on_control("     inactivated", 31);
		std::cout << "\n" << std::endl;

		// 2 User permissions
		pf::printf_color_on_control("> 2 PhaseMorpher: user permissions ?");
		std::cout << "\n" << std::endl;

		// 3 cpu ID of this computer
		pf::printf_color_on_control("> 3 PhaseMorpher: cpu ID ?");
		std::cout << "\n" << std::endl;

		// 4 about copyright, developers & e-mail
		pf::printf_color_on_control("> 4 PhaseMorpher: about ?");
		std::cout << "\n" << std::endl;

		// 5 back to previous menu
		pf::printf_color_on_control("> 5 back to previous menu");
		std::cout << "\n" << std::endl;

		// 6 quit MInDes
		pf::printf_color_on_control("> 6 Quit PhaseMorpher");
		std::cout << "\n" << std::endl;


		bool end_selection{ false };
		while (!end_selection) {
			char option{ input_single_char() };
			switch (option) {
			case '1': {
				license_info();
				end_selection = true;
				break;
			}
			case '2': {
				permissions_info();
				end_selection = true;
				break;
			}
			case '3': {
				cpu_id_info();
				end_selection = true;
				break;
			}
			case '4': {
				about_info();
				end_selection = true;
				break;
			}
			case '5': {
				end_selection = true;
				break;
			}
			case '6': {
				exit(0);
			}
			default: {
				std::cout << "invalid option, please re-enter an option." << std::endl;
				end_selection = false;
				break;
			}
			}
		}
	}
	inline void pm_info() {
		license::init_license();
		printf_color_on_control("------------------------------------------------------------------------------------------------------\n", 34);
		printf_color_on_control(">>> >  >    >      >          >        PhaseMorpher   INFORMATION         <         <      <    <  < <<<", 34);
		cout << endl;
		pm_selection();
	}
}