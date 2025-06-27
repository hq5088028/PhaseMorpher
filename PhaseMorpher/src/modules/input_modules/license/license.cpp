#include "license.h"
namespace pf{
	enum ACCOUNT { ROOT };
	static vector<string> _account = { "PhaseMorpher" };
	static vector<string> _password = { "PhaseMorpher" };
	static bool is_license_init = false;
	const int header_num = 1234;
	static char* license_header;
	namespace license_private {
		bool check_developer_counts() {
			string account, password;
			printf_color_on_control("> (PhaseMorpher) Account(press enter to confirm):");
			std::cout << endl;
			account = get_string_from_consol();
			printf_color_on_control("> (PhaseMorpher) Password(press enter to confirm):");
			std::cout << endl;
			password = get_string_from_consol(false, '*');
			bool check = false;
			for (int index = 0; index < _account.size(); index++)
				if (_account[index] == account && _password[index] == password)
					check = true;
			return check;
		}
		bool read_license(string name, string& cpu_id, ints_time& time1, ints_time& time2, ints_time& time3) {
			std::fstream fin(name, std::ios::binary | ios::in);
			if (!fin) {
				fin.close();
				return false;
			}
			char header[header_num];
			fin.read((char*)header, header_num);
			// generate base board code
			int cpu_idSize;
			fin.read((char*)&cpu_idSize, sizeof(int));
			char buff_str; cpu_id = "";
			for (int i = 0; i < cpu_idSize; i++) {
				fin.read((char*)&buff_str, sizeof(buff_str));
				cpu_id += buff_str;
			}
			fin.read((char*)&time1.year, sizeof(int));
			fin.read((char*)&time1.month, sizeof(int));
			fin.read((char*)&time1.day, sizeof(int));
			fin.read((char*)&time1.hour, sizeof(int));
			fin.read((char*)&time1.minute, sizeof(int));
			fin.read((char*)&time1.second, sizeof(int));

			fin.read((char*)&time2.year, sizeof(int));
			fin.read((char*)&time2.month, sizeof(int));
			fin.read((char*)&time2.day, sizeof(int));
			fin.read((char*)&time2.hour, sizeof(int));
			fin.read((char*)&time2.minute, sizeof(int));
			fin.read((char*)&time2.second, sizeof(int));

			fin.read((char*)&time3.year, sizeof(int));
			fin.read((char*)&time3.month, sizeof(int));
			fin.read((char*)&time3.day, sizeof(int));
			fin.read((char*)&time3.hour, sizeof(int));
			fin.read((char*)&time3.minute, sizeof(int));
			fin.read((char*)&time3.second, sizeof(int));
			return true;
		}
		bool write_license(string name, string& cpu_id, ints_time& time1, ints_time& time2, ints_time& time3) {
			ofstream fout(name, std::ios::binary);
			if (!fout) {
				fout.close();
				return false;
			}
			fout.write((const char*)license_header, header_num);
			// generate base board code
			int cpu_idSize = int(cpu_id.size());
			fout.write((char*)&cpu_idSize, sizeof(int));
			fout.write(cpu_id.data(), cpu_idSize);
			// generate license.dat time
			fout.write((const char*)&time1.year, sizeof(int));
			fout.write((const char*)&time1.month, sizeof(int));
			fout.write((const char*)&time1.day, sizeof(int));
			fout.write((const char*)&time1.hour, sizeof(int));
			fout.write((const char*)&time1.minute, sizeof(int));
			fout.write((const char*)&time1.second, sizeof(int));
			// current read/write license.dat time
			fout.write((const char*)&time2.year, sizeof(int));
			fout.write((const char*)&time2.month, sizeof(int));
			fout.write((const char*)&time2.day, sizeof(int));
			fout.write((const char*)&time2.hour, sizeof(int));
			fout.write((const char*)&time2.minute, sizeof(int));
			fout.write((const char*)&time2.second, sizeof(int));
			// valid time (month)
			fout.write((const char*)&time3.year, sizeof(int));
			fout.write((const char*)&time3.month, sizeof(int));
			fout.write((const char*)&time3.day, sizeof(int));
			fout.write((const char*)&time3.hour, sizeof(int));
			fout.write((const char*)&time3.minute, sizeof(int));
			fout.write((const char*)&time3.second, sizeof(int));
			fout.close();
			return true;
		}
		void write_current_time_into_license() {
			ints_time license_time, last_use_time, current_time, valid_time;
			string safe_base_board_code;
			bool find_license = read_license(license_path, safe_base_board_code, license_time, last_use_time, valid_time);
			if (find_license) {
				timer::get_cunrrent_time_by_ints_time(current_time);
				if (current_time > last_use_time || current_time == last_use_time) {
					write_license(license_path, safe_base_board_code, license_time, current_time, valid_time);
				}
			}
		}
		bool second_protection_read(string name, string& cpu_id, ints_time& time) {
			std::fstream fin(name, std::ios::binary | ios::in);
			if (!fin) {
				fin.close();
				return false;
			}
			char header[header_num];
			fin.read((char*)header, header_num);
			int cpu_idSize;
			fin.read((char*)&cpu_idSize, sizeof(int));
			char buff_str; cpu_id = "";
			for (int i = 0; i < cpu_idSize; i++) {
				fin.read((char*)&buff_str, sizeof(buff_str));
				cpu_id += buff_str;
			}
			fin.read((char*)&time.year, sizeof(int));
			fin.read((char*)&time.month, sizeof(int));
			fin.read((char*)&time.day, sizeof(int));
			fin.read((char*)&time.hour, sizeof(int));
			fin.read((char*)&time.minute, sizeof(int));
			fin.read((char*)&time.second, sizeof(int));
			return true;
		}
		bool second_protection_write(string name, string& cpu_id, ints_time& time) {
			ofstream fout(name, std::ios::binary);
			if (!fout) {
				fout.close();
				return false;
			}
			fout.write((const char*)license_header, header_num);
			int cpu_idSize = int(cpu_id.size());
			fout.write((char*)&cpu_idSize, sizeof(int));
			fout.write(cpu_id.data(), cpu_idSize);
			fout.write((const char*)&time.year, sizeof(int));
			fout.write((const char*)&time.month, sizeof(int));
			fout.write((const char*)&time.day, sizeof(int));
			fout.write((const char*)&time.hour, sizeof(int));
			fout.write((const char*)&time.minute, sizeof(int));
			fout.write((const char*)&time.second, sizeof(int));
			fout.close();
			return true;
		}
	}
	namespace license {
		bool is_license() {
			return is_license_init;
		}
		bool check_license(bool _debug) {
			is_license_init = false;
			//#ifdef _WIN32
					// try to open license
			ints_time license_time, last_use_time, current_time, valid_time;
			string safe_cpu_id, cpu_id;
			if (license_private::read_license(license_path, safe_cpu_id, license_time, last_use_time, valid_time)) {
				timer::get_cunrrent_time_by_ints_time(current_time);
				cpu_id = get_cpuId();
				if (safe_cpu_id.compare(cpu_id) == 0) {
					if (last_use_time > current_time || license_time > current_time) {
						if (_debug) {
							printf_color_on_control("> (license) Time mismatch ! Don't change the local time of your computer/server !");
							std::cout << endl;
						}
						return false;
					}
				}
				else {
					if (_debug) {
						printf_color_on_control("> (license) Cpu ID mismatch ! Don't use the license of another computer/server !");
						std::cout << endl;
					}
					return false;
				}
				license_time = license_time + valid_time;
				if (license_time > current_time) {
					license_private::write_current_time_into_license();
					is_license_init = true;
					return true;
				}
				else {
					return false;
				}
			}
			return false;
		}
		void print_license_info() {
			stringstream info;
			ints_time license_time, last_use_time, current_time, valid_time;
			string safe_base_board_code;
			if (license_private::read_license(license_path, safe_base_board_code, license_time, last_use_time, valid_time)) {
				printf_color_on_control("> (license) license has been generated on " + timer::trans_int_time_to_string(license_time));
				cout << endl;
				printf_color_on_control("> (license) license was last used on " + timer::trans_int_time_to_string(last_use_time));
				cout << endl;
				printf_color_on_control("> (license) license valid until " + timer::trans_int_time_to_string(license_time + valid_time));
				cout << endl;
				printf_color_on_control("> (license) valid cpu ID of this license is " + safe_base_board_code);
				cout << endl;
			}
			else {
				printf_color_on_control("> (license) license isn't there !");
				cout << endl;
			}
		}
		bool generate_a_license_for_this_computer() {
			ints_time valid_months;
			if (license_private::check_developer_counts()) {
				printf_color_on_control("> (license) The length of time for the license is active (integer type, and press enter to confirm):");
				cout << endl;
				printf_color_on_control("> (license) Years:");
				cin >> valid_months.year;
				printf_color_on_control("> (license) Months:");
				cin >> valid_months.month;
				printf_color_on_control("> (license) Days:");
				cin >> valid_months.day;
				printf_color_on_control("> (license) Hours:");
				cin >> valid_months.hour;
				valid_months.minute = 0;
				valid_months.second = 0;
				ints_time t;
				timer::get_cunrrent_time_by_ints_time(t);
				string safe_cpu_id = get_cpuId();
				bool write = license_private::write_license(license_path, safe_cpu_id, t, t, valid_months);
				if (write) {
					printf_color_on_control("> (license) A license is generated successfully !");
					is_license_init = true;
					cout << endl;
				}
				return write;
			}
			printf_color_on_control("> (license) Account/password error ! Generate license failed !");
			cout << endl;
			return false;
		}
		bool generate_a_license_for_other_computer(string cpu_id) {
			ints_time valid_months;
			if (license_private::check_developer_counts()) {
				printf_color_on_control("> (license) The length of time for the license is active (integer type, and press enter to confirm):");
				cout << endl;
				printf_color_on_control("> (license) Years:");
				cin >> valid_months.year;
				printf_color_on_control("> (license) Months:");
				cin >> valid_months.month;
				printf_color_on_control("> (license) Days:");
				cin >> valid_months.day;
				printf_color_on_control("> (license) Hours:");
				cin >> valid_months.hour;
				valid_months.minute = 0;
				valid_months.second = 0;
				ints_time t;
				timer::get_cunrrent_time_by_ints_time(t);
				bool write = license_private::write_license(license_path, cpu_id, t, t, valid_months);
				if (write) {
					printf_color_on_control("> (license) A license is generated successfully !");
					cout << endl;
				}
				return write;
			}
			printf_color_on_control("> (license) Account/password error ! Generate license failed !");
			cout << endl;
			return false;
		}
		void init_license() {
			license_header = new char[header_num];
			std::string head = "## PM PROGRAM FILE ##";
			memcpy(license_header, head.c_str(), header_num);
			if (!check_license()) {
				printf_color_on_control("> (license) Activate license failed !");
				cout << endl;
			}
			else {
				printf_color_on_control("> (license) Activate license success !");
				cout << endl;
			}
		}
		void deinit_license() {
			license_private::write_current_time_into_license();
		}
	}
}