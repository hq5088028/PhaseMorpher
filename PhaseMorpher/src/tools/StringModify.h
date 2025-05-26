#pragma once
#include "../base/sysfiles.h"
namespace pf {
	inline void str_char_delete(std::string& str, char _word) {
		std::size_t start_position{ str.find(_word) };
		while (start_position != std::string::npos) {
			str.erase(start_position, 1);
			start_position = str.find(_word);
		}
	}
	inline string GetFolderOfPath(std::string file_path) {
		while (file_path.size() != 0)
		{
			char _c = *(file_path.end() - 1);
			if (_c == '/' || _c == '\\') {
				file_path.erase(file_path.end() - 1);
				break;
			}
			else {
				file_path.erase(file_path.end() - 1);
			}
		}
		return file_path;
	}
	inline string GetFileNameOfPath(std::string file_path) {
		string name = "";
		int index = 1;
		while (index <= file_path.size())
		{
			char _c = *(file_path.end() - index);
			if (_c == '/' || _c == '\\') {
				break;
			}
			else {
				index++;
				name = _c + name;
			}
		}
		return name;
	}
	inline string erase_tail_of_infile(string input_file_name) {
		string tail = ".pm", name_without_tail = input_file_name;
		int index = 0;
		bool is_name_correct = true;
		while (name_without_tail.size() != 0 && index < tail.size())
		{
			index++;
			if (*(name_without_tail.end() - 1) != *(tail.end() - index)) {
				is_name_correct = false;
				break;
			}
			name_without_tail.erase(name_without_tail.end() - 1);
		}
		if (name_without_tail.size() == 0)
			is_name_correct = false;
		if (!is_name_correct) {
			cout << "> input file name error, file name = " << input_file_name << ", aim tail is " << tail << endl;
			SYS_PROGRAM_STOP;
		}
		return name_without_tail;
	}
}