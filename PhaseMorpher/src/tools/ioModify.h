/*

Created:     Qi Huang 2025.03

Modified:    Qi Huang 2025.03;

*/
#pragma once
// #include "../base/sysfiles.h"
#include <conio.h>
#include <iostream>
#include <fstream>
namespace pf {
	//printf("\x1b[%d;%dm%s\x1b[%dm", backcolor, frountcolor, str, control);
	/*��һ��% d:backcolor��ʾ��ʾ�ַ����ı�����ɫ, ��ֵ���£�
	40 : ��
	41 : ��
	42 : ��
	43 : ��
	44 : ��
	45 : ��
	46 : ����
	47 : ��ɫ
	�ڶ���% d : frountcolor��ʾ������ɫ, ��ֵ���£�
	30 : ��
	31 : ��
	32 : ��
	33 : ��
	34 : ��
	35 : ��
	36 : ����
	37 : ��ɫ

	������ % s : str ��ʾ��Ҫ��ʾ���ַ���

	���ĸ� % d : control��ʾANSI������, ��ֵ���±���ʾ��

	ANSI������:

	QUOTE:
	\x1b[0m     �ر���������
	\x1b[1m     ���ø�����
	\x1b[4m     �»���
	\x1b[5m     ��˸
	\x1b[7m     ����
	\x1b[8m     ����
	\x1b[30m--  \x1b[37m   ����ǰ��ɫ
	\x1b[40m--  \x1b[47m   ���ñ���ɫ
	\x1b[nA    �������n��
	\x1b[nB    �������n��
	\x1b[nC    �������n��
	\x1b[nD    �������n��
	\x1b[y; xH  ���ù��λ��
	\x1b[2J    ����
	\x1b[K     ����ӹ�굽��β������
	\x1b[s     ������λ��
	\x1b[u     �ָ����λ��
	\x1b[? 25l  ���ع��
	\x1b[? 25h  ��ʾ���
	����
	printf("\x1b[%d;%dmhello world\n\x1b[0m",i, j);
	*/
	inline void printf_color_on_control(std::string str, int front_color = 30, int back_color = 43) {
		printf("\x1b[%d;%dm%s\x1b[0m", back_color, front_color, str.c_str());
	}

	inline char get_char_not_show() {
		char c;
#ifdef _WIN32
		c = _getch();
#else
		system("stty -echo");
		c = getchar();
		system("stty echo");
#endif
		return c;
	}

	inline std::string get_string_from_consol(bool is_show = true, char replace_char = '*') {
		std::string str = "";
		int i = 0;
		if (is_show) {
			while (true) {
				char ch = getchar();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
			}
		}
		else {
			while (true) {
				char ch = get_char_not_show();
				if (ch == '\r' || ch == '\n') {
					break;
				}
				str.push_back(ch);
				putchar(replace_char);
			}
		}
		putchar('\n');
		return str;
	};

	inline void write_string_to_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path);
		if (!fout) {
			std::cout << "Failed to write the txt file!" << std::endl;
			fout.close();
			return;
		}
		fout << content << std::endl;
		fout.close();
	}

	inline void add_string_to_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path, std::ios::app);
		if (!fout) {
			std::cout << "Failed to add the string to txt file!" << std::endl;
			fout.close();
			return;
		}
		fout << content;
		fout.close();
	}

	inline void add_string_to_screen_and_file(std::string content, std::string file_path) {
		std::ofstream fout(file_path, std::ios::app);
		if (!fout) {
			std::cout << "Failed to add the string to txt file!" << std::endl;
			fout.close();
			return;
		}
		std::cout << content;
		fout << content;
		fout.close();
	}

}