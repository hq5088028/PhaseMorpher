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
	/*第一个% d:backcolor表示显示字符串的背景颜色, 其值如下：
	40 : 黑
	41 : 红
	42 : 绿
	43 : 黄
	44 : 蓝
	45 : 紫
	46 : 深绿
	47 : 白色
	第二个% d : frountcolor表示字体颜色, 其值如下：
	30 : 黑
	31 : 红
	32 : 绿
	33 : 黄
	34 : 蓝
	35 : 紫
	36 : 深绿
	37 : 白色

	第三个 % s : str 表示需要显示的字符串

	第四个 % d : control表示ANSI控制码, 其值如下表所示：

	ANSI控制码:

	QUOTE:
	\x1b[0m     关闭所有属性
	\x1b[1m     设置高亮度
	\x1b[4m     下划线
	\x1b[5m     闪烁
	\x1b[7m     反显
	\x1b[8m     消隐
	\x1b[30m--  \x1b[37m   设置前景色
	\x1b[40m--  \x1b[47m   设置背景色
	\x1b[nA    光标上移n行
	\x1b[nB    光标下移n行
	\x1b[nC    光标右移n行
	\x1b[nD    光标左移n行
	\x1b[y; xH  设置光标位置
	\x1b[2J    清屏
	\x1b[K     清除从光标到行尾的内容
	\x1b[s     保存光标位置
	\x1b[u     恢复光标位置
	\x1b[? 25l  隐藏光标
	\x1b[? 25h  显示光标
	例：
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