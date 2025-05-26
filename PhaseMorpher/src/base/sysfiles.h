/*

Created:     Qi Huang 2025.03

Modified:    Qi Huang 2025.03

Email:       qihuang0908@163.com

*/

#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <omp.h>
#include <vector>
#include <cmath>
#include <string>
#include <time.h>
#include <iostream>
#include <cstring>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <initializer_list>
#include <complex>
#include <stdio.h>
#include <float.h>
#include <thread>
#include <chrono>
#include <vector>
#include <random>
#include <array>
#ifdef _WIN32
#define SYS_PROGRAM_STOP while(true){(void)getchar();};
//#define SYS_PROGRAM_STOP std::exit(0);
#include <tchar.h>
#include <shlobj.h>
#include <direct.h>
#include <Windows.h>
#undef max // 取消 max 宏定义
#undef min // 取消 min 宏定义
#include <io.h>
#include <conio.h>
#include <psapi.h>
#include <process.h>
#include <intrin.h>
#define dirSeparator std::string("\\")                                     //< Windows style directory separator
#else
#include <sys/io.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <unistd.h>
#include <cpuid.h>
#include <sys/types.h>
#define SYS_PROGRAM_STOP std::exit(0);
#define dirSeparator std::string("/")                                      //< Unix/Linux style directory separator
#endif

// control the REAL data
// #define USE_DOUBLE

#ifdef USE_DOUBLE
using REAL = double; // use double

#define SYS_EPSILON   (0.000001)
#define SYS_EPSILON_R (0.999999)

#define Phi_Num_Cut_Off   (0.001)
#define Phi_Num_Cut_Off_R (0.999)

#define PhiCon_Num_Cut_Off   (0.01)
#define PhiCon_Num_Cut_Off_R (0.99)

inline double NaN() { return std::numeric_limits<double>::max(); };
inline double REAL_MAX() { return std::numeric_limits<double>::max(); };

#define PI (3.1415926535897932)

#define AngleToRadians(angle) double(angle/180.0*PI)

#else
using REAL = float;  // use float

#define SYS_EPSILON   (0.000001f)
#define SYS_EPSILON_R (0.999999f)

#define Phi_Num_Cut_Off   (0.001f)
#define Phi_Num_Cut_Off_R (0.999f)

#define PhiCon_Num_Cut_Off   (0.01f)
#define PhiCon_Num_Cut_Off_R (0.99f)

inline float NaN() { return std::numeric_limits<float>::max(); };

inline float REAL_MAX() { return std::numeric_limits<float>::max(); };

#define PI (3.1415927f)

#define AngleToRadians(angle) float(angle/180.0f*PI)

#endif

#define MESH_INDEX(x,y,z,Nx,Ny) (x+y*Nx+z*Nx*Ny)
#define MESH_BC_X_INDEX(y,z,Ny) (y+z*Ny)
#define MESH_BC_Y_INDEX(x,z,Nx) (x+z*(Nx+1))
#define MESH_BC_Z_INDEX(x,y,Nx) (x+y*(Nx+1))




//   std::random_device rd; // 高质量随机数种子生成器
//   // std::mt19937 gen(static_cast<unsigned int>(std::time(nullptr))); // 时间戳种子：static_cast<unsigned int>(std::time(nullptr))，或固定的值 int
//   std::mt19937 gen(rd()); // 使用 Mersenne Twister 引擎初始化随机数生成器
//   std::uniform_real_distribution<> real_dist(0.0, 1.0); // [0.0, 1.0] 范围内的浮点数
//   // std::uniform_int_distribution<> int_dist(1, 100); // [1, 100] 范围内的整数
//   // std::normal_distribution<> normal_dist(50.0, 10.0); // 正态分布，均值 50，标准差 10
//   REAL rand = real_dist(gen);  // random