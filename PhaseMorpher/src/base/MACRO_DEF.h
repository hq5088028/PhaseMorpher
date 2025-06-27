/*

Created:     Qi Huang 2025.03

Modified:    Qi Huang 2025.03

Email:       qihuang0908@163.com

*/
#pragma once
#include <limits>
#define _CRT_SECURE_NO_WARNINGS

// control the REAL data
// #define USE_DOUBLE

#undef max // 取消 max 宏定义
#undef min // 取消 min 宏定义


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


#define SYS_PROGRAM_STOP while(true){(void)getchar();};

#define dirSeparator std::string("\\")                                     //< Windows style directory separator