#pragma once
#include <windows.h>
#include <intrin.h>
namespace pf {
#ifdef _WIN32
	// char[128]
	inline BOOL GetBaseBoardByCmd(char* lpszBaseBoard, int len/*=128*/)
	{
		const long MAX_COMMAND_SIZE = 10000; // 命令行输出缓冲大小	
		WCHAR szFetCmd[] = L"wmic BaseBoard get SerialNumber"; // 获取主板序列号命令行	
		const std::string strEnSearch = "SerialNumber"; // 主板序列号的前导信息

		BOOL   bret = FALSE;
		HANDLE hReadPipe = NULL; //读取管道
		HANDLE hWritePipe = NULL; //写入管道	
		PROCESS_INFORMATION pi;   //进程信息	
		STARTUPINFO			si;	  //控制命令行窗口信息
		SECURITY_ATTRIBUTES sa;   //安全属性

		char			szBuffer[MAX_COMMAND_SIZE + 1] = { 0 }; // 放置命令行结果的输出缓冲区
		std::string			strBuffer;
		unsigned long	count = 0;
		long			ipos = 0;

		memset(&pi, 0, sizeof(pi));
		memset(&si, 0, sizeof(si));
		memset(&sa, 0, sizeof(sa));

		pi.hProcess = NULL;
		pi.hThread = NULL;
		si.cb = sizeof(STARTUPINFO);
		sa.nLength = sizeof(SECURITY_ATTRIBUTES);
		sa.lpSecurityDescriptor = NULL;
		sa.bInheritHandle = TRUE;

		//1.0 创建管道
		bret = CreatePipe(&hReadPipe, &hWritePipe, &sa, 0);
		if (!bret)
		{
			//goto END;
			return false;
		}

		//2.0 设置命令行窗口的信息为指定的读写管道
		GetStartupInfo(&si);
		si.hStdError = hWritePipe;
		si.hStdOutput = hWritePipe;
		si.wShowWindow = SW_HIDE; //隐藏命令行窗口
		si.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;

		//3.0 创建获取命令行的进程
		bret = CreateProcess(NULL, szFetCmd, NULL, NULL, TRUE, 0, NULL, NULL, &si, &pi);
		if (!bret)
		{
			//goto END;
			return false;
		}

		//4.0 读取返回的数据
		WaitForSingleObject(pi.hProcess, 500/*INFINITE*/);
		bret = ReadFile(hReadPipe, szBuffer, MAX_COMMAND_SIZE, &count, 0);
		if (!bret)
		{
			//goto END;
			return false;
		}

		//5.0 查找主板序列号
		bret = FALSE;
		strBuffer = szBuffer;
		ipos = long(strBuffer.find(strEnSearch));

		if (ipos < 0) // 没有找到
		{
			//goto END;
			return false;
		}
		else
		{
			strBuffer = strBuffer.substr(ipos + strEnSearch.length());
		}

		memset(szBuffer, 0x00, sizeof(szBuffer));
		strcpy_s(szBuffer, strBuffer.c_str());

		//去掉中间的空格 \r \n \0 \t
		int j = 0;
		for (int i = 0; i < strlen(szBuffer); i++)
		{
			if (szBuffer[i] != ' ' && szBuffer[i] != '\n' && szBuffer[i] != '\r' && szBuffer[i] != '\t' && szBuffer[i] != '\0')
			{
				lpszBaseBoard[j] = szBuffer[i];
				j++;
			}
		}

		bret = TRUE;

		//关闭所有的句柄
		CloseHandle(hWritePipe);
		CloseHandle(hReadPipe);
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);

		return(bret);
	}

	inline void getcpuidex(unsigned int* CPUInfo, unsigned int InfoType, unsigned int ECXValue)
	{
#if defined(_MSC_VER) // MSVC  
#if defined(_WIN64) // 64位下不支持内联汇编. 1600: VS2010, 据说VC2008 SP1之后才支持__cpuidex.  
		__cpuidex((int*)(void*)CPUInfo, (int)InfoType, (int)ECXValue);
#else  
		if (NULL == CPUInfo)
			return;
		_asm {
			// load. 读取参数到寄存器.  
			mov edi, CPUInfo;
			mov eax, InfoType;
			mov ecx, ECXValue;
			// CPUID  
			cpuid;
			// save. 将寄存器保存到CPUInfo  
			mov[edi], eax;
			mov[edi + 4], ebx;
			mov[edi + 8], ecx;
			mov[edi + 12], edx;
		}
#endif  
#endif  
	}
#endif
	inline void getcpuid(unsigned int* CPUInfo, unsigned int InfoType)
	{
#if defined(__GNUC__)// GCC  
		__cpuid(InfoType, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
#elif defined(_MSC_VER)// MSVC  
		getcpuidex(CPUInfo, InfoType, 0);
#endif  
	}
	//static void get_cpuId(char* pCpuId)
	//{

	//	int dwBuf[4];
	//	getcpuid((unsigned int*)dwBuf, 1);
	//	sprintf(pCpuId, "%08X", dwBuf[3]);
	//	sprintf(pCpuId + 8, "%08X", dwBuf[0]);
	//	return;
	//}

	inline std::string get_cpuId() {
		int dwBuf[4];

		// 获取CPU ID
#ifdef _WIN32
		__cpuid(dwBuf, 1); // Windows平台使用__cpuid
#else
		__cpuid(1, dwBuf[0], dwBuf[1], dwBuf[2], dwBuf[3]); // Linux平台使用__cpuid
#endif

		// 格式化CPU ID
		char cpuId[17] = { 0 }; // 16字符 + 1终止符
#ifdef _WIN32
		sprintf_s(cpuId, sizeof(cpuId), "%08X", dwBuf[3]); // Windows平台使用sprintf_s
		sprintf_s(cpuId + 8, sizeof(cpuId) - 8, "%08X", dwBuf[0]);
#else
		snprintf(cpuId, sizeof(cpuId), "%08X", dwBuf[3]); // Linux平台使用snprintf
		snprintf(cpuId + 8, sizeof(cpuId) - 8, "%08X", dwBuf[0]);
#endif

		return std::string(cpuId); // 返回std::string
	}

}