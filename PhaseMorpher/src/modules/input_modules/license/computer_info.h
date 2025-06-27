#pragma once
#include <windows.h>
#include <intrin.h>
namespace pf {
#ifdef _WIN32
	// char[128]
	inline BOOL GetBaseBoardByCmd(char* lpszBaseBoard, int len/*=128*/)
	{
		const long MAX_COMMAND_SIZE = 10000; // ��������������С	
		WCHAR szFetCmd[] = L"wmic BaseBoard get SerialNumber"; // ��ȡ�������к�������	
		const std::string strEnSearch = "SerialNumber"; // �������кŵ�ǰ����Ϣ

		BOOL   bret = FALSE;
		HANDLE hReadPipe = NULL; //��ȡ�ܵ�
		HANDLE hWritePipe = NULL; //д��ܵ�	
		PROCESS_INFORMATION pi;   //������Ϣ	
		STARTUPINFO			si;	  //���������д�����Ϣ
		SECURITY_ATTRIBUTES sa;   //��ȫ����

		char			szBuffer[MAX_COMMAND_SIZE + 1] = { 0 }; // ���������н�������������
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

		//1.0 �����ܵ�
		bret = CreatePipe(&hReadPipe, &hWritePipe, &sa, 0);
		if (!bret)
		{
			//goto END;
			return false;
		}

		//2.0 ���������д��ڵ���ϢΪָ���Ķ�д�ܵ�
		GetStartupInfo(&si);
		si.hStdError = hWritePipe;
		si.hStdOutput = hWritePipe;
		si.wShowWindow = SW_HIDE; //���������д���
		si.dwFlags = STARTF_USESHOWWINDOW | STARTF_USESTDHANDLES;

		//3.0 ������ȡ�����еĽ���
		bret = CreateProcess(NULL, szFetCmd, NULL, NULL, TRUE, 0, NULL, NULL, &si, &pi);
		if (!bret)
		{
			//goto END;
			return false;
		}

		//4.0 ��ȡ���ص�����
		WaitForSingleObject(pi.hProcess, 500/*INFINITE*/);
		bret = ReadFile(hReadPipe, szBuffer, MAX_COMMAND_SIZE, &count, 0);
		if (!bret)
		{
			//goto END;
			return false;
		}

		//5.0 �����������к�
		bret = FALSE;
		strBuffer = szBuffer;
		ipos = long(strBuffer.find(strEnSearch));

		if (ipos < 0) // û���ҵ�
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

		//ȥ���м�Ŀո� \r \n \0 \t
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

		//�ر����еľ��
		CloseHandle(hWritePipe);
		CloseHandle(hReadPipe);
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);

		return(bret);
	}

	inline void getcpuidex(unsigned int* CPUInfo, unsigned int InfoType, unsigned int ECXValue)
	{
#if defined(_MSC_VER) // MSVC  
#if defined(_WIN64) // 64λ�²�֧���������. 1600: VS2010, ��˵VC2008 SP1֮���֧��__cpuidex.  
		__cpuidex((int*)(void*)CPUInfo, (int)InfoType, (int)ECXValue);
#else  
		if (NULL == CPUInfo)
			return;
		_asm {
			// load. ��ȡ�������Ĵ���.  
			mov edi, CPUInfo;
			mov eax, InfoType;
			mov ecx, ECXValue;
			// CPUID  
			cpuid;
			// save. ���Ĵ������浽CPUInfo  
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

		// ��ȡCPU ID
#ifdef _WIN32
		__cpuid(dwBuf, 1); // Windowsƽ̨ʹ��__cpuid
#else
		__cpuid(1, dwBuf[0], dwBuf[1], dwBuf[2], dwBuf[3]); // Linuxƽ̨ʹ��__cpuid
#endif

		// ��ʽ��CPU ID
		char cpuId[17] = { 0 }; // 16�ַ� + 1��ֹ��
#ifdef _WIN32
		sprintf_s(cpuId, sizeof(cpuId), "%08X", dwBuf[3]); // Windowsƽ̨ʹ��sprintf_s
		sprintf_s(cpuId + 8, sizeof(cpuId) - 8, "%08X", dwBuf[0]);
#else
		snprintf(cpuId, sizeof(cpuId), "%08X", dwBuf[3]); // Linuxƽ̨ʹ��snprintf
		snprintf(cpuId + 8, sizeof(cpuId) - 8, "%08X", dwBuf[0]);
#endif

		return std::string(cpuId); // ����std::string
	}

}