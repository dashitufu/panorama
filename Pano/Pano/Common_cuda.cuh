#include "cuda_runtime.h"
#include <iostream>
#include "stdio.h"

extern "C"
{
#include "Buddy_System.h"
}
template<typename _T, int iSize>struct Data_Block {
    _T Data[iSize];
};

extern Mem_Mgr oMem_Mgr_GPU;

void Init_Env_GPU();
void Free_Env_GPU();
void Init_Env_All();
void Free_Env_All();
void* pMalloc_GPU(unsigned int iSize);
void Free_GPU(void* p);
void Disp_Mem_GPU();
void Shrink_GPU(void* p, unsigned int iSize);
void Disp_Cuda_Error();
template<typename _T>void Disp_Part_GPU(_T* M, int iStride, int x, int y, int w, int h, const char* pcCaption = NULL);
template<typename _T>void Disp_GPU(_T* M, int iHeight, int iWidth, const char* pcCaption = NULL);
template<typename _T>void Disp_Sum_GPU(_T M[], int iSize);
int bLoad_Raw_Data_GPU(const char* pcFile, unsigned char** ppBuffer, int* piSize);
int bSave_Raw_Data_GPU(const char* pcFile, unsigned char* pBuffer, int iSize);
