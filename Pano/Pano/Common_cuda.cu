#include "Common.h"
#include "Common_cuda.cuh"
#include "cuda_runtime.h"

extern "C"
{
#include "Buddy_System.h"
}
using namespace std;
Mem_Mgr oMem_Mgr_GPU{};

void SB_Common_Cuda()
{//所有的模板函数实例化
	Disp_Part_GPU<char>(NULL, 0, 0, 0, 0, 0);
    Disp_Part_GPU<unsigned char>(NULL, 0, 0, 0, 0, 0);
	Disp_Part_GPU<int>(NULL, 0, 0, 0, 0, 0);
	Disp_Part_GPU<unsigned int>(NULL, 0, 0, 0, 0, 0);
	Disp_Part_GPU<float>(NULL, 0, 0, 0, 0, 0);
	Disp_Part_GPU<double>(NULL, 0, 0, 0, 0, 0);
	Disp_Part_GPU<short>(NULL, 0, 0, 0, 0, 0);
	Disp_Part_GPU<unsigned short>(NULL, 0, 0, 0, 0, 0);

	Disp_GPU<char>(NULL, 0, 0);
	Disp_GPU<unsigned char>(NULL, 0, 0);
	Disp_GPU<int>(NULL, 0, 0);
	Disp_GPU<unsigned int>(NULL, 0, 0);
	Disp_GPU<float>(NULL, 0, 0);
	Disp_GPU<double>(NULL, 0, 0);
	Disp_GPU<short>(NULL, 0, 0);
	Disp_GPU<unsigned short > (NULL, 0, 0);

	Disp_Sum_GPU<double> (NULL, 0);
	Disp_Sum_GPU<float>(NULL, 0);
	Disp_Sum_GPU<int>(NULL, 0);
	Disp_Sum_GPU<unsigned int>(NULL, 0);
	Disp_Sum_GPU<char>(NULL, 0);
	Disp_Sum_GPU<unsigned char>(NULL, 0);
	Disp_Sum_GPU<short>(NULL, 0);
	Disp_Sum_GPU<unsigned short>(NULL, 0);
}

/*****************一组显存管理函数**********************/
void Init_Env()
{//初始化个环境，分配些内存供一切函数临时使用
	//这个项目用统一内存
	unsigned long long iSize = 1000000000;
	const int iBlock_Size = 2048;
	void* pBuffer = NULL;
	cudaError cudaStatus=cudaMallocHost(&pBuffer, iSize);
	if (cudaStatus != cudaSuccess)
	{
		printf("Fail to allocate memory in Init_Env");
		return;
	}
	Init_Mem_Mgr_GPU(&oMem_Mgr, iSize, iBlock_Size, 997, pBuffer);
	return;
}

void Init_Env_GPU()
{//初始化个环境，分配些内存供一切函数临时使用
	//初始化GPU的内存管理
	unsigned long long iSize = 2000000000;
	const int iBlock_Size = 2048;
	void* pBuffer = NULL;
	cudaError cudaStatus;
	cudaStatus = cudaMalloc(&pBuffer, iSize + iBlock_Size);
	if (cudaStatus != cudaSuccess)
	{
		printf("Fail to allocate memory in Init_Env");
		return;
	}
	Init_Mem_Mgr_GPU(&oMem_Mgr_GPU, iSize, iBlock_Size, 997, pBuffer);
	return;
}

void Free_Env()
{
	printf("CPU Memory\n");
	if (oMem_Mgr.m_iPiece_Count)
		Disp_Mem(&oMem_Mgr, 0);

	if (oMem_Mgr.m_pBuffer)
	{
		cudaFree(oMem_Mgr.m_pOrg_Buffer);
		oMem_Mgr.m_pOrg_Buffer = NULL;
		Free_Mem_Mgr(&oMem_Mgr);
	}
	oMem_Mgr = {};
}

void Free_Env_GPU()
{
	printf("GPU Memory\n");
	//释放咸淳管理器
	if (oMem_Mgr_GPU.m_iPiece_Count)
		Disp_Mem(&oMem_Mgr_GPU, 0);

	if (oMem_Mgr_GPU.m_pBuffer)
	{
		cudaFree(oMem_Mgr_GPU.m_pOrg_Buffer);
		oMem_Mgr_GPU.m_pOrg_Buffer = NULL;
		Free_Mem_Mgr(&oMem_Mgr_GPU);
	}
	oMem_Mgr_GPU = {};
}
void* pMalloc_GPU(unsigned int iSize)
{//原来的pMalloc多了个参数太麻烦，干脆包一层更少
	return pMalloc(&oMem_Mgr_GPU, iSize);
}
void Free_GPU(void* p)
{
	Free(&oMem_Mgr_GPU, p);
}

void Init_Env_All()
{
	Init_Env();
	Init_Env_GPU();
}
void Free_Env_All()
{
	Free_Env();
	Free_Env_GPU();
}
void Disp_Mem_GPU()
{
	Disp_Mem(&oMem_Mgr_GPU, 0);
}
void Shrink_GPU(void* p, unsigned int iSize)
{
	Shrink(&oMem_Mgr_GPU, p, iSize);
	return;
}
/*****************一组显存管理函数**********************/
void Disp_Cuda_Error()
{
	cudaDeviceSynchronize();
	cudaError_t cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess)
		printf("%s\n", cudaGetErrorString(cudaStatus));
}

template<typename _T>__global__ void _Disp_Part_GPU(_T* M, int iStride, int x, int y, int w, int h)
{
	int i, j;
	for (i = 0; i < h; i++)
	{
		for (j = 0; j < w; j++)
		{
			int iPos = (y + i) * iStride + (x + j);
			if (std::is_same_v<_T, float> || 
				std::is_same_v<_T, double>)
			{
				if(M[iPos] - (int)(M[iPos]))
					printf("%.8f\t", (float)M[iPos]);
				else
					printf("%d ", (int)M[iPos]);
			}
			else if (std::is_same_v<_T, int> || std::is_same_v<_T, unsigned int> ||
				std::is_same_v<_T, short> || std::is_same_v<_T, unsigned short> ||
				std::is_same_v<_T, char> || std::is_same_v<_T, unsigned char>)
				printf("%d,", (int)M[iPos]);

			//printf("Pos:%d\n", iPos);
		}
		printf("\n");
	}
}

template<typename _T>void Disp_Part_GPU(_T* M, int iStride, int x, int y, int w, int h, const char* pcCaption)
{
	if (w + x > iStride)
	{
		printf("x position exceeds width\n");
		return;
	}
	if (pcCaption)
		printf("%s\n", (char*)pcCaption);
	Disp_Cuda_Error();
	_Disp_Part_GPU << <1, 1 >> > (M, iStride, x, y, w, h);
	Disp_Cuda_Error();
}

template<typename _T>__global__ void _Disp_GPU(_T* M, int iHeight, int iWidth)
{
	for (int i = 0; i < iHeight; i++)
	{
		for (int j = 0; j < iWidth; j++)
		{
			if (std::is_same_v<_T, float>)
			{
				if (M[i * iWidth + j] - (int)M[i * iWidth + j])
					//printf("%.10ef, ", (double)M[i * iWidth + j]);
					printf("%f, ", (double)M[i * iWidth + j]);
				else
					printf("%d, ", (int)M[i * iWidth + j]);
			}
			else if (std::is_same_v<_T, double>)
			{
				if (M[i * iWidth + j] - (int)M[i * iWidth + j])
					printf("%f, ", (double)M[i * iWidth + j]);
				//printf("%.10ef, ", (double)M[i * iWidth + j]);
				else
					printf("%d, ", (int)M[i * iWidth + j]);
				//printf("%f,", (double)M[i * iWidth + j]);
			}
			else if ( std::is_same_v<_T, unsigned int> ||
				std::is_same_v<_T, int> ||
				std::is_same_v<_T, short> ||
				std::is_same_v<_T, unsigned short> ||
				std::is_same_v<_T, char> ||
				std::is_same_v<_T, unsigned char>)
				printf("%d   ", (int)(M[i * iWidth + j]));
		}
		printf("\n");
	}
	return;
}

template<typename _T>void Disp_GPU(_T* M, int iHeight, int iWidth, const char* pcCaption)
{
	if (pcCaption)
		printf("%s\n", (char*)pcCaption);
	Disp_Cuda_Error();
	_Disp_GPU << <1, 1 >> > (M, iHeight,iWidth);
	Disp_Cuda_Error();
}
template<typename _T>__global__ void _Disp_Sum_GPU(_T M[], int iSize)
{
	double fTotal = 0;
	for (int i = 0; i < iSize; i++)
		fTotal += M[i];
	printf("Sum: %lf\n", fTotal);
}
template<typename _T>void Disp_Sum_GPU(_T M[], int iSize)
{
	_Disp_Sum_GPU << <1, 1 >> > (M, iSize);
	Disp_Cuda_Error();
}
int bSave_Raw_Data_GPU(const char* pcFile, unsigned char* pBuffer, int iSize)
{
	unsigned char* pBuffer_1 = (unsigned char*)pMalloc(iSize);
	cudaMemcpy(pBuffer_1, pBuffer, iSize, cudaMemcpyDeviceToHost);
	Disp_Cuda_Error();
	int bResult=bSave_Raw_Data(pcFile, pBuffer_1, iSize);
	Free(pBuffer_1);
	return bResult;
}

int bLoad_Raw_Data_GPU(const char* pcFile, unsigned char**ppBuffer, int* piSize)
{
	unsigned char* pBuffer = NULL, * pBuffer_GPU = *ppBuffer;
	if (!bLoad_Raw_Data(pcFile, &pBuffer, piSize))
		return 0;
	
	if (!pBuffer_GPU)
		pBuffer_GPU= (unsigned char*)pMalloc(*piSize);
	cudaMemcpy(pBuffer_GPU, pBuffer, *piSize, cudaMemcpyHostToDevice);
	Disp_Cuda_Error();
	return 0;
}