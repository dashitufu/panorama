//#undef WIN32

#include "Common.h"
#ifdef WIN32
#include "io.h"
#include "windows.h"
#endif
#ifndef WIN32
#include "sys/time.h"
#endif

extern "C"
{
#include "Buddy_System.h"
}

Mem_Mgr oMem_Mgr = {};

#ifdef WIN32

void Init_BitPtr(BitPtr* poBitPtr, unsigned char* pBuffer, int iSize)
{
	poBitPtr->m_iBitPtr = poBitPtr->m_iCur = 0;
	poBitPtr->m_iEnd = iSize;
	poBitPtr->m_pBuffer = pBuffer;
}

int iGetBits(BitPtr* poBitPtr, int iLen)
{
	unsigned int iValue, iCur = poBitPtr->m_iCur, iShiftBits;
	int iBitLeft = iLen;
	//先搞定第一字节
	iValue = ((poBitPtr->m_pBuffer[poBitPtr->m_iCur] << poBitPtr->m_iBitPtr) & 0xFF) << 24;
	iBitLeft -= (8 - poBitPtr->m_iBitPtr);
	iCur++;
	iShiftBits = 16 + poBitPtr->m_iBitPtr;
	while (iBitLeft > 0)
	{
		if (iShiftBits > 0)
			iValue |= (poBitPtr->m_pBuffer[iCur] & 0xFF) << iShiftBits;
		iBitLeft -= 8;
		iShiftBits -= 8;
		iCur++;
	}
	iValue = iLen ? iValue >> (32 - iLen) : 0;
	poBitPtr->m_iCur += (poBitPtr->m_iBitPtr + iLen) >> 3;
	poBitPtr->m_iBitPtr = (iLen + poBitPtr->m_iBitPtr) & 7;
	return iValue;
}

void Get_All_File(const char* pcPath,char *pBuffer)
{//取一个目录所有文件 例: c:\\tmp\\*.bin
	intptr_t handle;
	_finddata_t findData;
	int iCount = 0;
	char* pCur = pBuffer;
	handle = _findfirst(pcPath, &findData);    // 查找目录中的第一个文件
	if (handle == -1)
	{
		printf("File not found\n");
		return;
	}
	do
	{
		if (findData.attrib & _A_ARCH)
		{
			strcpy(pCur, findData.name);
			pCur += strlen(pCur)+1;
			iCount++;
		}	
	} while (_findnext(handle, &findData) == 0);    // 查找目录中的下一个文件
	_findclose(handle);    // 关闭搜索句柄
}
int iGet_File_Count(const char* pcPath)
{//给定文件路径，求此路径下得文件数量
	intptr_t handle;
	_finddata_t findData;
	int iCount = 0;
	handle = _findfirst(pcPath, &findData);    // 查找目录中的第一个文件
	if (handle == -1)
	{
		printf("File not found\n");
		return 0;
	}
	do
	{
		if (findData.attrib & _A_ARCH)
			iCount++;
	} while (_findnext(handle, &findData) == 0);    // 查找目录中的下一个文件
	_findclose(handle);    // 关闭搜索句柄
	return iCount;
}
#endif

int bSave_Bin(const char* pcFile, float* pData, int iSize)
{
	FILE* pFile = fopen(pcFile, "wb");
	if (!pFile)
	{
		printf("Fail to save:%s\n", pcFile);
		return 0;
	}
	int iResult = (int)fwrite(pData, 1, iSize, pFile);
	if (iResult != iSize)
	{
		printf("Fail to save:%s\n", pcFile);
		iResult = 0;
	}
	else
		iResult = 1;
	return iResult;
}
int iGet_Upper_Triangle_Size(int w)
{
	int iSize;
	iSize = w * (w - 1) / 2;
	return iSize;
}

int iUpper_Triangle_Cord_2_Index(int x, int y, int w)
{//一个上三角矩阵，给定(x,y)坐标，转换为索引值, 注意，w为上三角矩形的列数, 不是上三角第一行有效元数个数
	int iPos;
	//比如 w=4, y=2 x=3 Index=7
	if (x <= y || y >= w - 1)
		return -1;
	iPos = ((w - 1) + (w - y)) * y / 2 + x - y - 1;
	return iPos;
}

unsigned long long iGet_File_Length(char* pcFile)
{//return: >-0 if success; -1 if fail
	FILE* pFile = fopen(pcFile, "rb");
	long long iLen;
#ifdef WIN32
	if (!pFile)
	{

		int iResult = GetLastError();
		printf("Fail to get file length, error:%d\n", iResult);

		return -1;
	}		
	_fseeki64(pFile, 0, SEEK_END);
	iLen = _ftelli64(pFile);
#else
	if (!pFile)
	{
		printf("Fail to get file length\n");
		return -1;
	}

	fseeko64(pFile, 0, SEEK_END);
	iLen = ftello64(pFile);
#endif

	fclose(pFile);
	return iLen;
}

unsigned long long iGet_Tick_Count()
{//求当前的毫秒级Tick Count。之所以不搞微秒级是因为存在时间片问题，即使毫秒级也不准，很有可能Round到16毫秒一个时间片，聊胜于无
#ifndef WIN32
	timeval oTime;
	gettimeofday(&oTime,NULL);
	unsigned long long tTime = oTime.tv_sec * 1000000 + oTime.tv_usec;
	return tTime/1000;
#else
	 timeb tp;
	 ftime(&tp);
	 return (unsigned long long) ((unsigned long long)tp.time * 1000 + tp.millitm);
#endif

}
unsigned int iGet_Random_No()
{//伪随机	(789, 123)	(389,621) 已通过，可以自定义自己的种子
#define b 389
#define c 621
	static unsigned int iNum = 0xFFFFFFF;	//GetTickCount();	//种子
	return iNum = iNum * b + c;
#undef b
#undef c
}
template<typename _T>void Get_Random_Norm_Vec(_T V[], int n)
{//生成一个归一化向量
	double fTotal=0;
	int i;
	for (i = 0; i < n; i++)
	{
		unsigned int iValue =  iGet_Random_No();
		V[i] = (_T)iValue;
		fTotal += (unsigned long long)iValue * iValue;
	}
	_T fMod = (_T)sqrt((double)fTotal);
	for (i = 0; i < n; i++)
		V[i] /= fMod;
}

int iRandom(int iStart, int iEnd)
{//从iStart到iEnd之间随机出一个数字，由于RandomInteger太傻逼，没有必要把时间浪费在傻逼身上
	return iStart + iGet_Random_No() % (iEnd - iStart + 1);
}
int iRandom()
{//尝试搞一个无参数获得随机数的函数，计算办法为 iRandom_No= iRandom_No*a + c; 最后取模 0x7FFFFFFF
#define m 1999999973			//基于模为素数能令散列更均匀的理论
	static unsigned int a = (unsigned int)iGet_Tick_Count(); //1103515245;		//1103515245为很好的初始值，但选取静态值变成了确定问题
	static unsigned int c = 2347;	// 12345;
	static unsigned long long iRandom_No = 1;
	iRandom_No = (iRandom_No * a + c) % m;	//求和不溢出，再求模
	return (int)iRandom_No;
#undef m
}
template<typename _T>int iQuick_Sort_Partition(_T pBuffer[], int left, int right)
{//小到大的顺序
 //Real fRef;
	_T iValue, oTemp;
	int pos = right;

	//试一下将中间元素交换到头，多一步可能挽救了极端差形态
	oTemp = pBuffer[right];
	pBuffer[right] = pBuffer[(left + right) >> 1];
	pBuffer[(left + right) >> 1] = oTemp;

	right--;
	iValue = pBuffer[pos];
	while (left <= right)
	{
		while (left < pos && pBuffer[left] <= iValue)
			left++;
		while (right >= 0 && pBuffer[right] > iValue)
			right--;
		if (left >= right)
			break;
		oTemp = pBuffer[left];
		pBuffer[left] = pBuffer[right];
		pBuffer[right] = oTemp;
	}

	oTemp = pBuffer[left];
	pBuffer[left] = pBuffer[pos];
	pBuffer[pos] = oTemp;

	return left;
}
template<typename _T>int iAdjust_Left(_T* pStart, _T* pEnd)
{
	_T oTemp, * pCur_Left = pEnd - 1,
		* pCur_Right;
	_T oRef = *pEnd;

	//为了减少一次判断，此处先扫过去
	while (pCur_Left >= pStart && *pCur_Left == oRef)
		pCur_Left--;
	pCur_Right = pCur_Left;
	pCur_Left--;

	while (pCur_Left >= pStart)
	{
		if (*pCur_Left == oRef)
		{
			oTemp = *pCur_Left;
			*pCur_Left = *pCur_Right;
			*pCur_Right = oTemp;
			pCur_Right--;
		}
		pCur_Left--;
	}
	return (int)(pCur_Right - pStart);
}

template<typename _T> _T oGet_Nth_Elem(_T Seq[],int iCount, int iStart, int iEnd, int iNth,int bAdjust_Pos=0)
{//取第n大元素，用Quick_Sort的Partition做
	int iPos;
	if (iStart < iEnd)
	{
		iPos = iQuick_Sort_Partition(Seq, iStart, iEnd);
		if (iPos > iNth)//那么第n大元素在左分区内
			return oGet_Nth_Elem(Seq,iCount, iStart, iPos - 1, iNth);
		else if (iPos < iNth)//第n大元素在右分区内
			return oGet_Nth_Elem(Seq,iCount, iPos + 1, iEnd, iNth);
		else //找到了，正好在iPos中
			return Seq[iPos];
	}else
	{//此时又分奇偶两种情况
		if(iCount&1 || !bAdjust_Pos)
			return Seq[iStart];	//奇数好办，返回便是
		else //偶数的还要往前找最大值
		{
			_T fMax = Seq[iStart - 1];
			for (int i = iStart - 2; i >= 0; i--)
				if (Seq[i] > fMax)
					fMax = Seq[i];
			return (_T)((Seq[iStart] + fMax) / 2.f);
		}
	}	
}

template<typename _T> _T oGet_Nth_Elem(_T Seq[], int iCount, int iNth)
{
	return oGet_Nth_Elem(Seq, iCount, 0, iCount - 1, iNth);
}
template<typename _T> void Quick_Sort(_T Seq[], int iStart, int iEnd)
{
	int iPos, iLeft;
	if (iStart < iEnd)
	{
		iPos = iQuick_Sort_Partition(Seq, iStart, iEnd);
		//如果iPos>=iEnd, 则遇上极差形态的组了，这时做一个推进，找出所有与Seq[iEnd]相等的项目，赶到右边
		if (iPos >= iEnd)
		{
			iLeft = iAdjust_Left(Seq + iStart, Seq + iPos);
			iLeft = iStart + iLeft;
		}
		else
			iLeft = iPos - 1;

		if (iStart < iLeft)
			Quick_Sort(Seq, iStart, iLeft);

		Quick_Sort(Seq, iPos + 1, iEnd);
	}
}
void SB_Common()
{//template实例化，只对vc有效
	Temp_Load_File_2(NULL, NULL, NULL, (Point_2D<double>**)NULL,(double(**)[3])NULL,(double(**)[9])NULL);
	Temp_Load_File_2(NULL, NULL, NULL, (Point_2D<float>**)NULL, (float(**)[3])NULL, (float(**)[9])NULL);

	bSave_PLY(NULL, (double(*)[3])NULL, 0);
	bSave_PLY(NULL, (float(*)[3])NULL, 0);

	oGet_Nth_Elem((double*)NULL, 0, 0);
	oGet_Nth_Elem((float*)NULL, 0, 0);
	oGet_Nth_Elem((int*)NULL, 0, 0);

	Quick_Sort((double*)NULL, 0, 0);
	Quick_Sort((float*)NULL, 0, 0);
	Quick_Sort((int*)NULL, 0, 0);

	Get_Random_Norm_Vec((double*)NULL, 0);
	Get_Random_Norm_Vec((float*)NULL, 0);

	Init_Point_Cloud((Point_Cloud<float>*)NULL, 0, 0);
	Init_Point_Cloud((Point_Cloud<double>*)NULL, 0, 0);

	Free_Point_Cloud((Point_Cloud<float>*)NULL);
	Free_Point_Cloud((Point_Cloud<double>*)NULL);

	bSave_PLY(NULL, Point_Cloud < double>{});
	bSave_PLY(NULL, Point_Cloud < float>{});

	Draw_Line((Point_Cloud<double>*)NULL, (double)0, (double)0, (double)0, (double)0, (double)0, (double)0);
	Draw_Line((Point_Cloud<float>*)NULL, (float)0, (float)0, (float)0, (float)0, (float)0, (float)0);

	Draw_Sphere((Point_Cloud<double>*)NULL, (double)0, (double)0, (double)0);
	Draw_Sphere((Point_Cloud<float>*)NULL, (float)0, (float)0, (float)0);

	Draw_Rect<float>(NULL, 0, 0, 0, 0, 0);
}
int bLoad_Text_File(const char* pcFile, char** ppBuffer, int* piSize)
{
	FILE* pFile = fopen(pcFile, "rb");
	int bRet = 0, iResult, iSize;
	char* pBuffer;
	iSize = (int)iGet_File_Length((char*)pcFile);
	pBuffer = (char*)malloc(iSize+1);

	if (!pFile)
	{
		printf("Fail to open file:%s\n", pcFile);
		goto END;
	}
	if (!pBuffer)
	{
		printf("Fail to allocate memory\n");
		goto END;
	}

	iResult = (int)fread(pBuffer, 1, iSize, pFile);
	if (iResult != iSize)
	{
		if (pBuffer)
			free(pBuffer);
		*ppBuffer = NULL;
		printf("Fail to read data\n");
		goto END;
	}
	pBuffer[iSize] = '\0';
	*ppBuffer = pBuffer;
	if (piSize)
		*piSize = iSize;
	bRet = 1;
END:
	if (pFile)
		fclose(pFile);
	if (!bRet)
	{
		if (pBuffer)
			free(pBuffer);
	}
	return bRet;
}

int bSave_Raw_Data(const char* pcFile, unsigned char* pBuffer, int iSize)
{
	FILE* pFile = fopen(pcFile, "wb");
	if (!pFile)
	{
		printf("Fail to save file:%s\n", pcFile);
		return 0;
	}
	int iResult=(int)fwrite(pBuffer, 1, iSize, pFile);
	if (iResult != iSize)
	{
		//printf("Fail to save file:%s %d\n", pcFile,GetLastError());
		fclose(pFile);
		return 0;
	}
	fclose(pFile);
	return 1;
}
int bLoad_Raw_Data(const char* pcFile, unsigned char** ppBuffer, int iSize, int bNeed_Malloc, int iFrame_No)
{
	FILE* pFile = fopen(pcFile, "rb");
	unsigned long long iPos;
	int bRet = 0, iResult;
	iPos = (unsigned long long)iSize * iFrame_No;
	unsigned char* pBuffer;

	if (!iSize)
		iSize = (int)iGet_File_Length((char*)pcFile);

	if (bNeed_Malloc)
		pBuffer = (unsigned char*)malloc(iSize);
	else
		pBuffer = *ppBuffer;

	if (!pFile)
	{
		printf("Fail to open file:%s\n", pcFile);
		goto END;
	}
	if (!pBuffer)
	{
		printf("Fail to allocate memory\n");
		goto END;
	}
#ifdef WIN32
	_fseeki64(pFile, iPos, SEEK_SET);
#else
	fseeko64(pFile, (unsigned long long)iPos, SEEK_SET);
#endif

	iResult = (int)fread(pBuffer, 1, iSize, pFile);
	if (iResult != iSize)
	{
		if (pBuffer)
			free(pBuffer);
		*ppBuffer = NULL;
		printf("Fail to read data\n");
		bRet = 0;
		goto END;
	}
	*ppBuffer = pBuffer;
	bRet = 1;
END:
	if (pFile)
		fclose(pFile);
	if (!bRet)
	{
		if (pBuffer && bNeed_Malloc)
			free(pBuffer);
	}
	return bRet;
}
int bLoad_Raw_Data(const char* pcFile, unsigned char** ppBuffer,int *piSize)
{
	FILE* pFile = fopen(pcFile, "rb");
	int bRet = 0, iResult,iSize;
	unsigned char* pBuffer;
	iSize = (int)iGet_File_Length((char*)pcFile);
	pBuffer = (unsigned char*)malloc(iSize);

	if (!pFile)
	{
		printf("Fail to open file:%s\n", pcFile);
		goto END;
	}
	if (!pBuffer)
	{
		printf("Fail to allocate memory\n");
		goto END;
	}

	iResult = (int)fread(pBuffer, 1, iSize, pFile);
	if (iResult != iSize)
	{
		if (pBuffer)
			free(pBuffer);
		*ppBuffer = NULL;
		printf("Fail to read data\n");
		goto END;
	}
	*ppBuffer = pBuffer;
	if (piSize)
		*piSize = iSize;
	bRet = 1;
END:
	if (pFile)
		fclose(pFile);
	if (!bRet)
	{
		if (pBuffer)
			free(pBuffer);
	}
	return bRet;
}

template<typename _T>int bSave_PLY(const char* pcFile, _T Point[][3], int iPoint_Count, unsigned char Color[][3], int bText)
{//存点云，最简形式，用于实验，连结构都不要
	FILE* pFile = fopen(pcFile, "wb");
	char Header[512];
	int i;
	_T* pPos;

	if (!pFile)
	{
		printf("Fail to open file:%s\n", pcFile);
		return 0;
	}
	if (!iPoint_Count)
	{
		printf("No point to save\n");
		return 0;
	}

	//先写入Header
	sprintf(Header, "ply\r\n");
	if (bText)
		sprintf(Header + strlen(Header), "format ascii 1.0\r\n");
	else
		sprintf(Header + strlen(Header), "format binary_little_endian 1.0\r\n");
	sprintf(Header + strlen(Header), "comment HQYT generated\r\n");
	sprintf(Header + strlen(Header), "element vertex %d\r\n", iPoint_Count);
	sprintf(Header + strlen(Header), "property float x\r\n");
	sprintf(Header + strlen(Header), "property float y\r\n");
	sprintf(Header + strlen(Header), "property float z\r\n");

	if (Color)
	{
		sprintf(Header + strlen(Header), "property uchar red\r\n");
		sprintf(Header + strlen(Header), "property uchar green\r\n");
		sprintf(Header + strlen(Header), "property uchar blue\r\n");
	}

	sprintf(Header + strlen(Header), "end_header\r\n");
	fwrite(Header, 1, strlen(Header), pFile);

	for (i = 0; i < iPoint_Count; i++)
	{
		pPos = Point[i];
		if (bText)
		{
			fprintf(pFile, "%f %f %f ", pPos[0], pPos[1], pPos[2]);
			if (pPos[0] >= 1000000.f || isnan(pPos[0]))
				printf("here");
			if (Color)
				fprintf(pFile, "%d %d %d\r\n", Color[i][0], Color[i][1], Color[i][2]);
			else
				fprintf(pFile, "\r\n");
		}
	}
	fclose(pFile);
	return 1;
}

template<typename _T>void Temp_Load_File_2(int* piCameta_Count, int* piPoint_Count, int* piObservation_Count, Point_2D<_T>** ppPoint_2D, _T(**ppPoint_3D)[3], _T(**ppCamera)[3 * 3])
{//这次装入problem-16-22106-pre.txt。搞个好点的数据结构指出匹配关系

	int i, iCamera_Count, iPoint_Count, iObservation_Count;
	_T(*pPoint_3D)[3], (*pCamera)[3 * 3];
	Point_2D<_T>* pPoint_2D;
	char* pcFile = (char*)"c:\\tmp\\temp\\problem-16-22106-pre.txt";	//"Sample\\problem-16-22106-pre.txt";
	FILE* pFile = fopen(pcFile, "rb");
	i = fscanf(pFile, "%d %d %d\n", &iCamera_Count, &iPoint_Count, &iObservation_Count);
	pPoint_2D = (Point_2D<_T>*)malloc(iObservation_Count * 2 * sizeof(Point_2D<_T>));
	pPoint_3D = (_T(*)[3])malloc(iPoint_Count * 3 * sizeof(_T));
	pCamera = (_T(*)[3 * 3])malloc(iCamera_Count * 16 * sizeof(_T));

	for (i = 0; i < iObservation_Count; i++)
	{
		float xy[2];
		fscanf(pFile, "%d %d %f %f", &pPoint_2D[i].m_iCamera_Index, &pPoint_2D[i].m_iPoint_Index, &xy[0], &xy[1]);
		pPoint_2D[i].m_Pos[0] = xy[0], pPoint_2D[i].m_Pos[1] = xy[1];
		//fscanf(pFile, "%d %d %f %f", &pPoint_2D[i].m_iCamera_Index, &pPoint_2D[i].m_iPoint_Index, &pPoint_2D[i].m_Pos[0], &pPoint_2D[i].m_Pos[1]);
	}

	int iResult;
	for (i = 0; i < iCamera_Count; i++)
		for (int j = 0; j < 9; j++)
			if(sizeof(_T)==4)
				iResult=fscanf(pFile, "%f", (float*)&pCamera[i][j]);
			else
				iResult=fscanf(pFile, "%lf", (double*)&pCamera[i][j]);
	for (i = 0; i < iPoint_Count; i++)
		if(sizeof(_T)==4)
			fscanf(pFile, "%f %f %f ",(float*)&pPoint_3D[i][0], (float*)&pPoint_3D[i][1], (float*)&pPoint_3D[i][2]);
		else
			fscanf(pFile, "%lf %lf %lf ", (double*)&pPoint_3D[i][0], (double*)&pPoint_3D[i][1], (double*)&pPoint_3D[i][2]);
	fclose(pFile);
	*piCameta_Count = iCamera_Count;
	*piPoint_Count = iPoint_Count;
	*piObservation_Count = iObservation_Count;
	*ppPoint_2D = pPoint_2D;
	*ppPoint_3D = pPoint_3D;
	*ppCamera = pCamera;
	//Disp((_T*)pCamera, 16, 9, "Camera");
	return;
}

void* pMalloc(Light_Ptr* poPtr, int iSize)
{//轻量级的内存分配
	unsigned char* pBuffer;
	Malloc(*poPtr, iSize, pBuffer);
	return (void*)pBuffer;
}

void *pMalloc(unsigned int iSize)
{//原来的pMalloc多了个参数太麻烦，干脆包一层更少
	return pMalloc(&oMem_Mgr, iSize);
}
void Free(void* p)
{
	Free(&oMem_Mgr, p);
}
void Shrink(void* p, unsigned int iSize)
{
	Shrink(&oMem_Mgr, p, iSize);
	return;
}
void Disp_Mem()
{
	Disp_Mem(&oMem_Mgr, 0);
}

void Init_Env_CPU()
{//初始化个环境，分配些内存供一切函数临时使用
	//这个项目用统一内存
	unsigned long long iSize = 1000000000;
	//void* p;
	Init_Mem_Mgr(&oMem_Mgr, iSize, 2048, 997);
	return;
}

void Free_Env_CPU()
{
	if (oMem_Mgr.m_iPiece_Count)
		Disp_Mem(&oMem_Mgr, 0);

	if (oMem_Mgr.m_pBuffer)
		Free_Mem_Mgr(&oMem_Mgr);
	oMem_Mgr = {};
}

void Cal_Line(Line_1* poLine, float x0, float y0, float x1, float y1)
{//过两点决定直线方程
	poLine->x0 = x0, poLine->y0 = y0, poLine->x1 = x1, poLine->y1 = y1;
	//求ax+bx+c=0
	poLine->a = (float)(y0 - y1);
	poLine->b = (float)(x1 - x0);
	poLine->c = (float)(x0 * y1 - x1 * y0);

	if (x0 != x1)
	{//此时才有y=kx+m的表示
		poLine->k = (float)(y0 - y1) / (x0 - x1);
		poLine->m = y0 - poLine->k * x0;
		poLine->theta = (float)atan(poLine->k);
		//还要决定直线的方向
		if (x1 < x0)
			poLine->theta += (float)PI;

		if (y0 == y1)
			poLine->m_iMode = Line_1::Mode::Degree_180;
		else
			poLine->m_iMode = Line_1::Mode::Degree_Other;
	}
	else
	{
		poLine->k = NAN;
		poLine->m = x0;	//x=x0
		if (y0 > y1)
			poLine->theta =(float)( 3 * PI / 2);
		else
			poLine->theta = (float)(PI / 2);
		poLine->m_iMode = Line_1::Mode::Degree_90;
	}
	return;
}
float fGet_Line_y(Line_1* poLine, float x)
{//代入x，求y
	if (poLine->m_iMode == Line_1::Mode::Degree_90)
	{
		if (x != poLine->m)
			return NAN;	//此处不适用于cuda，得另想办法
		else
			return 0;	//随便返回一个值即可
	}
	else
		return poLine->k * x + poLine->m;
}
float fGet_Line_x(Line_1* poLine, float y)
{//代入y,求x= (y-m)/k
	if (poLine->k == 0)
	{//水平线，难求
		if (poLine->m == y)
			return 0;	//随便返回一个了事
		else
			return NAN;
	}
	else if (poLine->m_iMode == Line_1::Mode::Degree_90)
		return poLine->x0;
	else
		return (y - poLine->m) / poLine->k;
}
/**********************一组点云函数****************************/
template<typename _T>void bSave_PLY(const char* pcFile, Point_Cloud<_T> oPC)
{
	bSave_PLY(pcFile, oPC.m_pPoint, oPC.m_iCount, oPC.m_pColor);
}
template<typename _T>void Init_Point_Cloud(Point_Cloud<_T>* poPC,int iMax_Count,int bHas_Color)
{
	int iSize;
	Point_Cloud<_T> oPC = {};
	oPC.m_iMax_Count = iMax_Count;
	oPC.m_bHas_Color = bHas_Color;
	iSize = iMax_Count * 3 * sizeof(_T);
	if (bHas_Color)
		iSize += ALIGN_SIZE_128(iMax_Count) + iMax_Count * 3;

	oPC.m_pBuffer = (unsigned char*)pMalloc(iSize);
	iSize = iMax_Count;
	oPC.m_pPoint = (_T(*)[3])oPC.m_pBuffer;
	oPC.m_pColor = (unsigned char(*)[3])(oPC.m_pPoint + iMax_Count);
	oPC.m_pColor = (unsigned char(*)[3])ALIGN_ADDR_128(oPC.m_pColor);
	memset(oPC.m_pPoint, 0, iMax_Count * 3 * sizeof(_T));
	memset(oPC.m_pColor, 255, iMax_Count * 3);

	*poPC = oPC;
	return;
}
template<typename _T>void Free_Point_Cloud(Point_Cloud<_T>* poPC)
{
	if (poPC->m_pBuffer)
		Free(poPC->m_pBuffer);
	*poPC = {};
}
template<typename _T>void Draw_Point(Point_Cloud<_T>* poPC, _T x, _T y, _T z, int R, int G, int B)
{
	Point_Cloud<_T> oPC = *poPC;
	if (oPC.m_iCount + 1 > oPC.m_iMax_Count)
	{
		printf("Insufficient memroy\n");
		return;
	}
	_T* pPoint = oPC.m_pPoint[oPC.m_iCount];
	pPoint[0] = x;
	pPoint[1] = y;
	pPoint[2] = z;

	if (oPC.m_pColor)
	{
		unsigned char* pColor = oPC.m_pColor[oPC.m_iCount];
		pColor[0] = R, pColor[1] = G, pColor[2] = B;
	}

	oPC.m_iCount++;
	*poPC = oPC;
	return;
}

template<typename _T>void Gen_Sphere(_T(**ppPoint_3D)[3], int* piCount, _T r ,int iStep_Count)
{//生成一个球，半径为1
 // x = r*sin(beta)*cos(alpha)
 // y = r*sin(beta)*sin(alpha)
 // z = r*cos(beta)
	_T alpha, beta, (*pPoint_3D)[3], x, y, z, fDelta_2, r_sin_beta, fDelta_2_Start;
	int i, j, iCur_Point = 0, iStep_Count_1;
	const _T fDelta_1 = (_T)(PI / iStep_Count);
	Line_1 oLine;
	pPoint_3D = (_T(*)[3])pMalloc(&oMem_Mgr, (iStep_Count * iStep_Count/2) * 3 * sizeof(_T));

	Cal_Line(&oLine, 0.f, 1.f, (float)iStep_Count, (float)iStep_Count);
	for (i = 0; i <= (iStep_Count >> 1); i++)
	{//外层是z
		beta = i * fDelta_1;
		z = r* (_T)cos(beta);
		iStep_Count_1 = (int)fGet_Line_y(&oLine, (float)i);
		fDelta_2 = (_T)(PI * 2.f / iStep_Count_1);
		r_sin_beta = r* (_T)sin(beta);
		fDelta_2_Start = (iRandom() % 100) * (3.14f / 100.f);
		for (j = 0; j < iStep_Count_1; j++)
		{
			alpha = j * fDelta_2 + fDelta_2_Start;
			x = r_sin_beta * (_T)cos(alpha);
			y = r_sin_beta * (_T)sin(alpha);
			pPoint_3D[iCur_Point][0] = x;
			pPoint_3D[iCur_Point][1] = y;
			pPoint_3D[iCur_Point][2] = z;
			iCur_Point++;
			if (z > 0)
			{
				pPoint_3D[iCur_Point][0] = x;
				pPoint_3D[iCur_Point][1] = y;
				pPoint_3D[iCur_Point][2] = -z;
				iCur_Point++;
			}
		}
		//printf("%d\n", iStep_Count_1);
	}
	Shrink(&oMem_Mgr, pPoint_3D, iCur_Point * 3 * sizeof(_T));
	if (ppPoint_3D)
		*ppPoint_3D = pPoint_3D;
	if (piCount)
		*piCount = iCur_Point;
	return;
}
template<typename _T>void Draw_Rect(Point_Cloud<_T>* poPC, _T x, _T y, _T z,_T w, _T h, int iStep, int R, int G, int B)
{
	Point_Cloud<_T> oPC = *poPC;
	Draw_Line(&oPC, x, y, z, x + w, y, z,(int)(w/ iStep), R, G, B);
	Draw_Line(&oPC, x + w, y, z, x + w, y + h, z, (int)(h / iStep), R, G, B);
	Draw_Line(&oPC, x + w, y + h, z, x, y + h, z, (int)(w / iStep), R, G, B);
	Draw_Line(&oPC, x, y + h, z, x, y, z, (int)(h / iStep), R, G, B);
	*poPC = oPC;
	//bSave_PLY("c:\\tmp\\1.ply", oPC);
	return;
}
template<typename _T>void Draw_Sphere(Point_Cloud<_T>* poPC, _T x, _T y, _T z,_T r, int iStep_Count,int R,int G,int B)
{
	Point_Cloud<_T> oPC = *poPC;
	_T (*pPoint_3D)[3];
	int i,iCount,iPos;
	Gen_Sphere(&pPoint_3D, &iCount, r,iStep_Count);
	if (oPC.m_iCount + iCount > oPC.m_iMax_Count)
	{
		printf("Insufficient memroy");
		return;
	}
	iPos = oPC.m_iCount;
	for (i = 0; i < iCount; i++,iPos++)
	{
		_T* pPoint_3D_1 = oPC.m_pPoint[iPos];
		unsigned char* pColor = oPC.m_pColor[iPos];
		memcpy(pPoint_3D_1, pPoint_3D[i], 3 * sizeof(_T));
		pPoint_3D_1[0] += x;
		pPoint_3D_1[1] += y;
		pPoint_3D_1[2] += z;
		pColor[0] = R, pColor[1] = G, pColor[2] = B;
	}
	//memcpy(oPC.m_pPoint + oPC.m_iCount, pPoint_3D, iCount * 3*sizeof(_T));
	oPC.m_iCount += iCount;
	*poPC = oPC;
	Free(pPoint_3D);
	return;
}
template<typename _T>void Draw_Line(Point_Cloud<_T>* poPC, _T x0, _T y0, _T z0, _T x1, _T y1, _T z1, int iCount, int R, int G, int B)
{
	int i;
	for (i = 0; i < iCount; i++)
	{
		_T Pos[3] = {x0+ (x1 - x0) * i / iCount,
			y0+ (y1 - y0) * i / iCount,
			z0+(z1 - z0) * i / iCount };
		Draw_Point(poPC, Pos[0], Pos[1], Pos[2], R, G, B);
	}
	return;
}
/**********************一组点云函数****************************/

int iRead_Line(FILE* pFile, char Line[], int iLine_Size)
{//从文件读入一行,返回一行大小
	int iResult;
	while (iResult=(int)fread((void*)Line, 1, 1, pFile))
		if (Line[0] != '\r' && Line[0] != '\n' && Line[0] != '\t' && Line[0] != ' ')
			break;
	if(!iResult)
		return 0;
	int i = 0;
	while (iResult=(int)fread((void*)&Line[++i], 1, 1, pFile))
	{
		if (Line[i] == '\r' || Line[i] == '\n')
			break;
		if (i > iLine_Size)
			return 0;
	}
	if(i<iLine_Size)
		Line[i] = 0;

	return iResult || i?i:0;
}
