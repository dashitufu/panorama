//这个文件管所有独立于特殊平台的通用函数
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "sys/timeb.h"
#include "math.h"
#include <mutex>

extern "C"
{
#include "Buddy_System.h"
}

#ifndef WIN32
	#include <unistd.h>	//给sleep用
#endif

#define PI 3.14159265358979323846

#define Abs(A) ((A)>=0?(A):(-(A)))
#ifndef Min 
#define Min(A,B)( A<=B?A:B)
#endif
#ifndef Max
#define Max(A,B)(A>=B?A:B)
#endif
#ifndef Clip3
#define Clip3(x,y,z)  ( (z)<(x)?(x): (z)>(y)?(y):(z))
#endif

#define ALIGN_SIZE_8(iSize) ((( (unsigned long long)(iSize)+7)>>3)<<3) 
#define ALIGN_SIZE_128(iSize) ((( (unsigned long long)(iSize)+127)>>7)<<7) 
#define ALIGN_SIZE_1024(iSize) ((( (unsigned long long)(iSize)+1023)>>10)<<10) 
#define ALIGN_ADDR_128(pAddr) (((((unsigned long long)(pAddr))+127)>>7)<<7)

#define bGet_Bit(pBuffer, iBit_Pos) (pBuffer)[(iBit_Pos) >> 3] & (1 << ((iBit_Pos) & 0x7))
//int bGet_Bit(unsigned char* pBuffer, int iBit_Pos)
//{	return pBuffer[iBit_Pos >> 3] & (1 << (iBit_Pos & 0x7));}
#define Set_Bit(pBuffer, iBit_Pos) \
{\
	(pBuffer)[(iBit_Pos) >> 3] |= (1 << ((iBit_Pos) & 0x7)); \
}

#define Attach_Light_Ptr(oPtr,pBuffer,iSize,iGPU_ID) oPtr = { (int)0,(int)(iSize),pBuffer,iGPU_ID }

////Light Ptr在此分配
//#define Malloc(oPtr, iSize,pBuffer) \
//{ \
//	(pBuffer)= (unsigned char*)(oPtr).m_pBuffer+(oPtr).m_iCur; \
//	(oPtr).m_iCur +=ALIGN_SIZE_128((iSize)); \
//	if ( (oPtr).m_iCur > (oPtr).m_iMax_Buffer_Size) \
//	{ \
//		(oPtr).m_iCur -= ALIGN_SIZE_128((iSize));\
//		(pBuffer)=NULL; \
//		printf("Fail to allocate memory in Malloc, Total:%d Remain:%d Need:%d\n",(int)(oPtr).m_iMax_Buffer_Size,(int)(oPtr).m_iMax_Buffer_Size-(oPtr).m_iCur,(int)(iSize)); \
//	} \
//}

template<typename _T> struct Point_2D {
	unsigned int m_iCamera_Index;
	unsigned int m_iPoint_Index;
	_T m_Pos[2];
};

typedef struct BitPtr {
	unsigned char* m_pBuffer;
	int m_iCur;		//当前字节,位于m_oBuffer[]中的第m_iCur Byte.
	int m_iBitPtr;	//当前字节中的当前位
	int m_iEnd;		//m_Buffer的合法数据是有范围的, 合法数据的最后一个字节位于m_iEnd-1, 如果m_iCur>=m_iEnd则为越界
}BitPtr;

//尝试做一个可以增加增加删除点的结构
template<typename _T> struct Point_Cloud {
	_T(*m_pPoint)[3];				//记录点的位置(x,y,z)
	unsigned char (*m_pColor)[3];	//RGB
	unsigned char *m_pBuffer;		//内存所在
	unsigned char m_bHas_Color : 1;	//标志是否有颜色
	int m_iMax_Count;				//m_pBuffer最多容纳多少点
	int m_iCount;					//目前有多少点
};

typedef struct Line_1 {		//重做直线结构，全部用浮点表示
	enum Mode {
		Degree_Other = 0,
		Degree_90,
		Degree_180
	};
	float x0, y0, x1, y1;	//未必有值，因为直线也可以表示为y=kx+m
	float a, b, c;
	float k, m;		//y=kx+m表示
	float theta;	//倾角，当theta=PI/2时，y=kx+m无效,只能表示为x=m
	Mode m_iMode;	//=0时， y=kx+m =1时，x=m;
}Line_1;

extern Mem_Mgr oMem_Mgr;;

//void Set_Bit(unsigned char* pBuffer, int iBit_Pos)
//{	pBuffer[iBit_Pos >> 3] |= (1 << (iBit_Pos & 0x7)); }
const int WriteBitsMask[] = { 0,0x80,0xC0,0xE0,0xF0,0xF8,0xFC,0xFE,0xFF };
const int WriteBitsMask2[] = { 0xFF,0x7F,0x3F,0x1F,0x0F,0x07,0x03,0x01 };

#define WriteBits2(oBitPtr,iLen,iValue1)	\
{\
	unsigned int iValue=(iValue1)<<(32- (oBitPtr).m_iBitPtr-(iLen));	\
	unsigned char *pCur=(oBitPtr).m_pBuffer+ (oBitPtr).m_iCur;	\
	unsigned int iTemp1=(oBitPtr).m_iBitPtr+iLen;\
	*pCur= (*pCur &WriteBitsMask[(oBitPtr).m_iBitPtr]) | ((iValue>>24)&WriteBitsMask2[(oBitPtr).m_iBitPtr]);	\
	if(iTemp1>=8)									\
	{											\
		pCur[1]=(iValue & 0x00FF0000)>>16;		\
		pCur[2]=(iValue & 0x0000FF00)>>8;		\
		pCur[3]=(iValue & 0x000000FF);			\
		(oBitPtr).m_iCur+=iTemp1>>3;			\
		(oBitPtr).m_iBitPtr=iTemp1 & 0x7;		\
	}else										\
		(oBitPtr).m_iBitPtr=iTemp1;			\
}

unsigned long long iGet_File_Length(char* pcFile);
unsigned long long iGet_Tick_Count();
int iGet_File_Count(const char* pcPath);	//获取一个目录所有文件
void Get_All_File(const char* pcPath, char* pBuffer);
int bSave_Bin(const char* pcFile, float* pData, int iSize);
int bLoad_Raw_Data(const char* pcFile, unsigned char** ppBuffer, int* piSize=NULL);
int bLoad_Raw_Data(const char* pcFile, unsigned char** ppBuffer, int iSize = 0, int bNeed_Malloc = 1, int iFrame_No = 0);
int bLoad_Text_File(const char* pcFile, char** ppBuffer, int* piSize=NULL);
int bSave_Raw_Data(const char* pcFile, unsigned char* pBuffer, int iSize);

//上三角坐标转换为索引值
int iUpper_Triangle_Cord_2_Index(int x, int y, int w);

//上三角有效元数个数
int iGet_Upper_Triangle_Size(int w);

//一组有的没的随机数生成
int iRandom(int iStart, int iEnd);
int iRandom();

//早晚得废
//以后废弃的临时函数
template<typename _T>
void Temp_Load_Match_Point(_T(**ppPoint_1)[2], _T(**ppPoint_2)[2], int* piCount)
{
	_T(*pPoint_1)[2], (*pPoint_2)[2];
	int i, iCount = (int)iGet_File_Length((char*)"c:\\tmp\\2.bin") / (4 * sizeof(float));
	pPoint_1 = (_T(*)[2])malloc(iCount * 2 * sizeof(_T));
	pPoint_2 = (_T(*)[2])malloc(iCount * 2 * sizeof(_T));
	FILE* pFile = fopen("c:\\tmp\\1.bin", "rb");
	for (i = 0; i < iCount; i++)
	{
		float Data[2];
		fread(Data, 1, 2 * sizeof(float), pFile);
		pPoint_1[i][0] = (_T)Data[0];
		pPoint_1[i][1] = (_T)Data[1];
		fread(Data, 1, 2 * sizeof(float), pFile);
		pPoint_2[i][0] = (_T)Data[0];
		pPoint_2[i][1] = (_T)Data[1];		
	}
	fclose(pFile);
	*ppPoint_1 = pPoint_1, * ppPoint_2 = pPoint_2;
	*piCount = iCount;
}

void Init_BitPtr(BitPtr* poBitPtr, unsigned char* pBuffer, int iSize);

int iGetBits(BitPtr* poBitPtr, int iLen);

//解决两组基本数据类型的第n大与及快速排序，待优化
template<typename _T> _T oGet_Nth_Elem(_T Seq[], int iCount, int iNth);
template<typename _T> void Quick_Sort(_T Seq[], int iStart, int iEnd);
template<typename _T>int bSave_PLY(const char* pcFile, _T Point[][3], int iPoint_Count,unsigned char Color[][3]=NULL, int bText=1);

//Temp code 
template<typename _T>void Temp_Load_File(const char* pcFile, _T(**ppPoint_3D_1)[3], _T(**ppPoint_3D_2)[3], int* piCount);
template<typename _T>void Temp_Load_File_1(const char* pcFile, _T(**ppPoint_3D_1)[3], _T(**ppPoint_3D_2)[3], int* piCount);
template<typename _T>void Temp_Load_File_2(int* piCameta_Count, int* piPoint_Count, int* piObservation_Count, Point_2D<_T>** ppPoint_2D, _T(**ppPoint_3D)[3], _T(**ppCamera)[3 * 3]);
//以下函数也没啥用，第九讲用
template<typename _T>void Normalize(_T Point_3D[][3], Point_2D<_T> Point_2D[], int iPoint_Count, _T Camera[][9], int iCamera_Count);

//分配内存
//Input： iSize: 要分配内存的字节数
//return: void 指针，用户自己转换为自己的类型
//void* pMalloc(Mem_Mgr* poMem_Mgr, unsigned int iSize);
//void* pMalloc_1(Mem_Mgr* poMem_Mgr, unsigned int iSize);
void Init_Env();
void Free_Env();
void Init_Env_CPU();
void Free_Env_CPU();
void* pMalloc(unsigned int iSize);

void Shrink(void* p, unsigned int iSize);
void Free(void* p);
void Disp_Mem();

//直线函数
void Cal_Line(Line_1* poLine, float x0, float y0, float x1, float y1);
float fGet_Line_x(Line_1* poLine, float y);
float fGet_Line_y(Line_1* poLine, float x);

//用于点云简单画图
template<typename _T>void bSave_PLY(const char* pcFile, Point_Cloud<_T> oPC);
template<typename _T>void Init_Point_Cloud(Point_Cloud<_T>* poPC, int iMax_Count, int bHas_Color=0);
template<typename _T>void Free_Point_Cloud(Point_Cloud<_T>* poPC);
template<typename _T>void Draw_Point(Point_Cloud<_T>* poPC, _T x, _T y, _T z, int R = 255, int G = 255, int B = 255);
template<typename _T>void Draw_Sphere(Point_Cloud<_T>* poPC, _T x, _T y, _T z, _T r = 1.f, int iStep_Count = 40, int R = 255, int G = 255, int B = 255);
template<typename _T>void Draw_Line(Point_Cloud<_T>* poPC, _T x0, _T y0, _T z0, _T x1, _T y1, _T z1, int iCount = 50, int R = 255, int G = 255, int B = 255);
template<typename _T>void Draw_Camera(Point_Cloud<_T>* poPC, _T T[4 * 4], int R = 255, int G = 255, int B = 255);
template<typename _T>void Draw_Rect(Point_Cloud<_T>* poPC, _T x, _T y, _T z, _T w, _T h, int iStep = 1, int R = 255, int G = 255, int B = 255);

int iRead_Line(FILE* pFile, char Line[], int iLine_Size);
int bStricmp(char* pStr_0, char* pStr_1);
