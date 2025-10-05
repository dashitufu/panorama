#pragma once
//#include "Common.h"

extern "C"
{
#include "Buddy_System.h"
}

#define Clip(iValue) (( (iValue)&0xFFFFFF00)==0?(iValue): (iValue)<0?0:255)	//��!
#define Clip_16(iValue) (( (iValue)&0xFFFF0000)==0?(iValue): (iValue)<0?0:65535)	//��!
#define Clip3(x,y,z)  ( (z)<(x)?(x): (z)>(y)?(y):(z))
//#define Min(a,b) ((a)<(b)?(a):(b))
#define SWAP(A,B,Temp) Temp=A,	A = B,	B=Temp;

#ifdef PI
#undef PI
#endif

#define PI 3.14159265358979323846
#define S_FUNC_MAP_SIZE 128

#define _RGB_2_YUV(R, G, B, Y, U, V) \
{ \
	Y = (int)(0.256788f * R + 0.504129f * G + 0.097906f * B + 16.01f); \
	U = (int)(-0.148223f * R - 0.290993f * G + 0.439216f * B + 128.01f); \
	V = (int)(0.439216f * R - 0.367788f * G - 0.071427f * B + 128.01f); \
	Y = Clip(Y); \
	U = Clip(U); \
	V = Clip(V); \
} \

#define _YUV_2_RGB( Y, U, V,R, G, B) \
{ \
	R = (int)(1.164383 * (Y - 16) + 1.596027 * (V - 128)); \
	G =(int)( 1.164383 * (Y - 16) - 0.391762 * (U - 128) - 0.812968 * (V - 128)); \
	B = (int)(1.164383 * (Y - 16) + 2.017232 * (U - 128));\
	R = Clip(R);	\
	G = Clip(G);	\
	B = Clip(B);	\
}

enum Border_Type {
	BORDER_CONSTANT = 0, //!< `iiiiii|abcdefgh|iiiiiii`  with some specified `i`
	BORDER_REPLICATE = 1, //!< `aaaaaa|abcdefgh|hhhhhhh`
	BORDER_REFLECT = 2, //!< `fedcba|abcdefgh|hgfedcb`
	BORDER_WRAP = 3, //!< `cdefgh|abcdefgh|abcdefg`
	BORDER_REFLECT_101 = 4, //!< `gfedcb|abcdefgh|gfedcba`
	BORDER_TRANSPARENT = 5, //!< `uvwxyz|abcdefgh|ijklmno`

	BORDER_REFLECT101 = BORDER_REFLECT_101, //!< same as BORDER_REFLECT_101
	BORDER_DEFAULT = BORDER_REFLECT_101, //!< same as BORDER_REFLECT_101
	BORDER_ISOLATED = 16 //!< do not look outside of ROI
};

enum Interpolation_Flag {
	INTER_NEAREST = 0,  //�����    
	INTER_LINEAR = 1,   //Bi_Linear��˫����
	INTER_AREA = 3,
	WARP_RELATIVE_MAP = 32
};

typedef struct Pixel_8 {
	union {
		unsigned char Data[8];
		char Data_c[8];
	};
}Pixel_8;
typedef struct Pixel_4 {
	union {
		unsigned char Data[4];
		char Data_c[4];
	};	
}Pixel_4;

typedef struct BitMap_8Bit {
#pragma pack(2)
	typedef struct tagBITMAPFILEHEADER {
		unsigned short    bfType;
		unsigned int   bfSize;
		unsigned short    bfReserved1;
		unsigned short    bfReserved2;
		unsigned int   bfOffBits;
	} BITMAPFILEHEADER;

	typedef struct tagBITMAPINFOHEADER {
		unsigned int      biSize;
		int       biWidth;
		int       biHeight;
		unsigned short       biPlanes;
		unsigned short       biBitCount;
		unsigned int      biCompression;
		unsigned int      biSizeImage;
		int       biXPelsPerMeter;
		int       biYPelsPerMeter;
		unsigned int      biClrUsed;
		unsigned int      biClrImportant;
	} BITMAPINFOHEADER;
#pragma pack()
	typedef struct Part_1 {
		unsigned char* m_pRGBA[4];
		int m_iWidth, m_iHeight;
	}Part_1;

	union {
		struct {
			unsigned char* m_pRGBA[4];
			int m_iWidth, m_iHeight;
		};
		Part_1 m_oPart_1;
	};
	unsigned int m_iMax_Buffer_Size : 28;	//����ڴ�����
	unsigned int m_iGPU_ID : 4;	//�������15�ſ�

	unsigned short m_iBit_Count : 6;	//8,16,24,32
	unsigned short m_iChannel_Count : 3;
	unsigned short m_iMem_Src : 1;	//��ͼ��������CPU�ֻ�����GPU��
	unsigned short m_bIs_Attached : 1;	//��ͼ���Ƿ�Attach�ⲿBuffer
	unsigned char Reserved[88];	//�ճ�136���ֽڲſ죬�ܹ��죬��ʱ������1650ti�����������
}BitMap_8Bit;

typedef struct YUV_Image_444_8Bit {

	typedef struct {
		unsigned char* m_pYUV[4];
		int m_iWidth, m_iHeight;
	}Part_1;

	union {
		struct {
			unsigned char* m_pYUV[4];
			int m_iWidth, m_iHeight;
		};
		Part_1 m_oPart_1;
	};

	unsigned int m_iMax_Buffer_Size : 28;	//����ڴ�����
	unsigned int m_iGPU_ID : 4;	//�������15�ſ�

	unsigned short m_Reserve : 6;	//8,16,24,32
	unsigned short m_iChannel_Count : 3;
	unsigned short m_iMem_Src : 1;	//��ͼ��������CPU�ֻ�����GPU��
	unsigned short m_bIs_Attached : 1;	//�Ƿ������ⲿBuffer
}YUV_Image_444_8Bit;

typedef struct Image {	//��һ����һͳ��Image
	typedef enum {
		IMAGE_TYPE_BMP,
		IMAGE_TYPE_YUV_444
	}Type;

	typedef struct {
		unsigned char* m_pChannel[4];
		int m_iWidth, m_iHeight;
	}Part_1;

	union {
		BitMap_8Bit m_oBitMap;
		YUV_Image_444_8Bit m_oYUV;
		struct {
			union {
				struct {
					unsigned char* m_pChannel[4];
					int m_iWidth, m_iHeight;
				};
				Part_1 m_oPart_1;
			};

			unsigned int m_iMax_Buffer_Size : 28;	//����ڴ�����
			unsigned int m_iGPU_ID : 4;	//�������15�ſ�

			unsigned short m_iBit_Count : 6;	//8,16,24,32
			unsigned short m_iChannel_Count : 3;
			unsigned short m_iMem_Src : 1;	//��ͼ��������CPU�ֻ�����GPU��
			unsigned short m_bIs_Attached : 1;	//�Ƿ������ⲿBuffer
			unsigned short m_iImage_Type : 2;		//0:BMP 1:YUV_444

			unsigned char* m_pBuffer;			//������Buffer��ʼ֮���������Attached,�����ΪNULL
		};
	};
}Image;

typedef struct YUV_Image_420_8Bit {
	typedef struct Part_1 {
		unsigned char* m_pYUV[3];
		int m_iWidth, m_iHeight;
	}Part_1;
	union {
		struct {
			unsigned char* m_pYUV[3];
			int m_iWidth, m_iHeight;
		};
		Part_1 m_oPart_1;
	};

	unsigned int m_iMax_Buffer_Size : 26;	//����ڴ�����
	unsigned int m_iGPU_ID : 4;	//�������15�ſ�
	unsigned int m_iMem_Src : 1;	//��ͼ��������CPU�ֻ�����GPU��
	unsigned int m_bIs_Attached : 1;	//�Ƿ������ⲿBuffer
}YUV_Image_420_8Bit;

enum Mem_Src {
	CPU,
	GPU
};

enum class Image_Format {
	YUV_420,
	YUVA_420,
	YUV_444,
	YUVA_444,
	RGBA,
	RGB,
	BGR,
	BGRA,
	Mono,
	None	//ɶҲ���ǣ�һ��������Ϊ�Ƿ�ֵ
};

//typedef struct Line_1 {		//����ֱ�߽ṹ��ȫ���ø����ʾ
//	enum Mode {
//		Degree_Other = 0,
//		Degree_90,
//		Degree_180
//	};
//	float x0, y0, x1, y1;	//δ����ֵ����Ϊֱ��Ҳ���Ա�ʾΪy=kx+m
//	float a, b, c;
//	float k, m;		//y=kx+m��ʾ
//	float theta;	//��ǣ���theta=PI/2ʱ��y=kx+m��Ч,ֻ�ܱ�ʾΪx=m
//	Mode m_iMode;	//=0ʱ�� y=kx+m =1ʱ��x=m;
//}Line_1;

//һЩͨ�ú���
void Get_Image_Info(const char* pcFile, Image* poImage);
void Init_Image(Image* poImage, int iWidth, int iHeight, Image::Type iType, int iBit_Count);
int bLoad_Image(const char* pcFile, Image* poImage, int iWidth = 0, int iHeight = 0, int iFrame = 0, int bNeed_Malloc = 1, Mem_Mgr* poMem_Mgr=NULL);
int bLoad_Image_GPU(const char* pcFile, Image* poImage, int iWidth=0, int iHeight=0, int iFrame=0);
int bLoad_Comp_GPU(const char* pcFile, Image oDest, int iComp);

int bSave_Image(const char* pcFile, Image oImage, int iFormat = -1);
int bSave_Image(const char* pcFile, float* pImage, int iWidth, int iHeight);
int bSave_Comp(const char* pcFile, Image oImage, int iComp);

void Free_Image(Image* poImage);
void Free_Image(Mem_Mgr* poMem_Mgr, Image oImage);
void Free_Image(BitMap_8Bit* poBitMap);
void Attach_Buffer(Image* poImage, unsigned char* pBuffer, int iWidth, int iHeight, int iChannel_Count, int iImage_Type);

void Set_Color(Image oImage, int R = 0, int G = 0, int B = 0, int A = 255);
//void Cal_Line(Line_1* poLine, float x0, float y0, float x1, float y1);
//float fGet_Line_x(Line_1* poLine, float y);
//float fGet_Line_y(Line_1* poLine, float x);
void Fill_Region(Image oImage, int iSeed_x, int iSeed_y, int iBack_Color = 0);
void Draw_Line(Image oImage, int x0, int y0, int x1, int y1, int R = 255, int G = 255, int B = 255);
void Draw_Point(Image oImage, int x, int y, int r = 10, int R = 255, int G = 255, int B = 255);
void Draw_Arc(Image oImage, int r, int iCenter_x, int iCenter_y, float fAngle_Start=0, float fAngle_End=PI*2, int R=255, int G=255, int B=255);
void Draw_Rect(Image oImage, int x, int y, int w, int h, int iThickness=1, int r = 255, int g = 255, int b = 255);
template<typename _T>void Get_Bounding_Box(_T Point[][2], int iPoint_Count, _T Bounding_Box[2][2]);

void RGB_2_Gray(Image oImage, float* pImage);	//RGBת�Ҷ�
void Float_2_Image(float* pImage, Image oImage);	//������������ת��ΪImage
void Crop_Image(Image oSource, int x, int y, int w, int h, Image oDest);

//��˹ģ��
void Gen_Gauss_Filter(int r, float fSigma, float** ppFilter);
void Gen_Gauss_Filter(int r, float** ppFilter);		//����r�Ĵ�Сȷ��sigma
void Gauss_Filter_Ref(float* pSource, int iWidth, int iHeight, float* pDest, int r, float* pFilter);
void Gauss_Filter_AVX512(float* pSource, int iWidth, int iHeight, float* pDest, float* pAux, int r, float* pFilter);
void Transpose_AVX512(float* pSource, int iWidth, int iHeight, float* pDest);
void Box_Filter_Ref(Image oSource, Image oDest, int r);	//�����˲��ο�����

//�����˲�
void Box_Filter_Ref(Image oSource, Image oDest, int r);
void Box_Filter(Image oSource, Image oDest, int r);		//����������
void Dilate_Bin_Ref(Image oSource, Image oDest, int r = 1);	//����
void Dilate_Bin(Image oSource, int r = 1, Image oDest = {});			//���ٷ���
int Compare_Image(Image oSource, Image oDest, int iDiff_Threshold=0);
int Compare_Image(const char* pcSource, const char* pcDest, int iDiff_Threshold = 0);
void Place_Image(Image oTile, Image oScreen, int x, int y);
void Box_Filter_1(Image oSource, Image oDest, int r);	//������з���

//ƴ��ͼƬ
void Concat_Image(Image oA, Image oB, Image* poC, int iFlag = 0);
template<typename _T>void Draw_Match_Point(_T Point_1[][2], _T Point_2[][2], unsigned char Mask[], int iCount, Image oImage);
template<typename _T>void Draw_Match_Point(_T Point_1[][2], _T Point_2[][2], int iCount, Image oImage);

//��ֵ����
int iGet_Border_x(int x, int iWidth, Border_Type iBorder_Type);
int iGet_Border_y(int y, int iHeight, Border_Type iBorder_Type);
void Bi_Linear_Ref(Image oSource, Image oDest, Border_Type iBorder_Type = BORDER_REFLECT_101);
void Bi_Linear_cv(Image oSource, Image oDest, float fScale_x, float fScale_y, Border_Type iBorder_Type = BORDER_REFLECT_101);

template<typename _T>void Bi_Linear(_T Source[], int w_s, int h_s, int w_d, int h_d, _T Dest[]);
//���з�������padding �þ���
template<typename _T>void Sep_Filter_2D(_T A[], int w, int h, _T Ker_x[], _T Ker_y[], int iKernel_Size, _T B[]);

//��һ���ſ���������
template<typename Source_Type, typename Dest_Type, typename Kernel_Type>void Sep_Filter_2D_1
(Source_Type A[], int w, int h, Kernel_Type Ker_x[], Kernel_Type Ker_y[], int iKernel_Size, Dest_Type B[]);

//�������ϲ���
void Pyr_Up_Ref(Image oSource, Image oDest);
void Pyr_Down_Ref(Image oSource, Image oDest, Border_Type iBorder_Type = BORDER_REFLECT101);