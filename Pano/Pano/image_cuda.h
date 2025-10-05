#pragma once    //这个少了会编译出错

#include "image.h"
#include "cuda_runtime.h"

extern "C"
{
#include "Buddy_System.h"
}

#define GET_THREAD_ID() (blockDim.x * blockIdx.x + threadIdx.x)

void Init_Image_GPU(Image* poImage, int iWidth, int iHeight, Image::Type iType, int iBit_Count, Light_Ptr* poPtr=NULL, int iGPU_ID=0);
void Free_Image_GPU(Image* poImage);
void Copy_Image_To_CPU(Image oOrg, Image oNew);
void Copy_Image_To_GPU(Image oOrg, Image oNew);
int bSave_Image_GPU(const char* pcFile, Image oImage);
int bSave_Image_GPU(const char* pcFile, Image* poHeader_GPU);
int bSave_Comp_GPU(const char* pcFile, Image oImage, int iComp);
int bSave_Comp_GPU(const char* pcFile, Image* poHeader_GPU, int iComp);

//双线性组处理,接口很差
void Bi_Linear_cv_GPU(Image Source[], Image Dest[], int iCount, 
    int w_s, int h_s, int w_d, int h_d, int iChannel,
    float fScale_x, float fScale_y, Border_Type iBorder_Type);
void Set_Color_GPU(Image oImage, int R = 0, int G = 0, int B = 0, unsigned long long iStream = 0);

//边扩展
void Copy_Make_Border_GPU(Image oSource, Image oDest, short iLeft, short iTop, short iRight, short iBottom);

void Pyr_Down_GPU(Image oSource, Image oDest, unsigned short* pAux=NULL);
void Pyr_Up_GPU(Image oSource, Image oDest);
