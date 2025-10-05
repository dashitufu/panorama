#include "cuda_runtime.h"
#include "Common.h"
#include "Common_Cuda.cuh"
//#include "Image.h"
#include "image_cuda.h"
#include "Matrix.h"

void SB_Image_Cuda()
{

}
//__device__ static int iGet_Border_y_GPU(int y, int iHeight, Border_Type iBorder_Type)
//{
//    if (y < 0)
//    {
//        switch (iBorder_Type)
//        {
//        case Border_Type::BORDER_CONSTANT:
//            return -1;
//        case Border_Type::BORDER_REFLECT:
//            return -y - 1;
//        case Border_Type::BORDER_REFLECT_101:
//            return -y;
//        }
//    }
//    else if (y >= iHeight)
//    {
//        switch (iBorder_Type)
//        {
//        case Border_Type::BORDER_CONSTANT:
//            return -1;
//        case Border_Type::BORDER_REFLECT:
//            return iHeight - (y - iHeight + 1);
//        case Border_Type::BORDER_REFLECT_101:
//            return iHeight - (y - iHeight + 1);
//        }
//    }
//    return y;
//}
__device__ static int iGet_Border_y_GPU(int y, int iHeight, Border_Type iBorder_Type, int iThread_ID = 0)
{
    if (y < 0)
    {
        switch (iBorder_Type)
        {
        case Border_Type::BORDER_REFLECT:
            return -y - 1;
        case Border_Type::BORDER_REFLECT_101:
            return -y;
        case Border_Type::BORDER_CONSTANT:
        case Border_Type::BORDER_REPLICATE:
            return 0;
        }
    }
    else if (y >= iHeight)
    {
        /*if (iThread_ID == 77803)
            printf("Here");*/
        switch (iBorder_Type)
        {

        case Border_Type::BORDER_REFLECT:
            return iHeight - (y - iHeight + 1);
        case Border_Type::BORDER_REFLECT_101:
            return iHeight - (y - iHeight + 2);
        case Border_Type::BORDER_CONSTANT:
        case Border_Type::BORDER_REPLICATE:
            return iHeight - 1;
        }
    }
    return y;
}

//__device__ static int iGet_Border_x_GPU(int x, int iWidth, Border_Type iBorder_Type)
//{
//    if (x < 0)
//    {
//        switch (iBorder_Type)
//        {
//        case Border_Type::BORDER_REFLECT:
//            return -x - 1;
//        case Border_Type::BORDER_REFLECT_101:
//            return -x;
//        }
//    }
//    else if (x >= iWidth)
//    {
//        switch (iBorder_Type)
//        {
//        case Border_Type::BORDER_REFLECT:
//            return iWidth - (x - iWidth + 1);
//        case Border_Type::BORDER_REFLECT_101:
//            return iWidth - (x - iWidth + 2);
//        }
//    }
//    return x;
//}

__device__ static int iGet_Border_x_GPU(int x, int iWidth, Border_Type iBorder_Type)
{
    if (x < 0)
    {
        switch (iBorder_Type)
        {
        case Border_Type::BORDER_REFLECT:
            return -x - 1;
        case Border_Type::BORDER_REFLECT_101:
            return -x;
        case Border_Type::BORDER_CONSTANT:
        case Border_Type::BORDER_REPLICATE:
            return 0;
        }
    }
    else if (x >= iWidth)
    {
        switch (iBorder_Type)
        {
        case Border_Type::BORDER_REFLECT:
            return iWidth - (x - iWidth + 1);
        case Border_Type::BORDER_REFLECT_101:
            return iWidth - (x - iWidth + 2);
        case Border_Type::BORDER_CONSTANT:
        case Border_Type::BORDER_REPLICATE:
            return iWidth - 1;
        }
    }
    return x;
}
void Init_Image_GPU(Image* poImage, int iWidth, int iHeight, Image::Type iType, int iBit_Count,Light_Ptr *poPtr,int iGPU_ID)
{
    Light_Ptr oPtr;
	int i,iSize, iSize_With_Remain;
	if (!poImage || iWidth == 0 || iHeight == 0)
	{
		printf("Invalid parameter in Init_Image_GPU\n");
		*poImage = { {0 } };
		return;
	}

	poImage->m_iWidth = iWidth;
	poImage->m_iHeight = iHeight;
	poImage->m_iBit_Count = iBit_Count;
	poImage->m_iChannel_Count = iBit_Count >> 3;
	poImage->m_iMem_Src = Mem_Src::GPU;
	poImage->m_iGPU_ID = iGPU_ID;

	iSize = iWidth * iHeight;
	poImage->m_iMax_Buffer_Size = iSize * poImage->m_iChannel_Count;

	//留有余地，1，起码留64字节以上；2，128字节倍数，以便GPU对齐
	iSize_With_Remain = ((poImage->m_iMax_Buffer_Size + 127) / 128) * 128;
	if (iSize_With_Remain - poImage->m_iMax_Buffer_Size < 64)
		iSize_With_Remain += 128;

	if (poPtr)
	{
        oPtr = *poPtr;
		Malloc(oPtr, iSize_With_Remain, poImage->m_pBuffer);
		poImage->m_pChannel[0]= poImage->m_pBuffer;
	}
	else
		poImage->m_pChannel[0] = poImage->m_pBuffer = (unsigned char*)pMalloc_GPU(iSize_With_Remain);
	
	for (i = 1; i < poImage->m_iChannel_Count; i++)
		poImage->m_pChannel[i] = poImage->m_pChannel[i - 1] + iSize;
	for (; i < 4; i++)
		poImage->m_pChannel[i] = NULL;

	if (poPtr)
		*poPtr = oPtr;
}
void Copy_Image_To_CPU(Image oOrg, Image oNew)
{
	int iSize;
	iSize = oOrg.m_iWidth * oOrg.m_iHeight * Min(oOrg.m_iChannel_Count, oNew.m_iChannel_Count);
	if (iSize <= 0 || !oOrg.m_pChannel[0] || !oNew.m_pChannel[0])
	{
		printf("Invalid parameter in Copy_Image_To_GPU\n");
		return;
	}
    /*for(int i=0;i<oNew.m_iChannel_Count;i++)
	    cudaMemcpyAsync(oNew.m_pChannel[i], oOrg.m_pChannel[i], oNew.m_iWidth * oNew.m_iHeight, cudaMemcpyDeviceToHost);*/
    cudaMemcpyAsync(oNew.m_pChannel[0], oOrg.m_pChannel[0], iSize, cudaMemcpyDeviceToHost);
    //Disp_Cuda_Error();
}
void Copy_Image_To_GPU(Image oOrg, Image oNew)
{
	int iSize;
	iSize = oOrg.m_iWidth * oOrg.m_iHeight * Min(oOrg.m_iChannel_Count, oNew.m_iChannel_Count);
	if (iSize <= 0 || !oOrg.m_pChannel[0] || !oNew.m_pChannel[0] ||
		oNew.m_iMem_Src != Mem_Src::GPU)
	{
		printf("Invalid parameter in Copy_Image_To_GPU\n");
		return;
	}
	cudaMemcpyAsync(oNew.m_pChannel[0], oOrg.m_pChannel[0], iSize, oOrg.m_iMem_Src == Mem_Src::CPU ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToDevice);
}
void Free_Image_GPU(Image* poImage)
{
	if (poImage && poImage->m_pBuffer)
	{
		Free_GPU(poImage->m_pBuffer);
		*poImage = {};
	}
	return;
}
int bLoad_Comp_GPU(const char* pcFile, Image oDest, int iComp)
{//暂时将单色图抄到目标图某通道上
    Image oImage;
    if (!bLoad_Image(pcFile, &oImage))
        return 0;
    if ((oDest.m_iWidth != oImage.m_iWidth) ||
        (oDest.m_iHeight != oImage.m_iHeight) ||
        (!oDest.m_pChannel[iComp]))
    {
        printf("Invalid image parameter in bLoad_Comp_GPU\n");
        return 0;
    }
    cudaMemcpy(oDest.m_pChannel[iComp], oImage.m_pChannel[0], oImage.m_iWidth * oImage.m_iHeight, cudaMemcpyHostToDevice);
    Disp_Cuda_Error();
    Free_Image(&oImage);
    return 0;
}
int bLoad_Image_GPU(const char* pcFile, Image* poImage, int iWidth, int iHeight, int iFrame)
{//简化操作，如果poImage->m_pBuffer已经有东西，则不用开内存
	Image oImage;
	if (!bLoad_Image(pcFile, &oImage))
		return 0;

	if (!poImage->m_pBuffer)
		Init_Image_GPU(poImage, oImage.m_iWidth, oImage.m_iHeight, Image::IMAGE_TYPE_BMP, oImage.m_iBit_Count);
	else if (oImage.m_iWidth - poImage->m_iWidth +
		oImage.m_iHeight - poImage->m_iHeight /*+
		oImage.m_iChannel_Count - poImage->m_iChannel_Count*/
        )
	{
		printf("Invalid image parameter in bLoad_Image_GPU\n");
		return 0;
	}
	Copy_Image_To_GPU(oImage, *poImage);
	Disp_Cuda_Error();
    Free_Image(&oImage);
	return 1;
}

int bSave_Image_GPU(const char* pcFile, Image oImage)
{//头在内存，数据在显存
	Image oDest;
	Init_Image(&oDest, oImage.m_iWidth, oImage.m_iHeight, Image::IMAGE_TYPE_BMP, oImage.m_iBit_Count);
	Copy_Image_To_CPU(oImage, oDest);
    Disp_Cuda_Error();
	int bResult=bSave_Image(pcFile, oDest);
	Free_Image(&oDest);
	return bResult;
}

int bSave_Image_GPU(const char* pcFile, Image *poHeader_GPU)
{//连头都在GPU
    Image oImage;
    cudaMemcpy(&oImage, poHeader_GPU, sizeof(Image), cudaMemcpyDeviceToHost);
    return bSave_Image_GPU(pcFile, oImage);
}
int bSave_Comp_GPU(const char* pcFile, Image * poHeader_GPU, int iComp)
{
    Image oImage;
    cudaMemcpy(&oImage, poHeader_GPU, sizeof(Image), cudaMemcpyDeviceToHost);
    return bSave_Comp_GPU(pcFile, oImage,iComp);
}
int bSave_Comp_GPU(const char* pcFile, Image oImage,int iComp)
{
    Image oImage_1;
    Attach_Buffer(&oImage_1, oImage.m_pChannel[iComp], oImage.m_iWidth, oImage.m_iHeight, 1, Image::IMAGE_TYPE_BMP);
    return bSave_Image_GPU(pcFile, oImage_1);
}


//__global__ void Bi_Linear_cv_Reflect_GPU(Image Source[],Image Dest[], 
//    int iCount, float f_x, float f_y)
//{
//    //threadIdx: threadId in a block
//    //blockDim.x = Thread_Per_Block
//    //gridDim.x = Size/blockDim.x
//    //gridDim.y = Channel_Count
//    //gridDim.z = Image_Count
//    //blockIdx.x = Block_ID 
//    //blockIdx.y = Channel in Image
//    //blockIdx.z = Image ID in Group
//
//    int iThread_ID = blockIdx.x * blockDim.x + threadIdx.x;
//    //if (iThread_ID == 0)
//    //{
//    //    //printf("%d %d %d\n", blockDim.x, blockDim.y, blockDim.z);
//    //    //printf("%d %d %d\n", gridDim.x, gridDim.y, gridDim.z);
//    //    printf("blockIdx: %d %d %d\n", blockIdx.x, blockIdx.y, blockIdx.z);
//    //}
//    
//    unsigned short w_d = Dest[blockIdx.z].m_iWidth,
//        h_d = Dest[blockIdx.z].m_iHeight,
//        w_s = Source[blockIdx.z].m_iWidth,
//        h_s = Source[blockIdx.z].m_iHeight;
//
//    if (iThread_ID > w_d * h_d)
//    {
//        //printf("exit");
//        return;
//    }
//      
//    int x_d = iThread_ID % w_d,
//        y_d = iThread_ID / w_d;
//        
//    float x_s_f = (x_d + 0.5f) * f_x - 0.5f;
//    float y_s_f = (y_d + 0.5f) * f_y - 0.5f,
//
//    int y_s_0 = (int)floor(y_s_f);
//    int y_s_1 = y_s_0 + 1;
//    float w2 = (y_s_1 - y_s_f), w3 = 1.f - w2;
//
//    unsigned char* pCur_Line = &Source[blockIdx.z].m_pChannel[blockIdx.y][iGet_Border_y(y_s_0, h_s, BORDER_REFLECT) * w_s];
//    unsigned char* pNext_Line = &Source[blockIdx.z].m_pChannel[blockIdx.y][iGet_Border_y(y_s_1, h_s, BORDER_REFLECT) * w_s];
//
//    int x_s_0 = (int)floor(x_s_f);
//    int x_s_1 = x_s_0 + 1;
//    float w0 = (x_s_1 - x_s_f), w1 = 1.f - w0;
//
//    int x_s_0_r = iGet_Border_x(x_s_0, w_s, BORDER_REFLECT),
//        x_s_1_r = iGet_Border_x(x_s_1, w_s, BORDER_REFLECT);
//
//    float fValue_0, fValue_1;
//    fValue_0 = w0 * pCur_Line[x_s_0_r] + w1 * pCur_Line[x_s_1_r];
//    fValue_1 = w0 * pNext_Line[x_s_0_r] + w1 * pNext_Line[x_s_1_r];
//
//    Dest[blockIdx.z].m_pChannel[blockIdx.y][iThread_ID] = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5);
//    
//    return;
//}

__global__ void Bi_Linear_cv_Reflect_GPU(Image Source[], Image Dest[],
    int iCount, float f_x, float f_y,
    short w_s, short h_s, short w_d, short h_d)
{
    //threadIdx: threadId in a block
    //blockDim.x = Thread_Per_Block
    //gridDim.x = Size/blockDim.x
    //gridDim.y = Channel_Count
    //gridDim.z = Image_Count
    //blockIdx.x = Block_ID 
    //blockIdx.y = Channel in Image
    //blockIdx.z = Image ID in Group

    int iThread_ID = blockIdx.x * blockDim.x + threadIdx.x;
    /*unsigned short w_d = Dest[blockIdx.z].m_iWidth,
        h_d = Dest[blockIdx.z].m_iHeight,
        w_s = Source[blockIdx.z].m_iWidth,
        h_s = Source[blockIdx.z].m_iHeight;*/

    if (iThread_ID >= w_d * h_d)
        return;

    int x_d = iThread_ID % w_d,
        y_d = iThread_ID / w_d;

    float w2, w3, x_s_f = (x_d + 0.5f) * f_x - 0.5f;

    /*__shared__ unsigned char* pSource, * pDest;
    if (threadIdx.x == 0)
    {
        pSource = Source[blockIdx.z].m_pChannel[blockIdx.y];
        pDest = Dest[blockIdx.z].m_pChannel[blockIdx.y];
    }*/

    __syncthreads();

    unsigned char* pCur_Line, * pNext_Line;
    {
        float y_s_f = (y_d + 0.5f) * f_y - 0.5f;
        int y_s_0 = (int)floor(y_s_f);
        int y_s_1 = y_s_0 + 1;
        w2 = (y_s_1 - y_s_f), w3 = 1.f - w2;
        //pCur_Line = &pSource[iGet_Border_y_GPU(y_s_0, h_s, BORDER_REFLECT) * w_s];
        //pNext_Line = &pSource[iGet_Border_y_GPU(y_s_1, h_s, BORDER_REFLECT) * w_s];
        pCur_Line = &Source[blockIdx.z].m_pChannel[blockIdx.y][iGet_Border_y_GPU(y_s_0, h_s, BORDER_REFLECT) * w_s];
        pNext_Line = &Source[blockIdx.z].m_pChannel[blockIdx.y][iGet_Border_y_GPU(y_s_1, h_s, BORDER_REFLECT) * w_s];
    }

    float w0, w1;
    int x_s_0_r, x_s_1_r;
    {
        int x_s_0 = (int)floor(x_s_f);
        int x_s_1 = x_s_0 + 1;
        w0 = (x_s_1 - x_s_f), w1 = 1.f - w0;
        x_s_0_r = iGet_Border_x_GPU(x_s_0, w_s, BORDER_REFLECT);
        x_s_1_r = iGet_Border_x_GPU(x_s_1, w_s, BORDER_REFLECT);
    }

    float fValue_0, fValue_1;
    fValue_0 = w0 * pCur_Line[x_s_0_r] + w1 * pCur_Line[x_s_1_r];
    fValue_1 = w0 * pNext_Line[x_s_0_r] + w1 * pNext_Line[x_s_1_r];

    Dest[blockIdx.z].m_pChannel[blockIdx.y][iThread_ID] = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5);
    //pDest[iThread_ID] = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5);
    return;
}

void Bi_Linear_cv_GPU(Image Source[], Image Dest[], int iCount, 
    int w_s, int h_s, int w_d, int h_d, int iChannel,
    float fScale_x, float fScale_y, Border_Type iBorder_Type)
{//将opencv的双线性插值写成GPU版本
//无可置疑这个接口很垃圾，在探索怎样才能快点
    int iSize;  // = iCount * 2 * sizeof(Image);
    //先对Seam_est进行Bi_Linear
    int iThread_Per_Block = 256;
    //iSize = Dest[0].m_iWidth * Dest[0].m_iHeight;
    iSize = w_d * h_d;
    dim3 oGrid;
    oGrid.x = (iSize + iThread_Per_Block - 1) / iThread_Per_Block;
    oGrid.y = iChannel;
    oGrid.z = iCount;

    //Disp_Cuda_Error();
    Bi_Linear_cv_Reflect_GPU << <oGrid, iThread_Per_Block >> > (Source, Dest, iCount, 1.f / fScale_x, 1.f / fScale_y, w_s, h_s, w_d, h_d);
    //Disp_Cuda_Error();
    return;
}

__global__ void _Set_Color_GPU(Image oBitMap, unsigned int C0, unsigned int C1, unsigned int C2)
{//这个比异步拷贝cudaMemcpyAsync快很多
    typedef struct Pixel_4 {
        char Data[4];
    }Pixel_4;

    int iThread_ID = GET_THREAD_ID() << 2;
    union {
        Pixel_4 oValue;
        unsigned int iValue;
    };
    if (iThread_ID + 4 <= oBitMap.m_iWidth * oBitMap.m_iHeight)
    {
        iValue = C0;
        *(Pixel_4*)&oBitMap.m_pChannel[0][iThread_ID] = oValue;
        if(oBitMap.m_iChannel_Count==3)
        {
            iValue = C1;
            *(Pixel_4*)&oBitMap.m_pChannel[1][iThread_ID] = oValue;
            iValue = C2;
            *(Pixel_4*)&oBitMap.m_pChannel[2][iThread_ID] = oValue;
        }
    }else
    {
        while (iThread_ID < oBitMap.m_iWidth * oBitMap.m_iHeight)
        {
            oBitMap.m_pChannel[0][iThread_ID] = C0 & 0xFF;
            if (oBitMap.m_iChannel_Count == 3)
            {
                oBitMap.m_pChannel[1][iThread_ID] = C1 & 0xFF;
                oBitMap.m_pChannel[2][iThread_ID] = C2 & 0xFF;
            }
            iThread_ID++;
        }
    }
}

void Set_Color_GPU(Image oImage, int R, int G, int B, unsigned long long iStream)
{
    int iThread_Per_Block = 1024, iSize = oImage.m_iWidth * oImage.m_iHeight,
        iBlock_Count;
    int Color[3];
    if (oImage.m_iImage_Type == Image::IMAGE_TYPE_BMP)
    {
        Color[0] = (R << 24) + (R << 16) + (R << 8) + R;
        Color[1] = (G << 24) + (G << 16) + (G << 8) + G;
        Color[2] = (B << 24) + (B << 16) + (B << 8) + B;
    }
    else
    {
        _RGB_2_YUV(R, G, B, Color[0], Color[1], Color[2]);
        Color[0] = (Color[0] << 24) + (Color[0] << 16) + (Color[0] << 8) + Color[0];
        Color[1] = (Color[1] << 24) + (Color[1] << 16) + (Color[1] << 8) + Color[1];
        Color[2] = (Color[2] << 24) + (Color[2] << 16) + (Color[2] << 8) + Color[2];
    }

    iBlock_Count = (((iSize + 3) >> 2) + iThread_Per_Block - 1) / iThread_Per_Block;
    _Set_Color_GPU << <iBlock_Count, iThread_Per_Block >> > (oImage, Color[0], Color[1], Color[2]);
    return;
}

__global__ void _Pyr_Down_col_GPU(Data_Block<unsigned short*, 3>oMid, int iMid_Height, Image oDest, Border_Type iBorder_Type = BORDER_REFLECT101)
{//列方向
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oDest.m_iHeight * oDest.m_iWidth)
        return;
    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;

    short y1 = y << 1;
    int iMid_Pos = y1 * oDest.m_iWidth;
    int iMid_Size = iMid_Height * oDest.m_iWidth;
    unsigned short* pSource = &oMid.Data[0][x];
    unsigned int Mid_Pos[4] = { (unsigned int) iGet_Border_y_GPU(y1 - 2, iMid_Height, iBorder_Type) * oDest.m_iWidth,
        (unsigned int)iGet_Border_y_GPU(y1 - 1, iMid_Height, iBorder_Type)* oDest.m_iWidth,
        (unsigned int)iGet_Border_y_GPU(y1 + 1, iMid_Height, iBorder_Type)* oDest.m_iWidth,
        (unsigned int)iGet_Border_y_GPU(y1 + 2, iMid_Height, iBorder_Type)* oDest.m_iWidth };

    //for(int i=0;i<oDest.m_iChannel_Count;i++, pSource += iMid_Size)
    //{
        unsigned short iValue = pSource[iMid_Pos] * 6 +
            ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) << 2) +
            pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];
        oDest.m_pChannel[0][iThread_ID] = (iValue + 128) >> 8;
    //}
        if(oDest.m_iChannel_Count >1)
        {
            pSource += iMid_Size;
            iValue = pSource[iMid_Pos] * 6 +
                ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) << 2) +
                pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];
            oDest.m_pChannel[1][iThread_ID] = (iValue + 128) >> 8;

            if (oDest.m_iChannel_Count > 2)
            {
                pSource += iMid_Size;
                    iValue = pSource[iMid_Pos] * 6 +
                    ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) << 2) +
                    pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];
                oDest.m_pChannel[2][iThread_ID] = (iValue + 128) >> 8;
            }
        }

    /*for (int i = 0; i < oDest.m_iChannel_Count; i++)
    {
        unsigned short* pSource = &oMid.Data[i][x];
        unsigned short iValue = pSource[y1 * oDest.m_iWidth] * 6 +
            (pSource[iGet_Border_y_GPU(y1 - 1, iMid_Height, iBorder_Type) * oDest.m_iWidth] + pSource[iGet_Border_y_GPU(y1 + 1, iMid_Height, iBorder_Type) * oDest.m_iWidth]) * 4 +
            pSource[iGet_Border_y_GPU(y1 - 2, iMid_Height, iBorder_Type) * oDest.m_iWidth] + pSource[iGet_Border_y_GPU(y1 + 2, iMid_Height, iBorder_Type) * oDest.m_iWidth];
        oDest.m_pChannel[i][iThread_ID] = (iValue + 128) >> 8;
    }*/
    return;
}

__global__ void _Pyr_Down_row_Ref_GPU(Image oSource, int iMid_Width, Data_Block<unsigned short*, 3>oMid, Border_Type iBorder_Type = BORDER_REFLECT101)
{//行处理，这个更快
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oSource.m_iHeight * iMid_Width)
        return;

    short x = iThread_ID % iMid_Width,
        y = iThread_ID / iMid_Width;
    short x1 = x << 1;
    int iSize_s = oSource.m_iWidth * oSource.m_iHeight;
    unsigned char* pSource = &oSource.m_pChannel[0][y * oSource.m_iWidth];
    unsigned short Source_Pos[4] = { (unsigned short)iGet_Border_x_GPU(x1 - 2, oSource.m_iWidth, iBorder_Type),
        (unsigned short)iGet_Border_x_GPU(x1 - 1, oSource.m_iWidth, iBorder_Type),
       (unsigned short)iGet_Border_x_GPU(x1 + 1, oSource.m_iWidth, iBorder_Type),
       (unsigned short)iGet_Border_x_GPU(x1 + 2, oSource.m_iWidth, iBorder_Type) };

    oMid.Data[0][iThread_ID] = pSource[x1] * 6 +
        ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) << 2) +
        pSource[Source_Pos[0]] + pSource[Source_Pos[3]];
    pSource += iSize_s;

    if (oSource.m_iChannel_Count > 1)
    {
        oMid.Data[1][iThread_ID] = pSource[x1] * 6 +
            ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) << 2) +
            pSource[Source_Pos[0]] + pSource[Source_Pos[3]];
        pSource += iSize_s;

        if (oSource.m_iChannel_Count > 2)
            oMid.Data[2][iThread_ID] = pSource[x1] * 6 +
            ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) << 2) +
            pSource[Source_Pos[0]] + pSource[Source_Pos[3]];
    }

    return;
}

void Pyr_Down_GPU(Image oSource, Image oDest,unsigned short *pAux)
{//拉普拉斯金字塔下采样
    Data_Block<unsigned short*, 3>oMid;
    int iSize = oSource.m_iHeight * oDest.m_iWidth;

    if (pAux)
        oMid.Data[0] = pAux;
    else
        oMid.Data[0] = (unsigned short*)pMalloc_GPU(iSize * oSource.m_iChannel_Count * sizeof(unsigned short));
    
    oMid.Data[1] = oMid.Data[0] + iSize;
    oMid.Data[2] = oMid.Data[1] + iSize;

    dim3 oThread, oGrid;
    //先搞行方向
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    //oGrid.y = oDest.m_iChannel_Count;
    //Disp_GPU(oSource.m_pChannel[0], oSource.m_iHeight, oSource.m_iWidth, "Source");
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    _Pyr_Down_row_Ref_GPU << <oGrid, oThread >> > (oSource, oDest.m_iWidth, oMid);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);
    //Disp_GPU(oMid.Data[0], 1, 4);

    //再到列方向
    iSize = oDest.m_iWidth * oDest.m_iHeight;
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    _Pyr_Down_col_GPU << <oGrid, oThread >> > (oMid, oSource.m_iHeight, oDest);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //Disp_GPU(oDest.m_pChannel[0], oDest.m_iHeight, oDest.m_iWidth, "Dest");
    ////以下Commit似乎不能省，否则后续会因为Mid地提早释放出错，待测
    ////Disp_Cuda_Error();
    //bSave_Image_GPU("c:\\tmp\\4.bmp", oDest);
    //Compare_Image("c:\\tmp\\3.bmp", "c:\\tmp\\4.bmp");
    if (!pAux && oMid.Data[0])
        Free_GPU(oMid.Data[0]);
    return;
}
__global__ void _Pyr_Up_row_GPU(Image oSource, int iMid_Width, Data_Block<unsigned short*, 3> oMid)
{
    int iThread_ID = GET_THREAD_ID();
    int iSize_s = oSource.m_iWidth * oSource.m_iHeight;
    if (iThread_ID >= iSize_s)
        return;

    short x = iThread_ID % oSource.m_iWidth,
        y = iThread_ID / oSource.m_iWidth;
    unsigned char bHas_Remain_x = iMid_Width > (oSource.m_iWidth << 1) && (x == oSource.m_iWidth - 1) ? 1 : 0;
    unsigned char bIs_Source_Border = (x == oSource.m_iWidth - 1);
    unsigned char bEven = bIs_Source_Border && ((x << 1) + 1 < iMid_Width);

    int iPos_m = y * iMid_Width,
        iPos_s = y * oSource.m_iWidth;
    int iSize_m = iMid_Width * oSource.m_iHeight;
    
    unsigned short* pMid = &oMid.Data[0][iPos_m + (x << 1)];
    for (short i = 0; i < oSource.m_iChannel_Count; i++,pMid+=iSize_m)
    {
        //优化方法，快不少
        
        unsigned short Mid[3];
        //中间点
        unsigned char iPix = oSource.m_pChannel[i][iPos_s + x];
        if (x == 0)
        {//左边两点
            Mid[0] = iPix * 6 + (oSource.m_pChannel[i][iPos_s + 1] << 1);
            Mid[1] = (iPix + oSource.m_pChannel[i][iPos_s + 1]) << 2;
        }
        else if (bIs_Source_Border)
        {//右边两点，2对齐情况下
            Mid[0] = oSource.m_pChannel[i][iPos_s + x - 1] + iPix * 7;
            if (bEven)  //源行最后一点乘以二，对应目标也有两点
                Mid[1] = iPix << 3;
        }
        else
        {//一般情况
            Mid[0] = oSource.m_pChannel[i][iPos_s + x - 1] + iPix * 6 + oSource.m_pChannel[i][iPos_s + x + 1];
            Mid[1] = (iPix + oSource.m_pChannel[i][iPos_s + x + 1]) << 2;
        }

        pMid[0] = Mid[0];
        if (!bIs_Source_Border || bEven)
        {
            pMid[1] = Mid[1];
            if (bHas_Remain_x)  //后改由余数
                pMid[2] = Mid[1];
        }

        //if (x == 0)
        //{//左边两点
        //    oMid.Data[i][iPos_m] = oSource.m_pChannel[i][iPos_s] * 6 + oSource.m_pChannel[i][iPos_s + 1] * 2;
        //    oMid.Data[i][iPos_m + 1] = (oSource.m_pChannel[i][iPos_s] + oSource.m_pChannel[i][iPos_s + 1]) * 4;
        //}else if (x == oSource.m_iWidth - 1)
        //{//右边两点
        //    oMid.Data[i][iPos_m + (x << 1)] = oSource.m_pChannel[i][iPos_s + x - 1] + oSource.m_pChannel[i][iPos_s + x] * 7;
        //    if ((x << 1) + 1 < iMid_Width)
        //    {
        //        oMid.Data[i][iPos_m + (x << 1) + 1] = oSource.m_pChannel[i][iPos_s + x] * 8;
        //        if (bHas_Remain_x)  //后改由余数
        //            oMid.Data[i][iPos_m + (x << 1) + 2] = oMid.Data[i][iPos_m + (x << 1) + 1];
        //    }
        //}else
        //{//中间点
        //    oMid.Data[i][iPos_m + (x << 1)] = oSource.m_pChannel[i][iPos_s + x - 1] + oSource.m_pChannel[i][iPos_s + x] * 6 + oSource.m_pChannel[i][iPos_s + x + 1];
        //    oMid.Data[i][iPos_m + (x << 1) + 1] = (oSource.m_pChannel[i][iPos_s + x] + oSource.m_pChannel[i][iPos_s + x + 1]) * 4;
        //}
    }
}
__global__ void _Pyr_Up_col_Subtract_GPU(Data_Block<unsigned short*, 3>oMid, int iMid_Height, Image oDest)
{//行操作以后再加上减法操作
    int iThread_ID = GET_THREAD_ID();
    int iSize = iMid_Height * oDest.m_iWidth;

    if (iThread_ID >= iSize)
        return;

    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;
    unsigned char bHas_Remain_y = oDest.m_iHeight > iMid_Height * 2 && y == iMid_Height - 1 ? 1 : 0;
    unsigned char bEven = (y * 2 + 1 < oDest.m_iHeight);

    //重复计算的全部提取出来，快不少
    unsigned short* r0 = &oMid.Data[0][iGet_Border_y_GPU(y - 1, iMid_Height, BORDER_REFLECT101) * oDest.m_iWidth + x],
        * r1 = &oMid.Data[0][iGet_Border_y_GPU(y, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x],
        * r2 = &oMid.Data[0][iGet_Border_y_GPU(y + 1, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x];

    //int iPos_d = (y * 2) * oDest.m_iWidth + x;
    //for (int i = 0; i < oDest.m_iChannel_Count; i++, r0 += iSize, r1 += iSize, r2 += iSize)
    //{
    //    short iValue = oDest.m_pChannel[i][iPos_d] - ((*r0 + *r1 * 6 + *r2 + 32) >> 6);
    //    oDest.m_pChannel[i][iPos_d] = Clip3(-128, 127, iValue);
    //    if (bEven)
    //    {
    //        //d1行 =   [(r1 + r2)*4 + 32]>>6
    //        iValue = oDest.m_pChannel[i][iPos_d + oDest.m_iWidth] - ((((*r1 + *r2) << 2) + 32) >> 6);
    //        oDest.m_pChannel[i][iPos_d + oDest.m_iWidth] = Clip3(-128, 127, iValue);
    //        if (bHas_Remain_y)
    //        {
    //            iValue = oDest.m_pChannel[i][iPos_d + (oDest.m_iWidth << 1)] - oDest.m_pChannel[i][iPos_d];
    //            oDest.m_pChannel[i][iPos_d + (oDest.m_iWidth << 1)] = Clip3(-128, 127, iValue);
    //        }
    //    }
    //}

    unsigned char* pDest = &oDest.m_pChannel[0][(y * 2) * oDest.m_iWidth + x];
    int iDest_Size = oDest.m_iWidth * oDest.m_iHeight;
    int iValue = *pDest - ((*r0 + *r1 * 6 + *r2 + 32) >> 6);
    *pDest = Clip3(-128, 127, iValue);
    if (bEven)
    {
        iValue = pDest[oDest.m_iWidth] - ((((*r1 + *r2) << 2) + 32) >> 6);
        pDest[oDest.m_iWidth] = Clip3(-128, 127, iValue);
        if (bHas_Remain_y)
        {
            iValue = pDest[oDest.m_iWidth << 1] - (*pDest);
            pDest[oDest.m_iWidth << 1] = Clip3(-128, 127, iValue);
        }
    }

    if (oDest.m_iChannel_Count > 1)
    {
        pDest += iDest_Size, r0 += iSize, r1 += iSize, r2 += iSize;
        iValue = *pDest - ((*r0 + *r1 * 6 + *r2 + 32) >> 6);
        *pDest = Clip3(-128, 127, iValue);
        if (bEven)
        {
            iValue = pDest[oDest.m_iWidth] - ((((*r1 + *r2) << 2) + 32) >> 6);
            pDest[oDest.m_iWidth] = Clip3(-128, 127, iValue);
            if (bHas_Remain_y)
            {
                iValue = pDest[oDest.m_iWidth << 1] - (*pDest);
                pDest[oDest.m_iWidth << 1] = Clip3(-128, 127, iValue);
            }
        }

        if (oDest.m_iChannel_Count > 2)
        {
            pDest += iDest_Size, r0 += iSize, r1 += iSize, r2 += iSize;
            iValue = *pDest - ((*r0 + *r1 * 6 + *r2 + 32) >> 6);
            *pDest = Clip3(-128, 127, iValue);
            if (bEven)
            {
                iValue = pDest[oDest.m_iWidth] - ((((*r1 + *r2) << 2) + 32) >> 6);
                pDest[oDest.m_iWidth] = Clip3(-128, 127, iValue);
                if (bHas_Remain_y)
                {
                    iValue = pDest[oDest.m_iWidth << 1] - (*pDest);
                    pDest[oDest.m_iWidth << 1] = Clip3(-128, 127, iValue);
                }
            }
        }
    }

}
__global__ void _Pyr_Up_col_GPU(Data_Block<unsigned short*, 3>oMid, int iMid_Height, Image oDest)
{
    int iThread_ID = GET_THREAD_ID();
    int iSize = iMid_Height * oDest.m_iWidth;

    if (iThread_ID >= iSize)
        return;

    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;
    unsigned char bHas_Remain_y = oDest.m_iHeight > iMid_Height * 2 && y == iMid_Height - 1 ? 1 : 0;
    unsigned char bEven = (y * 2 + 1 < oDest.m_iHeight);

    //重复计算的全部提取出来，快不少
    unsigned short* r0 = &oMid.Data[0][iGet_Border_y_GPU(y - 1, iMid_Height, BORDER_REFLECT101) * oDest.m_iWidth + x],
        * r1 = &oMid.Data[0][iGet_Border_y_GPU(y, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x],
        * r2 = &oMid.Data[0][iGet_Border_y_GPU(y + 1, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x];
    
    //int iPos_d = (y * 2) * oDest.m_iWidth + x;
    //for (int i = 0; i < oDest.m_iChannel_Count; i++, r0 += iSize, r1 += iSize, r2 += iSize)
    //{
    //    oDest.m_pChannel[i][iPos_d] = (*r0 + *r1 * 6 + *r2 + 32) >> 6;
    //    if (bEven)
    //    {
    //        //d1行 =   [(r1 + r2)*4 + 32]>>6
    //        oDest.m_pChannel[i][iPos_d + oDest.m_iWidth] = (((*r1 + *r2) << 2) + 32) >> 6;
    //        if (bHas_Remain_y)
    //            oDest.m_pChannel[i][iPos_d + (oDest.m_iWidth << 1)] = oDest.m_pChannel[i][iPos_d];
    //    }
    //}
    unsigned char* pDest = &oDest.m_pChannel[0][(y * 2) * oDest.m_iWidth + x];
    int iDest_Size = oDest.m_iWidth * oDest.m_iHeight;
    *pDest = (*r0 + *r1 * 6 + *r2 + 32) >> 6;
    if (bEven)
    {
        pDest[oDest.m_iWidth] = (((*r1 + *r2) << 2) + 32) >> 6;
        if (bHas_Remain_y)
            pDest[oDest.m_iWidth << 1] = *pDest;
    }

    if(oDest.m_iChannel_Count>1)
    {
        pDest += iDest_Size, r0 += iSize, r1 += iSize, r2 += iSize;
        *pDest = (*r0 + *r1 * 6 + *r2 + 32) >> 6;
        if (bEven)
        {
            pDest[oDest.m_iWidth] = (((*r1 + *r2) << 2) + 32) >> 6;
            if (bHas_Remain_y)
                pDest[oDest.m_iWidth << 1] = *pDest;
        }

        if (oDest.m_iChannel_Count > 2)
        {
            pDest += iDest_Size, r0 += iSize, r1 += iSize, r2 += iSize;
            *pDest = (*r0 + *r1 * 6 + *r2 + 32) >> 6;
            if (bEven)
            {
                pDest[oDest.m_iWidth] = (((*r1 + *r2) << 2) + 32) >> 6;
                if (bHas_Remain_y)
                    pDest[oDest.m_iWidth << 1] = *pDest;
            }
        }
    }
}

void Pyr_Up_GPU(Image oSource, Image oDest)
{
    Data_Block<unsigned short*, 3>oMid;
    int iSize = oSource.m_iHeight * oDest.m_iWidth;
    oMid.Data[0] = (unsigned short*)pMalloc_GPU(iSize * oDest.m_iChannel_Count * sizeof(unsigned short));
    oMid.Data[1] = oMid.Data[0] + iSize;
    oMid.Data[2] = oMid.Data[1] + iSize;

    dim3 oThread, oGrid;
    iSize = oSource.m_iWidth * oSource.m_iHeight;
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;

    /*Disp_Cuda_Error();
    unsigned long long tStart = iGet_Tick_Count();
    for(int i=0;i<10000;i++)*/
    _Pyr_Up_row_GPU << <oGrid, oThread >> > (oSource, oDest.m_iWidth, oMid);
    /*Disp_Cuda_Error();
    printf("%lld\n", iGet_Tick_Count() - tStart);*/

    //Disp_GPU(oMid.Data[0], oSource.m_iHeight, oDest.m_iWidth, "Mid");
    iSize = oDest.m_iWidth * oSource.m_iHeight;
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    //_Pyr_Up_col_GPU << <oGrid, oThread >> > (oMid, oSource.m_iHeight, oDest);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    _Pyr_Up_col_Subtract_GPU << <oGrid, oThread >> > (oMid, oSource.m_iHeight, oDest);

    ////Disp_GPU(oDest.m_pChannel[0], oDest.m_iHeight, oDest.m_iWidth, "Dest");
    //bSave_Image_GPU("c:\\tmp\\3.bmp", oDest);
    //Compare_Image("c:\\tmp\\2.bmp", "c:\\tmp\\3.bmp");
    if (oMid.Data[0])
        Free_GPU(oMid.Data[0]);
    return;
}

__global__ void Gen_Gauss_Filter_GPU(int r, float fSigma, float Kernel[])
{
    int i, j, d = r + r + 1;
    float* pFilter = Kernel;
    float fValue, fSum = 0;

    for (i = 0; i < r; i++)
    {
        fValue = (float)(i - r) / fSigma;
        pFilter[i] = (float)exp(-0.5f * fValue * fValue);
        fSum += pFilter[i];
    }

    fValue = (float)(i - r) / fSigma;
    pFilter[i] = (float)exp(-0.5f * fValue * fValue);
    fSum = fSum * 2 + pFilter[i];

    fValue = 1.f / fSum;
    j = d - 1;
    for (i = 0; i < r; i++, j--)
        pFilter[j] = (pFilter[i] *= fValue);
    pFilter[r] *= fValue;
    return;
}
__global__ void _Copy_Make_Border_GPU(Image oSource, Image oDest, short iLeft, short iTop, short iRight, short iBottom, char iPix_Group_Size)
{//有点乱，需要化简

    int iDest_Size = oDest.m_iWidth * oDest.m_iHeight;
    short x_s = threadIdx.x * iPix_Group_Size;
    short x_d = x_s + iLeft,
        y_d = blockIdx.x + iTop;

    {//先抄中间部分
        //int iPos_s = blockIdx.x * oSource.m_iWidth + x_s;
        //int iPos_d = y_d * oDest.m_iWidth + x_d;
        int iSource_Size = oSource.m_iWidth * oSource.m_iHeight;
        short iRemain_x = oSource.m_iWidth - x_s;
        unsigned char* pDest = &oDest.m_pChannel[0][y_d * oDest.m_iWidth + x_d],
            * pSource = &oSource.m_pChannel[0][blockIdx.x * oSource.m_iWidth + x_s];

        for (short j = 0; j < oDest.m_iChannel_Count; j++, pSource += iSource_Size, pDest += iDest_Size)
        {
            if (iRemain_x < iPix_Group_Size)
            {
                for (short i = 0; i < iRemain_x; i++)
                    pDest[i] = pSource[i];
            }
            else
                *(Pixel_4*)pDest = *(Pixel_4*)pSource;
        }

        //if (iRemain_x < iPix_Group_Size)
        //{//最后一个,并且不满，要搞特殊
        //    /*for (short i = 0; i < iRemain_x; i++)
        //        for (short j = 0; j < oDest.m_iChannel_Count; j++)
        //            oDest.m_pChannel[j][iPos_d + i] = oSource.m_pChannel[j][iPos_s + i];*/
        //    for (short j = 0; j < oDest.m_iChannel_Count; j++)
        //        for (short i = 0; i < iRemain_x; i++)
        //            oDest.m_pChannel[j][iPos_d + i] = oSource.m_pChannel[j][iPos_s + i];
        //}else
        //{
        //    for (short j = 0; j < oDest.m_iChannel_Count; j++)
        //        *(Pixel_4*)&oDest.m_pChannel[j][iPos_d] = *(Pixel_4*)&oSource.m_pChannel[j][iPos_s];
        //}
    }

    __syncthreads();

    //将中间部分抄到左边界
    if (threadIdx.x < iLeft)
    {
        /* unsigned char* pDest = &oDest.m_pChannel[0][y_d * oDest.m_iWidth + threadIdx.x],
             * pSource = &oDest.m_pChannel[0][y_d * oDest.m_iWidth + iLeft * 2 - threadIdx.x - 1];
         for (int j = 0; j < oDest.m_iChannel_Count; j++,pSource+= iDest_Size,pDest+= iDest_Size)
             *pDest = *pSource;*/
        for (int j = 0; j < oDest.m_iChannel_Count; j++)
            oDest.m_pChannel[j][y_d * oDest.m_iWidth + threadIdx.x] = oDest.m_pChannel[j][y_d * oDest.m_iWidth + iLeft * 2 - threadIdx.x - 1];
    }

    if (threadIdx.x < iRight)
    {
        unsigned char* pDest = &oDest.m_pChannel[0][y_d * oDest.m_iWidth + oDest.m_iWidth - threadIdx.x - 1],
            * pSource = &oDest.m_pChannel[0][y_d * oDest.m_iWidth + oDest.m_iWidth - iRight * 2 + threadIdx.x];
        for (int j = 0; j < oDest.m_iChannel_Count; j++, pSource += iDest_Size, pDest += iDest_Size)
            *pDest = *pSource;
        //for (int j = 0; j < oDest.m_iChannel_Count; j++)
            //oDest.m_pChannel[j][ y_d * oDest.m_iWidth + oDest.m_iWidth - threadIdx.x - 1] = oDest.m_pChannel[j][y_d * oDest.m_iWidth  + oDest.m_iWidth - iRight*2 + threadIdx.x];
    }
    __syncthreads();

    short iWidth_Div_4 = oDest.m_iWidth >> 2;
    //把中间的内容抄到上面去
    if (blockIdx.x < iTop)
    {
        unsigned char* pDest = &oDest.m_pChannel[0][blockIdx.x * oDest.m_iWidth],
            * pSource = &oDest.m_pChannel[0][((iTop << 1) - blockIdx.x - 1) * oDest.m_iWidth];
        int iDist = oDest.m_iWidth - (iWidth_Div_4 << 2);

        for (short j = 0; j < oDest.m_iChannel_Count; j++, pSource += iDest_Size, pDest += iDest_Size)
        {
            for (short x = threadIdx.x; x < iWidth_Div_4; x += blockDim.x)
                *(Pixel_4*)&pDest[x * 4] = *(Pixel_4*)&pSource[x * 4];

            //for (short x = threadIdx.x; x < iWidth_Div_4; x += blockDim.x)
            //    *(Pixel_4*)&oDest.m_pChannel[j][blockIdx.x * oDest.m_iWidth + x * 4] = *(Pixel_4*)&oDest.m_pChannel[j][(iTop * 2 - blockIdx.x - 1) * oDest.m_iWidth + x * 4];

            //行尾余数
            if (threadIdx.x < iDist)
                pDest[oDest.m_iWidth - threadIdx.x - 1] = pDest[oDest.m_iWidth - iRight * 2 + threadIdx.x];
            //oDest.m_pChannel[j][blockIdx.x * oDest.m_iWidth + oDest.m_iWidth - threadIdx.x - 1] = oDest.m_pChannel[j][blockIdx.x * oDest.m_iWidth + oDest.m_iWidth - iRight * 2 + threadIdx.x];
        }
    }

    //将中间内容抄到下边去
    if (blockIdx.x >= oSource.m_iHeight - iBottom)
    {
        int iDist_y = oSource.m_iHeight - 1 - blockIdx.x;
        int iDist = oDest.m_iWidth - (iWidth_Div_4 << 2);
        unsigned char* pDest = &oDest.m_pChannel[0][(oDest.m_iHeight - iBottom + iDist_y) * oDest.m_iWidth],
            * pSource = &oDest.m_pChannel[0][(oDest.m_iHeight - iBottom - iDist_y - 1) * oDest.m_iWidth];

        for (int j = 0; j < oDest.m_iChannel_Count; j++, pSource += iDest_Size, pDest += iDest_Size)
        {
            for (short x = threadIdx.x; x < iWidth_Div_4; x += blockDim.x)
            {
                *(Pixel_4*)&pDest[x * 4] = *(Pixel_4*)&pSource[x * 4];
                //*(Pixel_4*)&oDest.m_pChannel[j][(oDest.m_iHeight - iBottom + iDist_y) * oDest.m_iWidth + x * 4] =
                    //*(Pixel_4*)&oDest.m_pChannel[j][(oDest.m_iHeight - iBottom - iDist_y - 1) * oDest.m_iWidth + x * 4];
            }
            //行尾余数
            if (threadIdx.x < iDist)
                pDest[oDest.m_iWidth - threadIdx.x - 1] = pDest[oDest.m_iWidth - iRight * 2 + threadIdx.x];
            //oDest.m_pChannel[j][(oDest.m_iHeight - iBottom + iDist_y) * oDest.m_iWidth + oDest.m_iWidth - threadIdx.x - 1] = oDest.m_pChannel[j][(oDest.m_iHeight - iBottom + iDist_y) * oDest.m_iWidth + oDest.m_iWidth - iRight * 2 + threadIdx.x];
        }
    }
}

void Copy_Make_Border_GPU(Image oSource, Image oDest, short iLeft, short iTop, short iRight, short iBottom)
{//尝试优化扩边拷贝
    dim3 oThread, oGrid;
    const int iPix_Group_Size = 4;
    oThread.x = (oSource.m_iWidth + iPix_Group_Size - 1) / iPix_Group_Size;
    oGrid.x = oSource.m_iHeight;
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();

    //最快427，最慢480
    //for (int i = 0; i < 10000; i++)
        _Copy_Make_Border_GPU << <oGrid, oThread >> > (oSource, oDest, iLeft, iTop, iRight, iBottom, iPix_Group_Size);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);
    //bSave_Image_GPU("c:\\tmp\\3.bmp", oDest);
    //Compare_Image("c:\\tmp\\2.bmp", "c:\\tmp\\3.bmp");

    return;
}