#include <iostream>
#include "cuda_runtime.h"
#include "Common.h"
#include "Common_cuda.cuh"
#include "Matrix.h"
#include "image_cuda.h"
#include "Panorama.h"

void SB_Pano()
{
    Init_Stitch<double>(NULL, 0, 0);
    Free_Stitch<double>(NULL);
    Resize_Seam_Image<double>(NULL);
    Warp_2<double>({}, {}, NULL, NULL, 0, NULL, NULL, INTER_LINEAR, INTER_NEAREST, NULL);
    Warp_3<double>({}, NULL, NULL, 0, NULL, INTER_AREA, INTER_AREA, NULL, BORDER_CONSTANT, BORDER_CONSTANT);

    Feed<double>(NULL);
    Find<double>(NULL);
    Re_Map_3_GPU<double>({}, {}, 0, 0, NULL, INTER_LINEAR, INTER_LINEAR, BORDER_CONSTANT, BORDER_CONSTANT);
    Block_Compensate<double>(NULL, NULL, NULL);
    Resize_Bitwise_And<double>(NULL, NULL, NULL);

}

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

template<typename _T>void Free_Stitch(Stitch<_T>* poStitch)
{
    int i;
    if (poStitch->m_pCamera)
    {
        Free(poStitch->m_pCamera);
        poStitch->m_pCamera = NULL;
    }

    if (poStitch->m_pImage_Source)
    {
        //已经改进，一个指针搞定
        for(int i=0;i<poStitch->m_iImage_Count;i++)
            Free_Image_GPU(&poStitch->m_pImage_Source[0]);
        Free(poStitch->m_pImage_Source);
        poStitch->m_pImage_Source = NULL;
    }
    if (poStitch->m_pSeam_Est)
    {
        //for (i = 0; i < poStitch->m_iImage_Count; i++)
        Free_Image_GPU(&poStitch->m_pSeam_Est[0]);
        //Free_GPU(poStitch->m_pSeam_Est);
        poStitch->m_pSeam_Est = NULL;
    }
    /*if (poStitch->m_pMask)
    {
        for (i = 0; i < poStitch->m_iImage_Count; i++)
            Free_Image(&poStitch->m_pMask[i]);
        Free(poStitch->m_pMask);
        poStitch->m_pMask = NULL;
    }*/
    if (poStitch->m_pImage_Match)
    {
        for (i = 0; i < poStitch->m_iImage_Count; i++)
        {
            for (int j = 0; j < poStitch->m_iImage_Count; j++)
            {
                if (poStitch->m_pImage_Match[i * poStitch->m_iImage_Count + j].m_iMatch_Count)
                    Free(poStitch->m_pImage_Match[i & poStitch->m_iImage_Count].m_pPoint_1);
            }
        }
        Free(poStitch->m_pImage_Match);
        poStitch->m_pImage_Match = NULL;
    }

    if (poStitch->m_pImage_Warp)
    {
        for (i = 0; i < poStitch->m_iImage_Count; i++)
            Free_Image_GPU(&poStitch->m_pImage_Warp[i]);
        //Free(poStitch->m_pImage_Warped);
        poStitch->m_pImage_Warp = NULL;
    }

    if (poStitch->m_pMasks_Warped)
    {
        for (i = 0; i < poStitch->m_iImage_Count; i++)
            Free_Image_GPU(&poStitch->m_pMasks_Warped[i]);
        //Free(poStitch->m_pMasks_Warped);
        //poStitch->m_pMasks_Warped = NULL;
    }

    //if (poStitch->m_pBlock_Images_GPU)//每一个Image分出多少个Title
    //{
    //    for (i = 0; i < poStitch->m_iBlock_Count; i++)
    //        Free_Image(&poStitch->m_pBlock_images[i]);
    //    Free(poStitch->m_pBlock_images);
    //    poStitch->m_pBlock_images = NULL;
    //}

    if (poStitch->m_pBlock_Masks)//每个Image分出多少个mask
    {
        for (i = 0; i < poStitch->m_iBlock_Count; i++)
            Free_Image(&poStitch->m_pBlock_Masks[i]);
        Free(poStitch->m_pBlock_Masks);
        poStitch->m_pBlock_Masks = NULL;
    }

    /*if (poStitch->m_pBlock_Corners)
    {
        Free(poStitch->m_pBlock_corners);
        poStitch->m_pBlock_corners = NULL;
    }*/

    /*if (poStitch->m_pBlock_per_imgs)
    {
        Free(poStitch->m_pBlock_per_imgs);
        poStitch->m_pBlock_per_imgs = NULL;
    }*/
    /*if (poStitch->m_pCorner)
    {
        Free(poStitch->m_pCorner);
        poStitch->m_pCorner = NULL;
    }*/
    if (poStitch->m_pImage_Source_Header_GPU)
        Free_GPU(poStitch->m_pImage_Source_Header_GPU);
    if (poStitch->m_pBlock_Corner)
        Free(poStitch->m_pBlock_Corner);
    if (poStitch->m_pBlock_Corner_GPU)
        Free_GPU(poStitch->m_pBlock_Corner_GPU);
    if (poStitch->m_pBlock_Image)
        Free(poStitch->m_pBlock_Image);
    if (poStitch->m_pBlock_Image_Header_GPU)
        Free_GPU(poStitch->m_pBlock_Image_Header_GPU);
    if (poStitch->m_pBlock_Image_Data_GPU)
        Free_GPU(poStitch->m_pBlock_Image_Data_GPU);
    if (poStitch->m_oComp.m_pGain)
        Free_GPU(poStitch->m_oComp.m_pGain);

    return;
}

template<typename _T>void Set_Scale(Stitch<_T>* poStitch)
{//拼接器中有大量scale，能过一个算一个

    //关于配准？
    poStitch->registr_resol_ = 0.3;

    //设workscale，由registr_resol_决定
    //work_scale_ = std::min(1.0, std::sqrt(registr_resol_ * 1e6 / full_img_sizes_[i].area()));
    int iArea = poStitch->m_Source_Size[0] * poStitch->m_Source_Size[1];
    _T fValue = poStitch->registr_resol_ * 1e6 / (_T)iArea;
    fValue = sqrt(fValue);  
    poStitch->work_scale = Min(1.0, fValue);
    poStitch->m_bWork_Scale_Set = 1;

    //设seam_scale_/ seam_work_aspect_
    //seam_scale_ = std::min(1.0, std::sqrt(seam_est_resol_ * 1e6 / full_img_sizes_[i].area()));
    poStitch->seam_est_resol_ = 0.1;    //拍脑袋值
    fValue = poStitch->seam_est_resol_ * 1e6 / (_T)iArea;
    fValue = sqrt(fValue);
    poStitch->seam_scale = Min(1.0, fValue);
    poStitch->seam_work_aspect = poStitch->seam_scale / poStitch->work_scale;
    poStitch->m_bSeam_Scale_set = 1;
    poStitch->m_Seam_Size[0] = (short)(poStitch->m_Source_Size[0] * poStitch->seam_scale + 0.5f);
    poStitch->m_Seam_Size[1] = (short)(poStitch->m_Source_Size[1] * poStitch->seam_scale + 0.5f);

    //设置warped_image_scale
    //可见，warped_image_scale并非一个和其它相比的scale, 而是各个相机的焦距参数均值
        
    //此处有个奇怪的算法，暂不知何意，原来的猜想对不上cv。待验证
    //_T fTotal = 0;
    //for (int i = 0; i < poStitch->m_iImage_Count; i++)
    //    fTotal += poStitch->m_pCamera[i].K[0];
    //poStitch->warped_image_scale = fTotal/poStitch->m_iImage_Count;
    if (poStitch->m_iImage_Count & 1)
    {//偶数时，取中间帧的焦距为warped_image_scale
        poStitch->warped_image_scale = poStitch->m_pCamera[poStitch->m_iImage_Count >> 1].K[0];
    }else
    {
        int iHalf = poStitch->m_iImage_Count >> 1;
        _T fTotal = poStitch->m_pCamera[iHalf - 1].K[0] + poStitch->m_pCamera[iHalf].K[0];
        poStitch->warped_image_scale = fTotal / 2;
    }

    return;
}
template<typename _T>void Init_Stitch(Stitch<_T>* poStitch, int iWidth, int iHeight)
{//iWidth,iHeight为原图的分辨率
    int i, iSize;
    Stitch<_T> oSt = *poStitch;
    Light_Ptr oPtr_GPU;
    unsigned char* p;

    //统一设置各种scale
    oSt.m_Source_Size[0] = iWidth, oSt.m_Source_Size[1] = iHeight;
    Set_Scale(&oSt);
   
    //先搞CPU
    //光是放在显存的图像头
    int iSize_Header = oSt.m_iImage_Count * sizeof(Image) * (1 +        //Source
        1 +     //Seam_est
        1 +     //Mask
        1 +     //Image_Warped
        1) +       //Mask_Warped
        oSt.m_iImage_Count * sizeof(Image_Size_In_Block) +    //Block_Per_Image
        oSt.m_iImage_Count * 2 * 2 * sizeof(int) + //Corner
        oSt.m_iImage_Count * 2 * sizeof(int) ;      //Size;

    p = (unsigned char*)pMalloc(iSize_Header);
    memset(p, 0, iSize_Header);
    
    oSt.m_pImage_Source = (Image*)p;
    oSt.m_pSeam_Est = oSt.m_pImage_Source + oSt.m_iImage_Count;
    oSt.m_pMask = oSt.m_pSeam_Est + oSt.m_iImage_Count;
    oSt.m_pImage_Warp = oSt.m_pMask + oSt.m_iImage_Count;
    oSt.m_pMasks_Warped = oSt.m_pImage_Warp + oSt.m_iImage_Count;
    oSt.m_pBlock_Per_Image = (Image_Size_In_Block*)(oSt.m_pMasks_Warped + oSt.m_iImage_Count);
    oSt.m_pCorner = (int(*)[2][2])(oSt.m_pBlock_Per_Image + oSt.m_iImage_Count);
    oSt.m_pSize = (int(*)[2])(oSt.m_pCorner + oSt.m_iImage_Count);

    //可分离卷积核
    oSt.m_pKer_GPU = (_T*)(oSt.m_pCorner + oSt.m_iImage_Count);
    _T Ker[] = { 0.25, 0.5, 0.25,0,0 };
    cudaMemcpy(oSt.m_pKer_GPU, Ker, 3 * sizeof(_T), cudaMemcpyHostToDevice);

    iSize = iSize_Header +
        oSt.m_iImage_Count * (iWidth * iHeight * 4 + 128);                //Source Image

    Attach_Light_Ptr(oPtr_GPU, (unsigned char*)pMalloc_GPU(iSize), iSize, 0);
    p = (unsigned char*)pMalloc_GPU(iSize_Header);
    oSt.m_pImage_Source_Header_GPU = (Image*)p;
    oSt.m_pSeam_Est_Header_GPU = oSt.m_pImage_Source_Header_GPU + oSt.m_iImage_Count;
    oSt.m_pMask_Header_GPU = oSt.m_pSeam_Est_Header_GPU + oSt.m_iImage_Count;
    oSt.m_pImage_Warp_Header_GPU = oSt.m_pMask_Header_GPU + oSt.m_iImage_Count;
    oSt.m_pMask_Warp_Header_GPU = oSt.m_pImage_Warp_Header_GPU + oSt.m_iImage_Count;
    oSt.m_pBlock_Per_Image_GPU = (Image_Size_In_Block*)(oSt.m_pMask_Warp_Header_GPU + oSt.m_iImage_Count);
    oSt.m_pCorner_GPU =(int (*)[2][2]) oSt.m_pBlock_Per_Image_GPU + oSt.m_iImage_Count;
    for (i = 0; i < poStitch->m_iImage_Count; i++)
    {
        Init_Image_GPU(&oSt.m_pImage_Source[i], iWidth, iHeight, Image::IMAGE_TYPE_BMP, 32, &oPtr_GPU);
        cudaMemset(oSt.m_pImage_Source[i].m_pChannel[3], 255, iWidth * iHeight * 3);
    }

    //这些在后面会释放，所以另开一块显存 Seam_est + Seam_Mask + Image_Warp
    iSize = oSt.m_iImage_Count * (oSt.m_Seam_Size[0] * oSt.m_Seam_Size[1] * 3 + 128) +      //Seam_est
        oSt.m_iImage_Count * (oSt.m_Seam_Size[0] * oSt.m_Seam_Size[1] + 128) ;               //Seam_Mask
    Attach_Light_Ptr(oPtr_GPU, (unsigned char*)pMalloc_GPU(iSize), iSize, 0);
   
    //Seam_est + Seam_Mask + Image_Warp
    for (i = 0; i < poStitch->m_iImage_Count; i++)
        Init_Image_GPU(&oSt.m_pSeam_Est[i], oSt.m_Seam_Size[0], oSt.m_Seam_Size[1], Image::IMAGE_TYPE_BMP, 24, &oPtr_GPU);
    for (i = 0; i < poStitch->m_iImage_Count; i++)
    {
        Init_Image_GPU(&oSt.m_pMask[i], oSt.m_Seam_Size[0], oSt.m_Seam_Size[1], Image::IMAGE_TYPE_BMP, 8, &oPtr_GPU);
        Set_Color_GPU(oSt.m_pMask[i],255);
        //Disp_Cuda_Error();
    }

    //继续对其他可以通过scale算出来的Image开辟空间
    cudaMemcpy(p, oSt.m_pImage_Source, iSize_Header,cudaMemcpyHostToDevice);
        

    *poStitch = oSt;
    //Disp_Cuda_Error();
    return;
}
template<typename _T>void Map_Forward(float x, float y, float* pu, float* pv, Sphere_Projector<_T> oProjector, Point_Cloud<float>* poPC = NULL)
{
    static int iCount = 0;
    //给定的(x,y,1),求经过投影后的坐标
    //(x_,y_,z_)' = R * K(-1) * (x,y,1)'

    //投影？
    float x_ = (float)(oProjector.r_kinv[0] * x + oProjector.r_kinv[1] * y + oProjector.r_kinv[2]);
    float y_ = (float)(oProjector.r_kinv[3] * x + oProjector.r_kinv[4] * y + oProjector.r_kinv[5]);
    float z_ = (float)(oProjector.r_kinv[6] * x + oProjector.r_kinv[7] * y + oProjector.r_kinv[8]);

    if (poPC)
    {
        if (iCount < 1278)
            Draw_Point(poPC, x_, y_, z_, 255, 0, 0);
        else
            Draw_Point(poPC, x_, y_, z_, 0, 255, 0);
    }

    *pu = (float)(oProjector.scale * atan2f(x_, z_));
    float w = y_ / sqrtf(x_ * x_ + y_ * y_ + z_ * z_);
    *pv = (float)(oProjector.scale * (PI - acosf(w == w ? w : 0)));

    iCount++;
    return;
}
template<typename _T>void Detect_Result_Roi_By_Border(int w, int h, int Dest_tl[2], int Dest_br[2],
    Sphere_Projector<_T> oProjector, Point_Cloud<float>* poPC = NULL)
{//tl: Top left     br: Bottom Right
    static int iCount = 0;
    float tl_uf = (std::numeric_limits<float>::max)();
    float tl_vf = (std::numeric_limits<float>::max)();
    float br_uf = -(std::numeric_limits<float>::max)();
    float br_vf = -(std::numeric_limits<float>::max)();
    float u, v;

    //此处也许是扫过一段弧，看看最左在哪，最右又在哪的意思
    for (int x = 0; x < w; x++)
    {
        Map_Forward((float)x, 0, &u, &v, oProjector, poPC);
        tl_uf = Min(tl_uf, u);
        tl_vf = Min(tl_vf, v);
        br_uf = Max(br_uf, u);
        br_vf = Max(br_vf, v);

        Map_Forward((float)x, (float)(h - 1), &u, &v, oProjector, poPC);
        tl_uf = Min(tl_uf, u);
        tl_vf = Min(tl_vf, v);
        br_uf = Max(br_uf, u);
        br_vf = Max(br_vf, v);
        //printf("x:%f %f %f %f %f\n", (float)x, tl_uf, tl_vf, br_uf, br_vf);
    }
    for (int y = 0; y < h; y++)
    {
        Map_Forward(0, (float)y, &u, &v, oProjector, poPC);
        tl_uf = Min(tl_uf, u); tl_vf = Min(tl_vf, v);
        br_uf = Max(br_uf, u); br_vf = Max(br_vf, v);

        Map_Forward((float)(w - 1), (float)y, &u, &v, oProjector, poPC);
        tl_uf = Min(tl_uf, u); tl_vf = Min(tl_vf, v);
        br_uf = Max(br_uf, u); br_vf = Max(br_vf, v);
    }
    //bSave_PLY<float>("c:\\tmp\\1.ply", oPC);

    Dest_tl[0] = (int)tl_uf;
    Dest_tl[1] = (int)tl_vf;
    Dest_br[0] = (int)br_uf;
    Dest_br[1] = (int)br_vf;
    iCount++;
    return;
}
template<typename _T>void Detect_Result_Roi(int w, int h, int Dest_tl[2], int Dest_br[2],
    Sphere_Projector<_T> oProjector, Point_Cloud<float>* poPC = NULL)
{//tl: Top left     br: Bottom Right
    static int iCount = 0;

    Detect_Result_Roi_By_Border(w, h, Dest_tl, Dest_br, oProjector, poPC);
    float tl_uf = (float)Dest_tl[0];
    float tl_vf = (float)Dest_tl[1];
    float br_uf = (float)Dest_br[0];
    float br_vf = (float)Dest_br[1];

    float x = (float)oProjector.rinv[1];
    float y = (float)oProjector.rinv[4];
    float z = (float)oProjector.rinv[7];
    if (y > 0.f)
    {
        float x_ = (float)((oProjector.K[0] * x + oProjector.K[1] * y) / z + oProjector.K[2]);
        float y_ = (float)(oProjector.K[4] * y / z + oProjector.K[5]);

        if (x_ > 0.f && x_ < w && y_ > 0.f && y_ < h)
        {
            tl_uf = Min(tl_uf, 0.f); tl_vf = std::min(tl_vf, float(PI * oProjector.scale));
            br_uf = Max(br_uf, 0.f); br_vf = std::max(br_vf, float(PI * oProjector.scale));
        }
    }

    x = (float)oProjector.rinv[1];
    y = -(float)oProjector.rinv[4];
    z = (float)oProjector.rinv[7];

    if (y > 0.f)
    {
        float x_ = (float)((oProjector.K[0] * x + oProjector.K[1] * y) / z + oProjector.K[2]);
        float y_ = (float)(oProjector.K[4] * y / z + oProjector.K[5]);
        if (x_ > 0.f && x_ < w && y_ > 0.f && y_ < h)
        {
            tl_uf = Min(tl_uf, 0.f); tl_vf = Min(tl_vf, 0);
            br_uf = Max(br_uf, 0.f); br_vf = Max(br_vf, 0);
        }
    }
    Dest_tl[0] = (int)tl_uf;
    Dest_tl[1] = (int)tl_vf;
    Dest_br[0] = (int)br_uf;
    Dest_br[1] = (int)br_vf;
    iCount++;
    return;
}

template<typename _T>void Set_Camera_Params(Sphere_Projector<_T>* poProjector, _T _K[3 * 3], _T _R[3 * 3], _T t[3])
{
    static int iCount = 0;
    memcpy(poProjector->K, _K, 3 * 3 * sizeof(_T));
    //Disp(_K, 3, 3, "K");
    //Disp(_R, 3, 3, "R");*/
     //Disp(t, 1, 3, "t");

     //对R做一个转置
    Matrix_Transpose(_R, 3, 3, poProjector->rinv);
    //R_Kinv = R * K.inv();
    _T K_inv[3 * 3];
    int iResult = 0;
    Get_Inv_Matrix_Row_Op(_K, K_inv, 3, &iResult);
    /*if (iCount == 4)
        Disp(poProjector->rinv, 3, 3, "rinv");*/

    Matrix_Multiply_3x3(_R, K_inv, poProjector->r_kinv);
    /*if (iCount == 4)
        Disp(_K, 3, 3, "K");*/

        //K_Rinv = K * Rinv;
    Matrix_Multiply_3x3(_K, poProjector->rinv, poProjector->k_rinv);

    if (t)
        memcpy(poProjector->t, t, 3 * sizeof(_T));
    else
        memset(poProjector->t, 0, 3 * sizeof(_T));
    iCount++;
}

template<typename _T> void Build_Map(int w, int h, _T K[], _T R[], _T fScale, 
    /*_T** ppx_Map, _T** ppy_Map,*/
    int Dest_roi[2][2], Sphere_Projector<_T>* poProjector, Point_Cloud<float>* poPC)
{//Dest_roi: [0] 为左上角， [1]为 w,h
    static int iCount = 0;

    Sphere_Projector<_T> oProjector = {};

    //每次调用前必须将scale带入
    oProjector.scale = fScale;  // 530.474915;   //529.692383;  //529.69240946121784;
    Set_Camera_Params(&oProjector, K, R, (_T*)NULL);

    int Dest_tl[2], Dest_br[2];
    Detect_Result_Roi(w, h, Dest_tl, Dest_br, oProjector, poPC);

    Dest_roi[0][0] = Dest_tl[0];
    Dest_roi[0][1] = Dest_tl[1];

    Dest_roi[1][0] = Dest_br[0] - Dest_tl[0];
    Dest_roi[1][1] = Dest_br[1] - Dest_tl[1];

    *poProjector = oProjector;
    iCount++;
    return;
}

__device__ static void Pix_Inter_2_GPU(Image::Part_1 oImage, float x, float y, unsigned char Pix[4],
    Border_Type iImage_Border_Type, Border_Type iMask_Border_Type, int iThread_ID = 0)
{//Bi_Linear + Nearest 二合一
    short x1 = floor(x),
        y1 = floor(y);
    //h_1 = oImage.m_iHeight - 1,
    //w_1 = oImage.m_iWidth - 1;
    unsigned char A, B, C, D;
    int iCur_Line_Pos, iNext_Line_Pos;
    short xl_Pos, xr_Pos;

    //先搞Image_Warp的插值
    if (iImage_Border_Type == BORDER_CONSTANT)
    {//单领出来做
        iCur_Line_Pos = iNext_Line_Pos = -1;
        short y2 = y1;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iCur_Line_Pos = y2 * oImage.m_iWidth;
        y2++;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iNext_Line_Pos = y2 * oImage.m_iWidth;

        xl_Pos = xr_Pos = -1;
        short x2 = x1;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xl_Pos = x2;
        x2++;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xr_Pos = x2;
    }
    else
    {
        iCur_Line_Pos = iGet_Border_y_GPU(y1, oImage.m_iHeight, iImage_Border_Type) * oImage.m_iWidth;
        iNext_Line_Pos = iGet_Border_y_GPU(y1 + 1, oImage.m_iHeight, iImage_Border_Type) * oImage.m_iWidth;
        xl_Pos = iGet_Border_x_GPU(x1, oImage.m_iWidth, iImage_Border_Type);
        xr_Pos = iGet_Border_x_GPU(x1 + 1, oImage.m_iWidth, iImage_Border_Type);
    }

    for (int i = 0; i < 3; i++)
    {
        if (oImage.m_pChannel[i])
        {
            if (iCur_Line_Pos >= 0)
            {
                A = xl_Pos >= 0 ? oImage.m_pChannel[i][iCur_Line_Pos + xl_Pos] : 0;
                B = xr_Pos >= 0 ? oImage.m_pChannel[i][iCur_Line_Pos + xr_Pos] : 0;
            }
            else
                A = B = 0;
            if (iNext_Line_Pos >= 0)
            {
                C = xl_Pos >= 0 ? oImage.m_pChannel[i][iNext_Line_Pos + xl_Pos] : 0;
                D = xr_Pos >= 0 ? oImage.m_pChannel[i][iNext_Line_Pos + xr_Pos] : 0;
            }
            else
                C = D = 0;

            float fValue_0, fValue_1;
            {
                float w1 = x - x1, w0 = 1.f - w1;
                fValue_0 = w0 * A + w1 * B;
                fValue_1 = w0 * C + w1 * D;
                //if (iThread_ID == 2466)
                    //printf("%d %d %d %d\n", iCur_Line_Pos, iNext_Line_Pos, xl_Pos, xr_Pos);
            }

            {
                float w3 = y - y1, w2 = 1.f - w3;
                Pix[i] = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5f);
            }
        }
    }

    //再搞Mask的插值
    if (iMask_Border_Type == BORDER_CONSTANT)
    {//单领出来做
        iCur_Line_Pos = iNext_Line_Pos = -1;
        short y2 = y1;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iCur_Line_Pos = y2 * oImage.m_iWidth;
        y2++;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iNext_Line_Pos = y2 * oImage.m_iWidth;

        xl_Pos = xr_Pos = -1;
        short x2 = x1;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xl_Pos = x2;
        x2++;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xr_Pos = x2;
    }
    else
    {
        iCur_Line_Pos = iGet_Border_y_GPU(y1, oImage.m_iHeight, iMask_Border_Type) * oImage.m_iWidth;
        iNext_Line_Pos = iGet_Border_y_GPU(y1 + 1, oImage.m_iHeight, iMask_Border_Type) * oImage.m_iWidth;
        xl_Pos = iGet_Border_x_GPU(x1, oImage.m_iWidth, iMask_Border_Type);
        xr_Pos = iGet_Border_x_GPU(x1 + 1, oImage.m_iWidth, iMask_Border_Type);
    }

    //判断取那个位置更近
    int iLine_Pos;
    short x_Pos;
    if (y - y1 <= 0.5f)
        iLine_Pos = iCur_Line_Pos;
    else
        iLine_Pos = iNext_Line_Pos;

    if (x - x1 <= 0.5f)
        x_Pos = xl_Pos;
    else
        x_Pos = xr_Pos;

    if (iLine_Pos >= 0 && x_Pos >= 0)
        Pix[3] = oImage.m_pChannel[3][iLine_Pos + x_Pos];
    else
        Pix[3] = 0;
}

template<typename _T> __global__ void _Re_Map_2_GPU(Image::Part_1 oImage, Image::Part_1 oImage_Warped,
    short roi_x, short roi_y, float fScale, K_Rinv<_T> oK_Rinv,
    Border_Type iImage_Border_Type, Border_Type iMask_Border_Type = BORDER_CONSTANT)
{
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oImage_Warped.m_iWidth * oImage_Warped.m_iHeight)
        return;
    float x, y;
    {
        //假设u,v是从Dest的x,y 恢复到虚拟的坐标上
        float u = fScale * (int)(iThread_ID % oImage_Warped.m_iWidth + roi_x),
            v = fScale * (int)(iThread_ID / oImage_Warped.m_iWidth + roi_y);

        float sinv = sinf((float)(PI - v));
        float x_ = sinv * sinf(u);
        float y_ = (float)cosf((float)(PI - v));
        float z_ = sinv * cosf(u);

        float z;
        x = (float)(oK_Rinv.M[0] * x_ + oK_Rinv.M[1] * y_ + oK_Rinv.M[2] * z_);
        y = (float)(oK_Rinv.M[3] * x_ + oK_Rinv.M[4] * y_ + oK_Rinv.M[5] * z_);
        z = (float)(oK_Rinv.M[6] * x_ + oK_Rinv.M[7] * y_ + oK_Rinv.M[8] * z_);

        if (z > 0)
            x /= z, y /= z;
        else
            x = y = -1;
    }

    unsigned char Pix[4];
    Pix_Inter_2_GPU(oImage, x, y, Pix, iImage_Border_Type, iMask_Border_Type, iThread_ID);
    oImage_Warped.m_pChannel[0][iThread_ID] = Pix[0];
    oImage_Warped.m_pChannel[1][iThread_ID] = Pix[1];
    oImage_Warped.m_pChannel[2][iThread_ID] = Pix[2];
    oImage_Warped.m_pChannel[3][iThread_ID] = Pix[3];

    /*if (iThread_ID == 100 * oImage_Warped.m_iWidth + 100)
        printf("%d\n", Pix[3]);*/

}

template<typename _T>void Re_Map_3_GPU(Image oImage, Image oImage_Warped,
    int roi_x, int roi_y, Sphere_Projector<_T>* poProjector,
    Interpolation_Flag iImage_Inter_Type, Interpolation_Flag iMask_Inter_Type,
    Border_Type iImage_Border_Type, Border_Type iMask_Border_Type)
{//Image, Mask一起搞了
    static int iCount = 0;
    int iThread_Per_Block = 512;
    dim3 oGrid;
    oGrid.x = (oImage_Warped.m_iWidth * oImage_Warped.m_iHeight + iThread_Per_Block - 1) / iThread_Per_Block;
    oGrid.y = 1;    //oDest.m_iChannel_Count;
    oGrid.z = 1;

    _Re_Map_2_GPU<_T> << <oGrid, iThread_Per_Block >> > (oImage.m_oPart_1, oImage_Warped.m_oPart_1, roi_x, roi_y,
        (float)(1.f / poProjector->scale), poProjector->oK_Rinv, iImage_Border_Type, iMask_Border_Type);

    iCount++;
    return;
}

template<typename _T>void Re_Map_2_GPU(Image oImage, Image oMask,
    Image oImage_Warped, Image oMask_Warped,
    int roi_x, int roi_y, Sphere_Projector<_T>* poProjector,
    Interpolation_Flag iImage_Inter_Type, Interpolation_Flag iMask_Inter_Type,
    Border_Type iImage_Border_Type, Border_Type iMask_Border_Type)
{//Image, Mask一起搞了
    static int iCount = 0;
    int iThread_Per_Block = 512;
    dim3 oGrid;
    oGrid.x = (oImage_Warped.m_iWidth * oImage_Warped.m_iHeight + iThread_Per_Block - 1) / iThread_Per_Block;
    oGrid.y = 1;    //oDest.m_iChannel_Count;
    oGrid.z = 1;

    //小技巧，将Mask作为Image的Alpha通道
    oImage.m_pChannel[3] = oMask.m_pChannel[0];
    oImage_Warped.m_pChannel[3] = oMask_Warped.m_pChannel[0];
    _Re_Map_2_GPU<_T> << <oGrid, iThread_Per_Block >> > (oImage.m_oPart_1, oImage_Warped.m_oPart_1, roi_x, roi_y,
        (float)(1.f / poProjector->scale), poProjector->oK_Rinv, iImage_Border_Type, iMask_Border_Type);
    //oImage.m_pChannel[3] = oImage_Warped.m_pChannel[3] = NULL;
    //Disp_Cuda_Error();
    /*Disp_Part_GPU(oImage_Warped.m_pChannel[3], oImage_Warped.m_iWidth, 100, 100, 2, 2);
    Disp_Part_GPU(oMask_Warped.m_pChannel[0], oMask_Warped.m_iWidth, 100, 100, 2, 2);*/

    /*bSave_Image_GPU("c:\\tmp\\1.bmp", oImage_Warped);
    bSave_Image_GPU("c:\\tmp\\2.bmp", oMask_Warped);*/
    iCount++;
    return;
}

template<typename _T>void Warp_2(Image oImage, Image oMask, _T K[3 * 3], _T R[3 * 3], _T fScale,
    Image* poImage_Warped, Image* poMask_Warped,Interpolation_Flag iImage_Inter_Type, Interpolation_Flag iMask_Inter_Type,
    int Corner[2][2], Border_Type iImage_Border_Type, Border_Type iMask_Border_Type, Point_Cloud<float>* poPC)
{//将图与Mask的投影合二为一，实际上，这种话深度优化已经牺牲了程序的灵活性，
    //属于没有任何重用价值的死代码，唯独快！
    static int iCount = 0;
    //_T* puxmap, * puymap;    //暂时未知
    int Dest_roi[2][2];   //[0][0-1]: x,y [1][0-1]: w,h
    Sphere_Projector<_T> oProjector = {};
    Build_Map<_T>(oImage.m_iWidth, oImage.m_iHeight, K, R, fScale, /*&puxmap, &puymap,*/ Dest_roi, &oProjector, poPC);
    Init_Image_GPU(poImage_Warped, Dest_roi[1][0] + 1, Dest_roi[1][1] + 1, Image::IMAGE_TYPE_BMP, oImage.m_iBit_Count);
    Init_Image_GPU(poMask_Warped, Dest_roi[1][0] + 1, Dest_roi[1][1] + 1, Image::IMAGE_TYPE_BMP, 8);
    Re_Map_2_GPU<_T>(oImage, oMask, *poImage_Warped, *poMask_Warped,
        Dest_roi[0][0], Dest_roi[0][1], &oProjector,
        iImage_Inter_Type, iMask_Inter_Type, iImage_Border_Type, iMask_Border_Type);

    poImage_Warped->m_pChannel[3] = poMask_Warped->m_pChannel[0];
    if (Corner)
    {
        Corner[0][0] = Dest_roi[0][0];
        Corner[0][1] = Dest_roi[0][1];
        Corner[1][0] = Dest_roi[0][0] + Dest_roi[1][0];
        Corner[1][1] = Dest_roi[0][1] + Dest_roi[1][1];
    }
    iCount++;
    return;
}

__device__ void Pix_Bi_Linear_GPU(Image::Part_1 oImage, float x, float y, unsigned char Pix[3], Border_Type iBorder_Type, int iThread_ID = 0)
{
    short x1 = floor(x),
        y1 = floor(y);
    //h_1 = oImage.m_iHeight - 1,
    //w_1 = oImage.m_iWidth - 1;
    unsigned char A, B, C, D;
    int iCur_Line_Pos, iNext_Line_Pos;
    short xl_Pos, xr_Pos;
    //int iSize = oImage.m_iWidth * oImage.m_iHeight;

    if (iBorder_Type == BORDER_CONSTANT)
    {//单领出来做
        iCur_Line_Pos = iNext_Line_Pos = -1;
        short y2 = y1;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iCur_Line_Pos = y2 * oImage.m_iWidth;
        y2++;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iNext_Line_Pos = y2 * oImage.m_iWidth;

        xl_Pos = xr_Pos = -1;
        short x2 = x1;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xl_Pos = x2;
        x2++;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xr_Pos = x2;
    }
    else
    {
        iCur_Line_Pos = iGet_Border_y_GPU(y1, oImage.m_iHeight, iBorder_Type) * oImage.m_iWidth;
        iNext_Line_Pos = iGet_Border_y_GPU(y1 + 1, oImage.m_iHeight, iBorder_Type) * oImage.m_iWidth;
        xl_Pos = iGet_Border_x_GPU(x1, oImage.m_iWidth, iBorder_Type);
        xr_Pos = iGet_Border_x_GPU(x1 + 1, oImage.m_iWidth, iBorder_Type);
    }

    for (int i = 0; i < 3; i++)
    {
        if (oImage.m_pChannel[i])
        {
            if (iCur_Line_Pos >= 0)
            {
                A = xl_Pos >= 0 ? oImage.m_pChannel[i][iCur_Line_Pos + xl_Pos] : 0;
                B = xr_Pos >= 0 ? oImage.m_pChannel[i][iCur_Line_Pos + xr_Pos] : 0;
            }
            else
                A = B = 0;
            if (iNext_Line_Pos >= 0)
            {
                C = xl_Pos >= 0 ? oImage.m_pChannel[i][iNext_Line_Pos + xl_Pos] : 0;
                D = xr_Pos >= 0 ? oImage.m_pChannel[i][iNext_Line_Pos + xr_Pos] : 0;
            }
            else
                C = D = 0;

            float fValue_0, fValue_1;
            {
                float w1 = x - x1, w0 = 1.f - w1;
                fValue_0 = w0 * A + w1 * B;
                fValue_1 = w0 * C + w1 * D;
                //if (iThread_ID == 2466)
                    //printf("%d %d %d %d\n", iCur_Line_Pos, iNext_Line_Pos, xl_Pos, xr_Pos);
            }

            {
                float w3 = y - y1, w2 = 1.f - w3;
                Pix[i] = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5f);
            }
        }
    }
    return;
}

template<typename _T> __global__ void _Re_Map_Bi_Binear_GPU(Image::Part_1 oSource, Image::Part_1 oDest, short roi_x,
    short roi_y, float fScale, K_Rinv<_T> oK_Rinv, Border_Type iBorder_Type = BORDER_REFLECT)
{
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oDest.m_iWidth * oDest.m_iHeight)
        return;
    float x, y;
    {
        //假设u,v是从Dest的x,y 恢复到虚拟的坐标上
        float u = fScale * (int)(iThread_ID % oDest.m_iWidth + roi_x),
            v = fScale * (int)(iThread_ID / oDest.m_iWidth + roi_y);

        float sinv = sinf((float)(PI - v));
        float x_ = sinv * sinf(u);
        float y_ = (float)cosf((float)(PI - v));
        float z_ = sinv * cosf(u);

        float z;
        x = (float)(oK_Rinv.M[0] * x_ + oK_Rinv.M[1] * y_ + oK_Rinv.M[2] * z_);
        y = (float)(oK_Rinv.M[3] * x_ + oK_Rinv.M[4] * y_ + oK_Rinv.M[5] * z_);
        z = (float)(oK_Rinv.M[6] * x_ + oK_Rinv.M[7] * y_ + oK_Rinv.M[8] * z_);

        if (z > 0)
            x /= z, y /= z;
        else
            x = y = -1;
    }

    unsigned char Pix[3];
    Pix_Bi_Linear_GPU(oSource, x, y, Pix, iBorder_Type, iThread_ID);
    oDest.m_pChannel[0][iThread_ID] = Pix[0];
    oDest.m_pChannel[1][iThread_ID] = Pix[1];
    oDest.m_pChannel[2][iThread_ID] = Pix[2];

    //oDest.m_pChannel[blockIdx.y][iThread_ID] = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5);
    //if (iThread_ID == 6 * oDest.m_iWidth + 312)
    //{
    //    //printf("x:%f y:%f\n", x, y);
    //    //printf("%d %d %d %d\n", Pix[0],Pix[1],Pix[2],iThread_ID);
    //    //printf("%f %f %d\n", fValue_0, fValue_1, oDest.m_pChannel[blockIdx.y][iThread_ID]);
    //}

    return;
}

__device__ void Pix_Nearest_GPU(Image::Part_1 oImage, float x, float y, unsigned char* Pix, Border_Type iBorder_Type, int iThread_ID = 0)
{
    short x1 = floor(x),
        y1 = floor(y);

    int iCur_Line_Pos, iNext_Line_Pos;
    short xl_Pos, xr_Pos;

    if (iBorder_Type == BORDER_CONSTANT)
    {//单领出来做
        iCur_Line_Pos = iNext_Line_Pos = -1;
        short y2 = y1;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iCur_Line_Pos = y2 * oImage.m_iWidth;
        y2++;
        if (y2 >= 0 && y2 < oImage.m_iHeight)
            iNext_Line_Pos = y2 * oImage.m_iWidth;

        xl_Pos = xr_Pos = -1;
        short x2 = x1;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xl_Pos = x2;
        x2++;
        if (x2 >= 0 && x2 < oImage.m_iWidth)
            xr_Pos = x2;
    }
    else
    {
        iCur_Line_Pos = iGet_Border_y_GPU(y1, oImage.m_iHeight, iBorder_Type) * oImage.m_iWidth;
        iNext_Line_Pos = iGet_Border_y_GPU(y1 + 1, oImage.m_iHeight, iBorder_Type) * oImage.m_iWidth;
        xl_Pos = iGet_Border_x_GPU(x1, oImage.m_iWidth, iBorder_Type);
        xr_Pos = iGet_Border_x_GPU(x1 + 1, oImage.m_iWidth, iBorder_Type);
    }

    //判断取那个位置更近
    int iLine_Pos;
    short x_Pos;
    if (y - y1 <= 0.5f)
        iLine_Pos = iCur_Line_Pos;
    else
        iLine_Pos = iNext_Line_Pos;

    /*if (iThread_ID == 35300)
    {
        printf("%d %d %f\n", iLine_Pos, x_Pos, y - y1);
    }*/

    if (x - x1 <= 0.5f)
        x_Pos = xl_Pos;
    else
        x_Pos = xr_Pos;

    if (iLine_Pos >= 0 && x_Pos >= 0)
    {
        int iPos = iLine_Pos + x_Pos;
        for (int i = 0; i < 3; i++)
            if (oImage.m_pChannel[i])
                Pix[i] = oImage.m_pChannel[i][iPos];
    }
    else
        Pix[0] = Pix[1] = Pix[2] = 0;

}
template<typename _T> __global__ void _Re_Map_Nearest_GPU(Image::Part_1 oSource, Image::Part_1 oDest, short roi_x,
    short roi_y, float fScale, K_Rinv<_T> oK_Rinv, Border_Type iBorder_Type = BORDER_REFLECT)
{
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oDest.m_iWidth * oDest.m_iHeight)
        return;
    float x, y;
    {
        //假设u,v是从Dest的x,y 恢复到虚拟的坐标上
        float u = fScale * (int)(iThread_ID % oDest.m_iWidth + roi_x),
            v = fScale * (int)(iThread_ID / oDest.m_iWidth + roi_y);

        float sinv = sinf((float)(PI - v));
        float x_ = sinv * sinf(u);
        float y_ = (float)cosf((float)(PI - v));
        float z_ = sinv * cosf(u);

        float z;
        x = (float)(oK_Rinv.M[0] * x_ + oK_Rinv.M[1] * y_ + oK_Rinv.M[2] * z_);
        y = (float)(oK_Rinv.M[3] * x_ + oK_Rinv.M[4] * y_ + oK_Rinv.M[5] * z_);
        z = (float)(oK_Rinv.M[6] * x_ + oK_Rinv.M[7] * y_ + oK_Rinv.M[8] * z_);

        if (z > 0)
            x /= z, y /= z;
        else
            x = y = -1;
    }
    unsigned char Pix[3];
    Pix_Nearest_GPU(oSource, x, y, Pix, iBorder_Type, iThread_ID);
    for (int i = 0; i < 3; i++)
        if (oDest.m_pChannel[i])
            oDest.m_pChannel[i][iThread_ID] = Pix[i];

    /*if (iThread_ID == 98 * oDest.m_iWidth + 118)
    {
        printf("%f %f thread:%d\n", x, y,iThread_ID);
    }*/

    return;
}

template<typename _T>void Re_Map_GPU(Image oSource, Image oDest, int roi_x, int roi_y,
    Sphere_Projector<_T>* poProjector, Interpolation_Flag iInter_Type, Border_Type iBorder_Type)
{//打算连Map_Backward一起做
    static int iCount = 0;
    int iThread_Per_Block = 512;
    dim3 oGrid;
    oGrid.x = (oDest.m_iWidth * oDest.m_iHeight + iThread_Per_Block - 1) / iThread_Per_Block;
    oGrid.y = 1;    //oDest.m_iChannel_Count;
    oGrid.z = 1;
    //unsigned long long tStart = iGet_Tick_Count();

    //for(int i=0;i<10000;i++)
    if (iInter_Type == INTER_LINEAR)
        _Re_Map_Bi_Binear_GPU<_T> << <oGrid, iThread_Per_Block >> > (oSource.m_oPart_1, oDest.m_oPart_1, roi_x, roi_y,
            (float)(1.f / poProjector->scale), poProjector->oK_Rinv, iBorder_Type);
    else
        _Re_Map_Nearest_GPU<_T> << <oGrid, iThread_Per_Block >> > (oSource.m_oPart_1, oDest.m_oPart_1, roi_x, roi_y,
            (float)(1.f / poProjector->scale), poProjector->oK_Rinv, iBorder_Type);

    ////if (iCount == 1)
    //{
    //    Disp_Cuda_Error();
    //    //Disp_Part_GPU(oSource.m_pChannel[2], oSource.m_iWidth, 121, 99, 2, 2);
    //    //bSave_Image_GPU("c:\\tmp\\3.bmp", oDest);
    //    Compare_Image("c:\\tmp\\2.bmp", "c:\\tmp\\3.bmp",0);
    //}

    iCount++;
    return;
}

//template<typename _T>void Warp(Image oImage, _T K[3 * 3], _T R[3 * 3], _T fScale, Image* poImages_Warped,
//    Interpolation_Flag iInter_Type, int Corner[2][2], Point_Cloud<float>* poPC = NULL, Border_Type iBorder_Type = BORDER_REFLECT)
//{//投影函数，将图按照原来的K,T重新还原到归一化平面，再投影到另一个更大些的平面上
//    //具体的数学推导还需慢慢做一遍
//    static int iCount = 0;
//    _T* puxmap, * puymap;    //暂时未知
//    int Dest_roi[2][2];   //[0][0-1]: x,y [1][0-1]: w,h
//    Sphere_Projector<_T> oProjector = {};
//    //unsigned long long tStart;
//
//    //可能是利用K，将原图oImage长款，虚拟打到一个空间上
//    //Build_Map<_T>(oImage.m_iWidth, oImage.m_iHeight, K, R, fScale, &puxmap, &puymap, Dest_roi, &oProjector, poPC);
//    //Init_Image_GPU(poImages_Warped, Dest_roi[1][0] + 1, Dest_roi[1][1] + 1, Image::IMAGE_TYPE_BMP, oImage.m_iBit_Count);
//    //Re_Map_GPU(oImage, *poImages_Warped, Dest_roi[0][0], Dest_roi[0][1], &oProjector, iInter_Type, iBorder_Type);
//
//    iCount++;
//    return;
//}

template<typename _T>void Resize_Seam_Image(Stitch<_T>* poStitch, Light_Ptr oPtr)
{
    Stitch<_T> oSt = *poStitch;
    //bLoad_Image_GPU("c:\\tmp\\1.bmp", &oSt.m_pImage_Source[0]);

    //BORDER_REFLECT
    //unsigned long long tStart = iGet_Tick_Count();

    /*Bi_Linear_cv_GPU(oSt.m_pImage_Source, oSt.m_pSeam_Est,oSt.m_iImage_Count,
        (float)oSt.seam_scale, (float)oSt.seam_scale, oPtr, BORDER_REFLECT);*/
    //for(int i=0;i<10000;i++)
    Bi_Linear_cv_GPU(oSt.m_pImage_Source_Header_GPU,
        oSt.m_pSeam_Est_Header_GPU, oSt.m_iImage_Count,
        oSt.m_Source_Size[0], oSt.m_Source_Size[1],
        oSt.m_Seam_Size[0], oSt.m_Seam_Size[1], 3,
        (float)oSt.seam_scale, (float)oSt.seam_scale, BORDER_REFLECT);
    /*Bi_Linear_cv_GPU(oSt.m_pImage_Source,
            oSt.m_pSeam_Est, oSt.m_iImage_Count,
            oSt.m_Source_Size[0], oSt.m_Source_Size[1],
            oSt.m_Seam_Size[0], oSt.m_Seam_Size[1], 3,
            (float)oSt.seam_scale, (float)oSt.seam_scale, BORDER_REFLECT);*/
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //for(int i=0;i<4;i++)
    //    bSave_Image_GPU("c:\\tmp\\1.bmp", oSt.m_pSeam_Est[i]);
    //Compare_Image("c:\\tmp\\2.bmp", "c:\\tmp\\3.bmp");
    return;
}

__global__ void Feed_Split_Image(Image Warp[],
    Image_Size_In_Block Block_Per_Image[],      // unsigned char Block_Per_Image[][2],
    /*int bw,int bh,*/
    Image Block_Image[], unsigned char Space[],
    int Corner[][2][2], int Block_Corner[][2])
{
    int iThread_ID = threadIdx.y * blockDim.x + threadIdx.x;

    //要特别清晰哪些变量放在共享内存，哪些是私有变量，否则很容易出现内存存取错误
    __shared__ Image::Part_1 oImage;
    __shared__ Image oBlock;
    __shared__ int block_tl[2], block_br[2];
    __shared__ unsigned char Block_Size[2];
    __shared__ int iBlock_ID;
    if (iThread_ID == 0)
    {
        Image_Size_In_Block oBlock_Per_Image = Block_Per_Image[blockIdx.z];
        Block_Size[0] = oBlock_Per_Image.m_Block_Per_Image[0];
        Block_Size[1] = oBlock_Per_Image.m_Block_Per_Image[1];

        if (blockIdx.x < Block_Size[0] && blockIdx.y < Block_Size[1])
        {
            //int iBlock_ID = blockIdx.z * bw * bh + blockIdx.y * bw + blockIdx.x;
            iBlock_ID = oBlock_Per_Image.m_iBlock_Start + blockIdx.y * Block_Size[0] + blockIdx.x;
            oImage = Warp[blockIdx.z].m_oPart_1;
            unsigned char block_width = (oImage.m_iWidth + Block_Size[0] - 1) / Block_Size[0];
            unsigned char block_height = (oImage.m_iHeight + Block_Size[1] - 1) / Block_Size[1];
            block_tl[0] = blockIdx.x * block_width;
            block_tl[1] = blockIdx.y * block_height;
            block_br[0] = min(block_tl[0] + block_width, oImage.m_iWidth);
            block_br[1] = min(block_tl[1] + block_height, oImage.m_iHeight);

            Block_Corner[iBlock_ID][0] = Corner[blockIdx.z][0][0] + block_tl[0];
            Block_Corner[iBlock_ID][1] = Corner[blockIdx.z][0][1] + block_tl[1];

            oBlock.m_iWidth = block_br[0] - block_tl[0];
            oBlock.m_iHeight = block_br[1] - block_tl[1];
            oBlock.m_iBit_Count = 32;
            oBlock.m_iImage_Type = Image::IMAGE_TYPE_BMP;
            oBlock.m_iChannel_Count = oBlock.m_iBit_Count >> 3;

            short iBlock_Size = oBlock.m_iWidth * oBlock.m_iHeight;
            oBlock.m_pChannel[0] = &Space[iBlock_ID * 32 * 32 * 4];
            oBlock.m_pChannel[1] = oBlock.m_pChannel[0] + iBlock_Size;
            oBlock.m_pChannel[2] = oBlock.m_pChannel[1] + iBlock_Size;
            oBlock.m_pChannel[3] = oBlock.m_pChannel[2] + iBlock_Size;

            Block_Image[iBlock_ID] = oBlock;
        }
    }
    __syncthreads();

    if (blockIdx.x >= Block_Size[0] || blockIdx.y >= Block_Size[1] ||
        threadIdx.x >= oBlock.m_iWidth || threadIdx.y >= oBlock.m_iHeight)
        return;

    int iSource_Pos = (block_tl[1] + threadIdx.y) * oImage.m_iWidth + block_tl[0] + threadIdx.x;
    int iDest_Pos = threadIdx.y * oBlock.m_iWidth + threadIdx.x;

    oBlock.m_pChannel[0][iDest_Pos] = oImage.m_pChannel[0][iSource_Pos];
    oBlock.m_pChannel[1][iDest_Pos] = oImage.m_pChannel[1][iSource_Pos];
    oBlock.m_pChannel[2][iDest_Pos] = oImage.m_pChannel[2][iSource_Pos];
    oBlock.m_pChannel[3][iDest_Pos] = oImage.m_pChannel[3][iSource_Pos];
    //if (oImage.m_pChannel[0][iSource_Pos] == 123)
    //{
    //}
}
template<typename _T>__device__ _T _Sep_Filter_2D_GPU(_T A[], short w, short h, short x, short y, _T Kernel[], int iKernel_Size)
{
    float fTotal;
    int x1, i;  //,iPadding = iKernel_Size >> 1;
    _T* pLine = &A[y * w];

    for (fTotal = 0, i = 0, x1 = x - (iKernel_Size >> 1); i < iKernel_Size; x1++, i++)
    {
        float fValue;
        if (x1 < 0)
        {
            if (-x1 >= w)
                fValue = pLine[w - 1];
            else
                fValue = pLine[-x1];
        }
        else if (x1 >= w)
        {
            if (x - (x1 - x) < 0)
                fValue = pLine[0];
            else
                fValue = pLine[x - (x1 - x)];
        }
        else
            fValue = pLine[x1];
        /*if (blockIdx.y == 0 && y == 0 && x == 9)
            printf("x1:%d fValue: %f\n", x1, fValue);*/
        fTotal += fValue * Kernel[i];
    }
    return fTotal;  // / iKernel_Size;
}
template<typename _T>__global__ void Batch_Sep_Filter_2D(_T Gain[], Image_Size_In_Block Block_Per_Image[],
    _T Kernel[], int iKernel_Size, int iImage_Count, short* pN_2_A_Map, int iIter_Count)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image_Size_In_Block oBlock_Per_Image;
    __shared__ _T Share_Mid[32 * 32];
    __shared__ _T Share_Gain[32 * 32];

    if (iThread_ID == 0)
        oBlock_Per_Image = Block_Per_Image[blockIdx.y];
    __syncthreads();

    if (iThread_ID >= oBlock_Per_Image.m_Block_Per_Image[0] * oBlock_Per_Image.m_Block_Per_Image[1])
        return;

    unsigned short /*w = oBlock_Per_Image.m_Block_Per_Image[0],*/
        /*h = oBlock_Per_Image.m_Block_Per_Image[1],*/
        y = iThread_ID / oBlock_Per_Image.m_Block_Per_Image[0],
        x = iThread_ID % oBlock_Per_Image.m_Block_Per_Image[0];

    //抄数据到Buffer中
    Share_Gain[iThread_ID] = Gain[oBlock_Per_Image.m_iBlock_Start + iThread_ID];
    __syncthreads();

    for (short i = 0; i < iIter_Count; i++)
    {
        Share_Mid[x * oBlock_Per_Image.m_Block_Per_Image[1] + y] =
            _Sep_Filter_2D_GPU(Share_Gain, oBlock_Per_Image.m_Block_Per_Image[0],
                oBlock_Per_Image.m_Block_Per_Image[1], x, y, Kernel, iKernel_Size);
        __syncthreads();

        Share_Gain[y * oBlock_Per_Image.m_Block_Per_Image[0] + x] =
            _Sep_Filter_2D_GPU(Share_Mid, oBlock_Per_Image.m_Block_Per_Image[1],
                oBlock_Per_Image.m_Block_Per_Image[0], y, x, Kernel, iKernel_Size);
        __syncthreads();
    }
    Gain[oBlock_Per_Image.m_iBlock_Start + iThread_ID] = Share_Gain[iThread_ID];
}

__device__ static int bGet_Matrix_xy_GPU(int iPos, int n, int* px, int* py)
{//给定一个矩阵的阶n,和一个顺序位置，求其在矩阵中行列
    //最后不等为  (n + n-y-1)*y/2 <=N
    //先设有h, 第h-1行的个数为 n-(h-1)
    //总数为 [n + n-(h-1)]*h/2 = (2n+1 -h)*h/2 <=N
    // (2n+1)*h - h^2 <=2N
    //h^2 - (2n+1)h >=-2N
    //h^2 - (2n+1)h + 2N >=0

    //曲线开口向上
    int/* a = 1,*/
        b = -(1 + 2 * n), c = 2 * iPos;
    float fDelta = b * b - 4 * c;
    float sqrt_b_sqr_4ac = sqrt(fDelta);
    if (fDelta < 0)
    {
        printf("错误数据");
        return 0;
    }

    float h[2] = { (-b + sqrt_b_sqr_4ac) / 2.f,	//两个根
        (-b - sqrt_b_sqr_4ac) / 2.f };

    //验根，具体情况，y1[0]为大根，y1[1]为小根
    //开口向上
    int h1 = floor(h[1]);
    if (h1 > n || h1 < 0)
    {
        h1 = std::floor(h[1]);
        if (h1 > n || h1 < 0)
            return 0;
    }

    //y为什么等于高？因为在下一行开始
    int x, x_Start = h1, y = h1;
    int iLine_Start = (n + (n - y + 1)) * y / 2;
    if (iLine_Start == iPos)
        x = y;
    else
        x = iPos - iLine_Start + x_Start;
    *px = x, * py = y;
    return 1;
}

__device__ int bIs_Overlap_GPU(int tl1[2], int tl2[2], int w1, int h1, int w2, int h2, short roi[2][2])
{
    int x_tl = max(tl1[0], tl2[0]);
    int y_tl = max(tl1[1], tl2[1]);
    int x_br = min(tl1[0] + w1, tl2[0] + w2);
    int y_br = min(tl1[1] + h1, tl2[1] + h2);
    if (x_tl < x_br && y_tl < y_br)
    {
        roi[0][0] = x_tl, roi[0][1] = y_tl;
        roi[1][0] = x_br - x_tl, roi[1][1] = y_br - y_tl;
        return 1;
    }
    return 0;
}

__global__ void _Get_Overlap(Image Block_Image[], int Block_Corner[][2],
    int iBlock_Count, Overlap_Pair Pair[], int* piCount)
{
    __shared__ Overlap_Pair Pair_Buffer[1024];
    __shared__ int iPair_Count, iPos_1;

    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= (iBlock_Count + 1) * iBlock_Count / 2)
        return;
    if (threadIdx.x == 0)
        iPair_Count = 0;
    __syncthreads();

    int x, y;
    if (!bGet_Matrix_xy_GPU(iThread_ID, iBlock_Count, &x, &y))
    {
        printf("err");
        return;
    }
    /*if(x>=iBlock_Count || y>=iBlock_Count || x<y)
        printf("%d %d\n",x,y);*/
    Overlap_Pair oPair;
    if (bIs_Overlap_GPU(Block_Corner[x], Block_Corner[y],
        Block_Image[x].m_iWidth, Block_Image[x].m_iHeight,
        Block_Image[y].m_iWidth, Block_Image[y].m_iHeight,
        oPair.roi))
    {
        oPair.m_Pair[0] = x, oPair.m_Pair[1] = y;
        int iPos = atomicAdd(&iPair_Count, 1);
        Pair_Buffer[iPos] = oPair;
    }
    __syncthreads();

    //加完自己以后，加到内存中

    if (threadIdx.x == 0)
        iPos_1 = atomicAdd(piCount, iPair_Count);
    __syncthreads();
    for (int i = threadIdx.x; i < iPair_Count; i++)
        Pair[iPos_1 + i] = Pair_Buffer[i];

    return;
}

void Single_Feed_Get_Overlap(Image Block_Image[], int Block_Corner[][2], int iBlock_Count,
    Overlap_Pair** ppPair, int* piCount)
{//给定一些Block和它的位置，求其是否有交
    int iPair_Count, iSize = iBlock_Count * 10;
    int* piCount_GPU = (int*)pMalloc_GPU(iBlock_Count * 10 * sizeof(Overlap_Pair));
    cudaMemset(piCount_GPU, 0, sizeof(int));
    Overlap_Pair* pPair = (Overlap_Pair*)(piCount_GPU + 1);

    //算上三角的大小
    iSize = (iBlock_Count + 1) * iBlock_Count / 2;
    int iThread_Per_Block = 1024;
    int iBlock_Count_1 = (iSize + iThread_Per_Block - 1) / iThread_Per_Block;

    unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    {
        _Get_Overlap << <iBlock_Count_1, iThread_Per_Block >> > (Block_Image,
            Block_Corner, iBlock_Count, pPair, piCount_GPU);
    }

    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    cudaMemcpy(&iPair_Count, piCount_GPU, 4, cudaMemcpyDeviceToHost);
    if (iPair_Count > iBlock_Count * 10)
    {
        printf("iPair_Count exceeds max count\n");
        exit(0);
    }

    iSize = sizeof(int) + iPair_Count * sizeof(Overlap_Pair);
    Shrink_GPU(piCount_GPU, iSize);
    *piCount = iPair_Count;
    *ppPair = pPair;
    return;
}

__device__ static void Get_Overlap_GPU(short roi[2][2], int Corner[2], short Overlap[2][2])
{
    Overlap[0][0] = roi[0][0] - Corner[0];
    Overlap[0][1] = roi[0][1] - Corner[1];
    //roi.br() - corners[i]
    Overlap[1][0] = roi[0][0] + roi[1][0] - Corner[0] - Overlap[0][0];
    Overlap[1][1] = roi[0][1] + roi[1][1] - Corner[1] - Overlap[0][1];
}

template<typename _T> __device__ float fGet_Mod_3(_T x, _T y, _T z)
{//就是求个模
    return (_T)sqrt(x * x + y * y + z * z);
}

__global__ void Get_Intersect_1(Image Block_Image[], int iBlock_Count,
    Overlap_Pair Pair[], int iPair_Count, int Block_Corner[][2],
    int N[], double I[], unsigned char Skip[])
{//该函数内存排列有问题，即将废弃
    __shared__ unsigned int iShare_Intersect_Count;     //纯用这个 1200ms
    __shared__ unsigned short iSub_Image_Width, iSub_Image_Height,            //Crop出来的块都是一样大小
        iSub_Image_Size,                                //子图的大小
        iChannel_Count;                                 //颜色通道数，不算alpha

    __shared__ short Overlap_1[2][2], Overlap_2[2][2];
    //重点是在次换个排列，不用plannar，用rgba    
    __shared__ unsigned char Sub_Image_1[32 * 32][4], Sub_Image_2[32 * 32][4], Intersect[32 * 32];
    __shared__  Overlap_Pair oPair;
    __shared__ Image::Part_1  oBlock_Image_1, oBlock_Image_2;

    //__shared__ int Share_Intersect_Count[32];
    __shared__ float Share_Sum[2][32];  //放Isum1, Isum2

    if (threadIdx.x == 0)
    {//显然要比原来的清爽
        iShare_Intersect_Count = 0;     //交集数量设初值0
        oPair = Pair[blockIdx.x];       //注意，对于一对Block_Image, 假象[1]为y向, [0]为x向

        //根据roi与Block_Corner计算重叠区域，形成Sub_Image
        //此处需要整理数学原理
        Get_Overlap_GPU(oPair.roi, Block_Corner[oPair.m_Pair[1]], Overlap_1);
        iSub_Image_Width = Overlap_1[1][0];
        iSub_Image_Height = Overlap_1[1][1];
        iSub_Image_Size = iSub_Image_Width * iSub_Image_Height;

        Get_Overlap_GPU(oPair.roi, Block_Corner[oPair.m_Pair[0]], Overlap_2);
        iChannel_Count = Block_Image[oPair.m_Pair[0]].m_iChannel_Count;
        oBlock_Image_1 = Block_Image[oPair.m_Pair[1]].m_oPart_1;
        oBlock_Image_2 = Block_Image[oPair.m_Pair[0]].m_oPart_1;
    }
    for (int i = threadIdx.x; i < 32; i += blockDim.x)
        Share_Sum[0][i] = Share_Sum[1][i] = 0;
    __syncthreads();
    if (threadIdx.x >= iSub_Image_Size) //凡是超出Sub_Image大小的线程都可以退出了
        return;
    //Share_Intersect_Count[i] = 0;

    unsigned short Intersect_Count = 0;
    for (unsigned int iPos_d = threadIdx.x; iPos_d < iSub_Image_Size; iPos_d += blockDim.x)
    {
        //经典计算方法，给定线程求(x,y)
        unsigned short x_d = iPos_d % iSub_Image_Width;
        unsigned short y_d = iPos_d / iSub_Image_Width;
        //求原图Block_Image的位置
        unsigned short iPos_s = (Overlap_1[0][1] + y_d) * oBlock_Image_1.m_iWidth +
            Overlap_1[0][0] + x_d;
        //像素抄到Sub_Image_1中，快不快看这了
        /*Pixel_4 oPixel = { oBlock_Image_1.m_pChannel[0][iPos_s],
         oBlock_Image_1.m_pChannel[0][iPos_s],
         oBlock_Image_1.m_pChannel[0][iPos_s] ,
            oBlock_Image_1.m_pChannel[0][iPos_s] };*/

        *(Pixel_4*)Sub_Image_1[iPos_d] =
        { oBlock_Image_1.m_pChannel[0][iPos_s],
         oBlock_Image_1.m_pChannel[1][iPos_s],
         oBlock_Image_1.m_pChannel[2][iPos_s] ,
            oBlock_Image_1.m_pChannel[3][iPos_s] };

        /*Sub_Image_1[iPos_d][0] = oBlock_Image_1.m_pChannel[0][iPos_s];
        Sub_Image_1[iPos_d][1] = oBlock_Image_1.m_pChannel[1][iPos_s];
        Sub_Image_1[iPos_d][2] = oBlock_Image_1.m_pChannel[2][iPos_s];
        Sub_Image_1[iPos_d][3] = oBlock_Image_1.m_pChannel[3][iPos_s];*/

        iPos_s = (Overlap_2[0][1] + y_d) * oBlock_Image_2.m_iWidth +
            Overlap_2[0][0] + x_d;
        *(Pixel_4*)Sub_Image_2[iPos_d] = { oBlock_Image_2.m_pChannel[0][iPos_s],
            oBlock_Image_2.m_pChannel[1][iPos_s],
            oBlock_Image_2.m_pChannel[2][iPos_s],
            oBlock_Image_2.m_pChannel[3][iPos_s] };
        /*Sub_Image_2[iPos_d][0] = oBlock_Image_2.m_pChannel[0][iPos_s];
        Sub_Image_2[iPos_d][1] = oBlock_Image_2.m_pChannel[1][iPos_s];
        Sub_Image_2[iPos_d][2] = oBlock_Image_2.m_pChannel[2][iPos_s];
        Sub_Image_2[iPos_d][3] = oBlock_Image_2.m_pChannel[3][iPos_s];*/

        if (Sub_Image_1[iPos_d][3] == 255 && Sub_Image_2[iPos_d][3] == 255)
        {//此处就是交集的判断依据，两个Mask同为255则为交集
            Intersect[iPos_d] = 255;
            Intersect_Count++;
        }
        else
            Intersect[iPos_d] = 0;
    }
    atomicAdd(&iShare_Intersect_Count, Intersect_Count);
    __syncthreads();
    /*if (oPair.m_Pair[0] == 0 && oPair.m_Pair[1] == 0)
        printf("%d\n", iShare_Intersect_Count);*/

        //设置N矩阵的值，设置Skip对应的值
    if (threadIdx.x == 0)
    {
        //可见，N矩阵就是装着自己与其他图片交集的大小
        N[oPair.m_Pair[0] * iBlock_Count + oPair.m_Pair[1]] =
            N[oPair.m_Pair[1] * iBlock_Count + oPair.m_Pair[0]] =
            Max(1, iShare_Intersect_Count);

        //如果自己与其他块有交集，不许跳过
        if (oPair.m_Pair[0] != oPair.m_Pair[1])
            Skip[oPair.m_Pair[0]] = Skip[oPair.m_Pair[1]] = 0;    //表示i,j 元素有东西
        //Isum1 = Isum2 = 0;
    }
    __syncthreads();

    //再算I矩阵
    for (unsigned short iPos_d = threadIdx.x; iPos_d < iSub_Image_Size; iPos_d += blockDim.x)
    {
        //unsigned short x_d = iPos_d % iSub_Image_Width;
        //unsigned short y_d = iPos_d / iSub_Image_Width;
        if (iChannel_Count >= 3)
        {
            if (Intersect[iPos_d])
            {
                atomicAdd(&Share_Sum[0][iPos_d % 32], fGet_Mod_3<float>(Sub_Image_1[iPos_d][0],
                    Sub_Image_1[iPos_d][1],
                    Sub_Image_1[iPos_d][2]));

                atomicAdd(&Share_Sum[1][iPos_d % 32], fGet_Mod_3<float>(Sub_Image_2[iPos_d][0],
                    Sub_Image_2[iPos_d][1],
                    Sub_Image_2[iPos_d][2]));
            }
        }
        else if (iChannel_Count == 1)
        {
            printf("Not implemented yet\n");
            /*atomicAdd(&Isum1, oSub_Image_1.m_pChannel[3][iPos_d]);
            atomicAdd(&Isum2, oSub_Image_2.m_pChannel[3][iPos_d]);*/
        }
        else
        {
            printf("Not implemented yet\n");
        }
    }
    __syncthreads();

    if (threadIdx.x == 0)
    {
        for (int i = 1; i < 32; i++)
        {//直接用这个来存Isum1, Isum2
            Share_Sum[0][0] += Share_Sum[0][i];
            Share_Sum[1][0] += Share_Sum[1][i];
            /*if (oPair.m_Pair[1] == 2 && oPair.m_Pair[0] == 91)
                printf("%f %f\n", Share_Sum[0][i], Share_Sum[0][0]);*/
        }

        int iPos_1 = oPair.m_Pair[1] * iBlock_Count + oPair.m_Pair[0];
        int iPos_2 = oPair.m_Pair[0] * iBlock_Count + oPair.m_Pair[1];
        I[iPos_1] = Share_Sum[0][0] / N[iPos_1];
        I[iPos_2] = Share_Sum[1][0] / N[iPos_2];
        /*if (oPair.m_Pair[1] == 2 && oPair.m_Pair[0] == 91)
            printf("I[2,9]:%f N:%d Isam1:%f\n",
                I[iPos_1],N[iPos_1], Share_Sum[0][0]);*/
    }
    return;
}

__global__ void Get_num_eq(unsigned char Skip[], int iBlock_Count, int* piNum_eq,
    short A_2_N_Map[], short N_2_A_Map[])
{//一块搞定   
    int iNum_eq = 0;
    for (int i = 0; i < iBlock_Count; i++)
    {
        if (Skip[i] == 0)
        {
            A_2_N_Map[iNum_eq] = i;
            N_2_A_Map[i] = iNum_eq;
            iNum_eq++;
        }
    }
    *piNum_eq = iNum_eq;
    /*for(int i=0;i<iNum_eq;i++)
        printf("%d\n", A_2_N_Map[i]);*/
}

__global__ void Gen_Eq(double A[], int iOrder, /*double x[], */double b[],
    short pA_2_N_Map[], int N[], double I[], int iBlock_Count,
    double alpha, double beta)
{
    int iThread_ID = GET_THREAD_ID();
    int x1, y1;

    if (iThread_ID >= iOrder * iOrder)
        return;
    y1 = iThread_ID / iOrder;
    x1 = iThread_ID % iOrder;

    /*if (iThread_ID >= (iOrder + 1) * iOrder / 2)
        return;

    if (!bGet_Matrix_xy_GPU(iThread_ID, iOrder, &x1, &y1))
    {
        printf("error");
        return;
    }*/

    int iPos_N = pA_2_N_Map[y1] * iBlock_Count + pA_2_N_Map[x1];
    //显然这是对角线元素
    if (N[iPos_N])
    {
        atomicAdd(&b[y1], beta * N[iPos_N]);
        __threadfence();
        atomicAdd(&A[y1 * iOrder + y1], beta * N[iPos_N]);
        //if (y1 == 0 /*&& x1 == 0*/)
            //printf("%f %f\n", beta * N[iPos_N], A[y1 * iOrder + y1]);
        __threadfence();
        if (x1 != y1)
        {//非对角线元素
            atomicAdd(&A[y1 * iOrder + y1], 2 * alpha * I[iPos_N] * I[iPos_N] * N[iPos_N]);
            //if (y1 == 0)
                //printf("y_s:%d x_s:%d %f %f\n", pA_2_N_Map[y1], pA_2_N_Map[x1], I[iPos_N], A[y1 * iOrder + y1]);
            __threadfence();
            A[y1 * iOrder + x1] -= 2 * alpha * I[iPos_N] * I[pA_2_N_Map[x1] * iBlock_Count + pA_2_N_Map[y1]] * N[iPos_N];
            //if (y1 == 0 /*&& x1 == 64*/)
                //printf("%d, %d, %lf\n", y1, x1, I[pA_2_N_Map[x1] * iBlock_Count + pA_2_N_Map[y1]]);
        }
    }
    return;
}

__global__ void Set_Gain(unsigned char Skip[], double x[], int iBlock_Count,
    double Gain[], short N_2_A_Map[])
{
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= iBlock_Count)
        return;
    Gain[iThread_ID] = !Skip[iThread_ID] ? x[N_2_A_Map[iThread_ID]] : 1;
}

template<typename _T>void Single_Feed(Stitch<_T>* poStitch,
    Compensator* poComp)
{//对所有的Block求交，交集处做一个权值计算，大概这样
//内存，除了产生pGain一块外，其余全部可以释放
    Overlap_Pair* pPair;
    int iSize, iPair_Count, iBlock_Count = poStitch->m_iBlock_Count;
    Single_Feed_Get_Overlap(poStitch->m_pBlock_Image_Header_GPU,
        poStitch->m_pBlock_Corner_GPU, poStitch->m_iBlock_Count, &pPair, &iPair_Count);

    //放交集的像素个数,pN是个对称矩阵
    //int* piOverlap_Count = (int*)pMalloc_GPU(4);
    int* pN = (int*)pMalloc_GPU(iBlock_Count * iBlock_Count * sizeof(int));
    //I[i] = Isum/pN[i]  Isum1 是累加像素值之和，再除以像素个数，就是像素值平均值
    //就是这块总体是什么颜色
    double* pI = (double*)pMalloc_GPU(iBlock_Count * iBlock_Count * sizeof(double));
    cudaMemset(pN, 0, iBlock_Count * iBlock_Count * sizeof(int));
    cudaMemset(pI, 0, iBlock_Count * iBlock_Count * sizeof(double));

    iSize = iPair_Count * 2 * sizeof(Image) +
        iPair_Count * 2 * 32 * 32 * 4;
    //Image* pSub_Image_1_GPU, * pSub_Image_2_GPU;
    //unsigned char* pSpace;
    //pSub_Image_1_GPU = (Image*)pMalloc_GPU(iSize);
    //pSub_Image_2_GPU = pSub_Image_1_GPU + iPair_Count;
    //pSpace = (unsigned char*)(pSub_Image_2_GPU + iPair_Count);

    //所有的块都有可能与别人形成交集，若和别人形成不了交集
    //则跳过这个块
    unsigned char* pSkip = (unsigned char*)pMalloc_GPU(iBlock_Count);
    cudaMemset(pSkip, 1, iBlock_Count); //缺省下，所有的块都与别有交集

    dim3 oThread, oGrid;
    //oThread.x = oThread.y = 32;
    oThread.x = 512;
    oGrid.x = iPair_Count;
    Disp_Cuda_Error();
    //bSave_Image_GPU("c:\\tmp\\2.bmp", &poStitch->m_pBlock_Image_Header_GPU[2]);
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Get_Intersect_1 << <oGrid, oThread >> >
        (poStitch->m_pBlock_Image_Header_GPU, iBlock_Count, pPair, iPair_Count, poStitch->m_pBlock_Corner_GPU,
            /* pSub_Image_1_GPU, pSub_Image_2_GPU, pSpace,*/ pN, pI, pSkip);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);;

    //unsigned short* pN_2_A_Map = (unsigned short*)pMalloc_GPU(iBlock_Count);
    int* piNum_eq = (int*)pMalloc(4);  //系数矩阵A的大小
    short* pA_2_N_Map = (short*)pMalloc_GPU(iBlock_Count * 2 * sizeof(short));
    short* pN_2_A_Map = pA_2_N_Map + iBlock_Count;

    cudaMemset(pA_2_N_Map, -1, iBlock_Count * 2 * sizeof(short));
    /* tStart = iGet_Tick_Count();
     for(int i=0;i<10000;i++)*/
    Get_num_eq << <1, 1 >> > (pSkip, iBlock_Count, piNum_eq, pA_2_N_Map, pN_2_A_Map);
    Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //接着需要构造 矩阵
    int num_eq = *piNum_eq;
    iSize = num_eq * num_eq * sizeof(double) + 128 +
        num_eq * sizeof(double) + 128;
    Light_Ptr oPtr;
    Attach_Light_Ptr(oPtr, (unsigned char*)pMalloc_GPU(iSize), iSize, 0);
    double* pA_GPU, * pb_GPU, * pA, * pb;
    unsigned char* p;
    Malloc(oPtr, num_eq * num_eq * sizeof(double), p);
    pA_GPU = (double*)p;
    Malloc(oPtr, num_eq * sizeof(double), p);
    pb_GPU = (double*)p;

    Attach_Light_Ptr(oPtr, (unsigned char*)pMalloc(iSize), iSize, 0);
    Malloc(oPtr, num_eq * num_eq * sizeof(double), p);
    pA = (double*)p;
    Malloc(oPtr, num_eq * sizeof(double), p);
    pb = (double*)p;

    double* px = (double*)pMalloc(num_eq * sizeof(double));
    cudaMemset(pA_GPU, 0, num_eq * num_eq * sizeof(double));
    cudaMemset(pb_GPU, 0, num_eq * sizeof(double));
    const double beta = 100, alpha = 0.01;

    int iThread_Count = num_eq * num_eq;    //(num_eq + 1) * num_eq / 2;
    oThread.x = 1024;
    oGrid.x = (iThread_Count + oThread.x - 1) / oThread.x;
    Disp_Cuda_Error();

    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    {
        Gen_Eq << <oGrid, oThread >> > (pA_GPU, num_eq, /*px,*/ pb_GPU, pA_2_N_Map, pN, pI,
            iBlock_Count, alpha, beta);
        cudaMemcpy(pA, pA_GPU, iSize, cudaMemcpyDeviceToHost);
    }
    Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //解方程
    int iResult;


    //tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    //CPU解 1400ms, 可以忍受
    Solve_Linear_Gause_AAt(pA, num_eq, pb, px, &iResult);   //|px| = 16.933621

    //printf("%lld\n", iGet_Tick_Count() - tStart);
    printf("Sum of x:%f\n", fGet_Mod(px, num_eq));

    double* pGain = (double*)pMalloc_GPU(iBlock_Count * sizeof(double));
    oThread.x = 512;
    oGrid.x = (iBlock_Count + oThread.x - 1) / oThread.x;

    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Set_Gain << <oGrid, oThread >> > (pSkip, px, iBlock_Count, pGain, pN_2_A_Map);    //78ms, 没有优化必要
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //Disp_GPU(pGain, 1, iBlock_Count);
    poComp->m_pGain = pGain;

    //然后将有的没得全释放
    Free_GPU(((unsigned int*)pPair) - 1);   //此处有些奇怪，以后再收拾它
    Free_GPU(pN);
    Free_GPU(pI);
    Free_GPU(pSkip);
    Free(piNum_eq);
    Free_GPU(pA_2_N_Map);
    Free_GPU(pA_GPU);
    Free(pA);
    Free(px);
    return;
}

template<typename _T>void Feed(Stitch<_T>* poStitch, Compensator* poComp)
{
    short* pN_2_A_Map = NULL;
    for (int n = 0; n < poComp->nr_feed; n++)
    {
        Single_Feed(poStitch, poComp);
    }
    //poComp->m_pGain_Map = (double**)pMalloc(poStitch->m_iImage_Count * sizeof(double*));
    //_T* pKer = poStitch->m_pKer_GPU;
    int iMax_Size = 0;
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        int iSize = poStitch->m_pBlock_Per_Image[i].m_Block_Per_Image[0] * poStitch->m_pBlock_Per_Image[i].m_Block_Per_Image[1];
        if (iSize > iMax_Size)
            iMax_Size = iSize;
    }
    if (iMax_Size >= 32 * 32)
    {
        printf("Exceed max size in Feed\n");
        exit(0); //无可挽回
    }
    //此处做一个复杂的可分离卷积
    //计算Mid的大小
    //_T* pMid = (_T*)pMalloc_GPU(poStitch->m_iBlock_Count * sizeof(_T));
    dim3 oThread, oGrid;
    oThread.x = Min(iMax_Size, 512);
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = poStitch->m_iImage_Count;

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Batch_Sep_Filter_2D<_T> << <oGrid, oThread >> > (poComp->m_pGain, poStitch->m_pBlock_Per_Image_GPU,
        poStitch->m_pKer_GPU, 3, /*pMid,*/ poStitch->m_iImage_Count, pN_2_A_Map, 2);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    ////可能没什么用，只能不要
    //for (int i = 0; i < poStitch->m_iImage_Count; i++)
    //{
    //    Image_Block oBlock = poStitch->m_pBlock_Per_Image[i];
    //    poComp->m_pGain_Map[i] = &poComp->m_pGain[oBlock.m_iBlock_Start];
    //}
    return;
}

template<typename _T>__device__ static _T Pix_Inter_2_GPU(_T Image[], short w, short h, float x, float y,
    Border_Type iImage_Border_Type, int iThread_ID = 0)
{
    short x1 = floor(x),
        y1 = floor(y);
    //h_1 = oImage.m_iHeight - 1,
    //w_1 = oImage.m_iWidth - 1;
    _T A, B, C, D;
    int iCur_Line_Pos, iNext_Line_Pos;
    short xl_Pos, xr_Pos;

    //先搞Image_Warp的插值
    if (iImage_Border_Type == BORDER_CONSTANT)
    {//单领出来做
        iCur_Line_Pos = iNext_Line_Pos = -1;
        short y2 = y1;
        if (y2 >= 0 && y2 < h)
            iCur_Line_Pos = y2 * w;
        y2++;
        if (y2 >= 0 && y2 < h)
            iNext_Line_Pos = y2 * w;

        xl_Pos = xr_Pos = -1;
        short x2 = x1;
        if (x2 >= 0 && x2 < w)
            xl_Pos = x2;
        x2++;
        if (x2 >= 0 && x2 < w)
            xr_Pos = x2;
    }
    else
    {
        iCur_Line_Pos = iGet_Border_y_GPU(y1, h, iImage_Border_Type) * w;
        iNext_Line_Pos = iGet_Border_y_GPU(y1 + 1, h, iImage_Border_Type) * w;
        xl_Pos = iGet_Border_x_GPU(x1, w, iImage_Border_Type);
        xr_Pos = iGet_Border_x_GPU(x1 + 1, w, iImage_Border_Type);
    }

    if (iCur_Line_Pos >= 0)
    {
        A = xl_Pos >= 0 ? Image[iCur_Line_Pos + xl_Pos] : 0;
        B = xr_Pos >= 0 ? Image[iCur_Line_Pos + xr_Pos] : 0;
    }
    else
        A = B = 0;
    if (iNext_Line_Pos >= 0)
    {
        C = xl_Pos >= 0 ? Image[iNext_Line_Pos + xl_Pos] : 0;
        D = xr_Pos >= 0 ? Image[iNext_Line_Pos + xr_Pos] : 0;
    }
    else
        C = D = 0;


    float fValue_0, fValue_1;
    {
        float w1 = x - x1, w0 = 1.f - w1;
        fValue_0 = w0 * A + w1 * B;
        fValue_1 = w0 * C + w1 * D;
    }

    {
        float w3 = y - y1, w2 = 1.f - w3;
        return w2 * fValue_0 + w3 * fValue_1;   //由于不用舍入，故此不必加上0.5
    }
}

template<typename _T>__global__ void _Block_Compensate(Image Warp[], _T Gain[],
    Image_Size_In_Block Block_Per_Image[], _T Aux[] = NULL)
{//这个函数已经完全对齐opencv，除了个别像素因为前面的浮点运算有轻微的差别意外
    //属于前面生成增益的问题，不属于这个函数本身的问题
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image_Size_In_Block oBlock_Per_Image;
    __shared__ Image::Part_1 oImage;
    __shared__ _T Share_Gain[32 * 32];

    if (threadIdx.x == 0)
    {
        oImage = Warp[blockIdx.y].m_oPart_1;
        oBlock_Per_Image = Block_Per_Image[blockIdx.y];
    }
    __syncthreads();

    if (iThread_ID >= oImage.m_iHeight * oImage.m_iWidth)
        return;
    //printf("Thread:%d %d\n", iThread_ID,oImage.m_iWidth*oImage.m_iHeight );
    
        
    short w = oBlock_Per_Image.m_Block_Per_Image[0],
        h = oBlock_Per_Image.m_Block_Per_Image[1],
        y = iThread_ID / oImage.m_iWidth,
        x = iThread_ID % oImage.m_iWidth;

    //Gain足够小，装入共享内存了事
    int iSize = w * h;

    for (int i = threadIdx.x; i < iSize; i += blockDim.x)
        Share_Gain[i] = Gain[oBlock_Per_Image.m_iBlock_Start + i];
    __syncthreads();

    float f_y = (float)h / oImage.m_iHeight,
        f_x = (float)w / oImage.m_iWidth;

    //opencv方案
    //y_s_f = (y + 0.5f) * f_y - 0.5f;
    //x_s_f = (x + 0.5f) * f_x - 0.5f;
    _T fPix = Pix_Inter_2_GPU(Share_Gain, w, h, (x + 0.5f) * f_x - 0.5f, (y + 0.5f) * f_y - 0.5f, BORDER_REFLECT, iThread_ID);
    /*unsigned*/ short iPix;
    iPix = oImage.m_pChannel[0][iThread_ID] * fPix + 0.5f;
    oImage.m_pChannel[0][iThread_ID] = Clip(iPix);

    iPix = oImage.m_pChannel[1][iThread_ID] * fPix + 0.5f;
    oImage.m_pChannel[1][iThread_ID] = Clip(iPix);

    iPix = oImage.m_pChannel[2][iThread_ID] * fPix + 0.5f;
    oImage.m_pChannel[2][iThread_ID] = Clip(iPix);
    /*if (Aux && blockIdx.y == 0 && y == 1089 && x == 257)
    {
        printf("B:%d\n",iPix);
    }*/
    return;
}

template<typename _T>void Block_Compensate(Stitch<_T>* poStitch, 
    Image Image_Warp[],Image Image_Warp_Header_GPU[])
{//对原来的Warp 投影图进行补光？
    //此处不搞扩大为Warp Image那套，直接算
    static int iCount = 0;
    int iMax_Size = 0;
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        /*int iSize = poStitch->m_pImage_Warp[i].m_iWidth *
            poStitch->m_pImage_Warp[i].m_iHeight;*/
        int iSize = Image_Warp[i].m_iWidth * Image_Warp[i].m_iHeight;
        if (iSize > iMax_Size)
            iMax_Size = iSize;
        /*if (iCount == 1)
            printf("%d\n", iSize);*/
    }
    dim3 oThread, oGrid;

    oThread.x = 512;
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = poStitch->m_iImage_Count;
    //Disp_Part_GPU<unsigned char>(poStitch->m_pImage_Warp[0].m_pChannel[1], poStitch->m_pImage_Warp[0].m_iWidth, 129, 0, 3, 2);
    //double* pGain_Resize = (double*)pMalloc(iMax_Size * sizeof(double*));
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    //_Block_Compensate << <oGrid, oThread >> > (poStitch->m_pImage_Warp_Header_GPU,
    //    poStitch->m_oComp.m_pGain, poStitch->m_pBlock_Per_Image_GPU/*, pGain_Resize*/);

    /*_T* pAux = NULL;
    if (iCount == 1)
        pAux = (_T*)pMalloc(10);*/

    _Block_Compensate << <oGrid, oThread >> > (Image_Warp_Header_GPU,
        poStitch->m_oComp.m_pGain, poStitch->m_pBlock_Per_Image_GPU);
    //Disp_Cuda_Error();
    //if (iCount == 1)
    //{
    //    //bSave_Image_GPU("c:\\tmp\\1.bmp", &Image_Warp_Header_GPU[0]);
    //    //Compare_Image("c:\\tmp\\Warp_Comp_0.bmp", "c:\\tmp\\1.bmp",1);
    //    //printf("Here");
    //}
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //Image oImage = poStitch->m_pImage_Warp[0];
    //Disp_Part<double>(pGain_Resize, oImage.m_iWidth, 0, 0, oImage.m_iWidth, 1);
    //Disp_Part<double>(pGain_Resize, oImage.m_iWidth, 129,0, 3, 2);
    //Disp_Part_GPU<unsigned char>(oImage.m_pChannel[1], oImage.m_iWidth, 130, 0, 1, 1);

    /*bSave_Image_GPU("c:\\tmp\\2.bmp", &poStitch->m_pImage_Warp_Header_GPU[3]);
    Compare_Image("c:\\tmp\\1.bmp", "c:\\tmp\\2.bmp");*/
    iCount++;
    return;
}

template<typename _T>void Feed(Stitch<_T>* poStitch)
{//这个是否有普遍性？不好说。如果有，则做成通用接口
    //在次开辟也无妨
    int (*Corner)[2][2] = poStitch->m_pCorner;
    Image* Images_Warped = poStitch->m_pImage_Warp;
    //*Masks_Warped = poStitch->m_pMasks_Warped;
    int iBlock_Count = 0;

    //先统计数据量
    for (int iImage_Index = 0; iImage_Index < poStitch->m_iImage_Count; iImage_Index++)
    {
        int block_per_img[2] = { (Images_Warped[iImage_Index].m_iWidth + 32 - 1) / 32,
                                (Images_Warped[iImage_Index].m_iHeight + 32 - 1) / 32 };

        poStitch->m_pBlock_Per_Image[iImage_Index].m_Block_Per_Image[0] = block_per_img[0];
        poStitch->m_pBlock_Per_Image[iImage_Index].m_Block_Per_Image[1] = block_per_img[1];
        if (iImage_Index == 0)
            poStitch->m_pBlock_Per_Image[iImage_Index].m_iBlock_Start = 0;
        else
        {
            poStitch->m_pBlock_Per_Image[iImage_Index].m_iBlock_Start =
                poStitch->m_pBlock_Per_Image[iImage_Index - 1].m_iBlock_Start +
                poStitch->m_pBlock_Per_Image[iImage_Index - 1].m_Block_Per_Image[0] *
                poStitch->m_pBlock_Per_Image[iImage_Index - 1].m_Block_Per_Image[1];
        }
        iBlock_Count += block_per_img[0] * block_per_img[1];
    }
    poStitch->m_iBlock_Count = iBlock_Count;
    cudaMemcpy(poStitch->m_pBlock_Per_Image_GPU, poStitch->m_pBlock_Per_Image,
        poStitch->m_iImage_Count * sizeof(Image_Size_In_Block), cudaMemcpyHostToDevice);

    /********************大图分块*********************************/
    dim3 oBlock, oGrid;
    int i, iMax_bx = 0, iMax_by = 0;
    for (i = 0; i < poStitch->m_iImage_Count; i++)
    {
        if (poStitch->m_pBlock_Per_Image[i].m_Block_Per_Image[0] > iMax_bx)
            iMax_bx = poStitch->m_pBlock_Per_Image[i].m_Block_Per_Image[0];
        if (poStitch->m_pBlock_Per_Image[i].m_Block_Per_Image[1] > iMax_by)
            iMax_by = poStitch->m_pBlock_Per_Image[i].m_Block_Per_Image[1];
    }

    //此处浪费点内存了事
    oBlock.x = 32;
    oBlock.y = 32;

    oGrid.x = iMax_bx;
    oGrid.y = iMax_by;
    oGrid.z = 4;    // poStitch->m_iImage_Count;

    poStitch->m_pBlock_Image = (Image*)pMalloc(poStitch->m_iImage_Count * iMax_bx * iMax_by * sizeof(Image));
    poStitch->m_pBlock_Image_Header_GPU = (Image*)pMalloc_GPU(poStitch->m_iImage_Count * iMax_bx * iMax_by * sizeof(Image));
    poStitch->m_pBlock_Image_Data_GPU = (unsigned char*)pMalloc_GPU(poStitch->m_iImage_Count * iMax_bx * iMax_by * 32 * 32 * 4);
    //cudaMemset(poStitch->m_pBlock_Image_Data_GPU, 0, poStitch->m_iImage_Count * iMax_bx * iMax_by * 32 * 32 * 4);

    poStitch->m_pBlock_Corner = (int(*)[2])pMalloc(poStitch->m_iImage_Count * iMax_bx * iMax_by * 2 * sizeof(int));
    poStitch->m_pBlock_Corner_GPU = (int(*)[2])pMalloc_GPU(poStitch->m_iImage_Count * iMax_bx * iMax_by * 2 * sizeof(int));

    //此处是个祸害，应该在Warp的时候就顺所做了
    cudaMemcpy(poStitch->m_pImage_Warp_Header_GPU, poStitch->m_pImage_Warp,
        poStitch->m_iImage_Count * sizeof(Image), cudaMemcpyHostToDevice);

    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Feed_Split_Image << <oGrid, oBlock >> > (poStitch->m_pImage_Warp_Header_GPU,
        poStitch->m_pBlock_Per_Image,/*iMax_bx,iMax_by,*/
        poStitch->m_pBlock_Image_Header_GPU, poStitch->m_pBlock_Image_Data_GPU,
        Corner, poStitch->m_pBlock_Corner_GPU);

    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);
    /********************大图分块*********************************/

    //for (i = 0; i < iBlock_Count; i++)
    //{
    //    char File_2[256];
    //    sprintf(File_2, "c:\\tmp\\2\\%04d.bmp", i);
    //    bSave_Image_GPU(File_2, &poStitch->m_pBlock_Image_Header_GPU[i]);
    //    //bSave_Comp_GPU(File_2, &poStitch->m_pBlock_Image_Header_GPU[i], 3);
    //    ////printf("%s\n", File_2);
    //    //char File_1[256];
    //    //sprintf(File_1, "c:\\tmp\\1\\%04d.bmp", i);
    //    //if (!Compare_Image(File_1, File_2))
    //    //    break;
    //}

    Compensator* poComp = &poStitch->m_oComp;
    Feed<_T>(poStitch, poComp);

    Block_Compensate(poStitch,poStitch->m_pImage_Warp,poStitch->m_pImage_Warp_Header_GPU);
    return;
}

int bIs_Overlap(int tl1[2], int tl2[2], int w1, int h1, int w2, int h2, short roi[2][2])
{
    int x_tl = max(tl1[0], tl2[0]);
    int y_tl = max(tl1[1], tl2[1]);
    int x_br = min(tl1[0] + w1, tl2[0] + w2);
    int y_br = min(tl1[1] + h1, tl2[1] + h2);
    if (x_tl < x_br && y_tl < y_br)
    {
        if (roi)
        {
            roi[0][0] = x_tl, roi[0][1] = y_tl;
            roi[1][0] = x_br - x_tl, roi[1][1] = y_br - y_tl;
        }
        return 1;
    }
    return 0;
}

__global__ void Find_Pair_Copy_Image_1(Image Warp[], int Corner[][2][2],
    int iGap, Image_Pair oPair)
{//同时计算dx,dy, Pair的两图
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image::Part_1 oImage_A, oImage_B;
    __shared__ int Share_Corner[2][2][2];

    if (threadIdx.x == 0)
    {
        oImage_A = Warp[oPair.m_iImage_A].m_oPart_1;
        oImage_B = Warp[oPair.m_iImage_B].m_oPart_1;
        Share_Corner[0][0][0] = Corner[oPair.m_iImage_A][0][0];
        Share_Corner[0][0][1] = Corner[oPair.m_iImage_A][0][1];
        Share_Corner[0][1][0] = Corner[oPair.m_iImage_A][1][0];
        Share_Corner[0][1][1] = Corner[oPair.m_iImage_A][1][1];
        Share_Corner[1][0][0] = Corner[oPair.m_iImage_B][0][0];
        Share_Corner[1][0][1] = Corner[oPair.m_iImage_B][0][1];
        Share_Corner[1][1][0] = Corner[oPair.m_iImage_B][1][0];
        Share_Corner[1][1][1] = Corner[oPair.m_iImage_B][1][1];
    }
    __syncthreads();

    if (iThread_ID >= oPair.m_oImage_A.m_iWidth * oPair.m_oImage_A.m_iHeight)
        return;

    short y = iThread_ID / oPair.m_oImage_A.m_iWidth;
    short x = iThread_ID % oPair.m_oImage_A.m_iWidth;
    //if (blockIdx.y == 0 && y == 100 /*&& x == 100*/)
    //    printf("Here");

    {
        //表示的是原图的位置
        short y_A = oPair.roi[0][1] - Share_Corner[0][0][1] + y - iGap;
        short x_A = oPair.roi[0][0] - Share_Corner[0][0][0] + x - iGap;
        //int iPos = (y + iGap) * oPair.m_oImage_A.m_iWidth + x + iGap

        if (y_A >= 0 && x_A >= 0 && y_A < oImage_A.m_iHeight && x_A < oImage_A.m_iWidth)
        {
            int iPos_1 = y_A * oImage_A.m_iWidth + x_A;
            oPair.m_oImage_A.m_pChannel[0][iThread_ID] = oImage_A.m_pChannel[0][iPos_1];
            oPair.m_oImage_A.m_pChannel[1][iThread_ID] = oImage_A.m_pChannel[1][iPos_1];
            oPair.m_oImage_A.m_pChannel[2][iThread_ID] = oImage_A.m_pChannel[2][iPos_1];
            oPair.m_oImage_A.m_pChannel[3][iThread_ID] = oImage_A.m_pChannel[3][iPos_1];
        }
        else
        {
            oPair.m_oImage_A.m_pChannel[0][iThread_ID] =
                oPair.m_oImage_A.m_pChannel[1][iThread_ID] =
                oPair.m_oImage_A.m_pChannel[2][iThread_ID] =
                oPair.m_oImage_A.m_pChannel[3][iThread_ID] = 0;
        }
    }

    {
        //表示的是原图的位置
        short y_B = oPair.roi[0][1] - Share_Corner[1][0][1] + y - iGap;
        short x_B = oPair.roi[0][0] - Share_Corner[1][0][0] + x - iGap;

        if (y_B >= 0 && x_B >= 0 && y_B < oImage_B.m_iHeight && x_B < oImage_A.m_iWidth)
        {
            int iPos_1 = y_B * oImage_B.m_iWidth + x_B;
            oPair.m_oImage_B.m_pChannel[0][iThread_ID] = oImage_B.m_pChannel[0][iPos_1];
            oPair.m_oImage_B.m_pChannel[1][iThread_ID] = oImage_B.m_pChannel[1][iPos_1];
            oPair.m_oImage_B.m_pChannel[2][iThread_ID] = oImage_B.m_pChannel[2][iPos_1];
            oPair.m_oImage_B.m_pChannel[3][iThread_ID] = oImage_B.m_pChannel[3][iPos_1];
        }
        else
        {
            oPair.m_oImage_B.m_pChannel[0][iThread_ID] =
                oPair.m_oImage_B.m_pChannel[1][iThread_ID] =
                oPair.m_oImage_B.m_pChannel[2][iThread_ID] =
                oPair.m_oImage_B.m_pChannel[3][iThread_ID] = 0;
        }
    }
    return;
}

__global__ void Find_Pair_Copy_Image(Image Warp[], int Corner[][2][2],
    int iGap, int* dx[], int* dy[], Image_Pair Pair[])
{//同时计算dx,dy, Pair的两图
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image_Pair oPair;
    __shared__ Image::Part_1 oImage_A, oImage_B;
    __shared__ int Share_Corner[2][2][2];

    if (threadIdx.x == 0)
    {
        oPair = Pair[blockIdx.y];
        oImage_A = Warp[oPair.m_iImage_A].m_oPart_1;
        oImage_B = Warp[oPair.m_iImage_B].m_oPart_1;
        Share_Corner[0][0][0] = Corner[oPair.m_iImage_A][0][0];
        Share_Corner[0][0][1] = Corner[oPair.m_iImage_A][0][1];
        Share_Corner[0][1][0] = Corner[oPair.m_iImage_A][1][0];
        Share_Corner[0][1][1] = Corner[oPair.m_iImage_A][1][1];
        Share_Corner[1][0][0] = Corner[oPair.m_iImage_B][0][0];
        Share_Corner[1][0][1] = Corner[oPair.m_iImage_B][0][1];
        Share_Corner[1][1][0] = Corner[oPair.m_iImage_B][1][0];
        Share_Corner[1][1][1] = Corner[oPair.m_iImage_B][1][1];
    }
    __syncthreads();

    if (iThread_ID >= oPair.m_oImage_A_GPU.m_iWidth * oPair.m_oImage_A_GPU.m_iHeight)
        return;

    short y = iThread_ID / oPair.m_oImage_A_GPU.m_iWidth;
    short x = iThread_ID % oPair.m_oImage_A_GPU.m_iWidth;
    {
        //表示的是原图的位置
        short y_A = oPair.roi[0][1] - Share_Corner[0][0][1] + y - iGap;
        short x_A = oPair.roi[0][0] - Share_Corner[0][0][0] + x - iGap;
        //int iPos = (y + iGap) * oPair.m_oImage_A.m_iWidth + x + iGap

        if (y_A >= 0 && x_A >= 0 && y_A < oImage_A.m_iHeight && x_A < oImage_A.m_iWidth)
        {
            int iPos_1 = y_A * oImage_A.m_iWidth + x_A;
            oPair.m_oImage_A_GPU.m_pChannel[0][iThread_ID] = oImage_A.m_pChannel[0][iPos_1];
            oPair.m_oImage_A_GPU.m_pChannel[1][iThread_ID] = oImage_A.m_pChannel[1][iPos_1];
            oPair.m_oImage_A_GPU.m_pChannel[2][iThread_ID] = oImage_A.m_pChannel[2][iPos_1];
            oPair.m_oImage_A_GPU.m_pChannel[3][iThread_ID] = oImage_A.m_pChannel[3][iPos_1];
            oPair.m_pSub_dx_A[iThread_ID] = dx[oPair.m_iImage_A][iPos_1];
            oPair.m_pSub_dy_A[iThread_ID] = dy[oPair.m_iImage_A][iPos_1];
            /*if (blockIdx.y == 0 && y == 50 && x == 0)
                printf("%d Thread:%d \n", dx[oPair.m_iImage_A][iPos_1],iThread_ID);*/
        }
        else
        {
            oPair.m_oImage_A_GPU.m_pChannel[0][iThread_ID] =
                oPair.m_oImage_A_GPU.m_pChannel[1][iThread_ID] =
                oPair.m_oImage_A_GPU.m_pChannel[2][iThread_ID] =
                oPair.m_oImage_A_GPU.m_pChannel[3][iThread_ID] = 0;
            oPair.m_pSub_dx_A[iThread_ID] =
                oPair.m_pSub_dy_A[iThread_ID] = 0;
        }
    }

    {
        //表示的是原图的位置
        short y_B = oPair.roi[0][1] - Share_Corner[1][0][1] + y - iGap;
        short x_B = oPair.roi[0][0] - Share_Corner[1][0][0] + x - iGap;

        if (y_B >= 0 && x_B >= 0 && y_B < oImage_B.m_iHeight && x_B < oImage_A.m_iWidth)
        {
            int iPos_1 = y_B * oImage_B.m_iWidth + x_B;
            oPair.m_oImage_B_GPU.m_pChannel[0][iThread_ID] = oImage_B.m_pChannel[0][iPos_1];
            oPair.m_oImage_B_GPU.m_pChannel[1][iThread_ID] = oImage_B.m_pChannel[1][iPos_1];
            oPair.m_oImage_B_GPU.m_pChannel[2][iThread_ID] = oImage_B.m_pChannel[2][iPos_1];
            oPair.m_oImage_B_GPU.m_pChannel[3][iThread_ID] = oImage_B.m_pChannel[3][iPos_1];
            oPair.m_pSub_dx_B[iThread_ID] = dx[oPair.m_iImage_B][iPos_1];
            oPair.m_pSub_dy_B[iThread_ID] = dy[oPair.m_iImage_B][iPos_1];
            /*if (blockIdx.y == 0 && y == 50 && x == oPair.m_oImage_B.m_iWidth - 1)
                printf("%d %d\n", y_B, x_B);*/
        }
        else
        {
            oPair.m_oImage_B_GPU.m_pChannel[0][iThread_ID] =
                oPair.m_oImage_B_GPU.m_pChannel[1][iThread_ID] =
                oPair.m_oImage_B_GPU.m_pChannel[2][iThread_ID] =
                oPair.m_oImage_B_GPU.m_pChannel[3][iThread_ID] = 0;
            oPair.m_pSub_dx_B[iThread_ID] =
                oPair.m_pSub_dy_B[iThread_ID] = 0;
        }
    }

    return;
}
template<typename _T>float fGet_Distance_Sqr(_T V_1[], _T V_2[], int n)
{//取欧几里得距离的平方
    float  fTotal = 0;
    for (int i = 0; i < n; i++)
        fTotal += ((float)V_1[i] - V_2[i]) * (V_1[i] - V_2[i]);

    return fTotal;
}
void Add_Edge(GCGraph* poGraph, int i, int j, float w, float revw)
{//网络流加边
    if (!poGraph->m_iCur_Edge)
    {
        poGraph->m_Edge[0] = poGraph->m_Edge[1] = {};
        poGraph->m_iCur_Edge = 2;
    }

    Graph_Edge fromI, toI;
    //雕虫小技，类似建树，将点放在Parent的第一个孩子
    fromI.dst = j;  //表示该点流向哪个点？
    fromI.next = poGraph->m_Vertex[i].first;
    fromI.weight = w;
    poGraph->m_Vertex[i].first = poGraph->m_iCur_Edge;
    poGraph->m_Edge[poGraph->m_iCur_Edge++] = fromI;
    /*if (poGraph->m_iCur_Edge >= 164775)
        printf("Here");*/

    toI.dst = i;
    toI.next = poGraph->m_Vertex[j].first;
    toI.weight = revw;
    poGraph->m_Vertex[j].first = poGraph->m_iCur_Edge;
    poGraph->m_Edge[poGraph->m_iCur_Edge++] = toI;

}
void Set_Graph_Weights_Color(Image oImage_A, Image oImage_B, GCGraph* poGraph)
{//没有必要搞个Mask，实际上Alpha Channel就是Mask
    int i, v, x, y, iSize = oImage_A.m_iWidth * oImage_B.m_iHeight;// ,
        //iCur_Edge = 0;
    *poGraph = {};
    Graph_Vertex* pVertex = (Graph_Vertex*)pMalloc(iSize * sizeof(Graph_Vertex));
    poGraph->m_iMax_Vertex_Count = iSize;
    memset(pVertex, 0, iSize * sizeof(Graph_Vertex));

    //此处应该可以算出来，后面不用Shrink
    //暂时先跟opencv的算法，边数量尽可能加
    poGraph->m_iMax_Edge_Count = (iSize * 2 - oImage_A.m_iWidth - oImage_A.m_iHeight) * 2 + 2; //加2是为了一个源点一个汇点？

    Graph_Edge* pEdge = (Graph_Edge*)pMalloc(poGraph->m_iMax_Edge_Count * sizeof(Graph_Edge));
    const float terminal_cost_ = 10000;
    for (i = 0; i < iSize; i++)
    {
        Graph_Vertex* pV = &pVertex[i];
        float dw = pV->weight;
        //注意，source表示源点，sink表示汇点
        /*if (oMask_A.m_pChannel[0][i])
            printf("Here");*/
        float sourceW = oImage_A.m_pChannel[3][i] ? terminal_cost_ : 0;
        float sinkW = oImage_B.m_pChannel[3][i] ? terminal_cost_ : 0;
        if (dw > 0)
            sourceW += dw;
        else
            sinkW -= dw;

        poGraph->flow += (sourceW < sinkW) ? sourceW : sinkW;
        pV->weight = sourceW - sinkW;
        //printf("%f %f\n", poGraph->flow, pV->weight);
    }
    poGraph->m_iMax_Vertex_Count = iSize;
    poGraph->m_Edge = pEdge;
    poGraph->m_Vertex = pVertex;

    const float weight_eps = 1.f, bad_region_penalty_ = 1000;
    int w_1 = oImage_A.m_iWidth - 1,
        h_1 = oImage_A.m_iHeight - 1;
    i = 0;
    for (y = 0; y < oImage_A.m_iHeight; y++)
    {
        for (x = 0; x < oImage_A.m_iWidth; x++, i++)
        {
            //if (y == 271 && x == 136)
            /*if (y == 272 && x == 135)
                printf("Here");*/
            if (x < w_1)
            {
                float weight;
                {//取两图对应像素间的距离平方，反映对应像素间得差异程度
                    unsigned char Pix_3_A[] = { oImage_A.m_pChannel[0][i],oImage_A.m_pChannel[1][i],oImage_A.m_pChannel[2][i] },
                        Pix_3_B[] = { oImage_B.m_pChannel[0][i],oImage_B.m_pChannel[1][i],oImage_B.m_pChannel[2][i] };
                    weight = fGet_Distance_Sqr(Pix_3_A, Pix_3_B, 3);
                }

                {//取两Mssk间对应位置的差的平方
                    v = i + 1;
                    unsigned char Pix_3_A[] = { oImage_A.m_pChannel[0][v],oImage_A.m_pChannel[1][v],oImage_A.m_pChannel[2][v] },
                        Pix_3_B[] = { oImage_B.m_pChannel[0][v],oImage_B.m_pChannel[1][v],oImage_B.m_pChannel[2][v] };
                    weight += fGet_Distance_Sqr(Pix_3_A, Pix_3_B, 3) + weight_eps;
                }
                if (!oImage_A.m_pChannel[3][i] || !oImage_A.m_pChannel[3][v] ||
                    !oImage_B.m_pChannel[3][i] || !oImage_B.m_pChannel[3][v])
                {
                    weight += bad_region_penalty_;
                }

                //(i,v)为边，边上有两个权值，不知是不是容量与流量，待考
                //pEdge[iCur_Edge++] = { i,v,weight,weight };
                Add_Edge(poGraph, i, v, weight, weight);
                //printf("(%d %d)\n", i, v);
            }

            if (y < h_1)
            {
                float weight;
                {
                    unsigned char Pix_3_A[] = { oImage_A.m_pChannel[0][i],oImage_A.m_pChannel[1][i],oImage_A.m_pChannel[2][i] },
                        Pix_3_B[] = { oImage_B.m_pChannel[0][i],oImage_B.m_pChannel[1][i],oImage_B.m_pChannel[2][i] };
                    weight = fGet_Distance_Sqr(Pix_3_A, Pix_3_B, 3);
                }

                {
                    v = i + oImage_A.m_iWidth;
                    unsigned char Pix_3_A[] = { oImage_A.m_pChannel[0][v],oImage_A.m_pChannel[1][v],oImage_A.m_pChannel[2][v] },
                        Pix_3_B[] = { oImage_B.m_pChannel[0][v],oImage_B.m_pChannel[1][v],oImage_B.m_pChannel[2][v] };
                    weight += fGet_Distance_Sqr(Pix_3_A, Pix_3_B, 3) + weight_eps;
                }
                if (!oImage_A.m_pChannel[3][i] || !oImage_A.m_pChannel[3][v] ||
                    !oImage_B.m_pChannel[3][i] || !oImage_B.m_pChannel[3][v])
                    weight += bad_region_penalty_;
                //加边
                Add_Edge(poGraph, i, v, weight, weight);
                //printf("(%d %d)\n", i, v);
                //pEdge[iCur_Edge++] = { i,v,weight,weight };
            }
        }
    }
    return;
}
void Max_Flow(GCGraph* poGraph)
{//先干一遍，还不知如何与图像匹配有关
    const int TERMINAL = -1, ORPHAN = -2;
    GCGraph oG = *poGraph;
    Graph_Vertex stub = {}, * nilNode = &stub, * first = nilNode, * last = nilNode;
    stub.next = nilNode;

    //Vtx* vtxPtr = &vtcs[0];
    Graph_Edge* edgePtr = oG.m_Edge;
    Graph_Vertex* vtxPtr = oG.m_Vertex;

    int i, curr_ts = 0;
    for (i = 0; i < oG.m_iMax_Vertex_Count; i++)
    {
        Graph_Vertex* v = &oG.m_Vertex[i];
        v->ts = 0;
        if (v->weight != 0)
        {
            last = last->next = v;
            v->dist = 1;
            v->parent = TERMINAL;
            v->t = v->weight < 0;
        }
        else
            v->parent = 0;
    }

    first = first->next;
    last->next = nilNode;
    nilNode->next = 0;


    Graph_Vertex** orphans = (Graph_Vertex**)pMalloc(oG.m_iMax_Vertex_Count * sizeof(Graph_Vertex*));
    int iCur_Orphan = 0;
    int iCount = 0;

    for (;;)
    {
        Graph_Vertex* v, * u;
        int e0 = -1, ei = 0, ej = 0;
        float minWeight, weight;
        unsigned char vt;
        iCount++;
        while (first != nilNode)
        {
            v = first;
            if (v->parent)
            {
                vt = v->t;
                for (ei = v->first; ei != 0; ei = edgePtr[ei].next)
                {
                    if (edgePtr[ei ^ vt].weight == 0)
                        continue;
                    u = vtxPtr + edgePtr[ei].dst;
                    if (!u->parent)
                    {
                        u->t = vt;
                        u->parent = ei ^ 1;
                        u->ts = v->ts;
                        u->dist = v->dist + 1;
                        if (!u->next)
                        {
                            u->next = nilNode;
                            last = last->next = u;
                        }
                        continue;
                    }

                    if (u->t != vt)
                    {
                        e0 = ei ^ vt;
                        break;
                    }

                    if (u->dist > v->dist + 1 && u->ts <= v->ts)
                    {
                        // reassign the parent
                        u->parent = ei ^ 1;
                        u->ts = v->ts;
                        u->dist = v->dist + 1;
                    }
                }
                if (e0 > 0)
                    break;
            }

            // exclude the vertex from the active list
            first = first->next;
            v->next = 0;
        }//while

        if (e0 <= 0)
            break;

        // find the minimum edge weight along the path
        minWeight = edgePtr[e0].weight;

        // k = 1: source tree, k = 0: destination tree
        for (int k = 1; k >= 0; k--)
        {
            for (v = vtxPtr + edgePtr[e0 ^ k].dst;; v = vtxPtr + edgePtr[ei].dst)
            {
                if ((ei = v->parent) < 0)
                    break;
                weight = edgePtr[ei ^ k].weight;
                minWeight = Min(minWeight, weight);
            }
            weight = abs(v->weight);
            minWeight = Min(minWeight, weight);
        }

        // modify weights of the edges along the path and collect orphans
        edgePtr[e0].weight -= minWeight;
        edgePtr[e0 ^ 1].weight += minWeight;
        oG.flow += minWeight;

        // k = 1: source tree, k = 0: destination tree
        for (int k = 1; k >= 0; k--)
        {
            for (v = vtxPtr + edgePtr[e0 ^ k].dst;; v = vtxPtr + edgePtr[ei].dst)
            {
                if ((ei = v->parent) < 0)
                    break;
                edgePtr[ei ^ (k ^ 1)].weight += minWeight;
                if ((edgePtr[ei ^ k].weight -= minWeight) == 0)
                {
                    orphans[iCur_Orphan++] = v;
                    v->parent = ORPHAN;
                }
            }

            v->weight = v->weight + minWeight * (1 - k * 2);
            if (v->weight == 0)
            {
                orphans[iCur_Orphan++] = v;
                v->parent = ORPHAN;
            }
        }//for

        // restore the search trees by finding new parents for the orphans
        curr_ts++;
        while (iCur_Orphan)
        {
            Graph_Vertex* v2 = orphans[--iCur_Orphan];

            int d, minDist = INT_MAX;
            e0 = 0;
            vt = v2->t;

            for (ei = v2->first; ei != 0; ei = edgePtr[ei].next)
            {
                if (edgePtr[ei ^ (vt ^ 1)].weight == 0)
                    continue;
                u = vtxPtr + edgePtr[ei].dst;
                if (u->t != vt || u->parent == 0)
                    continue;
                // compute the distance to the tree root
                for (d = 0;; )
                {
                    if (u->ts == curr_ts)
                    {
                        d += u->dist;
                        break;
                    }
                    ej = u->parent;
                    d++;
                    if (ej < 0)
                    {
                        if (ej == ORPHAN)
                            d = INT_MAX - 1;
                        else
                        {
                            u->ts = curr_ts;
                            u->dist = 1;
                        }
                        break;
                    }
                    u = vtxPtr + edgePtr[ej].dst;
                }

                // update the distance
                if (++d < INT_MAX)
                {
                    if (d < minDist)
                    {
                        minDist = d;
                        e0 = ei;
                    }
                    for (u = vtxPtr + edgePtr[ei].dst; u->ts != curr_ts; u = vtxPtr + edgePtr[u->parent].dst)
                    {
                        u->ts = curr_ts;
                        u->dist = --d;
                    }
                }
            }

            if ((v2->parent = e0) > 0)
            {
                v2->ts = curr_ts;
                v2->dist = minDist;
                continue;
            }

            /* no parent is found */
            v2->ts = 0;
            for (ei = v2->first; ei != 0; ei = edgePtr[ei].next)
            {
                u = vtxPtr + edgePtr[ei].dst;
                ej = u->parent;
                if (u->t != vt || !ej)
                    continue;
                if (edgePtr[ei ^ (vt ^ 1)].weight && !u->next)
                {
                    u->next = nilNode;
                    last = last->next = u;
                }
                if (ej > 0 && vtxPtr + edgePtr[ej].dst == v2)
                {
                    orphans[iCur_Orphan++] = u;
                    u->parent = ORPHAN;
                }
            }

        }//while
    }
    Free(orphans);
    return;
}

int bIn_Source_Segment(GCGraph* poGraph, int i)
{
    return poGraph->m_Vertex[i].t == 0;
}

template<typename _T>void Find_Pair_1(Stitch<_T>* poStitch, int* dx[], int* dy[], Image_Pair** ppPair, Image_Pair** ppImage_Pair_GPU, int* piPair_Count)
{//用内存搞，数据源为 Image_Warp 的图像与 Mask数据
    int iSize, iMax_Size, iPair_Count = 0, iMax_Pair_Count = 10,
        iImage_Count = poStitch->m_iImage_Count;
    const int gap = 10;
    unsigned char* p;
    //现在CPU上找到交集数
    Light_Ptr oPtr;
    iSize = 100000000;
    Attach_Light_Ptr(oPtr, (unsigned char*)pMalloc(iSize), iSize,0);
    Image_Pair *pImage_Pair;

    iSize = iMax_Pair_Count * sizeof(Image_Pair);
    Malloc(oPtr, iSize, p);
    pImage_Pair = (Image_Pair*)p;

    //先将Image_Warp从GPU  拷贝到内存中
    iSize = iImage_Count * sizeof(Image);
    Malloc(oPtr, iSize, p);
    Image* pImage_Warp = (Image*)p;

    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_0.bmp", &poStitch->m_pImage_Warp[0]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_1.bmp", &poStitch->m_pImage_Warp[1]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_2.bmp", &poStitch->m_pImage_Warp[2]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_3.bmp", &poStitch->m_pImage_Warp[3]);
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        int w = poStitch->m_pImage_Warp[i].m_iWidth,
            h = poStitch->m_pImage_Warp[i].m_iHeight;
        Malloc(oPtr, w * h * 4, p);
        Attach_Buffer(&pImage_Warp[i], p, w, h, 4, Image::IMAGE_TYPE_BMP);
        //Copy_Image_To_CPU(poStitch->m_pImage_Warp[i],pImage_Warp[i]);
        //bSave_Image_GPU("c:\\tmp\\1.bmp", poStitch->m_pImage_Warp[0]);
        pImage_Warp[i].m_pChannel[0] = poStitch->m_pImage_Warp[i].m_pChannel[0];
        pImage_Warp[i].m_pChannel[1] = poStitch->m_pImage_Warp[i].m_pChannel[1];
        pImage_Warp[i].m_pChannel[2] = poStitch->m_pImage_Warp[i].m_pChannel[2];
        cudaMemcpy(pImage_Warp[i].m_pChannel[3], poStitch->m_pImage_Warp[i].m_pChannel[3],
            pImage_Warp[i].m_iWidth * pImage_Warp[i].m_iHeight, cudaMemcpyDeviceToHost);
        //Copy_Image_To_CPU(poStitch->m_pMasks_Warped[i], pImage_Mask[i]);
    }
    
    //此处可以装入Warp Image以便对数据
    
    iMax_Size = 0;
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        for (int j = i + 1; j < poStitch->m_iImage_Count; j++)
        {
            Image_Pair oPair;
            if (bIs_Overlap(poStitch->m_pCorner[i][0], poStitch->m_pCorner[j][0],
                poStitch->m_pImage_Warp[i].m_iWidth, poStitch->m_pImage_Warp[i].m_iHeight,
                poStitch->m_pImage_Warp[j].m_iWidth, poStitch->m_pImage_Warp[j].m_iHeight, oPair.roi))
            {
                if (iPair_Count >= iMax_Pair_Count)
                {
                    printf("Pair count exceed:%d in Find_Pair\n", iMax_Pair_Count);
                    return;
                }

                oPair.m_iImage_A = i, oPair.m_iImage_B = j;
                //Init_Image_GPU(&oPair.m_oImage_A_GPU, oPair.roi[1][0] + 2 * gap, oPair.roi[1][1] + 2 * gap, Image::IMAGE_TYPE_BMP, 32, &oPtr_GPU);
                //Init_Image_GPU(&oPair.m_oImage_B_GPU, oPair.m_oImage_A_GPU.m_iWidth, oPair.m_oImage_A_GPU.m_iHeight, Image::IMAGE_TYPE_BMP, oPair.m_oImage_A_GPU.m_iBit_Count, &oPtr_GPU);
                int w = oPair.roi[1][0] + 2 * gap,
                    h = oPair.roi[1][1] + 2 * gap;
                Malloc(oPtr, w * h*4, p);
                Attach_Buffer(&oPair.m_oImage_A, p, w, h, 4, Image::IMAGE_TYPE_BMP);
                Malloc(oPtr, w * h*4, p);
                Attach_Buffer(&oPair.m_oImage_B, p, w, h, 4, Image::IMAGE_TYPE_BMP);

                if (oPair.m_oImage_A.m_iWidth * oPair.m_oImage_A.m_iHeight > iMax_Size)
                    iMax_Size = oPair.m_oImage_A.m_iWidth * oPair.m_oImage_A.m_iHeight;

                pImage_Pair[iPair_Count] = oPair;
                iPair_Count++;
            }
        }
    }
    //Disp_Cuda_Error();

    /*dim3 oThread, oGrid;
    oThread.x = 512;
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = iPair_Count;*/

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    //Find_Pair_Copy_Image_1 << <oGrid, oThread >> > (poStitch->m_pImage_Warp_Header_GPU, poStitch->m_pCorner_GPU, gap, pImage_Pair);

    //Disp_Cuda_Error();
    //Image_Pair oPair = pImage_Pair[0];
    //bSave_Image_GPU("c:\\tmp\\1.bmp", oPair.m_oImage_A);
    //bSave_Image_GPU("c:\\tmp\\2.bmp", oPair.m_oImage_B);

    for (int i = 0; i < iPair_Count; i++)
    {
        Image_Pair oPair = pImage_Pair[i];
        dim3 oThread, oGrid;
        iSize = oPair.m_oImage_A.m_iHeight * oPair.m_oImage_A.m_iWidth;
        oThread.x = Min(512, iSize);
        oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;

        Find_Pair_Copy_Image_1 << <oGrid, oThread >> > (pImage_Warp, poStitch->m_pCorner_GPU, gap, oPair);
        Disp_Cuda_Error();

        GCGraph oGraph;
        Set_Graph_Weights_Color(oPair.m_oImage_A, oPair.m_oImage_B, &oGraph);
        //unsigned long long tStart = iGet_Tick_Count();
        Max_Flow(&oGraph);
        //printf("%lld\n", iGet_Tick_Count() - tStart);
        //printf("flow:%f\n", oGraph.flow);

        int* tl_A = poStitch->m_pCorner[oPair.m_iImage_A][0],
            * tl_B = poStitch->m_pCorner[oPair.m_iImage_B][0];
        Image oImage_A = pImage_Warp[oPair.m_iImage_A],
            oImage_B = pImage_Warp[oPair.m_iImage_B];
        unsigned char* pMask_A = oImage_A.m_pChannel[3],
            * pMask_B = oImage_B.m_pChannel[3];
        
        for (int y = 0; y < oPair.roi[1][1]; y++)
        {
            for (int x = 0; x < oPair.roi[1][0]; x++)
            {
                if (bIn_Source_Segment(&oGraph, (y + gap) * (oPair.roi[1][0] + 2 * gap) + x + gap))
                {
                    if (pMask_A[(oPair.roi[0][1] - tl_A[1] + y) * oImage_A.m_iWidth + oPair.roi[0][0] - tl_A[0] + x])
                        pMask_B[(oPair.roi[0][1] - tl_B[1] + y) * oImage_B.m_iWidth + oPair.roi[0][0] - tl_B[0] + x] = 0;
                }
                else
                {
                    if (pMask_B[(oPair.roi[0][1] - tl_B[1] + y) * oImage_B.m_iWidth + oPair.roi[0][0] - tl_B[0] + x])
                        pMask_A[(oPair.roi[0][1] - tl_A[1] + y) * oImage_A.m_iWidth + oPair.roi[0][0] - tl_A[0] + x] = 0;
                }
            }
        }

        Free(oGraph.m_Edge);
        Free(oGraph.m_Vertex);
        //bSave_Comp("c:\\tmp\\1.bmp", pImage_Warp[1],3);
    }
    /*bSave_Comp("c:\\tmp\\1.bmp", pImage_Warp[0],3);
    bSave_Comp("c:\\tmp\\2.bmp", pImage_Warp[1],3);
    bSave_Comp("c:\\tmp\\3.bmp", pImage_Warp[2],3);
    bSave_Comp("c:\\tmp\\4.bmp", pImage_Warp[3],3);*/


    //将 Image_Warp的Mask
    for (int i = 0; i < iImage_Count; i++)
    {
        cudaMemcpy(poStitch->m_pImage_Warp[i].m_pChannel[3], pImage_Warp[i].m_pChannel[3],
            pImage_Warp[i].m_iWidth * pImage_Warp[i].m_iHeight, cudaMemcpyHostToDevice);
    }
    Free(oPtr.m_pBuffer);
    return;
}
template<typename _T>void Find_Pair(Stitch<_T>* poStitch, int* dx[], int* dy[], Image_Pair** ppPair, Image_Pair** ppImage_Pair_GPU, int* piPair_Count)
{
    int iSize, iMax_Size, iPair_Count = 0, iMax_Pair_Count = 10;
    const int gap = 10;
    unsigned char* p;

    //现在CPU上找到交集数
    Light_Ptr oPtr, oPtr_GPU;
    iSize = 100000000;
    Attach_Light_Ptr(oPtr_GPU, (unsigned char*)pMalloc_GPU(iSize), iSize,0);
    Image_Pair* pImage_Pair_GPU, * pImage_Pair = (Image_Pair*)pMalloc(iMax_Pair_Count * sizeof(Image_Pair));

    iMax_Size = 0;
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        for (int j = i + 1; j < poStitch->m_iImage_Count; j++)
        {
            Image_Pair oPair;
            if (bIs_Overlap(poStitch->m_pCorner[i][0], poStitch->m_pCorner[j][0],
                poStitch->m_pImage_Warp[i].m_iWidth, poStitch->m_pImage_Warp[i].m_iHeight,
                poStitch->m_pImage_Warp[j].m_iWidth, poStitch->m_pImage_Warp[j].m_iHeight, oPair.roi))
            {
                if (iPair_Count >= iMax_Pair_Count)
                {
                    printf("Pair count exceed:%d in Find_Pair\n", iMax_Pair_Count);
                    return;
                }

                oPair.m_iImage_A = i, oPair.m_iImage_B = j;
                Init_Image_GPU(&oPair.m_oImage_A_GPU, oPair.roi[1][0] + 2 * gap, oPair.roi[1][1] + 2 * gap, Image::IMAGE_TYPE_BMP, 32, &oPtr_GPU);
                Init_Image_GPU(&oPair.m_oImage_B_GPU, oPair.m_oImage_A_GPU.m_iWidth, oPair.m_oImage_A_GPU.m_iHeight, Image::IMAGE_TYPE_BMP, oPair.m_oImage_A_GPU.m_iBit_Count, &oPtr_GPU);

                iSize = oPair.m_oImage_A_GPU.m_iWidth * oPair.m_oImage_A_GPU.m_iHeight * sizeof(int);
                Malloc(oPtr_GPU, iSize, p);
                oPair.m_pSub_dx_A = (int*)p;
                Malloc(oPtr_GPU, iSize, p);
                oPair.m_pSub_dy_A = (int*)p;
                Malloc(oPtr_GPU, iSize, p);
                oPair.m_pSub_dx_B = (int*)p;
                Malloc(oPtr_GPU, iSize, p);
                oPair.m_pSub_dy_B = (int*)p;

                if (oPair.m_oImage_A_GPU.m_iWidth * oPair.m_oImage_A_GPU.m_iHeight > iMax_Size)
                    iMax_Size = oPair.m_oImage_A_GPU.m_iWidth * oPair.m_oImage_A_GPU.m_iHeight;

                pImage_Pair[iPair_Count] = oPair;
                iPair_Count++;
                //打印一下可以看出是否重叠过度，一张图同时和同一方向上多图重叠
                //printf("%d %d\n", i, j);
            }
        }
    }

    Shrink(pImage_Pair, iPair_Count * sizeof(Image_Pair));

    iSize = iPair_Count * sizeof(Image_Pair);
    Malloc(oPtr_GPU, iSize, p);
    pImage_Pair_GPU = (Image_Pair*)p;
    cudaMemcpy(pImage_Pair_GPU, pImage_Pair, iSize, cudaMemcpyHostToDevice);

    dim3 oThread, oGrid;
    oThread.x = 512;
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = iPair_Count;

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Find_Pair_Copy_Image << <oGrid, oThread >> > (poStitch->m_pImage_Warp_Header_GPU, poStitch->m_pCorner_GPU, gap,
        dx, dy, pImage_Pair_GPU);

    Attach_Light_Ptr(oPtr, (unsigned char*)pMalloc(oPtr_GPU.m_iCur), oPtr_GPU.m_iCur, 0);
    cudaMemcpy(oPtr.m_pBuffer, oPtr_GPU.m_pBuffer, oPtr_GPU.m_iCur, cudaMemcpyDeviceToHost);
    
    
    Image* pImage_Mask = (Image*)pMalloc(poStitch->m_iImage_Count * sizeof(Image));
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        Init_Image(&pImage_Mask[i], poStitch->m_pMasks_Warped[i].m_iWidth,
            poStitch->m_pMasks_Warped[i].m_iHeight, Image::IMAGE_TYPE_BMP, 8);
        Copy_Image_To_CPU(poStitch->m_pMasks_Warped[i], pImage_Mask[i]);
    }
    Disp_Cuda_Error();
    //bSave_Image("c:\\tmp\\1.bmp", pImage_Mask[0]);

    //此处补充个最大流i就完活，暂时没工夫去改进，只能保留Opencv分离代码
    for (int i = 0; i < iPair_Count; i++)
    {
        GCGraph oGraph;
        Image_Pair oPair = pImage_Pair[i];
        int w = oPair.m_oImage_A_GPU.m_iWidth,
            h = oPair.m_oImage_B_GPU.m_iHeight;
        iSize = w * h * 4;
        p = oPtr.m_pBuffer + (oPair.m_oImage_A_GPU.m_pChannel[0] - oPtr_GPU.m_pBuffer);
        Attach_Buffer(&oPair.m_oImage_A, p, w, h, 4, Image::IMAGE_TYPE_BMP);
        p = oPtr.m_pBuffer + (oPair.m_oImage_B_GPU.m_pChannel[0] - oPtr_GPU.m_pBuffer);
        Attach_Buffer(&oPair.m_oImage_B, p, w, h, 4, Image::IMAGE_TYPE_BMP);

        /*bSave_Image_GPU("c:\\tmp\\1.bmp", oPair.m_oImage_A_GPU);
        bSave_Image("c:\\tmp\\2.bmp", oPair.m_oImage_A);
        Compare_Image("c:\\tmp\\1.bmp", "c:\\tmp\\2.bmp");

        bSave_Image_GPU("c:\\tmp\\1.bmp", oPair.m_oImage_B_GPU);
        bSave_Image("c:\\tmp\\2.bmp", oPair.m_oImage_B);
        Compare_Image("c:\\tmp\\1.bmp", "c:\\tmp\\2.bmp");*/
                
        Set_Graph_Weights_Color(oPair.m_oImage_A, oPair.m_oImage_B, &oGraph);
        Max_Flow(&oGraph);
        printf("flow:%f\n", oGraph.flow);

        //Free_Image(&oPair.m_oImage_A);
        //Free_Image(&oPair.m_oImage_B);

        //Image oMask_A = poStitch->m_pMasks_Warped[i];

        //最后求得两个Image的Mask
        Image oMask_A = pImage_Mask[oPair.m_iImage_A],
            oMask_B = pImage_Mask[oPair.m_iImage_B];
        for (int y = 0; y < oPair.roi[1][1]; y++)
        {
            for (int x = 0; x < oPair.roi[1][0]; x++)
            {
                /*if (bIn_Source_Segment(&oGraph, (y + gap) * (roi[1][0] + 2 * gap) + x + gap))
                {
                    if (oMask_A.m_pChannel[0][(roi[0][1] - tl_A[1] + y) * oMask_A.m_iWidth + roi[0][0] - tl_A[0] + x])
                        oMask_B.m_pChannel[0][(roi[0][1] - tl_B[1] + y) * oMask_B.m_iWidth + roi[0][0] - tl_B[0] + x] = 0;
                }
                else
                {
                    if (oMask_B.m_pChannel[0][(roi[0][1] - tl_B[1] + y) * oMask_B.m_iWidth + roi[0][0] - tl_B[0] + x])
                        oMask_A.m_pChannel[0][(roi[0][1] - tl_A[1] + y) * oMask_A.m_iWidth + roi[0][0] - tl_A[0] + x] = 0;
                }*/
            }
        }

        Free(oGraph.m_Edge);
        Free(oGraph.m_Vertex);
    }

    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //Disp_Cuda_Error();
    //Image_Pair oPair = pImage_Pair[0];
    //Disp_Part(oPair.m_pSub_dx_B, oPair.m_oImage_B.m_iWidth, 0, 50, oPair.m_oImage_B.m_iWidth, 1, "dx");

    //for(int i=0;i<iPair_Count;i++)
    //{
    //    Image_Pair oPair = pImage_Pair[i];
    //    Disp_Sum_GPU(oPair.m_oImage_A.m_pChannel[3], oPair.m_oImage_A.m_iWidth * oPair.m_oImage_A.m_iHeight);
    //    bSave_Comp_GPU("c:\\tmp\\2.bmp", oPair.m_oImage_A, 3);
    //    //Disp_Sum_GPU(oPair.m_pSub_dy_B, oPair.m_oImage_B.m_iWidth * oPair.m_oImage_B.m_iHeight);
    //}
    *ppPair = pImage_Pair;
    *piPair_Count = iPair_Count;
    /*Light_Ptr oPtr;
    Light_Ptr oPtr_GPU;
    int iSize = iMax_Pair_Count* sizeof(Image_Pair) +
        iMax_Pair_Count
    Attach_Light_Ptr(oPtr,(unsigned char*)pMalloc())*/
    //Image(*pImage_Pair)[2] = (Image(*)[2])pMalloc(iPair_Count * 2 * sizeof(Image));
    //Image(*pImage_Pair_GPU)[2] = (Image(*)[2])pMalloc_GPU(iPair_Count * 2 * sizeof(Image));

    return;
}

template<typename _T>__global__ void Set_Sobel_Mid(Image Warp[], int iImage_Count,
    _T* (*Sobel)[2][3], _T* (*Mid)[2][3], int** dx, int** dy, _T* p)
{
    int i, j, k;
    int* p1 = (int*)p;

    for (i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oImage = Warp[i].m_oPart_1;
        int iSize = oImage.m_iHeight * oImage.m_iWidth;
        dx[i] = p1;
        p1 += iSize;
        p1 = (int*)ALIGN_SIZE_128(p1);
        dy[i] = p1;
        p1 += iSize;
        p1 = (int*)ALIGN_SIZE_128(p1);
    }
    _T* p2 = (_T*)ALIGN_SIZE_128(p1);

    for (i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oImage = Warp[i].m_oPart_1;
        int iSize = oImage.m_iHeight * oImage.m_iWidth;
        for (j = 0; j < 2; j++)
            for (k = 0; k < 3; k++)
                Sobel[i][j][k] = p2 + (j * 3 + k) * iSize;
        p2 += iSize * 6;

        for (j = 0; j < 2; j++)
            for (k = 0; k < 3; k++)
                Mid[i][j][k] = p2 + (j * 3 + k) * iSize;
        p2 += iSize * 6;
    }
}

template<typename Dest_Type>__global__ void _Get_dx_dy(Image Warp[],
    Pixel_4 oKer_y, Pixel_4 oKer_x, int iKernel_Size, Dest_Type* Mid[][2][3], Dest_Type* B[][2][3],
    int** dx, int** dy, Border_Type iBorder_Type)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image::Part_1 oImage;
    if (threadIdx.x == 0)
        oImage = Warp[blockIdx.y].m_oPart_1;

    __syncthreads();
    if (iThread_ID >= oImage.m_iWidth * oImage.m_iHeight)
        return;         //多余线程退出
    short x = iThread_ID % oImage.m_iWidth,
        y = iThread_ID / oImage.m_iWidth;

    short r = iKernel_Size >> 1;
    Dest_Type Total[2][3] = {};
    int iPos;
    for (short i = -r; i <= r; i++)
    {//此处恐怕要估计Border_Type
        short y1 = iGet_Border_y_GPU(y + i, oImage.m_iHeight, iBorder_Type, iThread_ID);
        iPos = y1 * oImage.m_iWidth + x;
        Total[0][0] += Mid[blockIdx.y][0][0][iPos] * oKer_y.Data_c[i + r];
        Total[0][1] += Mid[blockIdx.y][0][1][iPos] * oKer_y.Data_c[i + r];
        Total[0][2] += Mid[blockIdx.y][0][2][iPos] * oKer_y.Data_c[i + r];

        Total[1][0] += Mid[blockIdx.y][1][0][iPos] * oKer_x.Data_c[i + r];
        Total[1][1] += Mid[blockIdx.y][1][1][iPos] * oKer_x.Data_c[i + r];
        Total[1][2] += Mid[blockIdx.y][1][2][iPos] * oKer_x.Data_c[i + r];

        /*if (blockIdx.y == 0 && y == 0 && x == 3)
            printf("y+i:%d y1:%d x:%d Mid:%d\n",y+i, y1, x, Mid[blockIdx.y][0][2][iPos]);*/
    }
    /*B[blockIdx.y][0][0][iThread_ID] = Total[0][0];
    B[blockIdx.y][0][1][iThread_ID] = Total[0][1];
    B[blockIdx.y][0][2][iThread_ID] = Total[0][2];*/

    //捋一捋dx,dy的数值范围
    //Mid 为[-1,0,1], [1,2,1] 的卷积结果，Max(Mid)为 255*(1+2+1) = 1020
    //Max(dx) = 1028*1020 * 3 = 3,121,200   只能用int表示
    dx[blockIdx.y][iThread_ID] =
        Total[0][0] * Total[0][0] +
        Total[0][1] * Total[0][1] +
        Total[0][2] * Total[0][2];

    /*B[blockIdx.y][1][0][iThread_ID] = Total[1][0];
    B[blockIdx.y][1][1][iThread_ID] = Total[1][1];
    B[blockIdx.y][1][2][iThread_ID] = Total[1][2];*/

    dy[blockIdx.y][iThread_ID] =
        Total[1][0] * Total[1][0] +
        Total[1][1] * Total[1][1] +
        Total[1][2] * Total[1][2];

    return;
}

template<typename Dest_Type>__global__ void _Sep_Filter_row(Image Warp[],
    Pixel_4 oKer_x, Pixel_4 oKer_y, int iKernel_Size, Dest_Type* Mid[][2][3], Border_Type iBorder_Type)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image::Part_1 oImage;
    if (threadIdx.x == 0)
    {
        oImage = Warp[blockIdx.y].m_oPart_1;
    }
    __syncthreads();
    if (iThread_ID >= oImage.m_iWidth * oImage.m_iHeight)
        return;         //多余线程退出
    short x = iThread_ID % oImage.m_iWidth,
        y = iThread_ID / oImage.m_iWidth;

    /*if (blockIdx.y == 0 && y == oImage.m_iHeight - 1 && x == 1)
        printf("Here");*/

    short r = iKernel_Size >> 1;
    Dest_Type Total[2][3] = {};
    int iPos;
    for (short i = -r; i <= r; i++)
    {//此处恐怕要估计Border_Type
        short x1 = iGet_Border_x_GPU(x + i, oImage.m_iWidth, iBorder_Type);
        iPos = y * oImage.m_iWidth + x1;
        Total[0][0] += oImage.m_pChannel[0][iPos] * oKer_x.Data_c[i + r];
        Total[0][1] += oImage.m_pChannel[1][iPos] * oKer_x.Data_c[i + r];
        Total[0][2] += oImage.m_pChannel[2][iPos] * oKer_x.Data_c[i + r];

        Total[1][0] += oImage.m_pChannel[0][iPos] * oKer_y.Data_c[i + r];
        Total[1][1] += oImage.m_pChannel[1][iPos] * oKer_y.Data_c[i + r];
        Total[1][2] += oImage.m_pChannel[2][iPos] * oKer_y.Data_c[i + r];

        //if (blockIdx.y == 0 && y == 0 && x == 3)
            //printf("y:%d x:%d x1:%d Pix:%d Mid:%d\n", y, x, x1, oImage.m_pChannel[2][iPos], Total[0][2]);
    }

    Mid[blockIdx.y][0][0][iThread_ID] = Total[0][0];
    Mid[blockIdx.y][0][1][iThread_ID] = Total[0][1];
    Mid[blockIdx.y][0][2][iThread_ID] = Total[0][2];

    Mid[blockIdx.y][1][0][iThread_ID] = Total[1][0];
    Mid[blockIdx.y][1][1][iThread_ID] = Total[1][1];
    Mid[blockIdx.y][1][2][iThread_ID] = Total[1][2];

    /*if (Mid[0][0][2][3] != 7)
        printf("%d ", Mid[0][0][2][3]);*/

}

template<typename Dest_Type>static void Sep_Filter_2D_2(Image Warp[], int iImage_Count, int iMax_Size,
    Pixel_4 oKer_x, Pixel_4 oKer_y, int iKernel_Size, Dest_Type* B[][2][3], Dest_Type* Mid[][2][3], int** dx, int** dy)
{//此处已经完全对准opencv的数据，将多步合而为一

    dim3 oThread, oGrid;
    oThread.x = 512;
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = iImage_Count;

    _Sep_Filter_row << <oGrid, oThread >> > (Warp, oKer_x, oKer_y, iKernel_Size, Mid, BORDER_REFLECT101);
    _Get_dx_dy << <oGrid, oThread >> > (Warp, oKer_y, oKer_x, iKernel_Size, Mid, B, dx, dy, BORDER_REFLECT101);

    //Image oImage = Warp[0];
    ////Disp_Sum_GPU<short>(B[0][0][2], oImage.m_iWidth * oImage.m_iHeight);
    ////Disp_Part_GPU<short>(B[0][0][2], oImage.m_iWidth, 0, 0, oImage.m_iWidth, 1, "Sobel");
    //Disp_Part_GPU<short>(&Mid[0][0][2][0], oImage.m_iWidth, 0, 0, oImage.m_iWidth, 1, "Mid");
    //Disp_Part_GPU(oImage.m_pChannel[2], oImage.m_iWidth, 2, 0, 3, 2, "r");
    return;
}

//template<typename Source_Type, typename Dest_Type, typename Kernel_Type>static void Sep_Filter_2D_1
//(Source_Type A[], int w, int h, Kernel_Type Ker_x[], Kernel_Type Ker_y[], int iKernel_Size, Dest_Type B[])
//{
//    Dest_Type* pMid = (Dest_Type*)pMalloc(w * h * sizeof(Dest_Type));
//
//    Sep_Filter_2D_Line(A, w, h, Ker_x, iKernel_Size, pMid);
//    //Disp(pMid, w,h);
//
//    //Matrix_Transpose(pMid, w, h, pMid);
//    //Disp(pMid, 1,w);
//
//    //出来以后，pMid 是一个宽为h, 高为w的转置矩阵，继续行推进即可
//    Sep_Filter_2D_Line(pMid, h, w, Ker_y, iKernel_Size, B);
//    //Disp(B, 1, w);
//
//    Free(pMid);
//}

template<typename _T>void Find(Stitch<_T>* poStitch)
{//注意，其实这版并不需要dx,dy，在做得过程中顺手算了，而且耗时很少
    typedef short Sobel_Type;
    //第一步，给Sobel 分配内存
    int i, iMax_Size = 0, iSize, iSize_All = 0;
    Sobel_Type* (*pSobel)[2][3];   //就是放Warp图经过可分离卷积之后的16位像素值
    Sobel_Type* (*pMid)[2][3];    //中间结果
    int** dx, ** dy;

    iSize = poStitch->m_iImage_Count * 3 * 2 * sizeof(Sobel_Type*) +    //Sobel
        poStitch->m_iImage_Count * 3 * 2 * sizeof(int*) +               //mid
        poStitch->m_iImage_Count * sizeof(int*);                        //dx,dy

    dx = (int**)pMalloc_GPU(iSize);
    dy = dx + poStitch->m_iImage_Count;
    pSobel = (Sobel_Type * (*)[2][3])(dy + poStitch->m_iImage_Count);
    pMid = (Sobel_Type * (*)[2][3])pSobel + poStitch->m_iImage_Count;

    for (i = 0; i < poStitch->m_iImage_Count; i++)
    {
        Image::Part_1 oImage = poStitch->m_pImage_Warp[i].m_oPart_1;
        iSize = oImage.m_iHeight * oImage.m_iWidth;
        if (iSize > iMax_Size)
            iMax_Size = iSize;
        iSize_All += iSize;
    }

    iSize = iSize_All * 3 * 2 * sizeof(Sobel_Type) +    //Soble
        iSize_All * 3 * 2 * sizeof(Sobel_Type) +        //Mid
        +iSize_All * 2 * sizeof(int) + 128 * 12;          //dx,dy

    Sobel_Type* p = (Sobel_Type*)pMalloc_GPU(iSize);
    poStitch->m_pdx_dy_GPU = (int*)p;
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Set_Sobel_Mid << <1, 1 >> > (poStitch->m_pImage_Warp_Header_GPU,
        poStitch->m_iImage_Count, pSobel, pMid, dx, dy, p);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);
    //exit(0);


    //全部图一起做一个可分离卷积
    Pixel_4 oKernel_x;  // = { -1, 0, 1 };
    Pixel_4     oKernel_y = { 1, 2, 1 };
    oKernel_x.Data_c[0] = -1, oKernel_x.Data_c[1] = 0, oKernel_x.Data_c[2] = 1;
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Sep_Filter_2D_2<Sobel_Type>(poStitch->m_pImage_Warp_Header_GPU, poStitch->m_iImage_Count,
        iMax_Size, oKernel_x, oKernel_y, 3, pSobel, pMid, dx, dy);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //Disp_Mem_GPU();
    //Free_GPU(pSobel[0][0][0]);
    Shrink_GPU(poStitch->m_pdx_dy_GPU, iSize_All * 2 * sizeof(int) + 128 * 2);

    //Disp_Mem_GPU();
    //Disp_Cuda_Error();
    //int** dx_1, ** dy_1;
    //dx_1 = (int**)pMalloc(poStitch->m_iImage_Count * 2 * sizeof(int*));
    //dy_1 = dx_1 + poStitch->m_iImage_Count;
    //cudaMemcpy(dx_1, dx, poStitch->m_iImage_Count * 2 * sizeof(int*),cudaMemcpyDeviceToHost);
    //for(int i=0;i< poStitch->m_iImage_Count;i++)
    //{
    //    Image oImage;
    //    cudaMemcpy(&oImage, &poStitch->m_pImage_Warp[i], sizeof(Image), cudaMemcpyDeviceToHost);

    //    //for (int j = 0; j < 2; j++)
    //    //{
    //    //    /*for (int k = 0; k < 3; k++)
    //    //        Disp_Sum_GPU<short>(B[i][j][k], oImage.m_iWidth * oImage.m_iHeight);*/
    //    //    Disp_Sum_GPU<short>(B[i][j][0], oImage.m_iWidth * oImage.m_iHeight*3);
    //    //}
    //    Disp_Sum_GPU(dy_1[i], oImage.m_iWidth * oImage.m_iHeight);
    //    printf("\n");
    //}
    //bSave_Image_GPU("c:\\tmp\\1.bmp", poStitch->m_pImage_Warp[0]);
    //后面要find_Pair
    Image_Pair* pPair, * pPair_GPU;
    int iPair_Count;    Disp_Cuda_Error();
    unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<100;i++)
    {
        Find_Pair_1<_T>(poStitch, dx, dy, &pPair, &pPair_GPU, &iPair_Count);
        //Find_Pair<_T>(poStitch, dx, dy, &pPair, &pPair_GPU, &iPair_Count);
        //Free(pPair);
        //Free(pPair[0].m_oImage_A.m_pBuffer);
    }
    Disp_Cuda_Error();
    printf("Max Pool: %lld ms\n", iGet_Tick_Count() - tStart);
    //Free_GPU(pPair[0].m_oImage_A_GPU.m_pBuffer);
    //Free(pPair);
    Free_GPU(poStitch->m_pdx_dy_GPU);
    poStitch->m_pdx_dy_GPU = NULL;
    Free_GPU(dx);
    return;
}

template<typename _T>void Warp_3(Image oImage, _T K[3 * 3], _T R[3 * 3], _T fScale,
    Image* poImage_Warped, Interpolation_Flag iImage_Inter_Type, Interpolation_Flag iMask_Inter_Type,
    int Corner[2][2], Border_Type iImage_Border_Type, Border_Type iMask_Border_Type, Point_Cloud<float>* poPC)
{//还是尝试用原来得方法，一直做到Blender Prepare
    //_T* puxmap, * puymap;    //暂时未知
    int Dest_roi[2][2];   //[0][0-1]: x,y [1][0-1]: w,h
    Sphere_Projector<_T> oProjector = {};

    Build_Map<_T>(oImage.m_iWidth, oImage.m_iHeight, K, R, fScale, /*&puxmap, &puymap,*/ Dest_roi, &oProjector, poPC);
    //注意，由于前面得warp用的roi大小不是严格得，故此此处要修正
    Init_Image_GPU(poImage_Warped, Dest_roi[1][0] + 1, Dest_roi[1][1] + 1, Image::IMAGE_TYPE_BMP, 32);
    Re_Map_3_GPU<_T>(oImage, *poImage_Warped, Dest_roi[0][0], Dest_roi[0][1], &oProjector,
        iImage_Inter_Type, iMask_Inter_Type, iImage_Border_Type, iMask_Border_Type);
    //Disp((int*)Dest_roi,2,2,"roi");

    //bSave_Image_GPU("c:\\tmp\\1.bmp", oImage);
    //bSave_Comp_GPU("c:\\tmp\\2.bmp", *poImage_Warped,3);
    if (Corner)
    {
        Corner[0][0] = Dest_roi[0][0];
        Corner[0][1] = Dest_roi[0][1];
        Corner[1][0] = Dest_roi[0][0] + Dest_roi[1][0];
        Corner[1][1] = Dest_roi[0][1] + Dest_roi[1][1];
    }
    //Disp((int*)Dest_roi, 2, 2);
    return;
}

__global__ void _Dilate_col_GPU(Image Warp[])
{
    __shared__ Image::Part_1 oImage;
    int iThread_ID = GET_THREAD_ID();

    if (threadIdx.x == 0)
        oImage = Warp[blockIdx.y].m_oPart_1;
    __syncthreads();
    if (iThread_ID >= oImage.m_iWidth)
        return;
    short iMove_To, y, y1, r = 1;
    unsigned char* pCur = &oImage.m_pChannel[3][iThread_ID + oImage.m_iWidth];

    for (y = 1; y < oImage.m_iHeight; y++, pCur += oImage.m_iWidth)
    {
        if (*pCur != pCur[-oImage.m_iWidth])
        {
            if (*pCur)
            {//白色，向上
                iMove_To = y - r;
                if (iMove_To < 0)
                    iMove_To = 0;
                unsigned char* pPrev = pCur - oImage.m_iWidth;
                for (y1 = y - 1; y1 >= iMove_To && !*pPrev; y1--)
                    *pPrev = 0xFF, pPrev -= oImage.m_iWidth;
            }
            else
            {//黑色，向下
                iMove_To = y + r;
                if (iMove_To > oImage.m_iHeight)
                    iMove_To = oImage.m_iHeight;
                unsigned char* pNext = pCur;
                for (y1 = y; y1 < iMove_To && !(*pNext); y1++, pNext += oImage.m_iWidth)
                    *pNext = 0xFF;
                y = y1;
                pCur = &oImage.m_pChannel[3][y * oImage.m_iWidth + iThread_ID];
            }
        }
    }
}

__global__ void _Dilate_row_GPU(Image Warp[])
{//慢到死
    __shared__ Image::Part_1 oImage;
    int iThread_ID = GET_THREAD_ID();

    if (threadIdx.x == 0)
        oImage = Warp[blockIdx.y].m_oPart_1;
    __syncthreads();

    if (iThread_ID >= oImage.m_iHeight)
        return;

    short iMove_To, x1, r = 1;
    //short w_8 = (oImage.m_iWidth >> 3) << 3;
    unsigned char* pLine = &oImage.m_pChannel[3][iThread_ID * oImage.m_iWidth];

    for (int x = 1; x < oImage.m_iWidth;)
    {
        if (pLine[x] != pLine[x - 1])
        {//有变化，要处理
            if (pLine[x])
            {//白色，向后扩大
                iMove_To = x - r;
                if (iMove_To < 0)
                    iMove_To = 0;
                for (x1 = x - 1; x1 >= iMove_To && !pLine[x1]; x1--)
                    pLine[x1] = 0xFF;
            }
            else
            {//黑色，向前推进
                iMove_To = x + r;
                if (iMove_To > oImage.m_iWidth)
                    iMove_To = oImage.m_iWidth;
                for (x1 = x; x1 < iMove_To && !pLine[x1]; x1++)
                    pLine[x1] = 0xFF;
                x = x1;
            }
        }
        x++;
    }
}

__global__ void _Dilate_row_GPU_1(Image Warp[])
{
    __shared__ Image::Part_1 oImage;
    if (threadIdx.x == 0)
    {
        oImage = Warp[blockIdx.y].m_oPart_1;
    }
    __syncthreads();
    short x = threadIdx.x,
        y = blockIdx.x;

    if (y >= oImage.m_iHeight || x >= oImage.m_iWidth)
        return;

    Pixel_4 oPix;
    oPix = *(Pixel_4*)(&oImage.m_pChannel[3][y * oImage.m_iWidth + x - 1]);
    if (x == 0)
        oPix.Data[0] = oPix.Data[1];
    if (x == oImage.m_iWidth - 1)
        oPix.Data[2] = oPix.Data[1];
    oPix.Data[3] = oPix.Data[1];    //备份一下

    if (oPix.Data[1] == 0)
    {
        if (oPix.Data[0] == 0xFF || oPix.Data[2] == 0xFF)
            oPix.Data[1] = 0xFF;
    }
    __syncthreads();
    if (oPix.Data[1] != oPix.Data[3])
        oImage.m_pChannel[3][y * oImage.m_iWidth + x] = 0xFF;
    //if (y == 1 && x == 126 && blockIdx.y == 0)
        //printf("%d %d %d %d\n", oPix.Data[0], oPix.Data[1], oPix.Data[2], oPix.Data[3]);

}

void Dilate_GPU(Image Warp[], Image Warp_Header_GPU[], int iImage_Count)
{//只对Image的alpha 通道进行
//这个函数已经完全对准 opencv

    dim3 oGrid, oThread;
    const int iThread_Per_Block = 1024;
    int iMax_w = 0, iMax_h = 0;
    for (int i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oImage = Warp[i].m_oPart_1;
        if (oImage.m_iWidth > iMax_w)
            iMax_w = oImage.m_iWidth;
        if (oImage.m_iHeight > iMax_h)
            iMax_h = oImage.m_iHeight;
    }

    if (iMax_h > iThread_Per_Block)
    {
        //先搞行方向
        oThread.x = Min(iThread_Per_Block, iMax_h);
        oGrid.x = (iMax_h + oThread.x - 1) / oThread.x;
        oGrid.y = iImage_Count;
        _Dilate_row_GPU << <oGrid, oThread >> > (Warp_Header_GPU);
    }
    else
    {
        oThread.x = Min(iThread_Per_Block, iMax_w);
        oGrid.x = iMax_h;
        oGrid.y = iImage_Count;
        _Dilate_row_GPU_1 << <oGrid, oThread >> > (Warp_Header_GPU);
    }

    oThread.x = Min(iThread_Per_Block, iMax_w);
    oGrid.x = (iMax_w + oThread.x - 1) / oThread.x;
    oGrid.y = iImage_Count;
    _Dilate_col_GPU << <oGrid, oThread >> > (Warp_Header_GPU);
    return;
}

__global__ void _Resize_Bitwise_And(Image Source[], Image Dest[], Border_Type iBorder_Type)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image oSource, oDest;
    __shared__ float f_x, f_y;
    if (threadIdx.x == 0)
    {
        oSource = Source[blockIdx.y];
        oDest = Dest[blockIdx.y];
        f_x = (float)oSource.m_iWidth / oDest.m_iWidth;
        f_y = (float)oSource.m_iHeight / oDest.m_iHeight;
    }
    __syncthreads();

    if (iThread_ID >= oDest.m_iWidth * oDest.m_iHeight)
        return;
    short x_d = iThread_ID % oDest.m_iWidth,
        y_d = iThread_ID / oDest.m_iWidth;

    float w2, w3, x_s_f = (x_d + 0.5f) * f_x - 0.5f;
    unsigned char* pCur_Line, * pNext_Line;

    {
        float y_s_f = (y_d + 0.5f) * f_y - 0.5f;
        int y_s_0 = (int)floor(y_s_f);
        int y_s_1 = y_s_0 + 1;
        w2 = (y_s_1 - y_s_f), w3 = 1.f - w2;
        pCur_Line = &oSource.m_pChannel[3][iGet_Border_y_GPU(y_s_0, oSource.m_iHeight, iBorder_Type) * oSource.m_iWidth];
        pNext_Line = &oSource.m_pChannel[3][iGet_Border_y_GPU(y_s_1, oSource.m_iHeight, iBorder_Type) * oSource.m_iWidth];
    }

    float w0, w1;
    int x_s_0_r, x_s_1_r;
    {
        int x_s_0 = (int)floor(x_s_f);
        int x_s_1 = x_s_0 + 1;
        w0 = (x_s_1 - x_s_f), w1 = 1.f - w0;
        x_s_0_r = iGet_Border_x_GPU(x_s_0, oSource.m_iWidth, iBorder_Type);
        x_s_1_r = iGet_Border_x_GPU(x_s_1, oSource.m_iWidth, iBorder_Type);
    }

    float fValue_0, fValue_1;
    fValue_0 = w0 * pCur_Line[x_s_0_r] + w1 * pCur_Line[x_s_1_r];
    fValue_1 = w0 * pNext_Line[x_s_0_r] + w1 * pNext_Line[x_s_1_r];

    oDest.m_pChannel[3][iThread_ID] &= (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5);

    return;
}

template<typename _T>void Resize_Bitwise_And(Stitch<_T>* poStitch,
    Image Blend_Warp[], Image Blend_Warp_Header_GPU[])
{
    int iMax_Size = 0;
    Image::Part_1 oImage;
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        oImage = Blend_Warp[i].m_oPart_1;
        if (oImage.m_iWidth * oImage.m_iHeight > iMax_Size)
            iMax_Size = oImage.m_iWidth * oImage.m_iHeight;
    }

    dim3 oThread, oGrid;
    oThread.x = 512;
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = poStitch->m_iImage_Count;
    //bSave_Comp_GPU("c:\\tmp\\2.bmp", Blend_Warp[2], 3);
    //Compare_Image("c:\\tmp\\1.bmp", "c:\\tmp\\2.bmp", 1);
    _Resize_Bitwise_And << <oGrid, oThread >> > (poStitch->m_pImage_Warp_Header_GPU,
        Blend_Warp_Header_GPU, BORDER_REFLECT);

    /*Disp_Cuda_Error();
    bSave_Comp_GPU("c:\\tmp\\2.bmp", Blend_Warp[2], 3);
    Compare_Image("c:\\tmp\\1.bmp", "c:\\tmp\\2.bmp",1);*/
    return;
}

template<typename _T>void Temp_Compare(const char File[], _T* pBuffer, int w, int h, _T iDiff_Threshold = 0)
{
    FILE* pFile = fopen(File, "rb");
    if (!pFile)
    {
        printf("Fail to open:%s\n", File);
        return;
    }
    _T* pBuffer_1 = (_T*)pMalloc(w * h * sizeof(_T));
    Disp_Cuda_Error();
    cudaMemcpy(pBuffer_1, pBuffer, w * h * sizeof(_T), cudaMemcpyDeviceToHost);
    const float eps = 0.000001f;
    int iCount = 0;
    for (int y = 0; y < h; y++)
        for (int x = 0; x < w; x++)
        {
            //if (y == 1 && x == 384)
                //printf("here");
            if (typeid(_T) == typeid(float))
            {
                float fValue;
                int iResult = (int)fread(&fValue, 1, 4, pFile);
                if (abs(fValue - pBuffer_1[y * w + x]) > eps)
                {
                    printf("y:%d x:%d Source:%f Dest:%f\n", y, x, fValue, pBuffer_1[y * w + x]);
                    iCount++;
                    return;
                }
            }
            else if (typeid(_T) == typeid(short))
            {
                short iValue;
                int iResult = (int)fread(&iValue, 1, 2, pFile);
                if (abs(iValue - pBuffer_1[y * w + x]) > iDiff_Threshold)
                {
                    printf("y:%d x:%d Source:%d Dest:%d\n", y, x, iValue, pBuffer_1[y * w + x]);
                    iCount++;
                    //return;
                }
            }
        }

    Free(pBuffer_1);
    fclose(pFile);
    printf("Mismatched Count:%d\n", iCount);

    return;
}
template<typename _T>void Temp_Load_Camera(const char* pcFile, Stitch<_T>* poStitch)
{
    Stitch<_T> oStitch = *poStitch;
    Match_Item<_T> oMatch = {};
    int i, j, iResult, iSize, bRet = 0,
        iMatch_Count;

    FILE* pFile = fopen(pcFile, "rb");
    if (!pFile)
    {
        printf("Fail to load:%s in Temp_Load_Camera\n", pcFile);
        exit(0);
        goto END;
    }
       

    //读入iCamera_Count
    iResult = (int)fread(&oStitch.m_iImage_Count, 1, 4, pFile);
    if (!iResult || !oStitch.m_iImage_Count)
        goto END;

    iSize = oStitch.m_iImage_Count * sizeof(Camera<_T>);
    oStitch.m_pCamera = (Camera<_T>*)pMalloc(iSize);
    memset(oStitch.m_pCamera, 1, oStitch.m_iImage_Count * sizeof(Camera<_T>));
    for (i = 0; i < oStitch.m_iImage_Count; i++)
    {
        //先读K
        for (j = 0; j < 3 * 3; j++)
        {
            double fValue;
            iResult = (int)fread(&fValue, 1, sizeof(double), pFile);
            oStitch.m_pCamera[i].K[j] = fValue;
        }

        for (j = 0; j < 3 * 3; j++)
        {
            double fValue;
            iResult = (int)fread(&fValue, 1, sizeof(fValue), pFile);
            oStitch.m_pCamera[i].R[j] = fValue;
        }

        for (j = 0; j < 3; j++)
        {
            double fValue;
            iResult = (int)fread(&fValue, 1, sizeof(fValue), pFile);
            oStitch.m_pCamera[i].t[j] = fValue;
        }
    }

    iResult = (int)fread(&iMatch_Count, 1, 4, pFile);
    oStitch.m_pImage_Match = (Match_Item<_T>*)pMalloc(iMatch_Count * sizeof(Match_Item<_T>));
    for (i = 0; i < iMatch_Count; i++)
    {
        int A, B;
        iResult = (int)fread(&A, 1, sizeof(A), pFile);
        iResult = (int)fread(&B, 1, sizeof(B), pFile);
        oMatch.m_iImage_A = A;
        oMatch.m_iImage_B = B;
        oStitch.m_pImage_Match[i] = oMatch;
    }
    *poStitch = oStitch;
    bRet = 1;

END:
    if (!bRet)
        Free_Stitch(&oStitch);
    if (pFile)
        fclose(pFile);
    return;
}

template<typename _T>void bLoad_Image(Stitch<_T> oStitch)
{
    for (int i = 0; i < oStitch.m_iImage_Count; i++)
    {
        char File[256];
        sprintf(File, "data\\%d.bmp", i);
        if (!bLoad_Image_GPU(File, &oStitch.m_pImage_Source[i]))
            exit(0);
    }
    return;
}

//template<typename _T>void Temp_Load_Image(Stitch<_T>* poStitch)
//{
//    char File[256];
//    //FILE* pFile;
//    int i, iSize;
//    Stitch<_T> oStitch = *poStitch;
//}

template<typename _T>void Set_K(Stitch<_T>* poStitch,
    _T K_s[3 * 3], _T K[3 * 3], _T fAspect)
{
    /*_T fSeam_Work_Aspect = poStitch->seam_work_aspect;
    K[0] = K_s[0] * fSeam_Work_Aspect;
    K[2] = K_s[2] * fSeam_Work_Aspect;
    K[4] = K_s[4] * fSeam_Work_Aspect;
    K[5] = K_s[5] * fSeam_Work_Aspect;*/

    K[0] = K_s[0] * fAspect;
    K[2] = K_s[2] * fAspect;
    K[4] = K_s[4] * fAspect;
    K[5] = K_s[5] * fAspect;
    K[8] = 1.f;
    K[1] = K[3] = K[6] = K[7] = 0.f;
    return;
}


template<typename _T>void Free_Partial(Stitch<_T>* poStitch)
{//在Feed之后释放部分内存
    //理应释放三部分内存，Seam_Est，Image_Warped，Mask
    Free_Image_GPU(&poStitch->m_pSeam_Est[0]);  //此处连Mask也释放了
    //Free_Image_GPU(&poStitch->m_pImage_Warp[0]);
    //以下三句做不做两可，纯粹为了显式
    memset(poStitch->m_pSeam_Est, 0, poStitch->m_iImage_Count * sizeof(Image));
    //memset(poStitch->m_pImage_Warp, 0, poStitch->m_iImage_Count * sizeof(Image));
    //memset(poStitch->m_pMask, 0, poStitch->m_iImage_Count * sizeof(Image));
    //以下三句做不做两可，纯粹为了显式
    cudaMemset(poStitch->m_pSeam_Est_Header_GPU, 0, poStitch->m_iImage_Count * sizeof(Image));
    //cudaMemset(poStitch->m_pImage_Warp_Header_GPU, 0, poStitch->m_iImage_Count * sizeof(Image));
    //cudaMemset(poStitch->m_pMask_Header_GPU, 0, poStitch->m_iImage_Count * sizeof(Image));
}


void Result_Roi(int Corner[][2][2], int Size[][2], int iImage_Count, int roi[2][2])
{//感觉就是求个bounding box
    int tl[2] = { std::numeric_limits<int>::max(), std::numeric_limits<int>::max() },
        br[2] = { std::numeric_limits<int>::min(), std::numeric_limits<int>::min() };
    for (int i = 0; i < iImage_Count; i++)
    {
        if (Corner[i][0][0] < tl[0])
            tl[0] = Corner[i][0][0];
        if (Corner[i][0][1] < tl[1])
            tl[1] = Corner[i][0][1];

        if (Corner[i][0][0] + Size[i][0] > br[0])
            br[0] = Corner[i][0][0] + Size[i][0];
        if (Corner[i][0][1] + Size[i][1] > br[1])
            br[1] = Corner[i][0][1] + Size[i][1];
    }
    roi[0][0] = tl[0], roi[0][1] = tl[1],
        roi[1][0] = br[0] - tl[0], roi[1][1] = br[1] - tl[1];
    return;
}
__global__ void Disp_Image(Image::Part_1 Group_1[], Image::Part_1 Group_2[])
{
    int i;
    for (i = 0; i < 6; i++)
        printf("%d %d\n", Group_1[i].m_iHeight, Group_1[i].m_iWidth);
    for (i = 0; i < 6; i++)
        printf("%d %d\n", Group_2[i].m_iHeight, Group_2[i].m_iWidth);
    return;
}
template<typename _T>void Prepare_Blender_1(Stitch<_T>* poStitch, Blender* poBlender)
{
    int roi[2][2];
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        poStitch->m_pSize[i][0] = poStitch->m_pCorner[i][1][0] - poStitch->m_pCorner[i][0][0] + 1;
        poStitch->m_pSize[i][1] = poStitch->m_pCorner[i][1][1] - poStitch->m_pCorner[i][0][1] + 1;
    }
    //Disp((int*)poStitch->m_pSize, 4, 2);
    Result_Roi(poStitch->m_pCorner, poStitch->m_pSize, poStitch->m_iImage_Count, roi);
    memcpy(poBlender->dst_roi_final, roi, sizeof(roi));

    double max_len = Max(roi[1][0], roi[1][1]);
    //注意，log是以e为底
    //log(max_len)/log(2) = log2(max_len) 看这个max_len能占多少位
    //真正的目的是求划分block的最大值32，如果max_len还不到32，就一个block
    int num_bands_ = std::min(5, static_cast<int>(ceil(std::log(max_len) / std::log(2.0))));

    int iBlock_Size = 1 << num_bands_;
    roi[1][0] += (iBlock_Size - roi[1][0] % iBlock_Size) % iBlock_Size;
    roi[1][1] += (iBlock_Size - roi[1][1] % iBlock_Size) % iBlock_Size;

    //此时，roi的[1]装着目标图的大小。目标图还得建金字塔
    memcpy(poBlender->dst_roi, roi, sizeof(roi));
    poBlender->m_iNum_Band = num_bands_;

    //计算使用空间
    //首先Image头的大小
    //开辟内存
    int iTotal = 0, iSize = roi[1][0] * roi[1][1];      //Mask?
    //可以先算一下需要总量
    iTotal += 6 * sizeof(short*);       //dst_pyr_laplace_Header_GPU
    iTotal += 6 * sizeof(float*);        //dst_band_weights
    iTotal += iSize * 3 * sizeof(short) + 128; //m_pDst 23003136
    iTotal += iSize + 128;                    //oMask; 3833856
    iTotal += iSize * sizeof(float) + 128;    //dst_band_weights   15335424 
    //再算金字塔的空间需要
    {
        int Size_1[] = { roi[1][0],roi[1][1] };
        for (int i = 1; i <= num_bands_; ++i)
        {
            Size_1[0] = (Size_1[0] + 1) >> 1;
            Size_1[1] = (Size_1[1] + 1) >> 1;
            iSize = Size_1[0] * Size_1[1];
            iTotal += iSize * 3 * sizeof(short) + 128;    //dst_pyr_laplace[i]
            iTotal += iSize * sizeof(float) + 128;        //dst_band_weights
        }
    }

    Light_Ptr oPtr;
    unsigned char* p;
    Attach_Light_Ptr(oPtr, (unsigned char*)pMalloc_GPU(iTotal), iTotal, 0);
    //次吃必要吗？
    cudaMemset(oPtr.m_pBuffer, 0, oPtr.m_iMax_Buffer_Size);
    poBlender->m_pBuffer = oPtr.m_pBuffer;

    //先分配dst_pyr_laplace
    iSize = roi[1][0] * roi[1][1];
    Image::Part_1 oImage;
    oImage.m_iHeight = roi[1][1];
    oImage.m_iWidth = roi[1][0];
    Malloc(oPtr, iSize * 3 * sizeof(short), oImage.m_pChannel[0]);
    oImage.m_pChannel[1] = oImage.m_pChannel[0] + iSize * sizeof(short);
    oImage.m_pChannel[2] = oImage.m_pChannel[1] + iSize * sizeof(short);
    poBlender->dst_pyr_laplace[0] = oImage;
    for (int i = 1; i < 6; i++)
    {
        Image::Part_1 oUpper, oLower;
        oUpper = poBlender->dst_pyr_laplace[i - 1];
        oLower.m_iHeight = (oUpper.m_iHeight + 1) >> 1;
        oLower.m_iWidth = (oUpper.m_iWidth + 1) >> 1;
        iSize = oLower.m_iHeight * oLower.m_iWidth;
        Malloc(oPtr, iSize * 3 * sizeof(short), oLower.m_pChannel[0]);
        oLower.m_pChannel[1] = oLower.m_pChannel[0] + iSize * sizeof(short);
        oLower.m_pChannel[2] = oLower.m_pChannel[1] + iSize * sizeof(short);
        poBlender->dst_pyr_laplace[i] = oLower;
    }

    //再分配
    iSize = roi[1][0] * roi[1][1];
    oImage.m_iHeight = roi[1][1];
    oImage.m_iWidth = roi[1][0];
    Malloc(oPtr, iSize * sizeof(float), oImage.m_pChannel[0]);
    poBlender->dst_band_weights[0] = oImage;
    for (int i = 1; i < 6; i++)
    {
        Image::Part_1 oUpper, oLower;
        oUpper = poBlender->dst_band_weights[i - 1];
        oLower.m_iHeight = (oUpper.m_iHeight + 1) >> 1;
        oLower.m_iWidth = (oUpper.m_iWidth + 1) >> 1;
        iSize = oLower.m_iHeight * oLower.m_iWidth;
        Malloc(oPtr, iSize * sizeof(float), oLower.m_pChannel[0]);
        poBlender->dst_band_weights[i] = oLower;
    }

    Malloc(oPtr, 12 * sizeof(Image::Part_1), p);
    poBlender->dst_pyr_laplace_Header_GPU = (Image::Part_1*)p;
    poBlender->dst_band_weights_Header_GPU = poBlender->dst_pyr_laplace_Header_GPU + 6;
    cudaMemcpy(poBlender->dst_pyr_laplace_Header_GPU, poBlender->dst_pyr_laplace, 12 * sizeof(Image::Part_1), cudaMemcpyHostToDevice);

    return;
}

//template<typename _T>void Prepare_Blender(Stitch<_T> *poStitch,Blender* poBlender)
//{//这个废弃
//    int roi[2][2];
//    for (int i = 0; i < poStitch->m_iImage_Count; i++)
//    {
//        poStitch->m_pSize[i][0] = poStitch->m_pCorner[i][1][0] - poStitch->m_pCorner[i][0][0] + 1;
//        poStitch->m_pSize[i][1] = poStitch->m_pCorner[i][1][1] - poStitch->m_pCorner[i][0][1] + 1;
//    }
//    //Disp((int*)poStitch->m_pSize, 4, 2);
//    Result_Roi(poStitch->m_pCorner, poStitch->m_pSize, poStitch->m_iImage_Count, roi);
//    
//    double max_len = Max(roi[1][0], roi[1][1]);
//    //注意，log是以e为底
//    //log(max_len)/log(2) = log2(max_len) 看这个max_len能占多少位
//    //真正的目的是求划分block的最大值32，如果max_len还不到32，就一个block
//    int num_bands_ = std::min(5, static_cast<int>(ceil(std::log(max_len) / std::log(2.0))));
//
//    int iBlock_Size = 1 << num_bands_;
//    roi[1][0] += (iBlock_Size - roi[1][0] % iBlock_Size) % iBlock_Size;
//    roi[1][1] += (iBlock_Size - roi[1][1] % iBlock_Size) % iBlock_Size;
//
//    memcpy(poBlender->dst_roi, roi, sizeof(roi));
//    poBlender->m_iNum_Band = num_bands_;
//
//    //开辟内存
//    int iTotal=0,iSize = roi[1][0] * roi[1][1];
//    //可以先算一下需要总量
//    iTotal += 6 * sizeof(short*);       //dst_pyr_laplace_Header_GPU
//    iTotal += 6 * sizeof(float*);        //dst_band_weights
//    iTotal += iSize * 3 * sizeof(short) + 128; //m_pDst 23003136
//    iTotal += iSize + 128;                    //oMask; 3833856
//    iTotal += iSize * sizeof(float) + 128;    //dst_band_weights   15335424 
//    {
//        int Size_1[] = { roi[1][0],roi[1][1] };
//        for (int i = 1; i <= num_bands_; ++i)
//        {
//            Size_1[0] = (Size_1[0] + 1) >> 1;
//            Size_1[1] = (Size_1[1] + 1) >> 1;
//            iSize = Size_1[0] * Size_1[1];
//            iTotal += iSize * 3 * sizeof(short) + 128;    //dst_pyr_laplace[i]
//            iTotal += iSize * sizeof(float) + 128;        //dst_band_weights
//        }
//    }
//    Light_Ptr oPtr;
//    unsigned char* p;
//    Attach_Light_Ptr(oPtr, (unsigned char*)pMalloc_GPU(iTotal), iTotal, 0);
//    cudaMemset(oPtr.m_pBuffer, 0, oPtr.m_iMax_Buffer_Size);
//    poBlender->m_pBuffer = oPtr.m_pBuffer;
//
//    //正式开始分配
//    iSize = 6 * sizeof(short*) + 6 * sizeof(float*);
//    Malloc(oPtr, iSize, p);
//    poBlender->dst_pyr_laplace_Header_GPU = (short**)p;
//    poBlender->dst_band_weights_Header_GPU =(float**)(poBlender->dst_pyr_laplace_Header_GPU + 6);
//
//    //m_pDst
//    iSize = roi[1][0] * roi[1][1];
//    Malloc(oPtr, iSize * 3 * sizeof(short), p);    
//    poBlender->m_pDst = (short(*)[3])p;
//    poBlender->dst_pyr_laplace[0] = (short*)poBlender->m_pDst;
//    //cudaMemset(poBlender->m_pDst, 0, iSize * 3 * sizeof(short));
//        
//    Init_Image_GPU(&poBlender->m_oMask, roi[1][0], roi[1][1], Image::IMAGE_TYPE_BMP, 8, &oPtr);
//    //Set_Color_GPU(poBlender->m_oMask);
//
//    //拉普拉斯金字塔开辟空间
//    Malloc(oPtr, iSize * sizeof(float), p);
//    poBlender->dst_band_weights[0] = (float*)p;
//    //cudaMemset(poBlender->dst_band_weights[0], 0, iSize * sizeof(float));
//    {
//        int Size_1[] = { roi[1][0],roi[1][1] };
//        for (int i = 1; i <= num_bands_; ++i)
//        {
//            Size_1[0] = (Size_1[0] + 1) >> 1;
//            Size_1[1] = (Size_1[1] + 1) >> 1;
//            iSize = Size_1[0] * Size_1[1];
//            Malloc(oPtr, iSize * 3 * sizeof(short), p);
//            poBlender->dst_pyr_laplace[i] = (short*)p;
//            //cudaMemset(poBlender->dst_pyr_laplace[i], 0, iSize * 3 * sizeof(short));
//
//            Malloc(oPtr, iSize * sizeof(float), p);
//            poBlender->dst_band_weights[i] = (float*)p;
//            //cudaMemset(poBlender->dst_band_weights[i], 0, iSize * sizeof(float));
//        }
//    }
//    cudaMemcpy(poBlender->dst_pyr_laplace_Header_GPU, poBlender->dst_pyr_laplace, 6 * 2 * 8, cudaMemcpyHostToDevice);
//
//    return;
//}

template<typename _T>void Feed_Blender_Get_Pos(Image Warp[], Stitch<_T>* poStitch, Blender* poBlender, int Pos[][2][2], int Size[][2], short LTRB_New[][2][2])
{
    int gap = 3 * (1 << poBlender->m_iNum_Band);
    for (int i = 0; i < poStitch->m_iImage_Count; i++)
    {
        Image::Part_1 oImage = Warp[i].m_oPart_1;
        int* tl = poStitch->m_pCorner[i][0];
        int tl_new[2] = { std::max(poBlender->dst_roi[0][0], tl[0] - gap),
                      std::max(poBlender->dst_roi[0][1], tl[1] - gap) };
        int dst_roi_br[2] = { poBlender->dst_roi[0][0] + poBlender->dst_roi[1][0],
                      poBlender->dst_roi[0][1] + poBlender->dst_roi[1][1] };
        int br_new[2] = { std::min(dst_roi_br[0] , tl[0] + oImage.m_iWidth + gap),
            std::min(dst_roi_br[1], tl[1] + oImage.m_iHeight + gap) };

        tl_new[0] = poBlender->dst_roi[0][0] + (((tl_new[0] - poBlender->dst_roi[0][0]) >> poBlender->m_iNum_Band) << poBlender->m_iNum_Band);
        tl_new[1] = poBlender->dst_roi[0][1] + (((tl_new[1] - poBlender->dst_roi[0][1]) >> poBlender->m_iNum_Band) << poBlender->m_iNum_Band);

        int width = br_new[0] - tl_new[0];
        int height = br_new[1] - tl_new[1];

        width += ((1 << poBlender->m_iNum_Band) - width % (1 << poBlender->m_iNum_Band)) % (1 << poBlender->m_iNum_Band);
        height += ((1 << poBlender->m_iNum_Band) - height % (1 << poBlender->m_iNum_Band)) % (1 << poBlender->m_iNum_Band);

        br_new[0] = tl_new[0] + width;
        br_new[1] = tl_new[1] + height;

        int dy = std::max(br_new[1] - dst_roi_br[1], 0);
        int dx = std::max(br_new[0] - dst_roi_br[0], 0);

        tl_new[0] -= dx; br_new[0] -= dx;
        tl_new[1] -= dy; br_new[1] -= dy;

        int top = tl[1] - tl_new[1];
        int left = tl[0] - tl_new[0];
        int bottom = br_new[1] - tl[1] - oImage.m_iHeight;
        int right = br_new[0] - tl[0] - oImage.m_iWidth;
        //printf("%d %d %d %d\n", top, left, bottom, right);
        Pos[i][0][0] = left + 5;
        Pos[i][0][1] = top;
        Pos[i][1][0] = right;
        Pos[i][1][1] = bottom;
        Size[i][0] = oImage.m_iWidth + Pos[i][0][0] + Pos[i][1][0];
        Size[i][1] = oImage.m_iHeight + Pos[i][0][1] + Pos[i][1][1];

        LTRB_New[i][0][0] = tl_new[0];
        LTRB_New[i][0][1] = tl_new[1];
        LTRB_New[i][1][0] = br_new[0];
        LTRB_New[i][1][1] = br_new[1];

        //if (i == 3)
            //printf("Here");
    }
    //Disp((int*)Pos, 4, 4, "Pos");
    return;
}
__global__ void Copy_Make_Border_1(Image Source[], Image Dest[], int Border[][2][2])
{//重写一个，妄图快点
    __shared__ Image oSource, oDest;
    __shared__ short iTop, iLeft, iRight, iBottom;;
    if (threadIdx.x == 0)
    {
        oSource = Source[blockIdx.y];
        oDest = Dest[blockIdx.y];
        iLeft = Border[blockIdx.y][0][0];
        iTop = Border[blockIdx.y][0][1];
        iRight = Border[blockIdx.y][1][0];
        iBottom = Border[blockIdx.y][1][1];
    }
    __syncthreads();
    short iWidth_Align_4 = (oSource.m_iWidth + 3) >> 2;
    if (threadIdx.x >= iWidth_Align_4 || blockIdx.x >= oSource.m_iHeight)
        return;             //没用的线程走人

    int iDest_Size = oDest.m_iWidth * oDest.m_iHeight;
    short x_s = threadIdx.x * 4;      //本线程对应的x坐标
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
            if (iRemain_x < 4)
            {
                for (short i = 0; i < iRemain_x; i++)
                    pDest[i] = pSource[i];
            }
            else
                *(Pixel_4*)pDest = *(Pixel_4*)pSource;
        }
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
        //此处有误，应该把自己这这行抄到目标上，而不是自己扮演了目标行
        //unsigned char* pDest = &oDest.m_pChannel[0][blockIdx.x * oDest.m_iWidth],
        //    * pSource = &oDest.m_pChannel[0][((iTop << 1) - blockIdx.x - 1) * oDest.m_iWidth];
        unsigned char* pSource = &oDest.m_pChannel[0][y_d * oDest.m_iWidth];
        unsigned char* pDest = &oDest.m_pChannel[0][(iTop - blockIdx.x - 1) * oDest.m_iWidth];

        /*if (blockIdx.y == 2 && threadIdx.x == 0 && blockIdx.x == 26)
            printf("%d\n", iTop - blockIdx.x - 1);*/

        int iDist = oDest.m_iWidth - (iWidth_Div_4 << 2);
        for (short j = 0; j < oDest.m_iChannel_Count; j++, pSource += iDest_Size, pDest += iDest_Size)
        {
            for (short x = threadIdx.x; x < iWidth_Div_4; x += iWidth_Align_4)
            {
                *(Pixel_4*)&pDest[x * 4] = *(Pixel_4*)&pSource[x * 4];
            }

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
            for (short x = threadIdx.x; x < iWidth_Div_4; x += iWidth_Align_4)
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
    return;
}

__global__ void Copy_Make_Border(Image Source[], Image Dest[],
    int Border[][2][2])   //顺手建立了Pytamid
{//Border_REFLECT

    int iThread_ID = GET_THREAD_ID();
    __shared__ Image::Part_1 oSource, oDest;
    __shared__ int iTop, iLeft;
    if (threadIdx.x == 0)
    {
        oSource = Source[blockIdx.y].m_oPart_1;
        oDest = Dest[blockIdx.y].m_oPart_1;
        iLeft = Border[blockIdx.y][0][0];
        iTop = Border[blockIdx.y][0][1];
        //iRight = Border[blockIdx.y][1][0];
        //iBottom = Border[blockIdx.y][1][1];
    }
    __syncthreads();

    if (iThread_ID >= oDest.m_iWidth * oDest.m_iHeight)
        return;
    short y_d = iThread_ID / oDest.m_iWidth,
        x_d = iThread_ID % oDest.m_iWidth;
    short x_s, y_s;
    if (x_d < iLeft)
        x_s = iLeft - x_d - 1;
    else if (x_d >= oSource.m_iWidth + iLeft)
        x_s = oSource.m_iWidth - 1 - (x_d - (oSource.m_iWidth + iLeft));
    else
        x_s = x_d - iLeft;

    if (y_d < iTop)
        y_s = iTop - y_d - 1;
    else if (y_d >= oSource.m_iHeight + iTop)   //注意，此处只能用 Reflect,不用101
        y_s = oSource.m_iHeight - 1 - (y_d - (oSource.m_iHeight + iTop));
    else
        y_s = y_d - iTop;

    //if (blockIdx.y == 0 && y_d == 0 && x_d == 82)
    //{
    //    printf("y_s:%d x_s:%d\n", y_s, x_s);
    //    //printf("%d\n", y_d >= oSource.m_iHeight + iTop);
    //}

    int iPos_d = y_d * oDest.m_iWidth + x_d,
        iPos_s = y_s * oSource.m_iWidth + x_s;
    oDest.m_pChannel[0][iPos_d] = oSource.m_pChannel[0][iPos_s];
    oDest.m_pChannel[1][iPos_d] = oSource.m_pChannel[1][iPos_s];
    oDest.m_pChannel[2][iPos_d] = oSource.m_pChannel[2][iPos_s];

    //if (blockIdx.y == 0 && iThread_ID < 6)
    //{
    //    //printf("%d Thread:%d \n", gridDim.y,iThread_ID);
    //    Pyramid[iThread_ID] = &Dest[gridDim.y * iThread_ID];
    //}
    return;
}
__global__ void _Pyr_Up_col_Add(Image::Part_1 oMid, Image::Part_1 oDest, int iChannel_Count = 3)
{
    int iThread_ID = GET_THREAD_ID();
    int iMid_Size = oMid.m_iHeight * oMid.m_iWidth;
    if (iThread_ID >= iMid_Size)
        return;

    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;
    unsigned char bHas_Remain_y = oDest.m_iHeight > oMid.m_iHeight * 2 && y == oMid.m_iHeight - 1 ? 1 : 0;
    unsigned char bEven = (y * 2 + 1 < oDest.m_iHeight);

    //重复计算的全部提取出来，快不少
    short* r0 = &((short*)oMid.m_pChannel[0])[iGet_Border_y_GPU(y - 1, oMid.m_iHeight, BORDER_REFLECT101) * oDest.m_iWidth + x],
        * r1 = &((short*)oMid.m_pChannel[0])[iGet_Border_y_GPU(y, oMid.m_iHeight, BORDER_REFLECT) * oDest.m_iWidth + x],
        * r2 = &((short*)oMid.m_pChannel[0])[iGet_Border_y_GPU(y + 1, oMid.m_iHeight, BORDER_REFLECT) * oDest.m_iWidth + x];

    int iPos_d = (y * 2) * oDest.m_iWidth + x;

    for (int i = 0; i < iChannel_Count; i++, r0 += iMid_Size, r1 += iMid_Size, r2 += iMid_Size)
    {
        int iValue = ((short*)oDest.m_pChannel[i])[iPos_d] + ((*r0 + *r1 * 6 + *r2 + 32) >> 6);
        //if (y * 2 == 16 && x == 0)
            //printf("%d %d %d\n", ((short*)oDest.m_pChannel[i])[iPos_d], ((*r0 + *r1 * 6 + *r2 + 32) >> 6),iValue);

        ((short*)oDest.m_pChannel[i])[iPos_d] = iValue;
        if (bEven)
        {
            iValue = ((short*)oDest.m_pChannel[i])[iPos_d + oDest.m_iWidth] + ((((*r1 + *r2) << 2) + 32) >> 6);
            ((short*)oDest.m_pChannel[i])[iPos_d + oDest.m_iWidth] = iValue;    // Clip3(-128, 127, iValue);
            if (bHas_Remain_y)
            {
                iValue = ((short*)oDest.m_pChannel[i])[iPos_d + (oDest.m_iWidth << 1)] + (((short*)oDest.m_pChannel[i])[iPos_d]);
                ((short*)oDest.m_pChannel[i])[iPos_d + (oDest.m_iWidth << 1)] = iValue; // Clip3(-128, 127, iValue);
            }
        }
    }
    return;
}
__global__ void _Pyr_Up_col_Subtract_Batch_GPU(Image Source[], unsigned short* Mid[][3], Image Dest[], Border_Type iBorder_Type = BORDER_REFLECT101)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image oDest;
    __shared__ short iMid_Height;
    __shared__ int iMid_Size;
    if (threadIdx.x == 0)
    {
        oDest = Dest[blockIdx.y];
        iMid_Height = Source[blockIdx.y].m_iHeight;
        iMid_Size = (int)(Mid[blockIdx.y][1] - Mid[blockIdx.y][0]);
    }
    __syncthreads();
    if (iThread_ID >= iMid_Height * oDest.m_iWidth)
        return;

    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;
    unsigned char bHas_Remain_y = oDest.m_iHeight > iMid_Height * 2 && y == iMid_Height - 1 ? 1 : 0;
    unsigned char bEven = (y * 2 + 1 < oDest.m_iHeight);

    //重复计算的全部提取出来，快不少
    unsigned short* r0 = &Mid[blockIdx.y][0][iGet_Border_y_GPU(y - 1, iMid_Height, BORDER_REFLECT101) * oDest.m_iWidth + x],
        * r1 = &Mid[blockIdx.y][0][iGet_Border_y_GPU(y, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x],
        * r2 = &Mid[blockIdx.y][0][iGet_Border_y_GPU(y + 1, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x];

    int iPos_d = (y * 2) * oDest.m_iWidth + x;
    for (int i = 0; i < oDest.m_iChannel_Count; i++, r0 += iMid_Size, r1 += iMid_Size, r2 += iMid_Size)
    {
        int iValue = oDest.m_pChannel[i][iPos_d] - ((*r0 + *r1 * 6 + *r2 + 32) >> 6);
        oDest.m_pChannel[i][iPos_d] = Clip3(-128, 127, iValue);
        if (bEven)
        {
            iValue = oDest.m_pChannel[i][iPos_d + oDest.m_iWidth] - ((((*r1 + *r2) << 2) + 32) >> 6);
            oDest.m_pChannel[i][iPos_d + oDest.m_iWidth] = Clip3(-128, 127, iValue);
            if (bHas_Remain_y)
            {
                iValue = oDest.m_pChannel[i][iPos_d + (oDest.m_iWidth << 1)] - (oDest.m_pChannel[i][iPos_d]);
                oDest.m_pChannel[i][iPos_d + (oDest.m_iWidth << 1)] = Clip3(-128, 127, iValue);
            }
        }
    }
}
__global__ void _Pyr_Up_col_Batch_GPU(Image Source[], unsigned short* Mid[][3], Image Dest[], Border_Type iBorder_Type = BORDER_REFLECT101)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image oDest;
    __shared__ short iMid_Height;
    __shared__ int iMid_Size;
    if (threadIdx.x == 0)
    {
        oDest = Dest[blockIdx.y];
        iMid_Height = Source[blockIdx.y].m_iHeight;
        iMid_Size = (int)(Mid[blockIdx.y][1] - Mid[blockIdx.y][0]);
    }
    __syncthreads();
    if (iThread_ID >= iMid_Height * oDest.m_iWidth)
        return;

    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;
    unsigned char bHas_Remain_y = oDest.m_iHeight > iMid_Height * 2 && y == iMid_Height - 1 ? 1 : 0;
    unsigned char bEven = (y * 2 + 1 < oDest.m_iHeight);

    //重复计算的全部提取出来，快不少
    unsigned short* r0 = &Mid[blockIdx.y][0][iGet_Border_y_GPU(y - 1, iMid_Height, BORDER_REFLECT101) * oDest.m_iWidth + x],
        * r1 = &Mid[blockIdx.y][0][iGet_Border_y_GPU(y, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x],
        * r2 = &Mid[blockIdx.y][0][iGet_Border_y_GPU(y + 1, iMid_Height, BORDER_REFLECT) * oDest.m_iWidth + x];

    int iPos_d = (y * 2) * oDest.m_iWidth + x;
    /*if (iThread_ID == 0 && blockIdx.y == 0)
        printf("%d\n", iMid_Size);*/

        //int i = 0;
    for (int i = 0; i < oDest.m_iChannel_Count; i++, r0 += iMid_Size, r1 += iMid_Size, r2 += iMid_Size)
    {
        oDest.m_pChannel[i][iPos_d] = (*r0 + *r1 * 6 + *r2 + 32) >> 6;
        if (bEven)
        {
            //d1行 =   [(r1 + r2)*4 + 32]>>6
            oDest.m_pChannel[i][iPos_d + oDest.m_iWidth] = (((*r1 + *r2) << 2) + 32) >> 6;
            if (bHas_Remain_y)
                oDest.m_pChannel[i][iPos_d + (oDest.m_iWidth << 1)] = oDest.m_pChannel[i][iPos_d];
        }
    }
}

__global__ void _Pyr_Down_col_Batch_GPU(Image Source[], unsigned short* Mid[][3], Image Dest[], Border_Type iBorder_Type = BORDER_REFLECT101)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image oDest;
    __shared__ short iMid_Height;
    if (threadIdx.x == 0)
    {
        oDest = Dest[blockIdx.y];
        iMid_Height = Source[blockIdx.y].m_iHeight;
    }
    __syncthreads();

    //后面照抄？
    if (iThread_ID >= oDest.m_iHeight * oDest.m_iWidth)
        return;
    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;

    short y1 = y << 1;
    int iMid_Pos = y1 * oDest.m_iWidth;
    //int iMid_Size = Mid[blockIdx]
    //unsigned short* pSource = &Mid[blockIdx.y][0][x];

    int Mid_Pos[4] = { iGet_Border_y_GPU(y1 - 2, iMid_Height, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 - 1, iMid_Height, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 + 1, iMid_Height, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 + 2, iMid_Height, iBorder_Type) * oDest.m_iWidth };

    //for (int i = 0; i < oDest.m_iChannel_Count; i++/*, pSource += iMid_Size*/)
    //{
    unsigned short* pSource = &Mid[blockIdx.y][0][x];
    unsigned short iValue = pSource[iMid_Pos] * 6 +
        ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) << 2) +
        pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];
    oDest.m_pChannel[0][iThread_ID] = (iValue + 128) >> 8;
    //}

    if (oDest.m_iChannel_Count > 1)
    {
        pSource = &Mid[blockIdx.y][1][x];
        iValue = pSource[iMid_Pos] * 6 +
            ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) << 2) +
            pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];
        oDest.m_pChannel[1][iThread_ID] = (iValue + 128) >> 8;

        if (oDest.m_iBit_Count > 2)
        {
            pSource = &Mid[blockIdx.y][2][x];
            iValue = pSource[iMid_Pos] * 6 +
                ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) << 2) +
                pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];
            oDest.m_pChannel[2][iThread_ID] = (iValue + 128) >> 8;
        }
    }
    return;
}

__global__ void _Pyr_Up_row_Batch_GPU(Image Source[], Image Dest[], unsigned short* Mid[][3], Border_Type iBorder_Type = BORDER_REFLECT101)
{//此处要用到Dest的Width，故此也要带入
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image oSource;
    __shared__ short iMid_Width;
    if (threadIdx.x == 0)
    {
        oSource = Source[blockIdx.y];
        iMid_Width = Dest[blockIdx.y].m_iWidth;
    }
    __syncthreads();

    int iSize_s = oSource.m_iWidth * oSource.m_iHeight;
    if (iThread_ID >= iSize_s)
        return;

    short x = iThread_ID % oSource.m_iWidth,
        y = iThread_ID / oSource.m_iWidth;
    unsigned char bHas_Remain_x = iMid_Width > (oSource.m_iWidth << 1) && (x == oSource.m_iWidth - 1) ? 1 : 0;
    unsigned char bIs_Source_Border = (x == oSource.m_iWidth - 1);
    unsigned char bEven = bIs_Source_Border && ((x << 1) + 1 < iMid_Width);

    int iPos_m = y * iMid_Width + (x << 1),
        iPos_s = y * oSource.m_iWidth;

    for (short i = 0; i < oSource.m_iChannel_Count; i++)
    {
        unsigned short* pMid = &Mid[blockIdx.y][i][iPos_m];
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
    }
}
__global__ void _Pyr_Down_row_Batch_GPU(Image Source[], unsigned short* Mid[][3], Border_Type iBorder_Type = BORDER_REFLECT101)
{
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image oSource;
    if (threadIdx.x == 0)
        oSource = Source[blockIdx.y];
    __syncthreads();

    int iMid_Width = (oSource.m_iWidth + 1) >> 1;
    //然后后面照抄？
    if (iThread_ID >= oSource.m_iHeight * iMid_Width)
        return;

    short x = iThread_ID % iMid_Width,
        y = iThread_ID / iMid_Width;
    short x1 = x << 1;
    int iSize_s = oSource.m_iWidth * oSource.m_iHeight;
    unsigned char* pSource = &oSource.m_pChannel[0][y * oSource.m_iWidth];
    short Source_Pos[4] = { (short)iGet_Border_x_GPU(x1 - 2, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 - 1, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 + 1, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 + 2, oSource.m_iWidth, iBorder_Type) };

    Mid[blockIdx.y][0][iThread_ID] = pSource[x1] * 6 +
        ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) << 2) +
        pSource[Source_Pos[0]] + pSource[Source_Pos[3]];
    pSource += iSize_s;

    if (oSource.m_iChannel_Count > 1)
    {
        Mid[blockIdx.y][1][iThread_ID] = pSource[x1] * 6 +
            ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) << 2) +
            pSource[Source_Pos[0]] + pSource[Source_Pos[3]];
        pSource += iSize_s;

        if (oSource.m_iChannel_Count > 2)
            Mid[blockIdx.y][2][iThread_ID] = pSource[x1] * 6 +
            ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) << 2) +
            pSource[Source_Pos[0]] + pSource[Source_Pos[3]];
    }
}

void Pyr_Down_Batch_GPU(Image Source[], Image Source_Header_GPU[], int iImage_Count,
    unsigned short* pMid_Header[][3], unsigned short* pMid_Header_GPU[][3],
    Image Dest[], Image Dest_Header_GPU[])
{
    //寻找最大Mid面积
    int iSize, iMax_Size = 0;
    for (int j = 0; j < iImage_Count; j++)
    {
        Image::Part_1 oSource = Source[j].m_oPart_1;
        iSize = ((oSource.m_iWidth + 1) / 2) * oSource.m_iHeight;
        if (iSize > iMax_Size)
            iMax_Size = iSize;
    }

    //开始批处理下采样
    dim3 oThread, oGrid;
    oThread.x = Min(512, iMax_Size);
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = iImage_Count;
    _Pyr_Down_row_Batch_GPU << <oGrid, oThread >> > (Source_Header_GPU, pMid_Header_GPU);

    //寻找最大图像面积
    iMax_Size = 0;
    for (int j = 0; j < iImage_Count; j++)
    {
        Image::Part_1 oDest = Dest[j].m_oPart_1;
        iSize = oDest.m_iWidth * oDest.m_iHeight;
        if (iSize > iMax_Size)
            iMax_Size = iSize;
    }

    oThread.x = Min(512, iMax_Size);
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = iImage_Count;
    _Pyr_Down_col_Batch_GPU << <oGrid, oThread >> > (Source_Header_GPU, pMid_Header_GPU, Dest_Header_GPU);
}
void Pyr_Up_Batch_GPU(Image Source[], Image Source_Header_GPU[], int iImage_Count,
    unsigned short* pMid_Header[][3], unsigned short* pMid_Header_GPU[][3],
    Image Dest[], Image Dest_Header_GPU[])
{
    //寻找最大Mid面积
    int iSize, iMax_Size = 0;
    for (int j = 0; j < iImage_Count; j++)
    {
        Image::Part_1 oSource = Source[j].m_oPart_1;
        iSize = oSource.m_iWidth * oSource.m_iHeight;
        if (iSize > iMax_Size)
            iMax_Size = iSize;
    }

    //开始批处理下采样
    dim3 oThread, oGrid;
    oThread.x = Min(512, iMax_Size);
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = iImage_Count;
    _Pyr_Up_row_Batch_GPU << <oGrid, oThread >> > (Source_Header_GPU, Dest_Header_GPU, pMid_Header_GPU);
    //Disp_Cuda_Error();

    iMax_Size = 0;
    for (int j = 0; j < iImage_Count; j++)
    {
        Image::Part_1 oSource = Source[j].m_oPart_1,
            oDest = Dest[j].m_oPart_1;
        iSize = oDest.m_iWidth * oSource.m_iHeight;
        if (iSize > iMax_Size)
            iMax_Size = iSize;
    }

    oThread.x = Min(512, iMax_Size);
    oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    oGrid.y = iImage_Count;

    //此处应该连减法都做了
    //_Pyr_Up_col_Batch_GPU << <oGrid, oThread >> > (Source_Header_GPU, pMid_Header_GPU, Dest_Header_GPU);
    _Pyr_Up_col_Subtract_Batch_GPU << <oGrid, oThread >> > (Source_Header_GPU, pMid_Header_GPU, Dest_Header_GPU);
    //Disp_Cuda_Error();
    //bSave_Image_GPU("c:\\tmp\\2.bmp", Dest[0]);
    return;
}

__global__ void _Pyr_Down_row_Batch_float_GPU(Image::Part_1 Source[], Image::Part_1 Mid[], Border_Type iBorder_Type = BORDER_REFLECT101)
{//行处理，这个更快
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image::Part_1 oSource, oMid;
    if (threadIdx.x == 0)
    {
        oSource = Source[blockIdx.y];
        oMid = Mid[blockIdx.y];
    }
    __syncthreads();
    if (iThread_ID >= oSource.m_iHeight * oMid.m_iWidth)
        return;

    short x = iThread_ID % oMid.m_iWidth,
        y = iThread_ID / oMid.m_iWidth;
    short x1 = x << 1;
    float* pSource = &((float*)(oSource.m_pChannel[0]))[y * oSource.m_iWidth];
    short Source_Pos[4] = { (short)iGet_Border_x_GPU(x1 - 2, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 - 1, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 + 1, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 + 2, oSource.m_iWidth, iBorder_Type) };

    ((float*)oMid.m_pChannel[0])[iThread_ID] = pSource[x1] * 6 +
        ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) * 4) +
        pSource[Source_Pos[0]] + pSource[Source_Pos[3]];

    //if (blockIdx.y == 1 && y == 1 && x == 384)
    //{
    //    //printf("%f\n", ((float*)oMid.m_pChannel[0])[iThread_ID]);
    //    printf("%f %f %f %f %f\n", pSource[Source_Pos[0]], pSource[Source_Pos[1]],
    //        pSource[x1],
    //        pSource[Source_Pos[2]], pSource[Source_Pos[3]]);
    //}
}

__global__ void _Pyr_Down_row_float_GPU(Image::Part_1 oSource, int iMid_Width, float* pMid, Border_Type iBorder_Type = BORDER_REFLECT101)
{//行处理，这个更快
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oSource.m_iHeight * iMid_Width)
        return;

    short x = iThread_ID % iMid_Width,
        y = iThread_ID / iMid_Width;
    short x1 = x << 1;
    float* pSource = &((float*)(oSource.m_pChannel[0]))[y * oSource.m_iWidth];
    short Source_Pos[4] = { (short)iGet_Border_x_GPU(x1 - 2, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 - 1, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 + 1, oSource.m_iWidth, iBorder_Type),
        (short)iGet_Border_x_GPU(x1 + 2, oSource.m_iWidth, iBorder_Type) };

    pMid[iThread_ID] = pSource[x1] * 6 +
        ((pSource[Source_Pos[1]] + pSource[Source_Pos[2]]) * 4) +
        pSource[Source_Pos[0]] + pSource[Source_Pos[3]];

}

__global__ void _Pyr_Down_col_Batch_float_GPU(Image::Part_1 Mid[], Image::Part_1 Dest[], Border_Type iBorder_Type = BORDER_REFLECT101)
{//列方向
    int iThread_ID = GET_THREAD_ID();
    __shared__ Image::Part_1 oMid, oDest;
    if (threadIdx.x == 0)
    {
        oMid = Mid[blockIdx.y];
        oDest = Dest[blockIdx.y];
    }
    __syncthreads();

    if (iThread_ID >= oDest.m_iHeight * oDest.m_iWidth)
        return;

    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;

    short y1 = y << 1;
    int iMid_Pos = y1 * oDest.m_iWidth;
    //int iMid_Size = iMid_Height * oDest.m_iWidth;
    float* pSource = &((float*)oMid.m_pChannel[0])[x];

    int Mid_Pos[4] = { iGet_Border_y_GPU(y1 - 2, oMid.m_iHeight, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 - 1, oMid.m_iHeight, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 + 1, oMid.m_iHeight, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 + 2, oMid.m_iHeight, iBorder_Type) * oDest.m_iWidth };

    float fValue = pSource[iMid_Pos] * 6 +
        ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) * 4) +
        pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];
    ((float*)(oDest.m_pChannel[0]))[iThread_ID] = fValue * (1.f / 256.f);

    /*if (blockIdx.y == 1 && y == 0 && x == 271)
        printf("%f %f %f %f %f\n", pSource[Mid_Pos[0]], pSource[Mid_Pos[1]],
            pSource[iMid_Pos],
            pSource[Mid_Pos[2]], pSource[Mid_Pos[3]]);*/
            //printf("%f\n", ((float*)(oDest.m_pChannel[0]))[iThread_ID]);

    return;
}

__global__ void _Pyr_Down_col_float_GPU(float* pMid, int iMid_Height, Image::Part_1 oDest, Border_Type iBorder_Type = BORDER_REFLECT101)
{//列方向
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oDest.m_iHeight * oDest.m_iWidth)
        return;
    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;

    short y1 = y << 1;
    int iMid_Pos = y1 * oDest.m_iWidth;
    //int iMid_Size = iMid_Height * oDest.m_iWidth;
    float* pSource = &pMid[x];
    int Mid_Pos[4] = { iGet_Border_y_GPU(y1 - 2, iMid_Height, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 - 1, iMid_Height, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 + 1, iMid_Height, iBorder_Type) * oDest.m_iWidth,
        iGet_Border_y_GPU(y1 + 2, iMid_Height, iBorder_Type) * oDest.m_iWidth };

    float fValue = pSource[iMid_Pos] * 6 +
        ((pSource[Mid_Pos[1]] + pSource[Mid_Pos[2]]) * 4) +
        pSource[Mid_Pos[0]] + pSource[Mid_Pos[3]];

    ((float*)(oDest.m_pChannel[0]))[iThread_ID] = fValue * (1.f / 256.f);
}

void Pyr_Down_float_GPU(Image::Part_1 oSource, Image::Part_1 oDest, float* pAux)
{//单图拉普拉斯下采样
    Data_Block<float*, 1>oMid;
    int iSize = oSource.m_iHeight * oDest.m_iWidth;

    if (pAux)
        oMid.Data[0] = pAux;
    else
        oMid.Data[0] = (float*)pMalloc_GPU(iSize * sizeof(float));

    dim3 oThread, oGrid;
    //先搞行方向
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Pyr_Down_row_float_GPU << <oGrid, oThread >> > (oSource, oDest.m_iWidth, oMid.Data[0]);
    //Disp_Part_GPU(oMid.Data[0], oDest.m_iWidth, 342, 12, 10, 1);

    //再到列方向
    iSize = oDest.m_iWidth * oDest.m_iHeight;
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Pyr_Down_col_float_GPU << <oGrid, oThread >> > (oMid.Data[0], oSource.m_iHeight, oDest);

    //Disp_Part_GPU(oMid.Data[0], oDest.m_iWidth, 342, 12, 10, 1);
    //Disp_Part_GPU((float*)oDest.m_pChannel[0], oDest.m_iWidth, 342, 5, 10, 1);
    //Disp_Cuda_Error();
    //Temp_Compare("c:\\tmp\\Dest_1.bin", (float*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight);

    //Disp_Part_GPU((float*)oDest.m_pChannel[0], oDest.m_iWidth, 0, 100, oDest.m_iWidth, 1);
    return;
}
__global__ void Disp_Image_Header(Image::Part_1 Image_Header_GPU[], int iImage_Count)
{
    for (int i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oImage = Image_Header_GPU[i];
        printf("%d %d\n", oImage.m_iWidth, oImage.m_iHeight);
        for (int j = 0; j < 100; j++)
            printf("%d\n", oImage.m_pChannel[0][j]);
    }

}
void Pyr_Down_Batch_1_Leve_float(Image::Part_1 Source[], Image::Part_1 Dest[], Image::Part_1 Mid_Header[],
    Image::Part_1 Source_Header_GPU[], Image::Part_1 Dest_Header_GPU[], Image::Part_1 Mid_Header_GPU[],
    int iImage_Count, int iLevel = -1)
{
    /*Disp_Image_Header<<<1,1>>>(Source_Header_GPU,iImage_Count);
    Disp_Cuda_Error();*/

    static int iCount = 0;
    {
        int iMax_Size = 0;
        for (int j = 0; j < iImage_Count; j++)
        {
            Image::Part_1 oSource = Source[j],
                oMid = Mid_Header[j];
            oMid.m_iHeight = oSource.m_iHeight;
            oMid.m_iWidth = (oSource.m_iWidth + 1) / 2;
            if (oMid.m_iHeight * oMid.m_iWidth > iMax_Size)
                iMax_Size = oMid.m_iHeight * oMid.m_iWidth;
            Mid_Header[j] = oMid;
        }
        //cudaMemcpy(Mid_Header_GPU, Mid_Header, iImage_Count * sizeof(Image::Part_1), cudaMemcpyHostToDevice);

        dim3 oThread, oGrid;
        oThread.x = Min(iMax_Size, 512);
        oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
        oGrid.y = iImage_Count;
        _Pyr_Down_row_Batch_float_GPU << <oGrid, oThread >> > (Source_Header_GPU, Mid_Header);

        /*Disp_Cuda_Error();
        for (int j = 0; j < iImage_Count; j++)
        {
            Image::Part_1 oMid = Mid_Header[j];
            char File[256];
            sprintf(File, "c:\\tmp\\Weight_Mid_Level_%d_%d.bin", iLevel, j);
            Temp_Compare(File, (float*)oMid.m_pChannel[0], oMid.m_iWidth, oMid.m_iHeight);
        }*/

        //Disp_Part_GPU((float*)pMid_Header[0].m_pChannel[0], pMid_Header[0].m_iWidth, 0, 100, pMid_Header[0].m_iWidth, 1);
        //Disp_Cuda_Error();
        //Temp_Compare("c:\\tmp\\1.bin", (float*)pMid_Header[0].m_pChannel[0], pMid_Header[0].m_iWidth, Src_Pyr_Laplace[0][0].m_iHeight);
        iMax_Size = 0;
        for (int j = 0; j < iImage_Count; j++)
        {
            Image::Part_1 oDest = Dest[j];
            if (oDest.m_iHeight * oDest.m_iWidth > iMax_Size)
                iMax_Size = oDest.m_iHeight * oDest.m_iWidth;
        }
        oThread.x = Min(iMax_Size, 512);
        oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
        oGrid.y = iImage_Count;
        _Pyr_Down_col_Batch_float_GPU << <oGrid, oThread >> > (Mid_Header, Dest_Header_GPU);

        /*Disp_Cuda_Error();
        for (int j = 0; j < iImage_Count; j++)
        {
            Image::Part_1 oDest = Dest[j];
            char File[256];
            sprintf(File, "c:\\tmp\\Weight_Level_%d_%d.bin", iLevel, j);
            Temp_Compare(File, (float*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight);
        }*/
    }

    //Disp_Cuda_Error();
    iCount++;
    return;
}

void Pyr_Down_Batch_float(Image::Part_1* Src_Pyr_Laplace[6], Image::Part_1* Src_Pyr_Laplace_Header_GPU[6], int iImage_Count)
{
    Image::Part_1* pMid_Header, * pMid_Header_GPU;
    pMid_Header = (Image::Part_1*)pMalloc(iImage_Count * sizeof(Image::Part_1));
    pMid_Header_GPU = (Image::Part_1*)pMalloc_GPU(iImage_Count * sizeof(Image::Part_1));

    //计算Mid 用的空间
    int iMid_Size = 0;
    for (int i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oSource = Src_Pyr_Laplace[0][i],
            oMid;
        oMid.m_iHeight = oSource.m_iHeight;
        oMid.m_iWidth = (oSource.m_iWidth + 1) / 2;
        int iMid_1 = oMid.m_iWidth * oSource.m_iHeight;
        //中间Mid所需空间，由于先行后列处理
        //故此Mid所需空间时 源宽一半 * 源高
        //iMid_Size += (oSource.m_iWidth + 1) / 2 * oSource.m_iHeight;
        iMid_Size += iMid_1;
        //oMid.m_pChannel[0] = (unsigned char*)((unsigned long long)iMid_Size*4); //借用一下，玩个花
        pMid_Header[i] = oMid;
    }
    float* p, * pMid = (float*)pMalloc_GPU(iMid_Size * sizeof(float));
    p = pMid;
    for (int i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oMid = pMid_Header[i];
        oMid.m_pChannel[0] = (unsigned char*)p;
        p += oMid.m_iHeight * oMid.m_iWidth;
        pMid_Header[i] = oMid;
    }
    //cudaMemcpy(pMid_Header_GPU, pMid_Header, iImage_Count * sizeof(Image::Part_1), cudaMemcpyHostToDevice);

    int iDown_Sample_Thresold = 120000;
    for (int i = 0; i < 5; i++)
    {
        //Disp_Cuda_Error();
        //unsigned long long tStart = iGet_Tick_Count();
        //for (int k = 0; k < 10000; k++)
        {
            Image::Part_1 oImage = Src_Pyr_Laplace[i + 1][0];
            if (oImage.m_iHeight * oImage.m_iWidth < iDown_Sample_Thresold)
            {
                Pyr_Down_Batch_1_Leve_float(Src_Pyr_Laplace[i], Src_Pyr_Laplace[i + 1], pMid_Header,
                    Src_Pyr_Laplace_Header_GPU[i], Src_Pyr_Laplace_Header_GPU[i + 1], pMid_Header_GPU,
                    iImage_Count, i + 1);
                //此处有个问题，一定要加同步，怀疑是用了pMid_Header非pMid_Header_GPU的缘故
                Disp_Cuda_Error();
            }
            else
            {
                for (int j = 0; j < iImage_Count; j++)
                    Pyr_Down_float_GPU(Src_Pyr_Laplace[i][j], Src_Pyr_Laplace[i + 1][j], pMid);
            }
        }
        //Disp_Cuda_Error();
        //printf("%lld\n", iGet_Tick_Count() - tStart);
    }

    //exit(0);

    //for (int i = 0; i < 5; i++)
    //{
    //    for (int j = 0; j < iImage_Count; j++)
    //    {
    //        char File[256];
    //        sprintf(File, "c:\\tmp\\Weight_Level_%d_%d.bin", i + 1, j);
    //        Image::Part_1 oDest = Src_Pyr_Laplace[i+1][j];
    //        //printf("Level:%d Image:%d\n", i + 1, j);
    //        Temp_Compare(File, (float*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight);
    //    }
    //}

    //Image::Part_1 oDest = Src_Pyr_Laplace[1][0];
    //Temp_Compare("c:\\tmp\\Dest_1.bin", (float*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight);

    if (pMid_Header)
        Free(pMid_Header);
    if (pMid_Header_GPU)
        Free_GPU(pMid_Header_GPU);
    if (pMid)
        Free_GPU(pMid);
    return;
}
void Pyr_Down_1_float(Image::Part_1* Src_Pyr_Laplace[6], Image::Part_1* Src_Pyr_Laplace_Header_GPU[6], int iImage_Count)
{
    float** pMid_Header, ** pMid_Header_GPU;
    pMid_Header = (float**)pMalloc(iImage_Count * sizeof(float*));
    pMid_Header_GPU = (float**)pMalloc_GPU(iImage_Count * sizeof(float*));

    //计算Mid 用的空间
    int iMid_Size = 0;
    for (int i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oSource = Src_Pyr_Laplace[0][i],
            oDest = Src_Pyr_Laplace[1][i];;
        int iMid_1 = (oSource.m_iWidth + 1) / 2 * oSource.m_iHeight,
            iMid_2 = oDest.m_iHeight * oSource.m_iWidth;

        if (iMid_2 > iMid_1)
            printf("Odd value in Create_Laplace_Pyramid\n");
        //中间Mid所需空间，由于先行后列处理
        //故此Mid所需空间时 源宽一半 * 源高
        //iMid_Size += (oSource.m_iWidth + 1) / 2 * oSource.m_iHeight;
        iMid_Size += Max(iMid_1, iMid_2);
        pMid_Header[i] = (float*)((unsigned long long)iMid_Size); //借用一下，玩个花
    }
    float* p, * pMid = (float*)pMalloc_GPU(iMid_Size * sizeof(float));

    //然后构造pMid_Header 数组
    p = pMid;
    for (int i = 0; i < iImage_Count; i++)
    {
        int iMid_Size = (int)((unsigned long long)(pMid_Header[i]));
        pMid_Header[i] = p, p += iMid_Size;
        //printf("here");
    }
    cudaMemcpy(pMid_Header_GPU, pMid_Header, iImage_Count * sizeof(float*), cudaMemcpyHostToDevice);
    for (int i = 0; i < 5; i++)
    {// //首先下采样 总耗时 4900ms
        for (int j = 0; j < iImage_Count; j++)
        {
            //Disp_Cuda_Error();
            //unsigned long long tStart = iGet_Tick_Count();
            //for(int k=0;k<10000;k++)
            Pyr_Down_float_GPU(Src_Pyr_Laplace[i][j], Src_Pyr_Laplace[i + 1][j], pMid);
            //Disp_Cuda_Error();
            //printf("%lld\n", iGet_Tick_Count() - tStart);
            //exit(0);
        }
    }

    //Disp_Part_GPU((float*)Src_Pyr_Laplace[0][1].m_pChannel[0], Src_Pyr_Laplace[0][1].m_iWidth, 
    //    0, 100, Src_Pyr_Laplace[0][1].m_iWidth, 1);
    //Temp_Compare("c:\\tmp\\Dest_1.bin", (float*)(Src_Pyr_Laplace[1][1].m_pChannel[0]), Src_Pyr_Laplace[1][1].m_iWidth, Src_Pyr_Laplace[1][1].m_iHeight);

    ////Pyr_Down_float_GPU(Src_Pyr_Laplace[0][1], Src_Pyr_Laplace[1][1], pMid);
    //for (int i = 1; i < 6; i++)
    //{
    //    char File[256];
    //    sprintf(File, "c:\\tmp\\Dest_%d.bin", i);
    //    int j = 3;
    //    Temp_Compare(File, (float*)(Src_Pyr_Laplace[i][j].m_pChannel[0]), Src_Pyr_Laplace[i][j].m_iWidth, Src_Pyr_Laplace[i][j].m_iHeight);
    //}
    if (pMid_Header)
        Free(pMid_Header);
    if (pMid_Header_GPU)
        Free_GPU(pMid_Header_GPU);
    if (pMid)
        Free_GPU(pMid);
    return;
}
void Create_Laplace_Pyramid(Image* Src_Pyr_Laplace[6], Image* Src_Pyr_Laplace_Header_GPU[6], int iImage_Count)
{
    unsigned short* (*pMid_Header)[3], * (*pMid_Header_GPU)[3];
    pMid_Header = (unsigned short* (*)[3])pMalloc(iImage_Count * 3 * sizeof(unsigned short*));
    //pKernel_GPU = (float*)pMalloc_GPU((r * 2 + 1)*sizeof(float));
    pMid_Header_GPU = (unsigned short* (*)[3])pMalloc_GPU(iImage_Count * 3 * sizeof(unsigned short*));

    //计算Mid 用的空间
    int iMid_Size = 0;
    for (int i = 0; i < iImage_Count; i++)
    {
        Image::Part_1 oSource = Src_Pyr_Laplace[0][i].m_oPart_1,
            oDest = Src_Pyr_Laplace[1][i].m_oPart_1;
        int iMid_1 = (oSource.m_iWidth + 1) / 2 * oSource.m_iHeight,
            iMid_2 = oDest.m_iHeight * oSource.m_iWidth;

        if (iMid_2 > iMid_1)
            printf("Odd value in Create_Laplace_Pyramid\n");
        //中间Mid所需空间，由于先行后列处理
        //故此Mid所需空间时 源宽一半 * 源高
        //iMid_Size += (oSource.m_iWidth + 1) / 2 * oSource.m_iHeight;
        iMid_Size += Max(iMid_1, iMid_2);
        pMid_Header[i][0] = (unsigned short*)((unsigned long long)(iMid_Size)); //借用一下，玩个花
    }
    unsigned short* p, * pMid = (unsigned short*)pMalloc_GPU(iMid_Size * 3 * sizeof(unsigned short));

    //然后构造pMid_Header 数组
    p = pMid;
    for (int i = 0; i < iImage_Count; i++)
    {
        int iMid_Size = (int)((unsigned long long)(pMid_Header[i][0]));
        for (int j = 0; j < Src_Pyr_Laplace[0][i].m_iChannel_Count; j++)
            pMid_Header[i][j] = p, p += iMid_Size;
        //printf("here");
    }
    cudaMemcpy(pMid_Header_GPU, pMid_Header, iImage_Count * 3 * sizeof(unsigned short*), cudaMemcpyHostToDevice);

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();

    //小于这个阀值用批处理，多于这个阀值用串行
    //这个值只在4060上成立，其他GPU还得测

    int iDown_Sample_Thresold = 250000;
    //先做下采样，从第一层向下搞
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int k=0;k<10000;k++)
    for (int i = 0; i < 5; i++)
    {// //首先下采样 总耗时 4900ms
        //Disp_Cuda_Error();
        //unsigned long long tStart = iGet_Tick_Count();
        //for(int k=0;k<10000;k++)
        {
            Image::Part_1 oImage = Src_Pyr_Laplace[i + 1][0].m_oPart_1;
            if (oImage.m_iHeight * oImage.m_iWidth < iDown_Sample_Thresold)
                Pyr_Down_Batch_GPU(Src_Pyr_Laplace[i], Src_Pyr_Laplace_Header_GPU[i], iImage_Count,
                    pMid_Header, pMid_Header_GPU, Src_Pyr_Laplace[i + 1], Src_Pyr_Laplace_Header_GPU[i + 1]);
            else
            {
                for (int j = 0; j < iImage_Count; j++)
                    Pyr_Down_GPU(Src_Pyr_Laplace[i][j], Src_Pyr_Laplace[i + 1][j], pMid);
            }

            /*char File[256];
            for (int j = 0; j < iImage_Count; j++)
            {
                sprintf(File, "c:\\tmp\\Dest_%d.bmp", j);
                bSave_Image_GPU(File, &Src_Pyr_Laplace_Header_GPU[i + 1][j]);
            }*/

            //for (int j = 0; j < iImage_Count; j++)
            //{
            //    sprintf(File, "c:\\tmp\\Source_%d.bmp", j);
            //    bSave_Image_GPU(File, &Src_Pyr_Laplace_Header_GPU[i + 1][j]);
            //}
            //Compare_Image("c:\\tmp\\Source_0.bmp", "c:\\tmp\\Dest_0.bmp");
            //Compare_Image("c:\\tmp\\Source_1.bmp", "c:\\tmp\\Dest_1.bmp");
            //Compare_Image("c:\\tmp\\Source_2.bmp", "c:\\tmp\\Dest_2.bmp");
            //Compare_Image("c:\\tmp\\Source_3.bmp", "c:\\tmp\\Dest_3.bmp");
            //////printf("Here");
        }
        //Disp_Cuda_Error();
        //printf("%lld\n", iGet_Tick_Count() - tStart);
    }

    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //再到上采样
    iDown_Sample_Thresold = 120000;     //上采样策略阀值
    for (int i = 0; i < 5; i++)
    {
        //Disp_Cuda_Error();
        //unsigned long long tStart = iGet_Tick_Count();
        //for (int k = 0; k < 10000; k++)
        {
            Image::Part_1 oImage = Src_Pyr_Laplace[i + 1][0].m_oPart_1;
            if (oImage.m_iHeight * oImage.m_iWidth < iDown_Sample_Thresold)
                Pyr_Up_Batch_GPU(Src_Pyr_Laplace[i + 1], Src_Pyr_Laplace_Header_GPU[i + 1], iImage_Count,
                    pMid_Header, pMid_Header_GPU, Src_Pyr_Laplace[i], Src_Pyr_Laplace_Header_GPU[i]);
            else
            {
                for (int j = 0; j < iImage_Count; j++)
                    Pyr_Up_GPU(Src_Pyr_Laplace[i + 1][j], Src_Pyr_Laplace[i][j]);
            }
        }
        //Disp_GPU((char*)Src_Pyr_Laplace[0][0].m_pChannel[0], 1, 2);

        //Disp_Cuda_Error();
        //printf("%lld\n", iGet_Tick_Count() - tStart);

        /*for (int j = 0; j < iImage_Count; j++)
        {
            char File_2[256], File_1[256];
            sprintf(File_2, "c:\\tmp\\Dest_%d.bmp", j);
            bSave_Image_GPU(File_2, Src_Pyr_Laplace[j]);

            Set_Color_GPU(Src_Pyr_Laplace[i][j]);
            Pyr_Up_GPU(Src_Pyr_Laplace[i + 1][j], Src_Pyr_Laplace[i][j]);
            sprintf(File_1, "c:\\tmp\\Source_%d.bmp", j);
            bSave_Image_GPU(File_1, Src_Pyr_Laplace[i][j]);
            Compare_Image(File_1, File_2);
        }*/
    }

    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);
    //bSave_Image_GPU("c:\\tmp\\1.bmp", Src_Pyr_Laplace[5][0]);


    if (pMid)
        Free_GPU(pMid);
    if (pMid_Header)
        Free(pMid_Header);
    if (pMid_Header_GPU)
        Free_GPU(pMid_Header_GPU);
    return;
}
__global__ void Set_Pyramid_1(Image* Buffer, Image* Pyramid[6], int iImage_Count)
{
    Pyramid[threadIdx.x] = &Buffer[threadIdx.x * iImage_Count];
    return;
}

__global__ void _Copy_Make_Border_float_GPU(Image oSource, Image::Part_1 oDest, short iLeft, short iTop/*, short iRight, short iBottom*/)
{//源是unsigned char, 一般就是Mask, 目标是 float，借用 oDest.m_pChannel[0]
    //由于特殊方法，故此没有通用价值，不放Image中, 扩边方法用 BORDER_CONSTANT
#define FACTOR  1.f/255.f   //注意，此处日过不用float, 就慢6倍

    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oSource.m_iWidth * oSource.m_iHeight)
        return;

    int y = iThread_ID / oSource.m_iWidth,
        x = iThread_ID % oSource.m_iWidth;
    //float fValue = oSource.m_pChannel[3][iThread_ID];
    //((float*)(oDest.m_pChannel[0]))[(y + iTop) * oDest.m_iWidth + x + iLeft] = fValue * FACTOR;

    ((float*)(oDest.m_pChannel[0]))[(y + iTop) * oDest.m_iWidth + x + iLeft] = (float)oSource.m_pChannel[3][iThread_ID] * (float)FACTOR;
    return;
#undef FACTOR
}


void Copy_Make_Border_float_GPU(Image oSource, Image::Part_1 oDest, short iLeft, short iTop, short iRight, short iBottom)
{//毫无商用价值，直接用BORDER_CONSTANT

    int iSize = oSource.m_iHeight * oSource.m_iWidth;
    dim3 oThread, oGrid;
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x) / oThread.x;
    oGrid.y = 1;

    //先设置0等于补黑边
    //cudaMemset(oDest.m_pChannel[0], 0, oDest.m_iWidth * oDest.m_iHeight * sizeof(float));
    //简单方法，后面再收拾它
    _Copy_Make_Border_float_GPU << <oGrid, oThread >> > (oSource, oDest, iLeft, iTop);

    //Disp_Cuda_Error();
    //Disp_Part_GPU(oSource.m_pChannel[3], oSource.m_iWidth, 686, 1, 10, 1);
    //Disp_Part_GPU((float*)oDest.m_pChannel[0], oDest.m_iWidth, 0, 100, oDest.m_iWidth, 1);
    //Disp_Cuda_Error();
    //Temp_Compare("c:\\tmp\\1.bin", (float*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight);

    return;
}
__global__ void _Feed_Blender_1_GPU(Image oSource, Image::Part_1 oSource_Weight,
    Image::Part_1 oDest, Image::Part_1 oDest_Weight, short iLeft, short iTop, char bSigned)
{//以oWeight的大小为准
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oSource_Weight.m_iHeight * oSource_Weight.m_iWidth)
        return;
    unsigned short y = iThread_ID / oSource_Weight.m_iWidth,
        x = iThread_ID % oSource_Weight.m_iWidth;
    float fWeight = ((float*)oSource_Weight.m_pChannel[0])[y * oSource_Weight.m_iWidth + x];
    int iSource_Pos = y * oSource.m_iWidth + x;
    //int iWeight_Pos = y * oSource_Weight.m_iWidth + x;
    int iDest_Pos = (y + iTop) * oDest.m_iWidth + x + iLeft;

    if (bSigned)
    {
        ((short*)oDest.m_pChannel[0])[iDest_Pos] += (short)(((char*)oSource.m_pChannel[0])[iSource_Pos] * fWeight);
        ((short*)oDest.m_pChannel[1])[iDest_Pos] += (short)(((char*)oSource.m_pChannel[1])[iSource_Pos] * fWeight);
        ((short*)oDest.m_pChannel[2])[iDest_Pos] += (short)(((char*)oSource.m_pChannel[2])[iSource_Pos] * fWeight);
    }
    else
    {
        ((short*)oDest.m_pChannel[0])[iDest_Pos] += (short)(((unsigned char*)oSource.m_pChannel[0])[iSource_Pos] * fWeight);
        ((short*)oDest.m_pChannel[1])[iDest_Pos] += (short)(((unsigned char*)oSource.m_pChannel[1])[iSource_Pos] * fWeight);
        ((short*)oDest.m_pChannel[2])[iDest_Pos] += (short)(((unsigned char*)oSource.m_pChannel[2])[iSource_Pos] * fWeight);
    }
    ((float*)oDest_Weight.m_pChannel[0])[iDest_Pos] += fWeight;

    return;
}
void Feed_Blender_1_GPU(Image oSource, Image::Part_1 oSource_Weight,
    Image::Part_1 oDest, Image::Part_1 oDest_Weight, short iLeft, short iTop, char bSigned = 1)
{
    static int iCount = 0;
    dim3 oGrid, oThread;
    int iSize = oSource_Weight.m_iWidth * oSource_Weight.m_iHeight;
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Feed_Blender_1_GPU << <oGrid.x, oThread.x >> > (oSource, oSource_Weight, oDest, oDest_Weight, iLeft, iTop, bSigned);
    //Disp_Cuda_Error();
    //if (iCount == 5)
    //{
    //    //Disp_Part_GPU((unsigned char *)oSource.m_pChannel[0], oSource.m_iWidth, 18, 0, 1, 1);
    //    //Disp_Part_GPU((short*)oDest.m_pChannel[0], oDest.m_iWidth, 59, 5, 1, 1);
    //    //Temp_Compare<short>("c:\\tmp\\1.bin", (short*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight); 
    //    //Disp_Part_GPU((float*)oSource_Weight.m_pChannel[0], oSource_Weight.m_iWidth, 59, 5, 1, 1);
    //    //Disp_Part_GPU((float*)oDest_Weight.m_pChannel[0], oDest_Weight.m_iWidth, 882, 15, 1, 1);
    //    //Temp_Compare<float>("c:\\tmp\\0_float.bin", (float*)oDest_Weight.m_pChannel[0], oDest_Weight.m_iWidth, oDest_Weight.m_iHeight);
    //}
    //Disp_Cuda_Error();
    iCount++;
    return;
}

void _Feed_Blender(Image* Src_Pyr_Laplace[6], Image::Part_1* pPyr_Weight[6],
    Image::Part_1 dst_pyr_laplace[6], Image::Part_1 dst_band_weights[6],
    Image* Src_Pyr_Laplace_Header_GPU[6], Image::Part_1* Pyr_Weight_Header_GPU[6],
    Image::Part_1 dst_pyr_laplace_Header_GPU[6], Image::Part_1 dst_band_weights_Header_GPU[6],
    short (*pLTRB)[2][2], int dst_roi[2][2], int iImage_Count)
{
    const int MAX_IMAGE_COUNT = 8;
    if (iImage_Count > MAX_IMAGE_COUNT)
    {
        printf("exceeds MAX_IMAGE_COUNT in _Feed_Blender\n");
        exit(0);
    }
    short LTBR[8][2][2], roi[8][2][2];  // , (*roi_Header_GPU)[2][2];

    for (int j = 0; j < iImage_Count; j++)
    {
        LTBR[j][0][0] = pLTRB[j][0][0] - dst_roi[0][0];
        LTBR[j][0][1] = pLTRB[j][0][1] - dst_roi[0][1];
        LTBR[j][1][0] = pLTRB[j][1][0] - dst_roi[0][0];
        LTBR[j][1][1] = pLTRB[j][1][1] - dst_roi[0][1];
    }

    //诸葛执行
    for (int i = 0; i <= 5; i++)
    {
        for (int j = 0; j < iImage_Count; j++)
        {
            roi[j][0][0] = LTBR[j][0][0];
            roi[j][0][1] = LTBR[j][0][1];
            roi[j][1][0] = LTBR[j][1][0] - LTBR[j][0][0];
            roi[j][1][1] = LTBR[j][1][1] - LTBR[j][0][1];

            /* Disp_Cuda_Error();
             unsigned long long tStart = iGet_Tick_Count();
             for(int k=0;k<10000;k++)*/
            Feed_Blender_1_GPU(Src_Pyr_Laplace[i][j], pPyr_Weight[i][j],
                dst_pyr_laplace[i], dst_band_weights[i], roi[j][0][0], roi[j][0][1], i < 5);
            /* Disp_Cuda_Error();
             printf("%lld\n", iGet_Tick_Count() - tStart);
             exit(0);*/
             //Temp_Compare<short>("c:\\tmp\\1.bin", (short*)dst_pyr_laplace[0].m_pChannel[0], dst_pyr_laplace[0].m_iWidth, dst_pyr_laplace[0].m_iHeight);
            LTBR[j][0][0] >>= 1, LTBR[j][0][1] >>= 1, LTBR[j][1][0] >>= 1, LTBR[j][1][1] >>= 1;
        }
    }

    ////批处理，由于存在写重叠，此处是个难点，问题事无法对字节进行源自操作
    //for (int i = 0; i <= 5; i++)
    //{
    //    for (int j = 0; j < iImage_Count; j++)
    //    {
    //        roi[j][0][0] = LTBR[j][0][0];
    //        roi[j][0][1] = LTBR[j][0][1];
    //        roi[j][1][0] = LTBR[j][1][0] - LTBR[j][0][0];
    //        roi[j][1][1] = LTBR[j][1][1] - LTBR[j][0][1];
    //    }
    //}
    //roi_Header_GPU = (short(*)[2][2])pMalloc_GPU(iImage_Count * 2 * 2 * sizeof(short));
    //cudaMemcpy(roi_Header_GPU, roi, iImage_Count * 2 * 2 * sizeof(short), cudaMemcpyHostToDevice);

   /* char File[256];
    for (int i = 0; i < 6; i++)
    {
        printf("Level:%d\n", i);
        for (int j = 0; j < 3; j++)
        {
            sprintf(File, "c:\\tmp\\Level_%d_Channel_%d.bin", i, j);
            Temp_Compare<short>(File, (short*)dst_pyr_laplace[i].m_pChannel[j], dst_pyr_laplace[i].m_iWidth, dst_pyr_laplace[i].m_iHeight,0);
        }
        sprintf(File, "c:\\tmp\\Level_%d_Weight.bin", i);
        Temp_Compare<float>(File, (float*)dst_band_weights[i].m_pChannel[0], dst_band_weights[i].m_iWidth, dst_band_weights[i].m_iHeight, 0);
    }*/

    return;
}
template<typename _T>void Feed_Blender(Image Warp[], Image Warp_Header_GPU[], Stitch<_T>* poStitch, Blender* poBlender)
{
    int iMax_Size,                                  //iMax_Size是最大加了Border图像面积
        iMax_Source_Width, iMax_Source_Height,      //最大加了源图像的长宽
        iSize, iImage_Count = poStitch->m_iImage_Count;
    int (*pLTRB)[2][2], //用于装Left, Top, Right, Bottom 4个Margin
        (*pLTRB_GPU)[2][2], (*pSize)[2];
    short (*pLTRB_New)[2][2]; //tl_new, br_new

    //Image** Src_Pyr_Laplace;
    Image* Src_Pyr_Laplace[6], * Src_Pyr_Laplace_Header_GPU[6];
    Image* pImage_Buffer, * pImage_Buffer_GPU;

    //先分配内存
    Light_Ptr oPtr, oPtr_GPU;
    unsigned char* p;
    iSize = iImage_Count * 2 * sizeof(int) +                //Size
        iImage_Count * 2 * 2 * sizeof(int) +                //pLTRB
        iImage_Count * 2 * 2 * sizeof(short) +                //pLTRB_New
        iImage_Count * 6 * sizeof(Image) + 128 * 2;              //Src_Pyr_Laplace

    iSize = ALIGN_SIZE_128(iSize);
    Attach_Light_Ptr(oPtr, (unsigned char*)pMalloc(iSize), iSize, 0);

    iSize = iImage_Count * 2 * 2 * sizeof(int) +
        iImage_Count * 2 * sizeof(int) +
        iImage_Count * 2 * 2 * sizeof(int);
    Malloc(oPtr, iSize, p);
    pLTRB = (int(*)[2][2])p;
    pSize = (int(*)[2])(pLTRB + iImage_Count);
    pLTRB_New = (short(*)[2][2])(pSize + iImage_Count);

    iSize = iImage_Count * 6 * sizeof(Image);
    Malloc(oPtr, iSize, p);
    Src_Pyr_Laplace[0] = pImage_Buffer = (Image*)p;
    for (int i = 0; i < 6; i++)
        Src_Pyr_Laplace[i] = &pImage_Buffer[i * iImage_Count];

    Feed_Blender_Get_Pos(Warp, poStitch, poBlender, pLTRB, pSize, pLTRB_New);

    iMax_Size = iSize = iMax_Source_Width = iMax_Source_Height = 0;  //iMax_Size是最大加了Border图像面积
    for (int i = 0; i < iImage_Count; i++)
    {
        int iSize_1 = pSize[i][0] * pSize[i][1];
        iSize += iSize_1 * 3;             //src_pyr_laplace
        if (iSize_1 > iMax_Size)
            iMax_Size = iSize_1;
        if (Warp[i].m_iWidth > iMax_Source_Width)
            iMax_Source_Width = Warp[i].m_iWidth;
        if (Warp[i].m_iHeight > iMax_Source_Height)
            iMax_Source_Height = Warp[i].m_iHeight;

        //iSize += pSize[i][0] * pSize[i][1] * sizeof(float); //weight_pyr_gaus
        int Size_1[2] = { (pSize[i][0] + 1) / 2,(pSize[i][1] + 1) / 2 };
        for (int j = 0; j < 6; j++)
        {
            iSize += Size_1[0] * Size_1[1] * 3;             //src_pyr_laplace
            //iSize += Size_1[0] * Size_1[1] * sizeof(float); //weight_pyr_gaus
            Size_1[0] = (Size_1[0] + 1) / 2;
            Size_1[1] = (Size_1[1] + 1) / 2;
        }
    }

    iSize += iImage_Count * 2 * 2 * sizeof(int) +   //pSize_GPU,
        iImage_Count * 6 * sizeof(Image);       //Src_Pyr_Laplace_Header_GPU
    Attach_Light_Ptr(oPtr_GPU, (unsigned char*)pMalloc_GPU(iSize), iSize, 0);
    iSize = 0;

    for (int i = 0; i < iImage_Count; i++)
        Init_Image_GPU(&Src_Pyr_Laplace[0][i], pSize[i][0], pSize[i][1], Image::IMAGE_TYPE_BMP, 24, &oPtr_GPU);

    for (int i = 1; i < 6; i++)
    {
        for (int j = 0; j < iImage_Count; j++)
        {
            Image::Part_1 oPrev = Src_Pyr_Laplace[i - 1][j].m_oPart_1;
            Init_Image_GPU(&Src_Pyr_Laplace[i][j], (oPrev.m_iWidth + 1) / 2,
                (oPrev.m_iHeight + 1) / 2, Image::IMAGE_TYPE_BMP, 24, &oPtr_GPU);
        }
        //printf("%d %d\n", Src_Pyr_Laplace[i][0].m_iWidth, Src_Pyr_Laplace[i][0].m_iHeight);
    }

    //将pPos_4内容抄到GPU
    iSize = iImage_Count * 2 * 2 * sizeof(int) +    //Pos_4
        iImage_Count * 6 * sizeof(Image) +          //pImage_Buffer_GPU
        6 * sizeof(Image*);                         //Src_Pyr_Laplace_Header_GPU
    Malloc(oPtr_GPU, iSize, p);
    pLTRB_GPU = (int (*)[2][2])p;
    pImage_Buffer_GPU = (Image*)(pLTRB_GPU + iImage_Count);
    //Src_Pyr_Laplace_Header_GPU = (Image**)(pImage_Buffer_GPU + iImage_Count * 6);
    cudaMemcpy(pLTRB_GPU, pLTRB, iSize, cudaMemcpyHostToDevice);
    cudaMemcpy(pImage_Buffer_GPU, pImage_Buffer, iImage_Count * 6 * sizeof(Image), cudaMemcpyHostToDevice);
    for (int i = 0; i < 6; i++)
        Src_Pyr_Laplace_Header_GPU[i] = &pImage_Buffer_GPU[i * iImage_Count];

    ////旧版扩版，有点慢
    //dim3 oThread, oGrid;
    //oThread.x = 512;
    //oGrid.x = (iMax_Size + oThread.x - 1) / oThread.x;
    //oGrid.y = iImage_Count;
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for (int i = 0; i < 10000; i++)
    ////3600 ms 很慢
    //Copy_Make_Border<<<oGrid,oThread>>>(Warp_Header_GPU, pImage_Buffer_GPU, pLTRB_GPU);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

   /* dim3 oThread, oGrid;
    oThread.x = (iMax_Source_Width + 3) / 4;
    oGrid.x = iMax_Source_Height;
    oGrid.y = iImage_Count;*/

    //bLoad_Image_GPU("c:\\tmp\\1.bmp", &Warp[3]);
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int j=0;j<10000;j++)
    {
        //3100ms, 有点慢
        //Copy_Make_Border_1 << <oGrid, oThread >> > (Warp_Header_GPU, pImage_Buffer_GPU, pLTRB_GPU);

        //注意，单独运行不需要要GPU头 
        //2900ms，咯快
        for (int i = 0; i < iImage_Count; i++)
            Copy_Make_Border_GPU(Warp[i], pImage_Buffer[i], pLTRB[i][0][0], pLTRB[i][0][1], pLTRB[i][1][0], pLTRB[i][1][1]);
    }
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);
    //bSave_Image_GPU("c:\\tmp\\3.bmp", &pImage_Buffer_GPU[3]);
    //Compare_Image("c:\\tmp\\2.bmp", "c:\\tmp\\3.bmp");

    //60-70 ms
    ////Set_Pyramid_1 << <1, 6 >> > (pImage_Buffer_GPU, Src_Pyr_Laplace_Header_GPU, iImage_Count);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\img_with_border_0.bmp", &Src_Pyr_Laplace[0][0]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\img_with_border_1.bmp", &Src_Pyr_Laplace[0][1]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\img_with_border_2.bmp", &Src_Pyr_Laplace[0][2]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\img_with_border_3.bmp", &Src_Pyr_Laplace[0][3]);

    Create_Laplace_Pyramid(Src_Pyr_Laplace, Src_Pyr_Laplace_Header_GPU, iImage_Count);
    /*for (int i = 0; i < 6; i++)
    {
        char Dest_File[256];
        sprintf(Dest_File, "c:\\tmp\\Dest_Level_%d.bmp", i);
        Image oImage = Src_Pyr_Laplace[i][0];
        bSave_Image_GPU(Dest_File, oImage);

        char Source_File[256];
        sprintf(Source_File, "c:\\tmp\\Source_Level_%d.bmp", i);
        Compare_Image(Source_File, Dest_File);
    }*/

    //轮到给搞Weight了，对应Mask
    //用一个文件头装着Weight了事
    Image::Part_1* pPyr_Weight[6], * pPyr_Weight_Header_GPU[6];
    pPyr_Weight[0] = (Image::Part_1*)pMalloc(iImage_Count * 6 * sizeof(Image::m_oPart_1));
    pPyr_Weight_Header_GPU[0] = (Image::Part_1*)pMalloc_GPU(iImage_Count * 6 * sizeof(Image::m_oPart_1));

    for (int i = 1; i < 6; i++)
        pPyr_Weight[i] = pPyr_Weight[i - 1] + 4;
    iSize = 0;
    for (int i = 0; i < iImage_Count; i++)
    {
        pLTRB[i][0][0] -= 5;
        pSize[i][0] -= 5;
        Image::Part_1 oImage;
        oImage.m_iWidth = pSize[i][0];
        oImage.m_iHeight = pSize[i][1];
        oImage.m_pChannel[0] = (unsigned char*)((unsigned long long)(iSize * 4));
        iSize += oImage.m_iHeight * oImage.m_iWidth;
        pPyr_Weight[0][i] = oImage;
    }
    for (int i = 1; i < 6; i++)
    {
        for (int j = 0; j < iImage_Count; j++)
        {
            Image::Part_1 oUpper = pPyr_Weight[i - 1][j];
            Image::Part_1 oLower;
            oLower.m_iHeight = (oUpper.m_iHeight + 1) >> 1;
            oLower.m_iWidth = (oUpper.m_iWidth + 1) >> 1;
            oLower.m_pChannel[0] = (unsigned char*)((unsigned long long)(iSize * 4));
            pPyr_Weight[i][j] = oLower;
            iSize += oLower.m_iHeight * oLower.m_iWidth;
        }
        pPyr_Weight_Header_GPU[i] = pPyr_Weight_Header_GPU[i - 1] + iImage_Count;
    }

    float* pWeight = (float*)pMalloc_GPU(iSize * sizeof(float));
    //相当于先补黑边
    cudaMemset(pWeight, 0, iSize * sizeof(float));
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < iImage_Count; j++)
            pPyr_Weight[i][j].m_pChannel[0] += (unsigned long long)pWeight;

    cudaMemcpy(pPyr_Weight_Header_GPU[0], pPyr_Weight[0], iImage_Count * 6 * sizeof(Image::Part_1), cudaMemcpyHostToDevice);

    //bSave_Comp_GPU("c:\\tmp\\1.bmp", Warp[0],3);
    /*Disp_Cuda_Error();
    unsigned long long tStart = iGet_Tick_Count();
    for(int k=0;k<10000;k++)*/
    //1800ms
    for (int i = 0; i < iImage_Count; i++)
        Copy_Make_Border_float_GPU(Warp[i], pPyr_Weight[0][i], pLTRB[i][0][0], pLTRB[i][0][1], pLTRB[i][1][0], pLTRB[i][1][1]);
    /*Disp_Cuda_Error();
    printf("%lld\n", iGet_Tick_Count() - tStart);*/

    /*Temp_Compare("c:\\tmp\\Dest_0.bin", (float*)pPyr_Weight[0][0].m_pChannel[0], pPyr_Weight[0][0].m_iWidth, pPyr_Weight[0][0].m_iHeight);
    Temp_Compare("c:\\tmp\\Dest_1.bin", (float*)pPyr_Weight[0][1].m_pChannel[0], pPyr_Weight[0][1].m_iWidth, pPyr_Weight[0][1].m_iHeight);
    Temp_Compare("c:\\tmp\\Dest_2.bin", (float*)pPyr_Weight[0][2].m_pChannel[0], pPyr_Weight[0][2].m_iWidth, pPyr_Weight[0][2].m_iHeight);
    Temp_Compare("c:\\tmp\\Dest_3.bin", (float*)pPyr_Weight[0][3].m_pChannel[0], pPyr_Weight[0][3].m_iWidth, pPyr_Weight[0][3].m_iHeight);*/
    //Disp_Part_GPU((float*)pPyr_Weight[0][0].m_pChannel[0], pPyr_Weight[0][0].m_iWidth, 
        //342*2, 12, 10, 1);

    //一个函数干完所有图的金字塔，
    Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    {
        //Pyr_Down_1_float(pPyr_Weight, pPyr_Weight_Header_GPU,iImage_Count);
        //目前为 4100ms
        Pyr_Down_Batch_float(pPyr_Weight, pPyr_Weight_Header_GPU, iImage_Count);
    }
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);


    //////接着算Weight 此处需要  poBlender->dst_roi
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    //此处最慢，12000ms
    _Feed_Blender(Src_Pyr_Laplace, pPyr_Weight, poBlender->dst_pyr_laplace, poBlender->dst_band_weights,
        Src_Pyr_Laplace_Header_GPU, pPyr_Weight_Header_GPU,
        poBlender->dst_pyr_laplace_Header_GPU, poBlender->dst_band_weights_Header_GPU,
        pLTRB_New, poBlender->dst_roi, iImage_Count);

    ////出来酒有了两个大图的金字塔，一个原图的差，一个事权重
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);
    if (oPtr.m_pBuffer)
        Free(oPtr.m_pBuffer);
    if (pWeight)
        Free_GPU(pWeight);
    if (pPyr_Weight[0])
        Free(pPyr_Weight[0]);
    if (pPyr_Weight_Header_GPU[0])
        Free_GPU(pPyr_Weight_Header_GPU[0]);
    if (Src_Pyr_Laplace[0][0].m_pChannel[0])
        Free_GPU(Src_Pyr_Laplace[0][0].m_pChannel[0]);

    return;
}

__global__ void _Normalize_Using_Weight_Map(Image::Part_1 oSource, Image::Part_1 oWeight)
{
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oSource.m_iHeight * oSource.m_iWidth)
        return;
    const float eps = 1e-5f;
    //row[x].x = static_cast<short>(row[x].x / (weight_row[x] + WEIGHT_EPS));
    float fWeight = 1.f / (((float*)oWeight.m_pChannel[0])[iThread_ID] + eps);
    short Value_3[] = { ((short*)oSource.m_pChannel[0])[iThread_ID],
        ((short*)oSource.m_pChannel[1])[iThread_ID],
        ((short*)oSource.m_pChannel[2])[iThread_ID] };

    ((short*)oSource.m_pChannel[0])[iThread_ID] = Value_3[0] * fWeight;
    ((short*)oSource.m_pChannel[1])[iThread_ID] = Value_3[1] * fWeight;
    ((short*)oSource.m_pChannel[2])[iThread_ID] = Value_3[2] * fWeight;

    return;
}

void Normalize_Using_Weight_Map(Image::Part_1 oSource, Image::Part_1 oWeight)
{
    int iSize = oSource.m_iHeight * oSource.m_iWidth;
    dim3 oThread, oGrid;
    oThread.x = Min(iSize, 512);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Normalize_Using_Weight_Map << <oGrid, oThread >> > (oSource, oWeight);

    //不用想分开3片来做了，慢一倍
    return;
}
__global__ void _Py_Up_row_short(Image::Part_1 oSource, Image::Part_1 oMid, int iChannel_Count = 3)
{
    int iThread_ID = GET_THREAD_ID();
    if (iThread_ID >= oMid.m_iHeight * oMid.m_iWidth)
        return;
    short x = iThread_ID % oSource.m_iWidth,
        y = iThread_ID / oSource.m_iWidth;

    unsigned char bHas_Remain_x = oMid.m_iWidth > (oSource.m_iWidth << 1) && (x == oSource.m_iWidth - 1) ? 1 : 0;
    unsigned char bIs_Source_Border = (x == oSource.m_iWidth - 1);
    unsigned char bEven = bIs_Source_Border && ((x << 1) + 1 < oMid.m_iWidth);

    int iPos_m = y * oMid.m_iWidth,
        iPos_s = y * oSource.m_iWidth;
    int iSize_m = oMid.m_iWidth * oMid.m_iHeight;

    short* pMid = &((short*)oMid.m_pChannel[0])[iPos_m + (x << 1)];
    for (short i = 0; i < iChannel_Count; i++, pMid += iSize_m)
    {
        short Mid[3];
        //中间点
        short iPix = ((short*)oSource.m_pChannel[i])[iPos_s + x];
        if (x == 0)
        {//左边两点
            Mid[0] = iPix * 6 + (((short*)oSource.m_pChannel[i])[iPos_s + 1] << 1);
            Mid[1] = (iPix + ((short*)oSource.m_pChannel[i])[iPos_s + 1]) << 2;
        }
        else if (bIs_Source_Border)
        {//右边两点，2对齐情况下
            Mid[0] = ((short*)oSource.m_pChannel[i])[iPos_s + x - 1] + iPix * 7;
            if (bEven)  //源行最后一点乘以二，对应目标也有两点
                Mid[1] = iPix << 3;
        }
        else
        {//一般情况
            Mid[0] = ((short*)oSource.m_pChannel[i])[iPos_s + x - 1] + iPix * 6 + ((short*)oSource.m_pChannel[i])[iPos_s + x + 1];
            Mid[1] = (iPix + ((short*)oSource.m_pChannel[i])[iPos_s + x + 1]) << 2;
        }

        pMid[0] = Mid[0];
        if (!bIs_Source_Border || bEven)
        {
            pMid[1] = Mid[1];
            if (bHas_Remain_x)  //后改由余数
                pMid[2] = Mid[1];
        }
    }
}
__global__ void _Pyr_Up_col_Add_Crop(Image::Part_1 oMid, Image::Part_1 oDest,
    Image::Part_1 oCrop, int iChannel_Count = 3)
{//连Crop都做了
    int iThread_ID = GET_THREAD_ID();
    int iMid_Size = oMid.m_iHeight * oMid.m_iWidth;

    short x = iThread_ID % oDest.m_iWidth,
        y = iThread_ID / oDest.m_iWidth;
    if (iThread_ID >= iMid_Size || x >= oCrop.m_iWidth)
        return;

    unsigned char bHas_Remain_y = oDest.m_iHeight > oMid.m_iHeight * 2 && y == oMid.m_iHeight - 1 ? 1 : 0;
    unsigned char bEven = (y * 2 + 1 < oDest.m_iHeight);

    //重复计算的全部提取出来，快不少
    short* r0 = &((short*)oMid.m_pChannel[0])[iGet_Border_y_GPU(y - 1, oMid.m_iHeight, BORDER_REFLECT101) * oDest.m_iWidth + x],
        * r1 = &((short*)oMid.m_pChannel[0])[iGet_Border_y_GPU(y, oMid.m_iHeight, BORDER_REFLECT) * oDest.m_iWidth + x],
        * r2 = &((short*)oMid.m_pChannel[0])[iGet_Border_y_GPU(y + 1, oMid.m_iHeight, BORDER_REFLECT) * oDest.m_iWidth + x];
    short y1 = y * 2;

    if (y1 >= oCrop.m_iHeight)
        return;

    int iPos_d = y1 * oDest.m_iWidth + x;
    for (int i = 0; i < iChannel_Count; i++, r0 += iMid_Size, r1 += iMid_Size, r2 += iMid_Size)
    {
        int iValue = ((short*)oDest.m_pChannel[i])[iPos_d] + ((*r0 + *r1 * 6 + *r2 + 32) >> 6);
        //if (y * 2 == 0 && x == 35)
        //    printf("%d %d %d\n", ((short*)oDest.m_pChannel[i])[iPos_d], ((*r0 + *r1 * 6 + *r2 + 32) >> 6),iValue);

        //((short*)oDest.m_pChannel[i])[iPos_d] = iValue;
        //((short*)oCrop.m_pChannel[i])[y1 * oCrop.m_iWidth + x] = iValue;
        int iPos_Crop = y1 * oCrop.m_iWidth + x;
        oCrop.m_pChannel[i][iPos_Crop] = oCrop.m_pChannel[3][iPos_Crop] ? iValue : 0;
        //if (y1 == 16 && x == 0)
            //printf("%d %d %d %d\n", iPos_d, ((short*)oDest.m_pChannel[i])[iPos_d], ((*r0 + *r1 * 6 + *r2 + 32) >> 6), iValue);
        if (y1 + 1 >= oCrop.m_iHeight)
            continue;
        if (bEven)
        {
            iValue = ((short*)oDest.m_pChannel[i])[iPos_d + oDest.m_iWidth] + ((((*r1 + *r2) << 2) + 32) >> 6);
            //((short*)oDest.m_pChannel[i])[iPos_d + oDest.m_iWidth] = iValue;    // Clip3(-128, 127, iValue);
            //((short*)oCrop.m_pChannel[i])[(y1 +1) * oCrop.m_iWidth + x] = iValue;
            iPos_Crop = (y1 + 1) * oCrop.m_iWidth + x;
            oCrop.m_pChannel[i][iPos_Crop] = oCrop.m_pChannel[3][iPos_Crop] ? iValue : 0;
            if (y1 + 2 >= oCrop.m_iHeight)
                continue;
            if (bHas_Remain_y)
            {
                iValue = ((short*)oDest.m_pChannel[i])[iPos_d + (oDest.m_iWidth << 1)] + (((short*)oDest.m_pChannel[i])[iPos_d]);
                //((short*)oDest.m_pChannel[i])[iPos_d + (oDest.m_iWidth << 1)] = iValue; // Clip3(-128, 127, iValue);
                //((short*)oCrop.m_pChannel[i])[(y1 + 2) * oCrop.m_iWidth + x] = iValue;
                iPos_Crop = (y1 + 2) * oCrop.m_iWidth + x;
                oCrop.m_pChannel[i][iPos_Crop] = oCrop.m_pChannel[3][iPos_Crop] ? iValue : 0;
            }
        }
    }
    return;
}
void Pyr_Up_Add_Crop_short(Image::Part_1 oSource, Image::Part_1 oDest, short* pMid, Image::Part_1 oCrop)
{//最后一步连Crop都做了
    Image::Part_1 oMid;
    oMid.m_iWidth = oDest.m_iWidth;
    oMid.m_iHeight = oSource.m_iHeight;
    int iSize = oSource.m_iHeight * oSource.m_iWidth;
    oMid.m_pChannel[0] = (unsigned char*)pMid;
    oMid.m_pChannel[1] = oMid.m_pChannel[0] + iSize * 2;
    oMid.m_pChannel[2] = oMid.m_pChannel[1] + iSize * 2;

    dim3 oThread, oGrid;
    //开始分行列进行上采样
    oThread.x = Min(iSize, 256);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Py_Up_row_short << <oGrid, oThread >> > (oSource, oMid);

    iSize = oMid.m_iHeight * oMid.m_iWidth;
    oThread.x = Min(iSize, 512);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Pyr_Up_col_Add_Crop << <oGrid, oThread >> > (oMid, oDest, oCrop);

    return;
}

void Pyr_Up_Add_short(Image::Part_1 oSource, Image::Part_1 oDest, short* pMid)
{//卷积模板大概是3x3，和下采样不一样
    Image::Part_1 oMid;
    oMid.m_iWidth = oDest.m_iWidth;
    oMid.m_iHeight = oSource.m_iHeight;
    int iSize = oSource.m_iHeight * oSource.m_iWidth;
    oMid.m_pChannel[0] = (unsigned char*)pMid;
    oMid.m_pChannel[1] = oMid.m_pChannel[0] + iSize * 2;
    oMid.m_pChannel[2] = oMid.m_pChannel[1] + iSize * 2;

    //Temp_Compare("c:\\tmp\\Py_Up_5_0.bin", (short*)oSource.m_pChannel[0], oSource.m_iWidth, oSource.m_iHeight);
    //Temp_Compare("c:\\tmp\\Py_Up_4_0.bin", (short*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight);

    //Disp_Part_GPU((short*)oSource.m_pChannel[0], oSource.m_iWidth, 17, 0, 3, 2, "Org Source");
    //Disp_Part_GPU((short*)oDest.m_pChannel[0], oDest.m_iWidth, 34, 0, 3, 2,"Org Dest");

    dim3 oThread, oGrid;
    //开始分行列进行上采样
    oThread.x = Min(iSize, 512);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Py_Up_row_short << <oGrid, oThread >> > (oSource, oMid);
    //Disp_Cuda_Error();
    //Disp_Part_GPU((short*)oMid.m_pChannel[0], oMid.m_iWidth, 6, 0, 1, 2);

    iSize = oMid.m_iHeight * oMid.m_iWidth;
    oThread.x = Min(iSize, 512);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Pyr_Up_col_Add << <oGrid, oThread >> > (oMid, oDest);


    //Disp_Cuda_Error();
    //Disp_Part_GPU((short*)oDest.m_pChannel[0], oDest.m_iWidth, 35, 0, 10, 1);
    //Temp_Compare("c:\\tmp\\Py_Up_4_0.bin", (short*)oDest.m_pChannel[0], oDest.m_iWidth, oDest.m_iHeight);

    return;
}

void Restore_Image_From_Laplace(Image::Part_1 Pyr[], Image::Part_1 oDst)
{//将一个金字塔恢复，已经完全对准opencv

    //先分配一个Mid
    Image::Part_1 oUpper = Pyr[0],
        oLower = Pyr[1];    // , oMid;

    short* pMid = (short*)pMalloc_GPU(oLower.m_iHeight * oUpper.m_iWidth * 3 * sizeof(short));

    ////临时装入文件,装入数据就完全对齐了
    //for (int i = 0; i < 6; i++)
    //{
    //    char File[256];
    //    int iSize = Pyr[i].m_iHeight * Pyr[i].m_iWidth * sizeof(short);
    //    for (int j = 0; j < 3; j++)
    //    {
    //        sprintf(File, "c:\\tmp\\normalizeUsingWeightMap_%d_%d.bin", i, j);
    //        bLoad_Raw_Data_GPU(File, &Pyr[i].m_pChannel[j], &iSize);
    //    }
    //}

    for (int i = 5; i > 1; i--)
    {
        //先安排好oMid
        oUpper = Pyr[i - 1];
        oLower = Pyr[i];
        //上采样
        Pyr_Up_Add_short(oLower, oUpper, pMid);
    }
    //Disp_Part_GPU(((short*)Pyr[0].m_pChannel[0]), Pyr[0].m_iWidth, 0, 16, 1, 1);

    oUpper = Pyr[0];
    oLower = Pyr[1];
    Pyr_Up_Add_Crop_short(oLower, oUpper, pMid, oDst);


    //bSave_Image_GPU("c:\\tmp\\1.bmp", oDst);
    //Disp_Cuda_Error();
    //Disp_Part_GPU(((short*)oDst.m_pChannel[0]), oDst.m_iWidth, 0, 16, 1, 1);
    //Temp_Compare<short>("c:\\tmp\\dst_0.bin", (short*)oDst.m_pChannel[0], oDst.m_iWidth, oDst.m_iHeight);
    //Temp_Compare<short>("c:\\tmp\\dst_1.bin", (short*)oDst.m_pChannel[1], oDst.m_iWidth, oDst.m_iHeight);
    //Temp_Compare<short>("c:\\tmp\\dst_2.bin", (short*)oDst.m_pChannel[2], oDst.m_iWidth, oDst.m_iHeight);

    ////验算
    //for (int i = 0; i < 6; i++)
    //{
    //    char File[256];
    //    for (int j = 0; j < 3; j++)
    //    {
    //        sprintf(File, "c:\\tmp\\Py_Up_%d_%d.bin", i, j);
    //        Temp_Compare(File, (short*)Pyr[i].m_pChannel[j], Pyr[i].m_iWidth, Pyr[i].m_iHeight);
    //    }
    //}
    if (pMid)
        Free_GPU(pMid);
    return;
}

__global__ void _Weight_Compare_Crop(Image::Part_1 oWeight, float eps, Image::Part_1 oMask)
{
    int iThread_ID = GET_THREAD_ID();
    short y = iThread_ID / oMask.m_iWidth,
        x = iThread_ID % oMask.m_iWidth;

    oMask.m_pChannel[0][iThread_ID] = ((float*)oWeight.m_pChannel[0])[y * oWeight.m_iWidth + x] > eps ? 255 : 0;
    return;
}
void Weight_Compare_Crop(Image::Part_1 oWeight, float eps, Image oMask)
{
    int iSize = oMask.m_iWidth * oMask.m_iHeight;
    dim3 oThread, oGrid;
    oThread.x = Min(512, iSize);
    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
    _Weight_Compare_Crop << <oGrid, oThread >> > (oWeight, eps, oMask.m_oPart_1);
    /*Disp_Cuda_Error();
    bSave_Image_GPU("c:\\tmp\\2.bmp", oMask);
    Compare_Image("c:\\tmp\\mask.bmp", "c:\\tmp\\2.bmp");*/

    return;
}

//__global__ void _Blend_GPU(Image::Part_1 oDst, Image::Part_1 oMask, Image::Part_1 oResult)
//{
//    int iThread_ID = GET_THREAD_ID();
//    if (iThread_ID > oDst.m_iHeight * oDst.m_iWidth)
//        return;
//
//    int iMask = oMask.m_pChannel[0][iThread_ID] & 1;
//    //short Value[3] = {};
//    //if (iMask)
//    //{
//    //    Value[0] = ((short*)oDst.m_pChannel[0])[iThread_ID];
//    //    Value[1] = ((short*)oDst.m_pChannel[1])[iThread_ID];
//    //    Value[2] = ((short*)oDst.m_pChannel[2])[iThread_ID];
//    //    /*if (Value[0] < 0 || Value[0]>255 ||
//    //        Value[1] < 0 || Value[1]>255 ||
//    //        Value[2] < 0 || Value[2]>255)
//    //    {
//    //        printf("*");
//    //    }*/
//    //}
//    oResult.m_pChannel[0][iThread_ID] = ((short*)oDst.m_pChannel[0])[iThread_ID] *iMask;
//    oResult.m_pChannel[1][iThread_ID] = ((short*)oDst.m_pChannel[1])[iThread_ID] *iMask;
//    oResult.m_pChannel[2][iThread_ID] = ((short*)oDst.m_pChannel[2])[iThread_ID] *iMask;
//}
//void Blend(Image::Part_1 oDst, Image oMask, Image oResult)
//{//最后的混合
//    //1，找出Mask里 =0的点
//    // 2，将dst中对应Mask中的0点，置为0
//    //3，可能要Clip一下     
//    int iSize = oDst.m_iHeight * oDst.m_iWidth;
//    dim3 oThread, oGrid;
//    oThread.x = Min(iSize, 512);
//    oGrid.x = (iSize + oThread.x - 1) / oThread.x;
//    _Blend_GPU << <oGrid, oThread >> > (oDst, oMask.m_oPart_1, oResult.m_oPart_1);
//    //Disp_Cuda_Error();
//    //bSave_Image_GPU("c:\\tmp\\dst.bmp", oResult);
//    return;
//}

template<typename _T>void Blend(Stitch<_T>* poStitch, Blender* poBlender)
{
    //先规格化
    //Disp_Cuda_Error();

    for (int i = 0; i < 6; i++)
    {
        //unsigned long long tStart = iGet_Tick_Count();
        //for(int j=0;j<10000;j++)
        Normalize_Using_Weight_Map(poBlender->dst_pyr_laplace[i], poBlender->dst_band_weights[i]);
        //Disp_Cuda_Error();
        //printf("%lld\n", iGet_Tick_Count() - tStart);

        ////对数据，此处有些许误差是正常的，是浮点运算的问题
        //char File[256];
        //Image::Part_1 oImage = poBlender->dst_pyr_laplace[i];
        /*if(i==5)
        for (int j = 0; j < 3; j++)
        {
            sprintf(File, "c:\\tmp\\normalizeUsingWeightMap_%d_%d.bin", i, j);
            Temp_Compare<short>(File, (short*)oImage.m_pChannel[j], oImage.m_iWidth, oImage.m_iHeight);
        }*/
    }

    /*Image::Part_1 oDst;
    oDst.m_iHeight = poBlender->dst_roi_final[1][1];
    oDst.m_iWidth = poBlender->dst_roi_final[1][0];
    int iSize = oDst.m_iHeight * oDst.m_iWidth;
    oDst.m_pChannel[0] = (unsigned char*)pMalloc_GPU(iSize * 3 *sizeof(short));
    oDst.m_pChannel[1] = oDst.m_pChannel[0] + iSize * sizeof(short);
    oDst.m_pChannel[2] = oDst.m_pChannel[1] + iSize * sizeof(short);*/

    Image oDst;
    Init_Image_GPU(&oDst, poBlender->dst_roi_final[1][0], poBlender->dst_roi_final[1][1], Image::IMAGE_TYPE_BMP, 32);
    poBlender->m_oResult = oDst;

    //short Dst_rc[4] = { 0,0, (short)poBlender->dst_roi_final[1][0],(short)poBlender->dst_roi_final[1][1] };
    //Init_Image_GPU(&poBlender->m_oMask, poBlender->dst_roi_final[1][0], poBlender->dst_roi_final[1][1], Image::IMAGE_TYPE_BMP, 8);
    Attach_Buffer(&poBlender->m_oMask, oDst.m_pChannel[3], poBlender->dst_roi_final[1][0], poBlender->dst_roi_final[1][1],
        1, Image::IMAGE_TYPE_BMP);

    Weight_Compare_Crop(poBlender->dst_band_weights[0], 1e-5f, poBlender->m_oMask);
    //bSave_Image_GPU("c:\\tmp\\1.bmp", poBlender->m_oMask);

    //接着做上采样
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Restore_Image_From_Laplace(poBlender->dst_pyr_laplace, oDst.m_oPart_1);
    if (poBlender->m_pBuffer)
    {
        Free_GPU(poBlender->m_pBuffer);
        poBlender->m_pBuffer = NULL;
    }

    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //bSave_Image_GPU("c:\\tmp\\1.bmp", oDst);
    ////最后的合成
    //Image oResult;
    //Init_Image_GPU(&oResult, oDst.m_iWidth, oDst.m_iHeight, Image::IMAGE_TYPE_BMP, 24);
    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    //Blend(oDst, poBlender->m_oMask, oResult);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    return;
}
template<typename _T>void Compose_Panorama(Stitch<_T>* poStitch,Image *poResult)
{
    int i;
    Stitch<_T> oStitch = *poStitch;
    Point_Cloud<float> oPC;
    Init_Point_Cloud(&oPC, 100000, 1);

    //以下做小规模投影
    for (i = 0; i < oStitch.m_iImage_Count; i++)
    {
        _T K[3 * 3];
        Camera<_T> oCamera = oStitch.m_pCamera[i];

        //此处实际上是将K的位移与投影乘以一个scale，变成另一张图
        //问题是加上z投影到哪?
        Set_K(&oStitch, oCamera.K, K, poStitch->seam_work_aspect);
        _T fScale = oStitch.warped_image_scale * oStitch.seam_work_aspect;

        //unsigned long long tStart = iGet_Tick_Count();
        //for(int j=0;j<10000;j++)
        {//合二为一
            Warp_2<_T>(oStitch.m_pSeam_Est[i], oStitch.m_pMask[i], K, oCamera.R, fScale,
                &oStitch.m_pImage_Warp[i], &oStitch.m_pMasks_Warped[i],
                INTER_LINEAR, INTER_NEAREST, oStitch.m_pCorner[i]);
            //bSave_Image_GPU("c:\\tmp\\2.bmp", oStitch.m_pMasks_Warped[i]);
            //Free_Image_GPU(&oStitch.m_pImage_Warp[i]);
            //Free_Image_GPU(&oStitch.m_pMasks_Warped[i]);
        }
        //Disp_Cuda_Error();
        //printf("%lld\n", iGet_Tick_Count() - tStart);
    }
    //投影完了以后还要把数据抄到GPU上
    cudaMemcpy(oStitch.m_pCorner_GPU, oStitch.m_pCorner, oStitch.m_iImage_Count * 2 * 2 * sizeof(int), cudaMemcpyHostToDevice);

    //Disp_Cuda_Error();
    //测试代码，装入数据
    /*bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\images_warped_0.bmp", &poStitch->m_pImage_Warp[0]);
    bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\images_warped_1.bmp", &poStitch->m_pImage_Warp[1]);
    bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\images_warped_2.bmp", &poStitch->m_pImage_Warp[2]);
    bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\images_warped_3.bmp", &poStitch->m_pImage_Warp[3]);*/

    //将前面的投影图像全部拿出来分块，每一块与别人求交集
    //根据交集算个补光值，一块一值。对各图补光？
    Feed<_T>(&oStitch);

    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_0.bmp", &poStitch->m_pImage_Warp[0]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_1.bmp", &poStitch->m_pImage_Warp[1]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_2.bmp", &poStitch->m_pImage_Warp[2]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\image_warp_f_3.bmp", &poStitch->m_pImage_Warp[3]);
    //对各个投影图进行求交？
    //bSave_Image_GPU("c:\\tmp\\1.bmp", oStitch.m_pImage_Warp[0]);
    //此处理应做最大流，暂时未做，慢慢改良
    Find<_T>(&oStitch);

    //至此，按照opencv得观点，可以释放部分内存了
    Free_Partial(&oStitch);

    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\0.bmp", &oStitch.m_pImage_Source[0]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\1.bmp", &oStitch.m_pImage_Source[1]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\2.bmp", &oStitch.m_pImage_Source[2]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\3.bmp", &oStitch.m_pImage_Source[3]);
    //bSave_Image_GPU("c:\\tmp\\1.bmp", oStitch.m_pImage_Source[0]);
    //int bIs_Compose_Scale_Set = 0;
    //Camera<_T> Cameras_Scaled[2];
    _T compose_work_aspect = oStitch.compose_scale / oStitch.work_scale;
    _T warp_scale = oStitch.warped_image_scale * compose_work_aspect;
    Image* pBlend_Warp = (Image*)pMalloc(poStitch->m_iImage_Count * sizeof(Image));
    Image* pBlend_Warp_Header_GPU = (Image*)pMalloc_GPU(poStitch->m_iImage_Count * sizeof(Image));

    for (int i = 0; i < oStitch.m_iImage_Count; i++)
    {
        _T K[3 * 3];
        Camera<_T> oCamera = oStitch.m_pCamera[i];

        //此处实际上是将K的位移与投影乘以一个scale，变成另一张图
        //问题是加上z投影到哪?
        Set_K(&oStitch, oCamera.K, K, compose_work_aspect);
        //Disp(K, 3, 3, "K");


        //借用原来poStitch->m_pWarp
        Image oImg_Warp;
        //Disp_Cuda_Error();
        //unsigned long long tStart = iGet_Tick_Count();
        //for(int j=0;j<10000;j++)
        {
            Warp_3<_T>(poStitch->m_pImage_Source[i], K, oCamera.R, warp_scale,
                &oImg_Warp, INTER_LINEAR, INTER_NEAREST, poStitch->m_pCorner[i],
                BORDER_REFLECT, BORDER_CONSTANT, NULL);
            //Free_Image_GPU(&oImg_Warp);
        }
        //Disp_Cuda_Error();
        //printf("%lld\n", iGet_Tick_Count() - tStart);
        pBlend_Warp[i] = oImg_Warp;
        //轮到膨胀收缩了
    }
    cudaMemcpy(pBlend_Warp_Header_GPU, pBlend_Warp, poStitch->m_iImage_Count * sizeof(Image), cudaMemcpyHostToDevice);

    ////装入Warp图像，测补光
    //bLoad_Image_GPU("C:\\tmp\\Warp_0.bmp", &pBlend_Warp[0]);
    //bLoad_Image_GPU("C:\\tmp\\Warp_1.bmp", &pBlend_Warp[1]);
    //bLoad_Image_GPU("C:\\tmp\\Warp_2.bmp", &pBlend_Warp[2]);
    //bLoad_Image_GPU("C:\\tmp\\Warp_3.bmp", &pBlend_Warp[3]);
    //补光
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    Block_Compensate(&oStitch, pBlend_Warp, pBlend_Warp_Header_GPU);
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    //bSave_Image_GPU("c:\\tmp\\0.bmp", &pBlend_Warp_Header_GPU[0]);
    //Compare_Image("c:\\tmp\\Warp_Comp_0.bmp", "c:\\tmp\\0.bmp",1);

    ////装入Warp图像
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Warp_0.bmp", &poStitch->m_pImage_Warp[0]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Warp_1.bmp", &poStitch->m_pImage_Warp[1]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Warp_2.bmp", &poStitch->m_pImage_Warp[2]);
    //bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Warp_3.bmp", &poStitch->m_pImage_Warp[3]);
    //cudaMemcpy(poStitch->m_pImage_Warp[0].m_pChannel[3], poStitch->m_pImage_Warp[0].m_pChannel[0], 
    //    poStitch->m_pImage_Warp[0].m_iWidth * poStitch->m_pImage_Warp[0].m_iHeight, cudaMemcpyDeviceToDevice);
    //cudaMemcpy(poStitch->m_pImage_Warp[1].m_pChannel[3], poStitch->m_pImage_Warp[1].m_pChannel[0],
    //    poStitch->m_pImage_Warp[1].m_iWidth * poStitch->m_pImage_Warp[1].m_iHeight, cudaMemcpyDeviceToDevice);
    //cudaMemcpy(poStitch->m_pImage_Warp[2].m_pChannel[3], poStitch->m_pImage_Warp[2].m_pChannel[0],
    //    poStitch->m_pImage_Warp[2].m_iWidth * poStitch->m_pImage_Warp[2].m_iHeight, cudaMemcpyDeviceToDevice);
    //cudaMemcpy(poStitch->m_pImage_Warp[3].m_pChannel[3], poStitch->m_pImage_Warp[3].m_pChannel[0],
    //    poStitch->m_pImage_Warp[3].m_iWidth * poStitch->m_pImage_Warp[3].m_iHeight, cudaMemcpyDeviceToDevice);

    Dilate_GPU(poStitch->m_pImage_Warp, poStitch->m_pImage_Warp_Header_GPU, poStitch->m_iImage_Count);

    /*Draw_Point<float>(&oPC, 0, 0, 0);
    bSave_PLY("c:\\tmp\\1.ply", oPC);*/

    //接着将Image_Warp中的Mask通道缩放到Warp那么大，再与Warp的Mask进行 Bitwise And
    //此处很慢，3500ms
    Resize_Bitwise_And(poStitch, pBlend_Warp, pBlend_Warp_Header_GPU);
    Blender oBlender = {};

    //Disp_Cuda_Error();
    //unsigned long long tStart = iGet_Tick_Count();
    //for(int i=0;i<10000;i++)
    {
        //此处要修改，过早分配内存，忘了后面咋用，后面再分配聚合度更高
        Prepare_Blender_1(&oStitch, &oBlender);
    }
    //Disp_Cuda_Error();
    //printf("%lld\n", iGet_Tick_Count() - tStart);

    ////先造数据，装入Mask_Warp
    //bLoad_Comp_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Mask_Warp_0.bmp", pBlend_Warp[0], 3);
    //bLoad_Comp_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Mask_Warp_1.bmp", pBlend_Warp[1], 3);
    //bLoad_Comp_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Mask_Warp_2.bmp", pBlend_Warp[2], 3);
    //bLoad_Comp_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\Mask_Warp_3.bmp", pBlend_Warp[3], 3);

    //bSave_Image_GPU("c:\\tmp\\1.bmp", pBlend_Warp[0]);
    Feed_Blender(pBlend_Warp, pBlend_Warp_Header_GPU, poStitch, &oBlender);

    //最后做Blend
    Blend(poStitch, &oBlender);
    *poResult = oBlender.m_oResult;

    //最后释放
    if (pBlend_Warp)
    {
        for (int i = 0; i < oStitch.m_iImage_Count; i++)
            Free_Image_GPU(&pBlend_Warp[i]);
        Free(pBlend_Warp);
    }

    if (pBlend_Warp_Header_GPU)
        Free_GPU(pBlend_Warp_Header_GPU);
    if (oBlender.m_pBuffer)
        Free_GPU(oBlender.m_pBuffer);
    //Free_Image_GPU(&oBlender.m_oResult);

    *poStitch = oStitch;
    Free_Point_Cloud(&oPC);
    return;
}



void Pry_Down_Test()
{
    Image oSource = {}, oDest;
    //bLoad_Image_GPU("C:\\tmp\\1.bmp", &oSource);
    bLoad_Image_GPU("C:\\tmp\\temp\\stitch\\1706x1279\\1.bmp", &oSource);
    Init_Image_GPU(&oDest, (oSource.m_iWidth + 1) >> 1, (oSource.m_iHeight + 1) >> 1, Image::IMAGE_TYPE_BMP, oSource.m_iBit_Count);

    Disp_Cuda_Error();
    unsigned long long tStart = iGet_Tick_Count();
    for (int i = 0; i < 10000; i++)
        Pyr_Down_GPU(oSource, oDest);
    Disp_Cuda_Error();
    printf("%lld\n", iGet_Tick_Count() - tStart);

    bSave_Image_GPU("c:\\tmp\\3.bmp", oDest);
    Compare_Image("c:\\tmp\\2.bmp", "c:\\tmp\\3.bmp");

    Free_Image_GPU(&oSource);
    Free_Image_GPU(&oDest);
    return;
}

void Pyr_Up_Test()
{//金字塔上采样实验
    Image oSource = {}, oDest;
    bLoad_Image_GPU("c:\\tmp\\1.bmp", &oSource);
    Init_Image_GPU(&oDest, oSource.m_iWidth * 2 - 1, oSource.m_iHeight * 2 - 1, Image::IMAGE_TYPE_BMP, oSource.m_iBit_Count);
    Disp_Cuda_Error();
    unsigned long long tStart = iGet_Tick_Count();
    for (int i = 0; i < 10000; i++)
        Pyr_Up_GPU(oSource, oDest);
    Disp_Cuda_Error();
    printf("%lld\n", iGet_Tick_Count() - tStart);

    Free_Image_GPU(&oSource);
    Free_Image_GPU(&oDest);
    return;
}
void Copy_Make_Border_Test()
{
    int iLeft = 5, iTop = 11, iRight = 125 + 2, iBottom = 85;
    Image oSource = {}, oDest;
    bLoad_Image_GPU("c:\\tmp\\1.bmp", &oSource);
    Init_Image_GPU(&oDest, oSource.m_iWidth + iLeft + iRight, oSource.m_iHeight + iTop + iBottom, Image::IMAGE_TYPE_BMP, oSource.m_iBit_Count);
    Copy_Make_Border_GPU(oSource, oDest, iLeft, iTop, iRight, iBottom);
    return;
}

static void Test_1()
{
    typedef double _T;
    Stitch<_T> oStitch = {};
    //由于将视觉分离，此处装入相机参数
    Temp_Load_Camera("data\\camera.bin", &oStitch);

    //初始化拼接器
    Init_Stitch<double>(&oStitch, 1706, 1279);

    //装入原图，以后替代此处
    bLoad_Image(oStitch);
    Resize_Seam_Image(&oStitch);

    Disp_Cuda_Error();
    unsigned long long tStart = iGet_Tick_Count();
    Image oResult;
    Compose_Panorama(&oStitch, &oResult);
    Disp_Cuda_Error();
    Free_Stitch(&oStitch);
    printf("Overall %lld ms\n", iGet_Tick_Count() - tStart);

    bSave_Image_GPU("c:\\tmp\\1.bmp", oResult);
    Free_Image_GPU(&oResult);
}
int main()
{
    Init_Env_All();
    //Copy_Make_Border_Test();
    //Gauss_Test();
    Test_1();
    //Test_2();
    //Pry_Down_Test();
    //Pyr_Up_Test();
    Free_Env_All();
    return 0;
}

