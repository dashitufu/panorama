// Pano.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include <iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include "stdio.h"

#include "Common.h"
#include "Common_cuda.cuh"
#include "Matrix.h"
//#include "Image.h"
#include "image_cuda.h"
#include "Panorama.h"

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







