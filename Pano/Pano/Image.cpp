//��潫����BitMap������Image��ͨ���ͽṹ
#pragma once

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "sys/timeb.h"
#include "Common.h"
#include "Image.h"
#include "matrix.h"	//���Ҫ���ں��棬��ΪMem_MgrҪ����ǰ��
#include "immintrin.h"

extern "C"
{
#include "Buddy_System.h"
}

int bStricmp(char* pStr_0, char* pStr_1)
{//����Linux��Windows��һ�£����ȸ�һ��ͬһ�޴�Сд�ȽϺ���
#ifndef WIN32
#include "strings.h"
	return strcasecmp(pStr_0, pStr_1);
#else
	return _stricmp(pStr_0, pStr_1);
#endif
}
float fNaN()
{
	unsigned int iRet = 0x7FFFFFFF;
	return *(float*)(&iRet);
}

void Init_Image(Image* poImage, int iWidth, int iHeight, Image::Type iType, int iBit_Count)
{
	int i, iSize, iSize_With_Remain;;
	poImage->m_iWidth = iWidth;
	poImage->m_iHeight = iHeight;
	poImage->m_iBit_Count = iBit_Count;
	poImage->m_iChannel_Count = iBit_Count >> 3;
	poImage->m_iMem_Src = Mem_Src::CPU;
	poImage->m_bIs_Attached = 0;
	iSize = iWidth * iHeight;

	poImage->m_iMax_Buffer_Size = iSize * poImage->m_iChannel_Count;

	//������أ�1��������64�ֽ����ϣ�2��128�ֽڱ������Ա�GPU����
	iSize_With_Remain = ((poImage->m_iMax_Buffer_Size + 127) / 128) * 128;
	if (iSize_With_Remain - poImage->m_iMax_Buffer_Size < 64)
		iSize_With_Remain += 128;

	//poImage->m_pChannel[0] = (unsigned char*)malloc(iSize_With_Remain);
	poImage->m_pChannel[0] = (unsigned char*)pMalloc(iSize_With_Remain);
	if (!poImage->m_pChannel[0])
	{
		printf("Fail to malloc in Init_Image\n");
		*poImage = {};
		return;
	}

	for (i = 1; i < poImage->m_iChannel_Count; i++)
		poImage->m_pChannel[i] = poImage->m_pChannel[i - 1] + iSize;
	for (; i < 4; i++)
		poImage->m_pChannel[i] = NULL;

	poImage->m_iImage_Type = iType;
	poImage->m_pBuffer = poImage->m_oBitMap.m_pRGBA[0];

	Set_Color(*poImage);
}

int bLoad_Image(const char* pcFile, BitMap_8Bit* poBitMap, int bNeed_Malloc,Mem_Mgr *poMem_Mgr)
{//װ��BMPͼ��δ���Ż��������������Բ�������
	BitMap_8Bit::BITMAPFILEHEADER oFile_Header;
	BitMap_8Bit::BITMAPINFOHEADER oInfo_Header;
	unsigned char* pLine = NULL, * pCur, * RGBA[4], * RGBA_1[4];
	int i, j, k, iSize, iResult, iSize_Width_Remain, iRet = 1;
	int iCount;
	FILE* pFile = fopen(pcFile, "rb");
	if (!pFile)
	{
		printf("Fail to open:%s\n", pcFile);
		return 0;
	}
	iResult = (int)fread(&oFile_Header, 1, sizeof(BitMap_8Bit::BITMAPFILEHEADER), pFile);
	iResult = (int)fread(&oInfo_Header, 1, sizeof(BitMap_8Bit::BITMAPINFOHEADER), pFile);
	if (!iResult)
	{
		printf("Fail to read header in bLoad_Image\n");
		iRet = 0;
		goto END;
	}
	if (/*oInfo_Header.biBitCount!=1 &&*/ oInfo_Header.biBitCount != 32 && oInfo_Header.biBitCount != 24 && oInfo_Header.biBitCount != 8)
	{//��ʱ����������λͼ
		printf("Invalid Bit-Depth in bitmap:%d\n", oInfo_Header.biBitCount);
		iRet = 0;
		goto END;
	}
	poBitMap->m_iBit_Count = oInfo_Header.biBitCount;
	poBitMap->m_iChannel_Count = poBitMap->m_iBit_Count / 8;
	poBitMap->m_iWidth = oInfo_Header.biWidth;
	poBitMap->m_iHeight = oInfo_Header.biHeight;
	poBitMap->m_iMem_Src = Mem_Src::CPU;
	if (bNeed_Malloc)
		poBitMap->m_bIs_Attached = 0;

	iSize = oInfo_Header.biWidth * oInfo_Header.biHeight;
	poBitMap->m_iMax_Buffer_Size = iSize * poBitMap->m_iChannel_Count;

	//������أ�1��������64�ֽ����ϣ�2��128�ֽڱ������Ա�GPU����
	iSize_Width_Remain = ((poBitMap->m_iMax_Buffer_Size + 127) / 128) * 128;
	if (iSize_Width_Remain - poBitMap->m_iMax_Buffer_Size < 64)
		iSize_Width_Remain += 128;

	if (bNeed_Malloc)
	{
		if (poMem_Mgr)
			poBitMap->m_pRGBA[0] = (unsigned char*)pMalloc(poMem_Mgr, iSize_Width_Remain);
		else
			poBitMap->m_pRGBA[0] = (unsigned char*)pMalloc(iSize_Width_Remain);
	}

	for (i = 1; i < poBitMap->m_iChannel_Count; i++)
		poBitMap->m_pRGBA[i] = poBitMap->m_pRGBA[i - 1] + iSize;
	for (; i < 4; i++)
		poBitMap->m_pRGBA[i] = NULL;

	if(poBitMap->m_iBit_Count>=8)
		iSize = ((oInfo_Header.biWidth * poBitMap->m_iChannel_Count + 3) / 4) * 4;

	if(poMem_Mgr)
		pLine = (unsigned char*)pMalloc(poMem_Mgr,iSize + 4);
	else
		pLine = (unsigned char*)pMalloc(iSize + 4);
	if (!poBitMap->m_pRGBA[0] || !pLine)
	{
		printf("Fail to allocate memory in Init_Image\n");
		iRet = 0;
		goto END;
	}

	iCount = Min(3, poBitMap->m_iChannel_Count);
	for (i = 0; i < iCount; i++)//�ļ�����BGR����
		RGBA_1[i] = poBitMap->m_pRGBA[iCount - i - 1] + poBitMap->m_iWidth * poBitMap->m_iHeight - oInfo_Header.biWidth;
	if (poBitMap->m_iChannel_Count == 4)	//��Alpha
		RGBA_1[3] = poBitMap->m_pRGBA[3] + poBitMap->m_iWidth * poBitMap->m_iHeight - oInfo_Header.biWidth;

	fseek(pFile, oFile_Header.bfOffBits, SEEK_SET);
	for (i = 0; i < oInfo_Header.biHeight; i++)
	{
		iResult = (int)fread(pLine, 1, iSize, pFile);
		pCur = pLine;
		for (k = 0; k < poBitMap->m_iChannel_Count; k++)
			RGBA[k] = RGBA_1[k];

		for (j = 0; j < oInfo_Header.biWidth; j++)
		{
			for (k = 0; k < poBitMap->m_iChannel_Count; k++)
			{
				*RGBA[k] = *pCur++;
				RGBA[k]++;
			}
		}
		for (k = 0; k < poBitMap->m_iChannel_Count; k++)
			RGBA_1[k] -= oInfo_Header.biWidth;
	}

END:
	fclose(pFile);
	if (pLine)
	{
		if (poMem_Mgr)
			Free(poMem_Mgr, pLine);
		else
			Free(pLine);
	}
	if (!iRet)
		Free(poBitMap->m_pRGBA[0]);

	return iRet;
}

void Init_Image(YUV_Image_420_8Bit* poYUV_Image, int iWidth, int iHeight)
{
	int iSize, iSize_With_Remain;
	int iChroma_Width, iChroma_Height;

	poYUV_Image->m_iMem_Src = Mem_Src::CPU;
	poYUV_Image->m_iWidth = iWidth;
	poYUV_Image->m_iHeight = iHeight;
	iChroma_Width = (iWidth + 1) >> 1;
	iChroma_Height = (iHeight + 1) >> 1;

	iSize = poYUV_Image->m_iWidth * poYUV_Image->m_iHeight;
	poYUV_Image->m_iMax_Buffer_Size = iSize + iChroma_Width * iChroma_Height * 2;

	//������أ�1��������64�ֽ����ϣ�2��128�ֽڱ������Ա�GPU����
	iSize_With_Remain = ((poYUV_Image->m_iMax_Buffer_Size + 127) / 128) * 128;
	if (iSize_With_Remain - poYUV_Image->m_iMax_Buffer_Size < 64)
		iSize_With_Remain += 128;

	//�˴����ǲ���·������Ҫ�д��������ڴ�ռ�ķ�������Ӧ���������ڴ�
	poYUV_Image->m_pYUV[0] = (unsigned char*)malloc(iSize_With_Remain);

	if (!poYUV_Image->m_pYUV[0])
	{
		printf("Fail to allocate memory in Init_Image\n");
		poYUV_Image->m_pYUV[1] = poYUV_Image->m_pYUV[2] = NULL;
		return;
	}
	poYUV_Image->m_pYUV[1] = poYUV_Image->m_pYUV[0] + iSize;
	poYUV_Image->m_pYUV[2] = poYUV_Image->m_pYUV[1] + iChroma_Width * iChroma_Height;
}

int bLoad_Image(const char* pcFile, YUV_Image_420_8Bit* poImage, int iWidth, int iHeight, int bNeed_Malloc, int iFrame)
{
	FILE* pFile = fopen(pcFile, "rb");
	int iResult, bRet = 1, iSize = iWidth * iHeight + ((iWidth + 1) >> 1) * ((iHeight + 1) >> 1) * 2;
	if (!pFile)
	{
		printf("Fail to open:%s\n", pcFile);
		return 0;
	}
	if (bNeed_Malloc)
		Init_Image(poImage, iWidth, iHeight);

	if (!poImage->m_pYUV[0])
	{
		bRet = 0;
		goto END;
	}

#ifdef WIN32
	_fseeki64(pFile, (unsigned long long)iSize * iFrame, SEEK_SET);
#else
	fseeko64(pFile, (unsigned long long)iSize * iFrame, SEEK_SET);
#endif

	iResult = (int)fread(poImage->m_pYUV[0], 1, iSize, pFile);
	if (iResult < iSize)
	{
		printf("Error in reading file:%s iSize:%d iResult:%d\n", pcFile, iSize, iResult);
		bRet = 0;
	}

END:
	if (pFile)
		fclose(pFile);
	return bRet;
}
void Get_Image_Info(const char* pcFile, Image* poImage)
{//���BMP��ͼ������
	if (bStricmp((char*)pcFile + strlen(pcFile) - 3, (char*)"bmp") == 0)
	{
		BitMap_8Bit::BITMAPFILEHEADER oFile_Header;
		BitMap_8Bit::BITMAPINFOHEADER oInfo_Header;
		//unsigned char* pLine = NULL, * pCur, *RGBA[4], * RGBA_1[4];
		int iResult, iSize;
		FILE* pFile = fopen(pcFile, "rb");
		if (!pFile)
		{
			printf("Fail to open:%s\n", pcFile);
			return;
		}
		iResult = (int)fread(&oFile_Header, 1, sizeof(BitMap_8Bit::BITMAPFILEHEADER), pFile);
		iResult = (int)fread(&oInfo_Header, 1, sizeof(BitMap_8Bit::BITMAPINFOHEADER), pFile);
		fclose(pFile);
		if (!iResult)
		{
			printf("Fail to read header in bLoad_Image\n");			
			goto END;
		}
		if (oInfo_Header.biBitCount != 32 && oInfo_Header.biBitCount != 24 && oInfo_Header.biBitCount != 8)
		{//��ʱ����������λͼ
			printf("Invalid Bit-Depth in bitmap:%d\n", oInfo_Header.biBitCount);
			goto END;
		}
		poImage->m_iBit_Count = oInfo_Header.biBitCount;
		poImage->m_iChannel_Count = poImage->m_iBit_Count / 8;
		poImage->m_iWidth = oInfo_Header.biWidth;
		poImage->m_iHeight = oInfo_Header.biHeight;
		poImage->m_iMem_Src = Mem_Src::CPU;
		
		iSize = oInfo_Header.biWidth * oInfo_Header.biHeight;
		poImage->m_iMax_Buffer_Size = iSize * poImage->m_iChannel_Count;
		poImage->m_iImage_Type = Image::IMAGE_TYPE_BMP;
	}else
	{
		printf("Fail to open file:%s\n", pcFile);
		return;
	}
END:

	return;
}
int bLoad_Image(const char* pcFile, Image* poImage, int iWidth, int iHeight, int iFrame, int bNeed_Malloc, Mem_Mgr* poMem_Mgr)
{//ֱ��װ��Image,����ʲô����
	int iResult;
	YUV_Image_420_8Bit oImage;
	if(bNeed_Malloc)
		*poImage = { };

	if (bStricmp((char*)pcFile + strlen(pcFile) - 3, (char*)"bmp") == 0)
	{
		iResult = bLoad_Image(pcFile, &poImage->m_oBitMap, bNeed_Malloc, poMem_Mgr);
		poImage->m_iImage_Type = Image::IMAGE_TYPE_BMP;
		return iResult;
	}
	else if (bStricmp((char*)pcFile + strlen(pcFile) - 3, (char*)"yuv") == 0)
	{
		iResult = bLoad_Image(pcFile, &oImage, iWidth, iHeight, bNeed_Malloc, iFrame);
		poImage->m_iImage_Type = Image::IMAGE_TYPE_YUV_444;
	}
	else
	{
		printf("Fail to open file extension:%s\n", pcFile);
		iResult = 0;
	}

	if (!iResult || poImage->m_iImage_Type == Image::IMAGE_TYPE_BMP)
		return iResult;

	//YUV��Ҫ��һ������
	if(bNeed_Malloc)
		Init_Image(poImage, iWidth, iHeight, Image::IMAGE_TYPE_YUV_444, 24);

	if (!poImage->m_oYUV.m_pYUV[0])
	{
		iResult = 0;
		goto END;
	}
	

END:
	return iResult;
}

void Set_Color(Image oImage, int R, int G, int B, int A)
{
	int iSize = oImage.m_iWidth * oImage.m_iHeight;
	int Color[3];
	if (oImage.m_iImage_Type == Image::IMAGE_TYPE_BMP)
		Color[0] = R, Color[1] = G, Color[2] = B;
	else
		_RGB_2_YUV(R, G, B, Color[0], Color[1], Color[2]);

	if (oImage.m_pChannel[0])
		memset(oImage.m_pChannel[0], Color[0], iSize);
	if (oImage.m_pChannel[1])
		memset(oImage.m_pChannel[1], Color[1], iSize);
	if (oImage.m_pChannel[2])
		memset(oImage.m_pChannel[2], Color[2], iSize);
	
	if (oImage.m_pChannel[3])
		memset(oImage.m_pChannel[3], A, iSize);

	return;
}

void Set_BitMap_Header(BitMap_8Bit::BITMAPFILEHEADER* poFile_Header, BitMap_8Bit::BITMAPINFOHEADER* poInfo_Header, int iWidth, int iHeight, int iChannel_Count)
{//����BMPͷ�����ͷ����û�ã�ֻ���ļ������Žṹ˵�����ã�����ʱ��ȫ�˷ѿռ�
	//int iSize;
	BitMap_8Bit::BITMAPFILEHEADER oFile_Header;
	BitMap_8Bit::BITMAPINFOHEADER oInfo_Header;

	oFile_Header.bfReserved1 = oFile_Header.bfReserved2 = 0;
	oFile_Header.bfType = 19778;	//BM
	oFile_Header.bfOffBits = 54;

	oInfo_Header.biPlanes = 1;
	oInfo_Header.biBitCount = iChannel_Count * 8;
	oInfo_Header.biClrImportant = 0;
	oInfo_Header.biClrUsed = 0;
	oInfo_Header.biCompression = 0;
	oInfo_Header.biHeight = iHeight;
	oInfo_Header.biWidth = iWidth;
	oInfo_Header.biSize = 40;
	oInfo_Header.biYPelsPerMeter = oInfo_Header.biXPelsPerMeter = 0;
	oInfo_Header.biSizeImage = oInfo_Header.biWidth * oInfo_Header.biHeight * oInfo_Header.biBitCount / 8;
	oFile_Header.bfSize = oInfo_Header.biSizeImage + oFile_Header.bfOffBits;
	//iSize = oInfo_Header.biWidth * oInfo_Header.biHeight;
	if (iChannel_Count == 1)
	{
		oFile_Header.bfOffBits = sizeof(BitMap_8Bit::BITMAPINFOHEADER) + sizeof(BitMap_8Bit::BITMAPFILEHEADER) + 1024;	//�˴�Ҫ���ϵ�ɫ��
		oFile_Header.bfSize = oInfo_Header.biSizeImage + oFile_Header.bfOffBits;
	}

	*poFile_Header = oFile_Header;
	*poInfo_Header = oInfo_Header;
}

int bSave_Image(const char* pcFile, BitMap_8Bit oBitMap)
{//���̣�ϣ���ܴ�һͳ����1ͨ���Ķ�ͳһ����������ڣ�δ��ɣ���Ƿ1ͨ�������
	FILE* pFile = fopen(pcFile, "wb");
	int iSize, i, j, k, iResult, iRet = 1;
	unsigned char* pLine = NULL, * pCur;
	unsigned char* RGBA[4], * RGBA_1[4];
	unsigned char Palatte[1024];

	BitMap_8Bit::BITMAPFILEHEADER oFile_Header;
	BitMap_8Bit::BITMAPINFOHEADER oInfo_Header;

	if (!pFile)
		return 0;
	Set_BitMap_Header(&oFile_Header, &oInfo_Header, oBitMap.m_iWidth, oBitMap.m_iHeight, oBitMap.m_iChannel_Count);

	iSize = oInfo_Header.biWidth * oInfo_Header.biHeight;
	if (!iSize)
		return 0;

	if (oBitMap.m_iBit_Count < 32)
	{
		for (i = 0; i < oBitMap.m_iChannel_Count; i++)
			RGBA_1[i] = oBitMap.m_pRGBA[oBitMap.m_iChannel_Count - i - 1] + iSize - oInfo_Header.biWidth;
	}
	else
	{
		for (i = 0; i < 3; i++)
			RGBA_1[i] = oBitMap.m_pRGBA[3 - i - 1] + iSize - oInfo_Header.biWidth;
		RGBA_1[3] = oBitMap.m_pRGBA[3] + iSize - oInfo_Header.biWidth;
	}


	iSize = ((oBitMap.m_iWidth * oBitMap.m_iChannel_Count + 3) / 4) * 4;
	pLine = (unsigned char*)malloc(iSize);
	memset(pLine + iSize - 4, 0, 4);

	iResult = (int)fwrite(&oFile_Header, (size_t)1, (size_t)sizeof(oFile_Header), pFile);
	iResult = (int)fwrite(&oInfo_Header, (size_t)1, (size_t)sizeof(oInfo_Header), pFile);
	if (oBitMap.m_iChannel_Count == 1)
	{
		//��ɫ��
		for (pCur = Palatte, i = 0; i < 256; i++)
		{
			pCur[0] = pCur[1] = pCur[2] = i; pCur[3] = 0;
			pCur += 4;
		}
		iResult = (int)fwrite(Palatte, 1, 1024, pFile);
	}

	for (i = 0; i < oInfo_Header.biHeight; i++)
	{
		pCur = pLine;
		for (k = 0; k < oBitMap.m_iChannel_Count; k++)
			RGBA[k] = RGBA_1[k];

		for (j = 0; j < oInfo_Header.biWidth; j++)
		{
			for (k = 0; k < oBitMap.m_iChannel_Count; k++)
			{
				*pCur++ = *RGBA[k];
				RGBA[k]++;
			}
		}
		iResult = (int)fwrite(pLine, (size_t)1, (size_t)iSize, pFile);
		if (iResult != iSize)
		{
			iRet = 0;
			goto END;
		}
		for (k = 0; k < oBitMap.m_iChannel_Count; k++)
			RGBA_1[k] -= oInfo_Header.biWidth;
	}

END:
	if (pLine)
		free(pLine);
	if (pFile)
		fclose(pFile);
	return iRet;
}
void Init_Image(BitMap_8Bit* poBitMap, int iWidth, int iHeight, int iBit_Count)
{//��ʼ��BMP
	int i, iSize, iSize_With_Remain;;
	poBitMap->m_iWidth = iWidth;
	poBitMap->m_iHeight = iHeight;
	poBitMap->m_iBit_Count = iBit_Count;
	poBitMap->m_iChannel_Count = iBit_Count >> 3;
	poBitMap->m_iMem_Src = Mem_Src::CPU;
	poBitMap->m_bIs_Attached = 0;
	iSize = iWidth * iHeight;

	poBitMap->m_iMax_Buffer_Size = iSize * poBitMap->m_iChannel_Count;

	//������أ�1��������64�ֽ����ϣ�2��128�ֽڱ������Ա�GPU����
	iSize_With_Remain = ((poBitMap->m_iMax_Buffer_Size + 127) / 128) * 128;
	if (iSize_With_Remain - poBitMap->m_iMax_Buffer_Size < 64)
		iSize_With_Remain += 128;

	poBitMap->m_pRGBA[0] = (unsigned char*)malloc(iSize_With_Remain);
	if (!poBitMap->m_pRGBA[0])
	{
		printf("Fail to malloc in Init_Image\n");
		*poBitMap = { };
		return;
	}

	for (i = 1; i < poBitMap->m_iChannel_Count; i++)
		poBitMap->m_pRGBA[i] = poBitMap->m_pRGBA[i - 1] + iSize;
	for (; i < 4; i++)
		poBitMap->m_pRGBA[i] = NULL;

	return;
}
void Free_Image(Mem_Mgr* poMem_Mgr, Image oImage)
{
	Free(poMem_Mgr, oImage.m_pChannel[0]);
	return;
}
void Free_Image(Image* poImage)
{
	if (poImage && poImage->m_pChannel[0])
	{
		if (poImage->m_iMem_Src == Mem_Src::CPU)
		{
			//free(poImage->m_pChannel[0]);
			Free(poImage->m_pChannel[0]);
			*poImage = {};	//��֪�ĸ���׼��c��Ҫ������д����������Warning
		}else
			printf("Buffer of this image is in GPU\n");
	}/*else
		printf("Ptr is NULL\n");*/
}
void Free_Image(BitMap_8Bit* poBitMap)
{//�ͷ�BMPͼ
	if (poBitMap && poBitMap->m_pRGBA[0])
	{
		if (poBitMap->m_iMem_Src == Mem_Src::CPU)
		{
			free(poBitMap->m_pRGBA[0]);
			*poBitMap = {  };	//��֪�ĸ���׼��c��Ҫ������д����������Warning
		}else
			printf("Buffer of this image is in GPU\n");
	}else
		printf("Ptr is NULL\n");
}
void RGB_2_YUV(BitMap_8Bit oBMP, YUV_Image_420_8Bit oYUV)
{//����ת������Ϊ��׼����
	unsigned char* pR, * pG, * pB, * pY, * pU, * pV;
	int R[4], G[4], B[4], Y[4], U[4], V[4];
	int y, x, iSource_Height, iSource_Width, iSource_Width_x2, iDest_Width, iDest_Width_x2, iDest_Height, iDest_Width_Half;
	pR = oBMP.m_pRGBA[0], pG = oBMP.m_pRGBA[1], pB = oBMP.m_pRGBA[2];
	pY = oYUV.m_pYUV[0], pU = oYUV.m_pYUV[1], pV = oYUV.m_pYUV[2];
	iDest_Width = oYUV.m_iWidth;
	iDest_Width_x2 = iDest_Width << 1;
	iDest_Width_Half = (iDest_Width + 1) >> 1;
	iDest_Height = oYUV.m_iHeight;
	iSource_Width = oBMP.m_iWidth;
	iSource_Height = oBMP.m_iHeight;
	iSource_Width_x2 = iSource_Width << 1;

	for (y = 0; y < iDest_Height; y += 2, pR += iSource_Width_x2, pG += iSource_Width_x2, pB += iSource_Width_x2, pY += iDest_Width_x2, pU += iDest_Width_Half, pV += iDest_Width_Half)
	{
		for (x = 0; x < iDest_Width; x += 2)
		{
			if (x < iSource_Width - 1)
			{
				R[0] = pR[x], R[1] = pR[x + 1];
				G[0] = pG[x], G[1] = pG[x + 1];
				B[0] = pB[x], B[1] = pB[x + 1];
				if (y < iSource_Height - 1)
				{
					R[2] = pR[x + iSource_Width], R[3] = pR[x + iSource_Width + 1];
					G[2] = pG[x + iSource_Width], G[3] = pG[x + iSource_Width + 1];
					B[2] = pB[x + iSource_Width], B[3] = pB[x + iSource_Width + 1];
				}
				else
				{
					R[2] = pR[x], R[3] = pR[x + 1];
					G[2] = pG[x], G[3] = pG[x + 1];
					B[2] = pB[x], B[3] = pB[x + 1];
				}
			}
			else
			{
				R[0] = pR[x]; R[1] = pR[x];
				G[0] = pG[x]; G[1] = pG[x];
				B[0] = pB[x]; B[1] = pB[x];
				if (y < iSource_Height - 1)
				{
					R[2] = pR[x + iSource_Width], R[3] = pR[x + iSource_Width];
					G[2] = pG[x + iSource_Width], G[3] = pG[x + iSource_Width];
					B[2] = pB[x + iSource_Width], B[3] = pB[x + iSource_Width];
				}
				else
				{
					R[2] = pR[x], R[3] = pR[x];
					G[2] = pG[x], G[3] = pG[x];
					B[2] = pB[x], B[3] = pB[x];
				}

			}
			_RGB_2_YUV(R[0], G[0], B[0], Y[0], U[0], V[0]);
			_RGB_2_YUV(R[1], G[1], B[1], Y[1], U[1], V[1]);
			_RGB_2_YUV(R[2], G[2], B[2], Y[2], U[2], V[2]);
			_RGB_2_YUV(R[3], G[3], B[3], Y[3], U[3], V[3]);

			if (x >= iDest_Width - 1)
			{	//���߽��ˣ�������ȣ���д��ȥ��
			}
			else
			{
				pY[x + 1] = Y[1];
				if (y < iDest_Height - 1)
					pY[x + iDest_Width + 1] = Y[3];
			}

			pY[x] = Y[0];
			if (y < iDest_Height - 1)
				pY[x + iDest_Width] = Y[2];

			U[0] = (U[0] + U[1] + U[2] + U[3] + 2) >> 2;
			pU[x >> 1] = Clip(U[0]);
			V[0] = (V[0] + V[1] + V[2] + V[3] + 2) >> 2;
			pV[x >> 1] = Clip(V[0]);
		}
	}
	return;
}
void YUV_2_RGB(YUV_Image_444_8Bit oYUV, BitMap_8Bit oBMP)
{
	int iSize = oBMP.m_iWidth * oBMP.m_iHeight;
	unsigned char* pR, * pG, * pB, * pY, * pU, * pV;
	int i, Y, U, V;
	pR = oBMP.m_pRGBA[0], pG = oBMP.m_pRGBA[1], pB = oBMP.m_pRGBA[2];
	pY = oYUV.m_pYUV[0], pU = oYUV.m_pYUV[1], pV = oYUV.m_pYUV[2];
	
	for (i = 0; i < iSize; i++)
	{
		U = *pU;
		V = *pV;
		Y = *pY;
		_YUV_2_RGB(Y, U, V, *pR, *pG, *pB);
		pR++, pG++, pB++, pY++, pU++, pV++;
	}	
	return;
}
int bSave_Image(const char* pcFile, YUV_Image_420_8Bit oImage, int bSave_Append)
{
	int iResult, bRet = 1, iSize = oImage.m_iWidth * oImage.m_iHeight + ((oImage.m_iWidth + 1) >> 1) * ((oImage.m_iHeight + 1) >> 1) * 2;
	FILE* pFile;
	if (bSave_Append)
		pFile = fopen(pcFile, "ab");
	else
		pFile = fopen(pcFile, "wb");

	if (!pFile)
	{
		printf("Fail to save image\n");
		return 0;
	}
	iResult = (int)fwrite(oImage.m_pYUV[0], 1, iSize, pFile);
	if (iResult != iSize)
	{
		printf("Fail to save image\n");
		bRet = 0;
	}

	fclose(pFile);
	return bRet;
}
int bSave_Comp(const char* pcFile, Image oImage, int iComp)
{//����һ��
	Image oImage_1;
	//�൱��Attach
	Attach_Buffer(&oImage_1, oImage.m_pChannel[iComp], oImage.m_iWidth,
		oImage.m_iHeight, 1, Image::IMAGE_TYPE_BMP);
	return bSave_Image(pcFile, oImage_1);
}
int bSave_Image(const char* pcFile, Image oImage, int iFormat)
{//����̣�������������
	int bRet=1;
	int iTo_Save;

	if (bStricmp((char*)pcFile + strlen(pcFile) - 3, (char*)"bmp") == 0)
		iTo_Save = Image::IMAGE_TYPE_BMP;
	else if (bStricmp((char*)pcFile + strlen(pcFile) - 3, (char*)"yuv") == 0)
		iTo_Save = Image::IMAGE_TYPE_YUV_444;

	if (iFormat == -1)
	{//����ԭ��ʽ
		if (iTo_Save == Image::IMAGE_TYPE_BMP)
		{
			if (oImage.m_iImage_Type == Image::IMAGE_TYPE_BMP || oImage.m_iChannel_Count == 1)
				return bSave_Image(pcFile, oImage.m_oBitMap);
			else
			{
				BitMap_8Bit oNew;
				Init_Image(&oNew, oImage.m_iWidth, oImage.m_iHeight, oImage.m_iBit_Count);
				YUV_2_RGB(oImage.m_oYUV, oNew);
				bRet = bSave_Image(pcFile, oNew);
				Free_Image(&oNew);
				return bRet;
			}
		}
		else
		{
			/*if (oImage.m_iImage_Type == Image::IMAGE_TYPE_YUV_444)
				return bSave_Image_GPU(pcFile, &oImage.m_oYUV);
			else*/
			{
				YUV_Image_420_8Bit oNew;
				Init_Image(&oNew, oImage.m_iWidth, oImage.m_iHeight);
				RGB_2_YUV(oImage.m_oBitMap, oNew);
				return bSave_Image(pcFile, oNew,0);
			}
		}
	}
	else if (iFormat == (int)Image_Format::RGB)
	{
		if (oImage.m_iImage_Type == Image::IMAGE_TYPE_BMP)
			return bSave_Image(pcFile, oImage.m_oBitMap);

		BitMap_8Bit oNew;
		Init_Image(&oNew, oImage.m_iWidth, oImage.m_iHeight, oImage.m_iBit_Count);
		YUV_2_RGB(oImage.m_oYUV, oNew);
		bRet = bSave_Image(pcFile, oNew);
		Free_Image(&oNew);
		return bRet;
	}
	/*else if (iFormat == (int)Image_Format::YUV_444 || iFormat == (int)Image_Format::YUVA_444)
	{
		if (oImage.m_iImage_Type == Image::IMAGE_TYPE_YUV_444)
			return bSave_Image_GPU(pcFile, &oImage.m_oYUV);

		printf("No implemented yet\n");
		return 0;
	}*/
	return 0;
}
int bSave_Image(const char* pcFile, float* pImage, int iWidth, int iHeight)
{//���������Ҷ�ֵ���ͼ
	Image oImage;
	Init_Image(&oImage, iWidth, iHeight, Image::IMAGE_TYPE_BMP, 8);
	Float_2_Image(pImage, oImage);
	return bSave_Image(pcFile, oImage);
}
void Draw_Line(Image oImage, int x0, int y0, int x1, int y1, int R, int G, int B)
{
	int d, d_upper, d_lower;
	int a, b, x, y, iPos;
	int iTemp, iWidth = oImage.m_iWidth, iHeight = oImage.m_iHeight;
	unsigned char* pR = oImage.m_pChannel[0],
		* pG = oImage.m_pChannel[1], * pB = oImage.m_pChannel[2];

	if (abs(y1 - y0) <= abs(x1 - x0))
	{//45������
		if (x0 > x1)
		{//��x����С���󽻻�ͷβ
			iTemp = x0, x0 = x1, x1 = iTemp;
			iTemp = y0, y0 = y1; y1 = iTemp;
		}
		a = y0 - y1;
		b = x1 - x0;	//���㷨����Ҫc�����Ƶ�������Լȥ

		if (y1 >= y0)
		{
			d = a * 2 + b;//�����ȷҪ�����㣬��ѭ���ǲſ�������������Ϊ�������жϣ�����Ϊ�ڶ����ж�
			d_upper = (a + b) * 2;
		}
		else
		{
			d = 2 * a - b;
			d_upper = (a - b) * 2;
		}
		d_lower = a * 2;

		//����0��
		if (y0 >= 0 && y0 < iHeight && x0 >= 0 && x0 < iWidth)
		{
			iPos = y0 * iWidth + x0;
			pR[iPos] = R;
			if (oImage.m_pChannel[1])
				pG[iPos] = G;
			if (oImage.m_pChannel[2])
				pB[iPos] = B;
		}
		if (y1 >= y0)
		{//0-45��
			for (x = x0, y = y0; x < x1;)
			{
				if (d > 0)//M����ֱ���Ϸ���ȡ�����
					d += d_lower;	//y���䣬����һ��һ��
				else
				{
					d += d_upper;
					y++;
				}
				x++;
				if (y >= 0 && y < iHeight && x >= 0 && x < iWidth)
				{
					iPos = y * iWidth + x;
					pR[iPos] = R;
					if (oImage.m_pChannel[1])
						pG[iPos] = G;
					if (oImage.m_pChannel[2])
						pB[iPos] = B;
				}
			}
		}
		else
		{//-45�ȵ�0
			for (x = x0, y = y0; x < x1;)
			{
				if (d < 0)//M����ֱ���Ϸ���ȡ�����
					d += d_lower;	//y���䣬����һ��һ��
				else
				{
					d += d_upper;
					y--;
				}
				x++;
				if (y >= 0 && y < iHeight && x >= 0 && x < iWidth)
				{
					iPos = y * iWidth + x;
					pR[iPos] = R;
					if (oImage.m_pChannel[1])
						pG[iPos] = G;
					if(oImage.m_pChannel[2])
						pB[iPos] = B;
				}
			}
		}
	}
	else
	{//��yΪ�Ա��������ϼ�һ���� ����45��
		if (y0 > y1)
		{//��x����С���󽻻�ͷβ
			iTemp = x0, x0 = x1, x1 = iTemp;
			iTemp = y0, y0 = y1; y1 = iTemp;
		}
		a = x0 - x1;
		b = y1 - y0;	//���㷨����Ҫc�����Ƶ�������Լȥ
		if (x1 >= x0)
		{
			d = 2 * b + a;	//�����ȷҪ�����㣬��ѭ���ǲſ�������������Ϊ�������жϣ�����Ϊ�ڶ����ж�
			d_upper = (a + b) * 2;
		}
		else
		{
			d = 2 * a - b;
			d_upper = (a - b) * 2;
		}
		d_lower = a * 2;
		//����0��
		if (y0 >= 0 && y0 < iHeight && x0 >= 0 && x0 < iWidth)
		{
			iPos = y0 * iWidth + x0;
			pR[iPos] = R;
			if (oImage.m_pChannel[1])
				pG[iPos] = G;
			if (oImage.m_pChannel[2])
				pB[iPos] = B;
		}

		if (x1 >= x0)
		{//45��-90��
			for (x = x0, y = y0; y <= y1;)
			{
				if (d > 0)//M����ֱ���Ϸ���ȡ�����
					d += d_lower;	//y���䣬����һ��һ��
				else
				{
					d += d_upper;
					x++;
				}
				y++;
				if (x < 0 || y < 0)
					continue;
				else if (y < iHeight && x < iWidth)
				{
					iPos = y * iWidth + x;
					pR[iPos] = R;
					if (oImage.m_pChannel[1])
						pG[iPos] = G;
					if (oImage.m_pChannel[2])
						pB[iPos] = B;
				}
				else
					break;
			}
		}
		else
		{//90�ȵ�135��
			for (x = x0, y = y0; y < y1;)
			{
				if (d < 0)//M����ֱ���Ϸ���ȡ�����
					d += d_lower;	//y���䣬����һ��һ��
				else
				{
					d += d_upper;
					x--;
				}
				y++;
				if (y > 0 && y < iHeight && x >= 0 && x < iWidth)
				{
					iPos = y * iWidth + x;
					pR[iPos] = R;
					if (oImage.m_pChannel[1])
						pG[iPos] = G;
					if (oImage.m_pChannel[2])
						pB[iPos] = B;
				}
				//else
					//break;
			}
		}
	}
	return;
}

void Draw_Arc(Image oImage, int r, int iCenter_x, int iCenter_y, float fAngle_Start, float fAngle_End, int R, int G, int B)
{//��Բ������Ȼ���Ապ��˳�ΪԲ
	int Y, U, V;
	_RGB_2_YUV(R, G, B, Y, U, V);
	float fTheta, delta_theta = (float)(1.f * 0.99f / (r * 2.f));
	int y, x, iPos;
	float fy, fx;

	if (fAngle_End > fAngle_Start)
	{
		for (fTheta = fAngle_Start; fTheta <= fAngle_End; fTheta += delta_theta)
		{
			fy = (float)r * (float)sin(fTheta);
			fx = (float)r * (float)cos(fTheta);
			y = (fy >= 0 ? (int)(fy + 0.5) : (int)(fy - 0.5)) + iCenter_y;
			x = (fx >= 0 ? (int)(fx + 0.5) : (int)(fx - 0.5)) + iCenter_x;
			iPos = y * oImage.m_iWidth + x;
			if (y >= 0 && y < oImage.m_iHeight && x >= 0 && x < oImage.m_iWidth)
			{
				if(oImage.m_pChannel[0])
					oImage.m_pChannel[0][iPos] = R;
				if (oImage.m_pChannel[1])
					oImage.m_pChannel[1][iPos] = G;
				if (oImage.m_pChannel[2])
					oImage.m_pChannel[2][iPos] = B;
			}
		}
	}
	else
	{
		for (fTheta = fAngle_Start; fTheta >= fAngle_End; fTheta -= delta_theta)
		{
			fy = (float)r * (float)sin(fTheta);
			fx = (float)r * (float)cos(fTheta);
			y = (fy >= 0 ? (int)(fy + 0.5) : (int)(fy - 0.5)) + iCenter_y;
			x = (fx >= 0 ? (int)(fx + 0.5) : (int)(fx - 0.5)) + iCenter_x;
			iPos = y * oImage.m_iWidth + x;
			if (y >= 0 && y < oImage.m_iHeight && x >= 0 && x < oImage.m_iWidth)
			{
				oImage.m_pChannel[0][iPos] = R;
				oImage.m_pChannel[1][iPos] = G;
				oImage.m_pChannel[2][iPos] = B;
			}
		}
	}
	return;
}

void _Fill_Region(Image* poImage, int iSeed_x, int iSeed_y, int iBack_Color, int iFill_Color)
{
	int iPos = iSeed_y * poImage->m_iWidth,
		iWidth = poImage->m_iWidth;
	unsigned char* pCur_Next,
		* pCur_R = poImage->m_pChannel[0] ? &poImage->m_pChannel[0][iPos] : NULL,
		* pCur_G = poImage->m_pChannel[1] ? &poImage->m_pChannel[1][iPos] : NULL,
		* pCur_B = poImage->m_pChannel[2] ? &poImage->m_pChannel[2][iPos] : NULL;

	int x, x_Left, x_Right;	//����ߣ����ұ߿��������

	if (iSeed_y < 0)
		return;
	//Ѱ������
	for (x = iSeed_x; x >= 0 && pCur_R[x] == iBack_Color;)
	{
		if (pCur_R)
			pCur_R[x] = iFill_Color;
		if(pCur_G)
			pCur_G[x] = iFill_Color;
		if(pCur_G)
			pCur_B[x] = iFill_Color;
		x--;
	}

	x_Left = x + 1;
	x_Right = iSeed_x;

	//���Ͽ�
	if (iSeed_y > 0)
	{
		pCur_Next = pCur_R - iWidth;
		for (x = x_Left; x < x_Right;)
		{
			if (pCur_Next[x] == iBack_Color)
			{
				//����Ѱ�����ҵ�
				while (x < iWidth && pCur_Next[x] == iBack_Color)
					x++;
				_Fill_Region(poImage, x - 1, iSeed_y - 1, iBack_Color, iFill_Color);
			}
			else
				x++;
		}
	}

	//������
	if (iSeed_y < iWidth - 1)
	{
		pCur_Next = pCur_R + iWidth;
		for (x = x_Left; x < x_Right;)
		{
			if (pCur_Next[x] == iBack_Color)
			{
				//����Ѱ�����ҵ�
				while (x < iWidth && pCur_Next[x] == iBack_Color)
					x++;
				_Fill_Region(poImage, x - 1, iSeed_y + 1, iBack_Color, iFill_Color);
			}
			else
				x++;
		}
	}

	return;
}
void Fill_Region(Image oImage, int iSeed_x, int iSeed_y, int iBack_Color)
{//Ϊ�˼򻯣�ֻ��һ�ֱ���ɫ��ֻ��Y����
	unsigned char* pCur = &oImage.m_pChannel[0][iSeed_y * oImage.m_iWidth];
	int x;

	if (iSeed_y >= 0)
	{
		if (pCur[iSeed_x] != iBack_Color)
		{
			printf("The Seed is not in back color:%d\n", oImage.m_pChannel[0][iSeed_y * oImage.m_iWidth + iSeed_x]);
			return;
		}
		//Ѱ������
		for (x = iSeed_x; x < oImage.m_iWidth && pCur[x] == iBack_Color; x++);
		iSeed_x = x - 1;
	}
	_Fill_Region(&oImage, iSeed_x, iSeed_y, iBack_Color, 255);
}
void Draw_Point(Image oImage, int x, int y, int r, int R, int G, int B)
{//ƽ���ϻ�һ����
	if (x < 0 || x >= oImage.m_iWidth || y < 0 || y >= oImage.m_iHeight)
		return;

	Draw_Arc(oImage, r, x, y, 0, (float)(2 * PI), R, G, B);

	//Fill_Region(oImage, x, y);
	if(oImage.m_pChannel[0])
		oImage.m_pChannel[0][y * oImage.m_iWidth + x] = R;
	if (oImage.m_pChannel[1])
		oImage.m_pChannel[1][y * oImage.m_iWidth + x] = G;
	if (oImage.m_pChannel[2])
		oImage.m_pChannel[2][y * oImage.m_iWidth + x] = B;
	return;
}
void Gen_Gauss_Filter(int r, float** ppFilter)
{//����r�Ĵ�С����sigma
	float fValue, fSum = 0;
	int i, j, d = r + r + 1;
	float fSigma = 0.3f * ((d - 1.f) * 0.5f - 1.f) + 0.8f;
	float* pFilter = (float*)malloc(d * sizeof(float));

	for (i = 0; i < r; i++)
	{
		fValue = (float)(i - r) / fSigma;
		pFilter[i] = (float)exp(-0.5f * fValue * fValue);
		fSum += pFilter[i];
	}
	fValue = (float)(i - r) / fSigma;
	pFilter[i] = (float)exp(-0.5f * fValue * fValue);
	fSum = fSum * 2 + pFilter[i];

	//Temp Code
	//for (i = 0; i < r + 1; i++)
		//printf("%f\n", pFilter[i]);

	fValue = 1.f / fSum;
	j = d - 1;
	for (i = 0; i < r; i++, j--)
		pFilter[j] = (pFilter[i] *= fValue);
	pFilter[r] *= fValue;
	*ppFilter = pFilter;

	return;
}

void Gen_Gauss_Filter(int r, float fSigma,float **ppFilter)
{//fSigma��Ӧ��̬�ֲ��е�sigma, �������ܹ����һ�������ĸ�˹�ˣ����Ϊ1
//�˴���Ҫ��һ���㣬fSigma��ȡֵ�����˺˴�С
	int i,j, d = r + r + 1;
	float* pFilter = (float*)malloc(d * sizeof(float));
	float fValue,fSum=0;

	////�о����������̫���£���һ��͹���
	//for (i = 0; i < d; i++)
	//{
	//	fValue = (float)(i - r) / fSigma;
	//	pFilter[i] = exp(-0.5f * fValue * fValue);
	//	fSum += pFilter[i];
	//}
	//fValue = 1.f / fSum;
	//for (i = 0; i < d; i++)
	//	pFilter[i] *= fValue;

	for (i = 0; i < r; i++)
	{
		fValue = (float)(i - r) / fSigma;
		pFilter[i] = (float)exp(-0.5f * fValue * fValue);
		fSum += pFilter[i];
	}

	fValue = (float)(i - r) / fSigma;
	pFilter[i] =(float)exp(-0.5f * fValue * fValue);
	fSum = fSum*2+ pFilter[i];

	fValue = 1.f / fSum;
	j = d - 1;
	for (i = 0; i < r; i++,j--)
		pFilter[j] = (pFilter[i] *= fValue);
	pFilter[r] *= fValue;
	*ppFilter = pFilter;
	
	return;
}

void Transpose_AVX512(float* pSource, int iWidth, int iHeight, float* pDest)
{//����ת��
	int y, x,i,iPos, iEnd;
	int iWidth_Align_16 = (iWidth >> 4) << 4,
		iHeight_Align_16 = (iHeight >> 4) << 4;
	int x_Remain_Mask = (1 << (iHeight - iHeight_Align_16)) - 1;
		//y_Remain_Mask = (1 << (iWidth - iWidth_Align_16)) - 1;

	union {
		__m512i vIndex;
		int _vIndex[16];
	};

	float Buffer[16][16];

	for (x = 0; x < 16; x++)
		_vIndex[x] = x * 16;
		//vIndex.m512i_i32[x] = x * 16;	//����Intel ����������
		
	for (y = 0; y < iHeight_Align_16; y += 16)
	{//���
		for (x = 0; x < iWidth_Align_16; x += 16)
		{
			iPos = y * iWidth + x;
			iEnd = y<iHeight_Align_16?16:iHeight - y;
			for (i = 0; i < iEnd; i++, iPos += iWidth)
				_mm512_i32scatter_ps(&Buffer[0][i], vIndex, _mm512_loadu_ps(&pSource[iPos]), 4);

			//Ȼ��д��
			iPos = x * iHeight + y;
			iEnd = x<iWidth_Align_16?16:iWidth - x;
			for (i = 0; i < iEnd; i++, iPos += iHeight)
				_mm512_storeu_ps(&pDest[iPos], *(__m512*)Buffer[i]);
				//*(__m512*)&pDest[iPos] = *(__m512*)Buffer[i];
			//*(__m512*)&pDest[iPos] = *(__m512*)&Buffer[i][1];
		}
		if (x < iWidth)
		{//��
			iPos = y * iWidth + x;
			//_mm512_i32scatter_ps �Ǹ���vIndex���λ��������д�룬ת�þ����������
			//�������һ��16�����ݣ�д��Buffer��һ��
			for (i = 0; i < 16; i++, iPos += iWidth)
				_mm512_i32scatter_ps(&Buffer[0][i], vIndex, _mm512_loadu_ps(&pSource[iPos]), 4);

			//Ȼ��д��
			iPos = x * iHeight + y;
			iEnd = iWidth - x;
			for (i = 0; i < iEnd; i++, iPos += iHeight)
				*(__m512*)&pDest[iPos] = *(__m512*)Buffer[i];
		}
	}

	if (y < iHeight)
	{
		for (x = 0; x < iWidth_Align_16; x += 16)
		{
			iPos = y * iWidth + x;
			iEnd = y < iHeight_Align_16 ? 16 : iHeight - y;
			for (i = 0; i < iEnd; i++, iPos += iWidth)
				_mm512_i32scatter_ps(&Buffer[0][i], vIndex, _mm512_loadu_ps(&pSource[iPos]), 4);

			//Ȼ��д��
			iPos = x * iHeight + y;
			iEnd = x < iWidth_Align_16 ? 16 : iWidth - x;
			
			for (i = 0; i < iEnd; i++, iPos += iHeight)
				_mm512_mask_compressstoreu_ps(&pDest[iPos], x_Remain_Mask, *(__m512*)Buffer[i]);
		}
	}
	return;
}

void _Gause_Filter_AVX512(float* pSource, int iWidth, int iHeight, float* pDest, int r, float* pFilter)
{
	__m512 Pixel_0_16, Pixel_1_16, Sum_16;
	int iWidth_Align_16 = (iWidth >> 4) << 4;		//���ǰ���16������������һЩ
		//iHeight_Align_16 = (iHeight >> 4) << 4;
	//�����Ǹ����룬ֻ�����ֵ���˼
	int x_Remain_Mask = (1 << (iWidth - iWidth_Align_16)) - 1;
		//y_Remain_Mask = (1 << (iHeight - iHeight_Align_16)) - 1;
	int x, y, y1, i;	// , d = r + r + 1;
	//int d_minus_1 = d - 1;
	float* pFilter_1 = &pFilter[r];

	//���з�����һά��˹�˲�
	for (y = 0; y < iHeight; y++)
	{//�����ϣ����ش���Ĺ��������к���
		for (x = 0; x < iWidth; x += 16)
		{
			Sum_16 = _mm512_set1_ps(0.f);	//��ʾ��Sum_16��16���������
			for (i = 1; i <= r; i++)
			{//���Լ�Ϊ���ģ�����������������ɢ��ȥ
				y1 = y - i;	//y1��ʾ�Լ������i�У�Pixel_0_16��ʾ���ڴ�����������16�����ݽ��Ĵ���
				Pixel_0_16 = _mm512_loadu_ps(&pSource[x + iWidth * (y1 < 0 ? 0 : y1)]);
				y1 = y + i; //y1��ʾ�����i�У�Pixel_1_16��ʾ���ڴ�����������16�����ݽ��Ĵ���
				Pixel_1_16 = _mm512_loadu_ps(&pSource[x + iWidth * (y1 >= iHeight ? iHeight - 1 : y1)]);
				//Sum = (Pixel_0 + Pixel_1) * Filter[i]���˴�������һ�����ʣ����Լ�Ϊ���ģ���˹�����ſ�
				//�����i�������������i�����ݶ�Ӧͬһ��Ȩֵ�������ɸ�˹�˵ĶԳ��Ծ����ġ�
				Sum_16 = _mm512_add_ps(Sum_16, _mm512_mul_ps(_mm512_add_ps(Pixel_0_16, Pixel_1_16), _mm512_set1_ps(pFilter_1[i])));
			}
			//����������Լ�Ҳ����
			Pixel_0_16 = _mm512_loadu_ps(&pSource[x + iWidth * y]);
			Sum_16 = _mm512_add_ps(Sum_16, _mm512_mul_ps(Pixel_0_16, _mm512_set1_ps(*pFilter_1)));

			//�㵽���16����Ȩ����Ѿ��õ������Ŵ�ص�Ŀ���ڴ�λ����
			if (x + 16 > iWidth)	//�����ǲ���16������£���
				_mm512_mask_compressstoreu_ps(&pDest[y * iWidth + x], x_Remain_Mask, Sum_16);
			else					//���ﹻ16�����ݣ���д���ڴ漴��
				_mm512_storeu_ps(&pDest[y*iWidth+x],Sum_16);	//����ԭͼ˳��
		}
	}
}

void Gauss_Filter_AVX512(float* pSource, int iWidth, int iHeight, float* pDest, float *pAux, int r, float* pFilter)
{//����һ����ʱ������pAux���н���
	//���������з���ķ�������������Ҫ������,�˴��ǵ�һ�Σ���˹���з�������
	_Gause_Filter_AVX512(pSource, iWidth, iHeight, pAux, r, pFilter);
	//ת��һ�£��ڴ�ķ����ٶȾ����˾�����һ������
	Transpose_AVX512(pAux, iWidth, iHeight, pDest);
	//ת�������Ժ�ͼת��Ϊh*w����״����ͬ���ķ������ڶ��Σ���˹�˻������з���
	_Gause_Filter_AVX512(pDest, iHeight, iWidth,pAux, r, pFilter);
	//��ת�ã��ָ�ԭͼ
	Transpose_AVX512(pAux, iHeight, iWidth, pDest);
	return;
}

void Gauss_Filter_Ref(float* pSource, int iWidth, int iHeight,float *pDest,int r,float *pFilter)
{//��һ�Ÿ�������ʾ��ͼ���и�˹ģ���������Ϊ��׼�������Զ�����
	float *pFilter_Mid;
	float* pTemp = (float*)malloc(iWidth * iHeight * sizeof(float));
	float fSum;
	int x, y,i,x1,y1,iPos;
	pFilter_Mid = pFilter + r;
	
	//����
	for (x = 0; x < iWidth; x++)
	{
		for (y = 0; y < iHeight; y++)
		{
			for (fSum = 0, i = -r; i <= r; i++)
			{
				y1 = y + i;
				y1 = Clip3(0, iHeight - 1, y1);
				fSum += pSource[y1 * iWidth + x] * pFilter_Mid[i];
				//printf("pixel:%f weight:%f sum:%f\n", pSource[y1 * iWidth + x] * pFilter_Mid[i], pFilter_Mid[i], fSum);
			}
			pTemp[y * iWidth + x] = fSum;
		}
	}
	Disp(pTemp, iHeight, iWidth, "Mid");

	//����
	for (y = 0; y < iHeight; y++)
	{
		iPos = y * iWidth;
		for (x = 0; x < iWidth; x++)
		{
			for (fSum = 0, i = -r; i <= r; i++)
			{
				x1 = x + i;
				fSum += pTemp[iPos + Clip3(0, iWidth - 1, x1)] * pFilter_Mid[i];
				//printf("pixel:%f weight:%f sum:%f\n", pTemp[iPos + Clip3(0, iWidth - 1, x1)], pFilter_Mid[i], fSum);
			}
			pDest[iPos + x] = fSum;
		}
	}
	free(pTemp);
	return;
}
void Float_2_Image(float* pImage, Image oImage)
{
	int i,j,iSize = oImage.m_iWidth * oImage.m_iHeight;
	int iValue;
	for (i = 0; i < iSize; i++)
	{
		iValue= (int)(pImage[i] * 255.f);
		for(j=0;j<oImage.m_iChannel_Count;j++)
			oImage.m_pChannel[j][i] = iValue;
	}
}
void RGB_2_Gray(Image oImage, float* pImage)
{//RGB ת���Ҷ�ͼ��grey=0.299r+0.587g+0.114b, ��/255.f, �Բ���
//��һ��ת����ʽ�� grey=0.2126r +0.7152g +0.0722b
	int i,iSize = oImage.m_iWidth * oImage.m_iHeight;
	if (oImage.m_iChannel_Count >= 3)
		for (i = 0; i < iSize; i++)
			pImage[i] = ((float)round(0.2126f * oImage.m_pChannel[0][i] + 0.7152 * oImage.m_pChannel[1][i] + 0.0722f * oImage.m_pChannel[2][i])) / 255.f;
	else if (oImage.m_iChannel_Count == 1)
		for (i = 0; i < iSize; i++)
			pImage[i] = oImage.m_pChannel[0][i]/255.f;

	
		//pImage[i] = (0.2126f / 255.f) * oImage.m_pChannel[0][i] + (0.7152 / 255.f) * oImage.m_pChannel[1][i] + (0.0722f / 255.f) * oImage.m_pChannel[2][i];
		//pImage[i] = (0.299f/255.f) * oImage.m_pChannel[0][i] + (0.587 / 255.f)* oImage.m_pChannel[1][i] + (0.114f / 255.f)* oImage.m_pChannel[2][i];
	return;
}

int iGet_Image_Size(int iWidth, int iHeight, int iChannel_Count)
{//����һ�����������ͼ�����Ҫ�ö����ֽ�
	int iSize = iWidth * iHeight;
	int iMax_Buffer_Size = iSize * iChannel_Count;
	int iSize_With_Remain = ((iMax_Buffer_Size + 127) >> 7) << 7;
	if (iSize_With_Remain - iMax_Buffer_Size < 64)
		iSize_With_Remain += 128;
	return ALIGN_SIZE_1024(iSize_With_Remain);	//�������ڴ�������Ŀ�
}
void Attach_Buffer(Image* poImage, unsigned char* pBuffer, int iWidth, int iHeight, int iChannel_Count, int iImage_Type)
{
	if (!poImage || iWidth * iHeight == 0)
	{
		printf("Invalid parameter in Attach_Buffer_GPU\n");
		if (poImage)
			*poImage = { {0 } };
		return;
	}

	poImage->m_iWidth = iWidth;
	poImage->m_iHeight = iHeight;
	poImage->m_iBit_Count = iChannel_Count << 3;
	poImage->m_iChannel_Count = iChannel_Count;
	poImage->m_iMax_Buffer_Size = iWidth * iHeight * poImage->m_iChannel_Count;
	poImage->m_iMem_Src = Mem_Src::CPU;
	poImage->m_bIs_Attached = 1;
	poImage->m_iImage_Type = iImage_Type;

	int i;
	for (i = 0; i < poImage->m_iChannel_Count; i++)
		poImage->m_pChannel[i] = pBuffer + i * iWidth * iHeight;
	for (; i < 4; i++)
		poImage->m_pChannel[i] = NULL;
}

void Place_Image(Image oTile, Image oScreen, int x, int y)
{//x,y�����Ǹ�����Ҳ���Գ�Խ��Ļ����
	int i,w,w_Align_8, h,
		iTile_Cur, iScreen_Cur,
		x_Screen_Start, y_Screen_Start,x_Screen_End,y_Screen_End,
		x_Tile_Start,y_Tile_Start,x_Tile_End,y_Tile_End;	//����oScreen��ʼλ��
	unsigned char* pTile_Cur, * pScreen_Cur;
	if (x < 0)
		x_Screen_Start = 0,x_Tile_Start = -x;
	else
		x_Tile_Start = 0,	x_Screen_Start = x;
	w = oTile.m_iWidth - x_Tile_Start;
	if (x_Screen_Start + w > oScreen.m_iWidth)
		w -= x_Screen_Start + w - oScreen.m_iWidth;
	x_Tile_End = x_Tile_Start + w;
	x_Screen_End = x_Screen_Start + w;	
	w_Align_8 = (w >> 3) << 3;

	if (y < 0)
		y_Screen_Start = 0, y_Tile_Start = -y;
	else
		y_Tile_Start = 0,	y_Screen_Start = y;
	h = oTile.m_iHeight - y_Tile_Start;
	if (y_Screen_Start + h > oScreen.m_iHeight)
		h -= y_Screen_Start + h - oScreen.m_iHeight;
	y_Tile_End =y_Tile_Start+ h, y_Screen_End =y_Screen_Start+ h;
	
	
	iTile_Cur = y_Tile_Start * oTile.m_iWidth + x_Tile_Start;
	iScreen_Cur = y_Screen_Start * oScreen.m_iWidth + x_Screen_Start;
	for (y = y_Tile_Start; y < y_Tile_End; y++,iTile_Cur+=oTile.m_iWidth,iScreen_Cur+=oScreen.m_iWidth)
	{
		for (i = 0; i < 3; i++)
		{
			pTile_Cur = &oTile.m_pChannel[i][iTile_Cur];
			pScreen_Cur = &oScreen.m_pChannel[i][iScreen_Cur];
			for (x = 0; x < w_Align_8; x += 8, pTile_Cur+=8,pScreen_Cur +=8)
				*(unsigned long long*)pScreen_Cur = *(unsigned long long*)pTile_Cur;
			for (; x < w; x++)
				*pScreen_Cur++ = *pTile_Cur++;
		}
	}
	return;
}
void Concat_Image(Image oA, Image oB, Image* poC, int iFlag)
{//��ͼƴһͼC = A + B, iFlag =0, ˮƽ�� iFlat=1: ��ֱ
	Image oC = *poC;
	int i, iBit_Count;
	iBit_Count = Min(oA.m_iBit_Count, oB.m_iBit_Count);
	iBit_Count = Min(iBit_Count, 24);

	if (!oC.m_pBuffer)
		if (iFlag == 0)   //ˮƽ
			Init_Image(&oC, oA.m_iWidth + oB.m_iWidth, Max(oA.m_iHeight, oB.m_iHeight), Image::IMAGE_TYPE_BMP, iBit_Count);
		else
			Init_Image(&oC, Max(oA.m_iWidth, oB.m_iWidth), oA.m_iHeight + oB.m_iHeight, Image::IMAGE_TYPE_BMP, iBit_Count);

	if (iFlag == 0)
	{//ˮƽƴ��
		for (i = 0; i < iBit_Count / 8; i++)
		{
			int y, x, iDest_Pos, iSource_Pos;
			for (iSource_Pos = 0, y = 0; y < oA.m_iHeight; y++)
			{
				iDest_Pos = y * oC.m_iWidth;
				for (x = 0; x < oA.m_iWidth; x++, iDest_Pos++, iSource_Pos++)
					oC.m_pChannel[i][iDest_Pos] = oA.m_pChannel[i][iSource_Pos];
			}
			for (iSource_Pos = 0, y = 0; y < oB.m_iHeight; y++)
			{
				iDest_Pos = y * oC.m_iWidth + oA.m_iWidth;
				for (x = 0; x < oB.m_iWidth; x++, iDest_Pos++, iSource_Pos++)
					oC.m_pChannel[i][iDest_Pos] = oB.m_pChannel[i][iSource_Pos];
			}
		}
	}
	*poC = oC;
	return;
}
template<typename _T>void Draw_Match_Point(_T Point_1[][2], _T Point_2[][2], int iCount, Image oImage)
{
	int iOrg_Image_Width = oImage.m_iWidth / 2;
	int i;
	for (i = 0; i < iCount; i++)
	{
		Draw_Point(oImage, (int)Point_1[i][0], (int)Point_1[i][1], 5, 255, 0, 0);
		Draw_Point(oImage, (int)Point_2[i][0] + iOrg_Image_Width, (int)Point_1[i][1], 5, 0, 255, 0);
		Draw_Line(oImage, (int)Point_1[i][0], (int)Point_1[i][1],
			(int)Point_2[i][0] + iOrg_Image_Width, (int)Point_1[i][1]);
	}
	return;
}
template<typename _T>void Draw_Match_Point(_T Point_1[][2], _T Point_2[][2], unsigned char Mask[], int iCount, Image oImage)
{//oImageͨ��������ͼƴ�Ӷ���
	int iOrg_Image_Width = oImage.m_iWidth / 2;
	int i;
	for (i = 0; i < iCount; i++)
	{
		if (Mask[i])
		{
			Draw_Point(oImage, (int)Point_1[i][0], (int)Point_1[i][1], 5, 255, 0, 0);
			Draw_Point(oImage, (int)Point_2[i][0] + iOrg_Image_Width, (int)Point_1[i][1], 5, 0, 255, 0);
			Draw_Line(oImage, (int)Point_1[i][0], (int)Point_1[i][1],
				(int)Point_2[i][0] + iOrg_Image_Width, (int)Point_1[i][1]);
		}
	}
	return;
}
void SB_Image()
{
	Get_Bounding_Box((unsigned short(*)[2])NULL, 0, (unsigned short(*)[2])NULL);
	
	Draw_Match_Point((double(*)[2])NULL, (double(*)[2])NULL, 0, {});
	Draw_Match_Point((float(*)[2])NULL, (float(*)[2])NULL, 0, {});

	Draw_Match_Point((double(*)[2])NULL, (double(*)[2])NULL, NULL, 0, {});
	Draw_Match_Point((float(*)[2])NULL, (float(*)[2])NULL, NULL, 0, {});

	Bi_Linear<double>(NULL,0,0,0,0,NULL);
	Bi_Linear<float>(NULL, 0, 0, 0, 0, NULL);
	Bi_Linear<unsigned char>(NULL, 0, 0, 0, 0, NULL);


	Sep_Filter_2D<double>(NULL, 0, 0, NULL, NULL, 0, NULL);
	Sep_Filter_2D<float>(NULL, 0, 0, NULL, NULL, 0, NULL);

	Sep_Filter_2D_1<unsigned char, int, int>(NULL, 0, 0, NULL, NULL, 0, NULL);
	Sep_Filter_2D_1<unsigned char, short, int>(NULL, 0, 0, NULL, NULL, 0, NULL);

	Sep_Filter_2D_1<unsigned char, short, short>(NULL, 0, 0, NULL, NULL, 0, NULL);
	Sep_Filter_2D_1<unsigned char, unsigned char,char>(NULL, 0, 0, NULL, NULL, 0, NULL);
	Sep_Filter_2D_1<unsigned char, char, char>(NULL, 0, 0, NULL, NULL, 0, NULL);
}
template<typename _T>void Get_Bounding_Box(_T Point[][2], int iPoint_Count,_T Bounding_Box[2][2])
{
	int i;
	//�˴����ָ�һ��������Ѱ�������Сֵ�ķ���
	Bounding_Box[0][0] = Bounding_Box[1][0] = Point[0][0];
	Bounding_Box[0][1] = Bounding_Box[1][1] = Point[0][1];
	for (i = 1; i < iPoint_Count; i++)
	{
		if (Point[i][0] < Bounding_Box[0][0])
			Bounding_Box[0][0] = Point[i][0];
		else if (Point[i][0] > Bounding_Box[1][0])
			Bounding_Box[1][0] = Point[i][0];

		if (Point[i][1] < Bounding_Box[0][1])
			Bounding_Box[0][1] = Point[i][1];
		else if (Point[i][1] > Bounding_Box[1][1])
			Bounding_Box[1][1] = Point[i][1];
	}
	return;
}
void Box_Filter_Ref(Image oSource,Image oDest,int r)
{//�������ü�ֵ��ֻ���������壬�ʴ���Ϊ�ο�����
	int iHeight_Minus_1 = oSource.m_iHeight - 1,
		iWidth_Minus_1 = oSource.m_iWidth - 1,
		iStencil_Size = (2*r+1) * (2*r+1);
	for (int iChannel = 0; iChannel < oSource.m_iChannel_Count; iChannel++)
	{
		unsigned char* pSource = oSource.m_pChannel[iChannel],
			* pDest = oDest.m_pChannel[iChannel];
		for (int y = 0; y < oSource.m_iHeight; y++)
		{
			for (int x = 0; x < oSource.m_iWidth; x++)
			{//�Լ���ȷ��������(x,y)λ��
				/*if (y == 100 && x == 799)
					printf("here");*/
				int iTotal = 0;
				for (int i = -r; i <= r; i++)
				{
					int y_s = y + i;
					unsigned char *pSource_1 = &pSource[Clip3(0, iHeight_Minus_1, y_s) * oSource.m_iWidth];
					for (int j = -r; j <= r; j++)
					{
						int x_s = x + j;
						x_s = Clip3(0, iWidth_Minus_1, x_s);
						iTotal += pSource_1[x_s];
						/*if (y == 100 && x == 799 && pSource_1[x_s]!=255)
							printf("here");*/
					}
				}
				pDest[y*oDest.m_iWidth + x] = iTotal / iStencil_Size;
			}
		}
	}
	return;
}

void Box_Filter(Image oSource, Image oDest, int r)
{//�������з��룬�������������伫���о�c��Χ���Ѿ����޿���
	unsigned int* pSum_Buffer = (unsigned int*)pMalloc(oSource.m_iWidth * oSource.m_iHeight*sizeof(unsigned int));
	const int iWidth_Minus_1 = oSource.m_iWidth - 1,
		iHeight_Minus_1 = oSource.m_iHeight -1,
		iStencil_Size = (2 * r + 1) * (2 * r + 1);
	float iStencil_Size_Recip = 1.f / iStencil_Size;
	//memset(pSum_Buffer, 0, oSource.m_iWidth * oSource.m_iHeight * sizeof(unsigned int));
	int bFast_Row_Process = oSource.m_iWidth >= 2 * r + 1?1 : 0,
		bFast_Col_Process = oSource.m_iHeight >= 2 * r + 1 ? 1 : 0;

	for (int iChannel = 0; iChannel < oSource.m_iChannel_Count; iChannel++)
	{
		//�ȸ���
		for (int y = 0; y < oSource.m_iHeight; y++)
		{
			int x, iLeft = -r;	// , iRight = r;
			unsigned char* pSource_Line = &oSource.m_pChannel[iChannel][y * oSource.m_iWidth];

			//�˴�Dest��Ϊת�þ���
			unsigned int* pDest_Line = &pSum_Buffer[y];
			unsigned int iTotal = r * pSource_Line[0];
			if (oSource.m_iWidth < r + 1)
			{//������������,�ݲ�����
				for (x = 0; x < oSource.m_iWidth; x++)
					iTotal += pSource_Line[x];
				//ʣ�²���Ĳ���β
				iTotal += pSource_Line[iWidth_Minus_1] * (r + 1 - oSource.m_iWidth);
			}else
			{
				for (x = 0; x <= r; x++)
					iTotal += pSource_Line[x];
			}
			pDest_Line[0] = iTotal;
			pDest_Line += oDest.m_iHeight;
			int iRight = r + 1;

			if (bFast_Row_Process)
			{
				//��[0]�ƶ���[r]
				int iMove_To = r + 1;
				for (x = 1; x < iMove_To; x++)
				{
					iTotal += pSource_Line[iRight++] - pSource_Line[0];
					*pDest_Line = iTotal;
					pDest_Line += oDest.m_iHeight;
				}
				iLeft += r;

				//��[r+1]�ƶ���[width-1-r]
				iMove_To = iWidth_Minus_1 - r;
				for (; x < iMove_To; x++)
				{
					iTotal += pSource_Line[iRight++] - pSource_Line[iLeft++];
					*pDest_Line = iTotal;
					pDest_Line += oDest.m_iHeight;
				}

				//��[width-1-r]�ƶ���[width-1]
				for (; x < oSource.m_iWidth; x++)
				{
					iTotal += pSource_Line[iWidth_Minus_1] - pSource_Line[iLeft++];
					*pDest_Line = iTotal;
					pDest_Line += oDest.m_iHeight;
				}
			}else
			{
				for (x = 1; x < oSource.m_iWidth; x++)
				{
					iTotal -= iLeft < 0 ? pSource_Line[0] : pSource_Line[iLeft];
					iTotal += iRight >= oSource.m_iWidth ? pSource_Line[iWidth_Minus_1] : pSource_Line[iRight];
					*pDest_Line = iTotal;
					iLeft++, iRight++;
					pDest_Line += oDest.m_iHeight;
				}
			}
		}

		//����һ�Σ��������ƽ�
		for (int y = 0; y < oSource.m_iWidth; y++)
		{
			int x, iLeft = -r;	// , iRight = r;
			unsigned int* pSource_Line = &pSum_Buffer[y*oSource.m_iHeight];
			//�˴�Ҫ�ĸģ���Dest��Ϊת�þ���
			unsigned char* pDest_Line = &oDest.m_pChannel[iChannel][y];
			unsigned int iTotal = r * pSource_Line[0];
			if (oSource.m_iHeight < r + 1)
			{//������������,�ݲ�����
				for (x = 0; x < oSource.m_iHeight; x++)
					iTotal += pSource_Line[x];
				iTotal += pSource_Line[iHeight_Minus_1] * (r + 1 - oSource.m_iHeight);
			}else
			{
				for (x = 0; x <= r; x++)
					iTotal += pSource_Line[x];
			}
			pDest_Line[0] = (unsigned char)(iTotal * iStencil_Size_Recip + 0.5f);
			pDest_Line += oDest.m_iWidth;

			int iRight = r + 1;
			//����Ҫ�ǳ����ѣ��Ǵ�ת�õ���������
			if (bFast_Col_Process)
			{//�˴��������ÿ��ٷ���
				int iMove_To = r + 1;
				for (x = 1; x < iMove_To; x++)
				{
					iTotal += pSource_Line[iRight++] - pSource_Line[0];
					//ע�ⲻҪɾ�������У�����������
					//*pDest_Line = iTotal / iStencil_Size;

					*pDest_Line =(unsigned char)(iTotal * iStencil_Size_Recip + 0.5f);
					pDest_Line += oDest.m_iWidth;
				}
				iLeft += r;

				iMove_To = iHeight_Minus_1 - r;
				for (; x < iMove_To; x++)
				{
					iTotal += pSource_Line[iRight++] - pSource_Line[iLeft++];
					//ע�ⲻҪɾ�������У�����������
					//*pDest_Line = iTotal / iStencil_Size;

					*pDest_Line = (unsigned char)(iTotal * iStencil_Size_Recip +0.5f);
					pDest_Line += oDest.m_iWidth;
				}

				//��[Height-1-r]�ƶ���[Height-1]
				for (; x < oSource.m_iHeight; x++)
				{
					iTotal += pSource_Line[iHeight_Minus_1] - pSource_Line[iLeft++];
					//ע�ⲻҪɾ�������У�����������
					//*pDest_Line = iTotal / iStencil_Size;

					*pDest_Line = (unsigned char)(iTotal * iStencil_Size_Recip + 0.5f);
					pDest_Line += oDest.m_iWidth;
				}
			}else
			{
				for (x = 1; x < oSource.m_iHeight; x++)
				{
					//if (x == 599 && y == 150)
					//printf("here");
					iTotal -= iLeft < 0 ? pSource_Line[0] : pSource_Line[iLeft];
					iTotal += iRight >= oSource.m_iHeight ? pSource_Line[iHeight_Minus_1] : pSource_Line[iRight];
					//ע�ⲻҪɾ�������У�����������
					*pDest_Line = iTotal / iStencil_Size;

					//�Ը������ǲ�ȷ���ģ������и��ҵ�̬��
					//*pDest_Line = (unsigned char)(iTotal * iStencil_Size_Recip + 0.5f);
					iLeft++, iRight++;
					pDest_Line += oDest.m_iWidth;
				}
			}			
		}
	}
	Free(pSum_Buffer);
	return;
}

void Box_Filter_1(Image oSource, Image oDest,int r)
{//�������з��룬����ǱȽϼ򵥵��������������Ż��������㷨��֤
	unsigned int* pSum_Buffer = (unsigned int*)pMalloc(oSource.m_iWidth * oSource.m_iHeight*sizeof(unsigned int));
	const int iWidth_Minus_1 = oSource.m_iWidth - 1,
		iHeight_Minus_1 = oSource.m_iHeight -1,
		iStencil_Size = (2 * r + 1) * (2 * r + 1);

	for (int iChannel = 0; iChannel < oSource.m_iChannel_Count; iChannel++)
	{
		//�ȸ���
		for (int y = 0; y < oSource.m_iHeight; y++)
		{
			int x, iLeft = -r;	// , iRight = r;
			unsigned char* pSource_Line = &oSource.m_pChannel[iChannel][y*oSource.m_iWidth];
			unsigned int* pDest_Line = &pSum_Buffer[y * oSource.m_iWidth];
			/*if (y == 40)
			printf("here");*/
			//�ȸ��һ����
			unsigned int iTotal = r * pSource_Line[0];
			if (oSource.m_iWidth<r + 1 )
			{//������������,�ݲ�����
				for (x = 0; x < oSource.m_iWidth; x++)
					iTotal += pSource_Line[x];
				iTotal += pSource_Line[oSource.m_iWidth - 1] * (r + 1 - oSource.m_iWidth);
			}else
			{
				for (x = 0; x <= r; x++)
					iTotal += pSource_Line[x];
			}
			pDest_Line[0] = iTotal;
			int iRight = r + 1;

			for (x = 1; x < oSource.m_iWidth; x++)
			{
				/*if (y == 40 && x == 40)
				printf("here");*/
				iTotal -= iLeft < 0 ? pSource_Line[0] : pSource_Line[iLeft];
				iTotal += iRight >= oSource.m_iWidth ? pSource_Line[iWidth_Minus_1] : pSource_Line[iRight];
				pDest_Line[x] = iTotal;
				iLeft++, iRight++;	
			}
		}

		//�ٸ���
		for (int x = 0; x < oSource.m_iWidth; x++)
		{
			/*if (x == 40)
			printf("here");*/
			int y, iTop = -r;//
			unsigned int* pSource_Start = &pSum_Buffer[x];
			unsigned char* pDest_Line = &oDest.m_pChannel[iChannel][x];
			unsigned int iTotal = r * pSource_Start[0];
			if (oSource.m_iHeight < r + 1)
			{//������������,�ݲ�����
				for (y = 0; y < oSource.m_iHeight; y++)
					iTotal += pSource_Start[y*oSource.m_iWidth];
				iTotal += pSource_Start[iHeight_Minus_1*oSource.m_iWidth] * (r + 1 - oSource.m_iWidth);
			}else
			{
				for (y = 0; y <= r; y++)
					iTotal += pSource_Start[y*oSource.m_iWidth];
			}
			pDest_Line[0] = iTotal / iStencil_Size;
			int iBottom = r+1;

			for (y = 1; y < oSource.m_iHeight; y++)
			{
				/*if(y==40 && x==40)
				printf("here");*/
				iTotal -= iTop < 0 ? *pSource_Start : pSource_Start[iTop*oSource.m_iWidth];
				iTotal += iBottom >= oSource.m_iHeight ? pSource_Start[iHeight_Minus_1*oSource.m_iWidth] : pSource_Start[iBottom*oSource.m_iWidth];
				pDest_Line[y * oSource.m_iWidth] = iTotal / iStencil_Size;
				/*if (iTotal / iStencil_Size <100)
				printf("%d ",iTotal / iStencil_Size);*/
				iTop++, iBottom++;	
			}
		}
		//Disp(pSum_Buffer, oSource.m_iHeight, oSource.m_iWidth, "Data");
	}
	Free(pSum_Buffer);
	return;
}
int Compare_Image(Image oSource, Image oDest, int iDiff_Threshold)
{
	if (oSource.m_iWidth != oDest.m_iWidth || oSource.m_iHeight != oDest.m_iHeight /*|| oSource.m_iChannel_Count != oDest.m_iChannel_Count*/)
	{
		printf("Size Mismatched\n");
		return 0;
	}
	int i, y, x;
	int iCount = 0;
	for (i = 0; i < oSource.m_iChannel_Count; i++)
		for (y = 0; y < oSource.m_iHeight; y++)
			for (x = 0; x < oSource.m_iWidth; x++)
			{
				//printf("%d\n", oSource.m_pChannel[i][y * oSource.m_iWidth + x] - oDest.m_pChannel[i][y * oDest.m_iWidth + x]);
				if (abs(oSource.m_pChannel[i][y * oSource.m_iWidth + x] - oDest.m_pChannel[i][y * oDest.m_iWidth + x]) > iDiff_Threshold)
				{
					//if(x%32==0 &&i==0)
					printf("Mismatched: i:%d y:%d x:%d Source:%d Dest:%d Diff:%d\n", i, y, x, oSource.m_pChannel[i][y * oSource.m_iWidth + x], oDest.m_pChannel[i][y * oDest.m_iWidth + x], abs(oSource.m_pChannel[i][y * oSource.m_iWidth + x] - oDest.m_pChannel[i][y * oDest.m_iWidth + x]));
					return 0;
					iCount++;
				}
			}
	printf("Mismatched count:%d\n", iCount);

	return 1;
}
int Compare_Image(const char* pcSource, const char* pcDest,int iDiff_Threshold)
{
	Image oSource, oDest;
	int iResult = 0;
	bLoad_Image(pcSource, &oSource);
	bLoad_Image(pcDest, &oDest);
	iResult = Compare_Image(oSource, oDest, iDiff_Threshold);
	Free_Image(&oSource);
	Free_Image(&oDest);	
	return iResult;
}
void Dilate_Bin_Ref(Image oSource, Image oDest, int r)
{//�ð뾶Ϊr�ĽṹԪ��oSource����ɫΪ�׽�������
	int iHeight_Minus_1 = oSource.m_iHeight - 1,
		iWidth_Minus_1 = oSource.m_iWidth - 1;
	unsigned char* pSource = oSource.m_pChannel[0],
		* pDest = oDest.m_pChannel[0];
	for (int y = 0; y < oSource.m_iHeight; y++)
	{
		for (int x = 0; x < oSource.m_iWidth; x++)
		{
			int bFound = 0;
			for (int i = -r; i <= r; i++)
			{
				int y_s = y + i;
				unsigned char* pSource_1 = &pSource[Clip3(0, iHeight_Minus_1, y_s) * oSource.m_iWidth];
				for (int j = -r; j <= r; j++)
				{
					int x_s = x + j;
					x_s = Clip3(0, iWidth_Minus_1, x_s);
					if (pSource_1[x_s])
					{//��ɫ��Ϊ0
						bFound = 1;
						goto SUCCESS;
					}
				}
			}
		SUCCESS:
			pDest[y * oDest.m_iWidth + x] = bFound * 255;
		}
	}
	return;
}
void Dilate_Bin(Image oSource, int r, Image oDest)
{//�о���̫�鷳�ˣ�����������򵥵ģ���ȫ������Դ���Ե���Ŀ��
	if (!oDest.m_pChannel[0])
		oDest = oSource;
	else if(oSource.m_pChannel[0]!=oDest.m_pChannel[0])
		memcpy(oDest.m_pChannel[0], oSource.m_pChannel[0], oSource.m_iWidth * oDest.m_iHeight);

	const int bFast_Line_Forward = !(oDest.m_iWidth & 0x7);
	int iMove_To, x, x1, y, y1;
	unsigned char* pLine;
	//�ȸ���
	for (y = 0; y < oDest.m_iHeight; y++)
	{
		pLine = &oDest.m_pChannel[0][y*oDest.m_iWidth];
		for (x = 1; x < oDest.m_iWidth;/* x++*/)
		{
			if (bFast_Line_Forward && *(unsigned long long*) & pLine[x] == pLine[x - 1] * 0x0101010101010101)
			{
				x += 8;
				continue;
			}

			if (pLine[x] != pLine[x - 1])
			{//�б仯��Ҫ����
				if (pLine[x])
				{//��ɫ���������
					iMove_To = x - r;
					if (iMove_To < 0)
						iMove_To = 0;
					for (x1 = x-1; x1 >= iMove_To && !pLine[x1]; x1--)
						pLine[x1] = 0xFF;
				}else
				{//��ɫ����ǰ�ƽ�
					iMove_To = x + r;
					if (iMove_To > oDest.m_iWidth)
						iMove_To = oDest.m_iWidth;
					for(x1=x;x1<iMove_To && !pLine[x1];x1++)
						pLine[x1] = 0xFF;
					x = x1;
				}
			}
			x++;
		}
	}

	//bSave_Image("c:\\tmp\\1.bmp", oDest);
	//�ȸ���
	for (x = 0; x < oDest.m_iWidth; x++)
	{
		unsigned char* pCur = &oDest.m_pChannel[0][x+oDest.m_iWidth];
		for (y = 1; y < oDest.m_iHeight; y++,pCur+=oDest.m_iWidth)
		{
			//unsigned char* pCur = &oDest.m_pChannel[0][y*oDest.m_iWidth+x];
			if (*pCur != pCur[-oDest.m_iWidth])
			{
				if (*pCur)
				{//��ɫ������
					iMove_To = y - r;
					if (iMove_To < 0)
						iMove_To = 0;
					unsigned char* pPrev = pCur-oDest.m_iWidth;
					for (y1 = y - 1; y1 >= iMove_To && !*pPrev; y1--)
						*pPrev = 0xFF, pPrev -= oDest.m_iWidth;
				}else
				{//��ɫ������
					iMove_To = y + r;
					if (iMove_To > oDest.m_iHeight)
						iMove_To = oDest.m_iHeight;
					unsigned char* pNext = pCur;
					for(y1=y;y1<iMove_To && !(*pNext);y1++,pNext+=oDest.m_iWidth)
						*pNext = 0xFF;
					y = y1;
					pCur = &oDest.m_pChannel[0][y*oDest.m_iWidth+x];
				}
			}
		}
	}
	//bSave_Image("c:\\tmp\\1.bmp", oDest);
	return;
}

void Draw_Rect(Image oImage, int x, int y, int w, int h, int iThickness,int r,int g,int b)
{//��������׳��ûʱ���廨����ͷ����ʰ
	int x1, x1_End;
	unsigned char Color[] = { (unsigned char)r,(unsigned char)g,(unsigned char)b };

	int iWidth_Minus_1 = oImage.m_iWidth - 1,
		iHeight_Minus_1 = oImage.m_iHeight - 1;
	x1_End = x + w;
	if (x1_End >oImage.m_iWidth)
		x1_End = oImage.m_iWidth;

	int x2 = x, x2_End = x + iThickness;
	x2_End = Clip3(0, oImage.m_iWidth, x2_End);

	int x3 = x + w - iThickness, x3_End = x3 + iThickness;
	x3_End=Clip3(0, oImage.m_iWidth, x3_End);


	for (int iChannel = 0; iChannel < oImage.m_iChannel_Count; iChannel++)
	{
		unsigned char iColor = Color[iChannel];
		int y1, y_End;
		//�������ˮƽֱ��
		y1 = Clip3(0, iHeight_Minus_1, y);
		y_End = Clip3(0, oImage.m_iHeight, y + iThickness);
		for (; y1 <y_End; y1++)
		{
			unsigned char* pLine = &oImage.m_pChannel[iChannel][y1 * oImage.m_iWidth];
			for (x1 = x; x1 < x1_End; x1++)
				pLine[x1] = iColor;
		}

		y1 = y + iThickness;
		y_End = y + h -iThickness;
		y1 = Clip3(0, iHeight_Minus_1,y1);
		y_End = Clip3(0,  oImage.m_iHeight, y_End);
		for (; y1 < y_End; y1++)
		{
			unsigned char* pLine = &oImage.m_pChannel[iChannel][y1 * oImage.m_iWidth];
			for (x1 = x2; x1 < x2_End; x1++)
				pLine[x1] = iColor;

			for (x1 = x3; x1 < x3_End; x1++)
				pLine[x1] = iColor;
		}

		y1 = y + h - iThickness;
		y_End = y + h;
		y1 = Clip3(0, iHeight_Minus_1,y1);
		y_End = Clip3(0,  oImage.m_iHeight, y_End);
		for (; y1 <y_End; y1++)
		{
			unsigned char* pLine = &oImage.m_pChannel[iChannel][y1 * oImage.m_iWidth];
			for (x1 = x; x1 < x1_End; x1++)
				pLine[x1] = iColor;
		}
	}
	//bSave_Image("c:\\tmp\\1.bmp", oImage);
}

void Crop_Image(Image oSource, int x, int y, int w, int h, Image oDest)
{//��Source��һ�鿽����Dest
	if (oSource.m_iChannel_Count != oDest.m_iChannel_Count)
	{
		printf("Error in Crop_Image\n");
		return;
	}

	int n = oSource.m_iWidth,
		iDest_Stride = oDest.m_iWidth;

	for (int i = 0; i < oSource.m_iChannel_Count; i++)
	{
		unsigned char * pSource_Cur = &oSource.m_pChannel[i][y * n + x],
			* pSource_End = pSource_Cur + h * n,
			* pDest_Cur = oDest.m_pChannel[i];

		for (; pSource_Cur < pSource_End; pSource_Cur += n, pDest_Cur += iDest_Stride)
			for (int i = 0; i < w; i++)
				pDest_Cur[i] = pSource_Cur[i];
	}
}

template<typename _T>void Bi_Linear(_T Source[], int w_s, int h_s, int w_d, int h_d, _T Dest[])
{
	int x, y;
	float x_s_f, y_s_f, //Ŀ��x,y��Ӧ��Դ����
		w0, w1, w2, w3; //����Ȩ��  w1 = 1-w0   w3 = 1-w2
	int x_s_0, x_s_1,   //Ŀ��x�����Ӧ�����ҵ�����
		y_s_0, y_s_1;   //Ŀ��y�����Ӧ����������
	int w_1 = w_s - 1,             //w-1
		h_1 = h_s - 1;            //h-1

	//����ϵ�����ڼ����Ӧԭλ��
	float f_y = (float)w_s / w_d;
	float f_x = (float)w_s / w_d;

	_T* pDest_Cur = Dest;
	for (y = 0; y < h_d; y++)
	{
		//opencv����
		//y_s_f = (y + 0.5f) * f_y - 0.5f;

		//һ�㷽��
		y_s_f = y * f_y;

		y_s_0 = (int)y_s_f;
		y_s_1 = Min(y_s_0 + 1, h_1);
		w2 = (y_s_1 - y_s_f), w3 = 1.f - w2;
		_T* pCur_Line = &Source[y_s_0 * w_s];
		_T* pNext_Line = &Source[y_s_1 * w_s];

		//printf("%d %d\n",y_s_0,y_s_1);
		for (x = 0; x < w_d; x++, pDest_Cur++)
		{
			//opencv����
			//x_s_f = (x + 0.5f) * f_x - 0.5f;

			//һ�㷽��
			x_s_f = x  * f_x;

			x_s_0 = (int)x_s_f;
			x_s_1 = Min(x_s_0 + 1, w_1);
			w0 = (x_s_1 - x_s_f), w1 = 1.f - w0;
			float fValue_0 = (float)(w0 * pCur_Line[x_s_0] + w1 * pCur_Line[x_s_1]);
			float fValue_1 = (float)(w0 * pNext_Line[x_s_0] + w1 * pNext_Line[x_s_1]);
			*pDest_Cur =(_T)( w2 * fValue_0 + w3 * fValue_1);
		}
	}

	return;
}
int iGet_Border_x(int x, int iWidth,  Border_Type iBorder_Type)
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

	/*if (x < 0)
	{
		switch (iBorder_Type)
		{
		case Border_Type::BORDER_REFLECT:
			return -x - 1;
		case Border_Type::BORDER_REFLECT_101:
			return -x;
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
		}
	}
	return x;*/
}

int iGet_Border_y(int y, int iHeight,  Border_Type iBorder_Type)
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

	/*if (y < 0)
	{
		switch (iBorder_Type)
		{
		case Border_Type::BORDER_CONSTANT:
			return -1;
		case Border_Type::BORDER_REFLECT:
			return -y-1;
		case Border_Type::BORDER_REFLECT_101:
			return -y;
		}
	}else if(y>=iHeight)
	{
		switch (iBorder_Type)
		{
		case Border_Type::BORDER_CONSTANT:
			return -1;
		case Border_Type::BORDER_REFLECT:
			return iHeight - (y - iHeight + 1);
		case Border_Type::BORDER_REFLECT_101:
			return iHeight -(y - iHeight + 1);
		}
	}
	return y;*/
}
void Bi_Linear_cv(Image oSource, Image oDest, float fScale_x, float fScale_y, Border_Type iBorder_Type)
{//����Opencv,����1
	int i, x, y;
	float x_s_f, y_s_f, //Ŀ��x,y��Ӧ��Դ����
		w0, w1, w2, w3; //����Ȩ��  w1 = 1-w0   w3 = 1-w2
	int x_s_0, x_s_1,   //Ŀ��x�����Ӧ�����ҵ�����
		y_s_0, y_s_1;   //Ŀ��y�����Ӧ����������
	int w_1 = oSource.m_iWidth - 1,             //w-1
		h_1 = oSource.m_iHeight - 1;            //h-1

	if (oSource.m_iChannel_Count != oDest.m_iChannel_Count)
	{
		printf("Error in Resize_Bi_Linear\n");
		return;
	}

	//0.45643546458763845

	//����ϵ�����ڼ����Ӧԭλ��
	float f_y = 1.f / fScale_y;
	float f_x = 1.f / fScale_x;

	unsigned char* pDest_Cur = oDest.m_pChannel[0];
	for (i = 0; i < oDest.m_iChannel_Count; i++)
	{
		for (y = 0; y < oDest.m_iHeight; y++)
		{
			y_s_f = (y + 0.5f) * f_y - 0.5f;
			//y_s_f = y * f_y;

			y_s_0 = (int)floor(y_s_f);
			y_s_1 = Min(y_s_0 + 1, h_1);
			w2 = (y_s_1 - y_s_f), w3 = 1.f - w2;

			int iCur_y = iGet_Border_y(y_s_0, oSource.m_iHeight, iBorder_Type);
			int iNext_y = iGet_Border_y(y_s_1, oSource.m_iHeight, iBorder_Type);
			unsigned char* pCur_Line = &oSource.m_pChannel[i][iCur_y * oSource.m_iWidth];
			unsigned char* pNext_Line = &oSource.m_pChannel[i][iNext_y * oSource.m_iWidth];
			for (x = 0; x < oDest.m_iWidth; x++, pDest_Cur++)
			{
				x_s_f = (x + 0.5f) * f_x - 0.5f;
				//x_s_f = x  * f_x;
				x_s_0 = (int)floor(x_s_f);
				x_s_1 = Min(x_s_0 + 1, w_1);
				w0 = (x_s_1 - x_s_f), w1 = 1.f - w0;
				//����Ϊ��ʵȡֵλ�� 
				int x_s_0_r = iGet_Border_x(x_s_0, oSource.m_iWidth, iBorder_Type),
					x_s_1_r = iGet_Border_x(x_s_1, oSource.m_iWidth, iBorder_Type);
				float fValue_0, fValue_1;

				fValue_0 = w0 * pCur_Line[x_s_0_r] + w1 * pCur_Line[x_s_1_r];
				fValue_1 = w0 * pNext_Line[x_s_0_r] + w1 * pNext_Line[x_s_1_r];
				*pDest_Cur = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5);
			}

		}
	}
}
void Bi_Linear_Ref(Image oSource, Image oDest, Border_Type iBorder_Type)
{//�߽���һ���ο������ô������λ��
//��δ����Opencv����,�д�����
	int i, x, y;
	float x_s_f, y_s_f, //Ŀ��x,y��Ӧ��Դ����
		w0, w1, w2, w3; //����Ȩ��  w1 = 1-w0   w3 = 1-w2
	int x_s_0, x_s_1,   //Ŀ��x�����Ӧ�����ҵ�����
		y_s_0, y_s_1;   //Ŀ��y�����Ӧ����������
	int w_1 = oSource.m_iWidth - 1,             //w-1
		h_1 = oSource.m_iHeight - 1;            //h-1

	if (oSource.m_iChannel_Count != oDest.m_iChannel_Count)
	{
		printf("Error in Resize_Bi_Linear\n");
		return;
	}

	//0.45643546458763845

	//����ϵ�����ڼ����Ӧԭλ��
	float f_y = (float)oSource.m_iHeight / oDest.m_iHeight;
	float f_x = (float)oSource.m_iWidth / oDest.m_iWidth;

	unsigned char* pDest_Cur = oDest.m_pChannel[0];
	for (i = 0; i < oDest.m_iChannel_Count; i++)
	{
		for (y = 0; y < oDest.m_iHeight; y++)
		{
			/*if (y == 1)
				printf("here");*/
			y_s_f = (y + 0.5f) * f_y - 0.5f;
			//y_s_f = y * f_y;

			y_s_0 = (int)floor( y_s_f);
			y_s_1 = y_s_0 + 1;	//Min(y_s_0 + 1, h_1);
			w2 = (y_s_1 - y_s_f), w3 = 1.f - w2;

			unsigned char A=1, B=1, C=1, D=1;
			int iCur_y = iGet_Border_y(y_s_0, oSource.m_iHeight, iBorder_Type);
			int iNext_y = iGet_Border_y(y_s_1, oSource.m_iHeight, iBorder_Type);
			unsigned char* pCur_Line = &oSource.m_pChannel[i][iCur_y * oSource.m_iWidth];
			unsigned char* pNext_Line = &oSource.m_pChannel[i][iNext_y * oSource.m_iWidth];
			/*if (iCur_y < 0)
				A = B = 0;
			else if (iNext_y < 0)
				C = D = 0;*/
			for (x = 0; x < oDest.m_iWidth; x++, pDest_Cur++)
			{
				x_s_f = (x + 0.5f) * f_x - 0.5f;
				//x_s_f = x  * f_x;
				x_s_0 = (int)floor(x_s_f);
				x_s_1 = Min(x_s_0 + 1, w_1);
				w0 = (x_s_1 - x_s_f), w1 = 1.f - w0;
				//����Ϊ��ʵȡֵλ�� 
				int x_s_0_r = iGet_Border_x(x_s_0, oSource.m_iWidth,  iBorder_Type),
					x_s_1_r = iGet_Border_x(x_s_1, oSource.m_iWidth,  iBorder_Type);
				float fValue_0, fValue_1;
				
				/*if (y == 0 && x == 114)
					printf("Here");*/
				/*if (x_s_0_r < 0)
					A = C = 0;
				else if (x_s_1_r < 0)
					B = D = 0;
				else
					C = D = 1;*/
				fValue_0 = w0 * pCur_Line[x_s_0_r]*A + w1 * pCur_Line[x_s_1_r]*B;
				fValue_1 = w0 * pNext_Line[x_s_0_r]*C + w1 * pNext_Line[x_s_1_r]*D;
				*pDest_Cur = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 +0.5);
			}

		}
	}
}

void Bi_Linear_Ref_1(Image oSource, Image oDest)
{//�߽���һ���ο������ô������λ��
//��δ����Opencv����,�д�����
	int i, x, y;
	float x_s_f, y_s_f, //Ŀ��x,y��Ӧ��Դ����
		w0, w1, w2, w3; //����Ȩ��  w1 = 1-w0   w3 = 1-w2
	int x_s_0, x_s_1,   //Ŀ��x�����Ӧ�����ҵ�����
		y_s_0, y_s_1;   //Ŀ��y�����Ӧ����������
	int w_1 = oSource.m_iWidth - 1,             //w-1
		h_1 = oSource.m_iHeight - 1;            //h-1

	if (oSource.m_iChannel_Count != oDest.m_iChannel_Count)
	{
		printf("Error in Resize_Bi_Linear\n");
		return;
	}

	//����ϵ�����ڼ����Ӧԭλ��
	float f_y = (float)oSource.m_iHeight / oDest.m_iHeight;
	float f_x = (float)oSource.m_iWidth / oDest.m_iWidth;

	unsigned char* pDest_Cur = oDest.m_pChannel[0];
	for (i = 0; i < oDest.m_iChannel_Count; i++)
	{
		for (y = 0; y < oDest.m_iHeight; y++)
		{
			y_s_f = (y + 0.5f) * f_y - 0.5f;
			//y_s_f = y * f_y;

			y_s_0 = (int)y_s_f;
			y_s_1 = Min(y_s_0 + 1, h_1);
			w2 = (y_s_1 - y_s_f), w3 = 1.f - w2;
			unsigned char* pCur_Line = &oSource.m_pChannel[i][y_s_0 * oSource.m_iWidth];
			unsigned char* pNext_Line = &oSource.m_pChannel[i][y_s_1 * oSource.m_iWidth];
			//printf("%d %d\n",y_s_0,y_s_1);
			for (x = 0; x < oDest.m_iWidth; x++, pDest_Cur++)
			{
				x_s_f = (x + 0.5f) * f_x - 0.5f;
				//x_s_f = x  * f_x;
				x_s_0 = (int)x_s_f;
				x_s_1 = Min(x_s_0 + 1, w_1);
				w0 = (x_s_1 - x_s_f), w1 = 1.f - w0;
				float fValue_0 = w0 * pCur_Line[x_s_0] + w1 * pCur_Line[x_s_1];
				float fValue_1 = w0 * pNext_Line[x_s_0] + w1 * pNext_Line[x_s_1];
				*pDest_Cur = (unsigned char)(w2 * fValue_0 + w3 * fValue_1);
			}
		}
	}

	/*bSave_Image("c:\\tmp\\1.bmp", oDest);
	Compare_Image((char*)"C:\\tmp\\2.bmp", (char*)"C:\\tmp\\1.bmp");*/
	return;
}

template<typename _T>static void Sep_Filter_2D_Line(_T A[], int w, int h, _T Kernel[], int iKernel_Size, _T B[])
{
	_T* pLine;
	int x, y, x1, i, w_1 = w - 1,   //w-1
		h_1 = h - 1;                //h-1
	int iPadding = iKernel_Size / 2;
	int x_End_No_Padding = w - iPadding;

	_T fTotal;
	//Disp(A, h, w);
	for (y = 0; y < h; y++)
	{//���ϵ��£������߸�
		pLine = &A[y * w];
		//�ȸ�Padding
		for (x = 0; x < iPadding && x < w; x++)
		{
			for (fTotal = 0, i = 0, x1 = x - iPadding; i < iKernel_Size; x1++, i++)
			{
				if (x1 < 0)
				{
					if (-x1 >= w)
						fTotal += Kernel[i] * pLine[w_1];
					else
						fTotal += Kernel[i] * pLine[-x1];
				}
				else if (x1 >= w)
				{
					if (x - (x1 - x)<0)
						fTotal += Kernel[i] * pLine[0];
					else
						fTotal += Kernel[i] * pLine[w_1 - x1];
				}
				else
					fTotal += Kernel[i] * pLine[x1];
			}
			B[x * h + y] = fTotal;
		}

		//�˴���Padding, �ʺϿ����ƽ�
		for (; x < x_End_No_Padding; x++)
		{
			for (fTotal = 0, x1 = x - iPadding, i = 0; i < iKernel_Size; i++, x1++)
				fTotal += Kernel[i] * pLine[x1];
			B[x * h + y] = fTotal;
		}

		for (; x < w; x++)
		{
			for (fTotal = 0, i = 0, x1 = x - iPadding; i < iKernel_Size; x1++, i++)
			{
				if (x1 < 0)
				{
					if (-x1 >= w)
						fTotal += Kernel[i] * pLine[w_1];
					else
						fTotal += Kernel[i] * pLine[-x1];
				}
				else if (x1 >= w)
				{
					int iDup_Pos = w_1 - (x1 - w_1);
					if (iDup_Pos < 0)
						fTotal += Kernel[i] * pLine[0];
					else
						fTotal += Kernel[i] * pLine[iDup_Pos];
				}
				else
					fTotal += Kernel[i] * pLine[x1];
			}
			B[x * h + y] = fTotal;
		}
		//printf("Here");
	}
	return;
}

template<typename _T>void Sep_Filter_2D(_T A[], int w, int h, _T Ker_x[], _T Ker_y[], int iKernel_Size, _T B[])
{//���з��������������Ҫ��Դ��Ŀ�ı���ͬ����������
//Ϊ��ȫ������Դ��Ŀ�������ͬ
	_T* pMid = (_T*)pMalloc(w * h * sizeof(_T));

	//Disp(A, h, w);
	Sep_Filter_2D_Line(A, w, h, Ker_x, iKernel_Size, pMid);
	//Disp(pMid, w,h);

	//Matrix_Transpose(pMid, w, h, pMid);
	//Disp(pMid, h,w);

	//�����Ժ�pMid ��һ����Ϊh, ��Ϊw��ת�þ��󣬼������ƽ�����
	Sep_Filter_2D_Line(pMid, h, w, Ker_y, iKernel_Size, B);
	//Disp(B, h, w);

	Free(pMid);
	return;
}

template<typename Source_Type, typename Dest_Type, typename Kernel_Typ>static void Sep_Filter_2D_Line
(Source_Type A[], int w, int h, Kernel_Typ Kernel[], int iKernel_Size, Dest_Type B[])
{
	Source_Type* pLine;
	int x, y, x1, i, w_1 = w - 1,   //w-1
		h_1 = h - 1;                //h-1
	int iPadding = iKernel_Size / 2;
	int x_End_No_Padding = w - iPadding;

	Dest_Type fTotal;
	for (y = 0; y < h; y++)
	{//���ϵ��£������߸�
		pLine = &A[y * w];

		//�ȸ�Padding
		for (x = 0; x < iPadding && x < w; x++)
		{
			for (fTotal = 0, i = 0, x1 = x - iPadding; i < iKernel_Size; x1++, i++)
			{
				if (x1 < 0)
				{
					if (-x1 >= w)
						fTotal += Kernel[i] * pLine[w_1];
					else
						fTotal += Kernel[i] * pLine[-x1];
				}
				else if (x1 >= w)
				{
					if (-x < 0)
						fTotal += Kernel[i] * pLine[0];
					else
						fTotal += Kernel[i] * pLine[w_1 - x1];
				}
				else
					fTotal += Kernel[i] * pLine[x1];
			}
			B[x * h + y] = fTotal;
		}

		//�˴���Padding, �ʺϿ����ƽ�
		for (; x < x_End_No_Padding; x++)
		{
			for (fTotal = 0, x1 = x - iPadding, i = 0; i < iKernel_Size; i++, x1++)
				fTotal += Kernel[i] * pLine[x1];
			B[x * h + y] = fTotal;
		}

		for (; x < w; x++)
		{
			for (fTotal = 0, i = 0, x1 = x - iPadding; i < iKernel_Size; x1++, i++)
			{
				if (x1 < 0)
				{
					if (-x1 >= w)
						fTotal += Kernel[i] * pLine[w_1];
					else
						fTotal += Kernel[i] * pLine[-x1];
				}
				else if (x1 >= w)
				{
					int iDup_Pos = w_1 - (x1 - w_1);
					if (iDup_Pos < 0)
						fTotal += Kernel[i] * pLine[0];
					else
						fTotal += Kernel[i] * pLine[iDup_Pos];
				}
				else
					fTotal += Kernel[i] * pLine[x1];
			}
			B[x * h + y] = fTotal;
		}
		//printf("Here");
	}
	return;
}



void Pix_Bi_Linear(Image::Part_1 oImage, float x, float y, unsigned char Pix[3], Border_Type iBorder_Type)
{
	short x1 = (short)floor(x),
		y1 = (short)floor(y);
	//h_1 = oImage.m_iHeight - 1,
	//w_1 = oImage.m_iWidth - 1;
	unsigned char A, B, C, D;
	int iCur_Line_Pos, iNext_Line_Pos;
	short xl_Pos, xr_Pos;
	//int iSize = oImage.m_iWidth * oImage.m_iHeight;

	if (iBorder_Type == BORDER_CONSTANT)
	{//���������
		iCur_Line_Pos = iNext_Line_Pos = -1;
		int y2 = y1;
		if (y2 >= 0 && y2 < oImage.m_iHeight)
			iCur_Line_Pos = y2 * oImage.m_iWidth;
		y2++;
		if (y2 >= 0 && y2 < oImage.m_iHeight)
			iNext_Line_Pos = y2 * oImage.m_iWidth;

		xl_Pos = xr_Pos = -1;
		if (x1 >= 0 && x1 < oImage.m_iHeight)
			xl_Pos = x1;
		x1++;
		if (x1 >= 0 && x1 < oImage.m_iHeight)
			xr_Pos = x1;
	}
	else
	{
		iCur_Line_Pos = iGet_Border_y(y1, oImage.m_iHeight, iBorder_Type) * oImage.m_iWidth;
		iNext_Line_Pos = iGet_Border_y(y1 + 1, oImage.m_iHeight, iBorder_Type) * oImage.m_iWidth;
		xl_Pos = iGet_Border_x(x1, oImage.m_iWidth, iBorder_Type);
		xr_Pos = iGet_Border_x(x1 + 1, oImage.m_iWidth, iBorder_Type);
	}

	for (int i = 0; i < 3; i++)
	{
		if (oImage.m_pChannel[i])
		{
			if (iCur_Line_Pos >= 0)
			{
				A = xl_Pos >= 0 ? oImage.m_pChannel[i][iCur_Line_Pos] : 0;
				B = xr_Pos >= 0 ? oImage.m_pChannel[i][iCur_Line_Pos + 1] : 0;
			}
			else
				A = B = 0;
			if (iNext_Line_Pos >= 0)
			{
				C = xl_Pos >= 0 ? oImage.m_pChannel[i][iNext_Line_Pos] : 0;
				D = xr_Pos >= 0 ? oImage.m_pChannel[i][iNext_Line_Pos] : 0;
			}
			else
				C = D = 0;

			float fValue_0, fValue_1;
			{
				float w0 = (xr_Pos - x), w1 = 1.f - w0;
				fValue_0 = w0 * A + w1 * B;
				fValue_1 = w0 * C + w1 * D;
			}

			{
				float w3 = y - y1, w2 = 1.f - w3;
				Pix[i] = (unsigned char)(w2 * fValue_0 + w3 * fValue_1 + 0.5);
			}
		}
	}
	return;
}
template<typename Source_Type, typename Dest_Type, typename Kernel_Type>void Sep_Filter_2D_1
(Source_Type A[], int w, int h, Kernel_Type Ker_x[], Kernel_Type Ker_y[], int iKernel_Size, Dest_Type B[])
{
	Dest_Type* pMid = (Dest_Type*)pMalloc(w * h * sizeof(Dest_Type));

	Sep_Filter_2D_Line(A, w, h, Ker_x, iKernel_Size, pMid);
	//Disp(pMid, w,h);

	//Matrix_Transpose(pMid, w, h, pMid);
	//Disp(pMid, 1,w);

	//�����Ժ�pMid ��һ����Ϊh, ��Ϊw��ת�þ��󣬼������ƽ�����
	Sep_Filter_2D_Line(pMid, h, w, Ker_y, iKernel_Size, B);
	//Disp(B, 1, w);

	Free(pMid);
}

void Pyr_Up_Ref(Image oSource, Image oDest)
{//�Ѿ���ȫ����opencv
	int i, x, y;
	int bHas_Remain_x = oDest.m_iWidth & 1,
		bHas_Remain_y;	//= oDest.m_iHeight & 1;

	if (oDest.m_iHeight > oSource.m_iHeight * 2)
		bHas_Remain_y = 1;
	else
		bHas_Remain_y = 0;

	//�˴�Ҫ�����м����ռ�
	short* pMid = (short*)pMalloc(oDest.m_iWidth * oSource.m_iHeight * sizeof(short));
	for (i = 0; i < oSource.m_iChannel_Count; i++)
	{
		for (y = 0; y < oSource.m_iHeight; y++)
		{
			int iPos_m = y * oDest.m_iWidth,
				iPos_s = y * oSource.m_iWidth;

			//�ȸ��������
			pMid[iPos_m] = oSource.m_pChannel[i][iPos_s] * 6 + oSource.m_pChannel[i][iPos_s + 1] * 2;
			pMid[iPos_m + 1] = (oSource.m_pChannel[i][iPos_s] + oSource.m_pChannel[i][iPos_s + 1]) * 4;

			//���м�
			for (x = 1; x < oSource.m_iWidth - 1; x++)
			{
				pMid[iPos_m + (x << 1)] = oSource.m_pChannel[i][iPos_s + x - 1] + oSource.m_pChannel[i][iPos_s + x] * 6 + oSource.m_pChannel[i][iPos_s + x + 1];
				pMid[iPos_m + (x << 1) + 1] = (oSource.m_pChannel[i][iPos_s + x] + oSource.m_pChannel[i][iPos_s + x + 1]) * 4;
			}

			//��������
			pMid[iPos_m + (x << 1)] = oSource.m_pChannel[i][iPos_s + x - 1] + oSource.m_pChannel[i][iPos_s + x] * 7;
			pMid[iPos_m + (x << 1) + 1] = oSource.m_pChannel[i][iPos_s + x] * 8;

			//������
			if (bHas_Remain_x)
			{
				pMid[iPos_m + (x << 1) + 2] = pMid[iPos_m + (x << 1) + 1];
			}
		}

		//Disp(pMid, oSource.m_iHeight, oDest.m_iWidth, "Mid");

		int y0 = -1, y1 = 0, y2 = 1;
		//int h = std::min(oSource.m_iHeight * 2, oDest.m_iHeight) / 2;
		for (y = 0; y < oSource.m_iHeight; y++, y0++, y1++, y2++)
		{
			short* r0 = &pMid[iGet_Border_y(y0, oSource.m_iHeight, BORDER_REFLECT101) * oDest.m_iWidth];
			short* r1 = &pMid[iGet_Border_y(y1, oSource.m_iHeight, BORDER_REFLECT) * oDest.m_iWidth];
			short* r2 = &pMid[iGet_Border_y(y2, oSource.m_iHeight, BORDER_REFLECT) * oDest.m_iWidth];
			/*if (y == oSource.m_iHeight-1)
				printf("%d %d %d\n",r0[1163],r1[1163],r2[1163]);*/

			for (x = 0; x < oDest.m_iWidth; x++)
			{
				/*if (y == oSource.m_iHeight - 1 && x == 1416)
					printf("here");*/

				//d0�� = ��r0 + r1*6 + r2 + 32)>>6
				oDest.m_pChannel[i][(y * 2) * oDest.m_iWidth + x] = (r0[x] + r1[x] * 6 + r2[x] + 32) >> 6;
				
				if(y*2+1<oDest.m_iHeight)
					//d1�� =   [(r1 + r2)*4 + 32]>>6
					oDest.m_pChannel[i][((y * 2) + 1) * oDest.m_iWidth + x] = ((r1[x] + r2[x]) * 4 + 32) >> 6;
			}
		}

		//���Ŀ���Ϊ��������Ҫ��һ��
		if (bHas_Remain_y)
		{
			unsigned char* r0 = &oDest.m_pChannel[i][(oDest.m_iHeight - 3) * oDest.m_iWidth];
			unsigned char* r1 = &oDest.m_pChannel[i][(oDest.m_iHeight - 1) * oDest.m_iWidth];
			memcpy(r1, r0, oDest.m_iWidth);
			//printf("Here");
		}
		//Disp(oDest.m_pChannel[i], oDest.m_iHeight, oDest.m_iWidth, "Dest");
	}
	if (pMid)
		Free(pMid);
	return;
}

void Pyr_Down_Ref(Image oSource, Image oDest, Border_Type iBorder_Type)
{//ע�⣬�߽������� BORDER_REFLECT_101
//�Ѿ���ȫ��׼opencv
	int i, x, y;
	//�˴�Ҫ�����м����ռ�
	unsigned short* pMid = (unsigned short*)pMalloc(oDest.m_iWidth * oSource.m_iHeight * sizeof(unsigned short));
	//����ˮƽ��ֵ
	for (i = 0; i < oSource.m_iChannel_Count; i++)
	{
		//Disp(oSource.m_pChannel[0], 1, 3);
		for (y = 0; y < oSource.m_iHeight; y++)
		{
			int iPos_m = y * oDest.m_iWidth,
				iPos_s = y * oSource.m_iWidth;

			for (x = 0; x < oDest.m_iWidth; x++)
			{
				int x1 = x * 2;
				pMid[iPos_m + x] = oSource.m_pChannel[i][iPos_s + x1] * 6 +  //C*6
					(oSource.m_pChannel[i][iPos_s + iGet_Border_x(x1 - 1, oSource.m_iWidth, iBorder_Type)] +            //B
						oSource.m_pChannel[i][iPos_s + iGet_Border_x(x1 + 1, oSource.m_iWidth, iBorder_Type)]) * 4 +     //D
					oSource.m_pChannel[i][iPos_s + iGet_Border_x(x1 - 2, oSource.m_iWidth, iBorder_Type)] +             //A
					oSource.m_pChannel[i][iPos_s + iGet_Border_x(x1 + 2, oSource.m_iWidth, iBorder_Type)];               //E
			}
		}
		//Disp(pMid, oSource.m_iHeight, oDest.m_iWidth, "Mid");

		//y�����ϵľ��
		unsigned short row_ID[5];
		int row_Pos[5];
		for (int j = 0; j < 3; j++)
		{
			row_ID[j] = iGet_Border_y(j - 2, oSource.m_iHeight, iBorder_Type);
			row_Pos[j] = row_ID[j] * oDest.m_iWidth;
		}

		int iPos_d = 0;
		for (y = 0; y < oDest.m_iHeight; y++)
		{
			int y1 = y * 2;
			row_Pos[3] = iGet_Border_y(y1 + 1, oSource.m_iHeight, iBorder_Type) * oDest.m_iWidth;
			row_Pos[4] = iGet_Border_y(y1 + 2, oSource.m_iHeight, iBorder_Type) * oDest.m_iWidth;
			for (x = 0; x < oDest.m_iWidth; x++, iPos_d++)
			{
				unsigned short iValue;
				iValue = pMid[row_Pos[2] + x] * 6 +
					(pMid[row_Pos[1] + x] + pMid[row_Pos[3] + x]) * 4 +
					pMid[row_Pos[0] + x] + pMid[row_Pos[4] + x];
				oDest.m_pChannel[i][iPos_d] = (iValue + 128) >> 8;
			}

			//ƽ������
			row_Pos[0] = row_Pos[2];
			row_Pos[1] = row_Pos[3];
			row_Pos[2] = row_Pos[4];
		}
	}

	//Disp(oDest.m_pChannel[0], oDest.m_iHeight, oDest.m_iWidth, "Dest");
	Free(pMid);
	return;
}