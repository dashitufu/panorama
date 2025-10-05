//�����ڴ�ص�����
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "Buddy_System.h"

int iGet_Prime_No(int iMax_Piece_Count)
{//Ѱ��һ��>=iMax_Piece_Count����С����
	return iMax_Piece_Count;	//û����Ū���Ժ���Ū
}

//static unsigned int iGet_Random_No()
//{//α���	(789, 123)	(389,621) ��ͨ���������Զ����Լ�������
//#define b 389
//#define c 621
//	static unsigned int iNum = 0xFFFFFFF;	//GetTickCount();	//����
//	return iNum = iNum * b + c;
//#undef b
//#undef c
//}

void Init_Hash_Table(Mem_Mgr* poMem_Mgr, int iMax_Piece_Count)
{//��ʼ��ɢ�б�iMax_Piece_Count�������ڴ���ܷ��������
	int bRet = 1;
	poMem_Mgr->m_iCur_Item = 1;
	poMem_Mgr->m_iHash_Size = iGet_Prime_No(iMax_Piece_Count);
	poMem_Mgr->m_pHash_Item = (struct Mem_Mgr_Hash_Item*)malloc((poMem_Mgr->m_iHash_Size + 1) * sizeof(struct Mem_Mgr_Hash_Item));
	poMem_Mgr->m_pHash_Table = (int*)malloc(poMem_Mgr->m_iHash_Size * sizeof(int));
	if (!poMem_Mgr->m_pHash_Item || !poMem_Mgr->m_pHash_Table)
	{
		bRet = 0;
		goto END;
	}
	memset(poMem_Mgr->m_pHash_Table, 0, poMem_Mgr->m_iHash_Size * sizeof(int));
	memset(poMem_Mgr->m_pHash_Item, 0, (poMem_Mgr->m_iHash_Size + 1) * sizeof(struct Mem_Mgr_Hash_Item));
END:
	if (!bRet)
	{
		if (poMem_Mgr->m_pHash_Item)
			free(poMem_Mgr->m_pHash_Item);
		if (poMem_Mgr->m_pHash_Table)
			free(poMem_Mgr->m_pHash_Table);
		poMem_Mgr->m_pHash_Item = NULL;
		poMem_Mgr->m_pHash_Table = NULL;
		poMem_Mgr->m_iHash_Size = 0;
	}
	return;
}

Mem_Mgr_Hash_Item* poFind(Mem_Mgr_Hash_Item Item[],int iSub_Block_Start, int Hash_Table[], int iHash_Size)
{//��ɢ�б����ҵ�iSub_Block_Start��Ӧ����Ŀ���������Hash_Item�ĵ�ַ
	struct Mem_Mgr_Hash_Item *poItem_Exist;
	int iPos,iHash_Pos = iSub_Block_Start % iHash_Size;
	
	if ( (iPos = Hash_Table[iHash_Pos])!=0)
	{//ɢ�б��ж���
		do {
			poItem_Exist = &Item[iPos];
			if (poItem_Exist->m_iPos == iSub_Block_Start)
				return poItem_Exist;//�ҵ���
		} while ( (iPos = (int)poItem_Exist->m_iNext)!=0);
	}
	return NULL;
}
void Remove_From_Hash_Table(Mem_Mgr_Hash_Item Item[], int* piCur_Item, int iSub_Block_Start, int Hash_Table[], int iHash_Size, int* piSize)
{//��ɢ�б���ɾ��һ��Ŀ��˳�㷵�����С
	unsigned long long iPos, iPos_1;
	int bFound = 0, iHash_Pos, iCur_Item = *piCur_Item - 1;
	struct Mem_Mgr_Hash_Item oItem_Exist, * poLast, * poPrevious = NULL;

	iHash_Pos = iSub_Block_Start % iHash_Size;
	if ( (iPos = Hash_Table[iHash_Pos])!=0)
	{//ɢ�б��ж���
		do {
			oItem_Exist = Item[iPos];
			if (oItem_Exist.m_iPos == iSub_Block_Start)
			{//�ҵ���
				*piSize = (int)oItem_Exist.m_iSub_Block_Count;
				bFound = 1;
				break;;
			}
			poPrevious = &Item[iPos];
		} while ( (iPos = oItem_Exist.m_iNext)!=0);
	}

	if (bFound)
	{
		if (poPrevious)	//�ҵ��Ľڵ�Ϊ�м�ڵ㣬�ð�
			poPrevious->m_iNext = oItem_Exist.m_iNext;
		else		//�ҵ��Ľڵ�Ϊͷ�ڵ㣬�ð�
			Hash_Table[iHash_Pos] = (int)oItem_Exist.m_iNext;

		if (iPos != iCur_Item)
		{//��β��Item �����ղ�ɾ����λ���ϣ�����ɢ�б�ˮ��
			//if (iCur_Item == 684)
				//printf("Here");
			//��Hash_Item���һ���ᵽ��ɾ�����λ����
			poLast = &Item[iPos];
			*poLast = Item[iCur_Item];	//��ʱ�������Ŀ�Ѿ��ƶ�����ɾ����λ�ã���Ҫ��һ���޸�����
			iHash_Pos = poLast->m_iPos % iHash_Size;
			poPrevious = NULL;

			if ( (iPos_1 = Hash_Table[iHash_Pos])!=0)
			{
				do {
					oItem_Exist = Item[iPos_1];
					if (oItem_Exist.m_iPos == poLast->m_iPos)
					{//
						if (poPrevious)
							poPrevious->m_iNext = iPos;
						else
							Hash_Table[iHash_Pos] = (int)iPos;
					}
					poPrevious = &Item[iPos_1];
				} while ( (iPos_1 = oItem_Exist.m_iNext)!=0);
			}
		}
		(*piCur_Item)--;
	}
	else
	{
		*piSize = 0;		//�൱�ڲ��޴���
		printf("Fail to find the Item\n");
	}
	return;
}

void Add_To_Hash_Table(struct Mem_Mgr_Hash_Item Item[], int* piCur_Item, struct Mem_Mgr_Hash_Item oItem, int Hash_Table[], int iHash_Size)
{//����ɢ�б�
	unsigned long long iPos;
	int iHash_Pos, iCur_Item = *piCur_Item;
	struct Mem_Mgr_Hash_Item oItem_Exist;
	//if (oItem.m_iPos == 13421479)
		//printf("Here");

	iHash_Pos = oItem.m_iPos % iHash_Size;

	if ( (iPos = Hash_Table[iHash_Pos])!=0)
	{//ɢ�б��ж���
		do {
			oItem_Exist = Item[iPos];
			if (oItem_Exist.m_iPos == oItem.m_iPos)
			{//�ظ��ˣ��˵㲻����ɢ�б�
				goto HASH_END;
			}
		} while ((iPos = oItem_Exist.m_iNext)!=0);
	}

	oItem.m_iNext = Hash_Table[iHash_Pos];
	Item[iCur_Item] = oItem;
	Hash_Table[iHash_Pos] = iCur_Item++;
HASH_END:
	;

	*piCur_Item = iCur_Item;
}

//��ʼ���ڴ��
//Input:	iSize��			�����ڴ�ع�����ڴ��С��С�ڵ�������ܷ�����ڴ�
//			iBlock_Size:	�û��۵���û�����С
//			iMax_Piece_Count: ����ܷ������Ƭ
//return:	���ɹ�poMem_Mgr->m_pBuffer�����ݣ�����ΪNULL���Դ��ж��Ƿ�ɹ�����
void Init_Mem_Mgr(Mem_Mgr* poMem_Mgr, unsigned long long iSize, int iBlock_Size, int iMax_Piece_Count)
{
	unsigned long long iSize_To_Allocate;
	struct Mem_Mgr_Layer* poLayer;
	int iRemain, iCur_Layer, iIndex_Size, bRet = 1;

	iSize_To_Allocate = ((iSize + iBlock_Size - 1) / iBlock_Size ) * iBlock_Size;
	memset(poMem_Mgr, 0, sizeof(Mem_Mgr));
	
	//�˴��ٸĸģ��ĳ�iMax_Piece_Count�ֽڶ������������
	poMem_Mgr->m_pOrg_Buffer = (unsigned char*)malloc(iSize_To_Allocate + iBlock_Size);	//ͻ��һ�㣬���ڶ���
	poMem_Mgr->m_pBuffer =(unsigned char*)(((unsigned long long)poMem_Mgr->m_pOrg_Buffer / iBlock_Size +1) * iBlock_Size);

	poMem_Mgr->m_iSize = iSize_To_Allocate;
	if (!poMem_Mgr->m_pBuffer)
	{
		printf("�ڴ����ʧ��\n");
		bRet = 0;
		goto END;
	}
	
	poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block = iBlock_Size;

	//�Ƚ��0��
	poMem_Mgr->m_Layer[0].m_iSub_Block_Count = (unsigned int)(iSize_To_Allocate / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
	if (poMem_Mgr->m_Layer[0].m_iSub_Block_Count > (1 << 24))
	{
		printf("������˼��û��ô�����������ӻ�����ÿ��%d�ֽڵ�����£������ֻ�ܷ���:%lld �ֽڣ�Ҫô������ÿ��Ĵ�С\n", poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block, (1ll << 24) * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
		bRet = 0;
		goto END;
	}
	iCur_Layer = 0;
	poMem_Mgr->m_Layer[0].m_iBlock_Count = (poMem_Mgr->m_Layer[0].m_iSub_Block_Count + 63) / 64;
	poMem_Mgr->m_Layer[0].m_iBytes_Per_Sub_Block = poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
	poMem_Mgr->m_Layer[0].m_iBytes_Per_Block = poMem_Mgr->m_Layer[0].m_iBytes_Per_Sub_Block * 64;
	poMem_Mgr->m_Layer[0].m_iCur_Layer = 0;
	poMem_Mgr->m_Layer[0].m_iSub_Block_Per_Block = 64;

	if (poMem_Mgr->m_Layer[0].m_iBlock_Count > 32)
	{
		iCur_Layer = 1;
		while (1)
		{
			poMem_Mgr->m_Layer[iCur_Layer].m_iSub_Block_Count = poMem_Mgr->m_Layer[iCur_Layer - 1].m_iBlock_Count;
			poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count = (poMem_Mgr->m_Layer[iCur_Layer].m_iSub_Block_Count + 31) / 32;
			poMem_Mgr->m_Layer[iCur_Layer].m_iBytes_Per_Sub_Block = poMem_Mgr->m_Layer[iCur_Layer - 1].m_iBytes_Per_Sub_Block * poMem_Mgr->m_Layer[iCur_Layer - 1].m_iSub_Block_Per_Block;
			poMem_Mgr->m_Layer[iCur_Layer].m_iBytes_Per_Block = (unsigned long long)poMem_Mgr->m_Layer[iCur_Layer].m_iBytes_Per_Sub_Block * 32;
			poMem_Mgr->m_Layer[iCur_Layer].m_iCur_Layer = iCur_Layer;
			poMem_Mgr->m_Layer[iCur_Layer].m_iSub_Block_Per_Block = 32;
			if (poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count <= 32)
				break;
			iCur_Layer++;
		}
	}
	
	poMem_Mgr->m_iLayer_Count = iCur_Layer + 1;
	
	iIndex_Size = 0;	//Ҫ64λ����
	for (iCur_Layer = 0; iCur_Layer < poMem_Mgr->m_iLayer_Count; iCur_Layer++)
		iIndex_Size += poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count * 8;
	poMem_Mgr->m_pIndex_Buffer = (unsigned char*)malloc(iIndex_Size);
	if (!poMem_Mgr->m_pIndex_Buffer)
	{
		bRet = 0;
		goto END;

	}
	memset(poMem_Mgr->m_pIndex_Buffer, 0, iIndex_Size);

	for (iIndex_Size = 0, iCur_Layer = 0; iCur_Layer < poMem_Mgr->m_iLayer_Count; iCur_Layer++)
	{
		poLayer = &poMem_Mgr->m_Layer[iCur_Layer];
		poLayer->m_pIndex = &poMem_Mgr->m_pIndex_Buffer[iIndex_Size];
		//�����һ�������������ã��൱���ڱ�����
		iRemain = poLayer->m_iSub_Block_Count % poLayer->m_iSub_Block_Per_Block;
		if (iRemain)
		{
			if (poLayer->m_iSub_Block_Per_Block == 64)
				((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] = ~(((1ll << iRemain) - 1))/*<<(64-iRemain)*/;
			else
			{
				((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] = ~((1ll << (iRemain << 1)) - 1) /*<< ((32 - iRemain) << 1)*/;
				if ((poMem_Mgr->m_iSize) % poLayer->m_iBytes_Per_Sub_Block)
				{
					((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] &= ~(3ll << ((iRemain - 1) << 1));
					((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] |= (((unsigned long long)Half_Full) << ((iRemain - 1) << 1));
				}
			}
		}
		iIndex_Size += poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count * 8;
	}

	Init_Hash_Table(poMem_Mgr, iMax_Piece_Count);
	if (!poMem_Mgr->m_pHash_Table)
	{
		bRet = 0;
		goto END;
	}

END:
	if (!bRet)
	{
		if (poMem_Mgr->m_pBuffer)
			free(poMem_Mgr->m_pBuffer);
		if (poMem_Mgr->m_pIndex_Buffer)
			free(poMem_Mgr->m_pIndex_Buffer);
		memset(poMem_Mgr, 0, sizeof(Mem_Mgr));
	}
	return;
}

//��ʼ���ڴ��
//Input:	iSize��			�����ڴ�ع�����ڴ��С��С�ڵ�������ܷ�����ڴ�
//			iBlock_Size:	�û��۵���û�����С
//			iMax_Piece_Count: ����ܷ������Ƭ
//			pBuffer_GPU�൱��Attach���Buffer
//return:	���ɹ�poMem_Mgr->m_pBuffer�����ݣ�����ΪNULL���Դ��ж��Ƿ�ɹ�����
void Init_Mem_Mgr_GPU(Mem_Mgr* poMem_Mgr, unsigned long long iSize, int iBlock_Size, int iMax_Piece_Count,void *pBuffer_GPU)
{//���Ǹ�Attach Memory
	unsigned long long iSize_To_Allocate;
	struct Mem_Mgr_Layer* poLayer;
	int iRemain, iCur_Layer, iIndex_Size, bRet = 1;

	iSize_To_Allocate = ((iSize + iBlock_Size - 1) / iBlock_Size) * iBlock_Size;
	memset(poMem_Mgr, 0, sizeof(Mem_Mgr));
		
	if (!pBuffer_GPU)
	{
		printf("NULL GPU buffer ptr\n");
		return;
	}

	//�˴��ٸĸģ��ĳ�iMax_Piece_Count�ֽڶ������������
	poMem_Mgr->m_pOrg_Buffer = (unsigned char*)pBuffer_GPU;		//malloc(iSize_To_Allocate + iBlock_Size);	//ͻ��һ�㣬���ڶ���
	poMem_Mgr->m_pBuffer = (unsigned char*)(((unsigned long long)poMem_Mgr->m_pOrg_Buffer / iBlock_Size + 1) * iBlock_Size);

	poMem_Mgr->m_iSize = iSize_To_Allocate;
	if (!poMem_Mgr->m_pBuffer)
	{
		printf("�ڴ����ʧ��\n");
		bRet = 0;
		goto END;
	}

	poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block = iBlock_Size;

	//�Ƚ��0��
	poMem_Mgr->m_Layer[0].m_iSub_Block_Count = (unsigned int)(iSize_To_Allocate / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
	if (poMem_Mgr->m_Layer[0].m_iSub_Block_Count > (1 << 24))
	{
		printf("������˼��û��ô�����������ӻ�����ÿ��%d�ֽڵ�����£������ֻ�ܷ���:%lld �ֽڣ�Ҫô������ÿ��Ĵ�С\n", poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block, (1ll << 24) * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
		bRet = 0;
		goto END;
	}
	iCur_Layer = 0;
	poMem_Mgr->m_Layer[0].m_iBlock_Count = (poMem_Mgr->m_Layer[0].m_iSub_Block_Count + 63) / 64;
	poMem_Mgr->m_Layer[0].m_iBytes_Per_Sub_Block = poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
	poMem_Mgr->m_Layer[0].m_iBytes_Per_Block = poMem_Mgr->m_Layer[0].m_iBytes_Per_Sub_Block * 64;
	poMem_Mgr->m_Layer[0].m_iCur_Layer = 0;
	poMem_Mgr->m_Layer[0].m_iSub_Block_Per_Block = 64;

	if (poMem_Mgr->m_Layer[0].m_iBlock_Count > 32)
	{
		iCur_Layer = 1;
		while (1)
		{
			poMem_Mgr->m_Layer[iCur_Layer].m_iSub_Block_Count = poMem_Mgr->m_Layer[iCur_Layer - 1].m_iBlock_Count;
			poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count = (poMem_Mgr->m_Layer[iCur_Layer].m_iSub_Block_Count + 31) / 32;
			poMem_Mgr->m_Layer[iCur_Layer].m_iBytes_Per_Sub_Block = poMem_Mgr->m_Layer[iCur_Layer - 1].m_iBytes_Per_Sub_Block * poMem_Mgr->m_Layer[iCur_Layer - 1].m_iSub_Block_Per_Block;
			poMem_Mgr->m_Layer[iCur_Layer].m_iBytes_Per_Block = (unsigned long long)poMem_Mgr->m_Layer[iCur_Layer].m_iBytes_Per_Sub_Block * 32;
			poMem_Mgr->m_Layer[iCur_Layer].m_iCur_Layer = iCur_Layer;
			poMem_Mgr->m_Layer[iCur_Layer].m_iSub_Block_Per_Block = 32;
			if (poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count <= 32)
				break;
			iCur_Layer++;
		}
	}

	poMem_Mgr->m_iLayer_Count = iCur_Layer + 1;

	iIndex_Size = 0;	//Ҫ64λ����
	for (iCur_Layer = 0; iCur_Layer < poMem_Mgr->m_iLayer_Count; iCur_Layer++)
		iIndex_Size += poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count * 8;
	poMem_Mgr->m_pIndex_Buffer = (unsigned char*)malloc(iIndex_Size);
	if (!poMem_Mgr->m_pIndex_Buffer)
	{
		bRet = 0;
		goto END;

	}
	memset(poMem_Mgr->m_pIndex_Buffer, 0, iIndex_Size);

	for (iIndex_Size = 0, iCur_Layer = 0; iCur_Layer < poMem_Mgr->m_iLayer_Count; iCur_Layer++)
	{
		poLayer = &poMem_Mgr->m_Layer[iCur_Layer];
		poLayer->m_pIndex = &poMem_Mgr->m_pIndex_Buffer[iIndex_Size];
		//�����һ�������������ã��൱���ڱ�����
		iRemain = poLayer->m_iSub_Block_Count % poLayer->m_iSub_Block_Per_Block;
		if (iRemain)
		{
			if (poLayer->m_iSub_Block_Per_Block == 64)
				((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] = ~(((1ll << iRemain) - 1))/*<<(64-iRemain)*/;
			else
			{
				((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] = ~((1ll << (iRemain << 1)) - 1) /*<< ((32 - iRemain) << 1)*/;
				if ((poMem_Mgr->m_iSize) % poLayer->m_iBytes_Per_Sub_Block)
				{
					((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] &= ~(3ll << ((iRemain - 1) << 1));
					((unsigned long long*)poLayer->m_pIndex)[poLayer->m_iBlock_Count - 1] |= (((unsigned long long)Half_Full) << ((iRemain - 1) << 1));
				}
			}
		}
		iIndex_Size += poMem_Mgr->m_Layer[iCur_Layer].m_iBlock_Count * 8;
	}

	Init_Hash_Table(poMem_Mgr, iMax_Piece_Count);
	if (!poMem_Mgr->m_pHash_Table)
	{
		bRet = 0;
		goto END;
	}

END:
	if (!bRet)
	{
		if (poMem_Mgr->m_pBuffer)
			free(poMem_Mgr->m_pBuffer);
		if (poMem_Mgr->m_pIndex_Buffer)
			free(poMem_Mgr->m_pIndex_Buffer);
		memset(poMem_Mgr, 0, sizeof(Mem_Mgr));
	}
	return;
}

void Free_Mem_Mgr(Mem_Mgr* poMem_Mgr)
{
	if (poMem_Mgr->m_pOrg_Buffer)
		free(poMem_Mgr->m_pOrg_Buffer);
	if (poMem_Mgr->m_pIndex_Buffer)
		free(poMem_Mgr->m_pIndex_Buffer);
	if (poMem_Mgr->m_pHash_Item)
		free(poMem_Mgr->m_pHash_Item);
	if (poMem_Mgr->m_pHash_Table)
		free(poMem_Mgr->m_pHash_Table);
	memset(poMem_Mgr, 0, sizeof(Mem_Mgr));
	return;
}
int iGet_Left_Free_Sub_Block_Count_32(unsigned long long iOccupancy)
{//Ѱ������ߵĿ����ӿ顣�͵�ַ��iOccupancy��λ���ߵ�ַ�ڸ�λ
	int i;// , iSub_Block_Count = 0;
	/*for (i = 0; i < 32; i++)
	{
		if (iOccupancy >> 56)
			break;
		else
			iOccupancy <<= 8;
	}
	for (; i < 32; i++)*/
	for (i = 0; i < 32; i++)
	{
		if ((iOccupancy >> 62) == Empty)
		{
			//iSub_Block_Count++;
			iOccupancy <<= 2;
		}
		else
			break;
	}
	//return iSub_Block_Count;
	return i;
}

void Get_Left_Free_Sub_Block_Count_32(unsigned long long iOccupancy, int* piSub_Block_Count_Free, Occupancy_Type* piType)
{//Return	piType: ��Block��ռ�����
	int i/*, iSub_Block_Count = 0*/;
	/*for (i = 0; i < 32; i += 4)
	{
		if (iOccupancy >> 56)
			break;
		else
			iOccupancy <<= 8;
	}
	for (; i < 32; i++)*/

	for (i = 0; i < 32; i++)
	{
		if ((iOccupancy >> 62) == Empty)
			iOccupancy <<= 2;
		else
		{
			*piType = (Occupancy_Type)(iOccupancy >> 62);
			break;
		}
	}

	if (i == 32)
		*piType = Empty;

	*piSub_Block_Count_Free = i;
	//*piSub_Block_Count_Free = iSub_Block_Count;
	return;
}

void Get_Right_Free_Sub_Block_Count_32(unsigned long long iOccupancy, int* piSub_Block_Count_Free, Occupancy_Type* piType)
{
	int i;
	/*for (i = 0; i < 32; i += 4)
	{
		if (iOccupancy & 0xFF)
			break;
		else
			iOccupancy >>= 8;
	}
	for (; i < 32; i++)*/
	for (i = 0; i < 32; i++)
	{
		if ((iOccupancy & 3) == Empty)
		{
			iOccupancy >>= 2;
		}
		else
			break;
	}
	if (i == 32)
		*piType = Empty;
	else
		*piType = (Occupancy_Type)(iOccupancy & 3);

	*piSub_Block_Count_Free = i;	//iSub_Block_Count;
}
int iGet_Right_Free_Sub_Block_Count_64(unsigned long long iOccupancy)
{
	int i;
	/*for (i = 0; i < 64; i += 8)
	{
		if (iOccupancy & 0xFF)
			break;
		else
			iOccupancy >>= 8;
	}
	for (; i < 64; i++)*/
	for (i = 0; i < 64; i++)
	{
		if (iOccupancy & 1)
			break;
		else
			iOccupancy >>= 1;
	}
	return i;
}
int iGet_Left_Free_Sub_Block_Count_64(unsigned long long iOccupancy)
{
	int i;
	/*for (i = 0; i < 64; i += 8)
	{
		if (iOccupancy >> 56)
			break;
		else
			iOccupancy <<= 8;
	}
	for (; i < 64; i++)*/
	for (i = 0; i < 64; i++)
	{
		if (iOccupancy >> 63)
			break;
		else
			iOccupancy <<= 1;
	}
	return i;
}

//int iGet_Left_Free_Sub_Block_Count_32(unsigned long long iOccupancy)
//{//Ѱ������ߵĿ����ӿ顣�͵�ַ��iOccupancy��λ���ߵ�ַ�ڸ�λ
//	int i, iSub_Block_Count = 0;
//	for (i = 0; i < 32; i++)
//	{
//		if ((iOccupancy >> 62) == Empty)
//		{
//			iSub_Block_Count++;
//			iOccupancy <<= 2;
//		}
//		else
//			break;
//	}
//	return iSub_Block_Count;
//}
//
//void Get_Left_Free_Sub_Block_Count_32(unsigned long long iOccupancy, int* piSub_Block_Count_Free, Occupancy_Type* piType)
//{//Return	piType: ��Block��ռ�����
//	int i, iSub_Block_Count = 0;
//	for (i = 0; i < 32; i++)
//	{
//		if ((iOccupancy >> 62) == Empty)
//		{
//			iSub_Block_Count++;
//			iOccupancy <<= 2;
//		}
//		else {
//			*piType = (Occupancy_Type)(iOccupancy >> 62);
//			break;
//		}
//	}
//	if (i == 32)
//		*piType = Empty;
//	*piSub_Block_Count_Free = iSub_Block_Count;
//	return;
//}
//
//void Get_Right_Free_Sub_Block_Count_32(unsigned long long iOccupancy, int* piSub_Block_Count_Free, Occupancy_Type* piType)
//{
//	int i, iSub_Block_Count = 0;
//	for (i = 0; i < 32; i++)
//	{
//		if ((iOccupancy & 3) == Empty)
//		{
//			iSub_Block_Count++;
//			iOccupancy >>= 2;
//		}
//		else
//			break;
//	}
//	if (i == 32)
//		*piType = Empty;
//	else
//		*piType = (Occupancy_Type)(iOccupancy & 3);
//
//	*piSub_Block_Count_Free = iSub_Block_Count;
//}
//int iGet_Right_Free_Sub_Block_Count_64(unsigned long long iOccupancy)
//{
//	int i, iSub_Block_Count = 0;
//	for (i = 0; i < 64; i++)
//	{
//		if (!(iOccupancy & 1))
//		{
//			iSub_Block_Count++;
//			iOccupancy >>= 1;
//		}
//	}
//	return iSub_Block_Count;
//}
//
//int iGet_Left_Free_Sub_Block_Count_64(unsigned long long iOccupancy)
//{
//	int i, iSub_Block_Count = 0;
//	for (i = 0; i < 64; i++)
//	{
//		if (!(iOccupancy >> 63))
//		{
//			iSub_Block_Count++;
//			iOccupancy <<= 1;
//		}
//	}
//	return iSub_Block_Count;
//}

int iGet_Left_Block_Start_of_Layer_0(Mem_Mgr* poMem_Mgr, int iCur_Layer, int iBlock_Start_of_Seeond_Layer)
{//һֱ��̽��0��Ŀ�ʼ����Sub_Blockλ�ã��������λ��(Sub_Block_Pos)
//iBlock_of_Second_Layer: �������������ڶ����Blockλ��
	struct Mem_Mgr_Layer* poLayer;
	int iLayer, iSub_Block_Remain;
	int iBlock_Start, iSub_Block_Start;
	unsigned long long iOccupancy;
	iBlock_Start = iBlock_Start_of_Seeond_Layer;
	iSub_Block_Start = 31;
	for (iLayer = iCur_Layer;; iLayer--)
	{
		poLayer = &poMem_Mgr->m_Layer[iLayer];
		iOccupancy = ((unsigned long long*)poLayer->m_pIndex)[iBlock_Start];

		if (iLayer)
			iSub_Block_Remain = iGet_Left_Free_Sub_Block_Count_32(iOccupancy);
		else
			iSub_Block_Remain = iGet_Left_Free_Sub_Block_Count_64(iOccupancy);
		if (iLayer == 0)
			break;

		iSub_Block_Start = 32 - iSub_Block_Remain;
		if (iSub_Block_Start)
		{//�����Լ��������ж�
			if ((iOccupancy >> ((iSub_Block_Start - 1) << 1) & 0x3) != Full)
				iSub_Block_Start--;
		}

		//�ٻ���Ϊ��һ���Block��ʼλ��
		iBlock_Start = iBlock_Start * 32 + iSub_Block_Start;
	}

	iSub_Block_Start = 64 - iSub_Block_Remain;
	return iBlock_Start * 64 + iSub_Block_Start;
}

Free_Status iGet_Free_Status(int iBlock_Free, int bHas_Neighbour, int iBlock_Need, int bHas_Remain)
{//ͨ����Ҫ�Ŀ��������������еĿ����������ٽ�����������Ƿ����㹻�Ŀռ�
	if (iBlock_Free > iBlock_Need)
		return Free_Status_Enough;

	if (bHas_Neighbour)
	{
		if (bHas_Remain)
		{
			if (iBlock_Free + 1 < iBlock_Need)
				return Free_Status_Not_Enough;
			else if (iBlock_Free > iBlock_Need + 1)
				return Free_Status_Enough;
			else
				return Free_Status_Need_To_Try;
		}
		else
		{
			if (iBlock_Free >= iBlock_Need)
				return Free_Status_Enough;
			else if (iBlock_Free + 1 >= iBlock_Need)
				return Free_Status_Need_To_Try;
			else
				return Free_Status_Not_Enough;
		}
	}
	else
	{//ֻ��������п�
		if (bHas_Remain)
			return Free_Status_Not_Enough;
		else
		{//��������
			if (iBlock_Free >= iBlock_Need)
				return Free_Status_Enough;
			else
				return Free_Status_Not_Enough;
		}
	}
	return Free_Status_Not_Enough;
}

int iGet_Right_Block_End_of_Layer_0(Mem_Mgr* poMem_Mgr, int iCur_Layer, int iBlock_Start_of_Seeond_Layer)
{//һֱ��̽��0��Ľ�������Sub_Blockλ�ã��������λ��(Sub_Block_Pos)
//iBlock_of_Second_Layer: �������������ڶ����Blockλ��
	struct Mem_Mgr_Layer* poLayer;
	int iLayer, iSub_Block_Remain;
	int iBlock_Start, iSub_Block_Start;
	unsigned long long iOccupancy;
	Occupancy_Type iOccupancy_Type;
	iBlock_Start = iBlock_Start_of_Seeond_Layer;
	iSub_Block_Start = 0;
	for (iLayer = iCur_Layer;; iLayer--)
	{
		poLayer = &poMem_Mgr->m_Layer[iLayer];
		iOccupancy = ((unsigned long long*)poLayer->m_pIndex)[iBlock_Start];

		if (iLayer)
			Get_Right_Free_Sub_Block_Count_32(iOccupancy, &iSub_Block_Remain, &iOccupancy_Type);
		else
			iSub_Block_Remain = iGet_Right_Free_Sub_Block_Count_64(iOccupancy);
		if (iLayer == 0)
			break;

		if (iOccupancy_Type == Half_Full)
			iSub_Block_Remain++;
		if (iSub_Block_Remain)
			iSub_Block_Start = iSub_Block_Remain - 1;
		else
		{
			iSub_Block_Start = 31;
			iBlock_Start--;
		}
		//�ٻ���Ϊ��һ���Block��ʼλ��
		iBlock_Start = iBlock_Start * 32 + iSub_Block_Start;
	}

	return iBlock_Start * 64 + iSub_Block_Remain - 1;
}

int bTest_Space_Layer_0(Mem_Mgr* poMem_Mgr, int iBlock_Start_Limit, int iBlock_End_Limit, unsigned long long iSize, int* piSub_Block_Start_All, int* piSub_Block_Count, int iCur_Piece)
{//Layer 0���Ĳ㲻һ�����ƺ�û��ô���ӣ������ܷ�һ�����㶨
	struct Mem_Mgr_Layer* poLayer = &poMem_Mgr->m_Layer[0];
	int iBlock_Count_Need = (int)(iSize / poLayer->m_iBytes_Per_Block),
		iBlock_Count_Free, iSub_Block_Count_Free,
		iSub_Block_Start,
		iSub_Block_Count_Need = (int)((iSize + poLayer->m_iBytes_Per_Sub_Block - 1) / poLayer->m_iBytes_Per_Sub_Block),
		iSub_Block_Layer_0_Start_All;
	int bHas_Left_Half_Full_Block, bHas_Right_Half_Full_Block,
		iLeft_Block_Start, iRight_Block_End, iSub_Block_Remain;
	int i, j, iZero_Bit_Count;
	unsigned long long iOccupancy, * pIndex_64 = (unsigned long long*)poLayer->m_pIndex;
	Free_Status iFree_Status = Free_Status_Not_Enough;

	if (iBlock_End_Limit >= (int)poLayer->m_iBlock_Count)
		iBlock_End_Limit = poLayer->m_iBlock_Count - 1;

	//if (iBlock_Count_Need<=1)
	if(iSub_Block_Count_Need<=64)
	{//������ռ���һ��Block���ò��˵�ʱ��ֱ�ӴӴӿ�ʼλ�����Ҳ���
		for (i = iBlock_Start_Limit; i <= iBlock_End_Limit; i++)
		{
			if (pIndex_64[i] == 0xFFFFFFFFFFFFFFFF)
				continue;
			//�ȹ���ǰ��Block����λ�ж��Ƿ��
			iOccupancy = pIndex_64[i];
			iZero_Bit_Count = 0;
			for (iSub_Block_Start = j = 0; j < 64 && iZero_Bit_Count < iSub_Block_Count_Need; j++, iOccupancy >>= 1)
			{
				if (iOccupancy & 1)
				{
					iSub_Block_Start = j + 1;
					iZero_Bit_Count = 0;
				}
				else
					iZero_Bit_Count++;
			}
			if (iZero_Bit_Count == iSub_Block_Count_Need)
			{
				iFree_Status = Free_Status_Enough;
				iSub_Block_Layer_0_Start_All = (i << 6) + iSub_Block_Start;
				break;
			}
			//�ٹ��ұ�Block
			if (iZero_Bit_Count && i < iBlock_End_Limit)
			{
				iOccupancy = pIndex_64[i + 1];
				for (j = 0; j < 64 && iZero_Bit_Count < iSub_Block_Count_Need; j++, iOccupancy >>= 1)
				{
					if (iOccupancy & 1)
					{//�ұ�����λ1��û�б�Ҫ����ȥ��
						iZero_Bit_Count = 0;
						break;
					}
					else
						iZero_Bit_Count++;
				}
				if (iZero_Bit_Count == iSub_Block_Count_Need)
				{
					iFree_Status = Free_Status_Enough;
					iSub_Block_Layer_0_Start_All = (i << 6) + iSub_Block_Start;
					break;
				}
			}
		}
	}else
	{//�˴���������bTest_Space_1�Ĳ��ԣ��������ҷֱ�����
		for (i = iBlock_Start_Limit; i <= iBlock_End_Limit;)
		{
			while (pIndex_64[i] == 0xFFFFFFFFFFFFFFFF && i < (int)poLayer->m_iBlock_Count && i <= iBlock_End_Limit)
				i++;
			if (i > iBlock_End_Limit)
			{
				iFree_Status = Free_Status_Not_Enough;
				break;
			}

			//��ߵ�Block
			if (pIndex_64[i])
			{
				bHas_Left_Half_Full_Block = 1;
				i++;
			}else
				bHas_Left_Half_Full_Block = 0;

			iLeft_Block_Start = i;	//��ʱ��iLeft_Block_Start������Left Half Full�飬���м䲿�ֵĿ�ʼ

			//�м��Block
			iSub_Block_Count_Free = 0;
			iBlock_Count_Free = 0;
			while (pIndex_64[i] == 0 && i <= iBlock_End_Limit && iBlock_Count_Free <= iBlock_Count_Need)
				iBlock_Count_Free++, i++;

			//�ұ�˳�����
			if (pIndex_64[i] != 0xFFFFFFFFFFFFFFFF && i < (int)poLayer->m_iBlock_Count && iBlock_Count_Free <= iBlock_Count_Need)
			{
				bHas_Right_Half_Full_Block = 1;
				iRight_Block_End = i;
				i++;
			}else
				bHas_Right_Half_Full_Block = 0;

			//�����Ѿ������Blockһ�����̽�⣬��������Sub_Blockһ��
			if (bHas_Left_Half_Full_Block)
				iSub_Block_Remain = iGet_Left_Free_Sub_Block_Count_64(pIndex_64[iLeft_Block_Start - 1]);
			else
				iSub_Block_Remain = 0;

			if (iSub_Block_Remain)
			{
				iLeft_Block_Start--;
				iSub_Block_Start = 64 - iSub_Block_Remain;
			}else
				iSub_Block_Start = 0;

			iSub_Block_Count_Free = iBlock_Count_Free * poLayer->m_iSub_Block_Per_Block + iSub_Block_Remain;
			if (bHas_Right_Half_Full_Block)
			{//˳����ұߵ�BlockҲ������
				iSub_Block_Remain = iGet_Right_Free_Sub_Block_Count_64(pIndex_64[iRight_Block_End]);
				iSub_Block_Count_Free += iSub_Block_Remain;
			}else
			{
				iRight_Block_End = i - 1;
				iSub_Block_Remain = 0;
			}
			if (iSub_Block_Count_Free >= iSub_Block_Count_Need)
			{
				iSub_Block_Layer_0_Start_All = (iLeft_Block_Start << 6) + iSub_Block_Start;
				iFree_Status = Free_Status_Enough;
				break;
			}
		}
	}
	if (iFree_Status == Free_Status_Enough)
	{
		*piSub_Block_Start_All = iSub_Block_Layer_0_Start_All;
		*piSub_Block_Count = (int)((iSize + poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block - 1) / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
		return 1;
	}else
		return 0;
}

static void Set_Index(Mem_Mgr* poMem_Mgr, unsigned int iSub_Block_Start, unsigned int iSub_Block_Count, int iBit, int iCur_Piece)
{//��0���iSub_Block_Start��ʼ����������iSub_Block_Count��ֵ�������������1:ռ�ã�0������
	unsigned int i, iBlock_Start,
		iSub_Block_Start_All, iSub_Block_End_All, iSub_Block_Remain, iSub_Block_To_Set;
	unsigned long long iSize, iSpace_Start, * pIndex_64, iOccupancy;
	struct Mem_Mgr_Layer* poLayer;

	iSpace_Start = (unsigned long long)iSub_Block_Start * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
	iBlock_Start = iSub_Block_Start >> 6;
	iSub_Block_Start = iSub_Block_Start & 63;
	iSize = iSub_Block_Count * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
	pIndex_64 = (unsigned long long*)poMem_Mgr->m_Layer[0].m_pIndex;

	if (iSub_Block_Start)
	{
		iSub_Block_Remain = 64 - iSub_Block_Start;
		iSub_Block_To_Set = Min(iSub_Block_Count, iSub_Block_Remain);

		if (iBit == 1)
			pIndex_64[iBlock_Start] |= ((1ll << iSub_Block_To_Set) - 1) << iSub_Block_Start;
		//pIndex_64[iBlock_Start] |= (0xFFFFFFFFFFFFFFFF >> iSub_Block_Start) << iSub_Block_Start;
		else
			pIndex_64[iBlock_Start] &= ~(((1ll << iSub_Block_To_Set) - 1) << iSub_Block_Start);

		iSub_Block_Count -= iSub_Block_To_Set;
		iBlock_Start++;
	}
	//if (iCur_Piece == 1)
		//printf("Here");
	if (iBit == 1)
		for (i = iBlock_Start; iSub_Block_Count > 63; iSub_Block_Count -= 64, i++)
			pIndex_64[i] = 0xFFFFFFFFFFFFFFFF;
	else
		for (i = iBlock_Start; iSub_Block_Count > 63; iSub_Block_Count -= 64, i++)
			pIndex_64[i] = 0x0;

	//�ٸ����µ�����
	if (iSub_Block_Count)
	{
		if (iBit)
			pIndex_64[i] |= ((1ll << iSub_Block_Count) - 1);
		else
			pIndex_64[i] &= ~((1ll << iSub_Block_Count) - 1);
	}

	//���ˣ��Ѿ����0����޸�
	for (i = 1; (int)i < poMem_Mgr->m_iLayer_Count; i++)
	{//����������޸�
		//�Ƚ����ߵ�
		poLayer = &poMem_Mgr->m_Layer[i];
		pIndex_64 = (unsigned long long*)poLayer->m_pIndex;
		iSub_Block_Count = (unsigned int)(iSize / poLayer->m_iBytes_Per_Sub_Block);	//�ӿ飨λ��ʼ����λ��
		iSub_Block_Start_All = iSub_Block_Start = (int)(iSpace_Start / poLayer->m_iBytes_Per_Sub_Block);
		iOccupancy = ((unsigned long long*)poMem_Mgr->m_Layer[i - 1].m_pIndex)[iSub_Block_Start];
		iBlock_Start = iSub_Block_Start >> 5;
		iSub_Block_Start &= 0x1F;
		//if (iCur_Piece == 1)
			//printf("Here");
		pIndex_64[iBlock_Start] &= ~(3ll << (iSub_Block_Start << 1));
		if (iOccupancy == 0xFFFFFFFFFFFFFFFF)
			pIndex_64[iBlock_Start] |= (unsigned long long)Full << (iSub_Block_Start << 1);
		else if (iOccupancy)
			pIndex_64[iBlock_Start] |= (unsigned long long)Half_Full << (iSub_Block_Start << 1);

		iSub_Block_Start_All++;	//��ʼ�����ƽ�һ��

		//�����һ��
		iSub_Block_End_All = iSub_Block_Start = (int)((iSpace_Start + iSize) / poLayer->m_iBytes_Per_Sub_Block);
		if (iSub_Block_Start_All <= iSub_Block_End_All)
		{
			iOccupancy = ((unsigned long long*)poMem_Mgr->m_Layer[i - 1].m_pIndex)[iSub_Block_Start];
			iBlock_Start = iSub_Block_Start >> 5;
			iSub_Block_Start &= 0x1F;
			pIndex_64[iBlock_Start] &= ~(3ll << (iSub_Block_Start << 1));
			if (iOccupancy == 0xFFFFFFFFFFFFFFFF)
				pIndex_64[iBlock_Start] |= (unsigned long long)Full << (iSub_Block_Start << 1);
			else if (iOccupancy)
				pIndex_64[iBlock_Start] |= (unsigned long long)Half_Full << (iSub_Block_Start << 1);

			iSub_Block_End_All--;
		}
		else
			continue;

		//����iSub_Block_Start_All��iSub_Block_End_All����Ͷ���
		//�Ƚ����ߵ�һ��Block������
		iBlock_Start = iSub_Block_Start_All >> 5;
		iSub_Block_Start = iSub_Block_Start_All & 0x1F;
		if (iSub_Block_Start)
		{
			iSub_Block_Remain = 32 - iSub_Block_Start;
			iSub_Block_Count = iSub_Block_End_All - iSub_Block_Start_All + 1;

			iSub_Block_To_Set = Min(iSub_Block_Count, iSub_Block_Remain);
			iSub_Block_Start_All += iSub_Block_To_Set;
			if (iBit)
				pIndex_64[iBlock_Start] |= ((1ll << (iSub_Block_To_Set << 1)) - 1) << (iSub_Block_Start << 1);
			else
				pIndex_64[iBlock_Start] &= ~(((1ll << (iSub_Block_To_Set << 1)) - 1) << (iSub_Block_Start << 1));
		}
		//�м����е�����������ȫ0
		iBlock_Start = iSub_Block_Start_All >> 5;
		if (iBit)
			for (; iSub_Block_Start_All <= ((long long)iSub_Block_End_All - 31); iSub_Block_Start_All += 32)
				pIndex_64[iBlock_Start++] = 0xFFFFFFFFFFFFFFFF;
		else
			for (; iSub_Block_Start_All <= ((long long)iSub_Block_End_All - 31); iSub_Block_Start_All += 32)
				pIndex_64[iBlock_Start++] = 0;

		//���ұ�һ������
		if (iSub_Block_Start_All <= iSub_Block_End_All)
		{
			//iSub_Block_Start_All -= 32;
			iSub_Block_Remain = iSub_Block_End_All - iSub_Block_Start_All + 1;
			iSub_Block_Start = iSub_Block_Start_All & 0x1F;
			if (iBit)
				pIndex_64[iBlock_Start] ^= ((1ll << (iSub_Block_Remain << 1)) - 1);
			else
				pIndex_64[iBlock_Start] &= ~((1ll << (iSub_Block_Remain << 1)) - 1);
		}
	}
	iCur_Piece++;
	return;
}

int bTest_Space_1(Mem_Mgr* poMem_Mgr, int iBlock_Start_Limit, int iBlock_End_Limit, unsigned int iSize, int iLayer, int* piSub_Block_Start_All, int* piSub_Block_Count, int iCur_Piece)
{//���������׼�ݹ飬ֻ�ܴ�iBlock_Start_Limit�ҵ�iBlock_End_Limit����������Χ
//Return: piSub_Block_Start_All: 0�����߿�ʼλ�ã�piSub_Block_Count��һ���ҵ����ٸ�0�����飨λ��
	struct Mem_Mgr_Layer* poLayer = &poMem_Mgr->m_Layer[iLayer];
	int //bBlock_Need_Half_Full = iSize % poLayer->m_iBytes_Per_Block ? 1 : 0,	//����Block���Ƿ�������
		iBlock_Count_Need = (int)(iSize / poLayer->m_iBytes_Per_Block),			//�����������������
		bSub_Block_Need_Half_Full = iSize % poLayer->m_iBytes_Per_Sub_Block ? 1 : 0,
		iSub_Block_Layer_0_Count_Need = (iSize + poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block - 1) / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block,
		iSub_Block_Count_Need = iSize / poLayer->m_iBytes_Per_Sub_Block;
	unsigned long long* pIndex_64 = (unsigned long long*)poLayer->m_pIndex;
	int i, iLeft_Block_Start, iSub_Block_Count_Free, iBlock_Count_Free, iSub_Block_Remain,
		iRight_Block_End, iRight_Sub_Block_End,	//�ҿ�������ţ����ӿ���32λ�е�������
		bHas_Left_Half_Full_Block, bHas_Right_Half_Full_Block, bSub_Block_Start_Is_Half_Full,
		iSub_Block_Start,
		iSub_Block_Layer_0_Start_All, iSub_Block_Layer_0_End_All;
	int iResult;
	Free_Status iFree_Status = Free_Status_Not_Enough;
	Occupancy_Type iOccupancy_Type;

	
	for (i = iBlock_Start_Limit; i <= iBlock_End_Limit && i < (int)poLayer->m_iBlock_Count;)
	{//�ȸ㶨Block�����ٸ�Sub Block����ÿһ���ٷ�����������
		while (pIndex_64[i] == 0xFFFFFFFFFFFFFFFF && i < (int)poLayer->m_iBlock_Count)
			i++;
		//if (iCur_Piece == 5133 && iLayer==1 && i>=1282)
			//printf("Here");
		if (i > iBlock_End_Limit)
			break;
		if (iSub_Block_Count_Need < 32 && iLayer>0)
		{
			//�����ǰBlock��ֱ����һ��Block
			if (iLayer > 1)
				iResult = bTest_Space_1(poMem_Mgr, i << 5, (i << 5) + 63, iSize, iLayer - 1, &iSub_Block_Layer_0_Start_All, piSub_Block_Count, iCur_Piece);
			else
				iResult = bTest_Space_Layer_0(poMem_Mgr, i << 5, (i << 5) + 63, iSize, &iSub_Block_Layer_0_Start_All, piSub_Block_Count, iCur_Piece);

			i++;
			if (iResult)
			{
				iFree_Status = Free_Status_Enough;
				break;
			}else
				iFree_Status = Free_Status_Not_Enough;
		}else
		{
			//��ߵ�Block
			if (pIndex_64[i])
			{
				bHas_Left_Half_Full_Block = 1;
				i++;
			}else
				bHas_Left_Half_Full_Block = 0;

			iLeft_Block_Start = i;	//��ʱ��iLeft_Block_Start������Left Half Full�飬���м䲿�ֵĿ�ʼ

			//�м��Block
			iSub_Block_Count_Free = 0;
			iBlock_Count_Free = 0;
			while (pIndex_64[i] == 0 && i <= iBlock_End_Limit && iBlock_Count_Free <= iBlock_Count_Need)
				iBlock_Count_Free++, i++;

			//�ұ�˳�����
			if (pIndex_64[i] != 0xFFFFFFFFFFFFFFFF && i < (int)poLayer->m_iBlock_Count && iBlock_Count_Free <= iBlock_Count_Need)
			{
				bHas_Right_Half_Full_Block = 1;
				iRight_Block_End = i;
				i++;
			}else
				bHas_Right_Half_Full_Block = 0;

			//�����Ѿ������Blockһ�����̽�⣬��������Sub_Blockһ��
			if (bHas_Left_Half_Full_Block)
				Get_Left_Free_Sub_Block_Count_32(pIndex_64[iLeft_Block_Start - 1], &iSub_Block_Remain, &iOccupancy_Type);
			else
			{
				iSub_Block_Remain = 0;
				iOccupancy_Type = Full;
			}

			if (iSub_Block_Remain || iOccupancy_Type == Half_Full)
			{
				bSub_Block_Start_Is_Half_Full = iOccupancy_Type == Half_Full ? 1 : 0;
				iLeft_Block_Start--;
				iSub_Block_Start = 32 - iSub_Block_Remain - bSub_Block_Start_Is_Half_Full;
			}else
			{
				bSub_Block_Start_Is_Half_Full = 0;
				iSub_Block_Start = 0;
			}

			iSub_Block_Count_Free = iBlock_Count_Free * poLayer->m_iSub_Block_Per_Block + iSub_Block_Remain;
			if (bHas_Right_Half_Full_Block)
			{//˳����ұߵ�BlockҲ������
				Get_Right_Free_Sub_Block_Count_32(pIndex_64[iRight_Block_End], &iSub_Block_Remain, &iOccupancy_Type);
				iSub_Block_Count_Free += iSub_Block_Remain;
				bHas_Right_Half_Full_Block = (iOccupancy_Type == Half_Full);
				if (bHas_Right_Half_Full_Block)
					iSub_Block_Remain++;
				iRight_Sub_Block_End = iSub_Block_Remain - 1;
				if (!iSub_Block_Remain && !bHas_Right_Half_Full_Block)
				{
					iRight_Block_End--;
					iRight_Sub_Block_End = 31;
				}
			}else
			{
				iRight_Block_End = i - 1;
				iRight_Sub_Block_End = 31;
				bHas_Right_Half_Full_Block = 0;
				iSub_Block_Remain = 0;
			}

			iFree_Status = iGet_Free_Status(iSub_Block_Count_Free, bHas_Left_Half_Full_Block | bHas_Right_Half_Full_Block, iSub_Block_Count_Need, bSub_Block_Need_Half_Full);
			if (iFree_Status == Free_Status_Not_Enough)
				continue;	//�����ռ䣬ת����һ����

			iSub_Block_Layer_0_Start_All = iGet_Left_Block_Start_of_Layer_0(poMem_Mgr, iLayer - 1, iLeft_Block_Start * 32 + iSub_Block_Start);
			if (iFree_Status == Free_Status_Need_To_Try)
			{//��Ҫ������̽
				iSub_Block_Layer_0_End_All = iGet_Right_Block_End_of_Layer_0(poMem_Mgr, iLayer - 1, iRight_Block_End * 32 + iRight_Sub_Block_End);
				if (iSub_Block_Layer_0_End_All - iSub_Block_Layer_0_Start_All + 1 >= iSub_Block_Layer_0_Count_Need)
					iFree_Status = Free_Status_Enough;//�㹻
				else
					iFree_Status = Free_Status_Not_Enough;
			}

			if (iFree_Status == Free_Status_Enough)
				break;
		}
	}

	//iCur_Piece++;
	if (iFree_Status == Free_Status_Enough)
	{
		*piSub_Block_Start_All = iSub_Block_Layer_0_Start_All;
		*piSub_Block_Count = (iSize + poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block - 1) / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
		return 1;
	}
	else
	{
		*piSub_Block_Start_All = -1;
		*piSub_Block_Count = 0;
		return 0;
	}
}

////�����ʼ�����ڴ�
////Input�� iSize: Ҫ�����ڴ���ֽ���
////return: void ָ�룬�û��Լ�ת��Ϊ�Լ�������
//void* pMalloc_At_End(Mem_Mgr* poMem_Mgr, unsigned int iSize)
//{//���ʵ�����鷳����0������ˡ���Ϊ�������淨��һ�㶼��С�ڴ棬��յĳ������棬����Codec
//	
//}

//�����ڴ�
//Input�� iSize: Ҫ�����ڴ���ֽ���
//return: void ָ�룬�û��Լ�ת��Ϊ�Լ�������
void* pMalloc_1(Mem_Mgr* poMem_Mgr, unsigned int iSize)
{//����Ƭ��������£�������졣��������ֻ�ߵײ�
	int iResult;
	int iLayer_0_Sub_Block_Start_All, iLayer_0_Sub_Block_Count;
	static int iCur_Piece = 0;	//�����ã�ûɶӪ��

	if (poMem_Mgr->m_iPiece_Count >= poMem_Mgr->m_iHash_Size)
	{
		printf("�����ɷ���Ƭ��\n");
		return NULL;
	}
	else if (iSize <= 0)
		return NULL;

	iResult = bTest_Space_Layer_0(poMem_Mgr, 0, poMem_Mgr->m_Layer[0].m_iBlock_Count - 1, iSize, &iLayer_0_Sub_Block_Start_All, &iLayer_0_Sub_Block_Count, iCur_Piece);
	if (!iResult)
		return NULL;
	//�ٽ�������ֵ��0���Ӧλ��1����������ʵ�������Full����Half Full
	Set_Index(poMem_Mgr, iLayer_0_Sub_Block_Start_All, iLayer_0_Sub_Block_Count, 1, iCur_Piece);

	//�����Piece�����������ɢ�б���ʵ��������������ʼλ�úͳ��ȣ������С��
	struct Mem_Mgr_Hash_Item oItem = { (unsigned long long)iLayer_0_Sub_Block_Start_All ,(unsigned long long)iLayer_0_Sub_Block_Count };
	Add_To_Hash_Table(poMem_Mgr->m_pHash_Item, &poMem_Mgr->m_iCur_Item, oItem, poMem_Mgr->m_pHash_Table, poMem_Mgr->m_iHash_Size);
	poMem_Mgr->m_iPiece_Count++;

	//if ((((unsigned char*)Allocate[i]) - oMem_Mgr.m_pBuffer) / 1024 == 13421479)
	//if (iLayer_0_Sub_Block_Start_All == 13421479)
		//printf("Here");

	iCur_Piece++;
	return poMem_Mgr->m_pBuffer + oItem.m_iPos * (unsigned long long)poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
}

//�����ڴ�
//Input�� iSize: Ҫ�����ڴ���ֽ���
//return: void ָ�룬�û��Լ�ת��Ϊ�Լ�������
void* pMalloc(Mem_Mgr* poMem_Mgr, unsigned int iSize)
{//�������뷨���ڵ�һ��Test_Space�н����еĿ�ʼ����λ����Ϣȫ���㶨
	int iResult;

	static int iCur_Piece = 0;	//�����ã�ûɶӪ��

	if (poMem_Mgr->m_iPiece_Count >= poMem_Mgr->m_iHash_Size)
	{
		printf("�����ɷ���Ƭ��\n");
		return NULL;
	}
	else if (iSize <= 0)
		return NULL;

	int iLayer_0_Sub_Block_Start_All, iLayer_0_Sub_Block_Count;

	//���ҿ���λ�õ���ʼλ��
	if(poMem_Mgr->m_iLayer_Count>1)
		iResult = bTest_Space_1(poMem_Mgr, 0, poMem_Mgr->m_Layer[poMem_Mgr->m_iLayer_Count - 1].m_iBlock_Count - 1, iSize, poMem_Mgr->m_iLayer_Count - 1, &iLayer_0_Sub_Block_Start_All, &iLayer_0_Sub_Block_Count, iCur_Piece);
	else
		iResult = bTest_Space_Layer_0(poMem_Mgr, 0, poMem_Mgr->m_Layer[poMem_Mgr->m_iLayer_Count - 1].m_iBlock_Count - 1, iSize, &iLayer_0_Sub_Block_Start_All, &iLayer_0_Sub_Block_Count, iCur_Piece);
	if (!iResult)
	{
		printf("Fail to allocate space\n");
		return NULL;
	}
		
	//if (iCur_Piece == 5132)
		//printf("Here");
	//if (iLayer_0_Sub_Block_Start_All >= 41048 * 64)
		//printf("Here");
	//�ٽ�������ֵ��0���Ӧλ��1����������ʵ�������Full����Half Full
	Set_Index(poMem_Mgr, iLayer_0_Sub_Block_Start_All, iLayer_0_Sub_Block_Count, 1, iCur_Piece);

	//�����Piece�����������ɢ�б���ʵ��������������ʼλ�úͳ��ȣ������С��
	struct Mem_Mgr_Hash_Item oItem = { (unsigned long long)iLayer_0_Sub_Block_Start_All ,(unsigned long long)iLayer_0_Sub_Block_Count };
	Add_To_Hash_Table(poMem_Mgr->m_pHash_Item, &poMem_Mgr->m_iCur_Item, oItem, poMem_Mgr->m_pHash_Table, poMem_Mgr->m_iHash_Size);
	poMem_Mgr->m_iPiece_Count++;

	iCur_Piece++;
	return poMem_Mgr->m_pBuffer + oItem.m_iPos * (unsigned long long)poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
}

//�ͷ�һ��
//Input:	p: ԭ�������ָ��
void Free(Mem_Mgr* poMem_Mgr, void* p)
{
	int iSub_Block_Start, iSub_Block_Count;
	unsigned long long /*iSize = 0, */iSpace_Start;

	if (!p || !poMem_Mgr || !poMem_Mgr->m_pBuffer)	return;
	if (((unsigned char*)p - poMem_Mgr->m_pBuffer) % poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block)
	{//��ʱ����Щ��ַ���ڸÿ��ٵ��ڴ���еĵ�һ�飬����ɼ�������Ϊ��ͷ
		printf("Pointer not found\n");
		return;
	}
	//�ȸ�0��
	iSpace_Start = ((unsigned char*)p) - poMem_Mgr->m_pBuffer;
	iSub_Block_Start = (int)(iSpace_Start / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
	Remove_From_Hash_Table(poMem_Mgr->m_pHash_Item, &poMem_Mgr->m_iCur_Item, iSub_Block_Start, poMem_Mgr->m_pHash_Table, poMem_Mgr->m_iHash_Size, &iSub_Block_Count);
	if (!iSub_Block_Count)
	{
		printf("Pointer not found\n");
		return;	//���޴���
	}
	poMem_Mgr->m_iPiece_Count--;
	Set_Index(poMem_Mgr, iSub_Block_Start, iSub_Block_Count, 0,0);
}

////�ͷ�һ��
////Input:	p: ԭ�������ָ��
//void Free(Mem_Mgr* poMem_Mgr, void* p, int iCur_Piece)
//{
//	int iSub_Block_Start, iSub_Block_Count;
//	unsigned long long iSize = 0, iSpace_Start;
//
//	if (!p)	return;
//	if (((unsigned char*)p - poMem_Mgr->m_pBuffer) % poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block)
//	{//��ʱ����Щ��ַ���ڸÿ��ٵ��ڴ���еĵ�һ�飬����ɼ�������Ϊ��ͷ
//		printf("Invalid address\n");
//		return;
//	}
//
//	//�ȸ�0��
//	iSpace_Start = ((unsigned char*)p) - poMem_Mgr->m_pBuffer;
//	iSub_Block_Start = (int)(iSpace_Start / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
//	Remove_From_Hash_Table(poMem_Mgr->m_pHash_Item, &poMem_Mgr->m_iCur_Item, iSub_Block_Start, poMem_Mgr->m_pHash_Table, poMem_Mgr->m_iHash_Size, &iSub_Block_Count);
//	poMem_Mgr->m_iPiece_Count--;
//
//	if (!iSub_Block_Count)
//	{
//		printf("���޴��ڴ�\n");
//		return;	//���޴���
//	}
//
//	Set_Index(poMem_Mgr, iSub_Block_Start, iSub_Block_Count, 0, iCur_Piece);
//}

static int iQuick_Sort_Partition(struct Mem_Mgr_Hash_Item* pBuffer, int left, int right)
{
	//Real fRef;
	unsigned long long iValue;
	int pos = right;
	struct Mem_Mgr_Hash_Item  oTemp;
	if (pBuffer == NULL)
		return -1;
	right--;
	iValue = pBuffer[pos].m_iPos;
	while (left <= right)
	{
		while (left < pos && pBuffer[left].m_iPos <= iValue)
			left++;
		while (right >= 0 && pBuffer[right].m_iPos > iValue)
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

static int iAdjust_Left(struct Mem_Mgr_Hash_Item* pStart, struct Mem_Mgr_Hash_Item* pEnd)
{
	struct Mem_Mgr_Hash_Item oTemp, * pCur_Left = pEnd - 1,
		* pCur_Right;
	struct Mem_Mgr_Hash_Item oRef = *pEnd;

	//Ϊ�˼���һ���жϣ��˴���ɨ��ȥ
	while (pCur_Left >= pStart && pCur_Left->m_iPos == oRef.m_iPos)
		pCur_Left--;
	pCur_Right = pCur_Left;
	pCur_Left--;

	while (pCur_Left >= pStart)
	{
		if (pCur_Left->m_iPos == oRef.m_iPos/* && pCur_Left!=pCur_Right*/)
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

static void Quick_Sort(struct Mem_Mgr_Hash_Item Seq[], int iStart, int iEnd)
{//ע�⣺iEnd������һ�����õ�������Ԫ�أ��������һ���Ѿ������ݵ�Ԫ�أ�������������һ��
	int iPos, iLeft;
	if (iStart < iEnd)
	{
		/*if (iStart == 9 && iEnd == 18)
		{
			for (int i = iStart; i <= iEnd; i++)
				printf("%d\n", Seq[i].m_iPos);

		}*/
		iPos = iQuick_Sort_Partition(Seq, iStart, iEnd);
		//���iPos>=iEnd, �����ϼ�����̬�����ˣ���ʱ��һ���ƽ����ҳ�������Seq[iEnd]��ȵ���Ŀ���ϵ��ұ�
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

void Disp_Mem(Mem_Mgr* poMem_Mgr,int iLayer_Count)
{//��ʾ�ڴ��ʹ������������е�Item��position�Ÿ���
//iLayer_Count: һ����ʾ���㣬���ϵ���
	int i, j;
	unsigned long long iTotal, iSpace_Used;
	struct Mem_Mgr_Hash_Item* pBuffer = (struct Mem_Mgr_Hash_Item*)malloc(poMem_Mgr->m_iPiece_Count * sizeof(struct Mem_Mgr_Hash_Item));
	memcpy(pBuffer, poMem_Mgr->m_pHash_Item + 1, poMem_Mgr->m_iPiece_Count * sizeof(struct Mem_Mgr_Hash_Item));
	
	//�Ʊ�Quick_Sort̫����Ӧ�Բ���������̬�����У�Ҫ��һ�ַ�ʽ
	Quick_Sort(pBuffer, 0, poMem_Mgr->m_iPiece_Count - 1);
	printf("End Ϊ��һ������λ��\n");
	for (i = poMem_Mgr->m_iLayer_Count - 1; i >poMem_Mgr->m_iLayer_Count- iLayer_Count-1 && i>=0; i--)
	{
		printf("Layer:%d\n", i);
		for (j = 0; j < (int)poMem_Mgr->m_Layer[i].m_iBlock_Count; j++)
			printf("%llX ", ((unsigned long long*)poMem_Mgr->m_Layer[i].m_pIndex)[j]);
		printf("\n");

		//˳����������£�Ӧ������ȫ0xFF
		////if (i == 2)
		//{	
		//	for (j = 0; j < (int)poMem_Mgr->m_Layer[i].m_iBlock_Count; j++)
		//		if (((unsigned long long*)poMem_Mgr->m_Layer[i].m_pIndex)[j] != 0xFFFFFFFFFFFFFFFF)
		//			printf("Here");
		//}
	}
	for (iSpace_Used = i = 0; i < poMem_Mgr->m_iPiece_Count; i++)
	{
		printf("Piece:%d Pos:%d block count:%d End:%d\n", i, (int)pBuffer[i].m_iPos, (int)pBuffer[i].m_iSub_Block_Count, (int)(pBuffer[i].m_iPos + pBuffer[i].m_iSub_Block_Count));
		iSpace_Used += (int)pBuffer[i].m_iSub_Block_Count;
	}
	iTotal = (unsigned long long)poMem_Mgr->m_Layer[0].m_iSub_Block_Count * poMem_Mgr->m_Layer[0].m_iBytes_Per_Sub_Block / 1024;
	iSpace_Used = (unsigned long long)iSpace_Used * poMem_Mgr->m_Layer[0].m_iBytes_Per_Sub_Block / 1024;
	printf("Total:%lldKB Used:%lldKB Remained:%lldKB\n", iTotal, iSpace_Used, iTotal - iSpace_Used);

	free(pBuffer);
	return;
}
//��ԭ������Ĵ�С��С��iSize��ԭ�������С������ڵ���iSize������ɶҲ����
//Input:	p: ԭ���ڴ�ָ��
//			iSize: �µĴ�С
void Shrink(Mem_Mgr* poMem_Mgr,void *p, unsigned int iSize)
{
	Mem_Mgr_Hash_Item* poItem;
	if (((unsigned char*)p - poMem_Mgr->m_pBuffer) % poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block)
	{
		printf("Invalid address\n");
		return;
	}
	if (iSize == 0)
	{
		Free(poMem_Mgr, p);
		return;
	}
	int iSub_Block_Start = (unsigned int)((unsigned char*)p - poMem_Mgr->m_pBuffer)/poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
	int iNew_Sub_Block_Count= (iSize + poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block - 1) / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;
	
	poItem = poFind(poMem_Mgr->m_pHash_Item, iSub_Block_Start, poMem_Mgr->m_pHash_Table, poMem_Mgr->m_iHash_Size);
	if (!poItem)
	{
		printf("Fail to find item in Shrink\n");
		return;
	}
	if (iSize >= poItem->m_iSub_Block_Count * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block)
	{
		//printf("iSize���ڵ���ԭ����ռ��С�������ı�\n");
		return;
	}

	Set_Index(poMem_Mgr, iSub_Block_Start + iNew_Sub_Block_Count, (int)(poItem->m_iSub_Block_Count - iNew_Sub_Block_Count), 0, 0);
	poItem->m_iSub_Block_Count = iNew_Sub_Block_Count;
	return;
}
//��ԭ��p���䵽���ڴ�����iSize���빬iSize<=ԭ����С����ɶҲ����
//Input:	p: ԭ���ڴ�ָ��
//			iSize: �µĴ�С
int bExpand(Mem_Mgr* poMem_Mgr, void* p, unsigned int iSize)
{
	Mem_Mgr_Hash_Item* poItem;
	int iSub_Block_Start_All = (unsigned int)((unsigned char*)p - poMem_Mgr->m_pBuffer) / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block,
		iSub_Block_Count_Need,
		iSub_Block_Count_Free,iBlock_Start,iSub_Block_Start;
	unsigned int iOrg_Size, iSize_Delta, iNew_Block_Count;
	int i, j, iZero_Bit_Count,iResult=1;
	unsigned long long iOccupancy;

	poItem = poFind(poMem_Mgr->m_pHash_Item, iSub_Block_Start_All, poMem_Mgr->m_pHash_Table, poMem_Mgr->m_iHash_Size);
	iOrg_Size = (unsigned int)(poItem->m_iSub_Block_Count * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block);
	iSize_Delta = iSize - iOrg_Size;
	if (iSize <= iOrg_Size)
	{
		printf("iSizeС���ڵ���ԭ����ռ��С�������ı�\n");
		return 0;
	}

	for (i = poMem_Mgr->m_iLayer_Count - 1; i >= 0; i--)
	{
		iNew_Block_Count =(unsigned int)(iSize_Delta / poMem_Mgr->m_Layer[i].m_iBytes_Per_Block);
		if (iNew_Block_Count)
			break;
	}
	iSub_Block_Count_Need = (iSize_Delta + poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block - 1) / poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block;

	if (!iNew_Block_Count)
	{//�������в�ûһ����ռ����һ��Block�ģ��϶���0������Ҳ���������Block
		iBlock_Start =(int) (((poItem->m_iPos + poItem->m_iSub_Block_Count) * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block) / poMem_Mgr->m_Layer[0].m_iBytes_Per_Block);
		iSub_Block_Start = (poItem->m_iPos + poItem->m_iSub_Block_Count) % 64;
		iOccupancy = ((unsigned long long*)poMem_Mgr->m_Layer[0].m_pIndex)[iBlock_Start];
		iOccupancy >>= iSub_Block_Start;
		iSub_Block_Start_All = (iBlock_Start << 6) + iSub_Block_Start;
		for (j=iSub_Block_Start, iZero_Bit_Count=0; j < 64 && iZero_Bit_Count < iSub_Block_Count_Need; j++, iOccupancy >>= 1, iZero_Bit_Count++)
			if (iOccupancy & 1)
				break;
		if (iZero_Bit_Count< iSub_Block_Count_Need && iBlock_Start+1 < (int)poMem_Mgr->m_Layer[0].m_iBlock_Count && !(iOccupancy&1) )
		{
			iBlock_Start++;
			iOccupancy = ((unsigned long long*)poMem_Mgr->m_Layer[0].m_pIndex)[iBlock_Start];
			for (j = 0; j < 64 && iZero_Bit_Count < iSub_Block_Count_Need; j++, iOccupancy >>= 1, iZero_Bit_Count++)
				if (iOccupancy & 1)
					break;
		}
		if (iZero_Bit_Count >= iSub_Block_Count_Need)
		{
			iResult = 1;
			iSub_Block_Count_Free = iSub_Block_Count_Need;
		}
	}else if (i == 0)
	{
		iBlock_Start =(int) (((poItem->m_iPos + poItem->m_iSub_Block_Count) * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block) / poMem_Mgr->m_Layer[0].m_iBytes_Per_Block);
		iResult = bTest_Space_Layer_0(poMem_Mgr, iBlock_Start, iBlock_Start + iNew_Block_Count + 1, iSize_Delta, &iSub_Block_Start_All, &iSub_Block_Count_Free, 0);
	}else
	{
		iBlock_Start = (int)(((poItem->m_iPos + poItem->m_iSub_Block_Count) * poMem_Mgr->m_iBytes_Per_Bottom_Sub_Block) / poMem_Mgr->m_Layer[i].m_iBytes_Per_Block);
		iResult = bTest_Space_1(poMem_Mgr, iBlock_Start, iBlock_Start + iNew_Block_Count + 1,iSize_Delta,i, &iSub_Block_Start_All, &iSub_Block_Count_Free, 0);
	}
	
	if (iResult && iSub_Block_Start_All == poItem->m_iPos + poItem->m_iSub_Block_Count)
	{
		Set_Index(poMem_Mgr, iSub_Block_Start_All, iSub_Block_Count_Free, 1, 0);
		poItem->m_iSub_Block_Count += iSub_Block_Count_Free;
		return 1;
	}else
		return 0;
}