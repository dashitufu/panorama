#pragma once

#ifndef Min
	#define Min(a,b) ((a)<(b)?(a):(b))
#endif
#ifndef Max
	#define Max(a,b) ((a)>(b)?(a):(b))
#endif

#define Attach_Light_Ptr(oPtr,pBuffer,iSize,iGPU_ID) oPtr = { (int)0,(int)(iSize),pBuffer,iGPU_ID }

//Light Ptr�ڴ˷���
#define Malloc(oPtr, iSize,pBuffer) \
{ \
	(pBuffer)= (unsigned char*)(oPtr).m_pBuffer+(oPtr).m_iCur; \
	(oPtr).m_iCur +=ALIGN_SIZE_128((iSize)); \
	if ( (oPtr).m_iCur > (oPtr).m_iMax_Buffer_Size) \
	{ \
		(oPtr).m_iCur -= ALIGN_SIZE_128((iSize));\
		(pBuffer)=NULL; \
		printf("Fail to allocate memory in Malloc, Total:%d Remain:%d Need:%zu\n",(oPtr).m_iMax_Buffer_Size,(oPtr).m_iMax_Buffer_Size-(oPtr).m_iCur,(size_t)(iSize)); \
	} \
}

//�������ӿ������״̬��0����ӿ�ֻ��һλ��ʾ�����Դ�ö�ٽṹû0��ʲô��
//1������ÿһ���ӿ�����λ��ɣ����Ա�ʾ����״̬
typedef enum Occupancy_Type {
	Empty = 0,				//�ӿ�Ϊ����״̬�����Է���
	Half_Full = 1,			//�ӿ�Ϊ�����״̬��������Щ������Щ�����ã�����̽һ�����
	Full = 3				//�ӿ��Ѿ���ռ��
}Occupancy_Type;

//��ĳһ���ƽ��Ĺ����У��õ�������������������֮�󣬿��Է������������
typedef enum Free_Status {	//���������
	Free_Status_Enough = 0,					//���㹻�Ŀռ�
	Free_Status_Need_To_Try = 1,			//��Ҫ��һ����̽�Ƿ����㹻���ӿ�
	Free_Status_Not_Enough = 3				//�����ռ�
}Free_Status;

typedef struct Mem_Mgr_Layer {	//�ýṹΪĳһ�������
	unsigned long long	m_iBytes_Per_Block;		//�ò�һ���������ʾ�����ֽڵ��ڴ�
	unsigned int		m_iBytes_Per_Sub_Block;	//�ò�һ�������ӿ��ʾ�����ֽ��ڴ�
	unsigned int		m_iBlock_Count;			//���㹲�ֶ���������
												//���Ƕ���0�����ϣ��鲻��Ҫ�������ӿ����λ 01,�����ӿ鶼������
	unsigned int m_iSub_Block_Count;			//���㹲�ж��ٸ��ӿ�
	unsigned int m_iCur_Layer : 3;				//��ǰ��ţ������ţ�
	unsigned int m_iSub_Block_Per_Block;		//һ��Block��Ϊ���ٸ�Sub Block��0����64��1��������32

	unsigned char* m_pIndex;					//Block����ʼλ�ã�64λ���룬����0�㣬һ��Sub_Block��һλ���ɣ�0���գ�1����ռ
												//����0�����ϵĲ㣬һ��Sub_Block����λ���ɣ�0:�գ� 11:��ռ�� 01����ռ
}Mem_Mgr_Layer;

//����λɢ�б��֡�һ���ڴ�ر���ָ������ܷ�������Piece)���������ָ��ֻ�ܷ���1023��Ŀ��
	//��ô�������1023���ʱ�򽫲�����䡣��Щ�ѷ��������һ��ɢ�б�����֯
typedef struct Mem_Mgr_Hash_Item {							//ɢ�б��ÿһ������
	unsigned long long m_iPos : 24;					//��0�����һ��Sub_Block�����ֻ�ܷ���1^24��0�������ӿ�
	unsigned long long m_iSub_Block_Count : 24;		//һ��ʹ���˵�0��Ķ��ٸ�Sub_Block
	unsigned long long m_iNext : 16;				//����
}Mem_Mgr_Hash_Item;

//ע�⣺��ĸ��������֣�һ�����û��۵�Ŀ飬���������ĸ��һ������Զ���Ϊ�����ֽڣ�����1024�ֽڡ���������Ϊ�û���
//��һ���ǳ����ڲ������õĿ飬���³�Ϊ�����顣��������64λ�޷���������ʾ��һ�������������ɸ������ӿ����
//0��һ����������64�������ӿ���ɣ�һ�������ӿ���һλ��ʾ��0λ���У�1Ϊռ��
//1��������������32�������ӿ���ɣ�һ�������ӿ�����λ��ʾ
typedef struct Mem_Mgr {
	//��ײ�Ϊ0�㣬һ��ɷ�Ϊ64���ӿ顣64λ����һ�顣ÿһ��Ϊm_iBytes_Per_Bottom_Sub_Block�ֽ�
	//�����ÿ��Ŀ��Ϊ32�ӿ�
	Mem_Mgr_Layer m_Layer[8];	//�������8��

	int m_iLayer_Count;					//���������ɶ��ٲ����
	int m_iBytes_Per_Bottom_Sub_Block;	//��ײ��һλ��Ӧ��������ֽ���
	unsigned char* m_pIndex_Buffer;		//���в����������
	unsigned long long m_iSize;			//����Mem_Mgr����������ڴ�
	unsigned char* m_pBuffer;			//����ķ�����ڴ�
	unsigned char* m_pOrg_Buffer;		//ϵͳ������ڴ濪ʼλ��
	
	Mem_Mgr_Hash_Item *m_pHash_Item;
	int m_iHash_Size;	//�����Է������Ƭ�ڴ档Ϊ����СIndex��ռ�ã���ɢ�б�����֯
	int* m_pHash_Table;	//��ɢ�б�
	int m_iPiece_Count;	//һ�������˶���Ƭ������Ϊm_iHash_Size
	int m_iCur_Item;	//m_pHash_Item��ǰ��ţ���һ������λ��

	//int m_iGPU_ID;		//���Ǹ�GPU�Ϸ���
}Mem_Mgr;

typedef struct Light_Ptr {	//˳����һ�����������Դ�ָ�룬����˳�����
	int m_iCur;					//���ܳ���m_iMax_Buffer_Size
	int m_iMax_Buffer_Size;
	unsigned char* m_pBuffer;	//�ڴ�ָ��
	int m_iGPU_ID;				//�ڴ�����GPU
}Light_Ptr;

//ȡ�����
//unsigned int iGet_Random_No();

//�����ڴ�
//Input�� iSize: Ҫ�����ڴ���ֽ���
//return: void ָ�룬�û��Լ�ת��Ϊ�Լ�������
void* pMalloc(Mem_Mgr* poMem_Mgr, unsigned int iSize);

void* pMalloc_1(Mem_Mgr* poMem_Mgr, unsigned int iSize);

//�ͷ�һ��
//Input:	p: ԭ�������ָ��
void Free(Mem_Mgr* poMem_Mgr, void* p);

//��ʼ���ڴ��
//Input:	iSize��			�����ڴ�ع�����ڴ��С��С�ڵ�������ܷ�����ڴ�
//			iBlock_Size:	�û��۵���û�����С
//			iMax_Piece_Count: ����ܷ������Ƭ
//return:	���ɹ�poMem_Mgr->m_pBuffer�����ݣ�����ΪNULL���Դ��ж��Ƿ�ɹ�����
void Init_Mem_Mgr(Mem_Mgr* poMem_Mgr, unsigned long long iSize, int iBlock_Size, int iMax_Piece_Count);
void Init_Mem_Mgr_GPU(Mem_Mgr* poMem_Mgr, unsigned long long iSize, int iBlock_Size, int iMax_Piece_Count, void* pBuffer_GPU);

//�ͷ������ڴ��
void Free_Mem_Mgr(Mem_Mgr* poMem_Mgr);
void Free_Mem_Mgr_GPU(Mem_Mgr* poMem_Mgr);

//��ԭ������Ĵ�С��С��iSize��ԭ�������С������ڵ���iSize������ɶҲ����
//Input: iSize: �µĴ�С
void Shrink(Mem_Mgr* poMem_Mgr, void* p, unsigned int iSize);

//��ԭ��p���䵽���ڴ�����iSize���빬iSize<=ԭ����С����ɶҲ����
//Input:	p: ԭ���ڴ�ָ��
//			iSize: �µĴ�С
int bExpand(Mem_Mgr* poMem_Mgr, void* p, unsigned int iSize);

//��ʾ�ڴ��ʹ������������е�Item��position�Ÿ���
void Disp_Mem(Mem_Mgr* poMem_Mgr, int iLayer_Count);
