//������ѧ��ͬʱӦ��float/double�������ͣ��ʴ˸ð濪ʼ��ȫ��Ϊģ�庯��
#pragma once 
//#include "Common.h"

#define sign(x) (x>=0?1:-1)

#define Get_Max_1(V, fMax,i) \
{ \
	double fAbs_Max; \
	fAbs_Max = Abs(fMax=V[0]); \
	i=0; \
	if (Abs(V[1]) > fAbs_Max) \
	{ \
		fAbs_Max = Abs(V[1]), fMax = V[1]; \
		i=1;	\
	} \
	if ((Abs(V[2]) > fAbs_Max)) \
	{ \
		fAbs_Max = Abs(V[2]), fMax = V[2]; \
		i=2;	\
	} \
}

//#define PI 3.14159265358979323846
#define Abs(A) ((A)>=0?(A):(-(A)))
#define MAX_FLOAT ((float)0xFFFFFFFFFFFFFFFF)	//������һ���㹻������֣�������iEEE�ĸ��������ֵ
#define ZERO_APPROCIATE	0.00001f

#define Get_Homo_Pos(Pos, Homo_Pos) Homo_Pos[0]=Pos[0], Homo_Pos[1]=Pos[1], Homo_Pos[2]=Pos[2],Homo_Pos[3]=1;

typedef struct Complex_d {//˫���ȸ���
	double real;	//ʵ��
	double im;		//�鲿 imagary part	
}Complex_d;
typedef struct Complex_f {//˫���ȸ���
	float real;	//ʵ��
	float im;		//�鲿 imagary part	
}Complex_f;

typedef struct Polynormial {	//��Ԫ����ʽ����������Ϊ������(>=0, =0��ʾû�и���
	unsigned char* m_pTerm;		//����ʽ�е���
	float* m_pCoeff;			//ÿ���ϵ��
	int m_iElem_Count;			//������Ԫ����
	int m_iTerm_Count;			//����
	int m_iMax_Term_Count;		//Ŀǰ���������Ϊ�պ���չ��
}Polynormial;

typedef struct SVD_Info {	//�б�Ҫ��SVD��һ�������Ĳ����ṹ
	void* A;		//ԭ����
	void* U;
	void* S;
	void* Vt;
	int h_A, w_A,
		h_Min_U, w_Min_U,
		w_Min_S,
		h_Min_Vt, w_Min_Vt;
	int m_bSuccess;
}SVD_Info;

template<typename _T> class Sparse_Matrix {
public:
	typedef struct Item {
	public:
		unsigned int x, y;	//0xFFFFFFFFδ�Ƿ�ֵ
		_T m_fValue;
		int m_iRow_Next;	//�з�����һ��������
		int m_iCol_Next;	//�з�����һ��������
	}Item;

	unsigned int* m_pRow=NULL;	//��Ͱ
	unsigned int* m_pCol;	//��Ͱ
	Item* m_pBuffer;		//���ݴ�[1]��ʼ [0]��ʾNULL
	int m_iCur_Item=1;		//��ǰ���ĸ�Item
	int m_iMax_Item_Count=0;		//����matrix���Ҷ��ٸ�item
	int m_iRow_Count;
	int m_iCol_Count;
};

unsigned long long iGet_Random_No_cv(unsigned long long* piState);
#define Malloc_1(oPtr, iSize,pBuffer) \
{ \
	(pBuffer)= (_T*)((oPtr).m_pBuffer+(oPtr).m_iCur); \
	(oPtr).m_iCur +=ALIGN_SIZE_128((iSize)); \
	if ( (oPtr).m_iCur > (oPtr).m_iMax_Buffer_Size) \
	{ \
		(oPtr).m_iCur -= ALIGN_SIZE_128((iSize));\
		(pBuffer)=NULL; \
		printf("Fail to allocate memory in Malloc, Total:%u Remain:%u Need:%d\n",(oPtr).m_iMax_Buffer_Size,(oPtr).m_iMax_Buffer_Size-(oPtr).m_iCur,(int)(iSize)); \
	} \
}
template<typename _T>void Disp(Sparse_Matrix<_T> oMatrix, const char Caption[] = NULL);
template<typename _T>void Disp_Fillness(Sparse_Matrix<_T> oMatrix);
template<typename _T>void Disp_Part(_T* M, int iStride, int x, int y, int w, int h, const char* pcCaption=NULL);
template<typename _T>void Disp(_T* M, int iHeight, int iWidth,const char* pcCaption=NULL)
{
	int i, j;
	if (pcCaption)
		printf("%s\n", (char*)pcCaption);

	for (i = 0; i < iHeight; i++)
	{
		for (j = 0; j < iWidth; j++)
		{
			if (typeid(_T) == typeid(float))
			{
				if (M[i * iWidth + j] - (int)M[i * iWidth + j])
					//printf("%.10ef, ", (double)M[i * iWidth + j]);
					printf("%f, ", (double)M[i * iWidth + j]);
				else
					printf("%d, ", (int)M[i * iWidth + j]);
			}else if (typeid(_T) == typeid(double))
			{
				if (M[i * iWidth + j] - (int)M[i * iWidth + j])
					printf("%f, ", (double)M[i * iWidth + j]);
					//printf("%.10ef, ", (double)M[i * iWidth + j]);
				else
					printf("%d, ", (int)M[i * iWidth + j]);
				//printf("%f,", (double)M[i * iWidth + j]);
			}				
			else if ( (typeid(_T) == typeid(unsigned int)) || 
				(typeid(_T) == typeid(int)) ||
				(typeid(_T) == typeid(short)) ||
				(typeid(_T) == typeid(unsigned short)) ||
				(typeid(_T) == typeid(unsigned char)) )
				printf("%d   ", (int)M[i * iWidth + j]);
			
		}
		printf("\n");
	}
	return;
}
template<typename _T>void Disp_Fillness(_T A[], int m, int n);

template<typename _T>void Matrix_Transpose(_T* A, int ma, int na, _T* At);

template<typename _T>void Matrix_Multiply_3x3(_T A[3*3], _T B[3*3], _T C[3*3])
{//����C=AxB
 //Light_Ptr oPtr = oMatrix_Mem;
	_T *C_1, C_2[9];
	C_1 = A == C || B == C ? C_2 : C;

	//Malloc_1(oPtr, 3 * 3 * sizeof(_T),C_1);
	C_1[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
	C_1[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
	C_1[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];

	C_1[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
	C_1[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
	C_1[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];

	C_1[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
	C_1[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
	C_1[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];

	if( A == C || B == C)
		memcpy(C, C_1, 9 * sizeof(_T));
	return;
}
template<typename _T>void Add_I_Matrix(_T A[], int iOrder, _T ramda=1);

//������Ϻ���
unsigned long long iFactorial(int n);
unsigned long long iFactorial(unsigned long long iPre_Value, int n);

unsigned int iGet_Random_No();
template<typename _T>void Get_Random_Norm_Vec(_T V[], int n);;
unsigned long long iGet_Random_No_cv();
unsigned long long iGet_Random_No_cv(unsigned long long* piState);
int iGet_Random_No_cv(int a, int b);
template<typename _T> _T fGet_Random_No(_T mean, _T sigma);

//void Schmidt_Orthogon(float* A, int m, int n, float* B);	//ʩ����������
template<typename _T>void Schmidt_Orthogon(_T* A, int m, int n, _T* B);

void Polor_2_Rect_Coordinate(float rho, float theta, float phi, float* px, float* py, float* pz);	//������תֱ������ϵ
void Polor_2_Rect_Coordinate(float rho, float theta, float* px, float* py);							//��ά
void Rect_2_Polor_Coordinate(float x, float y, float z, float* prho, float* ptheta, float* pphi);	//ֱ������ϵת������
void Rect_2_Screen_Coordinate(float x, float y, int* px_Screen, int* py_Screen, int iWidth = 1920, int iHeight = 1080);	//ֱ������ϵת��Ļ����
void Screen_2_Coordinate(int x_Screen, int y_Screen, float* px, float* py, int iWidth = 1920, int iHeight = 1080);			//��Ļ����תֱ������ϵ

//һ�������Է����йصĺ���
template<typename _T>void LLt_Decompose(_T A[], int iOrder, _T B[], int* pbSuccess=NULL);
template<typename _T>void Cholosky_Decompose(_T A[], int iOrder, _T B[]=NULL,int* pbSuccess=NULL);
void Conjugate_Gradient(float* A, const int n, float B[], float X[]);
void Get_Inv_Matrix(float* pM, float* pInv, int iOrder, int* pbSuccess);
void Solve_Linear_Cramer(float* A, int iOrder, float* B, float* X, int* pbSuccess); //����ķ��������Է�����
template<typename _T>void Test_Linear(_T A[], int iOrder, _T B[], _T X[],_T *pfError_Sum=NULL);

template<typename _T>void Solve_Linear_Gause_AAt(_T* A, int iOrder, _T* B, _T* X, int* pbSuccess);
template<typename _T>void Solve_Linear_Gause(_T* A, int iOrder, _T* B, _T* X, int* pbSuccess);	//��˹��������Է�����
//void Solve_Linear_Gause(float* A, int iOrder, float* B, float* X, int* pbSuccess=NULL);
template<typename _T>void Solve_Linear_LLt(_T A[], int iOrder, _T* B, _T* X, int* pbSuccess=NULL);

template<typename _T>void Solve_Linear_Gauss_Seidel(_T A[], _T B[], int iOrder, _T X[], int* pbResult);
template<typename _T>void Solve_Linear_Jocabi(_T A[], _T B[], int iOrder, _T X[], int* pbResult);				//�ſɱȵ�����������Է�����
void Solve_Linear_Gauss_Seidel(float A[], float B[], int iOrder, float X[], int* pbResult);
template<typename _T> void Solve_Linear_Contradictory(_T* A, int m, int n, _T* B, _T* X, int* pbSuccess);
//void Solve_Linear_Contradictory(float* A, int m, int n, float* B, float* X, int* pbSuccess);	//��ì�ܷ����飬����С���˽�
//void Solve_Linear_Solution_Construction(float* A, int m, int n, float B[], int* pbSuccess, float* pBasic_Solution = NULL, int* piBasic_Solution_Count = NULL, float* pSpecial_Solution = NULL);
template<typename _T> void Solve_Linear_Solution_Construction(_T* A, int m, int n, _T B[], int* pbSuccess, _T* pBasic_Solution=NULL, int* piBasic_Solution_Count=NULL, _T* pSpecial_Solution=NULL);
template<typename _T> void Solve_Linear_Solution_Construction_1(_T* A, int m, int n, _T B[], int* pbSuccess, _T* pBasic_Solution = NULL, int* piBasic_Solution_Count = NULL, _T* pSpecial_Solution = NULL);

void Get_Linear_Solution_Construction(float A[], int m, int n, float B[]);//����Ľṹ

//�������Է��̽�
template<typename _T>_T fLinear_Equation_Check(_T* A, int iOrder, _T* B, _T* X, _T eps=1e-9);

template<typename _T>void Elementary_Row_Operation(_T A[], int m, int n, _T A_1[], int* piRank=NULL, _T** ppBasic_Solution=NULL, _T** ppSpecial_Solution=NULL);
template<typename _T>void Elementary_Row_Operation_1(_T A[], int m, int n, _T A_1[], int* piRank=NULL, _T** ppBasic_Solution=NULL, _T** ppSpecial_Solution=NULL);

int bIs_Linearly_Dependent(float A[], int m, int n);	//�ж�������A�Ƿ��������
template<typename _T>int bIs_Orthogonal(_T* A, int h, int w=0);
template<typename _T>int bIs_R(_T R[9], _T eps = 1e-6);
template<typename _T>void Gen_I_Matrix(_T M[], int h, int w);	//���ɶԽ���

void Gen_Iterate_Matrix(float A[], int n, float B[], float B_1[], float C[]);	//���ɵ����е�B
int bIs_Contrative_Mapping(float B[], int n);					//�ж�һ�������Ƿ�Ϊѹ����������һ�����Է������Ƿ����õ��������ı�Ҫ����
int bIs_Diagonal_Dominant(float A[], int iOrder);			//�ж�һ�������Ƿ�Ϊ�Խ���ռ�Ŷ�������ռ�ž���
template<typename _T>int bIs_Symmetric_Matrix(_T A[], int iOrder, const _T eps = (_T)1e-4);
template<typename _T>int bIs_Unit_Vector(_T* V, int na, _T eps = ZERO_APPROCIATE);//���������Ƿ�Ϊ��λ����
void Diagonalize(float A[], int iOrder, float Diag[], float P[], int* pbSuccess);	//�Գƾ���Խǻ�

void Solve_Cubic(double A[3][3], Complex_d Root[3], int* piRoot_Count);
void Solve_Eigen(double A[3][3], Complex_d Root[3], int* piRoot_Count);
void Solve_Eigen_3x3(float A[], Complex_f Root[3]);
void Solve_Eigen_Vector(float A[], int iOrder, float fEigen_Value, float Q[], int* piCount);	//���ݸ���������ֵ����A�� ����������̵õ���������

void Solve_Poly(float Coeff[], int iCoeff_Count, Complex_f Root[], int* piRoot_Count = NULL);

//**************************ϡ�����*******************************
#define Get_Item(oMatrix, x1, y1, poItem) \
{ \
	int iCur; \
	iCur = (oMatrix).m_pRow[y1]; \
	if (!iCur) \
		poItem = NULL; \
	else \
	{\
		poItem = &(oMatrix).m_pBuffer[iCur];\
		while (poItem)\
		{\
			if (poItem->x == (unsigned int)x1)\
				break;\
			if (poItem->m_iRow_Next && poItem->x < (unsigned int)x1)\
				poItem = &(oMatrix).m_pBuffer[poItem->m_iRow_Next];\
			else\
			{\
				poItem = NULL;\
				break;\
			}\
		}\
	}\
}

template<typename _T> void Build_Link_Col_1(Sparse_Matrix<_T> oC)
{
	/*������һ�����Ͻ���Col Link*/
	int iItem_Index;
	typename Sparse_Matrix<_T>::Item* poItem;
	memset(oC.m_pCol, 0, oC.m_iCol_Count * sizeof(int));
	for (int i = oC.m_iRow_Count - 1; i >= 0; i--)
	{
		//if (i == 238)
		//printf("here");
		if (oC.m_pRow[i])
		{
			iItem_Index = oC.m_pRow[i];
			poItem = &oC.m_pBuffer[iItem_Index];
			while (1)
			{
				if (oC.m_pCol[poItem->x])/*��column��һ���ж���*/
					poItem->m_iCol_Next = oC.m_pCol[poItem->x];
				oC.m_pCol[poItem->x] = iItem_Index;
				if (poItem->m_iRow_Next)
				{
					iItem_Index = poItem->m_iRow_Next;
					poItem = &oC.m_pBuffer[poItem->m_iRow_Next];
				}
				else
					break;
			}
		}
	}
}

#define Build_Link_Col(oC) \
{ \
	/*������һ�����Ͻ���Col Link*/ \
	int iItem_Index; \
	typename Sparse_Matrix<_T>::Item* poItem; \
	memset(oC.m_pCol,0,oC.m_iCol_Count*sizeof(int)); \
	for (int i = oC.m_iRow_Count - 1; i >= 0; i--) \
	{ \
		if (oC.m_pRow[i]) \
		{ \
			iItem_Index = oC.m_pRow[i]; \
			poItem = &oC.m_pBuffer[iItem_Index]; \
			while (1) \
			{ \
				if (oC.m_pCol[poItem->x])/*��column��һ���ж���*/ \
					poItem->m_iCol_Next = oC.m_pCol[poItem->x]; \
				oC.m_pCol[poItem->x] = iItem_Index; \
				if (poItem->m_iRow_Next) \
				{ \
					iItem_Index = poItem->m_iRow_Next; \
					poItem = &oC.m_pBuffer[poItem->m_iRow_Next]; \
				}else \
					break; \
			} \
		} \
	} \
}
template<typename _T>int bIs_Symmetric_Matrix(Sparse_Matrix<_T> oA);
template<typename _T> _T fGet_Value(Sparse_Matrix<_T>* poMatrix, int x, int y);
template<typename _T>void Init_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix, int iItem_Count, int w, int h);
//template<typename _T>void Init_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix, int iItem_Count, int iMax_Order);
template<typename _T>void Compact_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix);
template<typename _T>void Free_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix);
template<typename _T>void Reset_Sparse_Matrix(Sparse_Matrix<_T>* poA);
template<typename _T>void Set_Value(Sparse_Matrix<_T>* poMatrix, int x, int y, _T fValue);
template<typename _T>void Matrix_Multiply(Sparse_Matrix<_T> A, Sparse_Matrix<_T> B, Sparse_Matrix<_T>* poC, int bCompact=1);
template<typename _T>void Matrix_Multiply(Sparse_Matrix<_T> A, _T a);
template<typename _T>void Resize_Matrix(Sparse_Matrix<_T>* poA, int iNew_Item_Count);
template<typename _T>void Matrix_Transpose_1(Sparse_Matrix<_T> A, Sparse_Matrix<_T>* poAt);
template<typename _T>void Matrix_Minus(Sparse_Matrix<_T> oA, Sparse_Matrix<_T> oB, Sparse_Matrix<_T>* poC, int bCompact=1);
template<typename _T>void Matrix_Add(Sparse_Matrix<_T> oA, Sparse_Matrix<_T> oB, Sparse_Matrix<_T>* poC, int bCompact=1);
template<typename _T>void Sparse_2_Dense(Sparse_Matrix<_T> oA, _T B[]);
template<typename _T>void Dense_2_Sparse(_T A[], int m, int n, Sparse_Matrix<_T>* poA, int bCompact=1);
template<typename _T>void Add_I_Matrix(Sparse_Matrix<_T>* poA, int* pbSuccess = NULL, _T ramda=1.f);
template<typename _T>void Solve_Linear_Gause(Sparse_Matrix<_T> oA, _T B[], _T X[], int* pbSuccess = NULL);
template<typename _T>void Solve_Linear_Gause_1(Sparse_Matrix<_T>oA1, _T B[], _T X[], int* pbSuccess = NULL);
template<typename _T>void Solve_Linear_Gause_2(Sparse_Matrix<_T>oA1, _T B[], _T X[], int* pbSuccess = NULL);
template<typename _T>void Cholosky_Decompose(Sparse_Matrix<_T> oA, Sparse_Matrix<_T> *poB,int *pbSuccess=NULL);

template<typename _T>void Re_Arrange_Sparse_Matrix(Sparse_Matrix<_T>* poA);
template<typename _T>void Get_Inv_Matrix_Row_Op(Sparse_Matrix<_T> oA, Sparse_Matrix<_T>* poA_Inv, int* pbSuccess=NULL);
//**************************ϡ�����*******************************

template<typename _T>void Vector_Add(_T A[], _T B[], int n, _T C[]);
template<typename _T>void Vector_Minus(_T A[], _T B[], int n, _T C[]);
template<typename _T>void Vector_Multiply(_T A[], int n,_T a,  _T B[]);
template<typename _T> void Matrix_Multiply(_T* A, int ma, int na, _T a, _T* C);
template<typename _T>void Matrix_Multiply_3x1(_T A[3 * 3], _T B[3], _T C[3]);

template<typename _T> void Matrix_Multiply(_T* A, int ma, int na, _T* B, int nb, _T* C);
template<typename _T> void Transpose_Multiply(_T A[], int m, int n, _T B[], int bAAt=1);
template<typename _T>void Matrix_Add(_T A[], _T B[], int iOrder, _T C[]);
template<typename _T>void Matrix_Minus(_T A[], _T B[], int iOrder, _T C[]);
template<typename _T>void Crop_Matrix(_T Source[], int m, int n, int x, int y, int w, int h, _T Dest[], int iDest_Stride=0);

void Matrix_x_Vector(double A[3][3], double X[3], double Y[3]);
void QR_Decompose(double A1[3][3], double Q[3][3], double R[3][3]);
void QR_Decompose(float* A, int ma, int na, float* Q, float* R);
void QR_Decompose(float* A, int na, float* R, float* Q, int* pbSuccess = NULL, int* pbDup_Root = NULL);

template<typename _T>int bInverse_Power(_T A[], int n, _T* pfEigen_Value, _T Eigen_Vector[], _T eps = 1e-9f);
template<typename _T>int bInverse_Power(_T A[3][3], _T* pfEigen_Value, _T Eigen_Vector[3], _T eps = 1e-9f);

template<typename _T>void LU_Decompose(_T A[], _T L[], _T U[], int n, int* pbSuccess=NULL);
template<typename _T>int bTest_LU(_T* A, _T* L, _T* U, int n);
template<typename _T>void Power_Method(_T A[], int n, _T Eigen_Vector[], _T* pfEigen_Value);
template<typename _T>int bTest_Eigen(_T A[], int n, _T fEigen_Value, _T Eigen_Vector[], _T eps = 1e-7);

template<typename _T> void Normalize(_T V[], int n, _T V_1[]);
template<typename _T>void Softmax(_T A[], int n, _T B[]);
template<typename _T>void Softmax(_T A[], int m, int n, int iPadding, _T B[]);
template<typename _T>void Homo_Normalize(_T V0[], int n, _T V1[]);	//����Ǹ�ɵ�Ʒ�����������ƭtemplate

template<typename _T> void Get_Inv_Matrix_Row_Op(_T* pM, _T* pInv, int iOrder, int* pbSuccess=NULL);	//��������
template<typename _T>void Get_Inv_AAt_Row_Op(_T* pM, _T* pInv, int iOrder, int* pbSuccess=NULL);
template<typename _T>void Get_Inv_AAt_3x3(_T M[3*3], _T Inv[3*3], int* pbSuccess);

void Get_Inv_Matrix(float* pM, float* pInv, int iOrder, int* pbSuccess);	//�ð����������
template<typename _T>_T fGet_Determinant(_T* A, int iOrder);			//������ʽ

template<typename _T> _T fGet_Mod(_T V[], int n);//��������ģ

template<typename _T> _T fGet_Distance(_T V_1[], _T V_2[], int n);	//������������
float fGet_Theta(float v0[], float v1[], float Axis[], int n);		//����������֮��ļн�
int iRank(double A[3][3]);								//��3x3�������
template<typename _T>int iGet_Rank(_T* A, int m, int n);				//��һ�������ȣ��ó����б任��
template<typename _T>void Cross_Product(_T V0[], _T V1[], _T V2[]);	//�����
template<typename _T>_T fDot(_T V0[], _T V1[], int iDim);		//���ڻ������

//����Ϊһ�����任���������
template<typename _T> void Gen_Roation_Matrix_2D(_T Rotation_Center[2], _T theta, _T R[3 * 3]);
template<typename _T>void Gen_Rotation_Matrix_2D(_T R[2 * 2], float fTheta);

template<typename _T>void Gen_Rotation_Matrix(_T Axis[3], _T fTheta, _T R[]);//������ת������ת�Ƕ�����һ����ת���󣬴˴���������������һ·���
template<typename _T>void Gen_Translation_Matrix(_T Offset[3], _T T[]);	//ƽ�Ʊ任
void Gen_Scale_Matrix(float Scale[3], float T[]);			//�����任

//��Ԫ��,��ת������ת������������
template<typename _T>void Quaternion_2_Rotation_Matrix(_T Q[4], _T R[]);	//��Ԫ��ת��ת����
template<typename _T>void Quaternion_2_Rotation_Vector(_T Q[4], _T V[4]);	//��Ԫ��ת��ת����
void Quaternion_Add(float Q_1[], float Q_2[], float Q_3[]);	//��Ԫ����
void Quaternion_Minus(float Q_1[], float Q_2[], float Q_3[]);	//��Ԫ��������ʵ������������������ͨ�����Ӽ�
void Quaternion_Conj(float Q_1[], float Q_2[]);				//���������
void Quaternion_Multiply(float Q_1[], float Q_2[], float Q_3[]);	//��Ԫ���˷�
void Quaternion_Inv(float Q_1[], float Q_2[]);	//��Ԫ������
template<typename _T>void Rotation_Matrix_2_Quaternion(_T R[], _T Q[]);		//��ת����ת����Ԫ��
template<typename _T>void Rotation_Matrix_2_Vector(_T R[3*3], _T V[4]);		//��ת����ת��ת����
template<typename _T>void Rotation_Matrix_2_Vector_3(_T R[3 * 3], _T V[3]);	
template<typename _T>void Rotation_Matrix_2_Vector_4(_T R[3 * 3], _T V[4]);

template<typename _T>void Rotation_Vector_4_2_Matrix(_T V[4], _T R[3 * 3]);		//��ת����ת����ת����
template<typename _T>void Rotation_Vector_3_2_Matrix(_T V[3], _T R[3 * 3]);		//��ת����ת����ת����
template<typename _T>void Rotation_Vector_2_Quaternion(_T V[4], _T Q[4]);	//��ת����ת����Ԫ��
template<typename _T>void Rotation_Vector_3_2_4(_T V[], _T V1[]);	
template<typename _T>void Rotation_Vector_4_2_3(_T V[], _T V1[]);
template<typename _T>void Hat(_T V[], _T M[]);				//���췴�Գƾ���
template<typename _T>void Vee(_T M[], _T V[3]);				//���Գƾ���ת��Ϊ����
template<typename _T>void Lie_Bracket(_T a[3], _T b[3], _T c[3]);	//������

template<typename _T>void se3_2_SE3(_T Ksi[6], _T T[]);		//��se3������SE3�����ת��
template<typename _T>void SE3_2_se3(_T T[4 * 4], _T ksi[6]);
template<typename _T>void SE3_2_se3(_T Rotation_Vector[4], _T t[3], _T Ksi[6]);	//һ����ת������һ��λ�����������һ��Ksi

void SIM3_2_sim3(float Rotation_Vector[], float t[], float s, float zeta[7]);
void sim3_2_SIM3(float zeta[7], float Rotation_Vector[4], float t[], float* ps);
template<typename _T>void Get_Jl_4(_T Rotation_Vector[4], _T J[]);
template<typename _T>void Get_Jl_3(_T Rotation_Vector[3], _T J[]);
template<typename _T>void Get_Jr_4(_T Rotation_Vector[4], _T J[]);
template<typename _T>void Get_Jr_3(_T Rotation_Vector[3], _T J[]);

template<typename _T> void Gen_Ksi_by_Rotation_Vector_t(_T Rotation_Vector[4], _T t[3], _T Ksi[6], int n=4);

template<typename _T> void Exp_Ref(_T A[], int n, _T B[],_T eps = (_T)1e-10);
template<typename _T>void Gen_Homo_Matrix_2D(_T R[2 * 2], _T t[2], _T T[3 * 3]);
template<typename _T>void Gen_Homo_Matrix(_T R[3*3], _T t[3], _T c2w[3*3]);			//����ת������λ�����깹��һ��SE3�任����
template<typename _T>void Gen_Homo_Matrix_1(_T Rotation_Vector[3], _T t[3], _T c2w[]);

template<typename _T>void Get_R_t(_T T[4 * 4], _T R[3 * 3]=NULL, _T t[3]=NULL);
template<typename _T>void Matrix_2_R(_T Hr[3 * 3], _T R[3 * 3]);

void Gen_Homo_Matrix(float R[], float t[], float s, float M[]);	//����ת,λ�ƣ����Ź���һ���任����

//ȡ��͸�ӱ任��һ������
void Perspective(float Pos_Source[3], float h[3], float Pos_Dest[3]);
void Perspective_Camera(float Pos_Source[3], float h[3], float Pos_Dest[3]);

template<typename _T> void svd_3(_T* A,SVD_Info oInfo, int* pbSuccess=NULL, double eps= 2.2204460492503131e-15);
template<typename _T>void Test_SVD(_T A[], SVD_Info oSVD, int* piResult=NULL, double eps= 2.2204460492503131e-15);
void Free_SVD(SVD_Info* poInfo);

template<typename _T>void Test_Inv_Matrix(_T A[], _T A_Inv[], int iOrder, _T* pfError_Sum=NULL);

//ע�⣬�÷�Ϊ SVD_Alloc<_T>(...)
template<typename _T>void SVD_Alloc(int h, int w, SVD_Info* poInfo, _T* A=NULL);

//һ���Ԫ����ʽ����
void Disp(Polynormial oPoly, int iElem_No=-1);
void Free_Polynormial(Polynormial* poPoly);				//�ͷŶ���ʽ��ռ�ڴ�
void Init_Polynormial(Polynormial* poPoly, int iElem_Count, int iMax_Term_Count);	//��ʼ��
void Add_Poly_Term(Polynormial* poPoly, float fCoeff, int first, ...);
void Get_Derivation(Polynormial* poSource, int iElem, Polynormial* poDest=NULL);	//����ʽ��xi��
float fGet_Polynormial_Value(Polynormial oPoly, float x[]);	//����x�����ʽֵ

template<typename _T>void Copy_Matrix_Partial(_T Source[], int iSource_Stride, int x, int y, _T Dest[], int m, int n);
template<typename _T>void Copy_Matrix_Partial(_T Source[], int m, int n, _T Dest[], int iDest_Stride, int x, int y);
template<typename _T>void Matrix_Add_Partial(_T Source[], int iSource_Stride, int x1, int y1, int w, int h, _T Dest[], int iDest_Stride, int x2, int y2);

/**********************һЩ������λ�û���*******************/
int bGet_Matrix_xy(int iPos, int n, int* px, int* py);
int iGet_Upper_Pos(int x, int y, int n);
/**********************һЩ������λ�û���*******************/

float fGet_e();
template<typename _T>_T fGet_Tr(_T M[], int iOrder);