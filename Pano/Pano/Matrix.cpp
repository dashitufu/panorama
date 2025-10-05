#include "stdio.h"
#include <iostream>
#include "assert.h"
#include "memory.h"
#include "math.h"
#include <stdarg.h>
#include "Common.h"
#include "Matrix.h"
#include "float.h"

extern "C"
{
#include "Buddy_System.h"
}

using namespace std;

#define  Init_LU(L, U) \
{/*��ʼ��L(�����ǣ� U�������ǣ�*/ \
	L[0][0] = L[1][1] = L[2][2] = 1; \
	L[0][1] = L[0][2] = L[1][2] = 0; \
	U[1][0] = U[2][0] = U[2][1] = 0; \
}

//Light_Ptr oMem_Mgr = { 0 };	//����ڴ�ר����Matrix�������ݣ��뿪�������ͷ�
//��ʱֻ��Ӧ�Ե��߳�

unsigned long long iGet_Random_No_cv(unsigned long long* piState)
{
	*piState = (unsigned long long)(unsigned) * piState * 4164903690U + (unsigned)(*piState >> 32);
	return *piState;
}

unsigned long long iGet_Random_No_cv()
{//α���������Opencv����
	static unsigned long long iNum = 0xFFFFFFFFFFFFFFFF;
	iNum = (unsigned long long)(unsigned)iNum * 4164903690U + (unsigned)(iNum >> 32);
	return iNum;
}
int iGet_Random_No_cv(int a, int b)
{//ȡ���a,b֮��
	unsigned iNum = (unsigned)iGet_Random_No_cv();
	return a == b ? a : (int)(iNum % (b - a) + a);
}
template<typename _T> _T fGet_Mod(_T V[], int n)
{
	_T fMod;
	int i;
	for (fMod = 0, i = 0; i < n; i++)
		fMod += V[i] * V[i];
	return sqrt(fMod);
}
template<typename _T> _T fGet_Distance(_T V_1[], _T V_2[], int n)
{//���������ľ��룬��ƽ����ʾ��������
	int i;
	_T fSum;
	for (i = 0, fSum = 0; i < n; i++)
		fSum += (V_1[i] - V_2[i]) * (V_1[i] - V_2[i]);
	return fSum;
}
void Normalize_Col(float M[], int m, int n, float M_1[])
{//��M���н��й��
	int y, x;
	float fTotal, fMod;
	for (x = 0; x < n; x++)
	{
		//��һ�е�ģ
		fTotal = 0;
		for (y = 0; y < m; y++)
			fTotal += M[y * n + x] * M[y * n + x];
		fMod = sqrt(fTotal);
		if (fMod != 0)
			for (y = 0; y < m; y++)
				M_1[y * n + x] = M[y * n + x] / fMod;
		else
			M_1[y * n + x] = 0;
	}
	return;
}
template<typename _T> void Normalize(_T V[], int n, _T V_1[])
{//���������й��
//����ÿ���������ܹ�񻯣���������0����ʱ���涨����0����
	_T fMod;
	int i;
#define eps 1e-10
	fMod = (_T)fGet_Mod(V, n);
	if (fMod > eps)
	{
		for (i = 0; i < n; i++)
			V_1[i] = V[i] / fMod;
	}
	else
		memset(V_1, 0, n * sizeof(_T));
#undef eps
}
void Disp(float* M, int iHeight, int iWidth, const char* pcCaption)
{
	int i, j;
	if (pcCaption)
		printf("%s\n", pcCaption);
	for (i = 0; i < iHeight; i++)
	{
		for (j = 0; j < iWidth; j++)
			printf("%.8f, ", (float)M[i * iWidth + j]);
		printf("\n");
	}
	return;
}
void Disp(double* M, int iHeight, int iWidth, const char* pcCaption)
{
	int i, j;
	if (pcCaption)
		printf("%s\n", pcCaption);
	for (i = 0; i < iHeight; i++)
	{
		for (j = 0; j < iWidth; j++)
			printf("%.8f, ", (double)M[i * iWidth + j]);
		printf("\n");
	}
	return;
}

template<typename _T>void Disp_Part(_T* M, int iStride, int x, int y, int w, int h, const char* pcCaption)
{
	int i, j;
	if (w + x > iStride)
		return;
	if (pcCaption)
		printf("%s\n", (char*)pcCaption);

	for (i = 0; i < h; i++)
	{
		for (j = 0; j < w; j++)
		{
			int iPos = (y + i) * iStride + (x + j);

			if (typeid(_T) == typeid(float))
			{
				if (M[i * iStride + j] - (int)M[iPos])
					printf("%.10ef, ", (double)M[iPos]);
				else
					printf("%f, ", (float)M[iPos]);
			}
			else if (typeid(_T) == typeid(double))
			{
				if (M[i * iStride + j] - (int)M[iPos])
					printf("%.10ef, ", (double)M[iPos]);
				else
					printf("%d, ", (int)M[iPos]);
			}
			else if ( (typeid(_T) == typeid(unsigned int)) ||
				(typeid(_T) == typeid(int)) ||
				(typeid(_T) == typeid(unsigned short)) ||
				(typeid(_T) == typeid(unsigned char)))
				printf("%d,", (int)M[iPos]);
		}
		printf("\n");
	}
}

void Disp(double V[3], const char* pcCaption)
{
	printf("%s\n", pcCaption);
	printf("%f %f %f\n", V[0], V[1], V[2]);
}
void Disp(float V[3], const char* pcCaption)
{
	printf("%s\n", pcCaption);
	printf("%f %f %f\n", V[0], V[1], V[2]);
}
void Disp(float A[3][3], const char* pcCaption = "")
{//��ʾһ������
	int y, x;
	printf("%s\n", pcCaption);
	for (y = 0; y < 3; y++)
	{
		for (x = 0; x < 3; x++)
			printf("%f\t", A[y][x]);
		printf("\n");
	}
	return;
}
void Disp(double A[3][3], const char* pcCaption)
{//��ʾһ������
	int y, x;
	printf("%s\n", pcCaption);
	for (y = 0; y < 3; y++)
	{
		for (x = 0; x < 3; x++)
			printf("%f\t", A[y][x]);
		printf("\n");
	}
	return;
}

void Matrix_x_Vector(double A[3][3], double X[3], double Y[3])
{//����Y=AX
	Y[0] = A[0][0] * X[0] + A[0][1] * X[1] + A[0][2] * X[2];
	Y[1] = A[1][0] * X[0] + A[1][1] * X[1] + A[1][2] * X[2];
	Y[2] = A[2][0] * X[0] + A[2][1] * X[1] + A[2][2] * X[2];
}
double fGet_Abs_Sum_3(double V[3])
{
	return Abs(V[0]) + Abs(V[1]) + Abs(V[2]);
}
float fGet_Next_Float(char** ppCur)
{
	char Value[256], * pCur = *ppCur;
	int i;
	for (i = 0; i < 256; i++, pCur++)
	{
		if (*pCur != '\t' && *pCur != '\0')
			Value[i] = *pCur;
		else
			break;
	}
	if (*pCur == '\0')
		*ppCur = pCur;
	else
		*ppCur = pCur + 1;
	return (float)atof(Value);
}
int bSave_Matrix(const char* pcFile, float* pMatrix, int iWidth, int iHeight)
{//����̣���\r\n����
	int i, j;
	FILE* pFile = fopen(pcFile, "wb");
	if (!pFile)
	{
		printf("Fail to open:%s\n", pcFile);
		return 0;
	}
	for (i = 0; i < iHeight; i++)
	{
		for (j = 0; j < iWidth; j++)
		{
			fprintf(pFile, "%f ", pMatrix[i * iWidth + j]);
		}
		fprintf(pFile, "\r\n");
	}
	fclose(pFile);
	return 1;
}
template<typename _T>void Cholosky_Decompose(Sparse_Matrix<_T> oA,Sparse_Matrix<_T> *poB, int* pbSuccess)
{	
	typename Sparse_Matrix<_T>::Item* poItem, * poCol_Item, * poRow_Item, * poPre, ** pCol_All;
	Sparse_Matrix<_T> oB;
	int x,y,iResult=1;
	_T* pCol_Value;
	_T d = 1.f, fValue;

	Init_Sparse_Matrix(&oB, oA.m_iCur_Item * 100, oA.m_iRow_Count, oA.m_iCol_Count);
	pCol_Value = (_T*)pMalloc(oB.m_iCol_Count * sizeof(_T));
	memset(pCol_Value, 0, oB.m_iCol_Count * sizeof(_T));
	pCol_All = (typename Sparse_Matrix<_T>::Item**)pMalloc(&oMem_Mgr, oA.m_iCol_Count * 8);
	memset(pCol_All, 0, oA.m_iCol_Count * 8);

	//��һ�������Ǿ���
	for (y = 0; y < oA.m_iRow_Count; y++)
	{
		if (!oA.m_pRow[y])
		{
			oB.m_pRow[y] = { 0 };
			continue;
		}
		//oRow_Index.
		poItem = &oA.m_pBuffer[oA.m_pRow[y]];
		if (poItem->x > (unsigned int)y)
		{//�������Ͳ��ǶԳƾ���
			printf("Not symmetric matrix in Cholosky_Decompose\n"); iResult = 0;
			goto END;
		}
		oB.m_pRow[y] = oB.m_iCur_Item;

		while (1)
		{//һ�мӽ�ȥ
			
			oB.m_pBuffer[oB.m_iCur_Item] = { poItem->x,poItem->y,poItem->m_fValue,oB.m_iCur_Item + 1 };
			oB.m_iCur_Item++;
			if (poItem->m_iRow_Next)
				poItem = &oA.m_pBuffer[poItem->m_iRow_Next];
			else
				break;
			if (poItem->x > (unsigned int)y)
				break;
		}
		oB.m_pBuffer[oB.m_iCur_Item-1].m_iRow_Next = 0;
	}
	Build_Link_Col_1(oB);
	//Disp(oB, "B");
	//��ʼ��
	for (y = 0; y < oB.m_iRow_Count; y++)
	{
		/*if (y == 1)
			printf("here");*/
		//�ȸ�Խ���Ԫ��
		if (!oB.m_pCol[y])
		{//����û��ֵ��ȫ0�����в����أ�
			continue;
		}
		poItem = &oB.m_pBuffer[oB.m_pCol[y]];
		if (!poItem)
		{//���У�
			printf("Err");
			continue;
		}
		d *= poItem->m_fValue;
		if (poItem->m_fValue < 0.f)
		{//ʵ�����ִ��ڶԳƵ������������
			iResult = 0;
			printf("Invalid value in Cholosky_Decompose\n");
			return;
		}
		fValue=poItem->m_fValue = sqrt(poItem->m_fValue);	//�Խ��߿������
		//ɨһ�У����о��Ǹղ������x�У��������������ĶԽ���Ԫ��
		if (!poItem->m_iCol_Next)
			continue;
		else
			poCol_Item = poItem = &oB.m_pBuffer[poItem->m_iCol_Next];

		//�˴������ǹ��ҽ�����������ĶԽ���Ԫ��Ϊ��Ԫ		
		//ɨ��Ԫ�У��������������ĶԽ���Ԫ��
		while (1)
		{
			pCol_Value[poCol_Item->y]=poCol_Item->m_fValue /= fValue;
			if (poCol_Item->m_iCol_Next)
				poCol_Item = &oB.m_pBuffer[poCol_Item->m_iCol_Next];
			else
				break;
		}
		
		//����pCol_All
		for (int i = y + 1; i < oB.m_iCol_Count; i++)
			pCol_All[i] = &oB.m_pBuffer[oB.m_pCol[i]];

		//�Ը���Ԫ���ұߵ�Ԫ�ؼ��㣬����ÿһ��Ԫ�أ�����һ��������x��������y
		//��Ԫ�ؼ�ȥ ��Ԫ�еĵ�i��Ԫ�س��Ե�j��Ԫ��
		
		poCol_Item = poItem;
		while (1)
		{
			if (poCol_Item->m_iRow_Next && poCol_Item->m_fValue!=0.f)
			{
				//�ٶ���һ������ֵ����һ��ÿ��Ԫ�ض�Ҫ�жϣ����鷳
				poPre = poCol_Item;
				fValue = poCol_Item->m_fValue;
				poRow_Item = &oB.m_pBuffer[poCol_Item->m_iRow_Next];
				for (x = poPre->x+1; (unsigned int)x <=poCol_Item->y; x++)
				{					
					if (!pCol_Value[x])
						continue;	//����Ҫ�仯
					if (poRow_Item->x != x)
					{//Ҫ��Ԫ��
						if (oB.m_iCur_Item >= oB.m_iMax_Item_Count)
						{
							printf("Insufficient space in Cholosky_Decompose\n");;
							iResult = 0;
							goto END;
						}
						oB.m_pBuffer[oB.m_iCur_Item] = { (unsigned int)x,poRow_Item->y, -fValue * pCol_Value[x], (int)(poRow_Item-oB.m_pBuffer) };
						if (x == 997)
							printf("Here");
						Add_Col_Link(oB, pCol_All, &oB.m_pBuffer[oB.m_iCur_Item]);
						poPre->m_iRow_Next = oB.m_iCur_Item;
						poPre = &oB.m_pBuffer[oB.m_iCur_Item];
						oB.m_iCur_Item++;
					}else
					{
						poRow_Item->m_fValue -= fValue * pCol_Value[poRow_Item->x];
						if (poRow_Item->m_iRow_Next)
							poRow_Item = &oB.m_pBuffer[poRow_Item->m_iRow_Next];
						poPre = poRow_Item;
					}
				}
			}
			if (poCol_Item->m_iCol_Next)
				poCol_Item = &oB.m_pBuffer[poCol_Item->m_iCol_Next];
			else
				break;
		}

		//����pCol_Value=0
		poCol_Item = &oB.m_pBuffer[oB.m_pCol[y]];
		while (poCol_Item->m_iCol_Next)
		{
			poCol_Item = &oB.m_pBuffer[poCol_Item->m_iCol_Next];
			pCol_Value[poCol_Item->y] = 0;
		}
		//Disp(oB, "oB");
	}
	Compact_Sparse_Matrix(&oB);
END:
	if (pCol_Value)
		Free(pCol_Value);
	if (pCol_All)
		Free(pCol_All);
	if (pbSuccess)
		*pbSuccess = iResult;
	if (poB)
		*poB = oB;
	return;
}
template<typename _T>void LLt_Decompose(_T A[], int iOrder, _T B[], int* pbSuccess)
{//������������LLt�ֽ⣬����Ҫ��ȫ��������һ������Ҫ��
//С��Χ����û����
	int i, j, k, iPos_kk, iPos_ik, iPos_ij, iPos_jk;
	_T* B1 = (_T*)pMalloc(&oMem_Mgr, iOrder * iOrder * sizeof(_T));
	_T d = 1;
	memcpy(B1, A, iOrder * iOrder * sizeof(_T));
	for (k = 0; k < iOrder; k++)
	{
		iPos_kk = k * iOrder + k;
		d *= B1[iPos_kk];
		if (B1[iPos_kk] < 0.f)
		{//ʵ�����ִ��ڶԳƵ������������
			if (pbSuccess)
				*pbSuccess = 0;
			printf("Invalid value in Cholosky_Decompose\n");
			return;
		}
		B1[iPos_kk] = sqrt(B1[iPos_kk]);	//�Խ��ߵ����俪�����Զ��׼�
		for (i = k + 1; i < iOrder; i++)
		{//ɨһ�У����о��Ǹղ������kk�У��������������ĶԽ���Ԫ��
			iPos_ik = i * iOrder + k;	//�ƺ��ڸ������ǵ�k��
			//���÷����
			if (B1[iPos_ik] == 0)
			{//�����ǿ��Բ����������Ԫ�Ƿ�Ϊ0�ģ����Լ����ֽ�
			}
			else if (B1[iPos_ik] != 0.f)
				B1[iPos_ik] /= B1[iPos_kk];
			else
			{
				if (pbSuccess)
					*pbSuccess = 0;
				printf("Divided by 0 in Cholosky_Decompose\n");
				return;
			}
		}

		for (j = k + 1; j < iOrder; j++)
		{//�Ը����ұߵ�Ԫ�ؼ���
			iPos_jk = j * iOrder + k;
			for (i = j; i < iOrder; i++)
			{
				iPos_ij = i * iOrder + j;
				//if (B1[iPos_ij] != 0.f)
				B1[iPos_ij] -= B1[i * iOrder + k] * B1[iPos_jk];
				/*else
				{
					if (pbSuccess)
						*pbSuccess = 0;
					printf("Divided by 0 in Cholosky_Decompose\n");
					return;
				}*/
			}
		}
		//Disp(B1, iOrder, iOrder, "B1");
	}

	//��Ȼ������������0
	if (B)
	{
		for (i = 0; i < iOrder; i++)
			for (j = i + 1; j < iOrder; j++)
				B1[i * iOrder + j] = 0;
		memcpy(B, B1, iOrder * iOrder * sizeof(_T));
	}
	//printf("d:%f\n", d);
	if (pbSuccess)
		*pbSuccess = 1;
	if (B1)
		Free(B1);
	return;
}

template<typename _T>void Cholosky_Decompose(_T A[], int iOrder, _T B[],int *pbSuccess)
{//������������Cholosky�ֽ⣬������ڰ�B A=BxB'
//С��Χ����û����
//ע�⣬A����Գ��������Գ�δ��һ��������������ʱ�ֽ�ʧ��
	int i, j, k, iPos_kk, iPos_ik, iPos_ij, iPos_jk;
	_T* B1 = (_T*)pMalloc(&oMem_Mgr, iOrder * iOrder * sizeof(_T));
	_T d = 1;
	memcpy(B1, A, iOrder * iOrder * sizeof(_T));
	for (k = 0; k < iOrder; k++)
	{
		iPos_kk = k * iOrder + k;
		d *= B1[iPos_kk];
		if (B1[iPos_kk] <= 0.f)
		{//ʵ�����ִ��ڶԳƵ������������
			if (pbSuccess)
				*pbSuccess = 0;
			printf("Not positive-definite Matrix in Cholosky_Decompose \n");
			return;
		}
		B1[iPos_kk] = sqrt(B1[iPos_kk]);	//�Խ��ߵ����俪�����Զ��׼�
		for (i = k + 1; i < iOrder; i++)
		{//ɨһ�У����о��Ǹղ������kk�У��������������ĶԽ���Ԫ��
			iPos_ik = i * iOrder + k;	//�ƺ��ڸ������ǵ�k��
			//���÷����
			//if (B1[iPos_ik] != 0.f)
				B1[iPos_ik] /= B1[iPos_kk];
			/*else
			{
				if(pbSuccess)
					*pbSuccess = 0;
				printf("Divided by 0 in Cholosky_Decompose \n");
				return;
			}*/
		}
		
		for (j = k + 1; j < iOrder; j++)
		{//�Ը����ұߵ�Ԫ�ؼ���
			iPos_jk = j * iOrder + k;
			for (i = j; i < iOrder; i++)
			{
				iPos_ij = i * iOrder + j;
				//if (B1[iPos_ij] != 0.f)
					B1[iPos_ij] -= B1[i * iOrder + k] * B1[iPos_jk];
				/*else
				{
					if (pbSuccess)
						*pbSuccess = 0;
					printf("Divided by 0 in Cholosky_Decompose\n");
					return;
				}*/
			}
		}
		//Disp(B1, iOrder, iOrder, "B1");
	}

	//��Ȼ������������0
	if (B)
	{
		for (i = 0; i < iOrder; i++)
			for (j = i + 1; j < iOrder; j++)
				B1[i * iOrder + j] = 0;
		memcpy(B, B1, iOrder * iOrder * sizeof(_T));
	}
	//printf("d:%f\n", d);
	if (pbSuccess)
		*pbSuccess = 1;
	if (B1)
		Free(B1);
	return;
}
void QR_Decompose(double A1[3][3], double Q[3][3], double R[3][3])
{//����һ��3x3����A������������Householder�����Q,R��R��������������Ǿ���
	//ע�⣺�п������������ݣ�����������������������ʵ���������²�����������
	double Q1[3][3], A2[3][3], Q2[3][3];
	double /*fMod,*/ fDelta, tau, v[3];
	//int x, y,i;

	//��һ������Q1����ʱ�������һ�����㣬��ȫ����������,��ʱ������I
	//�˴������Խ�һ��������ʱ����
	fDelta = sign(A1[0][0]) * sqrt(A1[0][0] * A1[0][0] + A1[1][0] * A1[1][0] + A1[2][0] * A1[2][0]);
	//v=(a11+delta,a21,a31...)
	v[0] = A1[0][0] + fDelta, v[1] = A1[1][0], v[2] = A1[2][0];
	tau = fDelta * (A1[0][0] + fDelta);
	tau = 1 / tau;	//�ܿ�������tau��Ϊ����
	//Q1=H1=I-vv'/tau
	Q1[0][0] = 1 - v[0] * v[0] * tau;
	Q1[1][0] = Q1[0][1] = -v[0] * v[1] * tau;
	Q1[2][0] = Q1[0][2] = -v[0] * v[2] * tau;
	Q1[1][1] = 1 - v[1] * v[1] * tau;
	Q1[2][1] = Q1[1][2] = -v[1] * v[2] * tau;
	Q1[2][2] = 1 - v[2] * v[2] * tau;

	//Disp(Q1,(char*)"Q1");
	//��A2=Q1*A
	A2[0][0] = Q1[0][0] * A1[0][0] + Q1[0][1] * A1[1][0] + Q1[0][2] * A1[2][0];
	A2[0][1] = Q1[0][0] * A1[0][1] + Q1[0][1] * A1[1][1] + Q1[0][2] * A1[2][1];
	A2[0][2] = Q1[0][0] * A1[0][2] + Q1[0][1] * A1[1][2] + Q1[0][2] * A1[2][2];
	A2[1][0] = A2[2][0] = 0;
	A2[1][1] = Q1[1][0] * A1[0][1] + Q1[1][1] * A1[1][1] + Q1[1][2] * A1[2][1];
	A2[1][2] = Q1[1][0] * A1[0][2] + Q1[1][1] * A1[1][2] + Q1[1][2] * A1[2][2];
	A2[2][1] = Q1[2][0] * A1[0][1] + Q1[2][1] * A1[1][1] + Q1[2][2] * A1[2][1];
	A2[2][2] = Q1[2][0] * A1[0][2] + Q1[2][1] * A1[1][2] + Q1[2][2] * A1[2][2];
	//Disp(A2, (char*)"A2");

	//�ڶ������˴��Ѿ�֤����Q[2][2]=-Q[1][1]��Ч��ţ�ƣ�
	//fDelta = sign(A2[1][1]) * sqrt(A2[1][1] * A2[1][1] + A2[2][1] * A2[2][1]);
	//v[0] = A2[1][1] + fDelta, v[1] = A2[2][1];
	//tau = 1/(fDelta * (A2[1][1] + fDelta));
	////Q=	I1	0
	////		0	H2
	//Q2[0][0] = 1, Q2[0][1] = Q2[0][2] = Q2[1][0] = Q2[2][0] = 0;
	//Q2[1][1] = 1 - tau * v[0] * v[0];
	//Q2[2][1]=Q2[1][2] = -tau * v[0] * v[1];
	//Q2[2][2] = 1 - tau * v[1] * v[1];
	//Disp(Q2,(char*)"Q2");

	//����Ϊ�ڶ����·���
	fDelta = 1 / (sign(A2[1][1]) * sqrt(A2[1][1] * A2[1][1] + A2[2][1] * A2[2][1]));
	//Q=	I1	0
	//		0	H2
	Q2[0][0] = 1, Q2[0][1] = Q2[0][2] = Q2[1][0] = Q2[2][0] = 0;
	Q2[1][1] = -(Q2[2][2] = A2[1][1] * fDelta);
	Q2[2][1] = Q2[1][2] = -A2[2][1] * fDelta;
	//Disp(Q2, (char*)"Q2");

	//Disp(Q2,(char*)"Q2");
	//A3����R���˴�����Q2,A2����Ԫ���򣬲��ң��������ۣ�R[1][0]=R[2]0]=R[2][1]=0
	R[0][0] =/* Q2[0][0] **/ A2[0][0];	//+Q2[0][1] * A2[1][0] + Q2[0][2] * A2[2][0];
	R[0][1] = /*Q2[0][0] **/ A2[0][1] /*+ Q2[0][1] * A2[1][1] + Q2[0][2] * A2[2][1]*/;
	R[0][2] =/* Q2[0][0] **/ A2[0][2] /*+ Q2[0][1] * A2[1][2] + Q2[0][2] * A2[2][2]*/;
	R[1][0] = R[2][0] = R[2][1] = 0;
	R[1][1] = /*Q2[1][0] * A2[0][1] +*/ Q2[1][1] * A2[1][1] + Q2[1][2] * A2[2][1];
	R[1][2] = /*Q2[1][0] * A2[0][2] +*/ Q2[1][1] * A2[1][2] + Q2[1][2] * A2[2][2];
	R[2][2] = /*Q2[2][0] * A2[0][2] +*/ Q2[2][1] * A2[1][2] + Q2[2][2] * A2[2][2];

	//Disp(R,(char*)"R=A3");
	//Q=(Q2*Q1)'
	//�˴���������Q2�����Լ��ټ������
	Q[0][0] = Q1[0][0];
	Q[1][0] = Q1[0][1];
	Q[2][0] = Q1[0][2];
	Q[0][1] = Q2[1][1] * Q1[1][0] + Q2[1][2] * Q1[2][0];
	Q[1][1] = Q2[1][1] * Q1[1][1] + Q2[1][2] * Q1[2][1];
	Q[2][1] = Q2[1][1] * Q1[1][2] + Q2[1][2] * Q1[2][2];
	Q[0][2] = Q2[2][1] * Q1[1][0] + Q2[2][2] * Q1[2][0];
	Q[1][2] = Q2[2][1] * Q1[1][1] + Q2[2][2] * Q1[2][1];
	Q[2][2] = Q2[2][1] * Q1[1][2] + Q2[2][2] * Q1[2][2];

	//Disp(Q,(char*)"Q");

	////������㿴���Ƿ�A=Q*R
	//for(y=0;y<3;y++)
	//	for(x=0;x<3;x++)
	//		for(A2[y][x]=0,i=0;i<3;i++)
	//			A2[y][x]+= Q[y][i] * R[i][x];
	return;
}

template<typename _T>void Vector_Minus(_T A[], _T B[], int n, _T C[])
{
	for (int i = 0; i < n; i++)
		C[i] = A[i] - B[i];
}
template<typename _T>void Vector_Add(_T A[], _T B[], int n, _T C[])
{
	for (int i = 0; i < n; i++)
		C[i] = A[i] + B[i];
}
template<typename _T>void Vector_Multiply(_T A[],int n, _T a, _T B[])
{//B = aA
	for (int i = 0; i < n; i++)
		B[i] = A[i] * a;
}
template<typename _T>int bIs_Symmetric_Matrix(Sparse_Matrix<_T> oA)
{//����һ��ϡ������Ƿ�Գ�
	typename Sparse_Matrix<_T>::Item* poRow_Item, * poCol_Item;
	int i;
	if (oA.m_iCol_Count != oA.m_iRow_Count)
		return 0;

	for(i=0;i<oA.m_iRow_Count;i++)
	{
		poRow_Item = &oA.m_pBuffer[oA.m_pRow[i]];
		poCol_Item = &oA.m_pBuffer[oA.m_pCol[i]];
		while (1)
		{
			if (poRow_Item->x != poCol_Item->y || poRow_Item->y != poCol_Item->x || poRow_Item->m_fValue != poCol_Item->m_fValue)
				return 0;
			if (!!poRow_Item->m_iRow_Next ^ !!poCol_Item->m_iCol_Next)
				return 0;

			if (poRow_Item->m_iRow_Next)
				poRow_Item = &oA.m_pBuffer[poRow_Item->m_iRow_Next];
			else
				break;
			if (poCol_Item->m_iCol_Next)
				poCol_Item = &oA.m_pBuffer[poCol_Item->m_iCol_Next];
			else
				break;
		}
	}
	return 1;
}
template<typename _T>int bIs_Symmetric_Matrix(_T A[], int iOrder, const _T eps)
{
	int y, x,bRet=1;
	float fA, fB, fMin, fDiff;
	for (y = 0; y < iOrder; y++)
	{
		for (x = 0; x < iOrder; x++)
		{
			fA = (float)abs(A[y * iOrder + x]);
			fB = (float)abs(A[x * iOrder + y]);
			fMin = Min(fA, fB);
			fDiff = abs(fA - fB);
			if (fDiff / fMin > eps)
			{
				if(typeid(_T)== typeid(float))
					printf("%f %f\n", (float)A[y * iOrder + x], (float)A[x * iOrder + y]);
				else if(typeid(_T)==typeid(double))
					printf("%lf %lf\n", (double)A[y * iOrder + x], (double)A[x * iOrder + y]);
				bRet = 0;
			}
		}
	}
	return bRet;
}
void Matrix_Multiply_Symmetric(float* A, int m, int n, float* AAt)
{//�� AxA', �����Ȼ��һ���Գƾ���A(mxn) x A'(n*m) = AAt(m*m)
	float* AAt_1 = (float*)malloc(m * m * sizeof(float));
	int y, x, i;
	float fValue;
	for (y = 0; y < m; y++)
	{
		for (x = y; x < m; x++)
		{
			//if (y == 0 && x == 1)
				//printf("here");
			for (fValue = 0, i = 0; i < n; i++)
				fValue += A[y * n + i] * A[x * n + i];
			AAt_1[x * m + y] = AAt_1[y * m + x] = fValue;
		}
	}
	//Disp(AAt_1, m, m, "AAt");
	memcpy(AAt, AAt_1, m * m * sizeof(float));
	free(AAt_1);
	return;
}
//void Matrix_Multiply(float* A,int ma,int na,float *B, int nb,float *C)
//{//Amn x Bno = Cmo
//	int y, x,i;
//	float fValue,*C_Dup=(float*)malloc(ma*nb*sizeof(float));
//	for (y = 0; y < ma; y++)
//	{
//		for (x = 0; x < nb; x++)
//		{
//			//if (y == 2 && x == 1)
//				//printf("Here");
//			for (fValue = 0, i = 0; i < na; i++)
//				fValue += A[y * na + i] * B[i * nb + x];
//			C_Dup[y * nb + x] = fValue;
//		}
//	}
//	memcpy(C, C_Dup, ma * nb * sizeof(float));
//	free(C_Dup);
//	return;
//}

void RQ_Multiply_3x3(double R[3][3], double Q[3][3], double A[3][3])
{//����A=RxQ
	//�˴�Ӧ�ð���R[1][0]=R[2][0]=R[2][1]=0���л���
	A[0][0] = R[0][0] * Q[0][0] + R[0][1] * Q[1][0] + R[0][2] * Q[2][0];
	A[0][1] = R[0][0] * Q[0][1] + R[0][1] * Q[1][1] + R[0][2] * Q[2][1];
	A[0][2] = R[0][0] * Q[0][2] + R[0][1] * Q[1][2] + R[0][2] * Q[2][2];
	A[1][0] = /*R[1][0] * Q[0][0] +*/ R[1][1] * Q[1][0] + R[1][2] * Q[2][0];
	A[1][1] = /*R[1][0] * Q[0][1] +*/ R[1][1] * Q[1][1] + R[1][2] * Q[2][1];
	A[1][2] = /*R[1][0] * Q[0][2] +*/ R[1][1] * Q[1][2] + R[1][2] * Q[2][2];
	A[2][0] = /*R[2][0] * Q[0][0] + R[2][1] * Q[1][0] +*/ R[2][2] * Q[2][0];
	A[2][1] = /*R[2][0] * Q[0][1] + R[2][1] * Q[1][1] +*/ R[2][2] * Q[2][1];
	A[2][2] = /*R[2][0] * Q[0][2] + R[2][1] * Q[1][2] +*/ R[2][2] * Q[2][2];
}
#define Scale_Matrix(A, A1,fMax) \
{ \
	double fRecep; \
	fMax=Abs(A[0][0]);	\
	if (Abs(A[0][1]) > fMax)fMax = Abs(A[0][1]); \
	if (Abs(A[0][2]) > fMax)fMax = Abs(A[0][2]); \
	if (Abs(A[1][0]) > fMax)fMax = Abs(A[1][0]); \
	if (Abs(A[1][1]) > fMax)fMax = Abs(A[1][1]); \
	if (Abs(A[1][2]) > fMax)fMax = Abs(A[1][2]); \
	if (Abs(A[2][0]) > fMax)fMax = Abs(A[2][0]); \
	if (Abs(A[2][1]) > fMax)fMax = Abs(A[2][1]); \
	if (Abs(A[2][2]) > fMax)fMax = Abs(A[2][2]); \
	fRecep = 1 / fMax; \
	A1[0][0] = A[0][0] * fRecep; \
	A1[0][1] = A[0][1] * fRecep; \
	A1[0][2] = A[0][2] * fRecep; \
	A1[1][0] = A[1][0] * fRecep; \
	A1[1][1] = A[1][1] * fRecep; \
	A1[1][2] = A[1][2] * fRecep; \
	A1[2][0] = A[2][0] * fRecep; \
	A1[2][1] = A[2][1] * fRecep; \
	A1[2][2] = A[2][2] * fRecep; \
}
#define Scale_Matrix_1(A1,fMax) \
{ \
	A1[0][0] *= fMax; \
	A1[0][1] *= fMax; \
	A1[0][2] *= fMax; \
	A1[1][0] *= fMax; \
	A1[1][1] *= fMax; \
	A1[1][2] *= fMax; \
	A1[2][0] *= fMax; \
	A1[2][1] *= fMax; \
	A1[2][2] *= fMax; \
}
void Solve_Homo_3x3(double A[3][3], double fEigen_Value, double x[3])
{//Ϊ3x3���������η���Ax=0
	double fMod, Sub[2];	//substitution;
	double A1[3][3] = { {A[0][0] - fEigen_Value,A[0][1],A[0][2]},
						{A[1][0],A[1][1] - fEigen_Value,A[1][2]},
						{A[2][0],A[2][1],A[2][2] - fEigen_Value} };
	//Disp(A1, "A1"); 
	//��һ����ȥA[1][0],A[2][0]
	//�ҳ�x1������Ԫ,������ʽ
	if (A1[0][0] != 0)
	{
		Sub[0] = -(A1[0][1] /= A1[0][0]);
		Sub[1] = -(A1[0][2] /= A1[0][0]);
		A1[1][1] += A1[1][0] * Sub[0];
		A1[1][2] += A1[1][0] * Sub[1];
		A1[2][1] += A1[2][0] * Sub[0];
		A1[2][2] += A1[2][0] * Sub[1];
		A1[1][0] = A1[2][0] = 0;
		A1[0][0] = 1;
	}
	else if (A1[1][0] != 0)
	{
		Sub[0] = -(A1[1][1] /= A1[1][0]);
		Sub[1] = -(A1[1][2] /= A1[1][0]);
		A1[2][1] += A1[2][0] * Sub[0];
		A1[2][2] += A1[2][0] * Sub[1];
		//A1[1][0] = 1;
		A1[2][0] = 0;
		//����1�����0�н���λ��
		A1[0][0] = 1;
		A1[1][0] = 0;
		std::swap(A1[0][1], A1[1][1]);
		std::swap(A1[0][2], A1[1][2]);
	}
	else if (A1[2][0] != 0)
	{
		A1[2][1] /= A1[2][0];
		A1[2][2] /= A1[2][0];

		//����2�����0�н���λ��
		A1[0][0] = 1;
		A1[2][0] = 0;
		std::swap(A1[0][1], A1[2][1]);
		std::swap(A1[0][2], A1[2][2]);
	}
	//Disp(A1,"A1");
	//�ڶ�������ȥA[2][1],��ȥA[0][1]
	if (A1[1][1] != 0)
	{//A1[1][1]��Ϊ0������������2��
		A1[1][2] /= A1[1][1];
		A1[1][1] = 1;
		Sub[0] = -A1[1][2];
		A1[2][2] += A1[2][1] * Sub[0];
		A1[2][1] = 0;
	}
	else if (A1[2][1] != 0)
	{
		A1[2][2] /= A1[2][1];

		A1[1][1] = 1;
		A1[2][1] = 0;
		//����2�����1�н���
		std::swap(A1[1][2], A1[2][2]);
	}

	//Disp(A1, "A1");
	//����Ҫ����һ�������ʴ˴˴�����ֱ���ж�Ϊ0
	//if (A1[2][2] !=0)
	//{//���ȣ���·�˳�
	//	x[0] = x[1] = x[2] = 0;
	//	return;
	//}

	//��3����Ԫ
	if (A1[2][2] != 0)
	{
		A1[2][2] = 1;
		A1[0][2] = A1[1][2] = 0;
	}
	//Disp(A1, "A1");

	//�����һ����������x, �˴����������⣬��δ�㶨����������
	if (A1[2][2] == 1)
		x[2] = 0;
	else
		x[2] = 1;	//����Ԫ

	if (A1[1][1] == 1)
		x[1] = -A[1][2] * x[2];
	else
		x[1] = 1;

	if (A1[0][0] == 1)
		x[0] = -A[0][1] * x[1] - A[0][2] * x[2];
	else
		x[0] = 1;

	//if (A1[1][1] == 1)
	//{
	//	x[2] = 1;	//��x[2]��Ϊ����Ԫ
	//	x[1] = -A1[1][2] * x[2];
	//}

	//if (A1[0][0] == 1)
	//{
	//	if (A1[1][1] == 0)
	//	{//����������Ԫ
	//		x[1] = 1;
	//		x[2] = 1;
	//	}
	//	//x[0]=-a01*x1-a02*x2
	//	x[0] = -A1[0][1] * x[1] - A1[0][2] * x[2];
	//}

	////��ʱ�������� Ax=0
	//double Temp[3];
	/*double A2[3][3] = { {A[0][0] - fEigen_Value,A[0][1],A[0][2]},
						{A[1][0],A[1][1] - fEigen_Value,A[1][2]},
						{A[2][0],A[2][1],A[2][2] - fEigen_Value} };*/

						//x[2] = 1;//ǿ�����㣬��Ϊ����Ԫ
						//Matrix_x_Vector(A2, x, Temp);
						//fMod = fGet_Abs_Sum_3(Temp);
						//if (fMod >= 0.5)
							//printf("Err:%f\n",fMod);

						//���һ����Normalize
	fMod = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	if (fMod < ZERO_APPROCIATE)
		return;
	else
		fMod = 1 / sqrt(fMod);

	x[0] *= fMod;
	x[1] *= fMod;
	x[2] *= fMod;

	//Disp(A2, "A2");
	//����Ϊʲôȫ������-1��������ǰ��ʲô���ˣ��о����ﲻӦ�ã������ǳ�ʼ��һ������
	/*x[0] *= -1;
	x[1] *= -1;
	x[2] *= -1;*/
	//Disp(x, "x");

	//double Temp[3];
	//Matrix_x_Vector(A2, x, Temp);
	//fMod = fGet_Abs_Sum_3(Temp);
	//if (fMod >= 0.5)
		//printf("Err:%f\n", fMod);
	//Disp(x, "x");
	return;
}
int bEigen_Value_3x3(double A[3][3], double Ramda[3], double Eigen_Vector[3][3])
{//�����������������ֵ����A���죬����һ������ֵ=0
//����ɹ�������1�����򣬷���0
#define MAX_ITER_COUNT 120
	double Q[3][3], A1[3][3], fMax;
	static double R[3][3] = {};
	//ʵ��֤�����ֽ�ǰ����һ������/����������ٵ�������

	Scale_Matrix(A, A1, fMax);
	//Disp(A, "A");
	//Disp(A1, "A1");
	int iCount = 0;
	float /*fPre_Sum=MAX_FLOAT,*/fSum;

	//�˴���ʼ����
	while (1)
	{
		QR_Decompose(A1, Q, R);
		//Disp(R, (char*)"R="); Disp(Q, (char*)"Q=");
		//A1 = RQ
		RQ_Multiply_3x3(R, Q, A1);
		//Disp(A1, (char*)"RQ=");
		//�ж��Ƿ�������
		fSum = (float)(Abs(A1[1][0]) + Abs(A1[2][0]) + Abs(A1[2][1]));
		//printf("Sum:%f\n", fSum);
		if (fSum < 0.000001)
			break;
		else if (iCount >= MAX_ITER_COUNT)
		{
			//printf("Sum of Err:%f\n", fSum);
			return 0;
		}
		//if (Abs(A1[1][0]) < ZERO_APPROCIATE && Abs(A1[2][0]) < ZERO_APPROCIATE && Abs(A1[2][1]) < ZERO_APPROCIATE)
			//break;
		iCount++;
	}
	Ramda[0] = A1[0][0] * fMax;
	Ramda[1] = A1[1][1] * fMax;
	Ramda[2] = A1[2][2] * fMax;
	//Disp(A1, "A1");
	//Disp(Ramda, (char*)"Eigen Value=");
	//��ʱ���룬��A1
	//printf("%f %f %f\n", Ramda[0] * fMax, Ramda[1] * fMax, Ramda[2] * fMax);
	//Scale_Matrix(A, A1, fMax);
	//Solve_Homo_3x3(A, Ramda[0]*fMax, Eigen_Vector[0]);

	//Solve_Homo_3x3(A, Ramda[0]*=fMax, Eigen_Vector[0]);
	//Solve_Homo_3x3(A, Ramda[1] *= fMax, Eigen_Vector[1]);
	//Solve_Homo_3x3(A, Ramda[2] *= fMax, Eigen_Vector[2]);
	//Scale_Matrix_1(A1, fMax);	//ʵ���ϣ���һ�����п���ʡ�ԣ���ΪMatlab����������һ��
	//Disp(Eigen_Vector, "Eigen_Vector");
	return 1;
#undef MAX_ITER_COUNT
}
#define LU_Decompose_1(A, L, U, bSucceed) \
{ \
	if (A[0][0] == 0) \
	{ \
		bSucceed = 0; \
		goto END; \
	} \
	U[0][0] = A[0][0]; \
	L[1][0] = A[1][0] / U[0][0]; \
	U[0][1] = A[0][1]; \
	U[1][1] = A[1][1] - L[1][0] * U[0][1]; \
	if (U[1][1] == 0) \
	{ \
		bSucceed = 0; \
		goto END; \
	} \
	L[2][0] = A[2][0] / U[0][0]; \
	L[2][1] = (A[2][1] - L[2][0] * U[0][1]) / U[1][1]; \
	U[0][2] = A[0][2]; \
	U[1][2] = A[1][2] - L[1][0] * U[0][2]; \
	U[2][2] = A[2][2] - L[2][0] * U[0][2] - L[2][1] * U[1][2]; \
	bSucceed = 1; \
END: \
	; \
}

template<typename _T>int bTest_LU(_T* A, _T* L, _T* U, int n)
{
	int bRet = 1;
	_T* A1 = (_T*)pMalloc(n * n * sizeof(_T));
	Matrix_Multiply(L, n, n, U, n, A1);
	const _T eps = 1e-7f;
	for (int i = 0; i < n * n; i++)
	{
		if (abs(A1[i] - A[i]) > eps)
		{
			printf("A:%f A1:%f\n", A[i], A1[i]);
			bRet = 0;
		}
	}
	Free(A1);
	return bRet;
}
template<typename _T>void LU_Decompose(_T A[], _T L[], _T U[], int n,int *pbSuccess) 
{//һ��LU�ֽ⣬L�������ϵ�λ�����ǣ�����Ψһ��
	int i, j, k, r;
	int iResult = 1;
	memset(L, 0, n * n * sizeof(_T));
	memset(U, 0, n * n * sizeof(_T));
	for (i = 0; i < n; i++)
		L[i * n + i] = 1;

	for (i = 0; i < n; i++/*, iFlag_rc = ~iFlag_rc*/)
	{
		_T fValue;

		//�ȸ��i��,���¸�
		for (j = 0; j <= i; j++)
		{
			fValue =/* aji =*/ A[j * n + i];
			//if (j == 1 && i == 1)
				//printf("here");
			//�Դ�aij��Ϊ���������һ��uij
			for (k = 0; k < j; k++)
				fValue -= L[j * n + k] * U[k * n + i];
			U[j * n + i] = fValue;
		}

		//�ٸ��i+1��,���Ҹ�
		r = i + 1;
		if (r >= n)
			break;  //�����ˣ�����û��

		for (j = 0; j < r; j++)
		{//����ֻ��LԪ�أ�
			fValue = /*arj =*/ A[r * n + j];
			/*if (r == 2 && j == 1)
				printf("here");*/
			for (k = 0; k < j; k++)
				fValue -= L[r * n + k] * U[k * n + j];
			if (U[j * n + j] == 0)
			{//��Ҫ�˴����γɷֽ�����
				iResult = 0;
				goto END;
			}
			L[r * n + j] = fValue / U[j * n + j];
		}
	}
END:
	if (pbSuccess)
		*pbSuccess = iResult;
	return;
}

template<typename _T>void LU_Decompose(_T A[3][3], _T L[3][3], _T U[3][3], int* pbSucceed)
{//���Խ�һ������A�ֽ�ΪLxU=A���ֽ������λ��ɧ�������Կ����
//���ƣ����Է�����������ȣ���Ч
	if (A[0][0] == 0)
	{
		*pbSucceed = 0;
		return;
	}
	Init_LU(L , U);
	U[0][0] = A[0][0];
	L[1][0] = A[1][0] / U[0][0];
	U[0][1] = A[0][1];
	U[1][1] = A[1][1] - L[1][0] * U[0][1];
	if (U[1][1] == 0)
	{
		*pbSucceed = 0;
		return;
	}
	L[2][0] = A[2][0] / U[0][0];
	L[2][1] = (A[2][1] - L[2][0] * U[0][1]) / U[1][1];
	U[0][2] = A[0][2];
	U[1][2] = A[1][2] - L[1][0] * U[0][2];
	U[2][2] = A[2][2] - L[2][0] * U[0][2] - L[2][1] * U[1][2];
	/*double Temp[3][3];
	Matrix_Multiply_3x3(L, U, Temp);
	Disp(Temp, "Temp");
	Disp(L, (char*)"L");
	Disp(U, (char*)"U");*/
	if(pbSucceed)
		*pbSucceed = 1;
	return;
}

template<typename _T>void Solve_Uk(_T L[3][3], _T Z[3], _T U[3])
{//���LU=Z�ֵ�U
	U[0] = Z[0];
	U[1] = Z[1] - L[1][0] * U[0];
	U[2] = Z[2] - L[2][0] * U[0] - L[2][1] * U[1];
	return;
}
template<typename _T>void Solve_yk(_T U[3][3], _T yk[3], _T uk[3])
{//���Uyk=uk�е�yk
	yk[2] = uk[2] / U[2][2];
	yk[1] = (uk[1] - U[1][2] * yk[2]) / U[1][1];
	yk[0] = (uk[0] - U[0][1] * yk[1] - U[0][2] * yk[2]) / U[0][0];
	return;
}
int iRank(double A[3][3])
{//�Ծ�������ȣ�ע�⣬δ���ԣ�ǧ�������
	double A1[3][3], Sub[2];
	//��һ����ȥA[1][0],A[2][0]
	//�ҳ�x1������Ԫ,������ʽ
	memcpy(A1, A, 3 * 3 * sizeof(double));
	if (A1[0][0] != 0)
	{
		Sub[0] = -(A1[0][1] /= A1[0][0]);
		Sub[1] = -(A1[0][2] /= A1[0][0]);
		A1[1][1] += A1[1][0] * Sub[0];
		A1[1][2] += A1[1][0] * Sub[1];
		A1[2][1] += A1[2][0] * Sub[0];
		A1[2][2] += A1[2][0] * Sub[1];
		A1[1][0] = A1[2][0] = 0;
		A1[0][0] = 1;
	}
	else if (A1[1][0] != 0)
	{
		Sub[0] = -(A1[1][1] /= A1[1][0]);
		Sub[1] = -(A1[1][2] /= A1[1][0]);
		A1[2][1] += A1[2][0] * Sub[0];
		A1[2][2] += A1[2][0] * Sub[1];
		//A1[1][0] = 1;
		A1[2][0] = 0;
		//����1�����0�н���λ��
		A1[0][0] = 1;
		A1[1][0] = 0;
		std::swap(A1[0][1], A1[1][1]);
		std::swap(A1[0][2], A1[1][2]);
	}
	else if (A1[2][0] != 0)
	{
		A1[2][1] /= A1[2][0];
		A1[2][2] /= A1[2][0];

		//����2�����0�н���λ��
		A1[0][0] = 1;
		A1[2][0] = 0;
		std::swap(A1[0][1], A1[2][1]);
		std::swap(A1[0][2], A1[2][2]);
	}
	//Disp(A1,"A1");
	//�ڶ�������ȥA[2][1],��ȥA[0][1]
	if (A1[1][1] != 0)
	{//A1[1][1]��Ϊ0������������2��
		A1[1][2] /= A1[1][1];
		A1[1][1] = 1;
		Sub[0] = -A1[1][2];
		A1[2][2] += A1[2][1] * Sub[0];
		A1[2][1] = 0;
	}
	else if (A1[2][1] != 0)
	{
		A1[2][2] /= A1[2][1];

		A1[1][1] = 1;
		A1[2][1] = 0;
		//����2�����1�н���
		std::swap(A1[1][2], A1[2][2]);
	}
	//Disp(A1, "A1");
	int i, iCount = 0;
	for (i = 0; i < 3; i++)
		if (A1[i][i] != 0)
			iCount++;
	return iCount;
}

template<typename _T>int bInverse_Power(_T A[3][3], _T* pfEigen_Value, _T Eigen_Vector[3], _T eps)
//int bInverse_Power(double A[3][3], double* pfEigen_Value, double Eigen_Vector[3])
{//�����С������ֵ����Ӧ������������С��Χ�����Ѿ��ɹ�����
//�ɹ�������1��ʧ��:����0
#define MAX_ITER_COUNT 100
	_T L[3][3], U[3][3], yk[3], zk[3] = { 1,1,1 }, z_pre[3] = { MAX_FLOAT,MAX_FLOAT,MAX_FLOAT }, uk[3];
	_T fDiff, fDiff_Total, tau, fSum;
	//Init_LU(L, U);
	int bSucceed = 1, i, iCount = 0;

	//double tm1[3][3],tv1[3];
	//iRank(A);

	//����㷨��ȱ�ݣ��������жϾ����Ƿ���棨���ȣ�
	LU_Decompose(A, L, U, &bSucceed);
	//�����Ժ�A=L*U

	Matrix_Multiply((double*)L, 3, 3, (double*)U, 3, (double*)A);
	//Disp((double*)A, 3, 3, "A");
	//Disp(A,"A");
	if (!bSucceed || U[0][0] == 0 || U[1][1] == 0 || U[2][2] == 0)
	{//���г��ֳ���0
		*pfEigen_Value = Eigen_Vector[0] = Eigen_Vector[1] = Eigen_Vector[2] = 0;
		return 0;
	}

	while (1)
	{//�������������Է����õ��Ǹ�˹��Ԫ��������ĳЩ����û��ͨ�õķ��ݷ����ȸߵ�ԭ��
		//�п�����ӦΪ��Щ�㷨�õ�������Ԫ������΢��Щ
		//�����Luk=Zk
		Solve_Uk(L, zk, uk);

		//���� Uyk=uk
		Solve_yk(U, yk, uk);

		//Get_Max(yk, tau);
		Get_Max_1(yk, tau, i);
		tau = zk[i] >= 0 ? 1 / tau : -1.f / tau;	//�˷��ȳ�����
		zk[0] = yk[0] * tau;
		zk[1] = yk[1] * tau;
		zk[2] = yk[2] * tau;

		//Disp(zk, "zk");
		fDiff_Total = 0;
		fDiff = zk[0] - z_pre[0];
		fDiff_Total = Abs(fDiff);
		fDiff = zk[1] - z_pre[1];
		fDiff_Total += Abs(fDiff);
		fDiff = zk[2] - z_pre[2];
		fDiff_Total += Abs(fDiff);
		//printf("Diff:%f\n", fDiff_Total);
		if (fDiff_Total < eps)
			//if(abs(zk[0] - z_pre[0])+abs(zk[1] - z_pre[1]) + abs(zk[2] - z_pre[2])< ZERO_APPROCIATE)
			break;
		else if ((iCount++) >= MAX_ITER_COUNT)
		{//��������֮�����������ף�������������
			//printf("Diff:%f\n", fDiff_Total);
			//*pfEigen_Value = Eigen_Vector[0] = Eigen_Vector[1] = Eigen_Vector[2] = 0;
			//return 0;
			break;
		}

		z_pre[0] = zk[0];
		z_pre[1] = zk[1];
		z_pre[2] = zk[2];
	}

	//���
	fSum = zk[0] * zk[0] + zk[1] * zk[1] + zk[2] * zk[2];
	if (fSum < ZERO_APPROCIATE)
		return 0;
	fSum = 1 / sqrt(fSum);
	Eigen_Vector[0] = zk[0] *= fSum;
	Eigen_Vector[1] = zk[1] *= fSum;
	Eigen_Vector[2] = zk[2] *= fSum;
	*pfEigen_Value = tau;

	////��ʱ��������
	//double A1[3][3],Temp[3],fErr;
	//memcpy(A1, A, 3 * 3 * sizeof(double));
	//A1[0][0] -= *pfEigen_Value;
	//A1[1][1] -= *pfEigen_Value;
	//A1[2][2] -= *pfEigen_Value;
	////Disp(A1, "A1");
	//Matrix_x_Vector(A1, Eigen_Vector,Temp);
	//if ((fErr=abs(Temp[0]) + abs(Temp[1]) + abs(Temp[2])) > 1)
	//	printf("err:%f\n",fErr);

	return 1;
#undef MAX_ITER_COUNT
}

template<typename _T>int bInverse_Power(_T A[], int n, _T* pfEigen_Value, _T Eigen_Vector[], _T eps)
{//�����С������ֵ����Ӧ������������С��Χ�����Ѿ��ɹ�����
//�ɹ�������1��ʧ��:����0
	//����Щ��·�Ż�
	if (n == 3)
		return bInverse_Power((_T(*)[3])A, pfEigen_Value, Eigen_Vector,eps);

#define MAX_ITER_COUNT 70
	_T* L, * U, * yk, * zk, * z_pre, * uk;
	_T* pBuffer = (_T*)pMalloc((2 * n + n*4)  * sizeof(_T));
	_T fDiff_Total;
	union {
		_T fMax;
		_T tau;
	};
	int bSucceed = 1, i, iCount = 0;
	
	L = pBuffer;
	U = L + n * n;
	yk = U + n * n;
	zk = yk + n;
	z_pre = zk + n;
	uk = z_pre + n;

	for (i = 0; i < n; i++)
		zk[i] = 1, z_pre[i] = MAX_FLOAT;

	//����㷨��ȱ�ݣ��������жϾ����Ƿ���棨���ȣ�
	LU_Decompose(A, L, U,n, &bSucceed);
	if (!bSucceed)
		goto END;
	//�����Ժ�A=L*U

	while (1)
	{
		//�����Luk=Zk
		//Solve_Uk(_T L[], _T Z[], _T U[], int n)//���LU=Z�ֵ�U
		Solve_Linear_Gause(L, n, zk, uk, &bSucceed);
		if (!bSucceed)
			break;
		//Disp(uk, 1, 3, "uk");
		//���� Uyk=uk	Solve_yk(U, yk, uk);
		Solve_Linear_Gause(U, n, uk, yk, &bSucceed);
		if (!bSucceed)
			break;

		/*if (iCount == 1)
			Disp(yk, 1, 3, "yk");*/
		//���ҳ�yk�о���ֵ�����
		int iMax = 0;
		for (i = 0; i < n; i++)
			if (abs(yk[i]) > fMax)
				fMax = abs(yk[i]), iMax = i;

		tau = zk[iMax] >= 0 ? 1 / yk[iMax] : -1.f / yk[iMax];	//�˷��ȳ�����
		for (i = 0; i < n; i++)
			zk[i] = yk[i] * tau;
		/*if(iCount==1)
			Disp(zk, 1, 3, "zk");*/
		fDiff_Total = 0;
		for (i = 0; i < n; i++)
			fDiff_Total += abs(zk[i] - z_pre[i]);

		//printf("iter:%d Diff:%f\n",iCount, fDiff_Total);
		if (fDiff_Total < eps)
			break;
		else if ((iCount++) >= MAX_ITER_COUNT)
		{//��������֮�����������ף�������������
			break;
		}
		memcpy(z_pre, zk, n * sizeof(_T));
		//Disp(z_pre, 1, 3, "z_pre");
	}
	if (iCount > MAX_ITER_COUNT && fDiff_Total > eps)
		bSucceed = 0;

	//���
	Normalize(zk, n, Eigen_Vector);
	if(pfEigen_Value)
		*pfEigen_Value = tau;
END:
	if (pBuffer)
		Free(pBuffer);
	return bSucceed;
#undef MAX_ITER_COUNT
}

void Power_Method_1(double A[3][3], double Eigen_Vector[3], double* pfEigen_Value)
{//��һ���ݷ�
	double fErr, v1[3], v[3] = { 1,1,1 };	//��ȷ��һ��v0��ʹ��(x1,v0)!=0,�ʴˣ�1,1,1)��������
	//������vk=A*vk-1

	while (1)
	{
		Matrix_x_Vector(A, v, v1);
		fErr = abs(v1[0] - v[0]) + abs(v1[1] - v[1]) + abs(v1[2] - v[2]);
		printf("v=(%f,%f,%f) v1=(%f,%f,%f) err=%f\n", v[0], v[1], v[2], v1[0], v1[1], v1[2], fErr);
		v[0] = v1[0], v[1] = v1[1], v[2] = v1[2];
	}
	return;
}

template<typename _T>void Power_Method(_T A[], int n, _T Eigen_Vector[], _T* pfEigen_Value)
{//�����ֵ����Ӧ����������
	_T fDiff_Total, fMod, tau, 
		*z_pre, *zk, *yk;
		//z_pre[3] = { MAX_FLOAT,MAX_FLOAT,MAX_FLOAT }, zk[3] = { 1,1,1 }, yk[3];
	int i, j;
	int iCount = 0;
	_T* pBuffer = (_T*)pMalloc(n * 3 * sizeof(_T));
	z_pre = pBuffer;
	zk = z_pre + n;
	yk = zk + n;
	const _T eps = 1e-7f;

	while (1)
	{//�˴�������֪������Ϊֹ
		Matrix_Multiply(A, n, n, zk, 1, yk);
		tau = 0;
		for (j = 0; j < n; j++)
		{
			if (abs(yk[j]) > tau)
			{
				tau = yk[j];
				i = j;
			}
		}
		tau = zk[i] >= 0 ? 1.f / tau : -1.f / tau;	//�˷��ȳ�����
		for(j=0;j<n;j++)
			zk[j] = yk[j] * tau;
		for (fDiff_Total=0, j = 0; j < n; j++)
			fDiff_Total += abs(zk[0] - z_pre[0]);

		if (fDiff_Total < eps)
			break;

		for (j = 0; j < n; j++)
			z_pre[j] = zk[j];		
	}

	*pfEigen_Value = 1 / tau;
	fMod = fGet_Mod(zk, n);
	for(j=0;j<n;j++)
		Eigen_Vector[j] = zk[j] / fMod;
	*pfEigen_Value = 1.f / tau;
	Free(pBuffer);
}
void Power_Method(double A[3][3], double Eigen_Vector[3], double* pfEigen_Value)
{//����������ֵ����Ӧ����������
	double fDiff, fDiff_Total, fMod, tau, z_pre[3] = { MAX_FLOAT,MAX_FLOAT,MAX_FLOAT }, zk[3] = { 1,1,1 }, yk[3];
	int i;
	int iCount = 0;
	while (1)
	{//�˴�������֪������Ϊֹ
		Matrix_x_Vector(A, zk, yk);
		//Get_Max(yk, tau);
		Get_Max_1(yk, tau, i);
		printf("count:%d zk=(%f %f %f) yk=(%f %f %f) tau=%f\n", iCount++, zk[0], zk[1], zk[2], yk[0], yk[1], yk[2], tau);
		//�˴��޸ĵ�����ʹ��zk����
		tau = zk[i] >= 0 ? 1.f / tau : -1.f / tau;	//�˷��ȳ�����
		zk[0] = yk[0] * tau;
		zk[1] = yk[1] * tau;
		zk[2] = yk[2] * tau;
		//Disp(zk, "zk");
		fDiff_Total = 0;
		fDiff = zk[0] - z_pre[0];
		fDiff_Total = Abs(fDiff);
		fDiff = zk[1] - z_pre[1];
		fDiff_Total += Abs(fDiff);
		fDiff = zk[2] - z_pre[2];
		fDiff_Total += Abs(fDiff);
		if (fDiff_Total < ZERO_APPROCIATE)
			break;
		z_pre[0] = zk[0];
		z_pre[1] = zk[1];
		z_pre[2] = zk[2];
	}
	//��ʱ��1/tau��Ϊ��������ֵ��zkΪ��Ӧ����������
	*pfEigen_Value = 1 / tau;
	//Normalize zk
	fMod = sqrt(zk[0] * zk[0] + zk[1] * zk[1] + zk[2] * zk[2]);
	Eigen_Vector[0] = zk[0] / fMod;
	Eigen_Vector[1] = zk[1] / fMod;
	Eigen_Vector[2] = zk[2] / fMod;
	*pfEigen_Value = 1 / tau;
	//printf("Eigen Value:%f\n", *pfEigen_Value);
	//Disp(Eigen_Vector, "Eigen_Vector");
	return;
}
void Solve_Cubic_Eq_PQ(double p, double q, double Root[3], int* piRoot_Count)
{//��� x^3 + px + q=0��ע�⣬����ʧ�ܣ����Ľ�
	//double w=-1.f/2.f;	//�����Ǹ���������ʱ���ܸ���
	double fTemp_1, fTemp_2, fDelta;	// , w_square = -1.f / 2.f;
	double fRoot_Part_1;	// delta^(1/2)
	double fRoot_Part_2, fRoot_Part_3;
	//(q / 2) ^ 2 + (p / 3) ^ 3;
	fTemp_1 = q / 2;
	fTemp_2 = p / 3;
	fDelta = fTemp_2 * fTemp_2 * fTemp_2 + fTemp_1 * fTemp_1;

	//��delta�����б�����
	if (fDelta < 0)
	{
		double r;
		r = -p / 3.f;
		r = sqrt(r * r * r);
		double R = 2 * sqrt(-p / 3.f);
		if (abs(-0.5 * (q / r)) > 1)
			printf("err");
		double theta = acos(-0.5 * (q / r)) / 3;
		Root[0] = R * cos(theta);
		Root[1] = R * cos(theta + 3.1415926 * 2 / 3);
		Root[2] = R * cos(theta + 3.1415926 * 4 / 3);
		*piRoot_Count = 3;
		//printf("����������ʵ��\n");
	}
	else
	{
		if (p == 0 && q == 0)
		{//fDelta=0
			Root[0] = Root[1] = Root[2] = 0;
			*piRoot_Count = 0;
			printf("������0��\n");
			return;
		}
		fRoot_Part_1 = sqrt(fDelta);
		fRoot_Part_2 = -fTemp_1 + fRoot_Part_1;
		fRoot_Part_3 = -fTemp_1 - fRoot_Part_1;
		fRoot_Part_2 = fRoot_Part_2 >= 0 ? pow(fRoot_Part_2, 1.f / 3.f) : -pow(-fRoot_Part_2, 1.f / 3.f);;
		fRoot_Part_3 = fRoot_Part_3 >= 0 ? pow(fRoot_Part_3, 1.f / 3.f) : -pow(-fRoot_Part_3, 1.f / 3.f);
		Root[0] = fRoot_Part_2 + fRoot_Part_3;
		if (fDelta > 0)
		{
			*piRoot_Count = 1;
			printf("��һ��ʵ����������");
			//return;
		}
		else
		{//��ôq�ز�����0
			Root[1] = Root[2] = Root[0] * (-0.5);
			*piRoot_Count = 3;
			printf("����ʵ�������������\n");
		}
	}

	//����
	float fErr;
	for (int i = 0; i < 3; i++)
	{
		fErr = (float)(Root[i] * Root[i] * Root[i] + p * Root[i] + q);
		if (abs(fErr) > 0.1)
			printf("err");
	}

	return;
}
void Solve_Cubic(double a, double b, double c, double d, Complex_d Root[3], int* piReal_Root_Count)
{//	��ʢ��ʽ��һԪ�����������̣���������ͨ��
#define SQRT_3 1.7320508075688772935274463415059
	//Complex_d Comp[2];
	double A = b * b - 3 * a * c;
	double B = b * c - 9 * a * d;
	double C = c * c - 3 * b * d;
	double Y1, Y2, K, T, theta, Part_1, Part_2;
	double delta;
	memset(Root, 0, 3 * sizeof(Complex_d));
	if (A == 0 && B == 0)
	{
		if (b != 0)
			Root[0].real = Root[1].real = Root[2].real = -c / b;
		else if (a != 0)
			Root[0].real = Root[1].real = Root[2].real = -b / (3 * a);
		else if (c != 0)
			Root[0].real = Root[1].real = Root[2].real = -(3 * d) / c;
		printf("��һ������ʵ��\n");
	}
	else
	{
		delta = B * B - 4 * A * C;
		if (delta > 0)
		{//����·���Ѿ�����
			//�Ȱ�ʵ�������
			delta = 1.5 * sqrt(delta);

			Y1 = A * b - 1.5 * a * B + delta;
			Y2 = Y1 - 2 * delta;
			//��Y1,Y2�����η�
			Y1 = Y1 >= 0 ? pow(Y1, 1.f / 3.f) : -pow(-Y1, 1.f / 3.f);
			Y2 = Y2 >= 0 ? pow(Y2, 1.f / 3.f) : -pow(-Y2, 1.f / 3.f);
			Root[0].real = (-b - Y1 - Y2) / (3 * a);
			*piReal_Root_Count = 1;
			//����һ�Թ����
			Root[1].real = Root[2].real = (-2 * b + Y1 + Y2) / (6 * a);
			Root[1].im = SQRT_3 * (Y1 - Y2) / (6 * a);
			Root[2].im = -Root[1].im;	//����
			printf("��һ��ʵ����һ�Թ����\n");
		}
		else if (delta == 0)
		{
			K = B / A;
			Root[0].real = -b / a + K;
			Root[1].real = Root[2].real = -0.5 * K;
			printf("������ʵ����������һ�����ظ�\n");
		}
		else //delta<0
		{
			T = (2 * A * b - 3 * a * B) / (2 * sqrt(A * A * A));
			theta = acos(T) / 3;
			A = sqrt(A);
			Part_1 = A * cos(theta);
			Part_2 = SQRT_3 * A * sin(theta);
			Root[0].real = -b - 2 * Part_1;
			Root[1].real = -b + Part_1 + Part_2;
			Root[2].real = -b + Part_1 - Part_2;
			a = 1.f / (3.f * a);
			Root[0].real *= a;
			Root[1].real *= a;
			Root[2].real *= a;
			*piReal_Root_Count = 3;
			//printf("����������ȵ�ʵ��\n");
		}
	}
}
void Solve_Eigen_3x3(float A[], Complex_f Root[3])
{//�����������
	double a, b, c, d;	//����ϵ��
	float(*A_1)[3] = (float(*)[3])A;
	Complex_d Root_1[3];
	int i, iCount;

	a = -1;
	b = A_1[0][0] + A_1[1][1] + A_1[2][2];//+ (A_100 + A_111 + A_122) * r ^ 2
	c = A_1[0][1] * A_1[1][0] + A_1[0][2] * A_1[2][0] - A_1[0][0] * A_1[1][1] - A_1[0][0] * A_1[2][2]
		- A_1[1][1] * A_1[2][2] + A_1[1][2] * A_1[2][1];
	d = A_1[0][0] * A_1[1][1] * A_1[2][2] + A_1[0][1] * A_1[1][2] * A_1[2][0] + A_1[0][2] * A_1[1][0] * A_1[2][1]
		- A_1[0][0] * A_1[1][2] * A_1[2][1] - A_1[0][1] * A_1[1][0] * A_1[2][2] - A_1[0][2] * A_1[1][1] * A_1[2][0];


	//Ȼ�����A_1x^3 + bx^2+ cx + d=0�����߳���A_1,��ȡ�෴������
	a = 1;
	b = -b;
	c = -c;
	d = -d;

	//ע�⣬�˴���õĽ����ڴ�����ֵ�����Բ���ȫ�š�Ӧ�ý�һ���Ե�������߾���
	Solve_Cubic(a, b, c, d, Root_1, &iCount);
	for (i = 0; i < 3; i++)
		Root[i].im = (float)Root_1->im, Root[i].real = (float)Root_1[i].real;
	return;
}
void Solve_Eigen_Vector(float A[], int iOrder, float fEigen_Value, float Q[], int* piCount)
{//���ݸ���������ֵ����A�� ����������̵õ����������� Q�� ��������  *piCount: ��������������ϵ��
//�˴�����һ�����⣬��һ����ֵ����һ���ľ������⣬�п�������ⷽ�̵�����
	float* A_1 = (float*)malloc((iOrder + 1) * iOrder * sizeof(float)),
		* A_2 = (float*)malloc((iOrder + 1) * iOrder * sizeof(float)),
		* Q_1;

	int y, x, iRank;
	for (y = 0; y < iOrder; y++)
	{
		for (x = 0; x < iOrder; x++)
		{
			if (x == y)
				A_1[y * (iOrder + 1) + x] = A[y * iOrder + x] - fEigen_Value;
			else
				A_1[y * (iOrder + 1) + x] = A[y * iOrder + x];
		}
		A_1[y * (iOrder + 1) + iOrder] = 0;
	}
	//Disp(A_1, iOrder, iOrder+1, "A_1");


	Elementary_Row_Operation(A_1, iOrder, iOrder + 1, A_2, &iRank, &Q_1);

	*piCount = iOrder - iRank;
	Disp(Q_1, 1, 3, "Q_1");
	memcpy(Q, Q_1, *piCount * iOrder * sizeof(float));
	free(A_1);
	free(A_2);
	free(Q_1);
	return;
}
void Solve_Eigen(double A[3][3], Complex_d Root[3], int* piRoot_Count)
{
	double a, b, c, d;	//����ϵ��
	a = -1;
	b = A[0][0] + A[1][1] + A[2][2];//+ (a00 + a11 + a22) * r ^ 2
	c = A[0][1] * A[1][0] + A[0][2] * A[2][0] - A[0][0] * A[1][1] - A[0][0] * A[2][2]
		- A[1][1] * A[2][2] + A[1][2] * A[2][1];
	d = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1]
		- A[0][0] * A[1][2] * A[2][1] - A[0][1] * A[1][0] * A[2][2] - A[0][2] * A[1][1] * A[2][0];

	//Ȼ�����ax^3 + bx^2+ cx + d=0�����߳���a,��ȡ�෴������
	a = 1;
	b = -b;
	c = -c;
	d = -d;

	//ע�⣬�˴���õĽ����ڴ�����ֵ�����Բ���ȫ�š�Ӧ�ý�һ���Ե�������߾���
	Solve_Cubic(a, b, c, d, Root, piRoot_Count);

	//����
	int i;
	float fResult;
	for (i = 0; i < 3; i++)
	{
		fResult = (float)(a * Root[i].real * Root[i].real * Root[i].real + b * Root[i].real * Root[i].real + c * Root[i].real + d);
		if (abs(fResult) > 0.1)
			printf("Err Result:%f\n", fResult);
	}

	return;
}
void Solve_Cubic(double A[3][3], double Root[3], int* piRoot_Count)
{//��һԪ���η��̣�δ�㶨
	double aa, a, b, bb, c, d, p, q/*, det*/;	//����ֵһԪ���η��̵�ϵ��
	//double r = 2;	//����ramda
	double fPart_1;
	Disp((double*)A, "A");
	// -r ^ 3 + (a00 + a11 + a22) * r ^ 2 + (a01 * a10 - a00 * a11 - a00 * a22 + a02 * a20 - a11 * a22 + a12 * a21) * r + a00 * a11 * a22 - a00 * a12 * a21 - a01 * a10 * a22 + a01 * a12 * a20 + a02 * a10 * a21 - a02 * a11 * a20
	a = -1; //-r ^ 3
	b = A[0][0] + A[1][1] + A[2][2];//+ (a00 + a11 + a22) * r ^ 2
	//(a01 * a10 - a00 * a11 - a00 * a22 + a02 * a20 - a11 * a22 + a12 * a21)* r
	c = A[0][1] * A[1][0] + A[0][2] * A[2][0] - A[0][0] * A[1][1] - A[0][0] * A[2][2]
		- A[1][1] * A[2][2] + A[1][2] * A[2][1];
	//a00 * a11 * a22 - a00 * a12 * a21 
	// - a01 * a10 * a22 + a01 * a12 * a20 
	//+ a02 * a10 * a21 - a02 * a11 * a20
	d = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[0][2] * A[1][0] * A[2][1]
		- A[0][0] * A[1][2] * A[2][1] - A[0][1] * A[1][0] * A[2][2] - A[0][2] * A[1][1] * A[2][0];
	aa = a * a;
	bb = b * b;

	p = (3 * a * c - bb) / (3 * aa);
	q = (27 * aa * d - 9 * a * b * c + 2 * bb * b) / (27 * aa * a);

	//Ȼ�����ax^3 + bx^2+ cx + d=0�����߳���a,��ȡ�෴������
	a = 1;
	b = -b;
	c = -c;
	d = -d;
	//�� x=y-b/3a, ���ɻ�Ϊ y^3 + py + q =0
	//det = a * r * r * r + b * r * r + c * r + d;

	//�˴�OK�ˣ����˷���
	Solve_Cubic_Eq_PQ(p, q, Root, piRoot_Count);
	//Solve_Cubic_Eq_PQ(-3, 2, Root,piRoot_Count);

	fPart_1 = b / (3 * a);
	Root[0] -= fPart_1;
	Root[1] -= fPart_1;
	Root[2] -= fPart_1;
	//Disp(Root, "Root");

	//����
	float fErr;
	for (int i = 0; i < 3; i++)
	{
		fErr = (float)(a * (Root[i] * Root[i] * Root[i]) + b * Root[i] * Root[i] + c * Root[i] + d);
		if (abs(fErr) > 0.1)
			printf("err");
	}
	return;
}


void Matrix_Test()
{
	//Test_7();
	//return;
	int i;
	double Ramda[3],/* V[3],*/ A[3][3] = { {1, 2, 3}, { 4,5,6 }, { 7,8,9 } };
	unsigned long long tStart = iGet_Tick_Count();
	double fTotal = 0,/*fEigen_Value,*/Eigen_Vector[3][3];

	double Root[3];
	int iCount;
	Solve_Cubic(A, Root, &iCount);
	//�ݷ����������ֵ����������
	//Power_Method(A, V,&fEigen_Value);
	//���ݷ�����С����ֵ����������
	//bInverse_Power(A, V);
	//Disp(V, "V");
	for (i = 0; i < 10000000; i++)
	{
		bEigen_Value_3x3(A, Ramda, Eigen_Vector);
		fTotal += Ramda[0];
		//Disp(Ramda, "Ramda");
	}
	printf("%lld %f\n", iGet_Tick_Count() - tStart, fTotal);
	Disp(Ramda, "Ramda");
	Disp(Eigen_Vector, "Eigen_Vector");
	return;
}
template<typename _T>_T fGet_Determinant(_T* A, int iOrder)
{//������ʽ��MΪ����iOrderΪ����ʽ�Ľ�,iStrikeΪM���д�С��һ��Ԫ�ظ����� Get n-order determinant
//���㷨�����ٶȣ�������֤��ѧԭ�����Բ��ð���չ������ D=ai1*Ai1+ai2*Ai2+...+ainAin
//�Ż�˼·��������ʽ������Ԫ�����г����б任��ת��Ϊ�����Ƿ���������Խ��߳˻�
	_T* pCur, * pCofactor; //����ʽ
	int i, j, k;
	_T fDeterminant, fTotal;
	if (iOrder == 2)//�ݹ鵽����һ�㣬 2������ʽ
		return A[0] * A[3] - A[1] * A[2];

	pCofactor = (_T*)malloc((iOrder - 1) * (iOrder - 1) * sizeof(_T)); //����ʽ
	A -= (iOrder + 1);	//Ϊ������ѧԭ�����Ӧ���˴��ı��ַ���±��1��ʼ
	//����һ��չ��
	for (fTotal = 0, i = 1; i <= iOrder; i++)
	{
		//���һ�е�i��Ԫ�ص�����ʽ
		for (pCur = pCofactor, j = 2; j <= iOrder; j++)
		{
			for (k = 1; k <= iOrder; k++)
			{
				if (k == i)
					continue;	//��ȥ
				*pCur++ = A[j * iOrder + k];
			}
		}
		fDeterminant = fGet_Determinant(pCofactor, iOrder - 1);
		fTotal += (_T)(A[1 * iOrder + i] * pow(-1, 1 + i) * fDeterminant);
		//Disp((float*)pCofactor, iOrder - 1, iOrder - 1);
	}
	free(pCofactor);
	return fTotal;
}
void Gen_Cofactor(float* pM, int iOrder, float* pCofactor, int i, int j)
{//����aij��Ӧ������ʽAij, iOrrderΪpM�Ľ�
	int x, y;
	float* pCur = pCofactor;
	for (y = 0; y < iOrder; y++)
	{
		for (x = 0; x < iOrder; x++)
		{
			if (y == i || x == j)
				continue;	//����i�У���j�л�ȥ
			*pCur++ = pM[y * iOrder + x];
		}
	}
	return;
}
void Get_Adjoint_Matrix(float* pA, int iOrder, float* pAdjoint_Matrix)
{//�������󡣰�������ɾ���A��Ԫ�صĴ������������
	float* pCofactor = (float*)malloc((iOrder - 1) * (iOrder - 1) * sizeof(float));
	int i, j, iFlag;
	float* pAdj = (float*)malloc(iOrder * iOrder * sizeof(float));
	//Disp(pA, iOrder, iOrder);
	//���������󣬼���iOrder*iOrder����������ʽ
	for (i = 0; i < iOrder; i++)
	{
		for (j = 0; j < iOrder; j++)
		{
			iFlag = (int)pow(-1, (i + j) % 2);
			if (iOrder >= 3)
			{
				Gen_Cofactor(pA, iOrder, pCofactor, i, j);
				//printf("%d %d\n", i, j);
				//Disp(pCofactor, iOrder - 1, iOrder - 1);
				//printf("Determinant:%f\n\n", fGet_Determinant(pCofactor, iOrder - 1));
				pAdj[i * iOrder + j] = iFlag * fGet_Determinant(pCofactor, iOrder - 1);
			}
			else
				pAdj[i * iOrder + j] = iFlag * pA[((i + 1) % 2) * 2 + (j + 1) % 2];

		}
	}

	//Disp(pAdj, iOrder, iOrder);
	memcpy(pAdjoint_Matrix, pAdj, iOrder * iOrder * sizeof(float));
	free(pAdj);
	return;
}
template<typename _T>void Get_Inv_Matrix_Row_Op_2(_T A[], _T A_Inv[], int iOrder, int* pbSuccess)
{//��ģ���δ�ҵ�����Ԫ���������ķ���
	//�˴�һ�λش������죬����Ҫ�����Ľ�ϡ�����Ľⷽ�̣��ø���ķ���
	int y, x, x_1, i, iRank = 0, iPos, iMax, bSuccess = 1;
	int m = iOrder, n = iOrder * 2;
	unsigned int* Q, iTemp;
	_T fValue, fMax, *pfMax_Row,* A_1;
	_T eps;	// = (_T)1e-5;
	printf("��ģ���δ�ҵ�����Ԫ���������ķ���\n");
	return;

	if(typeid(_T)==typeid(double))
		eps = (_T)1e-10;
	else
		eps = (_T)1e-5;
	Q = (unsigned int*)pMalloc(&oMem_Mgr, iOrder * sizeof(unsigned int));
	A_1 = (_T*)pMalloc(&oMem_Mgr, m * n * sizeof(_T));
	for (y = 0; y < m; y++)
	{
		for (x = 0; x < iOrder; x++)
		{
			A_1[y * n + x] = A[y * iOrder + x];
			A_1[y * n + iOrder + x] = (y == x ? 1.f : 0.f);
		}
		Q[y] = { (unsigned char)y };	//ÿ����Ԫ���ڵ���
	}
	iPos = 0;
	for (y = 0; y < m; y++)
	{//�������y��x�����ƽ����������
		iMax = y;
		x_1 = y;
		fMax = A_1[Q[iMax] * n + x_1];
		for (i = y + 1; i < m; i++)
		{
			if (abs(A_1[iPos = Q[i] * n + x_1]) > abs(fMax))
			{
				fMax = A_1[iPos];
				iMax = i;
			}
		}
		if (abs(fMax) <= eps)
		{
			//printf("������\n");
			bSuccess = 0;
			goto END;
		}

		//�����Ժ�iMax��Q�������������������кţ��ʴ˲���Reverse_Lookup, �����ԪSWAP��Q�ĵ�ǰλ����
		iTemp = Q[y];
		Q[y] = Q[iMax];
		Q[iMax] = iTemp;
		//Q[y].m_iCol_Index = x_1;
		iRank++;

		pfMax_Row = &A_1[Q[y] * n];
		pfMax_Row[x_1] = 1.f;
		for (x = x_1 + 1; x < n; x++)
			pfMax_Row[x] /= fMax;

		//��Q[y]���������д��룬�о������λش��ļ��������һ�λش��Ĵ���һ�����ڴ滹ɨ�����ߣ�
		//�ɴ�������ɻش����£�ʵ���1/3
		for (i = 0; i < m; i++)
		{//i��ʾ��i��
			if (i == Q[y])
				continue;
			iPos = i * n;
			if ((fValue = A_1[iPos + x_1]) != 0)
			{//���ڶ�ӦԪ��Ϊ0�����������
				A_1[iPos + x_1] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
				for (x = x_1; x < n; x++)
					A_1[iPos + x] -= fValue * pfMax_Row[x];
			}
		}
	}
END:
	if (bSuccess)
	{   //����Q��˳�򳭵�A_Inv��
		for (i = 0; i < m; i++)
			memcpy(&A_Inv[i * iOrder], &A_1[Q[i] * n + iOrder], iOrder * sizeof(_T));
	}
	if (pbSuccess)
		*pbSuccess = bSuccess;
	if (A_1)
		Free(&oMem_Mgr, A_1);
	if (Q)
		Free(&oMem_Mgr, Q);
	return;
}

//template<typename _T>void Get_Inv_Matrix_Row_Op_1(_T A[], _T A_Inv[], int iOrder, int* pbSuccess)
//{//������Ԫ���������, ���и���Ľ��ⷶΧ
//	typedef struct Q_Item {
//		unsigned int m_iRow_Index : 24;	//��ǰ�ж�Ӧ������Ԫ����������
//		unsigned int m_iCol_Index : 24;	//			����Ԫ��Ӧ����������x����
//	}Q_Item;
//
//	int y, x, x_1, i, iRank = 0, iPos, iMax, bSuccess=1;
//	int m = iOrder, n = iOrder * 2;
//	Q_Item* Q, iTemp;
//	_T fValue, fMax;
//	_T* A_1;
//	const _T eps = (_T)1e-10;
//	union {
//		_T* pfMax_Row;
//		_T* pfBottom_Row;
//		_T* pfCur_Row;
//	};
//
//	Q = (Q_Item*)pMalloc(&oMem_Mgr, iOrder * sizeof(Q_Item));
//	A_1 = (_T*)pMalloc(&oMem_Mgr, m * n * sizeof(_T));
//	//if (A_1)
//		//memcpy(A_1, A, m * n * sizeof(_T));
//	for (y = 0; y < m; y++)
//	{
//		for (x = 0; x < iOrder; x++)
//		{
//			A_1[y * n + x] = A[y * iOrder + x];
//			A_1[y * n + iOrder + x] = y == x ? 1 : 0;
//		}
//		Q[y] = { (unsigned char)y };	//ÿ����Ԫ���ڵ���
//	}
//
//	//Disp(A_1, m, n, "A_1");
//	iPos = 0;
//	for (x_1 = 0, y = 0; y < m; y++)
//	{//�������y��x�����ƽ����������
//		//while (1)
//		//{
//			iMax = y;
//			fMax = A_1[Q[iMax].m_iRow_Index * n + x_1];
//			for (i = y + 1; i < m; i++)
//			{
//				if (abs(A_1[iPos = Q[i].m_iRow_Index * n + x_1]) > abs(fMax))
//				{
//					fMax = A_1[iPos];
//					iMax = i;
//				}
//			}
//			//if (abs(fMax) <= eps && x_1 < n - 1)
//				//x_1++;
//			//else
//				//break;
//		//}
//		if (abs(fMax) < eps)//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
//		{
//			printf("������\n");
//			bSuccess = 0;
//			goto END;
//		}
//		//�����ԪSWAP��Q�ĵ�ǰλ����
//		iTemp = Q[y];
//		Q[y] = Q[iMax];
//		Q[iMax] = iTemp;
//		Q[y].m_iCol_Index = x_1;
//		iRank++;
//
//		pfMax_Row = &A_1[Q[y].m_iRow_Index * n];
//		pfMax_Row[x_1] = 1.f;
//		for (x = x_1 + 1; x < n; x++)
//			pfMax_Row[x] /= fMax;
//
//		//�Ժ��������д���
//		for (i = y + 1; i < m; i++)
//		{//i��ʾ��i��
//			iPos = Q[i].m_iRow_Index * n;
//			if (((fValue = A_1[iPos + x_1]) != 0) && i != y)
//			{//���ڶ�ӦԪ��Ϊ0�����������
//				for (x = x_1; x < n; x++)
//					A_1[iPos + x] -= fValue * pfMax_Row[x];
//				A_1[iPos + x_1] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
//			}
//		}
//		x_1++;
//	}
//
//	//Disp(A_1, m, n,"�����б任");
//	int y1;	//�Ѿ���֪�������
//	//Ȼ��˳������һ�����ϱ任����θ���˹�����Է��̲�һ��
//	//����Ժ�A_1����������
//	for (y = iRank - 1; y > 0; y--)
//	{//�߼��ϴ�����һ�����ϣ�ʵ������Qָ·
//		pfBottom_Row = &A_1[Q[y].m_iRow_Index * n];
//		x_1 = Q[y].m_iCol_Index;	//ǰ���Ѿ��õ���������Ԫλ��
//
//		for (y1 = y - 1; y1 >= 0; y1--)
//		{
//			iPos = Q[y1].m_iRow_Index * n;
//			x = x_1;	//�����е�xλ��
//			fValue = A_1[iPos + x];
//			A_1[iPos + x] = 0;
//			for (x++; x < n; x++)
//				A_1[iPos + x] -= fValue * pfBottom_Row[x];
//		}
//	}
//	/*for (i = 0; i < m; i++)
//		printf("i:%d Line_Index:%d\n", i, Q[i].m_iRow_Index);*/
//	//��A_1���Ұ벿�ְ���Q�����������������
//	for (i = 0; i < m; i++)
//		memcpy(&A_Inv[i * iOrder], &A_1[Q[i].m_iRow_Index * n + iOrder], iOrder * sizeof(_T));
//
//END:
//	if (pbSuccess)
//		*pbSuccess = bSuccess;
//	if (A_1)
//		Free(&oMem_Mgr, A_1);
//	if (Q)
//		Free(&oMem_Mgr, Q);
//}

template<typename _T>void Disp_Get_Inv_AAt_Row_Op(_T A[], int iOrder,const char *pcCaption)
{//��ʾ�������˫������
	int iPos, y, x;
	if(pcCaption)
		printf("%s\n", pcCaption);
	for (iPos = 0, y = 0; y < iOrder; y++)
	{
		for (x = 0; x < y; x++)
			printf("%f\t", 0);
		for (; x < y+iOrder + 1; x++)
			printf("%.8f\t", A[iPos++]);
		for (; x < iOrder * 2; x++)
			printf("%f\t", 0);
		printf("\n");
	}
	return;
}

template<typename _T>void Get_Inv_AAt_3x3(_T M[3*3], _T Inv[3*3], int* pbSuccess)
{//schur��Ԫ��������Ҫ��3x3�����, �Ѿ���ԭ����10��
	_T Aux[4 * 3];
	int bSuccess = 1;

	//��һ���������ԭ����+I
	Aux[0] = M[0], Aux[1] = M[1], Aux[2] = M[2];
	Aux[4] = M[4], Aux[5] = M[5];
	Aux[8] = M[8];
	Aux[6] = Aux[9] = Aux[10] = 0;
	Aux[3] = Aux[7] =Aux[11] = 1;
	//Disp_Get_Inv_AAt_Row_Op(Aux, 3, "Aux");

	const _T eps = (_T)1e-10;
	_T fPivot, fFactor;
	//�ȱ任����������ʽ
	//��һ��
	if (Abs(Aux[0]) < eps){	/*printf("����ԪΪ��%f\n", Aux[0]);*/	bSuccess = 0;		goto END;	}
	fPivot = 1.f/Aux[0];
		fFactor = Aux[1] * fPivot;
		Aux[4] -= fFactor * Aux[1];		Aux[5] -= fFactor * Aux[2];		Aux[6] -= fFactor * Aux[3];
		fFactor = Aux[2] * fPivot;
		Aux[8] -= fFactor * Aux[2];		Aux[9] -= fFactor * Aux[3];
		//����Ԫ�и���һ��
		Aux[1] *= fPivot; Aux[2] *= fPivot; Aux[3] *= fPivot;
		//Disp_Get_Inv_AAt_Row_Op(Aux, 3, "Aux");

	//�ڶ���
	if (Abs(Aux[4]) < eps){	/*printf("����ԪΪ��%f\n", Aux[4]);*/	bSuccess = 0;		goto END;	}
	fPivot = 1.f/Aux[4];
		fFactor = Aux[5] * fPivot;
		Aux[8] -= fFactor * Aux[5];	Aux[9] -= fFactor * Aux[6]; Aux[10] -= fFactor*Aux[7];
		//����Ԫ�и���һ��
		Aux[5] *= fPivot;Aux[6] *= fPivot;Aux[7] *= fPivot;

	//������
	if (Abs(Aux[8]) < eps){	/*printf("����ԪΪ��%f\n", Aux[8]);*/	bSuccess = 0;		goto END;	}
	fPivot = 1.f/Aux[8];
		Aux[9] *= fPivot; Aux[10] *= fPivot; Aux[11] *= fPivot; 
		//Disp_Get_Inv_AAt_Row_Op(Aux, 3, "Aux");
	
	//if (Aux[0] < 0)Aux[0] = -Aux[0];if (Aux[4] < 0)Aux[4] = -Aux[4];if (Aux[8] < 0)Aux[8] = -Aux[8];
	//if (Aux[0] < eps || Aux[4] < eps || Aux[8] < eps)
	//{
	//	//printf("����ԪΪ0\n");	
	//	bSuccess = 0;
	//	goto END;
	//}

	//�ش�
	//��һ��
	Inv[1] = Inv[3] = (Aux[6] -= Aux[5] * Aux[9]);	Inv[4] =(Aux[7] -= Aux[5] * Aux[10]);		//����һ����λ������Ŀ��
	Aux[3] -= Aux[2] * Aux[9];

	//Aux[3] -= Aux[1] * Aux[6];	//����һ����λ
	Inv[0] = Aux[3] - Aux[1] * Aux[6];
	Inv[6] = Inv[2] = Aux[9]; Inv[7] = Inv[5] = Aux[10]; Inv[8] = Aux[11];

	//Disp(Inv, 3, 3, "Inv");
	//Disp_Get_Inv_AAt_Row_Op(Aux, 3, "Aux");
END:
	if (pbSuccess)
		*pbSuccess = bSuccess;
	return;
}

template<typename _T>void Get_Inv_AAt_Row_Op(_T* pM, _T* pInv, int iOrder, int* pbSuccess)
{//���ڶԳƾ���, ������Ȼ�ǶԳƾ��󣬳��Կռ���٣�����Ҳ���ٿ���
	//Ȼ����û�п��ˣ����вο�����
	_T* pAux;
	int y, x, i,bSuccess = 1;
	int iOrder_Plus_1=iOrder+1, iOrder_Minus_1 = iOrder -1,iPos;	
	int iStep_1 = -1;

	pAux = (_T*)pMalloc(&oMem_Mgr, iOrder_Plus_1 * iOrder * sizeof(_T));
	if (!pAux)
	{
		printf("Fail to allocate memory in Get_Inv_AAt_Row_Op\n");
		return;
	}
	memset(pAux, 0, iOrder_Plus_1 * iOrder * sizeof(_T));

	//��һ���������ԭ����+I
	for (iPos=0,y = 0; y < iOrder; y++)
	{
		for (x = y; x < iOrder; x++,iPos++)
			pAux[iPos] = pM[y * iOrder + x];
		pAux[iPos+y] = 1;
		iPos += y+1;
	}
	//Disp_Get_Inv_AAt_Row_Op(pAux, iOrder, "Aux");

	const _T eps = (_T)1e-10;
	_T* pfPivot_Row = pAux, * pfCur_Row;
	for (y = 0; y < iOrder; y++)
	{
		//��Ԫ���ǶԽ���Ԫ��
		_T fPivot = *pfPivot_Row;
		if (Abs(fPivot) < eps)
		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
			printf("����ԪΪ��%f\n", fPivot);
			bSuccess = 0;
			goto END;
		}
		//Ϊ�˿�㣬�����õ���
		fPivot = 1.f / fPivot;
		pfCur_Row = pfPivot_Row + iOrder_Plus_1;
		int iOrder_Minus_y = iOrder - y,
			iOrder_Minus_i = iOrder_Minus_1;

		//����Ԫ�е���һ��һֱɨ�����һ��
		for (i = 1; i <iOrder_Minus_y; i++,iOrder_Minus_i--)
		{
			if ( _T fRow_Head =pfPivot_Row[i])
			{//���ڶ�ӦԪ��Ϊ0�����������
				//_T fFactor = fRow_Head / fPivot;
				_T fFactor = fRow_Head * fPivot;
				for (x = 0; x <=iOrder_Minus_i; x++)
					pfCur_Row[x] -= fFactor * pfPivot_Row[x + i];
			}
			pfCur_Row += iOrder_Plus_1;
		}

		//�ٽ���Ԫ��ɨһ��
		for (x = 1; x < iOrder_Plus_1; x++)
			//pfPivot_Row[x] /= pfPivot_Row[0];
			pfPivot_Row[x] *= fPivot;
		//pfPivot_Row[0] = 1;	//�Ǳ�Ҫ��
		pfPivot_Row += iOrder_Plus_1;
		//Disp_Get_Inv_AAt_Row_Op(pAux, iOrder, "Aux");
	}

	//Disp_Get_Inv_AAt_Row_Op(pAux, iOrder, "Aux");

	//����һ�����ϻش���Ԫ
	
	y = iOrder - 1;
	pfPivot_Row = &pAux[y * iOrder_Plus_1+1];
	for (y = iOrder - 1; y >= 1; y--,iStep_1--)
	{
		//printf("Pivot: %f\n", pfPivot_Row[0]);
		pfCur_Row = pfPivot_Row - iOrder;
		int x_End = y;
		for(x_End=y;x_End>=1;x_End--)
		{
			_T fPivot = pfCur_Row[iStep_1];
			//printf("\tcur_Row:%f Pivot:%f\n\t\t", *pfCur_Row,fPivot);
			for (x = 0; x < x_End; x++)
			{
				//printf("(%f %f) ", pfCur_Row[x],pfPivot_Row[x]);
				pfCur_Row[x] -= fPivot * pfPivot_Row[x];
			}
			//pfCur_Row[iStep_1] = 0;	//���Ǳ����
			//printf("\n");
			pfCur_Row -= iOrder;
		}
		//Disp_Get_Inv_AAt_Row_Op(pAux, iOrder, "Aux");
		pfPivot_Row -= iOrder;
	}
	//Disp_Get_Inv_AAt_Row_Op(pAux, iOrder, "Aux");

	//����pInv��
	pfCur_Row = &pAux[iOrder];
	for (y = 0; y < iOrder; y++)
	{
		//printf("%f ", *pfCur_Row);
		for (x = 0; x <y; x++,iPos++)
		{
			//printf("%f ", pfCur_Row[x]);
			pInv[y * iOrder + x] = pInv[x*iOrder + y] = pfCur_Row[x];
		}
		pInv[y * iOrder + x] = pfCur_Row[x];
		//printf("%\n");
		pfCur_Row += iOrder;
	}
	//Disp(pInv, iOrder, iOrder, "Inv");
END:
	if (pbSuccess)
		*pbSuccess = bSuccess;
	Free(&oMem_Mgr, pAux);
	return;
}
template<typename _T> void Get_Inv_Matrix_Row_Op(_T* pM, _T* pInv, int iOrder, int* pbSuccess)
{//�Ѿ�����
	//�����ð�������������̫�����˴��ó����б任�㣬��ʵ���Ǹ�˹�󷽳̵�һ��
	//��һ����Bufferװ����M���б任
	_T* pAux;
	int y, x, i, iRow_Size = iOrder * 2, iPos, bRet = 1;
	_T fMax, * pfCur_Row, fValue;
	/*if (iOrder > 256)
	{
		printf("Too large order:%d\n", iOrder);
		*pbSuccess = 0;
		return;
	}*/
	pAux = (_T*)pMalloc(&oMem_Mgr, iOrder * 2 * iOrder * sizeof(_T));
	if (!pAux)
		return;

	//��ʼ��ֵ
	for (y = 0; y < iOrder; y++)
	{
		iPos = y * iRow_Size;
		for (x = 0; x < iOrder; x++)
			pAux[iPos + x] = pM[y * iOrder + x];
		iPos += iOrder;
		for (x = 0; x < iOrder; x++)
			pAux[iPos + x] = y == x ? 1.f : 0.f;
	}

	//Disp(pAux, iOrder, iOrder * 2,"\n");
	const _T eps = (_T)1e-10;
	for (y = 0; y < iOrder; y++)
	{
		pfCur_Row = &pAux[y * iRow_Size];
		fMax = pAux[y * iRow_Size + y];
		if (Abs(fMax) < eps)
		{
			bRet = 0;
			printf("����ԪΪ0\n");
			goto END;
		}
		pfCur_Row[y] = 1.f;
		for (x = y + 1; x < iRow_Size; x++)
			pfCur_Row[x] /= fMax;
		//Disp(pAux, iOrder, iOrder * 2, "\n");

		//�Ժ��������д���
		for (i = y + 1; i < iOrder; i++)
		{//i��ʾ��i��
			//iPos = Q[i] * iRow_Size;
			iPos = i * iRow_Size;
			if ((fValue = pAux[iPos + y]) != 0)
			{//���ڶ�ӦԪ��Ϊ0�����������
				for (x = y + 1; x < iRow_Size; x++)
					pAux[iPos + x] -= fValue * pfCur_Row[x];
				pAux[iPos + y] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
			}
			//Disp(pAux, iOrder, iOrder * 2, "\n");
		}
	}
	//Disp(pAux, iOrder, iOrder * 2, "\n");

	_T* pfBottom_Row;
	int y_1;
	//Ȼ��˳������һ�����ϱ任����θ���˹�����Է��̲�һ��
	for (y = iOrder - 1; y > 0; y--)
	{//�߼��ϴ�����һ�����ϣ�ʵ������Qָ·
		pfBottom_Row = &pAux[y * iRow_Size];
		for (y_1 = y - 1; y_1 >= 0; y_1--)
		{
			//iPos = Q[y_1] * iRow_Size;
			iPos = y_1 * iRow_Size;
			x = y;	//�����е�xλ��
			fValue = pAux[iPos + x];
			pAux[iPos + x] = 0;
			for (x++; x < iRow_Size; x++)
				pAux[iPos + x] -= fValue * pfBottom_Row[x];
			//Disp(pAux, iOrder, iRow_Size, "\n");
		}
	}

	//��󳭵�Ŀ�����
	for (y = 0; y < iOrder; y++)
	{
		iPos = y * iRow_Size + iOrder;
		for (x = 0; x < iOrder; x++, iPos++, pInv++)
			*pInv = pAux[iPos];
	}

END:
	if (pbSuccess)
		*pbSuccess = bRet;
	Free(&oMem_Mgr, pAux);
	//free(pAux);
	return;
}

void Get_Inv_Matrix(float* pM, float* pInv, int iOrder, int* pbSuccess)
{//��pM��������ð������  Inv(A)= (A*)/|A| ��ʱֻ��3�����ϵľ���
	float fDet = fGet_Determinant(pM, iOrder);

	int i, j;
	if (abs(fDet) < ZERO_APPROCIATE)
	{
		*pbSuccess = 0;
		return;
	}
	float* pAdjoint_Matrix = (float*)malloc(iOrder * iOrder * sizeof(float));
	Get_Adjoint_Matrix(pM, iOrder, pAdjoint_Matrix);

	//���ˣ���������Ѿ�OK������pAdjoint_Matrix�Ȼ�������ս�����ܸ��̿���
	//׼ȷ���㷨Ӧ���� A(-1)= (A*)'/|A|	, �������Ҫת��һ�·���
	for (i = 0; i < iOrder; i++)
		for (j = 0; j < iOrder; j++)
			pInv[i * iOrder + j] = pAdjoint_Matrix[j * iOrder + i] / fDet;
	*pbSuccess = 1;

	free(pAdjoint_Matrix);
	return;
}

template<typename _T>_T fLinear_Equation_Check(_T* A, int iOrder, _T* B, _T* X, _T eps)
{//�����Է���Ax=B����, ��<epsΪ��
	_T* X1, fError;
	int i;
	X1 = (_T*)malloc(iOrder * sizeof(_T));

	Matrix_Multiply(A, iOrder, iOrder, X, 1, X1);
	for (i = 0, fError = 0; i < iOrder; i++)
		fError += abs(B[i] - X1[i]);

	if (fError < eps)
		printf("Correct:%f\n", fError);
	else
		printf("Incorrect, Error Sum:%f\n", fError);
	free(X1);
	return (_T)fError;
}
template<typename _T>void SVD_Get_Solution(SVD_Info oSVD, _T X[])
{//��svd�ֽ���vt�����з������������X��
	memcpy(X, &((_T*)oSVD.Vt)[(oSVD.w_A - 1) * oSVD.w_A], oSVD.w_A * sizeof(_T));

	/*for (int i = 0; i < oSVD.w_A; i++)
		X[i] = ((_T*)oSVD.Vt)[i * oSVD.w_A + oSVD.w_A - 1];*/
}
template<typename _T> void Solve_Linear_Contradictory(_T* A, int m, int n, _T* B, _T* X, int* pbSuccess)
{//��������С���˷���ì�ܷ����顣�ؼ����� A'Ax=A'B, ���н⣬��xΪ Ax=B����С���˽�
//����������⾡����ͨ�����߱��λ�Ϊ�������⣬Ȼ���ý�ì�ܷ�����ķ�������ͺð�
	_T* At, *AtA, *AtB;
	int bResult = 1;
	
	if (B)	//����η���
	{
		At = (_T*)pMalloc((m * n + n * n + n) * sizeof(_T));
		AtA = At + m * n;
		AtB = AtA + n * n;
		Matrix_Transpose(A, m, n, At);
		//Matrix_Multiply(At, n, m, A, n, AtA);
		Transpose_Multiply(At, n, m, AtA);
		//Disp(AtA, 12, 12, "AtA");
		Matrix_Multiply(At, n, m, B, 1, AtB);
		Solve_Linear_Gause(AtA, n, AtB, X, &bResult);
		Free(At);
	}else
	{
		//ì����η����飬��ʱ���м��ַ�����
		if (m == n)
		{
			int iCount;
			_T* pX = (_T*)pMalloc(n * (n - 1) * sizeof(_T));
			Solve_Linear_Solution_Construction_1(A, n, n,(_T*)NULL, &bResult,pX,&iCount);
			if (iCount == 1)
			{
				memcpy(X, pX, n * sizeof(_T));
				Free(pX);
				*pbSuccess = 1;
				return;
			}
			Free(pX);
		}
		

		//// ����1�����ٵķ��ݷ�
		//AtA = (_T*)pMalloc((n * n) * sizeof(_T));
		//Transpose_Multiply(A, m, n, AtA,0);
		//if (!bInverse_Power(AtA, n, (_T*)NULL, X))
		//{
		//	bResult = 0;
		//	memset(X, 0, 12 * sizeof(_T));
		//}
		//Free(AtA);
		//Test_Linear(A, n, (_T*)NULL, X);

		//����2�����ٵ�svd����
		SVD_Info oSVD;
		SVD_Alloc(m, n, &oSVD, A);
		svd_3(A, oSVD, &bResult);
		SVD_Get_Solution(oSVD, X);
		Free_SVD(&oSVD);
		Test_Linear(A, n, (_T*)NULL, X);

		//���ڶԷ��ݷ��˽ⲻ���ʱû���㹻�������������ָ���׳���ؼ�
		//����1������˭���ã�2���Ƿ����ĳЩsvd�ֽܷ⵫���ݷ��㲻������
		//��������ݷ�������������LU�ֽ�Ҫ��������һ������
	}
	
		
	/*Disp(AtA, n, n, "AtA");
	Disp(AtB, n, 1, "AtB");
	Disp(X, n, 1, "X");*/
	//if (bResult)
	//{//����
	//	_T fTotal = 0, * B_1 = (_T*)malloc(m * sizeof(_T));
	//	int i;
	//	Matrix_Multiply(A, m, n, X, 1, B_1);
	//	for (i = 0; i < m; i++)
	//		fTotal += (B[i] - B_1[i]) * (B[i] - B_1[i]);
	//	//printf("Error Sum:%f\n", fTotal);
	//	printf("%f\n", fTotal);
	//	free(B_1);
	//	//Disp(X, 1, n, "x");
	//}
	if(pbSuccess)
		*pbSuccess = bResult;	
}

//void Solve_Linear_Contradictory(float* A, int m, int n, float* B, float* X, int* pbSuccess)
//{//��������С���˷���ì�ܷ����顣�ؼ����� A'Ax=A'B, ���н⣬��xΪ Ax=B����С���˽�
////����������⾡����ͨ�����߱��λ�Ϊ�������⣬Ȼ���ý�ì�ܷ�����ķ�������ͺð�
//
//	float* At = (float*)malloc((m * n + n * n + m) * sizeof(float));
//	float* AtA = At + m * n;
//	float* AtB = AtA + n * n;
//	int bResult = 0;
//
//	Matrix_Transpose(A, m, n, At);
//	Matrix_Multiply(At, n, m, A, n, AtA);
//	Matrix_Multiply(At, n, m, B, 1, AtB);
//
//	Solve_Linear_Gause(AtA, n, AtB, X, &bResult);
//	/*Disp(AtA, n, n, "AtA");
//	Disp(AtB, n, 1, "AtB");
//	Disp(X, n, 1, "X");*/
//	if (bResult)
//	{//����
//		float fTotal = 0, * B_1 = (float*)malloc(m * sizeof(float));
//		int i;
//		Matrix_Multiply(A, m, n, X, 1, B_1);
//		for (i = 0; i < m; i++)
//			fTotal += (B[i] - B_1[i]) * (B[i] - B_1[i]);
//		printf("Error Sum:%f\n", fTotal);
//		free(B_1);
//		Disp(X, 1, n, "x");
//	}
//	*pbSuccess = bResult;
//	free(At);
//}

template<typename _T>int iGet_Rank(_T* A, int m, int n)
{//���ø�˹����Ԫ���������ȣ��ó����б任��
	int iRank = 0;
	int iMax, iTemp, iPos, * pQ;
	int y, x, i, iRow_To_Test;
	_T* Ai = (_T*)pMalloc(&oMem_Mgr, m * n * sizeof(_T));
	_T fMax, * pfMax_Row, fValue;
	iPos = 0;
	pQ = (int*)pMalloc(&oMem_Mgr, m * sizeof(int));
	if (Ai)
		memcpy(Ai, A, m * n * sizeof(_T));
	else
	{
		printf("Fail to malloc in fGet_Rank\n");
		return -1;
	}
	iRow_To_Test = Min(m, n);
	for (y = 0; y < m; y++)
		pQ[y] = y;	//ÿ����Ԫ���ڵ���

	for (y = 0; y < iRow_To_Test; y++)
	{
		iMax = y;
		fMax = Ai[pQ[iMax] * n + y];

		for (i = y + 1; i < m; i++)
		{//Ѱ������Ԫ
			if (abs(Ai[iPos = pQ[i] * n + y]) > abs(fMax))
			{
				fMax = Ai[iPos];
				iMax = i;
			}
		}
		if (abs(fMax) < ZERO_APPROCIATE)
		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
			printf("������,����ԪΪ��%f\n", fMax);
			Disp(Ai, m, n, "Ai");
			continue;
		}
		else
			iRank++;

		//�����ԪSWAP��Q�ĵ�ǰλ����
		iTemp = pQ[y];
		pQ[y] = pQ[iMax];
		pQ[iMax] = iTemp;

		//��iMax���ڵ��н���ϵ�����㣬��ϵ��/=A[y][y]
		pfMax_Row = &Ai[pQ[y] * n];
		pfMax_Row[y] = 1.f;
		for (x = y + 1; x < n; x++)
			pfMax_Row[x] /= fMax;

		//�Ժ��������д���
		for (i = y + 1; i < m; i++)
		{//i��ʾ��i��
			iPos = pQ[i] * n;
			if ((fValue = Ai[iPos + y]) != 0)
			{//���ڶ�ӦԪ��Ϊ0�����������
				for (x = y + 1; x < n; x++)
					Ai[iPos + x] -= fValue * pfMax_Row[x];
				Ai[iPos + y] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
			}
			//Disp(Ai, iOrder, iOrder + 1, "\n");
		}
	}
	Free(&oMem_Mgr, Ai);
	Free(&oMem_Mgr, pQ);
	return iRank;
}
template<typename _T>void Disp_Ai(_T A[], int iOrder)
{//��ʽһ���ԽǾ���������
	int i, j,iPos;
	for (iPos=i = 0; i < iOrder;i++)
	{
		for (j = 0; j < i; j++)
			printf("\t");
		for (; j < iOrder+1; j++, iPos++)
			printf("%.1f\t", A[iPos]);
		printf("\n");
	}
	return;
}
template<typename _T>void Solve_Linear_Gause_AAt(_T* A, int iOrder, _T* B, _T* X, int* pbSuccess)
{//��Գƾ��󷽳̣��������ǹ���������󣬿�ܶ�
//������Ԫ��ȣ���double�£������־����㹻�����£���������Ԫ����������
//����float�£����־��Ȳ����£�������Ԫ�����ܴ�����
	int y, x, i;
	union {
		int iRow_Size;
		int iPos;
	};

	_T fPivot, * pfPivot_Row;
	const _T eps = (_T)1e-10;
	int bSuccess = 1;
	_T* Ai = (_T*)pMalloc(&oMem_Mgr, (iOrder + 1 + 2)*iOrder/2 * sizeof(_T));
	if (!Ai)
	{
		bSuccess = 0;
		goto END;
	}
	iPos = 0;
	for (y = 0; y < iOrder; y++)
	{
		for (x = y; x < iOrder; x++, iPos++)
			Ai[iPos] = A[y * iOrder + x];
		Ai[iPos++] = B[y];
	}
	//Disp_Ai(Ai, iOrder);

	iRow_Size = iOrder + 1;
	//Ϊ�˱�����⣬����iMax��ʾQ�е�������������Ai�е��к�
	pfPivot_Row = Ai;

	_T* pfCur_Row;
	int iCur_Row_Size;
	for (y = 0; y < iOrder; y++)
	{
		//��Ԫ���ǶԽ���Ԫ��
		fPivot = *pfPivot_Row;
		if (abs(fPivot) < eps)
		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
			printf("����ԪΪ��%f\n", fPivot);
			bSuccess = 0;
			goto END;
		}
		//Ϊ�˿�㣬�����õ���
		fPivot = 1.f / fPivot;
		
		//��һ�´˴���һ������
		pfCur_Row = pfPivot_Row + iRow_Size;
		iCur_Row_Size = iRow_Size - 1;

		//�˴������е����⣬�ǲ���Ӧ�ô�y+1��iOrder�أ�
		//for (i = 1; i < iOrder; i++)
		for(i=1;iCur_Row_Size;i++)
		{//i��ʾ��i��
			union {
				_T fFactor;
				_T fRow_Head;
			};			
			if ( fRow_Head =pfPivot_Row[i])
			{//���ڶ�ӦԪ��Ϊ0�����������
				//fFactor = fRow_Head / fPivot;
				fFactor = fRow_Head * fPivot;
				for (x = 0; x < iCur_Row_Size; x++)
				{
					pfCur_Row[x] -= fFactor * pfPivot_Row[x+i];
				}
			}
			pfCur_Row += iCur_Row_Size;
			iCur_Row_Size--;
		}

		//pfPivot_Row[0] = 1.f;	//���Ҳ���ܲ���Ҫ
		for (x = 1; x < iRow_Size; x++)
			//pfPivot_Row[x] /= fPivot;
			pfPivot_Row[x] *= fPivot;
		pfPivot_Row += iRow_Size;
		iRow_Size--;
	}
	//Disp_Ai(Ai, iOrder);	
	//�ش�����Q[iOrder-1]��ʼ�ش���������һ�����ϻش�
	pfCur_Row = &Ai[(iOrder + 1 + 2) * iOrder / 2 - 1];
	//��һ����
	X[iOrder - 1] = pfCur_Row[0];

	//�������п�ʼ��iCur_Row_Size������ 1
	//pfCur_Row��1��ĵ�һ��Ԫ�ؿ�ʼ
	iCur_Row_Size = 1;
	pfCur_Row -= 3;
	for (y = iOrder - 2; y >= 0; y--)
	{
		//�������ϻش�
		//fValue=b
		_T fValue = pfCur_Row[iCur_Row_Size];
		for (x = 0; x < iCur_Row_Size; x++,iPos++)
			fValue -= pfCur_Row[x] * X[iOrder-iCur_Row_Size+x];
			//Ai[iPos + x] = 0;		//�˴����Ǳ���ģ�������0���ÿ�һЩ����
		X[y] = fValue;
		//Ai[iPos + iOrder] = fValue;	//�˴����Ǳ��룬�ÿ�����
		iCur_Row_Size++;
		pfCur_Row -= iCur_Row_Size+2;	//1ռһ����bռһ��
	}

	//Disp(X, iOrder, 1, "X");
	//����
	//fLinear_Equation_Check(A, iOrder, B, X, (_T)ZERO_APPROCIATE);
END:
	*pbSuccess = bSuccess;	
	if (Ai)
		Free(&oMem_Mgr, Ai);

	return;
}
template<typename _T>void Solve_Linear_Gause(_T* A, int iOrder, _T* B, _T* X, int* pbSuccess)
{//�ø�˹����Ԫ��������Է�����, Ҫ�㣺
	//1����˹����ȫ�ȼ������б任��ֻ����û��������Ĺ������������ǲ������ǽ���Ԫ��Ϊ1
	//2��ѡ����Ԫ��ԭ���Ǳ�֤�� ����ϵ��/aijʱ��ĸ������̫С����ĸС����
	//3,������Ԫ̫С��<eps)��������㲻���ȣ��˳���ʵ������Ԫ̫С�ᵼ�º���ĳ����������
	int y, x, i, iRow_Size;
	int iMax, iTemp, iPos, * pQ;
	_T fMax, * pfMax_Row, fValue;
	const _T eps = (_T)1e-10;
	int bSuccess = 1;
	pQ = (int*)pMalloc(&oMem_Mgr, iOrder * sizeof(int));
	_T* Ai = (_T*)pMalloc(&oMem_Mgr, (iOrder + 1) * iOrder * sizeof(_T));
	if (!pQ || !Ai)
	{
		bSuccess = 0;
		goto END;
	}

	iPos = 0;
	for (y = 0; y < iOrder; y++)
	{
		for (x = 0; x < iOrder; x++, iPos++)
			Ai[iPos] = A[y * iOrder + x];
		Ai[iPos++] = B[y];
		pQ[y] = y;	//ÿ����Ԫ���ڵ���
	}

	//Disp(Ai, iOrder, iOrder + 1,"\n");
	iRow_Size = iOrder + 1;
	//Ϊ�˱�����⣬����iMax��ʾQ�е�������������Ai�е��к�
	for (y = 0; y < iOrder; y++)
	{
		iMax = y;	//�о�������ˣ�Ӧ���� iMax = pQ[i]
		fMax = Ai[pQ[iMax] * iRow_Size + y];
		for (i = y + 1; i < iOrder; i++)
		{//Ѱ������Ԫ
			if (abs(Ai[iPos = pQ[i] * iRow_Size + y]) > abs(fMax))
			{
				fMax = Ai[iPos];
				iMax = i;
			}
		}
		/*if (iMax != y)
			printf("%d\n", iMax);*/

		if (abs(fMax) < eps)
		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
			printf("������,����ԪΪ��%f\n", fMax);
			bSuccess = 0;
			goto END;
		}

		//�����ԪSWAP��Q�ĵ�ǰλ����
		iTemp = pQ[y];
		pQ[y] = pQ[iMax];
		pQ[iMax] = iTemp;

		//��iMax���ڵ��н���ϵ�����㣬��ϵ��/=A[y][y]
		pfMax_Row = &Ai[pQ[y] * iRow_Size];
		pfMax_Row[y] = 1.f;
		for (x = y + 1; x < iRow_Size; x++)
			pfMax_Row[x] /= fMax;

		//Disp(Ai, iOrder, iOrder + 1, "\n");

		//�Ժ��������д���
		for (i = y + 1; i < iOrder; i++)
		{//i��ʾ��i��
			iPos = pQ[i] * iRow_Size;
			if ((fValue = Ai[iPos + y]) != 0)
			{//���ڶ�ӦԪ��Ϊ0�����������
				for (x = y + 1; x < iRow_Size; x++)
					Ai[iPos + x] -= fValue * pfMax_Row[x];
				Ai[iPos + y] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
			}
			//Disp(Ai, iOrder, iOrder + 1, "\n");
		}
	}

	//��һ����
	X[iOrder - 1] = Ai[pQ[iOrder - 1] * iRow_Size + iOrder];

	//�ش�����Q[iOrder-1]��ʼ�ش���������һ�����ϻش�
	for (y = iOrder - 2; y >= 0; y--)
	{
		//�������ϻش�
		iPos = pQ[y] * iRow_Size;
		//fValue=b
		fValue = Ai[iPos + iOrder];
		for (x = y + 1; x < iOrder; x++)
		{
			fValue -= Ai[iPos + x] * X[x];
			Ai[iPos + x] = 0;		//�˴����Ǳ���ģ�������0���ÿ�һЩ����
		}
		X[y] = fValue;
		Ai[iPos + iOrder] = fValue;	//�˴����Ǳ��룬�ÿ�����
	}

END:
	//Disp(Ai, iOrder, iOrder + 1, "\n");
	//*pbSuccess = 1;
	*pbSuccess = bSuccess;
	if (pQ)
		Free(&oMem_Mgr, pQ);
	if (Ai)
		Free(&oMem_Mgr, Ai);
	//����
	//Linear_Equation_Check(A, iOrder, B, X, (_T)ZERO_APPROCIATE);
	return;
}

//template<typename _T>void Solve_Linear_Gause(_T* A, int iOrder, _T* B, _T* X, int* pbSuccess)
//{//�ø�˹����Ԫ��������Է�����, Ҫ�㣺
//	//1����˹����ȫ�ȼ������б任��ֻ����û��������Ĺ������������ǲ������ǽ���Ԫ��Ϊ1
//	//2��ѡ����Ԫ��ԭ���Ǳ�֤�� ����ϵ��/aijʱ��ĸ������̫С����ĸС����
//	//3,������Ԫ̫С��<eps)��������㲻���ȣ��˳���ʵ������Ԫ̫С�ᵼ�º���ĳ����������
//	int y, x, i, iRow_Size;
//	int iMax, iTemp, iPos, Q[256];
//	_T fMax, * pfMax_Row, fValue;
//	if (iOrder > 256)
//	{
//		printf("Too large order:%d\n", iOrder);
//		*pbSuccess = 0;
//		return;
//	}
//
//	_T* Ai = (_T*)malloc((iOrder + 1) * iOrder * sizeof(_T));
//	iPos = 0;
//	for (y = 0; y < iOrder; y++)
//	{
//		for (x = 0; x < iOrder; x++, iPos++)
//			Ai[iPos] = A[y * iOrder + x];
//		Ai[iPos++] = B[y];
//		Q[y] = y;	//ÿ����Ԫ���ڵ���
//	}
//	//Disp(Ai, iOrder, iOrder + 1,"\n");
//	iRow_Size = iOrder + 1;
//	//Ϊ�˱�����⣬����iMax��ʾQ�е�������������Ai�е��к�
//	for (y = 0; y < iOrder; y++)
//	{
//		iMax = y;
//		fMax = Ai[Q[iMax] * iRow_Size + y];
//		for (i = y + 1; i < iOrder; i++)
//		{//Ѱ������Ԫ
//			if (abs(Ai[iPos = Q[i] * iRow_Size + y]) > abs(fMax))
//			{
//				fMax = Ai[iPos];
//				iMax = i;
//			}
//		}
//		if (abs(fMax) < ZERO_APPROCIATE)
//		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
//			printf("������,����ԪΪ��%f\n", fMax);
//			*pbSuccess = 0;
//			free(Ai);
//			return;
//		}
//
//		//�����ԪSWAP��Q�ĵ�ǰλ����
//		iTemp = Q[y];
//		Q[y] = Q[iMax];
//		Q[iMax] = iTemp;
//
//		//��iMax���ڵ��н���ϵ�����㣬��ϵ��/=A[y][y]
//		pfMax_Row = &Ai[Q[y] * iRow_Size];
//		pfMax_Row[y] = 1.f;
//		for (x = y + 1; x < iRow_Size; x++)
//			pfMax_Row[x] /= fMax;
//
//		//Disp(Ai, iOrder, iOrder + 1, "\n");
//
//		//�Ժ��������д���
//		for (i = y + 1; i < iOrder; i++)
//		{//i��ʾ��i��
//			iPos = Q[i] * iRow_Size;
//			if ((fValue = Ai[iPos + y]) != 0)
//			{//���ڶ�ӦԪ��Ϊ0�����������
//				for (x = y + 1; x < iRow_Size; x++)
//					Ai[iPos + x] -= fValue * pfMax_Row[x];
//				Ai[iPos + y] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
//			}
//			//Disp(Ai, iOrder, iOrder + 1, "\n");
//		}
//	}
//
//	//��һ����
//	X[iOrder - 1] = Ai[Q[iOrder - 1] * iRow_Size + iOrder];
//
//	//�ش�����Q[iOrder-1]��ʼ�ش���������һ�����ϻش�
//	for (y = iOrder - 2; y >= 0; y--)
//	{
//		//�������ϻش�
//		iPos = Q[y] * iRow_Size;
//		//fValue=b
//		fValue = Ai[iPos + iOrder];
//		for (x = y + 1; x < iOrder; x++)
//		{
//			fValue -= Ai[iPos + x] * X[x];
//			Ai[iPos + x] = 0;		//�˴����Ǳ���ģ�������0���ÿ�һЩ����
//		}
//		X[y] = fValue;
//		Ai[iPos + iOrder] = fValue;	//�˴����Ǳ��룬�ÿ�����
//	}
//
//	//Disp(Ai, iOrder, iOrder + 1, "\n");
//	*pbSuccess = 1;
//	free(Ai);
//	//free(X1);
//
//	//����
//	//Linear_Equation_Check(A, iOrder, B, X, (_T)ZERO_APPROCIATE);
//	return;
//}

//void Solve_Linear_Gause(float* A, int iOrder, float* B, float* X, int* pbSuccess)
//{//�ø�˹����Ԫ��������Է�����, Ҫ�㣺
//	//1����˹����ȫ�ȼ������б任��ֻ����û��������Ĺ������������ǲ������ǽ���Ԫ��Ϊ1
//	//2��ѡ����Ԫ��ԭ���Ǳ�֤�� ����ϵ��/aijʱ��ĸ������̫С����ĸС����
//	//3,������Ԫ̫С��<eps)��������㲻���ȣ��˳���ʵ������Ԫ̫С�ᵼ�º���ĳ����������
//	int y, x, i, iRow_Size;
//	int iMax, iTemp, iPos, Q[256];
//	float fMax, * pfMax_Row, fValue;
//	if (iOrder > 256)
//	{
//		printf("Too large order:%d\n", iOrder);
//		*pbSuccess = 0;
//		return;
//	}
//
//	float* Ai = (float*)malloc((iOrder + 1) * iOrder * sizeof(float));
//	iPos = 0;
//	for (y = 0; y < iOrder; y++)
//	{
//		for (x = 0; x < iOrder; x++, iPos++)
//			Ai[iPos] = A[y * iOrder + x];
//		Ai[iPos++] = B[y];
//		Q[y] = y;	//ÿ����Ԫ���ڵ���
//	}
//	//Disp(Ai, iOrder, iOrder + 1,"\n");
//	iRow_Size = iOrder + 1;
//	//Ϊ�˱�����⣬����iMax��ʾQ�е�������������Ai�е��к�
//	for (y = 0; y < iOrder; y++)
//	{
//		iMax = y;
//		fMax = Ai[Q[iMax] * iRow_Size + y];
//		for (i = y + 1; i < iOrder; i++)
//		{//Ѱ������Ԫ
//			if (abs(Ai[iPos = Q[i] * iRow_Size + y]) > abs(fMax))
//			{
//				fMax = Ai[iPos];
//				iMax = i;
//			}
//		}
//		if (abs(fMax) < ZERO_APPROCIATE)
//		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
//			printf("������,����ԪΪ��%f\n", fMax);
//			*pbSuccess = 0;
//			free(Ai);
//			return;
//		}
//
//		//�����ԪSWAP��Q�ĵ�ǰλ����
//		iTemp = Q[y];
//		Q[y] = Q[iMax];
//		Q[iMax] = iTemp;
//
//		//��iMax���ڵ��н���ϵ�����㣬��ϵ��/=A[y][y]
//		pfMax_Row = &Ai[Q[y] * iRow_Size];
//		pfMax_Row[y] = 1.f;
//		for (x = y + 1; x < iRow_Size; x++)
//			pfMax_Row[x] /= fMax;
//
//		//Disp(Ai, iOrder, iOrder + 1, "\n");
//
//		//�Ժ��������д���
//		for (i = y + 1; i < iOrder; i++)
//		{//i��ʾ��i��
//			iPos = Q[i] * iRow_Size;
//			if ( (fValue = Ai[iPos + y])!=0)
//			{//���ڶ�ӦԪ��Ϊ0�����������
//				for (x = y + 1; x < iRow_Size; x++)
//					Ai[iPos + x] -= fValue * pfMax_Row[x];
//				Ai[iPos + y] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
//			}
//			//Disp(Ai, iOrder, iOrder + 1, "\n");
//		}
//	}
//
//	//��һ����
//	X[iOrder - 1] = Ai[Q[iOrder - 1] * iRow_Size + iOrder];
//
//	//�ش�����Q[iOrder-1]��ʼ�ش���������һ�����ϻش�
//	for (y = iOrder - 2; y >= 0; y--)
//	{
//		//�������ϻش�
//		iPos = Q[y] * iRow_Size;
//		//fValue=b
//		fValue = Ai[iPos + iOrder];
//		for (x = y + 1; x < iOrder; x++)
//		{
//			fValue -= Ai[iPos + x] * X[x];
//			Ai[iPos + x] = 0;		//�˴����Ǳ���ģ�������0���ÿ�һЩ����
//		}
//		X[y] = fValue;
//		Ai[iPos + iOrder] = fValue;	//�˴����Ǳ��룬�ÿ�����
//	}
//
//	//Disp(Ai, iOrder, iOrder + 1, "\n");
//	*pbSuccess = 1;
//	free(Ai);
//	//free(X1);
//
//	//����
//	Linear_Equation_Check(A, iOrder, B, X, ZERO_APPROCIATE);
//	return;
//}
void Solve_Linear_Cramer(float* A, int iOrder, float* B, float* X, int* pbSuccess)
{//�ÿ���ķ��������Է���AX=B���˷��ȼ۸�˹����ֻ�ܽ����η�����Ψһ�⣬�����Աܿ�
//��˹���г���Ϊ0����������������ձ��ԡ���ʱ�临�ӶȺܲ�������ۼ�ֵ
	int i, j;
	float fDet = fGet_Determinant(A, iOrder);
	float* Ai = (float*)malloc(iOrder * iOrder * sizeof(float));
	if (abs(fDet) < ZERO_APPROCIATE)
	{
		free(Ai);
		*pbSuccess = 0;
		return;
	}
	for (i = 0; i < iOrder; i++)
	{//��A�ĵ�i�л��ɳ�������
		memcpy(Ai, A, iOrder * iOrder * sizeof(float));
		for (j = 0; j < iOrder; j++)
			Ai[j * iOrder + i] = B[j];
		X[i] = fGet_Determinant(Ai, iOrder) / fDet;
	}

	free(Ai);
	*pbSuccess = 1;
	return;
}
template<typename _T>_T fDot(_T V0[], _T V1[], int iDim)
{//���ڻ�
	_T fTotal = 0;
	for (int i = 0; i < iDim; i++)
		fTotal += V0[i] * V1[i];
	return fTotal;
}
template<typename _T>void Cross_Product(_T V0[], _T V1[], _T V2[])
{//�����ˣ�������������� V2=V0xV1��ֻ����ά a��b=��aybz-azby)i + (azbx-axbz)j + (axby-aybx)k
	_T Temp[3];
	Temp[0] = V0[1] * V1[2] - V0[2] * V1[1];
	Temp[1] = V0[2] * V1[0] - V0[0] * V1[2];
	Temp[2] = V0[0] * V1[1] - V0[1] * V1[0];
	V2[0] = Temp[0], V2[1] = Temp[1], V2[2] = Temp[2];
	return;
}

//void Schmidt_Orthogon(float* A, int m, int n, float* B)
//{//m��n�У�m������������
//	float fValue, * bi, * pB = (float*)malloc(m * n * sizeof(float));
//	int i, j, k;
//	//Disp(A, m, n);
//	for (i = 0; i < m; i++)
//	{//ÿ����һ��Bi
//		bi = &pB[i * n];
//		for (j = 0; j < n; j++)
//			bi[j] = A[i * n + j];	//bi=ai;
//		for (j = 0; j < i; j++)
//		{
//			fValue = fDot(&pB[j * n], &A[i * n], n);
//			fValue /= fDot(&pB[j * n], &pB[j * n], n);
//			for (k = 0; k < n; k++)
//				bi[k] -= fValue * pB[j * n + k];
//		}
//	}
//
//	//�ٵ�λ��
//	for (i = 0; i < m; i++)
//	{
//		bi = &pB[i * n];
//		Normalize(bi, n, bi);
//		//Disp(bi, 1, n);
//	}
//	memcpy(B, pB, m * n * sizeof(float));
//	free(pB);
//	return;
//}

template<typename _T>void Schmidt_Orthogon(_T* A, int m, int n, _T* B)
{//m��n�У�m������������
	_T fValue, * bi, * pB = (_T*)malloc(m * n * sizeof(_T));
	int i, j, k;
	//Disp(A, m, n);
	for (i = 0; i < m; i++)
	{//ÿ����һ��Bi
		bi = &pB[i * n];
		for (j = 0; j < n; j++)
			bi[j] = A[i * n + j];	//bi=ai;
		for (j = 0; j < i; j++)
		{
			fValue = fDot(&pB[j * n], &A[i * n], n);
			fValue /= fDot(&pB[j * n], &pB[j * n], n);
			for (k = 0; k < n; k++)
				bi[k] -= fValue * pB[j * n + k];
		}
	}

	//�ٵ�λ��
	for (i = 0; i < m; i++)
	{
		bi = &pB[i * n];
		Normalize(bi, n, bi);
		//Disp(bi, 1, n);
	}
	memcpy(B, pB, m * n * sizeof(_T));
	free(pB);
	return;
}

void Schmidt_Orthogon(float* A, int iOrder, float* B)
{//ʩ����������.Ϊ��򵥻������������Ծ�����ʽ���롣����������Ҳ�Ծ�����ʽ���
	float fValue, * bi, * pB = (float*)malloc(iOrder * iOrder * sizeof(float));
	int i, j, k;
	//b1=a1;
	//for (i = 0; i < iOrder; i++)
		//pTemp[i] = A[i];
	for (i = 0; i < iOrder; i++)
	{//ÿ����һ��Bi
		bi = &pB[i * iOrder];
		for (j = 0; j < iOrder; j++)
			bi[j] = A[i * iOrder + j];	//bi=ai;
		for (j = 0; j < i; j++)
		{
			fValue = fDot(&pB[j * iOrder], &A[i * iOrder], iOrder);
			fValue /= fDot(&pB[j * iOrder], &pB[j * iOrder], iOrder);
			for (k = 0; k < iOrder; k++)
				bi[k] -= fValue * pB[j * iOrder + k];
		}
		//Disp(bi, 1, 3);
	}
	//�ٵ�λ��
	for (i = 0; i < iOrder; i++)
	{
		bi = &pB[i * iOrder];
		fValue = sqrt(fDot(bi, bi, iOrder));
		for (j = 0; j < 3; j++)
			bi[j] /= fValue;
	}
	memcpy(B, pB, iOrder * iOrder * sizeof(float));
	free(pB);
	return;
}
template<typename _T>int bIs_Orthogonal(_T* A, int h, int w)
{//�ж�һ�������Ƿ�Ϊ�������������������Ƿ��󣬵��ǿ��Խ�һ���ſ���չ��������Ϊһ�����
	int y, x, iMin;
	_T* At, * S;
	//Light_Ptr oPtr = oMem_Mgr;
	if (w == 0)
		w = h;
	//Malloc_1(oPtr, w * h * sizeof(_T), At);
	At = (_T*)pMalloc(&oMem_Mgr, w * h * sizeof(_T));
	if (!At)
		return 0;
	Matrix_Transpose(A, h, w, At);
	iMin = Min(h, w);
	//Malloc_1(oPtr, iMin * iMin * sizeof(_T), S);
	S = (_T*)pMalloc(&oMem_Mgr, iMin * iMin * sizeof(_T));
	if (w > h)	//����ڸ�
		Matrix_Multiply(A, h, w, At, h, S);
	else
		Matrix_Multiply(At, w, h, A, w, S);
	//Disp(S, iMin, iMin, "S");
	for (y = 0; y < iMin; y++)
	{
		for (x = 0; x < iMin; x++)
		{
			if (y == x)
			{//�Խ���
				if (abs(S[y * iMin + x] - 1) > ZERO_APPROCIATE)
					return 0;
			}
			else
			{
				if (abs(S[y * iMin + x]) > ZERO_APPROCIATE)
					return 0;
			}
		}
	}
	Free(&oMem_Mgr, At);
	Free(&oMem_Mgr, S);
	return 1;
}

template<typename _T>int bIs_R(_T R[9], _T eps)
{
	_T Total[3];
	Total[0] = R[0] * R[1] + R[3] * R[4] + R[6] * R[7];
	Total[1] = R[0] * R[2] + R[3] * R[5] + R[6] * R[8];
	Total[2] = R[1] * R[2] + R[4] * R[5] + R[7] * R[8];
	if (!bIs_Orthogonal(R, 3, 3))
	{
		printf("����������\n");
		return 0;
	}
	if (abs(Total[0]) > eps || abs(Total[1]) > eps || abs(Total[2]) > eps)
	{
		Disp(Total, 1, 3, "Not regid transform");
		return 0;
	}
	return 1;
}

//template<typename _T>int bIs_Orthogonal(_T* A, int na)
//{//�жϾ���A�Ƿ�Ϊ�������жϷ���AA'=E
//	int y, x;
//	//_T* At = (_T*)malloc(na * na * sizeof(_T));
//	Light_Ptr oPtr = oMem_Mgr;
//	_T* At;
//	Malloc_1(oPtr, na * na * sizeof(_T), At);
//	if (!At)
//		return 0;
//
//	//float fTotal;
//	for (y = 0; y < na; y++)
//		for (x = 0; x < na; x++)
//			At[y * na + x] = A[x * na + y];	//ת��
//	Matrix_Multiply(A, na, na, At, na, At);
//
//	//Disp(At, na, na, "AAt");
//	for (y = 0; y < na; y++)
//	{
//		for (x = 0; x < na; x++)
//		{
//			if (y == x)
//			{
//				if (abs(At[y * na + x] - 1) > ZERO_APPROCIATE)
//					return 0;
//			}
//			else if (abs(At[y * na + x]) > ZERO_APPROCIATE)
//				return 0;
//		}
//	}
//
//	/* ʵ��֤�������������ÿһ�л�ÿ�в���Ϊ��λ����
//	for (y = 0; y < na; y++)
//	{
//		for (fTotal = 0, x = 0; x < na; x++)
//			fTotal += A[y * na + x];
//		if (abs(fTotal- 1.f)>ZERO_APPROCIATE)
//			printf("Here");
//	}
//	for (x = 0; x < na; x++)
//	{
//		for (fTotal = 0, y = 0; y < na; y++)
//			fTotal += A[y * na + x];
//		if (abs(fTotal - 1.f) > ZERO_APPROCIATE)
//			printf("Here");
//	}*/
//
//	//Disp(At, na, na, "AAt");
//	//free(At);
//	return 1;
//}

void Get_Householder(float* X, float* Y, int iDim, float* H, int* pbSuccess)
{//���� |X|==|Y|,��ģ��ȵ���������������Hx=y�ľ������
//���������һ���Գƾ�����������,H=H',ͬʱ��
	int i, j;
	float fMod_X, fMod_Y, fDenominator;
	float* u = (float*)malloc(iDim * sizeof(float));
	for (fDenominator = fMod_X = fMod_Y = 0, i = 0; i < iDim; i++)
	{
		fMod_X += X[i] * X[i], fMod_Y += Y[i] * Y[i];
		fDenominator += (X[i] - Y[i]) * (X[i] - Y[i]);
	}
	fDenominator = sqrt(fDenominator);
	if (abs(fMod_X - fMod_Y) > ZERO_APPROCIATE)
	{
		printf("x��y��ģ��һ��,|x|=%f |y|=%f", sqrt(fMod_X), sqrt(fMod_Y));
		*pbSuccess = 0;
		return;
	}

	//u�����˵�λ����
	for (i = 0; i < iDim; i++)
		u[i] = (X[i] - Y[i]) / fDenominator;

	//H=I-2uu', ����u�Ѿ��������������ʹ
	for (i = 0; i < iDim; i++)
		for (j = 0; j < iDim; j++)
			H[i * iDim + j] = (i == j ? 1 : 0) - 2 * u[i] * u[j];
	*pbSuccess = 1;

	//////����
	//float *y1 = (float*)malloc(iDim * sizeof(float));
	//memset(y1, 0, iDim * sizeof(float));
	////for (i = 0; i < iDim; i++)
	////	for (j = 0; j < iDim; j++)
	////		y1[i] += H[i*iDim+j] * X[j];
	//for (i = 0; i < iDim; i++)
	//	for (j = 0; j < iDim; j++)
	//		y1[i] += H[i * iDim + j];

	Disp(H, iDim, iDim, "H");
	//Disp(y1, 1, iDim, "y1");		//ȫ��1�Ŷ�

	free(u);
	return;
}

void Get_Householder(float* X, int iDim, int l, float* Q, int* pbSuccess)
{//�����Ѱ��Hʹ��Hx=y��һ�㣬ֻ����X�� Y��Xֱ���Ƶ�����,ʹ�� Qx=y
	//Q��Ȼ�Ǿ�����󣬵����Ǹ�����ľ������
	//l������������λ�ã��հ涨ΪX�������ٸ�����
	//�ȹ���һ��Y��ʹ��ǰl-1����x��ǰi-1���l�һ����ʹ��|X|=|y|
	float fMod, * Y = (float*)malloc(iDim * sizeof(float));
	int i;
	//for (fMod_Total = 0, i = 0; i < iDim; i++)
		//fMod_Total += X[i] * X[i];
	//if (bFront)
	{
		for (i = 0; i < l; i++)
			Y[i] = X[i];
		for (fMod = 0; i < iDim; i++)
			fMod += X[i] * X[i];
		Y[l] = sqrt(fMod);
		for (i = l + 1; i < iDim; i++)
			Y[i] = 0;
	}
	/*else
	{
		for (fMod=0,i = 0; i < iDim - l; i++)
		{
			fMod += X[i] * X[i];
			Y[i] = 0;
		}
		Y[i] = sqrt(fMod+X[i]*X[i]);
		for (i++; i < iDim; i++)
			Y[i] = X[i];
		Disp(Y, 1, 6);
	}*/

	//��������Y��������Qx=y
	Get_Householder(X, Y, iDim, Q, pbSuccess);
	//���ǣ�Q����һ������ľ������ ��I��H���ɣ���������һ���������
	Matrix_Multiply(Q, iDim, iDim, X, 1, Y);
	Disp(Y, 1, iDim);
	//Disp<float>(Y, 1, iDim);
	free(Y);
	return;
}
void Get_Householder_2(float X[], int iDim, int l, float* Q)
{//��һ��Householder�㷨
	float delta, tau;
	float* V = (float*)malloc((iDim - l) * sizeof(float));
	float* H_Dup = (float*)malloc((iDim - l) * (iDim - l) * sizeof(float));
	int i, j;
	//iH_Dim = iDim - l;

	//���delta
	for (delta = 0, i = l; i < iDim; i++)
		delta += X[i] * X[i];
	delta = sign(X[l]) * sqrt(delta);
	for (i = l, j = 0; i < iDim; i++, j++)
		V[j] = X[i];
	V[0] += delta;
	tau = V[0] * delta;
	//Disp(V, 1, iDim-l);
	Matrix_Multiply((float*)V, iDim - l, 1, (float*)V, iDim - l, (float*)H_Dup);

	/*if (abs(tau) < 2.2204460492503131e-15 && tau != 0)
		printf("Here");*/
	for (i = 0; i < iDim - l; i++)
		for (j = 0; j < iDim - l; j++)
			H_Dup[i * (iDim - l) + j] = (i == j ? 1 : 0) - (tau != 0 ? H_Dup[i * (iDim - l) + j] / tau : 0);
	//if (abs(tau) < 2.2204460492503131e-15 && tau != 0)
		//Disp(H_Dup, iDim - 1, iDim - 1);
	//���Թ���H��
	memset(Q, 0, iDim * iDim * sizeof(float));
	for (i = 0; i < l; i++)
		Q[i * iDim + i] = 1;		//�ȶ�I

	for (i = l; i < iDim; i++)
		for (j = l; j < iDim; j++)
			Q[i * iDim + j] = H_Dup[(i - l) * (iDim - l) + (j - l)];

	//for (i = 0; i < iDim * iDim; i++)
		//Q[i] *= 35;

	//Disp(Q, iDim, iDim);
	free(V);
	free(H_Dup);
	return;
}
void QR_Decompose(float* A, int ma, int na, float* Q, float* R)
{//�򵥷ֽ⣺ A=QnQn-1..Q1*R
	int i, j;
	float* An = (float*)malloc(ma * na * sizeof(float));
	float* v = (float*)malloc(ma * sizeof(float));
	float* Qn = (float*)malloc(ma * ma * sizeof(float));
	float* Q_Dup = (float*)malloc(ma * ma * sizeof(float));

	//�Ƚ�Q��ΪI
	memset(Q_Dup, 0, ma * ma * sizeof(float));
	for (i = 0; i < ma; i++)
		Q_Dup[i * ma + i] = 1;

	memcpy(An, A, ma * na * sizeof(float));
	for (i = 0; i < na; i++)
	{
		for (j = 0; j < ma; j++)
			v[j] = An[j * na + i];
		//if (iCur == 101 && i == 7)
		//{
		//	//Disp(Qn, na, na, "Qn");
		//	//Disp(Q_Dup, ma, na);
		//	Disp(v, ma, 1);
		//	//Disp(Qn, ma, ma);
		//}
		Get_Householder_2(v, ma, i, (float*)Qn);
		//��Ai=Qi-1*Ai-1
		Matrix_Multiply((float*)Qn, ma, ma, (float*)An, na, (float*)An);

		//����Q_Dup*Qn�ǣ��൱������Qn*...Q2*Q1
		Matrix_Multiply((float*)Q_Dup, ma, ma, (float*)Qn, ma, Q_Dup);
		//bIs_Orthogonal(Q_Dup, na);
		//Disp(Q_Dup, na, na, "Q");		
	}
	//��ʱ��R����һ�������Ǿ���
	if(An)
		memcpy(R, An, ma * ma * sizeof(float));
	//Q����һ���������
	memcpy(Q, Q_Dup, ma * ma * sizeof(float));
	//for (i = 0; i < ma * ma; i++)
		//Q[i] *= 21;
	//Disp(Q, ma, ma,"Q");
	//Disp(An, ma, na, "R");
	free(An);
	free(v);
	free(Qn);
	free(Q_Dup);
	return;
}
int bIs_Upper_Tri(float* A, int ma, int na)
{//����һ������A���ж����Ƿ�Ϊ�����Ǿ���
#define eps 2.2204460492503131e-15
	int y, x;
	for (x = 0; x < na; x++)
	{
		for (y = x + 1; y < ma; y++)
		{
			if (abs(A[y * na + x]) >= eps)	// ZERO_APPROCIATE)
				return 0;
		}
	}
	return 1;
#undef eps
}
template<typename _T>int bIs_Unit_Vector(_T* V, int na, _T eps)
{
	int i;
	_T fSum = 0;
	for (i = 0; i < na; i++)
		fSum += V[i] * V[i];
	if (fSum == 0 || abs(fSum - 1) <eps )
		return 1;
	else
		return 0;
}
int bIs_Unit_Vector(float* V, int na, int iStrike)
{
	int i;
	float fTotal;
	for (fTotal = 0, i = 0; i < na; i++)
	{
		fTotal += V[i * iStrike] * V[i * iStrike];
		//printf("%f ", V[i * iStrike]);
	}
	//printf(" Total:%f\n",fTotal);
	if (abs(fTotal - 1) < ZERO_APPROCIATE)
		return 1;
	else
		return 0;
}

void QR_Decompose(float* A, int na, float* R, float* Q, int* pbSuccess, int* pbDup_Root)
{//�������ֵR,��������Q. A=R*Q	
#define MAX_ITERATE_COUNT 200
//#define eps 2.2204460492503131e-15
	float* R_Dup = (float*)malloc(na * na * sizeof(float));
	float* Q_Dup = (float*)malloc(na * na * sizeof(float));
	float* An = (float*)malloc(na * na * sizeof(float));
	//float* Q_Dup_1 = (float*)malloc(na * na * sizeof(float));
	int i;
	memcpy(An, A, na * na * sizeof(float));

	//��Q_Dup_1��ΪI��������������
	//memset(Q_Dup_1, 0, na * na * sizeof(float));
	//for (i = 0; i < na; i++)
		//Q_Dup_1[i * na + i] = 1;

	float* Q1_2_Qn_Product = (float*)malloc(na * na * sizeof(float));
	//*Temp_3 = (float*)malloc(na * na * sizeof(float));
//int bSuccess;
	memset(Q1_2_Qn_Product, 0, na * na * sizeof(float));
	for (i = 0; i < na; i++)
		Q1_2_Qn_Product[i * na + i] = 1;

	i = 0;
	while (1)
	{
		/*for (int j = 0; j < na; j++)
			if (abs(An[j * na + j]) < ZERO_APPROCIATE)
				An[j * na + j] = 0;*/
				//�ֽ�ΪA=QR
		QR_Decompose(An, na, na, Q_Dup, R_Dup);
		//�ֽ��Q��������

		//printf("%d\n", bIs_Orthogonal(Q_Dup, na));
		//Matrix_Multiply(Q_Dup, na, na, R_Dup, na, An);
		//��An=RQ
		Matrix_Multiply(R_Dup, na, na, Q_Dup, na, An);

		//�۳�Qn*...*Q2*Q1,���ս��������������
		Matrix_Multiply(Q1_2_Qn_Product, na, na, Q_Dup, na, Q1_2_Qn_Product);

		////��һ�ѿ��� A= Q*An*Qt
		//Matrix_Multiply(Q1_2_Qn_Product, na, na, An, na, Temp_2);
		//Get_Inv_Matrix(Q1_2_Qn_Product, Temp_3, na, &bSuccess);
		//Matrix_Multiply(Temp_2, na, na, Temp_3,na, Temp_3);
		//Disp(A, na, na, "A");
		//Disp(Q_Dup, na, na, "Q");
		//Disp(Temp_3, na, na, "QAnQt");

		//����An�Ƿ�Ϊ�����Ǿ���
		//Disp(An, na, na, "An");
		if (bIs_Upper_Tri(An, na, na))
			break;
		if (i++ >= MAX_ITERATE_COUNT)
			break;
	}

	////�˴�������һ�����㿴���Ƿ� A= Q*An*Q'
	//float* Temp_1 = (float*)malloc(na * na * sizeof(float));
	//float* Temp_2 = (float*)malloc(na * na * sizeof(float));
	//printf("Orthogonal:%d\n", bIs_Orthogonal(Q1_2_Qn_Product,na));
	//Disp(Q1_2_Qn_Product, na, na, "Q");
	//Disp(An, na, na, "R");
	//Matrix_Multiply(Q1_2_Qn_Product, na, na, An, na, Temp_1);
	//Disp(Temp_1, na, na, "QR");
	//Matrix_Transpose(Q1_2_Qn_Product, na, na, Temp_2);
	//Matrix_Multiply(Temp_1, na, na, Temp_2, na, Temp_2);
	//Disp(A, na, na, "A");
	//Disp(Temp_2, na, na, "QRQt");

	////�����������Ƿ�Ϊ��λ����
	//printf("Is unit vector:\t");
	//for (i = 0; i < na; i++)
	//	printf("%d\t", bIs_Unit_Vector(&Q1_2_Qn_Product[i], na, na));
	//printf("\n");

	memcpy(R, An, na * na * sizeof(float));
	memcpy(Q, Q1_2_Qn_Product, na * na * sizeof(float));

	if (pbDup_Root)
	{
		*pbDup_Root = 0;
		for (i = 1; i < na; i++)
		{
			if (abs(R[(i - 1) * na + i - 1] - R[i * na + i]) < 0.01f)
				//if (abs(R[(i - 1) * na + i - 1] - R[i * na + i]) < eps)
				*pbDup_Root = 1;
		}
	}
	free(Q_Dup);
	free(R_Dup);
	free(An);
	free(Q1_2_Qn_Product);
	*pbSuccess = 1;
	//free(Q_Dup_1);
	return;
#undef MAX_ITERATE_COUNT
#undef eps
}

void svd_2(float* A, int ma, int na, float* U, float* S, float* Vt)
{//��������A'A��������������V, ����AA'��Ӧ�����������
	float* At = (float*)malloc(na * ma * sizeof(float));
	float* AtA = (float*)malloc(na * na * sizeof(float));
	float* V = (float*)malloc(na * na * sizeof(float));
	float* StS = (float*)malloc(na * na * sizeof(float));

	int i, bSuccess;
	memset(V, 0, na * na * sizeof(float));
	Matrix_Transpose(A, ma, na, At);
	Matrix_Multiply(At, na, ma, A, na, AtA);
	//Disp(AtA, na, na, "AtA");
	QR_Decompose(AtA, na, StS, V, &bSuccess);
	Disp(StS, na, na, "StS");
	Disp(V, na, na, "V");

	memset(S, 0, na * na * sizeof(float));
	//�������ɵõ�Sigma������ֵ
	for (i = 0; i < na; i++)
	{
		if (abs(StS[i * na + i]) < ZERO_APPROCIATE)
			StS[i * na + i] = 0;
		S[i * na + i] = sqrt(StS[i * na + i]);
	}

	//Disp(S, na, na, "S");
	//���²��ֿ��������ˣ�Ҫ�Ƶ���U= AV������Ҫ���壡
	//�����Ƶ�һ�Σ� �Ѿ�֪��A'A������������ɵľ���ΪV�� ��������yΪV�е�һ������������ 
	// rΪ��Ӧ������ֵ���� A'A y= ry,  ע�⣬����ҧ��Ҫ�� AA'����������������ǰ�澭�飬 ����ͬʱ���A
	//�� AA'Ay=rAy, ��ʱ�� ��Ay����һ�����壬 ����
	//	AA' ����������ΪAy, ������ֵ��r
	//��ô  U = AV, ֤�ϣ�

	float* U_1 = (float*)malloc(ma * na * sizeof(float));
	memset(U, 0, sizeof(ma * ma));

	//����AV, ��ΪU
	Matrix_Multiply(A, ma, na, V, na, U_1);
	//Ȼ�����ù�񻯣�����U_1����������ɣ��ʴ˱�����н��й��
	Normalize_Col(U_1, ma, na, U_1);
	//Disp(U_1, 9, 9);
	//�����Ժ�U_1�Ǹ�mxn���󣬸����ۻ���U����һ�������з���������
	//���Զ��ڵ�n����m-1�У�������0

	int y, x;
	for (y = 0; y < ma; y++)
	{
		for (x = 0; x < na; x++)
			U[y * ma + x] = U_1[y * na + x];
		for (; x < ma; x++)
			U[y * ma + x] = 0;
	}
	//Disp(U, ma, ma, "U");

	Matrix_Transpose(V, na, na, Vt);
	//Disp(Vt, na, na, "Vt");
	free(At);
	free(AtA);
	free(V);
	free(StS);
	return;
}

void svd_1(float* A, int ma, int na, float* U, float* S, float* Vt)
{//��mxn����A������ֵ�ֽ⣬�ֽ��U(mxm)Ϊ AA'������ֵ��V(nxn)ΪA'A������ֵ, S(mxn)ΪSigma
//��ʱֻ�㵽n>m�����
	float* AAt = (float*)malloc(ma * ma * sizeof(float));
	int y, i, bSuccess;

	Matrix_Multiply_Symmetric(A, ma, na, AAt);
	float* SSt = (float*)malloc(ma * ma * sizeof(float));
	memset(S, 0, ma * na * sizeof(float));
	QR_Decompose(AAt, ma, SSt, U, &bSuccess);
	//�����Ժ� AA' = U x SS' x U', �����￪ʼ���Ͳ��ܿ�������Щǧƪһ�ɵĶ��ѵ�

	//�������ɵõ�Sigma������ֵ
	for (i = 0; i < ma; i++)
		S[i * na + i] = sqrt(SSt[i * ma + i]);

	//����Ҫ��V, Ϊ�˼��Ƶ���ֱ�Ӹ���������֤���� V=A'U�� ��ô��ʱ��������Ҫ����ΪV��������������⣬���ַ����Ǵ��
	//֤��	A'A A'U = A' (AA')U = A' (r U) , ���� rΪ AA'������ֵ�� �����׺��沿��
	//����	A'A (A'U)= r (A'U)�� ��ʱ,�� A'U������������ɵľ�����ô��Ȼ�� (A'U)��Ϊһ�����壬��������������
	//����  V= A'U, ֤��
	//�ټ�һ��ת�� V'= (A'U)'= U'(mxm)A(mxn)

	//�˴�ֱ�������˷�
	int iSize = Max(ma, na) * Max(ma, na);
	float* Vt_1 = (float*)malloc(iSize * sizeof(float));
	memset(Vt_1, 0, iSize * sizeof(float));

	float* Ut = (float*)malloc(ma * ma * sizeof(float));
	Matrix_Transpose(U, ma, ma, Ut);
	Matrix_Multiply(Ut, ma, ma, A, na, Vt_1);

	for (y = 0; y < ma; y++)	//���
		Normalize(&Vt_1[y * na], na, &Vt_1[y * na]);

	//�������ͣ����ö�AtA����һ��QR�ֽ⣬����ò���һ��ä����ָ�⣬����Vt�����һ��
	float* AtA = (float*)malloc(na * na * sizeof(float));
	float* V_2 = (float*)malloc(na * na * sizeof(float)),
		* Vt_2 = (float*)malloc(na * na * sizeof(float));
	float* S_1 = (float*)malloc(na * na * sizeof(float));
	Matrix_Transpose(A, ma, na, AtA);
	Matrix_Multiply(AtA, na, ma, A, na, AtA);
	QR_Decompose(AtA, na, S_1, V_2, &bSuccess);
	Matrix_Transpose(V_2, na, na, Vt_2);
	Normalize(&Vt_2[(na - 1) * na], na, &Vt_1[(na - 1) * na]);
	//for (i = 0; i < na; i++)
		//Vt_1[i * na + na - 1] = Vt_2[(na - 1) * na + i];

	//Disp(Vt_1, na, na, "Vt");
	//Disp(Vt_2, na, na, "V_2");
	free(AtA);
	free(V_2);
	free(Vt_2);
	free(S_1);

	memcpy(Vt, Vt_1, na * na * sizeof(float));
	free(AAt);
	free(Vt_1);
	free(Ut);
	return;
}
void svd(float* A, int ma, int na, float* U, float* S, float* Vt)
{//�ӿڻ����ã��ڴ����˷�
//A = U * S * V' ������ U �� V������������ ������U���ܻ�ܴ󣬹ʴ˲����ж�������

	if (ma < na)
		svd_1(A, ma, na, U, S, Vt);
	else
		svd_2(A, ma, na, U, S, Vt);
}
//void svd(float* A, int ma, int na, float* U, float* S, float* Vt_1)
//{//��Amn��U,V �ֽ�
//	float* AAt = (float*)malloc(ma * ma * sizeof(float));
//	float* AtA = (float*)malloc(na * na * sizeof(float));
//	float fValue;
//	int x, y, i;
//
//	//�ȼ���AAt;AAt�Ǹ��Գƾ��󣬹ʴ˿��Ի���
//	for (y = 0; y < ma; y++)
//	{
//		for (x = 0; x < ma; x++)
//		{
//			//if (y == 0 && x == 1)
//				//printf("here");
//			for (fValue = 0, i = 0; i < na; i++)
//				fValue += A[y * na + i] * A[x * na + i];
//			AAt[y * ma + x] = fValue;
//		}
//	}
//
//	//�ټ���AtA��AtA�Ǹ��Գƾ��󣬹ʴ˿��Ի���
//	for (y = 0; y < na; y++)
//	{
//		for (x = 0; x < na; x++)
//		{
//			for (fValue = 0, i = 0; i < ma; i++)
//				fValue += A[y + i * na] * A[x + i * na];
//			AtA[y * na + x] = fValue;
//		}
//	}
//	//Disp(AtA, na, na,"AtA");
//	//Disp(AAt, ma, ma, "AAt");
//
//	int bSuccess;
//	float* SSt = (float*)malloc(ma * ma * sizeof(float));
//	float* StS = (float*)malloc(na * na * sizeof(float));
//	
//	memset(S, 0, ma * na * sizeof(float));
//	QR_Decompose(AAt, ma, SSt, U, &bSuccess);
//
//	//���·ֽ����ã���Ȼ�õ�������ֵ��SStһ�����ܲ�ס����������һ�������Էϵ�
//	//ԭ����ʹ��ڴ˴�
//	//QR_Decompose(AtA, na, StS, Vt_1, &bSuccess);	//����ֽ�û��
//	//Disp(Vt_1, na, na, "Incorrect V\n");
//	Disp(AAt, ma, ma, "AAt");
//
//	//��VҪ�ر���գ�
//	float* At = (float*)malloc(na * ma * sizeof(float));
//	float* Vt = (float*)malloc(na * na * sizeof(float));
//	Matrix_Transpose(A, ma, na, At);
//	Matrix_Multiply(At, na, ma, U, ma, Vt);
//	Matrix_Transpose(Vt, na, ma, Vt_1);
//	Disp(U, ma, ma, "U");
//
//	//��AtU�� ����AtΪ n*m����  UΪmxm����
//	
//	for (i = 0; i < ma; i++)
//		Normalize(&Vt_1[i * na], na, &Vt_1[i * na]);
//	Matrix_Transpose(Vt_1, na, na, Vt_1);
//	//Disp(Vt_1, na, na, "V");
//	Matrix_Transpose(Vt_1, ma, na, Vt);
//
//
//	//Disp(AAt, 3, 3, "AAt");
//	//Disp(SSt, 3, 3, "SSt");
//	//Disp(U, 3, 3, "U");
//	//Disp(AtA, 5, 5, "AtA");
//	//Disp(StS, 5, 5, "StS");
//	//Disp(Vt_1, 5, 5, "V");
//	for (i = 0; i < ma; i++)
//		S[i * na + i] = sqrt(SSt[i * ma + i]);
//	//for (i = 0; i < ma; i++)
//		//S[i * na + i] = sqrt(StS[i * na + i]);
//	
//
//	//Disp(U, ma, ma, "U");
//	//Disp(Vt_1, na, na, "V");
//	//printf("%d\n", bIs_Orthogonal((float*)V,na));
//	//Disp(SSt, ma, ma, "SSt");
//	//Disp(StS, na, na, "StS");
//	//Disp(S, ma, na, "Sigma");
//
//	//����һ��
//	float* Temp_1 = (float*)malloc(ma * na* sizeof(float));
//	Matrix_Transpose(Vt_1, na, na, Vt);
//	Matrix_Multiply(U, ma, ma, S, na, Temp_1);
//	//Disp(Temp_1, ma, na, "UxS");
//	Matrix_Multiply(Temp_1, ma, na, Vt, na, Temp_1);
//	
//	//Disp(Vt, na, na, "Vt");
//	//Disp(Temp_1, ma, na, "USVt");
//
//	memcpy(Vt_1, Vt, na * na * sizeof(float));
//
//	free(AAt);
//	free(AtA);
//	free(SSt);
//	free(StS);
//	free(Vt);
//	free(At);
//}

void Conjugate_Gradient(float* A, const int n, float B[], float X[])
{//�����ݶȷ������Է���Ax=b
#define eps 0.0001
	float* pBuffer, * r0, r0_dot_r0,	//r0.r0���ڻ�
		r1_dot_r1,	//r1.r1, �ڻ�
		* r1, * p0,/* * p1,*/* Aux, alpha, beta, fValue_0;

	int i, iCount = 0;
	pBuffer = (float*)malloc(n * 5 * sizeof(float));
	r0 = pBuffer;
	r1 = pBuffer + n;
	p0 = r1 + n;
	//p1 = p0 + n;
	Aux = r0 + n;

	for (i = 0; i < n; i++)
		X[i] = 0;		//��ʼx0Ϊ[0, 0, 0]

	//p0=r0= b-Ax0
	Matrix_Multiply(A, n, n, X, 1, Aux);
	for (i = 0; i < n; i++)
		p0[i] = r0[i] = B[i] - Aux[i];

	while (1)
	{
		r0_dot_r0 = fDot(r0, r0, n);

		//ak= rk'rk / pk'A pk;
		Matrix_Multiply(p0, 1, n, A, n, Aux);
		for (fValue_0 = 0, i = 0; i < n; i++)
			fValue_0 += Aux[i] * p0[i];
		alpha = r0_dot_r0 / fValue_0;

		//xk+1= xk+ak*pk;
		for (i = 0; i < n; i++)
			X[i] += alpha * p0[i];

		//rk+1= rk-ak.A.pk
		Matrix_Multiply(A, n, n, p0, 1, Aux);
		for (r1_dot_r1 = 0, i = 0; i < n; i++)
		{
			r1[i] = r0[i] - alpha * Aux[i];
			r1_dot_r1 += r1[i] * r1[i];
		}

		fValue_0 = fDot(r1, r1, n);
		if (fValue_0 < eps)
			return;
		printf("iterate:%d %f\n", iCount, fValue_0);

		//beta= r1.r1/r0.r0
		beta = r1_dot_r1 / r0_dot_r0;

		for (i = 0; i < n; i++)
		{
			p0[i] = r1[i] + beta * p0[i];	//pk+1= r1 + beta*pk
			r0[i] = r1[i];					//r0=r1
		}
		iCount++;
	}

	free(pBuffer);
	return;
#undef eps
}
template<typename _T>void Free_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix)
{
	/*if (poMatrix->m_pRow)
		Free(&oMem_Mgr,poMatrix->m_pRow);*/
	//if (poMatrix->m_pBuffer)
		//Free(&oMem_Mgr, poMatrix->m_pBuffer + 1);
	if (poMatrix->m_pRow)
		Free(&oMem_Mgr, poMatrix->m_pRow);
}

template<typename _T>void Compact_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix)
{//һ�����ں���Դ˾���ֻ�����������Item��ռ�Ŀռ�黹
	/*unsigned int* pNew_Addr = poMatrix->m_pRow + poMatrix->m_iRow_Count;
	memmove(pNew_Addr, poMatrix->m_pCol, poMatrix->m_iCol_Count * sizeof(unsigned int));
	poMatrix->m_pCol = pNew_Addr;
	poMatrix->m_pRow = (unsigned int*)realloc(poMatrix->m_pRow, (poMatrix->m_iRow_Count + poMatrix->m_iCol_Count) * sizeof(unsigned int));*/
	poMatrix->m_iMax_Item_Count = poMatrix->m_iCur_Item;
	int iSize = Max(poMatrix->m_iCol_Count,poMatrix->m_iRow_Count) * 2 * sizeof(unsigned int) + (poMatrix->m_iMax_Item_Count+1) * sizeof(typename Sparse_Matrix<_T>::Item);
	Shrink(&oMem_Mgr, poMatrix->m_pRow, iSize);
	return;
}
template<typename _T>void Disp_Link_Col(Sparse_Matrix<_T> oMatrix, int x)
{
	struct Sparse_Matrix<_T>::Item *poCur;
	printf("x:%d y:", x);
	if (!oMatrix.m_pCol[x])
	{
		printf("Null\n");
		return;
	}
	poCur = &oMatrix.m_pBuffer[oMatrix.m_pCol[x]];
	do
	{
		printf("%d ", poCur->y);
		poCur = poCur->m_iCol_Next ? &oMatrix.m_pBuffer[poCur->m_iCol_Next] : NULL;
	} while (poCur);
	printf("\n");
	return;
}

template<typename _T>void Disp_Link_Row(Sparse_Matrix<_T> oMatrix, int y)
{
	struct Sparse_Matrix<_T>::Item* poCur;
	printf("y:%d x:", y);
	if (!oMatrix.m_pRow[y])
	{
		printf("Null\n");
		return;
	}
	poCur = &oMatrix.m_pBuffer[oMatrix.m_pRow[y]];
	do
	{
		printf("%d ", poCur->x);
		poCur = poCur->m_iRow_Next ? &oMatrix.m_pBuffer[poCur->m_iRow_Next] : NULL;
	} while (poCur);
	printf("\n");
	return;
}

template<typename _T>void Resize_Matrix(Sparse_Matrix<_T>* poA, int iNew_Item_Count)
{
	//��С�ڴ�
	poA->m_pBuffer++;
	poA->m_pBuffer = (struct Sparse_Matrix<_T>::Item*)realloc(poA->m_pBuffer, iNew_Item_Count * sizeof(struct Sparse_Matrix<_T>::Item));
	poA->m_iMax_Item_Count = iNew_Item_Count;
	poA->m_pBuffer--;
}

template<typename _T>void Compare_Transpose(Sparse_Matrix<_T> A, Sparse_Matrix<_T> B)
{//��������ת���Ƿ��
	int x, y;
	struct Sparse_Matrix<_T>::Item oCur_Row, oCur_Col;
	for (y = 0; y < A.m_iRow_Count; y++)
	{
		oCur_Row = A.m_pBuffer[A.m_pRow[y]];
		x = y;
		oCur_Col = B.m_pBuffer[B.m_pCol[x]];
		while (1)
		{
			if (oCur_Row.m_fValue != oCur_Col.m_fValue || oCur_Row.y != oCur_Col.x ||
				oCur_Row.m_iCol_Next != oCur_Col.m_iRow_Next || oCur_Row.m_iRow_Next != oCur_Col.m_iCol_Next)
			{
				printf("err");
				break;
			}
			if (oCur_Row.m_iRow_Next)
				oCur_Row = A.m_pBuffer[oCur_Row.m_iRow_Next];
			else
				break;
			if (oCur_Col.m_iCol_Next)
				oCur_Col = B.m_pBuffer[oCur_Col.m_iCol_Next];
			else
				break;
		}
	}
	return;
}

template<typename _T>void Matrix_Add(Sparse_Matrix<_T>* poA, _T a, Sparse_Matrix<_T>* poB, _T b, _T* pC)
{//�������Ҫ��������Ϊֻ������ABͬһ��̬��̫����
	Sparse_Matrix<_T> A, B;
	struct Sparse_Matrix<_T>::Item oItem;
	_T* pCur;
	int y;
	A = *poA;
	B = *poB;
	memset(pC, 0, A.m_iRow_Count * A.m_iCol_Count * sizeof(_T));
	for (y = 0; y < A.m_iRow_Count; y++)
	{
		if (A.m_pRow[y])
		{
			oItem = A.m_pBuffer[A.m_pRow[y]];
			pCur = &pC[oItem.y * A.m_iCol_Count];
			while (1)
			{
				pCur[oItem.x] = a * oItem.m_fValue;
				if (oItem.m_iRow_Next)
					oItem = A.m_pBuffer[oItem.m_iRow_Next];
				else
					break;
			}
		}

		if (B.m_pRow[y])
		{
			oItem = B.m_pBuffer[B.m_pRow[y]];
			pCur = &pC[oItem.y * B.m_iCol_Count];
			while (1)
			{
				pCur[oItem.x] = b * oItem.m_fValue;
				if (oItem.m_iRow_Next)
					oItem = B.m_pBuffer[oItem.m_iRow_Next];
				else
					break;
			}
		}
	}
	return;
}
template<typename _T>void Matrix_Minus(_T A[], _T B[], int iOrder, _T C[])
{//����ӷ����˴��Ƿ���
	int i;
	for (i = 0; i < iOrder * iOrder; i++)
		C[i] = A[i] - B[i];
	return;
}
template<typename _T>void Matrix_Add(_T A[], _T B[], int iOrder, _T C[])
{//����ӷ����˴��Ƿ���
	int i;
	for (i = 0; i < iOrder * iOrder; i++)
		C[i] = A[i] + B[i];
	return;
}

//template<typename _T>void Matrix_Add(Sparse_Matrix<_T>* poA, _T a, Sparse_Matrix<_T>* poB, _T b)
//{//��ʱδ���á���������ӣ� aA+bB=>B�� �������B��
//	Sparse_Matrix<_T> A = *poA, B = *poB;
//	//pNew_Col_Preʵ������һ��Item, ���Դ���ϴθ��е����һ��Item,�Ա�ӿ��ٶ�
//	Sparse_Matrix<_T>::Item ** pNew_Col_Pre = (Sparse_Matrix<_T>::Item**)malloc(B.m_iCol_Count * sizeof(Sparse_Matrix<_T>::Item*));
//	Sparse_Matrix<_T>::Item oNew, * poNew, * poItem_A, * poItem_B, * poRow_Pre_B;	// , oHead = { };
//	int y, bAdd;
//	memset(pNew_Col_Pre, 0, B.m_iCol_Count * sizeof(Sparse_Matrix<_T>::Item*));
//	for (y = 0; y < A.m_iRow_Count; y++)
//	{
//		poItem_A = &A.m_pBuffer[A.m_pRow[y]];
//		poRow_Pre_B = poItem_B = &B.m_pBuffer[B.m_pRow[y]];
//
//		while (1)
//		{
//			bAdd = 0;
//			if (poItem_A->x == poItem_B->x)
//			{//���������
//				poItem_B->m_fValue = a * poItem_A->m_fValue + b * poItem_B->m_fValue;
//				if (poItem_A->m_iRow_Next)
//					poItem_A = &A.m_pBuffer[poItem_A->m_iRow_Next];
//				else
//					break;
//				poRow_Pre_B = poItem_B;
//				if (poItem_B->m_iRow_Next)
//					poItem_B = &B.m_pBuffer[poItem_B->m_iRow_Next];
//			}
//			else if (poItem_A->x < poItem_B->x)
//			{//��new_Item=aA[]���뵽B��������
//				oNew = *poItem_A;
//				oNew.m_iRow_Next = (int)(poItem_B - B.m_pBuffer);
//				B.m_pBuffer[++B.m_iMax_Item_Count] = oNew;
//				bAdd = 1;
//				if (poRow_Pre_B != poItem_B)
//				{//������ͷ
//					poRow_Pre_B->m_iRow_Next = B.m_iMax_Item_Count;
//				}
//				else
//				{//����ͷ
//					B.m_pRow[y] = B.m_iMax_Item_Count;
//				}
//				poRow_Pre_B = &B.m_pBuffer[B.m_iMax_Item_Count];
//				//poItem_A���ƽ�
//				if (!poItem_A->m_iRow_Next)
//					break;
//				else
//					poItem_A = &A.m_pBuffer[poItem_A->m_iRow_Next];
//			}
//			else
//			{//poItem_A->x > poItem_B->x, poItem_B����ƽ�һ��
//				poRow_Pre_B = poItem_B;
//				if (!poItem_B->m_iRow_Next)
//				{//��Ҫ������Ԫ��
//					oNew = *poItem_A;
//					oNew.m_iRow_Next = 0;
//					bAdd = 1;
//					B.m_pBuffer[++B.m_iMax_Item_Count] = oNew;
//					poItem_B->m_iRow_Next = B.m_iMax_Item_Count;
//					poItem_B = &B.m_pBuffer[B.m_iMax_Item_Count];
//					if (poItem_A->m_iRow_Next)
//						poItem_A = &A.m_pBuffer[poItem_A->m_iRow_Next];
//					else
//						break;
//				}
//				else
//					poItem_B = &B.m_pBuffer[poItem_B->m_iRow_Next];
//			}
//			if (bAdd)
//			{//�Ѿ����룬�ٸ��з���
//				Sparse_Matrix<_T>::Item* poCur, * poPrevious;
//				poNew = &B.m_pBuffer[B.m_iMax_Item_Count];
//				if (!pNew_Col_Pre[oNew.x])
//				{
//					if (!(B.m_pCol[oNew.x]))
//					{	//Ͱ��û��ָ��ֱ�Ӳ���
//						B.m_pCol[oNew.x] = B.m_iMax_Item_Count;
//					}
//					else
//					{
//						poCur = &B.m_pBuffer[B.m_pCol[oNew.x]];
//						poPrevious = NULL;
//						while (poCur->y < oNew.y && poCur->m_iCol_Next)
//						{
//							poPrevious = poCur;
//							poCur = &B.m_pBuffer[poCur->m_iCol_Next];
//						}
//						if (poCur->y < oNew.y)
//						{//�嵽ǰ��
//							poNew->m_iCol_Next = poCur->m_iCol_Next;
//							poCur->m_iCol_Next = B.m_iMax_Item_Count;
//						}
//						else
//						{
//							if (poPrevious)
//							{//����previous��cur֮��
//								poNew->m_iCol_Next = poPrevious->m_iCol_Next;
//								poPrevious->m_iCol_Next = B.m_iMax_Item_Count;
//							}
//							else
//							{//�鵽����ͷ
//								poNew->m_iCol_Next = B.m_pCol[oNew.x];
//								B.m_pCol[oNew.x] = B.m_iMax_Item_Count;
//							}
//						}
//					}
//				}
//				else
//				{//һֱ�����ң�ֱ���ҵ�����λ��
//					poCur = pNew_Col_Pre[oNew.x];
//					poPrevious = NULL;
//
//					while (poCur->y < oNew.y && poCur->m_iCol_Next)
//					{
//						poPrevious = poCur;
//						poCur = &B.m_pBuffer[poCur->m_iCol_Next];
//					}
//					if (poCur->y < oNew.y)
//					{//�嵽poCur֮��
//						poNew->m_iCol_Next = poCur->m_iCol_Next;
//						poCur->m_iCol_Next = B.m_iMax_Item_Count;
//					}
//					else
//					{
//						if (poPrevious)
//						{//����previous��cur֮��
//							poNew->m_iCol_Next = poPrevious->m_iCol_Next;
//							poPrevious->m_iCol_Next = B.m_iMax_Item_Count;
//						}
//					}
//				}
//				pNew_Col_Pre[oNew.x] = &B.m_pBuffer[B.m_iMax_Item_Count];
//			}
//		}
//	}
//	free(pNew_Col_Pre);
//	poB->m_iMax_Item_Count = B.m_iMax_Item_Count;
//	return;
//}
void Solve_abc(float x0, float x1, float x2, float y0, float y1, float y2, float abc[3])
{//�����ֵ�㷨
	float x0_Square = x0 * x0, x1_Square = x1 * x1, x2_Square = x2 * x2;
	float a, b, c;
	float fValue_1 = x1_Square * x2 - x1 * x2_Square - x0_Square * x1 + x0 * x1 * x2;
	a = (y1 * x2 - y2 * x1) / fValue_1;
	b = y0 / x0 - (y2 * x2 * x0 - y2 * x1 * x2) / fValue_1 - (y0 * x1 - y1 * x0) / (x0 * x1 - x0_Square) - (x1 * (y1 * x2 - y2 * x1)) / fValue_1;
	c = (y0 * x1 - y1 * x0) / (x1 - x0) + x0 * x1 * (y0 * x2 - y2 * x1) / fValue_1;
	abc[0] = a;
	abc[1] = b;
	abc[2] = c;
}

unsigned long long iFactorial(unsigned long long iPre_Value, int n)
{
	return n == 0 || n==1 ? 1 : iPre_Value * n;
}

unsigned long long iFactorial(int n)
{//��n��ȫ����,�׳�
	unsigned long long iResult = 1;
	if (n == 0)
		return 1;

	for (int i = 2; i <= n; i++)
		iResult *= i;
	return iResult;
}
unsigned int iGet_Perm(int n, int m)
{//��n��ȡm��������= n!/(n-m)! mСn��
	if (m > n)
		return 0;
	int iEnd = n - m + 1;
	unsigned int iResult = 1;
	while (n >= iEnd)
		iResult *= n--;
	return iResult;
}
unsigned int iGet_Combination(int n, int m)
{//������� n��ȡm = p(n,m)/m!
	if (m > n)
		return 0;
	return (unsigned int)(iGet_Perm(n, m) / iFactorial(m));
}
template<typename _T>_T fGet_Tr(_T M[], int iOrder)
{//��N�׷���ļ������Խ���֮��
	int i;
	_T fTotal = 0;
	for (i = 0; i < iOrder; i++)
		fTotal += M[i * iOrder + i];
	return fTotal;
}
template<typename _T>void Vee(_T M[], _T V[3])
{//���Գƾ�������
	V[0] = M[7];
	V[1] = M[2];
	V[2] = M[3];
	return;
}
template<typename _T>void Lie_Bracket(_T a[3],_T b[3], _T c[3])
{//so3������������ c = [a,b] = (a^b^ - b^a^)v
 //��ʱ�ò���
	_T Hat_a[3 * 3], Hat_b[3 * 3], Temp_1[3 * 3], Temp_2[3 * 3];
	Hat(a, Hat_a);
	Hat(b, Hat_b);
	Matrix_Multiply_3x3(Hat_a, Hat_b, Temp_1);
	Matrix_Multiply_3x3(Hat_b, Hat_a, Temp_2);
	Matrix_Minus(Temp_1, Temp_2, 3, Temp_1);
	Vee(Temp_1, c);	
	//Disp(c, 1, 3, "c");
}
template<typename _T>void Hat(_T V[], _T M[])
{//���ݸ������������췴�Գƾ��󣬸Ļ�������һ��
	if (V == M)
	{
		_T M1[3 * 3];
		M1[0] = M1[4] = M1[8] = 0;
		M1[1] = -V[2], M1[3] = V[2];
		M1[2] = V[1], M1[6] = -V[1];
		M1[5] = -V[0], M1[7] = V[0];
		memcpy(M, M1, 3 * 3 * sizeof(_T));
	}else
	{
		M[0] = M[4] = M[8] = 0;
		M[1] = -V[2], M[3] = V[2];
		M[2] = V[1], M[6] = -V[1];
		M[5] = -V[0], M[7] = V[0];
	}
	
	return;
}

template<typename _T>void Rotation_Matrix_2_Vector_3(_T R[3 * 3], _T V[3])
{//�������ϵĽⷽ�̲����ף���Ϊ������������������������ʹģ��һ��Ҳ�������������
	_T fCos_Theta = (_T)((R[0] + R[4] + R[8] - 1) * 0.5);
	_T fTheta = (_T)acos(fCos_Theta);
	_T fFactor= (_T)(0.5*fTheta/sin(fTheta));

	V[0] = (R[2 * 3 + 1] - R[1 * 3 + 2]) * fFactor;
	V[1] = (R[0 * 3 + 2] - R[2 * 3 + 0]) * fFactor;
	V[2] = (R[1 * 3 + 0] - R[0 * 3 + 1])*fFactor;

	return;
}

template<typename _T>void Rotation_Matrix_2_Vector_4(_T R[3 * 3], _T V[4])
{//�������ϵĽⷽ�̲����ף���Ϊ������������������������ʹģ��һ��Ҳ�������������
	_T fCos_Theta = (_T)((R[0] + R[4] + R[8] - 1) * 0.5);
	V[3] = (_T)acos(fCos_Theta);
	_T fFactor = (_T)(0.5 / sin(V[3]));

	V[0] = (R[2 * 3 + 1] - R[1 * 3 + 2]) * fFactor;
	V[1] = (R[0 * 3 + 2] - R[2 * 3 + 0]) * fFactor;
	V[2] = (R[1 * 3 + 0] - R[0 * 3 + 1])*fFactor;

	return;
}

template<typename _T>void Rotation_Matrix_2_Vector(_T R[3 * 3], _T V[4])
{//����ת������ת�������ǽ� Rn=n������n���Ǵ����ת�ᣬ��Ȼ����ֵΪ1�� �����������
//�������������
	//��������ֵ=1�� ����(A-rI)x=0, ���x����������������r=1,���Խ�(A-I)x=0����
	//_T I[3][3] = { 1,0,0,0,1,0,0,0,1 };
	printf("tended to be obsolete\n");

	_T R_1[3][3];// B[3] = { 0 }, V_1[3 * 3]
	_T B[3] = { 0 };
	int iResult;
	if (R == V)
	{
		printf("R connot be V\n");
		return;
	}

	//int y, x, iPos;
	////R_1=R-I
	//for (iPos = 0, y = 0; y < 3; y++)
	//	for (x = 0; x < 3; x++, iPos++)
	//		R_1[y][x] = R[iPos];
	////R_1[y][x] = y == x ? R[iPos] - 1.f: R[iPos];

	////Ȼ�����η��� R_1x=0
	//Solve_Homo_3x3(R_1, 1, V_1);

	memcpy(R_1, R, 9 * sizeof(_T));
	R_1[0][0] -= 1.f, R_1[1][1] -= 1.f, R_1[2][2] -= 1.f;
	
	//����Ӧ����svd�ֽ�
	SVD_Info oSVD;
	SVD_Alloc<_T>(3, 3, &oSVD);
	svd_3((_T*)R_1, oSVD, &iResult);

	//Vt�����һ�о��ǽ�
	memcpy(V, &((_T*)oSVD.Vt)[6], 3 * sizeof(_T));
	Free_SVD(&oSVD);
	//Matrix_Multiply((_T*)R_1, 3, 3, V, 1, V);	//������һ���Ƿ�Ax=0
		
	//����ˣ��˴���Ȼò�ƽ� Ax=0��׼ȷ��˵Ӧ���Ǹ���Լ�������⣬|x|=1
	//Solve_Linear_Solution_Construction((_T*)R_1, 3, 3, B, &iResult, V_1);
	/*V[0] = (_T)V_1[0];
	V[1] = (_T)V_1[1];
	V[2] = (_T)V_1[2];*/

	//������ת�Ƕȣ��о������������У����������廹����֤
	_T fTr = fGet_Tr(R, 3);
	const _T eps = (_T)1e-5;

	_T fTemp = (fTr - 1.f) / 2.f;
	fTemp = Clip3(-1.f, 1.f, fTemp);
	fTemp= -acos(fTemp);
	V[3] = fTemp;

	//���и�������δŪ��������ת���������������⡣������Ϊ�������������ɷ�����Ϊ
	//���������������ⳣ��������ԭ����������������ʴ�����������������Ƿ�Ӱ�����
	//����⣬�д��
	
	//if (abs(fTr - 3.f) <= eps)
	//{
	//	//printf("*****%f*****\n", fTr);
	//	V[3] = 0;
	//}
	//else
	//	V[3] = -acos((fTr - 1.f) / 2.f);
	
	
	return;
}
template<typename _T>void Rotation_Vector_2_Quaternion(_T V[4], _T Q[4])
{//��ת����ת��Ϊ��Ԫ�飬�������ն��壬һ����ת������һ����׼��������Ϊ��ת����һ����ת�Ƕȹ���
	_T fSin_Theta_Div_2;
	_T V_1[3];
	_T fTheta = V[3];
	Normalize(V, 3, V_1);
	Q[0] = cos(fTheta / 2.f);
	fSin_Theta_Div_2 = sin(fTheta / 2.f);
	Q[1] = V_1[0] * fSin_Theta_Div_2;
	Q[2] = V_1[1] * fSin_Theta_Div_2;
	Q[3] = V_1[2] * fSin_Theta_Div_2;
	return;
}
template<typename _T>void Rotation_Matrix_2_Quaternion(_T R[], _T Q[])
{//��ת������Ԫ�顣�˴�����ȫʵ�֣�����ȱ��ֱ���㷨���ʴ˿�һ����ת������Ϊ�м��̵��ڹ�ȥ
	_T V[4];
	Rotation_Matrix_2_Vector(R, V);
	Rotation_Vector_2_Quaternion(V, Q);
	return;
}
template<typename _T>void Rotation_Vector_3_2_Matrix(_T V[3], _T R[3 * 3])
{
	_T V1[4];
	Rotation_Vector_3_2_4(V, V1);
	Rotation_Vector_4_2_Matrix(V1, R);
}
template<typename _T>void Rotation_Vector_4_2_Matrix(_T V[4], _T R[3 * 3])
{//��ת��������ת������һ�ѿ���׼��׼, 3���Ѿ�һ��һ��
//������ǰ��������Ϊ��׼����ת�ᣬ���һ������Ϊ��ת��
//�ϸ�����ϵ��xyz��x����y��Զ,z����
	//��һ����ʽ��R = exp(V^)
	//�ڶ�����ʽ���޵����˹��ʽ���ܿ���������˴��õڶ���
	_T fCos_Theta, fSin_Theta;
	_T fTheta = V[3];
	_T V_1[3];
	//fTheta��Ҫ������ģ��
	fCos_Theta = cos(fTheta);
	fSin_Theta = sin(fTheta);
	_T nnt[3][3], I[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	_T Skew_Sym[3][3];

	//������������һ��
	Normalize(V, 3, V_1);

	Matrix_Multiply(V_1, 3, 1, V_1, 3, (_T*)nnt);
	Scale_Matrix_1(I, fCos_Theta);
	Scale_Matrix_1(nnt, (1 - fCos_Theta));
	Hat(V_1, (_T*)Skew_Sym);
	Scale_Matrix_1(Skew_Sym, fSin_Theta);

	Matrix_Add((_T*)I, (_T*)nnt, 3, (_T*)R);
	Matrix_Add((_T*)R, (_T*)Skew_Sym, 3, (_T*)R);

	return;
}

void Quaternion_Add(float Q_1[], float Q_2[], float Q_3[])
{
	for (int i = 0; i < 4; i++)
		Q_3[i] = Q_1[i] + Q_2[i];
	return;
}
void Quaternion_Minus(float Q_1[], float Q_2[], float Q_3[])
{
	for (int i = 0; i < 4; i++)
		Q_3[i] = Q_1[i] - Q_2[i];
	return;
}
void Quaternion_Conj(float Q_1[], float Q_2[])
{//���������
	Q_2[0] = Q_1[0];
	Q_2[1] = -Q_1[1];
	Q_2[2] = -Q_1[2];
	Q_2[3] = -Q_1[3];
}
void Quaternion_Multiply(float Q_1[], float Q_2[], float Q_3[])
{//�˷��Ȳ��ǵ��Ҳ������������䶨��
	Q_3[0] = Q_1[0] * Q_2[0] - Q_1[1] * Q_2[1] - Q_1[2] * Q_2[2] - Q_1[3] * Q_2[3];
	Q_3[1] = Q_1[0] * Q_2[1] + Q_1[1] * Q_2[0] + Q_1[2] * Q_2[3] - Q_1[3] * Q_2[2];
	Q_3[2] = Q_1[0] * Q_2[2] - Q_1[1] * Q_2[3] + Q_1[2] * Q_2[0] + Q_1[3] * Q_2[1];
	Q_3[3] = Q_1[0] * Q_2[3] + Q_1[1] * Q_2[2] - Q_1[2] * Q_2[1] + Q_1[3] * Q_2[0];
	return;
}
void Quaternion_Inv(float Q_1[], float Q_2[])
{//����Ԫ������
	float fMod = fGet_Mod(Q_1, 4);
	int i;
	Quaternion_Conj(Q_1, Q_2);
	fMod *= fMod;
	for (i = 0; i < 4; i++)
		Q_2[i] /= fMod;
	return;
}
template<typename _T>void Quaternion_2_Rotation_Matrix(_T Q[4], _T R[])
{//��Ԫ��ת��Ϊ��ת���� R= vv' + s^2*I + 2sv^ + (v^)^2
	_T fValue, M_2[3][3], M_1[3][3] = { {1,0,0},{0,1,0},{0,0,1} };	//��ʱ����	
	//�����vv' ֱ�ӷ�R����
	Matrix_Multiply(&Q[1], 3, 1, &Q[1], 3, R);

	//���� s^2*I
	fValue = Q[0] * Q[0];
	Scale_Matrix_1(M_1, fValue);
	Matrix_Add(R, (_T*)M_1, 3, R);

	//����2sv
	Hat(&Q[1], (_T*)M_1);
	fValue = 2.f * Q[0];
	memcpy(M_2, M_1, 3 * 3 * sizeof(_T));
	Scale_Matrix_1(M_2, fValue);
	Matrix_Add(R, (_T*)M_2, 3, R);

	//����(v^) ^ 2�������Ѿ��㶨��M_1= v^
	Matrix_Multiply((_T*)M_1, 3, 3, (_T*)M_1, 3, (_T*)M_2);
	Matrix_Add(R, (_T*)M_2, 3, R);

	//Disp((float*)R, 3, 3);
	return;
}

template<typename _T>void Quaternion_2_Rotation_Vector(_T Q[4], _T V[4])
{//����Ԫ��ת��Ϊ��ת����,ע���ˣ���Ԫ�������Ǳ�׼������|v|=1��������ת�ǶȾͲ���
	_T fSin_Theta_Div_2;
	V[3] = 2.f * acos(Q[0]);
	fSin_Theta_Div_2 = sin(V[3] / 2.f);
	V[0] = Q[1] / fSin_Theta_Div_2;
	V[1] = Q[2] / fSin_Theta_Div_2;
	V[2] = Q[3] / fSin_Theta_Div_2;
	return;
}

//void Gen_Homo_Matrix(float R[], float t[], float s, float M[])
//{//����һ��SIM3��������ת��ƽ�ƣ����Ź���
//	int y, x;
//	for (y = 0; y < 3; y++)
//		for (x = 0; x < 3; x++)
//			M[y * 4 + x] = s * R[y * 3 + x];
//	M[15] = 1;
//
//	M[3] = t[0];
//	M[7] = t[1];
//	M[11] = t[2];
//	M[12] = M[13] = M[14] = 0;
//	return;
//}

template<typename _T>void Gen_Homo_Matrix_1(_T Rotation_Vector[3], _T t[3], _T c2w[])
{//����ת������λ�ƹ�����α任����
	int y, x;
	_T M1[4 * 4], R1[3 * 3];
	if (Rotation_Vector)
	{
		_T V1[4];
		Normalize(Rotation_Vector, 3, V1);
		V1[3] = fGet_Mod(Rotation_Vector, 3);
		Rotation_Vector_4_2_Matrix(V1, R1);
	}
	else
	{
		_T V1[4] = { 0,0,1,0 };
		Rotation_Vector_4_2_Matrix(V1, R1);
	}

	for (y = 0; y < 3; y++)
		for (x = 0; x < 3; x++)
			M1[y * 4 + x] = R1[y * 3 + x];
	M1[15] = 1;
	if (t)
	{
		M1[3] = t[0];
		M1[7] = t[1];
		M1[11] = t[2];
	}
	else
		M1[3] = M1[7] = M1[11] = 0;
	M1[12] = M1[13] = M1[14] = 0;
	memcpy(c2w, M1, 4 * 4 * sizeof(_T));
	return;
}
template<typename _T>void Gen_Homo_Matrix(_T R[], _T t[], _T c2w[])
{//����ת������λ�����깹��һ��4x4 ��α任���󣬴˴�����ת��ƽ�ƹ���
//����������������Ӧ��������ת��ƽ�ƣ����ԣ��˴��õ�����һ��c2w����
//Ҫ��õ�w2c������c2w������󼴿�
	int y, x;
	_T M1[4 * 4];
	if (R)
	{
		for (y = 0; y < 3; y++)
			for (x = 0; x < 3; x++)
				M1[y * 4 + x] = R[y * 3 + x];
	}
	else
	{
		_T Rotation_Vector[4] = { 0,0,1,0 };
		_T R1[3 * 3];
		Rotation_Vector_4_2_Matrix(Rotation_Vector, R1);
		for (y = 0; y < 3; y++)
			for (x = 0; x < 3; x++)
				M1[y * 4 + x] = R1[y * 3 + x];
	}
	M1[15] = 1;

	if (t)
	{
		M1[3] = t[0];
		M1[7] = t[1];
		M1[11] = t[2];
	}
	else
		M1[3] = M1[7] = M1[11] = 0;

	M1[12] = M1[13] = M1[14] = 0;
	memcpy(c2w, M1, 4 * 4 * sizeof(_T));
	return;
}
void Roataion_Vector_2_Angle_Axis(float V_4[4], float V_3[3])
{
	float V_3_1[3], fTheta = V_4[3];
	//��4ά��ת��������Ϊ3ά�����᱾��Ϊ��λ���������Ƶ�ģ��Ϊ�Ƕ�
	//�ܼ򵥣���ÿ������*theta����
	V_3_1[0] = V_4[0] * fTheta;
	V_3_1[1] = V_4[1] * fTheta;
	V_3_1[2] = V_4[2] * fTheta;
	V_3[0] = V_3_1[0];
	V_3[1] = V_3_1[1];
	V_3[2] = V_3_1[2];
	return;
}
template<typename _T>void Get_Jr_4(_T Rotation_Vector_4[4], _T J[])
{//�Ŷ����ұ�
	//J[3] = -J[3];
	J[0] = -J[0], J[1] = -J[1], J[2] = -J[2];
	Get_Jl_4(Rotation_Vector_4, J);
}
template<typename _T>void Get_Jr_3(_T Rotation_Vector[3], _T J[])
{//�Ŷ����ұ�
 //J[3] = -J[3];
	_T v[4];
	Rotation_Vector_3_2_4(Rotation_Vector, v);
	v[0] = -v[0], v[1] = -v[1], v[2] = -v[2];
	Get_Jl_4(v, J);
}
template<typename _T>void Get_Jl_4(_T Rotation_Vector_4[4], _T J[])
{//��������ת���������J������ת����Ϊ4Ԫ��
	//����λ��t,�������J��aΪ��ת������ת��
	_T fValue, Temp_1[3][3], I[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	int i;
	memset(J, 0, 3 * 3 * sizeof(_T));

	//����J��һ���� (sin(theta)/theta)*I
	if (Rotation_Vector_4[3] != 0)
		fValue = sin(Rotation_Vector_4[3]) / Rotation_Vector_4[3];
	else
		fValue = 0;
	for (i = 0; i < 9; i++)
		((_T*)J)[i] = fValue * ((_T*)I)[i];

	//����J�ڶ����� (1-sin(theta)/theta) * axa'
	fValue = 1.f - fValue;
	Matrix_Multiply(Rotation_Vector_4, 3, 1, Rotation_Vector_4, 3, (_T*)Temp_1);
	//Disp((float*)Temp_1, 3, 3);
	for (i = 0; i < 9; i++)
		((_T*)J)[i] += fValue * ((_T*)Temp_1)[i];

	//����������� (1-cos(theta))/theta * a^
	if (Rotation_Vector_4[3] != 0)
		fValue = (1 - cos(Rotation_Vector_4[3])) / Rotation_Vector_4[3];
	else
		fValue = 0;
	Hat(Rotation_Vector_4, (_T*)Temp_1);
	//Disp((float*)Temp_1, 3, 3);
	for (i = 0; i < 9; i++)
		((_T*)J)[i] += fValue * ((_T*)Temp_1)[i];
	//Disp((float*)J, 3, 3);
}
template<typename _T>void Get_Jl_3(_T Rotation_Vector[3], _T J[])
{
	_T v[4];
	Rotation_Vector_3_2_4(Rotation_Vector, v);
	Get_Jl_4(v, J);
}
template<typename _T>void SE3_2_se3(_T T[4 * 4], _T ksi[6])
{//��һ��Homo matrix ת��Ϊksi
	_T R[3 * 3], t[3],V[4];
	Get_R_t(T, R, t);
	Rotation_Matrix_2_Vector(R, V);
	SE3_2_se3(V, t, ksi);
	return;
}

template<typename _T>void SE3_2_se3(_T Rotation_Vector[4], _T t[3], _T Ksi[6])
{//��һ����ת������һ��λ��ת��Ϊse3�ϵ�Ksi
//������SE3��һ��4x4���󣬰�����ת��λ�ơ� SE3->se3����4x4����ת��Ϊ6ά����
//Ȼ������ת������ʾ��ת��ֻ�����࣬�ʴ˴˴�����ת������λ������
//�ܽᣬ SE3�е���ת����se3�е���ά��������������ת��λ��
	_T J[3][3], J_Inv[3][3];
	int bResult;
	Get_Jl_4(Rotation_Vector, (_T*)J);
	//Disp(Rotation_Vector, 1, 4, "Rotation_Vector");
	//Disp((_T*)J, 3, 3, "J");
	//��J�����
	Get_Inv_Matrix_Row_Op_2((_T*)J, (_T*)J_Inv, 3, &bResult);
	if (bResult)
		Matrix_Multiply((_T*)J_Inv, 3, 3, t, 1, Ksi);
	else //����0��ת0�ȣ�t����
		memcpy(Ksi, t, 3 * sizeof(_T));

	//�ٰ�PhiҲ����һ�£���ʵ����Rotation_Vector 4άת��ά
	//Roataion_Vector_2_Angle_Axis(Rotation_Vector, &Ksi[3]);
	Rotation_Vector_4_2_3(Rotation_Vector, &Ksi[3]);
	//Disp(Ksi, 1, 6);	//Rho�������

	return;
}
//void se3_2_SE3(float Ksi[6], float T[])
//{//T��6ά se3����Ksi��Ӧ��4x4����, Ksiǰrho��phi
////ת����Ϻ�T��ȫ��ͼ��ѧ����άת������һ��
////�ܽᣬ SE3�е���ת����se3�е���ά��������������ת��λ��
//
//	//������R
//	float R[3][3];
//	float Rotation_Vector[4];
//
//	Normalize(&Ksi[3], 3, Rotation_Vector);
//	Rotation_Vector[3] = fGet_Mod(&Ksi[3], 3);	//�˴��Ѿ�����ת������Ϊ4ά��ʾ
//
//	Rotation_Vector_2_Matrix(Rotation_Vector, (float*)R);
//	//Disp((float*)R, 3, 3,"R");
//
//	float J[3][3], J_Rho[3];
//	//��Ȼ��J����޹أ�ֻ�Ӧ��Ƶ�����
//	Get_J_by_Rotation_Vector(Rotation_Vector, (float*)J);
//	Matrix_Multiply((float*)J, 3, 3, Ksi, 1, J_Rho);
//
//	//Ȼ�� R,J_Rho, 0', 1��ϳ�T
//	T[0] = R[0][0], T[1] = R[0][1], T[2] = R[0][2], T[3] = J_Rho[0];
//	T[4] = R[1][0], T[5] = R[1][1], T[6] = R[1][2], T[7] = J_Rho[1];
//	T[8] = R[2][0], T[9] = R[2][1], T[10] = R[2][2], T[11] = J_Rho[2];
//	T[12] = T[13] = T[14] = 0, T[15] = 1;
//
//	return;
//}

template<typename _T>void se3_2_SE3(_T Ksi[6], _T T[])
{//T��6ά se3����Ksi��Ӧ��4x4����, Ksiǰrho��phi
//ע�⣺�о����ת���Ǵ���ģ�������J��
//ת����Ϻ�T��ȫ��ͼ��ѧ����άת������һ��
//�ܽᣬ SE3�е���ת����se3�е���ά��������������ת��λ��
	//������R
	_T R[3][3];
	_T Rotation_Vector[4];

	Normalize(&Ksi[3], 3, Rotation_Vector);
	Rotation_Vector[3] = fGet_Mod(&Ksi[3], 3);	//�˴��Ѿ�����ת������Ϊ4ά��ʾ

	Rotation_Vector_4_2_Matrix(Rotation_Vector, (_T*)R);
	//Disp((_T*)R, 3, 3,"R");

	_T J[3][3], J_Rho[3];

	//������ܾ�������
	//��Ȼ��J����޹أ�ֻ�Ӧ��Ƶ�����
	Get_Jl_4(Rotation_Vector, (_T*)J);
	Matrix_Multiply((_T*)J, 3, 3, Ksi, 1, J_Rho);

	//Ȼ�� R,J_Rho, 0', 1��ϳ�T
	T[0] = R[0][0], T[1] = R[0][1], T[2] = R[0][2], T[3] = J_Rho[0];
	T[4] = R[1][0], T[5] = R[1][1], T[6] = R[1][2], T[7] = J_Rho[1];
	T[8] = R[2][0], T[9] = R[2][1], T[10] = R[2][2], T[11] = J_Rho[2];
	T[12] = T[13] = T[14] = 0, T[15] = 1;

	return;
}

//void Gen_Iterate_Matrix(float A[], int n, float B[], float B_1[], float C[])
//{//��Ax=B ��д�� x= Bx+C��BΪ��������cΪϵ������
//	//�ܼ򵥣�����i��/aii����
//	int y, x, iPos, iRow_Size = n + 1;
//	float aii;
//	//�ɴ��B,C�ո������������
//	float* pB_C = (float*)malloc(n * (n + 1) * sizeof(float));
//	//Disp(A, n, n, "\n");
//
//	for (iPos = 0, y = 0; y < n; y++)
//	{
//		aii = 1.f / A[y * n + y];
//		for (x = 0; x < n; x++, iPos++)
//			pB_C[y * iRow_Size + x] = y == x ? 0.f : A[iPos] * aii;
//		pB_C[y * iRow_Size + n] = B[y] * aii;
//		//Disp(pB_C, n, iRow_Size, "\n");
//	}
//	for (iPos = y = 0; y < n; y++)
//	{
//		for (x = 0; x < n; x++, iPos++)
//			B_1[iPos] = pB_C[y * iRow_Size + x];
//		C[y] = pB_C[y * iRow_Size + n];
//	}
//	free(pB_C);
//	return;
//}
template<typename _T>void Gen_Iterate_Matrix(_T A[], int n, _T B[], _T B_1[], _T C[])
{//��Ax=B ��д�� x= Bx+C��BΪ��������cΪϵ������
 //�ܼ򵥣�����i��/aii����
	int y, x, iPos, iRow_Size = n + 1;
	_T aii;
	//�ɴ��B,C�ո������������
	_T* pB_C = (_T*)pMalloc(n * (n + 1) * sizeof(_T));
	//Disp(A, n, n, "\n");

	for (iPos = 0, y = 0; y < n; y++)
	{
		aii = 1.f / A[y * n + y];
		for (x = 0; x < n; x++, iPos++)
			pB_C[y * iRow_Size + x] = y == x ? 0.f : A[iPos] * aii;
		pB_C[y * iRow_Size + n] = B[y] * aii;
		//Disp(pB_C, n, iRow_Size, "\n");
	}
	for (iPos = y = 0; y < n; y++)
	{
		for (x = 0; x < n; x++, iPos++)
			B_1[iPos] = pB_C[y * iRow_Size + x];
		C[y] = pB_C[y * iRow_Size + n];
	}
	Free(pB_C);
	return;
}
template<typename _T>int bIs_Contrative_Mapping(_T B[], int n)
{//����һ����������ͨ���������ļ��ֶ����ж����Ƿ�ѹ������
	_T fF, fRow_Max = -1.f, //�ж���
		fCol_Max = -1.f,		//�ж���
		fTotal;
	int y, x, iPos;
	fF = 0;	//F����
	for (iPos = y = 0; y < n; y++)
	{
		fTotal = 0;
		for (x = 0; x < n; x++, iPos++)
		{
			fTotal += abs(B[iPos]);
			fF += B[iPos] * B[iPos];
		}
		if (fTotal > fRow_Max)
			fRow_Max = fTotal;

	}
	fF = sqrt(fF);

	for (x = 0; x < n; x++)
	{
		fTotal = 0;
		for (y = 0; y < n; y++)
			fTotal += B[y * n + x];
		if (fTotal > fCol_Max)
			fCol_Max = fTotal;
	}
	if (fRow_Max < 1.f || fCol_Max < 1.f || fF < 1.f)
		return 1;	//������������ֻҪ��һ��С��1��BΪѹ��ӳ��
	else
		return 0;
}

template<typename _T>int bIs_Diagonal_Dominant(_T A[], int iOrder)
{//�ж�һ�������Ƿ�Ϊ�Խ���ռ�Ŷ�������ռ�ž���
//�ж϶Խ���ռ�ź�������ȸ��죬���Եÿ���������취ȡ���ж����ȵķ���
	int y, x;
	_T* pCur_Line, fSum;
	pCur_Line = A;
	for (y = 0; y < iOrder; y++, pCur_Line += iOrder)
	{
		for (fSum = 0, x = 0; x < iOrder; x++)
			if (y != x)
				fSum += abs(pCur_Line[x]);
		if (fSum > abs(pCur_Line[y]))
			return 0;
	}
	return 1;
}
template<typename _T>void Solve_Linear_Gauss_Seidel(_T A[], _T B[], int iOrder, _T X[], int* pbResult)
{//���ſɱȵ�����һ����ˣ���˼��������ϼ�������Ͷ����һ�����㣬����ȥ�ܺã�ʵ���ϲ���
//���ſɱȵ������������٣���Ϊ��ʹ�����ſɱ�����Ҳ�����Ҳ���ܷ�ɢ�����������ٶ�Ҳ�����ÿ���
#define ITERATE_COUNT 200
	int i;
	_T* pBuffer = (_T*)pMalloc((iOrder * iOrder + iOrder * 3) * sizeof(_T));
	_T* J = pBuffer;
	_T* C = J + iOrder * iOrder;
	_T* Xk_1 = C + iOrder;
	//float* Temp_1 = Xk_1 + iOrder;
	_T fValue;

	*pbResult = 0;
	//J��������ѹ������
	Gen_Iterate_Matrix(A, iOrder, B, J, C);
	if (!pBuffer)
		goto END;
	for (i = 0; i < iOrder; i++)
		Xk_1[i] = X[i] = 0;
	int y, x;
	for (i = 0; i < ITERATE_COUNT; i++)
	{
		for (y = 0; y < iOrder; y++)
		{
			for (fValue = C[y], x = 0; x < iOrder; x++)
				fValue += J[y * iOrder + x] * X[x];
			X[y] = fValue;	//ÿһ�����ϳ�Ϊ�⣬�˴��ǹؼ�
		}
		Disp(X, 1, 3, "Xk");
		if ((fValue = fGet_Distance(X, Xk_1, 3)) < 0.00001)
		{
			*pbResult = 1;
			goto END;
		}
		memcpy(Xk_1, X, 3 * sizeof(float));
	}
	//�ʹ�һ����������ſɱ��ţ��������ȷ
END:
	if (pBuffer)
		Free(pBuffer);
#undef ITERATE_COUNT
}

template<typename _T>void Solve_Linear_Jocabi(_T A[], _T B[], int iOrder, _T X[], int* pbResult)
{//�ſɱȵ����������Է����飬ͬ���о������⣬�����ⲻ��
#define ITERATE_COUNT 200
	int i;
	_T* pBuffer = (_T*)pMalloc((iOrder * iOrder + iOrder * 3) * sizeof(_T));
	_T* J = pBuffer;
	_T* C = J + iOrder * iOrder;
	_T* Xk_1 = C + iOrder;
	_T* Temp_1 = Xk_1 + iOrder;
	_T* Xk = X;
	_T fValue;

	*pbResult = 0;
	//J��������ѹ������
	Gen_Iterate_Matrix(A, iOrder, B, J, C);

	if (Xk_1)
		for (i = 0; i < iOrder; i++)
			Xk_1[i] = 0;
	for (i = 0; i < ITERATE_COUNT; i++)
	{
		Matrix_Multiply(J, 3, 3, Xk_1, 1, Temp_1);
		Vector_Add(Temp_1, C, 3, Xk);

		//Disp(Xk, 1, 3, "Xk");
		//������������飬����Ϊ����ԭʽ�ȽϽ������
		if ((fValue = fGet_Distance(Xk, Xk_1, 3)) < 0.00001)
		{
			*pbResult = 1;
			goto END;
		}
		//memcpy(Xk_1, Xk, 3 * sizeof(float));
		std::swap(Xk, Xk_1);
	}

	//������������ԭ�����ַ�����һ�ֲ�����������һ�����װ뾶��һ����ϵ�������Ƿ�Խ�ռ��
	if (!bIs_Contrative_Mapping(J, iOrder))
		printf("J������ѹ������\n");
	if (!bIs_Diagonal_Dominant(A, iOrder))
		printf("ϵ������ǶԽ�ռ��\n");
	//�װ뾶���ڻ�û�㶨QR�ֽ�ĸ�����ʾ����
END:
	if (Xk != X)
		memcpy(X, Xk, 3 * sizeof(float));
	Free(pBuffer);

#undef  ITERATE_COUNT
}

//void Solve_Linear_Jocabi(float A[], float B[], int iOrder, float X[], int* pbResult)
//{//�ſɱȵ����������Է����飬ͬ���о������⣬�����ⲻ��
//#define ITERATE_COUNT 200
//	int i;
//	float* pBuffer = (float*)malloc((iOrder * iOrder + iOrder * 3) * sizeof(float));
//	float* J = pBuffer;
//	float* C = J + iOrder * iOrder;
//	float* Xk_1 = C + iOrder;
//	float* Temp_1 = Xk_1 + iOrder;
//	float* Xk = X;
//	float fValue;
//
//	*pbResult = 0;
//	//J��������ѹ������
//	Gen_Iterate_Matrix(A, iOrder, B, J, C);
//
//	if (Xk_1)
//		for (i = 0; i < iOrder; i++)
//			Xk_1[i] = 0;
//	for (i = 0; i < ITERATE_COUNT; i++)
//	{
//		Matrix_Multiply((float*)J, 3, 3, Xk_1, 1, Temp_1);
//		Vector_Add(Temp_1, C, 3, Xk);
//
//		Disp(Xk, 1, 3, "Xk");
//		//������������飬����Ϊ����ԭʽ�ȽϽ������
//		if ((fValue = fGet_Distance(Xk, Xk_1, 3)) < 0.00001)
//		{
//			*pbResult = 1;
//			goto END;
//		}
//		//memcpy(Xk_1, Xk, 3 * sizeof(float));
//		std::swap(Xk, Xk_1);
//	}
//
//	//������������ԭ�����ַ�����һ�ֲ�����������һ�����װ뾶��һ����ϵ�������Ƿ�Խ�ռ��
//	if (!bIs_Contrative_Mapping(J, iOrder))
//		printf("J������ѹ������\n");
//	if (!bIs_Diagonal_Dominant(A, iOrder))
//		printf("ϵ������ǶԽ�ռ��\n");
//	//�װ뾶���ڻ�û�㶨QR�ֽ�ĸ�����ʾ����
//END:
//	if (Xk != X)
//		memcpy(X, Xk, 3 * sizeof(float));
//	free(pBuffer);
//
//#undef  ITERATE_COUNT
//}

void SIM3_Get_Js(float s, float sigma, float theta, float Phi[], float* Js)
{
	float Temp_1[3][3], Temp_2[3][3], I[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	float fValue;
	float hat[3][3];
	float a;	//=s*sin(theta)
	float b;	//=s*cos(theta)
	float c;	//=theta*theta + sigma*sigma
	float C;	//=(s-1)/sigma
	int i;

	a = s * sin(theta);
	b = s * cos(theta);
	c = theta * theta + sigma * sigma;

	C = (s - 1) / sigma;	//����
	for (i = 0; i < 9; i++)
		((float*)Temp_1)[i] = ((float*)I)[i] * C;	//����

	//ע�⣬�˴�����a^�� ������phi^
	//����и����ϣ�����ʹ
	//fValue = (sigma * s * sin(theta) + (1 - s * cos(theta)) * theta) / (sigma * sigma + theta * theta);
	fValue = (a * sigma + (1 - b) * theta) / (theta * c);
	Hat(Phi, (float*)hat);
	for (i = 0; i < 9; i++)
		((float*)Temp_1)[i] += fValue * ((float*)hat)[i];

	//fValue = (s - 1) / sigma - ((s * cos(theta) - 1) * sigma + s * sin(theta) * theta) / (sigma * sigma + theta * theta);
	fValue = (C - ((b - 1) * sigma + a * theta) / c) / (theta * theta);
	Matrix_Multiply((float*)hat, 3, 3, (float*)hat, 3, (float*)Temp_2);
	for (i = 0; i < 9; i++)
		((float*)Temp_2)[i] *= fValue;
	Matrix_Add((float*)Temp_1, (float*)Temp_2, 3, Js);
	return;
}
void sim3_2_SIM3(float zeta[7], float Rotation_Vector[4], float t[], float* ps)
{//��7άsim3������ԭ����ת������λ����������������s
//�ܽᣬsim3�а������ֱ任������׼ȷ˵���任��Ҫ����˳��R->s->t
	float s, theta, sigma;
	float Js[3][3];
	//��ת����ûʲô�ø��,��zeta����һ��
	Normalize(&zeta[3], 3, Rotation_Vector);
	Rotation_Vector[3] = theta = Rotation_Vector[3];
	sigma = zeta[6];
	//s����zigma��

	s = exp(zeta[6]);
	if (ps)
		*ps = s;
	//���Ѹ���t, �����Js����
	SIM3_Get_Js(s, sigma, theta, &zeta[3], (float*)Js);
	Disp((float*)Js, 3, 3, "Js");
	Matrix_Multiply((float*)Js, 3, 3, zeta, 1, t);

	return;
}
void SIM3_2_sim3(float Rotation_Vector[], float t[], float s, float zeta[7])
{//����SIM3Ӧ����һ��4x4����Ȼ�������4x4����ʵ������R,t,s���ɡ�
//����sΪ����ϵ����ת��Ϊһ��7ά�������ֱ�λ��rho,  ��תphi��sigma��Ӧs
//�ܽᣬsim3�а������ֱ任������׼ȷ˵���任��Ҫ����˳��R->s->t

	//s=e^sigma�������sigma
	float sigma = log(s);
	float rho[3];	//t= Js * rho. �Ӵ˴����Ƴ�rho
	float Js[3][3], Js_Inv[3][3];
	float Phi[3];
	int bResult;

	Roataion_Vector_2_Angle_Axis(Rotation_Vector, Phi);

	SIM3_Get_Js(s, sigma, Rotation_Vector[3], Phi, (float*)Js);
	Disp((float*)Js, 3, 3, "Js");
	//�����þ�������ķ���Ҳûë������������Ӧ������Sophus����ѧ�ķ�����ֱ��
	//Ӧ�������㷨
	Get_Inv_Matrix_Row_Op((float*)Js, (float*)Js_Inv, 3, &bResult);

	//Disp((float*)Js_Inv, 3, 3, "Js");
	Matrix_Multiply((float*)Js_Inv, 3, 3, t, 1, rho);
	zeta[0] = rho[0];
	zeta[1] = rho[1];
	zeta[2] = rho[2];
	zeta[3] = Phi[0];
	zeta[4] = Phi[1];
	zeta[5] = Phi[2];
	zeta[6] = sigma;
	Disp(zeta, 1, 7, "zeta");
	return;
}

//��������һ�η���任��ȫ�����б任��һ·���
template<typename _T>void Gen_Rotation_Matrix(_T Axis[3], _T fTheta, _T T[])
{//������ת������ת�Ƕ�����һ����ת���󣬴˴�T��4x4����
//��Rotation_Vector_2_Matrix_3D�ȼۣ����ǿ϶���ܶ࣬��Ϊ����Ҫ�漰����˷�

	_T V_1[3];	//��񻯷�������
	_T fSin_Theta = sin(fTheta),
		fCos_Theta = cos(fTheta),
		f1_Cos_Theta = 1 - fCos_Theta;
	_T fPart_1, fPart_2;

	Normalize(Axis, 3, V_1);

	memset(T, 0, 4 * 4 * sizeof(double));
	//memset(T, 0, 3 * 3 * sizeof(_T));
	T[0] = fCos_Theta + f1_Cos_Theta * V_1[0] * V_1[0];
	T[5] = fCos_Theta + f1_Cos_Theta * V_1[1] * V_1[1];
	T[10] = fCos_Theta + f1_Cos_Theta * V_1[2] * V_1[2];
	T[15] = 1;

	//���öԳ��Լ��㣬Ҳ����
	fPart_1 = f1_Cos_Theta * V_1[0] * V_1[1];
	fPart_2 = fSin_Theta * V_1[2];
	T[1] = fPart_1 - fPart_2;
	T[4] = fPart_1 + fPart_2;

	fPart_1 = f1_Cos_Theta * V_1[0] * V_1[2];
	fPart_2 = fSin_Theta * V_1[1];
	T[2] = fPart_1 + fPart_2;
	T[8] = fPart_1 - fPart_2;

	fPart_1 = f1_Cos_Theta * V_1[1] * V_1[2];
	fPart_2 = fSin_Theta * V_1[0];
	T[6] = fPart_1 - fPart_2;
	T[9] = fPart_1 + fPart_2;
	return;
}
template<typename _T>void Gen_Translation_Matrix(_T Offset[3], _T T[])
{//��ʾһ��λ�ƣ�����һ�����󣬴˾����������λ��߷Ǵ˴�
	memset(T, 0, 4 * 4 * sizeof(_T));
	T[0] = T[5] = T[10] = T[15] = 1;
	T[3] = Offset[0];
	T[7] = Offset[1];
	T[11] = Offset[2];
	return;
}
void Gen_Scale_Matrix(float Scale[3], float T[])
{//���������ϵı����任
	memset(T, 0, 4 * 4 * sizeof(float));
	T[0] = Scale[0];	//Scale_x
	T[5] = Scale[1];	//Scale_y
	T[10] = Scale[2];	//Scale_z
	T[15] = 1;
}
void Rect_2_Polor_Cordinate(float x, float y, float* prho, float* ptheta)
{//2d�µ�ֱ�����굽�����꣬�Ѿ��Ƚ�׼ȷ
	*prho = sqrt(x * x + y * y);
	if (*prho == 0)		//��r=0ʱ��ת���Ѿ�û�����壬������Ϊ0
		*ptheta = 0;
	else if (y > 0)
		*ptheta = acos(x / (*prho));
	else
		*ptheta = -acos(x / (*prho));
	return;
}
void Rect_2_Polor_Coordinate(float x, float y, float z, float* prho, float* ptheta, float* pphi)
{//��ֱ������ϵת��Ϊ�������ʾ���˴���δ���ã����ڸ����������⣬�����ٵ�
//�˴�δ���ƣ���δ��������
	*prho = sqrt(x * x + y * y + z * z);
	*pphi = acos(z / (*prho));
	float theta = atan(y / x);
	if (ptheta)
		*ptheta = theta;
	return;
}
void Screen_2_Coordinate(int x_Screen, int y_Screen, float* px, float* py, int iWidth, int iHeight)
{
	int iWidth_Half = iWidth >> 1,
		iHeight_Half = iHeight >> 1;
	*px = (float)(x_Screen - iWidth_Half);
	*py = (float)(-y_Screen + iHeight_Half);
}
void Rect_2_Screen_Coordinate(float x, float y, int* px_Screen, int* py_Screen, int iWidth, int iHeight)
{//ֱ������ϵ����ת��Ϊ��Ļ����
	int iWidth_Half = iWidth >> 1,
		iHeight_Half = iHeight >> 1;
	*px_Screen = (int)(x + iWidth_Half);
	*py_Screen = (int)(iHeight_Half - y);
	return;
}

void Polor_2_Rect_Coordinate(float rho, float theta, float* px, float* py)
{//��ά�µļ�����ת��Ϊֱ������
	*px = rho * cos(theta);
	*py = rho * sin(theta);
	return;
}
void Polor_2_Rect_Coordinate(float rho, float theta, float phi, float* px, float* py, float* pz)
{//��������(rho, theta, phi)��Ϊֱ������ϵ���˴��ϸ������϶��壬thetaΪrho��xyƽ����ͶӰ��x�н�
//phiΪrho��z��ļн�
	float fSin_Phi = sin(phi),
		fSin_Theta = sin(theta),
		fCos_Theta = cos(theta);
	*px = rho * fSin_Phi * fCos_Theta;
	*py = rho * fSin_Phi * fSin_Theta;
	*pz = rho * cos(phi);
	return;
}
template<typename _T>void Elementary_Row_Operation_1(_T A[], int m, int n, _T A_2[], int* piRank, _T** ppBasic_Solution, _T** ppSpecial_Solution)
{//�����Elementary_Row_Operation��ʲô�������� :)

	typedef struct Q_Item {
		unsigned short m_iRow_Index;	//��ǰ�ж�Ӧ������Ԫ����������
		unsigned short m_iCol_Index;	//����Ԫ��Ӧ����������x����
	}Q_Item;

	int y, x, x_1, i, iRank = 0, iPos, iMax;
	Q_Item* Q = NULL, iTemp;
	_T* pBasic_Solution = NULL, * pSpecial_Solution = NULL;
	short* pMap_Row_2_x_Index = NULL, * pMap_x_2_Basic_Solution_Index = NULL;
	int j, iRank_Basic_Solution;

	_T fValue, fMax, * A_1 = NULL;
	union {
		_T* pfMax_Row;
		_T* pfBottom_Row;
		_T* pfCur_Row;
	};
	if (piRank)
		*piRank = 0;
	if (m > 65535)
	{
		printf("Too large row count:%d\n", m);
		goto END;
	}
	A_1 = (_T*)pMalloc(&oMem_Mgr, m * n * sizeof(_T));
	Q = (Q_Item*)pMalloc(&oMem_Mgr, m * sizeof(Q_Item));
	if (A_1)
		memcpy(A_1, A, m * n * sizeof(_T));
	iPos = 0;
	for (y = 0; y < m; y++)
		Q[y] = { (unsigned short)y };	//ÿ����Ԫ���ڵ���

	//Disp(A_1, m, n, "\n");
	for (x_1 = 0, y = 0; y < m; y++)
	{//�������y��x�����ƽ����������
		while (1)
		{
			iMax = y;
			fMax = A_1[Q[iMax].m_iRow_Index * n + x_1];
			for (i = y + 1; i < m; i++)
			{
				if (abs(A_1[iPos = Q[i].m_iRow_Index * n + x_1]) > abs(fMax))
				{
					fMax = A_1[iPos];
					iMax = i;
				}
			}
			if (abs(fMax) <= ZERO_APPROCIATE && x_1 < n - 1)
				x_1++;
			else
				break;
		}

		if (abs(fMax) < ZERO_APPROCIATE)
		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
			//Disp(A_1, m, n,"\n");
			break;
		}

		//�����ԪSWAP��Q�ĵ�ǰλ����
		iTemp = Q[y];
		Q[y] = Q[iMax];
		Q[iMax] = iTemp;
		Q[y].m_iCol_Index = x_1;
		iRank++;

		//Disp(A_1, m, n, "\n");
		pfMax_Row = &A_1[Q[y].m_iRow_Index * n];
		pfMax_Row[x_1] = 1.f;
		for (x = x_1 + 1; x < n; x++)
			pfMax_Row[x] /= fMax;
		//Disp(A_1, m, n, "\n");

		//�Ժ��������д���
		for (i = y + 1; i < m; i++)
			//for (i = 0; i < m; i++)
		{//i��ʾ��i��
			iPos = Q[i].m_iRow_Index * n;
			if (((fValue = A_1[iPos + x_1]) != 0) && i != y)
			{//���ڶ�ӦԪ��Ϊ0�����������
				for (x = x_1; x < n; x++)
					A_1[iPos + x] -= fValue * pfMax_Row[x];
				A_1[iPos + x_1] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
			}
			//Disp(Ai, iOrder, iOrder + 1, "\n");
		}
		//Disp(A_1, m, n, "\n");
		x_1++;
	}

	//Disp(A_1, m, n,"�����б任");
	int y1;	//�Ѿ���֪�������
	//Ȼ��˳������һ�����ϱ任����θ���˹�����Է��̲�һ��
	//����Ժ�A_1����������
	for (y = iRank - 1; y > 0; y--)
	{//�߼��ϴ�����һ�����ϣ�ʵ������Qָ·
		pfBottom_Row = &A_1[Q[y].m_iRow_Index * n];
		x_1 = Q[y].m_iCol_Index;	//ǰ���Ѿ��õ���������Ԫλ��

		for (y1 = y - 1; y1 >= 0; y1--)
		{
			//iPos = Q[y_1] * iRow_Size;
			iPos = Q[y1].m_iRow_Index * n;
			x = x_1;	//�����е�xλ��
			fValue = A_1[iPos + x];
			A_1[iPos + x] = 0;
			for (x++; x < n; x++)
				A_1[iPos + x] -= fValue * pfBottom_Row[x];
			//Disp(A_1, m, n, "\n");
		}
	}
	//Disp(A_1, m, n, "A_1");
	if (piRank)
		*piRank = iRank;
	if (A_2)
		memcpy(A_2, A_1, m * n * sizeof(_T));
	iRank_Basic_Solution = n - 1 - iRank;	//������ϵ����
	if (!ppBasic_Solution || !ppSpecial_Solution)
		goto END;

	//���һ��������η����������ϵ��������������������Ϊn-1ά����������n-1- Rank��
	pBasic_Solution = (_T*)pMalloc(&oMem_Mgr, iRank_Basic_Solution * n * sizeof(_T));
	//����Ϊ�����Ľ�x��Ӧ�ĸ���������һ��Map
	pMap_Row_2_x_Index = (short*)pMalloc(&oMem_Mgr, (n - 1) * sizeof(short));
	pMap_x_2_Basic_Solution_Index = (short*)pMalloc(&oMem_Mgr, (n - 1) * sizeof(short));

	memset(pBasic_Solution, 0, iRank_Basic_Solution * n * sizeof(_T));
	memset(pMap_Row_2_x_Index, 0, (n - 1) * sizeof(short));
	memset(pMap_x_2_Basic_Solution_Index, 0, (n - 1) * sizeof(short));

	//����ǰ������Ԫ��xλ��Ϊ-1
	for (i = 0; i < iRank; i++)
		pMap_Row_2_x_Index[Q[i].m_iCol_Index] = -1;

	//ʣ�µ�ֵΪ0�ľ��ǻ�����ϵ��������Ӧ��xλ��
	for (j = 0, i = 0; i < n - 1; i++)
	{
		if (pMap_Row_2_x_Index[i] == 0)
		{
			pMap_x_2_Basic_Solution_Index[i] = j;
			pBasic_Solution[i * iRank_Basic_Solution + j] = 1;
			j++;
		}
	}

	//��󣬹�����η��̻�������
	for (y = 0; y < n; y++)
	{
		iPos = Q[y].m_iRow_Index * n;
		pfCur_Row = &A_1[iPos];
		x_1 = Q[y].m_iCol_Index + 1;
		for (; x_1 < n - 1; x_1++)
		{
			if (Abs(pfCur_Row[x_1]) > ZERO_APPROCIATE)
			{
				//�������кţ������к������к���أ���Ȼȡ�����кž����к�
				pBasic_Solution[Q[y].m_iCol_Index * iRank_Basic_Solution + pMap_x_2_Basic_Solution_Index[x_1]] = -A_1[iPos + x_1];
			}
		}
	}

	if (pSpecial_Solution)
	{
		memset(pSpecial_Solution, 0, (n - 1) * sizeof(_T));
		for (i = 0; i < iRank; i++)
			pSpecial_Solution[Q[i].m_iCol_Index] = A_1[Q[i].m_iRow_Index * n + (n - 1)];
	}
END:
	if (A_1)
		Free(&oMem_Mgr, A_1);
	if (pMap_Row_2_x_Index)
		Free(&oMem_Mgr, pMap_Row_2_x_Index);
	if (pMap_x_2_Basic_Solution_Index)
		Free(&oMem_Mgr, pMap_x_2_Basic_Solution_Index);
	if (ppBasic_Solution)
		*ppBasic_Solution = pBasic_Solution;
	if (pBasic_Solution && !ppBasic_Solution)
		Free(&oMem_Mgr, pBasic_Solution);
	if(ppSpecial_Solution)
		*ppSpecial_Solution = pSpecial_Solution;
	if (pSpecial_Solution && !ppSpecial_Solution)
			Free(&oMem_Mgr, pSpecial_Solution);
	if (Q)
		Free(&oMem_Mgr, Q);
	return;
}

template<typename _T>void Elementary_Row_Operation(_T A[], int m, int n, _T A_1[], int* piRank, _T** ppBasic_Solution, _T** ppSpecial_Solution)
{//��A�������б任����Ϊ����Σ�Ҫ������Ԫ��
//����A����
	typedef struct Q_Item {
		unsigned char m_iRow_Index;	//��ǰ�ж�Ӧ������Ԫ����������
		unsigned char m_iCol_Index;	//			����Ԫ��Ӧ����������x����
	}Q_Item;

	int y, x, x_1, i, iRank = 0, iPos, iMax;
	Q_Item Q[256], iTemp;
	_T fValue, fMax;
	union {
		_T* pfMax_Row;
		_T* pfBottom_Row;
		_T* pfCur_Row;
	};

	if (m > 256)
	{
		printf("Too large row count:%d\n", m);
		return;
	}
	if (A_1)
		memcpy(A_1, A, m * n * sizeof(_T));
	iPos = 0;
	for (y = 0; y < m; y++)
		Q[y] = { (unsigned char)y };	//ÿ����Ԫ���ڵ���

	//Disp(A_1, m, n, "\n");
	for (x_1 = 0, y = 0; y < m; y++)
	{//�������y��x�����ƽ����������
		while (1)
		{
			iMax = y;
			fMax = A_1[Q[iMax].m_iRow_Index * n + x_1];
			for (i = y + 1; i < m; i++)
			{
				if (abs(A_1[iPos = Q[i].m_iRow_Index * n + x_1]) > abs(fMax))
				{
					fMax = A_1[iPos];
					iMax = i;
				}
			}
			if (abs(fMax) <= ZERO_APPROCIATE && x_1 < n - 1)
				x_1++;
			else
				break;
		}

		if (abs(fMax) < ZERO_APPROCIATE)
		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
			//Disp(A_1, m, n,"\n");
			break;
		}

		//�����ԪSWAP��Q�ĵ�ǰλ����
		iTemp = Q[y];
		Q[y] = Q[iMax];
		Q[iMax] = iTemp;
		Q[y].m_iCol_Index = x_1;
		iRank++;

		//Disp(A_1, m, n, "\n");
		pfMax_Row = &A_1[Q[y].m_iRow_Index * n];
		pfMax_Row[x_1] = 1.f;
		for (x = x_1 + 1; x < n; x++)
			pfMax_Row[x] /= fMax;
		//Disp(A_1, m, n, "\n");

		//�Ժ��������д���
		for (i = y + 1; i < m; i++)
			//for (i = 0; i < m; i++)
		{//i��ʾ��i��
			iPos = Q[i].m_iRow_Index * n;
			if (((fValue = A_1[iPos + x_1]) != 0) && i != y)
			{//���ڶ�ӦԪ��Ϊ0�����������
				for (x = x_1; x < n; x++)
					A_1[iPos + x] -= fValue * pfMax_Row[x];
				A_1[iPos + x_1] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
			}
			//Disp(Ai, iOrder, iOrder + 1, "\n");
		}
		//Disp(A_1, m, n, "\n");
		x_1++;
	}

	//Disp(A_1, m, n, "A1");
	//Disp(A_1, m, n,"�����б任");
	int y1;	//�Ѿ���֪�������
	//Ȼ��˳������һ�����ϱ任����θ���˹�����Է��̲�һ��
	//����Ժ�A_1����������
	for (y = iRank - 1; y > 0; y--)
	{//�߼��ϴ�����һ�����ϣ�ʵ������Qָ·
		pfBottom_Row = &A_1[Q[y].m_iRow_Index * n];
		x_1 = Q[y].m_iCol_Index;	//ǰ���Ѿ��õ���������Ԫλ��

		for (y1 = y - 1; y1 >= 0; y1--)
		{
			//iPos = Q[y_1] * iRow_Size;
			iPos = Q[y1].m_iRow_Index * n;
			x = x_1;	//�����е�xλ��
			fValue = A_1[iPos + x];
			A_1[iPos + x] = 0;
			for (x++; x < n; x++)
				A_1[iPos + x] -= fValue * pfBottom_Row[x];
			//Disp(A_1, m, n, "\n");
		}
	}

	/*for (y = 0; y < m; y++)
		printf("y:%d %d\n", y, Q[y].m_iRow_Index);
	_T A2[18*18], Inv[18*18];
	for (y = 0; y < m; y++)
		for (x = 0; x < m; x++)
		{
			A2[y * m + x] = A[y * n + x];
			Inv[y * m + x] = A_1[Q[y].m_iRow_Index * n + m + x];
		}
	Matrix_Multiply(A2, m, m, Inv, m, A2);
	Disp(A2, m, m, "A2");	*/

	int j, iRank_Basic_Solution = n - 1 - iRank;	//������ϵ����
	//���һ��������η����������ϵ��������������������Ϊn-1ά����������n-1- Rank��
	_T* pBasic_Solution = (_T*)malloc(iRank_Basic_Solution * n * sizeof(_T));
	//����Ϊ�����Ľ�x��Ӧ�ĸ���������һ��Map
	short* pMap_Row_2_x_Index = (short*)malloc(n * sizeof(short));
	short* pMap_x_2_Basic_Solution_Index = (short*)malloc(n * sizeof(short));

	memset(pBasic_Solution, 0, iRank_Basic_Solution * n * sizeof(_T));
	memset(pMap_Row_2_x_Index, 0, (n - 1) * sizeof(short));
	memset(pMap_x_2_Basic_Solution_Index, 0, (n - 1) * sizeof(short));

	//����ǰ������Ԫ��xλ��Ϊ-1
	for (i = 0; i < iRank; i++)
		pMap_Row_2_x_Index[Q[i].m_iCol_Index] = -1;

	//ʣ�µ�ֵΪ0�ľ��ǻ�����ϵ��������Ӧ��xλ��
	for (j = 0, i = 0; i < n - 1; i++)
	{
		if (pMap_Row_2_x_Index[i] == 0)
		{
			pMap_x_2_Basic_Solution_Index[i] = j;
			pBasic_Solution[i * iRank_Basic_Solution + j] = 1;
			j++;
		}
	}

	//��󣬹�����η��̻�������
	for (y = 0; y < n; y++)
	{
		iPos = Q[y].m_iRow_Index * n;
		pfCur_Row = &A_1[iPos];
		x_1 = Q[y].m_iCol_Index + 1;
		for (; x_1 < n - 1; x_1++)
		{
			if (Abs(pfCur_Row[x_1]) > ZERO_APPROCIATE)
			{
				//�������кţ������к������к���أ���Ȼȡ�����кž����к�
				pBasic_Solution[Q[y].m_iCol_Index * iRank_Basic_Solution + pMap_x_2_Basic_Solution_Index[x_1]] = -A_1[iPos + x_1];
			}
		}
	}
	//Disp(pBasic_Solution, n - 1, iRank_Basic_Solution);
	//����һ���ؽ�
	_T* pSpecial_Solution = (_T*)malloc(n * sizeof(_T));
	memset(pSpecial_Solution, 0, (n - 1) * sizeof(_T));
	for (i = 0; i < iRank; i++)
		pSpecial_Solution[Q[i].m_iCol_Index] = A_1[Q[i].m_iRow_Index * n + (n - 1)];

	if (ppSpecial_Solution)
		*ppSpecial_Solution = pSpecial_Solution;
	else
		if(pSpecial_Solution)
			free(pSpecial_Solution);
	//Disp(A_1, m, n);
	if (piRank)
		*piRank = iRank;
	if (ppBasic_Solution)
		*ppBasic_Solution = pBasic_Solution;
	else
		if (pBasic_Solution)
			free(pBasic_Solution);

	/*if (pQ)
		for (y = 0; y < m; y++)
			pQ[y] = Q[y].m_iRow_Index;*/
	free(pMap_Row_2_x_Index);
	free(pMap_x_2_Basic_Solution_Index);
	return;
}

//void Elementary_Row_Operation(float A[], int m, int n, float A_1[], int* piRank, float** ppBasic_Solution, float** ppSpecial_Solution)
//{//��A�������б任����Ϊ����Σ�Ҫ������Ԫ��
////����A����
//	typedef struct Q_Item {
//		unsigned char m_iRow_Index;	//��ǰ�ж�Ӧ������Ԫ����������
//		unsigned char m_iCol_Index;	//			����Ԫ��Ӧ����������x����
//	}Q_Item;
//
//	int y, x, x_1, i, iRank = 0, iPos, iMax;
//	Q_Item Q[256], iTemp;
//	float fValue, fMax;
//	union {
//		float* pfMax_Row;
//		float* pfBottom_Row;
//		float* pfCur_Row;
//	};
//
//	if (m > 256)
//	{
//		printf("Too large row count:%d\n", m);
//		return;
//	}
//	if (A_1)
//		memcpy(A_1, A, m * n * sizeof(float));
//	iPos = 0;
//	for (y = 0; y < m; y++)
//		Q[y] = { (unsigned char)y };	//ÿ����Ԫ���ڵ���
//
//	//Disp(A_1, m, n, "\n");
//	for (x_1 = 0, y = 0; y < m; y++)
//	{//�������y��x�����ƽ����������
//		while (1)
//		{
//			iMax = y;
//			fMax = A_1[Q[iMax].m_iRow_Index * n + x_1];
//			for (i = y + 1; i < m; i++)
//			{
//				if (abs(A_1[iPos = Q[i].m_iRow_Index * n + x_1]) > abs(fMax))
//				{
//					fMax = A_1[iPos];
//					iMax = i;
//				}
//			}
//			if (abs(fMax) <= ZERO_APPROCIATE && x_1 < n - 1)
//				x_1++;
//			else
//				break;
//		}
//
//		if (abs(fMax) < ZERO_APPROCIATE)
//		{//����ԪΪ0����Ȼ�����ȣ��÷���û��Ψһ��
//			//Disp(A_1, m, n,"\n");
//			break;
//		}
//
//		//�����ԪSWAP��Q�ĵ�ǰλ����
//		iTemp = Q[y];
//		Q[y] = Q[iMax];
//		Q[iMax] = iTemp;
//		Q[y].m_iCol_Index = x_1;
//		iRank++;
//
//		//Disp(A_1, m, n, "\n");
//		pfMax_Row = &A_1[Q[y].m_iRow_Index * n];
//		pfMax_Row[x_1] = 1.f;
//		for (x = x_1 + 1; x < n; x++)
//			pfMax_Row[x] /= fMax;
//		//Disp(A_1, m, n, "\n");
//
//		//�Ժ��������д���
//		for (i = y + 1; i < m; i++)
//			//for (i = 0; i < m; i++)
//		{//i��ʾ��i��
//			iPos = Q[i].m_iRow_Index * n;
//			if (((fValue = A_1[iPos + x_1]) != 0) && i != y)
//			{//���ڶ�ӦԪ��Ϊ0�����������
//				for (x = x_1; x < n; x++)
//					A_1[iPos + x] -= fValue * pfMax_Row[x];
//				A_1[iPos + x_1] = 0;	//�˴�Ҳ���Ǳ���ģ�����ֻ�Ǻÿ�
//			}
//			//Disp(Ai, iOrder, iOrder + 1, "\n");
//		}
//		//Disp(A_1, m, n, "\n");
//		x_1++;
//	}
//
//	//Disp(A_1, m, n,"�����б任");
//	int y1;	//�Ѿ���֪�������
//	//Ȼ��˳������һ�����ϱ任����θ���˹�����Է��̲�һ��
//	//����Ժ�A_1����������
//	for (y = iRank - 1; y > 0; y--)
//	{//�߼��ϴ�����һ�����ϣ�ʵ������Qָ·
//		pfBottom_Row = &A_1[Q[y].m_iRow_Index * n];
//		x_1 = Q[y].m_iCol_Index;	//ǰ���Ѿ��õ���������Ԫλ��
//
//		for (y1 = y - 1; y1 >= 0; y1--)
//		{
//			//iPos = Q[y_1] * iRow_Size;
//			iPos = Q[y1].m_iRow_Index * n;
//			x = x_1;	//�����е�xλ��
//			fValue = A_1[iPos + x];
//			A_1[iPos + x] = 0;
//			for (x++; x < n; x++)
//				A_1[iPos + x] -= fValue * pfBottom_Row[x];
//			//Disp(A_1, m, n, "\n");
//		}
//	}
//
//	/*if (iRank == n-1)
//	{
//		*ppBasic_Solution = NULL;
//		return;
//	}*/
//
//
//	int j, iRank_Basic_Solution = n - 1 - iRank;	//������ϵ����
//	//���һ��������η����������ϵ��������������������Ϊn-1ά����������n-1- Rank��
//	float* pBasic_Solution = (float*)malloc(iRank_Basic_Solution * n * sizeof(float));
//	//����Ϊ�����Ľ�x��Ӧ�ĸ���������һ��Map
//	short* pMap_Row_2_x_Index = (short*)malloc((n - 1) * sizeof(short));
//	short* pMap_x_2_Basic_Solution_Index = (short*)malloc((n - 1) * sizeof(short));
//
//	memset(pBasic_Solution, 0, iRank_Basic_Solution * n * sizeof(float));
//	memset(pMap_Row_2_x_Index, 0, (n - 1) * sizeof(short));
//	memset(pMap_x_2_Basic_Solution_Index, 0, (n - 1) * sizeof(short));
//
//	//����ǰ������Ԫ��xλ��Ϊ-1
//	for (i = 0; i < iRank; i++)
//		pMap_Row_2_x_Index[Q[i].m_iCol_Index] = -1;
//
//	//ʣ�µ�ֵΪ0�ľ��ǻ�����ϵ��������Ӧ��xλ��
//	for (j = 0, i = 0; i < n - 1; i++)
//	{
//		if (pMap_Row_2_x_Index[i] == 0)
//		{
//			pMap_x_2_Basic_Solution_Index[i] = j;
//			pBasic_Solution[i * iRank_Basic_Solution + j] = 1;
//			j++;
//		}
//	}
//
//	//��󣬹�����η��̻�������
//	for (y = 0; y < n; y++)
//	{
//		iPos = Q[y].m_iRow_Index * n;
//		pfCur_Row = &A_1[iPos];
//		x_1 = Q[y].m_iCol_Index + 1;
//		for (; x_1 < n - 1; x_1++)
//		{
//			if (Abs(pfCur_Row[x_1]) > ZERO_APPROCIATE)
//			{
//				//�������кţ������к������к���أ���Ȼȡ�����кž����к�
//				pBasic_Solution[Q[y].m_iCol_Index * iRank_Basic_Solution + pMap_x_2_Basic_Solution_Index[x_1]] = -A_1[iPos + x_1];
//			}
//		}
//	}
//	//Disp(pBasic_Solution, n - 1, iRank_Basic_Solution);
//	//����һ���ؽ�
//	float* pSpecial_Solution = (float*)malloc((n - 1) * sizeof(float));
//	memset(pSpecial_Solution, 0, (n - 1) * sizeof(float));
//	for (i = 0; i < iRank; i++)
//		pSpecial_Solution[Q[i].m_iCol_Index] = A_1[Q[i].m_iRow_Index * n + (n - 1)];
//
//	if (ppSpecial_Solution)
//		*ppSpecial_Solution = pSpecial_Solution;
//	else
//		free(pSpecial_Solution);
//	//Disp(A_1, m, n);
//	if (piRank)
//		*piRank = iRank;
//	if (ppBasic_Solution)
//		*ppBasic_Solution = pBasic_Solution;
//	else
//		if (pBasic_Solution)
//			free(pBasic_Solution);
//
//	/*if (pQ)
//		for (y = 0; y < m; y++)
//			pQ[y] = Q[y].m_iRow_Index;*/
//	free(pMap_Row_2_x_Index);
//	free(pMap_x_2_Basic_Solution_Index);
//	return;
//}

//void Solve_Linear_Solution_Construction(float* A, int m, int n, float B[], int* pbSuccess, float* pBasic_Solution, int* piBasic_Solution_Count, float* pSpecial_Solution)
//{//�����Է��̣������������������ʽ
//	float* Ai = (float*)malloc(m * (n + 1) * sizeof(float));
//	float* pBasic_Solution_1, * pSpecial_Solution_1;
//	int y, x, iRow_Size = n + 1;
//	int iRank;	//ϵ���������
//
//	for (y = 0; y < m; y++)
//	{//AiΪ�������
//		for (x = 0; x < n; x++)
//			Ai[y * iRow_Size + x] = A[y * n + x];
//		Ai[y * iRow_Size + n] = B[y];
//	}
//	//�˴����г����б任�������ø����׵�����Ԫ
//	Elementary_Row_Operation(Ai, m, n + 1, Ai, &iRank, &pBasic_Solution_1, &pSpecial_Solution_1);
//	for (y = 0; y < m; y++)
//	{
//		int bIs_Zero = 1;
//		for (x = 0; x < n; x++)
//		{
//			if (abs(Ai[y * (n + 1) + x]) > ZERO_APPROCIATE)
//			{
//				bIs_Zero = 0;
//				break;
//			}
//		}
//		if (bIs_Zero && Ai[y * (n + 1) + n] > ZERO_APPROCIATE)
//		{
//			printf("�÷����޽⣬���������б任�Ժ󣬵�%d�еĳ�����Ϊ��%f\n", y, Ai[y * (n + 1) + n]);
//			Disp(Ai, m, n + 1);
//			*pbSuccess = 0;
//			return;
//		}
//	}
//	int bIs_Homo = 1;
//	for (y = 0; y < m; y++)
//	{
//		if (B[y])
//		{//���������b��Ϊ0����Ϊ����η���
//			bIs_Homo = 0;
//			break;
//		}
//	}
//	if (bIs_Homo)
//	{//��η�����Ľ⣬���� X= c0*v0 + c1*v1 + ... + c(n-r) * v(n-r)
//		if (iRank == n)
//		{
//			printf("ϵ���������ȣ�ֻ�����\n");
//			if (piBasic_Solution_Count)
//				*piBasic_Solution_Count = 0;
//			goto END;
//		}
//	}
//	if (piBasic_Solution_Count)
//		*piBasic_Solution_Count = n - iRank;
//	if (pBasic_Solution)
//	{
//		Matrix_Transpose(pBasic_Solution_1, n, n - iRank, pBasic_Solution_1);
//		Schmidt_Orthogon(pBasic_Solution_1, n - iRank, n, pBasic_Solution_1);
//		memcpy(pBasic_Solution, pBasic_Solution_1, (n - iRank) * n * sizeof(float));
//	}
//	if (!bIs_Homo && pSpecial_Solution)
//		memcpy(pSpecial_Solution, pSpecial_Solution_1, n * sizeof(float));
//
//END:
//	if (pSpecial_Solution_1)
//		free(pSpecial_Solution_1);
//	if (pBasic_Solution_1)
//		free(pBasic_Solution_1);
//	if (Ai)
//		free(Ai);
//	*pbSuccess = 1;
//}
template<typename _T> void Solve_Linear_Solution_Construction_1(_T* A, int m, int n, _T B[], int* pbSuccess, _T* pBasic_Solution, int* piBasic_Solution_Count, _T* pSpecial_Solution)
{//���һ����Ľṹ
	//�˴��������������ڴ�. m��n�е�������󣬼�һ�У�Ϊ m*(n+1)��Ԫ��
	//_T* Ai = (_T*)malloc(m * (n + 1) * sizeof(_T));
	_T* Ai = (_T*)pMalloc(m * (n + 1) * sizeof(_T));
	_T* pBasic_Solution_1=NULL, * pSpecial_Solution_1=NULL;
	int y, x, iRow_Size = n + 1;
	int iRank=0;	//ϵ���������
		
	for (y = 0; y < m; y++)
	{//AiΪ�������
		for (x = 0; x < n; x++)
			Ai[y * iRow_Size + x] = A[y * n + x];
		Ai[y * iRow_Size + n] = B?B[y]:0;
	}
	
	//�˴����г����б任�������ø����׵�����Ԫ
	Elementary_Row_Operation_1(Ai, m, n + 1, Ai, &iRank, &pBasic_Solution_1, &pSpecial_Solution_1);

	for (y = 0; y < m; y++)
	{
		int bIs_Zero = 1;
		for (x = 0; x < n; x++)
		{
			if (abs(Ai[y * (n + 1) + x]) > ZERO_APPROCIATE)
			{
				bIs_Zero = 0;
				break;
			}
		}
		if (bIs_Zero && Ai[y * (n + 1) + n] > ZERO_APPROCIATE)
		{
			printf("�÷����޽⣬���������б任�Ժ󣬵�%d�еĳ�����Ϊ��%f\n", y, Ai[y * (n + 1) + n]);
			Disp(Ai, m, n + 1);
			*pbSuccess = 0;
			return;
		}
	}
	int bIs_Homo = 1;
	int iBasic_Solution_Count;
	if (B)
	{
		for (y = 0; y < m; y++)
		{
			if (B[y])
			{//���������b��Ϊ0����Ϊ����η���
				bIs_Homo = 0;
				break;
			}
		}
	}

	if (bIs_Homo)
	{//��η�����Ľ⣬���� X= c0*v0 + c1*v1 + ... + c(n-r) * v(n-r)
		if (iRank == n)
		{
			printf("ϵ���������ȣ�ֻ�����\n");
			//if (piBasic_Solution_Count)
				//*piBasic_Solution_Count = 0;
			iBasic_Solution_Count = 0;
		}
		else
			iBasic_Solution_Count = n - iRank;
	}
	else
		iBasic_Solution_Count = n - iRank;

	if (piBasic_Solution_Count)
		*piBasic_Solution_Count = iBasic_Solution_Count;

	if (pBasic_Solution && pBasic_Solution_1)
	{
		Matrix_Transpose(pBasic_Solution_1, n, n - iRank, pBasic_Solution_1);
		//�˴���һ�²�����һ��Ҫ�������������������᲻��ı��˻��ļнǣ�
		Schmidt_Orthogon(pBasic_Solution_1, n - iRank, n, pBasic_Solution_1);
		if (iBasic_Solution_Count)
			memcpy(pBasic_Solution, pBasic_Solution_1, (n - iRank) * n * sizeof(_T));
		else
			memset(pBasic_Solution, 0, m * sizeof(_T));
	}
	if (!bIs_Homo && pSpecial_Solution)
		memcpy(pSpecial_Solution, pSpecial_Solution_1, n * sizeof(_T));

	if (pSpecial_Solution_1)
		Free(pSpecial_Solution_1);
	if (pBasic_Solution_1)
		Free(pBasic_Solution_1);
	if (Ai)
		Free(Ai);
	*pbSuccess = 1;
}
template<typename _T> void Solve_Linear_Solution_Construction(_T* A, int m, int n, _T B[], int* pbSuccess, _T* pBasic_Solution, int* piBasic_Solution_Count, _T* pSpecial_Solution)
{//���׼������������malloc
	//�����Է��̣������������������ʽ
	_T* Ai = (_T*)malloc(m * (n + 1) * sizeof(_T));
	_T* pBasic_Solution_1, * pSpecial_Solution_1;
	int y, x, iRow_Size = n + 1;
	int iRank;	//ϵ���������

	for (y = 0; y < m; y++)
	{//AiΪ�������
		for (x = 0; x < n; x++)
			Ai[y * iRow_Size + x] = A[y * n + x];
		Ai[y * iRow_Size + n] = B[y];
	}

	//�˴����г����б任�������ø����׵�����Ԫ
	Elementary_Row_Operation(Ai, m, n + 1, Ai, &iRank, &pBasic_Solution_1, &pSpecial_Solution_1);
	for (y = 0; y < m; y++)
	{
		int bIs_Zero = 1;
		for (x = 0; x < n; x++)
		{
			if (abs(Ai[y * (n + 1) + x]) > ZERO_APPROCIATE)
			{
				bIs_Zero = 0;
				break;
			}
		}
		if (bIs_Zero && Ai[y * (n + 1) + n] > ZERO_APPROCIATE)
		{
			printf("�÷����޽⣬���������б任�Ժ󣬵�%d�еĳ�����Ϊ��%f\n", y, Ai[y * (n + 1) + n]);
			Disp(Ai, m, n + 1);
			*pbSuccess = 0;
			return;
		}
	}
	int bIs_Homo = 1;
	int iBasic_Solution_Count;
	for (y = 0; y < m; y++)
	{
		if (B[y])
		{//���������b��Ϊ0����Ϊ����η���
			bIs_Homo = 0;
			break;
		}
	}
	if (bIs_Homo)
	{//��η�����Ľ⣬���� X= c0*v0 + c1*v1 + ... + c(n-r) * v(n-r)
		if (iRank == n)
		{
			printf("ϵ���������ȣ�ֻ�����\n");
			//if (piBasic_Solution_Count)
				//*piBasic_Solution_Count = 0;
			iBasic_Solution_Count = 0;
		}
		else
			iBasic_Solution_Count = n - iRank;
	}
	else
		iBasic_Solution_Count = n - iRank;

	if (piBasic_Solution_Count)
		*piBasic_Solution_Count = iBasic_Solution_Count;

	if (pBasic_Solution)
	{
		Matrix_Transpose(pBasic_Solution_1, n, n - iRank, pBasic_Solution_1);
		Schmidt_Orthogon(pBasic_Solution_1, n - iRank, n, pBasic_Solution_1);
		if (iBasic_Solution_Count)
			memcpy(pBasic_Solution, pBasic_Solution_1, (n - iRank) * n * sizeof(_T));
		else
			memset(pBasic_Solution, 0, m * sizeof(_T));
	}
	if (!bIs_Homo && pSpecial_Solution)
		memcpy(pSpecial_Solution, pSpecial_Solution_1, n * sizeof(_T));

	if (pSpecial_Solution_1)
		free(pSpecial_Solution_1);
	if (pBasic_Solution_1)
		free(pBasic_Solution_1);
	if (Ai)
		free(Ai);
	*pbSuccess = 1;
}

void Get_Linear_Solution_Construction(float A[], const int m, int n, float B[])
{//�����Է������Ľṹ
	float* Ai = (float*)malloc(m * (n + 1) * sizeof(float));
	int y, x, iRow_Size = n + 1;
	int iRank;	//ϵ���������
	for (y = 0; y < m; y++)
	{//AiΪ�������
		for (x = 0; x < n; x++)
			Ai[y * iRow_Size + x] = A[y * n + x];
		Ai[y * iRow_Size + n] = B[y];
	}
	Disp(Ai, m, n + 1, "Aumented Matrix");
	//�˴����г����б任�������ø����׵�����Ԫ
	float* pBasic_Solution = NULL,
		* pSpecial_Solution = NULL;

	Elementary_Row_Operation(Ai, m, n + 1, Ai, &iRank, &pBasic_Solution, &pSpecial_Solution);
	for (y = 0; y < m; y++)
	{
		int bIs_Zero = 1;
		for (x = 0; x < n; x++)
		{
			if (abs(Ai[y * (n + 1) + x]) > ZERO_APPROCIATE)
			{
				bIs_Zero = 0;
				break;
			}
		}
		if (bIs_Zero && Ai[y * (n + 1) + n] > ZERO_APPROCIATE)
		{
			printf("�÷����޽⣬���������б任�Ժ󣬵�%d�еĳ�����Ϊ��%f\n", y, Ai[y * (n + 1) + n]);
			Disp(Ai, m, n + 1);
			return;
		}
	}
	int bIs_Homo = 1;
	for (y = 0; y < m; y++)
	{
		if (B[y])
		{//���������b��Ϊ0����Ϊ����η���
			bIs_Homo = 0;
			break;
		}
	}
	//������������Ӳ�磬����Ԫ�������Ժ󲢲�������Ρ���զŪ
	Disp(Ai, m, n + 1, "����ξ�������");
	if (bIs_Homo)
	{//��η�����Ľ⣬���� X= c0*v0 + c1*v1 + ... + c(n-r) * v(n-r)
		if (iRank == n)
		{
			printf("ϵ���������ȣ�ֻ�����\n");
			goto END;
		}

		////����
		float* pB_1 = (float*)malloc(m * (n - iRank) * sizeof(float));
		Disp(pBasic_Solution, n, n - iRank, "������ϵ");
		Matrix_Multiply(A, m, n, pBasic_Solution, n - iRank, pB_1);

		float* pBasic_Solution_1 = (float*)malloc((n - iRank) * n * sizeof(float));
		Matrix_Transpose(pBasic_Solution, n, n - iRank, pBasic_Solution_1);
		Schmidt_Orthogon(pBasic_Solution_1, n - iRank, n, pBasic_Solution_1);
		Disp(pBasic_Solution_1, n - iRank, n, "������");
		free(pBasic_Solution_1);

		Disp(pB_1, m, n - iRank, "Axb");
		free(pB_1);

		////��һ��
		//float* pBasic_Solution_1 = (float*)malloc((n - iRank) * n * sizeof(float));
		//Matrix_Transpose(pBasic_Solution, n, n-iRank, pBasic_Solution_1);
		//
		//for (y = 0; y < n - iRank; y++)
		//	Normalize(&pBasic_Solution_1[y*n], n, &pBasic_Solution_1[y*n]);

		//Disp(pBasic_Solution_1, n - iRank, n);
		//free(pBasic_Solution_1);

		//float Matlab_Value[] = {0.3876f, -0.1577f ,0.3444f ,-0.1989f , -0.8165f};
		//printf("Mod:%f\n", fGet_Mod(Matlab_Value, n));
	}
	else
	{
		Disp(pSpecial_Solution, n, 1, "�ؽ�");
		if (iRank == n)
		{
			printf("ϵ���������ȣ�ֻ��Ψһ��\n");
			goto END;
		}
		Disp(pBasic_Solution, n, n - iRank, "������ϵ");
	}

	//����Ӹ������������
END:
	if (pSpecial_Solution)
		free(pSpecial_Solution);
	if (pBasic_Solution)
		free(pBasic_Solution);
	if (Ai)
		free(Ai);
	return;
}
int bIs_Linearly_Dependent(float A[], int m, int n)
{//��A��m�� nά�������Ƿ��������
	if (m > n)
		return 1;	//��������ά������Ȼ�������

	//float* B = (float*)malloc(m * sizeof(float));
	float* A1 = (float*)malloc(m * n * sizeof(float));

	//memset(B, 0, m * sizeof(float));
	int iRank;
	Elementary_Row_Operation(A, m, n, A1, &iRank);

	if (iRank == m)
		return 0;
	else
		return 1;
}
void Perspective(float Pos_Source[3], float h[3], float Pos_Dest[3])
{//h�����������ϵ���㣬��Ӧ͸��ԭ���е�h����
	//���ڿռ��е�һ��(x,y,z), ���ձ�����͸�ӱ任��Ϊ (x',y',z'). x'/x0= h/(h-z0) => x'= x * h/(h-z0)
	//͸��Ҫ�㣬1���ӵ㼴���
	//			2����λ��z��0��ʱ������ԭλ��
	//			3, ��zΪ��ʱ����ԭͼ�󣬵�zΪ��ʱ����ԭͼС
	//			4, ������û�п���������꣬�ʴ���֮��Щ����
	//			5�����ڱ任֮��
	//			6, ��zλ��Ϊ0ʱ�������Ѿ��任�����޴󣬸�ֵ������
	//			7����zλ�ô����ӵ�hλ��ʱ��ͬ�������壬��Ϊ����󿴲���
	if ((h[2] > 0 && Pos_Source[2] >= h[2]) || (h[2] < 0 && Pos_Source[2] <= h[2]))
	{
		printf("Invalid Pos: z>=h\n");
		return;
	}

	float H =/*h[0]/(h[0]-Pos_Source[0]) +*/ h[2] / (h[2] - Pos_Source[2]);
	Pos_Dest[0] = Pos_Source[0] * H;
	Pos_Dest[1] = Pos_Source[1] * H;
	Pos_Dest[2] = 0;
	return;
}
void Perspective_Camera(float Pos_Source[3], float h[3], float Pos_Dest[3])
{	//Pos_Source: ���λ�ã�fCamera_z������Լ���z���ϵ�λ�ã�Ϊ��
	//float z = h[2] + Pos_Source[2];
	//float H = h[2] / (h[2] - z);
	float /*p = 1.f / (h[2] * 2),*/
		/*q = 1.f / (h[2]),*/
		r = 1.f / (-h[2]);
	/*float p = 0.00215193816,
		r = -0.00186189834;*/
	float H = 1 / (/*p * abs(Pos_Source[0])*/ /*+ q * abs(Pos_Source[1])*/ +r * Pos_Source[2] + 1);
	//float H = 1.f / (p + r + 1);
	//float H =  ( h[2]/(h[2] + abs(Pos_Source[0]))) * (-h[2] / Pos_Source[2]);
	Pos_Dest[0] = Pos_Source[0] * H;
	Pos_Dest[1] = Pos_Source[1] * H;
	Pos_Dest[2] = 0;
}

float fGet_Theta(float v0[], float v1[], float Axis[], int n)
{//����������֮��ļнǣ�nΪά�� �õ������  v0.v1= |v0|*|v1|*cos(theta)
	//�ȵ��
	float fMod_v0 = fGet_Mod(v0, n),
		fMod_v1 = fGet_Mod(v1, n);
	float theta = acos(fDot(v0, v1, n) / (fMod_v0 * fMod_v1));
	//float v2[3];
	//������δ���theta ����������

	//�����ɵ�Ƶķ���
#define eps 0.0001f
	float v0_1[3], v1_1[3], R[3][3], Temp[4];
	float Axis_1[4] = { Axis[0],Axis[1],Axis[2],theta };
	Normalize(v0, 3, v0_1);
	Normalize(v1, 3, v1_1);
	Rotation_Vector_4_2_Matrix(Axis_1, (float*)R);
	Matrix_Multiply((float*)R, 3, 3, v0_1, 1, Temp);
	if (abs(Temp[0] - v1_1[0]) > eps || abs(Temp[1] - v1_1[1]) > eps || abs(Temp[2] - v1_1[2]) > eps)
		theta = -theta;
#undef eps
	//Cross_Product(v0, v1, v2);
	////���ݲ�˶��壬 |v2|=|v0|*|v1|*sin(theta)
	//theta = asin(fGet_Mod(v2, 3) / (fGet_Mod(v0, n) * fGet_Mod(v1, n)));

	return theta;
}

void Diagonalize(float A[], int iOrder, float Diag[], float P[], int* pbSuccess)
{//���һ���Գƾ���ĶԽǻ�������󽫶Խǻ��������Diag[]�У������ƾ���ı任�����P��
//ע�⣬��õĽ��������ϼ���˳�� A = P x Diag x P'
// ��P��һ���������󣬹ʴ�P'�������棬ת�ü���
	float* Diag_1 = (float*)malloc(iOrder * iOrder * sizeof(float)),
		* P_1 = (float*)malloc(iOrder * iOrder * sizeof(float));
	int y, x;
	QR_Decompose(A, iOrder, Diag_1, P_1, pbSuccess);
	if (!(*pbSuccess))
	{
		printf("�þ����޷�������ֵ����������\n");
		goto END;
	}

	//ֻ�����Խ��ߣ�˳�㿴��ԭ�����Ƿ�Գ�
	for (y = 0; y < iOrder; y++)
	{
		for (x = 0; x < iOrder; x++)
		{
			if (x != y)
				Diag_1[y * iOrder + x] = 0;
			else
			{
				if (A[y * iOrder + x] != A[x * iOrder + y])
				{
					printf("ԭ����ǶԳƾ��󣬲��ܶԽǻ�\n");
					*pbSuccess = 0;
					goto END;
				}
			}
		}
	}
	memcpy(Diag, Diag_1, iOrder * iOrder * sizeof(float));
	memcpy(P, P_1, iOrder * iOrder * sizeof(float));

END:
	free(Diag_1);
	free(P_1);
	return;
}
void Disp_Poly(float Coeff[], int iCoeff_Count, char* pCaption = NULL)
{
	int i, n = iCoeff_Count - 1;
	if (pCaption)
		printf("%s\n", pCaption);
	if (n >= 2)
		printf("%fx^%d ", Coeff[n], n);
	else if (n == 1)
		printf("%fx ", Coeff[n]);
	n--;

	for (i = n; i >= 0; i--, n--)
	{
		if (n >= 2)
			printf("+ %fx^%d ", Coeff[i], n);
		else if (n == 1)
			printf("+ %fx ", Coeff[i]);
		else
			printf("+ %f ", Coeff[i]);
	}
	printf("\n");
	return;
}
void Get_Poly_Derivative(float Coeff[], int iCoeff_Count, float Der_Coeff[])
{//�õ�һ������ʽ�ĵ�������ߴ�����[0]��
	int i;
	for (i = 1; i < iCoeff_Count; i++)
		Der_Coeff[i - 1] = i * Coeff[i];
}
void Solve_Poly(float Coeff[], int iCoeff_Count, Complex_f Root[], int* piRoot_Count)
{//�˴����Խ����һԪ����ʽ���ɵķ���
//iCount��ϵ�����������������Ϊn, ����ߴ�Ϊn-1
//Ϊ�˺��±�����ƥ�䣬Coeff��Ĵ�����С������
	Complex_f* pRoot_1 = (Complex_f*)malloc((iCoeff_Count - 1) * sizeof(Complex_f));
	memset(pRoot_1, 0, (iCoeff_Count - 1) * sizeof(Complex_f));
	if (iCoeff_Count == 3)
	{//���η��̣��н�
		float a = Coeff[2], b = Coeff[1], c = Coeff[0];
		float fDelta = b * b - 4 * a * c;
		if (fDelta > ZERO_APPROCIATE)
		{
			if (piRoot_Count)
				*piRoot_Count = 2;
			fDelta = sqrt(fDelta);
			pRoot_1[0] = { (-b + fDelta) / (2.f * a),0.f };
			pRoot_1[1] = { (-b - fDelta) / (2.f * a),0.f };
		}
		else if (fDelta < -ZERO_APPROCIATE)
		{
			if (piRoot_Count)
				*piRoot_Count = 2;
			fDelta = sqrt(-fDelta);
			pRoot_1[0] = { -b / (2.f * a), fDelta / (2.f * a) };
			pRoot_1[1] = { -b / (2.f * a), -fDelta / (2.f * a) };
		}
		else
		{
			if (piRoot_Count)
				*piRoot_Count = 1;
			pRoot_1[0] = pRoot_1[1] = { -b / (2.f * a),0 };
		}
	}
	else if (iCoeff_Count == 2)
	{//���� ax+b=0
		float a = Coeff[1], b = Coeff[0];
		*piRoot_Count = 1;
		pRoot_1[0] = { -b / a,0.f };
	}
	else
	{//�ߴη����ˣ���ʱҪ�õݹ��½������ǳ��鷳
		int n = iCoeff_Count - 1;
		//�ȶ�ԭ����ʽ�󵼣���ø���Ӷ����Ρ������ʽ��=0�����������̣��õ�һ���
		float* pDer_Coeff = (float*)malloc(n * sizeof(float));
		Complex_f* pDer_Root = (Complex_f*)malloc((n - 1) * sizeof(Complex_f));
		int iDer_Root_Count;
		//Disp_Poly(Coeff, iCoeff_Count);
		Get_Poly_Derivative(Coeff, iCoeff_Count, pDer_Coeff);
		//Disp_Poly(pDer_Coeff, iCoeff_Count - 1);
		Solve_Poly(pDer_Coeff, iCoeff_Count - 1, pDer_Root, &iDer_Root_Count);

		return;	//û��
	}

	memcpy(Root, pRoot_1, *piRoot_Count * sizeof(Complex_f));
	for (int i = 0; i < *piRoot_Count; i++)
	{
		if (Root[i].im == 0)
			printf("ʵ��%d:%f\n", i, Root[i].real);
		else
			if (Root[i].im > 0)
				printf("����%d:%f\t+%fi\n", i, Root[i].real, Root[i].im);
			else
				printf("����%d:%f\t%fi\n", i, Root[i].real, Root[i].im);
	}
	free(pRoot_1);
	return;
}

//void Decompose_E(float E[], float R_1[], float R_2[], float t_1[], float t_2[])
//{//��һ��E �ֽ�Ϊ һ������ R[3x3], ����t[3]
//	float R_z_t_1[3 * 3],		//Rz(pi/2)		�������ɳ���
//		R_z_t_2[3 * 3],		//Rz(-pi/2)
//		U[3 * 3], Vt[3 * 3], Ut[3 * 3],
//		S[3], Temp[3 * 3];
//	//int iResult;
//	//(E, 3, 3, U, S, Vt);
//	Matrix_Transpose(U, 3, 3, Ut);
//
//	float V[] = { 0,0,1,PI / 2.f };
//	float Sigma[3 * 3] = { S[0],0,0,
//						0,S[1],0,
//						0,0,S[2] };
//	Rotation_Vector_2_Matrix_3D(V, R_z_t_1);
//	V[3] = -V[3];
//	Rotation_Vector_2_Matrix_3D(V, R_z_t_2);
//
//	Matrix_Multiply(U, 3, 3, Sigma, 3, Temp);
//	Matrix_Multiply(Temp, 3, 3, Vt, 3, Temp);
//	/*Disp(Temp, 3, 3,"Temp");
//	Disp(U, 3, 3, "U");
//	Disp(S, 3, 1, "S");
//	Disp(Vt, 3, 3, "Vt");*/
//
//	//svd_3(E, 3, 3, U, S, Vt);
//	Disp(U, 3, 3, "U");
//	//Disp(Sigma, 3, 3, "Sigma");
//	Disp(Vt, 3, 3, "Vt");
//	Matrix_Multiply(U, 3, 3, Sigma, 3, Temp);
//	Matrix_Multiply(Temp, 3, 3, Vt, 3, Temp);
//	Disp(Temp, 3, 3, "Temp");
//
//	printf("%d", bIs_Orthogonal(U, 3));
//	printf("%d", bIs_Orthogonal(Vt, 3));
//
//	return;
//}

template<typename _T> void _svd_3(_T At[], int m, int n, int n1, _T Sigma[], _T Vt[], int* pbSuccess, double eps)
{
	int i, j, k, iter, max_iter = 30;
	_T sd;
	_T s, c;
	int iMax_Size = Max(m, n);
	memset(Sigma, 0, iMax_Size * sizeof(_T));
	memset(Vt, 0, m * n * sizeof(_T));

	int astep = m, vstep = n;

	for (i = 0; i < n; i++)
	{
		for (k = 0, sd = 0; k < m; k++)
		{
			_T t = At[i * astep + k];
			sd += (_T)t * t;
		}
		Sigma[i] = sd;

		if (Vt)
		{
			for (k = 0; k < n; k++)
				Vt[i * vstep + k] = 0;
			Vt[i * vstep + i] = 1;
		}
	}
	//Disp(Vt, 3, vstep, "Vt");
	for (iter = 0; iter < max_iter; iter++)
	{
		int changed = false;
		for (i = 0; i < n - 1; i++)
		{
			for (j = i + 1; j < n; j++)
			{
				_T* Ai = At + i * astep, * Aj = At + j * astep;
				_T a = Sigma[i], p = 0, b = Sigma[j];

				for (k = 0; k < m; k++)
					p += (_T)Ai[k] * Aj[k];

				if (std::abs(p) <= eps * std::sqrt((double)a * b))
					continue;

				p *= 2;
				double beta = a - b, gamma = hypot((double)p, beta);
				if (beta < 0)
				{
					double delta = (gamma - beta) * 0.5;
					s = (_T)std::sqrt(delta / gamma);
					c = (_T)(p / (gamma * s * 2));
				}
				else
				{
					c = (_T)std::sqrt((gamma + beta) / (gamma * 2));
					s = (_T)(p / (gamma * c * 2));
				}

				a = b = 0;
				for (k = 0; k < m; k++)
				{
					_T t0 = c * Ai[k] + s * Aj[k];
					_T t1 = -s * Ai[k] + c * Aj[k];
					Ai[k] = t0; Aj[k] = t1;
					a += (_T)t0 * t0; b += (_T)t1 * t1;
				}
				Sigma[i] = a; Sigma[j] = b;

				changed = true;

				if (Vt)
				{
					_T* Vi = Vt + i * vstep, * Vj = Vt + j * vstep;
					k = 0;	//vblas.givens(Vi, Vj, n, c, s);

					for (; k < n; k++)
					{
						_T t0 = c * Vi[k] + s * Vj[k];
						_T t1 = -s * Vi[k] + c * Vj[k];
						Vi[k] = t0; Vj[k] = t1;
					}
					//printf("iter:%d i:%d j:%d Vt[0]:%f\n", iter, i, j, Vt[0]);
				}
			}
		}
		if (!changed)
			break;
	}

	for (i = 0; i < n; i++)
	{
		for (k = 0, sd = 0; k < m; k++)
		{
			_T t = At[i * astep + k];
			sd += (_T)t * t;
		}
		Sigma[i] = std::sqrt(sd);
	}

	for (i = 0; i < n - 1; i++)
	{
		j = i;
		for (k = i + 1; k < n; k++)
		{
			if (Sigma[j] < Sigma[k])
				j = k;
		}
		if (i != j)
		{
			std::swap(Sigma[i], Sigma[j]);
			if (Vt)
			{
				for (k = 0; k < m; k++)
					std::swap(At[i * astep + k], At[j * astep + k]);

				for (k = 0; k < n; k++)
					std::swap(Vt[i * vstep + k], Vt[j * vstep + k]);
			}
		}
	}
	if (!Vt)
		return;
	//Disp(Vt, m, vstep,"Vt");

	unsigned long long iRandom_State = 0x12345678;
	for (i = 0; i < n1; i++)
	{
		sd = i < n ? Sigma[i] : 0;

		for (int ii = 0; ii < 100 && sd <= DBL_MIN; ii++)
		{
			// if we got a zero singular value, then in order to get the corresponding left singular vector
			// we generate a random vector, project it to the previously computed left singular vectors,
			// subtract the projection and normalize the difference.
			const _T val0 = (_T)(1. / m);
			for (k = 0; k < m; k++)
			{
				iGet_Random_No_cv(&iRandom_State);
				_T val = (iRandom_State & 256) != 0 ? val0 : -val0;
				At[i * astep + k] = val;
			}
			for (iter = 0; iter < 2; iter++)
			{
				for (j = 0; j < i; j++)
				{
					sd = 0;
					for (k = 0; k < m; k++)
						sd += At[i * astep + k] * At[j * astep + k];
					_T asum = 0;
					for (k = 0; k < m; k++)
					{
						_T t = (_T)(At[i * astep + k] - sd * At[j * astep + k]);
						At[i * astep + k] = t;
						asum += std::abs(t);
					}
					asum = asum > eps * 100 ? 1 / asum : 0;
					for (k = 0; k < m; k++)
						At[i * astep + k] *= asum;
				}
			}
			sd = 0;
			for (k = 0; k < m; k++)
			{
				_T t = At[i * astep + k];
				sd += (_T)t * t;
			}
			sd = std::sqrt(sd);
		}

		s = (_T)(sd > DBL_MIN ? 1 / sd : 0.);
		for (k = 0; k < m; k++)
			At[i * astep + k] *= s;
	}
}

template<typename _T> void svd_3(_T* A, SVD_Info oSVD, int* pbSuccess, double eps)
{/*��дһ��SVD������Ҫ�㣺
1�����еľ���ͨ��ת��ͳһΪ		nnnnn	���߱ȿ�����״
									nnnnn
									nnnnn
									nnnnn
									nnnnn
2, ����ԭ�����Ѿ��Ǹ߱ȿ�����״�������ϣ�����Ahxw, ��ʱVtֻ��wxw, U����wxh
3���������� mmmmmmmmmmmmmmm	����ȸߴ�ľ��󣬰�ʱVtֻ���� (h+1)*w, U����wxw
			mmmmmmmmmmmmmmm
			mmmmmmmmmmmmmmm
4,���ڷ�����ʱU Vt���Ƿ���
5,S����ʱֻ����һ��min(m,n)
*/
	_T* A_1, * Vt_1 = NULL, * Sigma;
	int bHor;
	int h = oSVD.h_A, w = oSVD.w_A;
	int iTemp, iMax_Size = Max(h, w), iMin_Size = Min(h, w);

	//���A����Ҫ���뵽��������������У���û�й����������ڴ氲��
	//A_1 = (_T*)pMalloc(&oMem_Mgr, iMax_Size * iMax_Size * sizeof(_T));
	A_1 = (_T*)pMalloc(&oMem_Mgr, iMax_Size * iMax_Size * 2 * sizeof(_T));
	Sigma = (_T*)pMalloc(&oMem_Mgr, iMax_Size * sizeof(_T));
	if (!A_1 || !Sigma)
	{
		oSVD.m_bSuccess = 0;
		goto END;
	}
	else
		oSVD.m_bSuccess = 1;
	memset(A_1, 0, iMax_Size * iMax_Size * sizeof(_T));
	if (h < w)
	{//Ҫת�ó�Ϊ�߱ȿ��
		bHor = 1;
		iTemp = h;
		h = w;
		w = iTemp;
		memcpy(A_1, A, h * w * sizeof(_T));
	}
	else
	{//���� ��ʱ��A_1=At
		Matrix_Transpose(A, h, w, A_1);
		bHor = 0;
	}

	Vt_1 = (_T*)pMalloc(&oMem_Mgr, w * (w + 1) * sizeof(_T));
	memset(Vt_1, 0, w * (w + 1) * sizeof(_T));
	//U�������ţ�����ֱ��д��U��
	//Disp(A, 3, 3, "A");
	//Disp(A_1, 3, 3, "A_1");
	_svd_3(A_1, h, w, w + 1, Sigma, Vt_1, &oSVD.m_bSuccess, eps);
	//Disp(A_1, h, h, "U");
	//Disp(Sigma,1,w,"S");
	//Disp(Vt_1, w, w, "Vt");
	//Disp(&A_1[2*3], 1, oSVD.h_A);
	memcpy(oSVD.S, Sigma, iMin_Size * sizeof(_T));
	if (!bHor)
	{//��A_1ת�þ���U��Ȼ������ʱUֻȡ������	mmmm
	//											mmmm
	//											mmmm
	//											mmmm
		if (oSVD.U)
			Matrix_Transpose(A_1, w, h, (_T*)oSVD.U);
		//Disp((_T*)oInfo.U, oInfo.h_Min_U, oInfo.w_Min_U, "U");
		if (oSVD.Vt)
			memcpy(oSVD.Vt, Vt_1, oSVD.h_Min_Vt * oSVD.h_Min_Vt * sizeof(_T));
	}
	else
	{//
		//��ʱ��Vt��ֵ��U
		//Disp((_T*)oInfo.Vt, h, h);
		Matrix_Transpose(Vt_1, w, w, (_T*)oSVD.U);
		memcpy(oSVD.Vt, A_1, (w + 1) * h * sizeof(_T));
	}
	if (pbSuccess)
		*pbSuccess = oSVD.m_bSuccess;
END:
	if (A_1)
		Free(&oMem_Mgr, A_1);
	if (Sigma)
		Free(&oMem_Mgr, Sigma);
	if (Vt_1)
		Free(&oMem_Mgr, Vt_1);
}
template<typename _T> void Matrix_Multiply(_T* A, int ma, int na, _T a, _T* C)
{//����A����һ������
	int i, iSize = ma * na;
	for (i = 0; i < iSize; i++)
		C[i] = A[i] * a;
	return;
}
template<typename _T>void Matrix_Multiply_3x1(_T A[3 * 3], _T B[3], _T C[3])
{//�����������������3x1������
	if (B == C)
	{
		_T C1[3];
		C1[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
		C1[1] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
		C1[2] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
		C[0] = C1[0], C[1] = C1[1], C[2] = C1[2];
	}else
	{
		C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
		C[1] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
		C[2] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
	}
	return;
}
template<typename _T> void Matrix_Multiply(_T* A, int ma, int na, _T* B, int nb, _T* C)
{//Amn x Bno = Cmo
	int y, x, i;
	_T fValue, * C_Dup;	// = (_T*)malloc(ma * nb * sizeof(_T));
	//Light_Ptr oPtr = oMem_Mgr;
	//Malloc_1(oPtr, ma * nb * sizeof(_T), C_Dup);
	if (C == A || C == B)
	{
		C_Dup = (_T*)pMalloc(&oMem_Mgr, ma * nb * sizeof(_T));
		if (!C_Dup)
		{
			printf("Fail to Malloc_1 in Matrix_Multiply\n");
			return;
		}
	}else
		C_Dup = C;
	
	for (y = 0; y < ma; y++)
	{
		for (x = 0; x < nb; x++)
		{
			/*if( (y==0 && x==1 ) || (y==1 && x==0))
				printf("here");*/
			for (fValue = 0, i = 0; i < na; i++)
				fValue += A[y * na + i] * B[i * nb + x];
			C_Dup[y * nb + x] = fValue;
		}
	}
	if (C == A || C == B)
	{
		memcpy(C, C_Dup, ma * nb * sizeof(_T));
		Free(&oMem_Mgr, C_Dup);
	}	
	return;
}

void SB_Matrix()
{//����Ǹ�ɵ�Ʒ�����������ƭtemplate
	fGet_Random_No((double)0, (double)0);
	fGet_Random_No((float)0, (float)0);

	bTest_Eigen((double*)NULL,0, (double)0.f, (double*)NULL);
	bTest_Eigen((float*)NULL, 0, (float)0.f, (float*)NULL);

	bInverse_Power((double(*)[3])NULL, (double*)NULL, (double*)NULL);
	bInverse_Power((float(*)[3])NULL, (float*)NULL, (float*)NULL);

	bInverse_Power((double*)NULL, 0, (double*)NULL, (double*)NULL);
	bInverse_Power((float*)NULL, 0, (float*)NULL, (float*)NULL);

	bTest_LU((double*)NULL, (double*)NULL, (double*)NULL, 0);
	bTest_LU((float*)NULL, (float*)NULL, (float*)NULL, 0);

	Test_Inv_Matrix((double*)NULL, (double*)NULL,0);
	Test_Inv_Matrix((float*)NULL, (float*)NULL,0);

	LU_Decompose((double*)NULL, (double*)NULL, (double*)NULL, 0);
	LU_Decompose((float*)NULL, (float*)NULL, (float*)NULL, 0);

	Power_Method((double*)NULL, 0, (double*)NULL, (double*)NULL);
	Power_Method((float*)NULL, 0, (float*)NULL, (float*)NULL);

	Matrix_2_R((double*)NULL, (double*)NULL);
	Matrix_2_R((float*)NULL, (float*)NULL);

	bIs_R((double*)NULL);
	bIs_R((float*)NULL);

	fGet_Tr((double*)NULL, 0);
	fGet_Tr((float*)NULL, 0);

	Gen_Rotation_Matrix((double*)NULL, (double)0.f, (double*)NULL);
	Gen_Rotation_Matrix((float*)NULL, (float)0.f, (float*)NULL);

	Softmax((double*)NULL, 0, (double*)NULL);
	Softmax((float*)NULL, 0, (float*)NULL);

	Cholosky_Decompose(Sparse_Matrix<double>{}, (Sparse_Matrix<double>*) NULL);
	Cholosky_Decompose(Sparse_Matrix<float>{}, (Sparse_Matrix<float>*) NULL);

	Test_Linear((double*)NULL, 0, (double*)NULL, (double*)NULL);
	Test_Linear((float*)NULL, 0, (float*)NULL, (float*)NULL);

	Solve_Linear_Gause_AAt((double*)NULL, 0, (double*)NULL, (double*)NULL, NULL);
	Solve_Linear_Gause_AAt((float*)NULL, 0, (float*)NULL, (float*)NULL, NULL);

	Solve_Linear_Jocabi((double*)NULL, (double*)NULL, 0, (double*)NULL, NULL);
	
	Solve_Linear_Gauss_Seidel((double*)NULL, (double*)NULL, 0, (double*)NULL, NULL);

	Solve_Linear_LLt((double*)NULL, 0, (double*)NULL, (double*)NULL);
	Solve_Linear_LLt((float*)NULL, 0, (float*)NULL, (float*)NULL);

	fLinear_Equation_Check((double*)NULL, (int)0, (double*)NULL, (double*)NULL,(double)0.f);
	fLinear_Equation_Check((float*)NULL, (int)0, (float*)NULL, (float*)NULL,(float)0.f);

	LLt_Decompose((double*)NULL, 0, (double*)NULL);
	LLt_Decompose((float*)NULL, 0, (float*)NULL);

	Cholosky_Decompose((double*)NULL,0, (double*)NULL);
	Cholosky_Decompose((float*)NULL, 0, (float*)NULL);

	bIs_Symmetric_Matrix(Sparse_Matrix<double>{});
	bIs_Symmetric_Matrix(Sparse_Matrix<float>{});

	Quaternion_2_Rotation_Matrix((double*)NULL, (double*)NULL);
	Quaternion_2_Rotation_Matrix((float*)NULL, (float*)NULL);

	Quaternion_2_Rotation_Vector((double*)NULL, (double*)NULL);
	Quaternion_2_Rotation_Vector((float*)NULL, (float*)NULL);

	bIs_Unit_Vector((double*)NULL, 0);
	bIs_Unit_Vector((float*)NULL, 0);

	bIs_Symmetric_Matrix<double>(NULL, 0);
	bIs_Symmetric_Matrix<float>(NULL, 0);
	bIs_Symmetric_Matrix<int>(NULL, 0);

	Crop_Matrix((double*)NULL, 0, 0, 0, 0, 0, 0, (double*)NULL, 0);
	Crop_Matrix((float*)NULL, 0, 0, 0, 0, 0, 0, (float*)NULL, 0);

	Add_I_Matrix((double*)NULL, 0);
	Add_I_Matrix((float*)NULL, 0);

	Disp_Fillness(Sparse_Matrix<double>{});
	Disp_Fillness(Sparse_Matrix<float>{});

	Disp_Fillness((double*)NULL, 0, 0);
	Disp_Fillness((float*)NULL, 0, 0);

	Matrix_Add_Partial((double*)NULL, 0, 0, 0,0,0, (double*)NULL, 0, 0, 0);
	Matrix_Add_Partial((float*)NULL, 0, 0, 0,0,0, (float*)NULL, 0, 0, 0);

	Copy_Matrix_Partial((double*)NULL, 0, 0, 0, (double*)NULL, 0, 0);
	Copy_Matrix_Partial((float*)NULL, 0, 0, 0, (float*)NULL, 0, 0);

	Copy_Matrix_Partial((double*)NULL, 0, 0, (double*)NULL, 0, 0, 0);
	Copy_Matrix_Partial((float*)NULL, 0, 0, (float*)NULL, 0, 0, 0);

	Rotation_Vector_3_2_4((double*)NULL, (double*)NULL);
	Rotation_Vector_3_2_4((float*)NULL, (float*)NULL);

	Rotation_Vector_4_2_3((double*)NULL, (double*)NULL);
	Rotation_Vector_4_2_3((float*)NULL, (float*)NULL);

	Re_Arrange_Sparse_Matrix((Sparse_Matrix<double>*)NULL);
	Re_Arrange_Sparse_Matrix((Sparse_Matrix<float>*)NULL);

	fGet_Value((Sparse_Matrix<double>*)NULL, 0, 0);
	fGet_Value((Sparse_Matrix<float>*)NULL, 0, 0);

	Disp(Sparse_Matrix<double>{});
	Disp(Sparse_Matrix<float>{});

	Disp_Part<float>(NULL, 0, 0, 0, 0, 0);
	Disp_Part<unsigned char>(NULL, 0, 0, 0, 0, 0);
	Disp_Part<float>(NULL, 0, 0, 0, 0, 0);
	Disp_Part<double>(NULL, 0, 0, 0, 0, 0);
	Disp_Part<int>(NULL, 0, 0, 0, 0, 0);
	Disp_Part<unsigned int>(NULL, 0, 0, 0, 0, 0);
	Disp_Part<short>(NULL, 0, 0, 0, 0, 0);
	Disp_Part<unsigned short>(NULL, 0, 0, 0, 0, 0);


	Reset_Sparse_Matrix((Sparse_Matrix<double>*)NULL);
	Reset_Sparse_Matrix((Sparse_Matrix<float>*)NULL);

	Matrix_Transpose_1(Sparse_Matrix<double>{}, (Sparse_Matrix<double>*)NULL);
	Matrix_Transpose_1(Sparse_Matrix<float>{}, (Sparse_Matrix<float>*)NULL);

	Matrix_Minus(Sparse_Matrix<double>{}, Sparse_Matrix<double>{}, (Sparse_Matrix<double>*)NULL);
	Matrix_Minus(Sparse_Matrix<float>{}, Sparse_Matrix<float>{}, (Sparse_Matrix<float>*)NULL);

	Matrix_Add(Sparse_Matrix<double>{}, Sparse_Matrix<double>{}, (Sparse_Matrix<double>*)NULL);
	Matrix_Add(Sparse_Matrix<float>{}, Sparse_Matrix<float>{}, (Sparse_Matrix<float>*)NULL);

	Add_I_Matrix((Sparse_Matrix<double>*)NULL);
	Add_I_Matrix((Sparse_Matrix<float>*)NULL);
	
	Matrix_Multiply(Sparse_Matrix<double>{}, Sparse_Matrix<double>{}, (Sparse_Matrix<double>*)NULL);
	Matrix_Multiply(Sparse_Matrix<float>{}, Sparse_Matrix<float>{}, (Sparse_Matrix<float>*)NULL);

	Matrix_Multiply(Sparse_Matrix<double>{},(double)0);
	Matrix_Multiply(Sparse_Matrix<float>{}, (float)0);

	Sparse_2_Dense(Sparse_Matrix<double>{}, (double*)NULL);
	Sparse_2_Dense(Sparse_Matrix<float>{}, (float*)NULL);

	Dense_2_Sparse((double*)NULL, 0, 0, (Sparse_Matrix<double>*)NULL);
	Dense_2_Sparse((float*)NULL, 0, 0, (Sparse_Matrix<float>*)NULL);

	Solve_Linear_Contradictory((double*)NULL, 0, 0, (double*)NULL, (double*)NULL, NULL);
	Solve_Linear_Contradictory((float*)NULL, 0, 0, (float*)NULL, (float*)NULL, NULL);

	Solve_Linear_Gause(Sparse_Matrix<double>{}, (double*)NULL, (double*)NULL);
	Solve_Linear_Gause(Sparse_Matrix<float>{}, (float*)NULL, (float*)NULL);
	
	Solve_Linear_Gause_2<double>(Sparse_Matrix<double>{}, (double*)NULL, (double*)NULL);
	Solve_Linear_Gause_2<float>(Sparse_Matrix<float>{}, (float*)NULL, (float*)NULL);

	Solve_Linear_Gause_1(Sparse_Matrix<double>{}, (double*)NULL, (double*)NULL);
	Solve_Linear_Gause_1(Sparse_Matrix<float>{}, (float*)NULL, (float*)NULL);

	Gen_Rotation_Matrix_2D((double*)NULL, 0);
	Gen_Rotation_Matrix_2D((float*)NULL, 0);

	Gen_Homo_Matrix_2D((double*)NULL, (double*)NULL, (double*)NULL);
	Gen_Homo_Matrix_2D((float*)NULL, (float*)NULL, (float*)NULL);

	//Init_Sparse_Matrix((Sparse_Matrix<double>*)NULL, 0, 0);
	//Init_Sparse_Matrix((Sparse_Matrix<float>*)NULL, 0, 0);

	Compact_Sparse_Matrix((Sparse_Matrix<double>*)NULL);
	Compact_Sparse_Matrix((Sparse_Matrix<float>*)NULL);

	Free_Sparse_Matrix((Sparse_Matrix<double>*)NULL);
	Free_Sparse_Matrix((Sparse_Matrix<float>*)NULL);

	Schmidt_Orthogon((double*)NULL, 0, 0, (double*)NULL);
	Schmidt_Orthogon((float*)NULL, 0, 0, (float*)NULL);

	Elementary_Row_Operation((double*)NULL, 0, 0, (double*)NULL);
	Elementary_Row_Operation((float*)NULL, 0, 0, (float*)NULL);

	Solve_Linear_Solution_Construction_1((double*)NULL, 0, 0, (double*)NULL, NULL);
	Solve_Linear_Solution_Construction_1((float*)NULL, 0, 0, (float*)NULL, NULL);

	Solve_Linear_Solution_Construction((double*)NULL, 0, 0, (double*)NULL, NULL);
	Solve_Linear_Solution_Construction((float*)NULL, 0, 0, (float*)NULL, NULL);

	Get_R_t((double*)NULL);
	Get_R_t((float*)NULL);

	Gen_Roation_Matrix_2D((double*)NULL, (double)0, (double*)NULL);
	Gen_Roation_Matrix_2D((float*)NULL, (float)0, (float*)NULL);

	se3_2_SE3((double*)NULL, (double*)NULL);
	se3_2_SE3((float*)NULL, (float*)NULL);

	SE3_2_se3((double*)NULL, (double*)NULL);
	SE3_2_se3((float*)NULL, (float*)NULL);

	SE3_2_se3((double*)NULL, (double*)NULL, (double*)NULL);
	SE3_2_se3((float*)NULL, (float*)NULL, (float*)NULL);

	Solve_Linear_Gause((double*)NULL, 0, (double*)NULL, (double*)NULL, NULL);
	Solve_Linear_Gause((float*)NULL, 0, (float*)NULL, (float*)NULL, NULL);

	Transpose_Multiply((double*)NULL, 0, 0, (double*)NULL);
	Transpose_Multiply((float*)NULL, 0, 0, (float*)NULL);

	Gen_Ksi_by_Rotation_Vector_t((double*)NULL, (double*)NULL, (double*)NULL);
	Gen_Ksi_by_Rotation_Vector_t((float*)NULL, (float*)NULL, (float*)NULL);

	Get_Jl_4((double*)NULL, (double*)NULL);
	Get_Jl_4((float*)NULL, (float*)NULL);

	Get_Jl_3((double*)NULL, (double*)NULL);
	Get_Jl_3((float*)NULL, (float*)NULL);

	Get_Jr_4((double*)NULL, (double*)NULL);
	Get_Jr_4((float*)NULL, (float*)NULL);

	Get_Jr_3((double*)NULL, (double*)NULL);
	Get_Jr_3((float*)NULL, (float*)NULL);

	Exp_Ref((double*)NULL, 0, (double*)NULL);
	Exp_Ref((float*)NULL, 0, (float*)NULL);

	fGet_Distance((double*)NULL, (double*)NULL, 0);
	fGet_Distance((float*)NULL, (float*)NULL, 0);

	Vector_Add((double*)NULL, (double*)NULL, 0, (double*)NULL);
	Vector_Add((float*)NULL, (float*)NULL, 0, (float*)NULL);

	Vector_Multiply((double*)NULL,0, (double)0, (double*)NULL);
	Vector_Multiply((float*)NULL,0, (float)0, (float*)NULL);

	Gen_Translation_Matrix((double*)NULL, (double*)NULL);
	Gen_Translation_Matrix((float*)NULL, (float*)NULL);

	iGet_Rank((double*)NULL, 0, 0);
	iGet_Rank((float*)NULL, 0, 0);

	Matrix_Add((double*)NULL, (double*)NULL, 0, (double*)NULL);
	Matrix_Add((float*)NULL, (float*)NULL, 0, (float*)NULL);

	Matrix_Minus((double*)NULL, (double*)NULL, 0, (double*)NULL);
	Matrix_Minus((float*)NULL, (float*)NULL, 0, (float*)NULL);

	Rotation_Vector_4_2_Matrix((double*)NULL, (double*)NULL);
	Rotation_Vector_4_2_Matrix((float*)NULL, (float*)NULL);

	Rotation_Vector_3_2_Matrix((double*)NULL, (double*)NULL);
	Rotation_Vector_3_2_Matrix((float*)NULL, (float*)NULL);

	Vector_Minus((double*)NULL, (double*)NULL, 0, (double*)NULL);
	Vector_Minus((float*)NULL, (float*)NULL, 0, (float*)NULL);

	Hat<double>(NULL, NULL);
	Hat<float>(NULL, NULL);

	Vee((double*)NULL, (double*)NULL);
	Vee((float*)NULL, (float*)NULL);

	Lie_Bracket((double*)NULL, (double*)NULL, (double*)NULL);
	Lie_Bracket((float*)NULL, (float*)NULL, (float*)NULL);

	Rotation_Vector_2_Quaternion((double*)NULL, (double*)NULL);
	Rotation_Vector_2_Quaternion((float*)NULL, (float*)NULL);

	Rotation_Matrix_2_Quaternion((double*)NULL, (double*)NULL);
	Rotation_Matrix_2_Quaternion((float*)NULL, (float*)NULL);

	fDot((double*)NULL, (double*)NULL, 0);
	fDot((float*)NULL, (float*)NULL, 0);

	Cross_Product((double*)NULL, (double*)NULL, (double*)NULL);
	Cross_Product((float*)NULL, (float*)NULL, (float*)NULL);

	Homo_Normalize((double*)NULL, 0, (double*)NULL);
	Homo_Normalize((float*)NULL, 0, (float*)NULL);

	Get_Inv_Matrix_Row_Op(Sparse_Matrix<double>{}, (Sparse_Matrix<double>*)NULL);
	Get_Inv_Matrix_Row_Op(Sparse_Matrix<float>{}, (Sparse_Matrix<float>*)NULL);

	Get_Inv_Matrix_Row_Op((float*)NULL, (float*)NULL, 0, NULL);
	Get_Inv_Matrix_Row_Op((double*)NULL, (double*)NULL, 0, NULL);

	Get_Inv_AAt_Row_Op((float*)NULL, (float*)NULL, 0, NULL);
	Get_Inv_AAt_Row_Op((double*)NULL, (double*)NULL, 0, NULL);

	Get_Inv_AAt_3x3((float*)NULL, (float*)NULL, NULL);
	Get_Inv_AAt_3x3((double*)NULL, (double*)NULL, NULL);

	Matrix_Multiply_3x1((float*)NULL, (float*)NULL, (float*)NULL);
	Matrix_Multiply_3x1((double*)NULL, (double*)NULL, (double*)NULL);

	Matrix_Multiply((float*)NULL, 0, 0, (float*)NULL, 0, (float*)NULL);
	Matrix_Multiply((double*)NULL, 0, 0, (double*)NULL, 0, (double*)NULL);
	Matrix_Multiply((int*)NULL, 0, 0, (int*)NULL, 0, (int*)NULL);

	Matrix_Multiply((float*)NULL, 0, 0, (float)0, (float*)NULL);
	Matrix_Multiply((double*)NULL, 0, 0, (double)0, (double*)NULL);

	bIs_Orthogonal((float*)NULL, 0);
	bIs_Orthogonal((double*)NULL, 0);

	svd_3((float*)NULL, {});
	svd_3((double*)NULL, {});

	Test_SVD((double*)NULL, {});
	Test_SVD((float*)NULL, {});

	//SVD_Allocate((double*)NULL, 0, 0, NULL);
	//SVD_Allocate((float*)NULL, 0, 0, NULL);

	SVD_Alloc(0, 0, NULL, (double*)NULL);
	SVD_Alloc(0, 0, NULL, (float*)NULL);

	Elementary_Row_Operation_1((double*)NULL, 0, 0, (double*)NULL);
	Elementary_Row_Operation_1((float*)NULL, 0, 0, (float*)NULL);

	fGet_Determinant((double*)NULL, 0);
	fGet_Determinant((float*)NULL, 0);

	Normalize((double*)NULL, 0, (double*)NULL);
	Normalize((float*)NULL, 0, (float*)NULL);

	fGet_Mod((double*)NULL, 0);
	fGet_Mod((float*)NULL, 0);

	Gen_Homo_Matrix((double*)NULL, (double*)NULL, (double*)NULL);
	Gen_Homo_Matrix((float*)NULL, (float*)NULL, (float*)NULL);

	Gen_Homo_Matrix_1((double*)NULL, (double*)NULL, (double*)NULL);
	Gen_Homo_Matrix_1((float*)NULL, (float*)NULL, (float*)NULL);

	Rotation_Matrix_2_Vector((double*)NULL, (double*)NULL);
	Rotation_Matrix_2_Vector((float*)NULL, (float*)NULL);

	Rotation_Matrix_2_Vector_3((double*)NULL, (double*)NULL);
	Rotation_Matrix_2_Vector_3((float*)NULL, (float*)NULL);

	Rotation_Matrix_2_Vector_4((double*)NULL, (double*)NULL);
	Rotation_Matrix_2_Vector_4((float*)NULL, (float*)NULL);

	Gen_I_Matrix((double*)NULL, 0, 0);
	Gen_I_Matrix((float*)NULL, 0, 0);

	Matrix_Transpose<double>(NULL, 0, 0, NULL);
	Matrix_Transpose<float>(NULL, 0, 0, NULL);
}
template<typename _T>void Test_SVD(_T A[], SVD_Info oSVD, int* piResult, double eps)
{//��һ�ּ򻯱�ʾ���˴�����һ�·ֽ����Ƿ����Ԥ��
	int y, x, iResult = 1;
	_T* A_1 = NULL;
	union {
		_T* S;
		_T* SVt;
		_T* US;
	};
	if (oSVD.h_A > oSVD.w_A)
	{//���Σ�A= U x (S x Vt) ������棬����ǰ��
		SVt = (_T*)pMalloc(&oMem_Mgr, oSVD.w_Min_S * oSVD.w_Min_S * sizeof(_T));
		memset(S, 0, oSVD.w_Min_S * oSVD.w_Min_S * sizeof(_T));
		for (y = 0; y < oSVD.w_Min_S; y++)
			S[y * oSVD.w_Min_S + y] = ((_T*)oSVD.S)[y];
		Disp(S, oSVD.w_Min_S, oSVD.w_Min_S, "S");
		Matrix_Multiply(S, oSVD.w_Min_S, oSVD.w_Min_S, (_T*)oSVD.Vt, oSVD.h_Min_Vt, SVt);
		Disp(SVt, oSVD.w_Min_S, oSVD.w_Min_S, "SVt");
		A_1 = (_T*)pMalloc(&oMem_Mgr, oSVD.h_A * oSVD.w_A * sizeof(_T));
		Matrix_Multiply((_T*)oSVD.U, oSVD.h_Min_U, oSVD.w_Min_U, SVt, oSVD.w_Min_Vt, A_1);
		Disp(A_1, oSVD.h_A, oSVD.w_A, "A_1");
		Free(&oMem_Mgr, SVt);
	}
	else
	{//����
		S = (_T*)pMalloc(&oMem_Mgr, oSVD.w_Min_S * oSVD.w_Min_S * sizeof(_T));
		memset(S, 0, oSVD.w_Min_S * oSVD.w_Min_S * sizeof(_T));
		for (y = 0; y < oSVD.w_Min_S; y++)
			((_T*)S)[y * oSVD.w_Min_S + y] = ((_T*)oSVD.S)[y];
		//Disp(S, oSVD.w_Min_S, oSVD.w_Min_S, "S");
		Matrix_Multiply((_T*)oSVD.U, oSVD.h_Min_U, oSVD.w_Min_U, S, oSVD.w_Min_S, US);
		Disp(US, oSVD.h_Min_U, oSVD.w_Min_S, "UxS");
		A_1 = (_T*)pMalloc(&oMem_Mgr, oSVD.h_A * oSVD.w_A * sizeof(_T));
		Matrix_Multiply(US, oSVD.h_Min_U, oSVD.w_Min_S, (_T*)oSVD.Vt, oSVD.w_Min_Vt, A_1);
		Disp(A_1, oSVD.h_A, oSVD.w_A, "A_1");
		Free(&oMem_Mgr, US);
	}
	if (bIs_Orthogonal((_T*)oSVD.U, oSVD.h_Min_U, oSVD.w_Min_U))
		printf("U ����\n");
	else
		printf("U ������\n");

	if (bIs_Orthogonal((_T*)oSVD.Vt, oSVD.h_Min_Vt, oSVD.w_Min_Vt))
		printf("Vt ����\n");
	else
		printf("Vt ������\n");

	for (y = 0; y < oSVD.h_A; y++)
	{
		for (x = 0; x < oSVD.w_A; x++)
		{
			if (abs(A_1[y * oSVD.w_A + x] - A[y * oSVD.w_A + x]) > eps)
			{
				printf("Correct:%f Error:%f\n", A[y * oSVD.w_A + x], A_1[y * oSVD.w_A + x]);
				iResult = 0;
			}
		}
	}

	if (A_1)
		Free(&oMem_Mgr, A_1);
	if (piResult)
		*piResult = iResult;
}
template<typename _T>void Test_SVD(_T A[], int h, int w, _T U[], _T S[], _T Vt[], int* piResult, double eps)
{//����SVD�ֽ���
	_T* pTemp = (_T*)pMalloc(&oMem_Mgr, h * w * sizeof(_T));
	int y, x, iResult;
	Disp(S, h, w, "S,���Ĺ۲�����ֵ,����0ֵ��λ�ÿ�������Vt��Ӧ����");
	Matrix_Multiply(U, h, h, S, w, pTemp);
	Disp(pTemp, h, w, "UxS");
	Matrix_Multiply(pTemp, h, w, Vt, w, pTemp);
	Disp(pTemp, h, w, "USVt");
	iResult = 1;
	if (bIs_Orthogonal(U, h))
		printf("U ����\n");
	else
		printf("U ������\n");

	if (bIs_Orthogonal(Vt, w))
		printf("Vt ����\n");
	else
		printf("Vt ������\n");

	for (y = 0; y < h; y++)
	{
		for (x = 0; x < w; x++)
		{
			if (abs(pTemp[y * w + x] - A[y * w + x]) > eps)
			{
				printf("Correct:%f Error:%f\n", A[y * w + x], pTemp[y * w + x]);
				iResult = 0;
				//goto END;
			}
		}
	}

	Free(&oMem_Mgr, pTemp);
	if (piResult)
		*piResult = iResult;
	return;
}
void Free_SVD(SVD_Info* poInfo)
{
	Free(&oMem_Mgr, poInfo->U);
	Free(&oMem_Mgr, poInfo->S);
	Free(&oMem_Mgr, poInfo->Vt);
}
template<typename _T>void SVD_Alloc(int h, int w, SVD_Info* poInfo, _T* A)
{//����ӿڣ����㻻�����, ע�⣬�÷�Ϊ SVD_Alloc<_T>(...)
	_T* U, * S, * Vt;
	SVD_Info oInfo;
	oInfo.A = A;
	oInfo.h_A = h;
	oInfo.w_A = w;
	if (h >= w)
	{//��׼��
		U = (_T*)pMalloc(&oMem_Mgr, w * h * sizeof(_T));
		oInfo.h_Min_U = h;
		oInfo.w_Min_U = w;

		Vt = (_T*)pMalloc(&oMem_Mgr, w * w * sizeof(_T));
		oInfo.h_Min_Vt = oInfo.w_Min_Vt = w;
	}
	else
	{//
		U = (_T*)pMalloc(&oMem_Mgr, h * h * sizeof(_T));
		oInfo.h_Min_U = oInfo.w_Min_U = h;
		Vt = (_T*)pMalloc(&oMem_Mgr, (h + 1) * w * sizeof(_T));
		oInfo.h_Min_Vt = h + 1;
		oInfo.w_Min_Vt = w;
	}
	oInfo.U = U;
	oInfo.Vt = Vt;
	S = (_T*)pMalloc(&oMem_Mgr, oInfo.w_Min_S = min(h, w));
	oInfo.S = S;
	*poInfo = oInfo;
}
template<typename _T>void SVD_Allocate(_T* A, int h, int w, SVD_Info* poInfo)
{
	_T* U, * S, * Vt;
	SVD_Info oInfo;
	oInfo.A = A;
	oInfo.h_A = h;
	oInfo.w_A = w;
	if (h >= w)
	{//��׼��
		U = (_T*)pMalloc(&oMem_Mgr, w * h * sizeof(_T));
		oInfo.h_Min_U = h;
		oInfo.w_Min_U = w;

		Vt = (_T*)pMalloc(&oMem_Mgr, w * w * sizeof(_T));
		oInfo.h_Min_Vt = oInfo.w_Min_Vt = w;
	}
	else
	{//
		U = (_T*)pMalloc(&oMem_Mgr, h * h * sizeof(_T));
		oInfo.h_Min_U = oInfo.w_Min_U = h;
		Vt = (_T*)pMalloc(&oMem_Mgr, (h + 1) * w * sizeof(_T));
		oInfo.h_Min_Vt = h + 1;
		oInfo.w_Min_Vt = w;
	}
	oInfo.U = U;
	oInfo.Vt = Vt;
	S = (_T*)pMalloc(&oMem_Mgr, oInfo.w_Min_S = min(h, w));
	oInfo.S = S;
	*poInfo = oInfo;
}
template<typename _T>void Homo_Normalize(_T V0[], int n, _T V1[])
{//homogeneous normalization, ��ʵ���ǹ�һ��Ϊ������꣬ƨӪ��Ҳ��Ӵiiu
	_T fCoeff = V0[n - 1];
	for (int i = 0; i < n - 1; i++)
		V1[i] = V0[i] / fCoeff;
	return;
}
template<typename _T>void Gen_I_Matrix(_T M[], int h, int w)
{//����һ�����Ƶ�λ����ľ��󣬶Խ���Ϊ1
	int iMin = Min(h, w);
	memset(M, 0, h * w * sizeof(_T));
	for (int i = 0; i < iMin; i++)
		M[i * w + i] = 1;
	return;
}

void Free_Polynormial(Polynormial* poPoly)
{//�ͷŶ���ʽ��ռ�ڴ�
	if (poPoly)
	{
		if (poPoly->m_pCoeff)
			Free(&oMem_Mgr, poPoly->m_pCoeff);
		if (poPoly->m_pTerm)
			Free(&oMem_Mgr, poPoly->m_pTerm);
		*poPoly = {};
	}
}
void Copy_Polynormial(Polynormial* poSource, Polynormial* poDest)
{
	poDest->m_iElem_Count = poSource->m_iElem_Count;
	poDest->m_iMax_Term_Count = poSource->m_iMax_Term_Count;
	poDest->m_iTerm_Count = poSource->m_iTerm_Count;
	memcpy(poDest->m_pCoeff, poSource->m_pCoeff, poDest->m_iTerm_Count * sizeof(float));
	memcpy(poDest->m_pTerm, poSource->m_pTerm, poDest->m_iTerm_Count * poDest->m_iElem_Count * sizeof(unsigned char));
	return;
}
void Init_Polynormial(Polynormial* poPoly, int iElem_Count, int iMax_Term_Count)
{//��ʼ��һ������ʽ��iElem_Count: Ԫ��	iMax_Term_Count:����ж�����
	if (!poPoly)
		return;
	Polynormial oPoly = {};
	oPoly.m_iMax_Term_Count = iMax_Term_Count;
	oPoly.m_iElem_Count = iElem_Count;
	oPoly.m_iTerm_Count = 0;
	oPoly.m_pCoeff = (float*)pMalloc(&oMem_Mgr, iMax_Term_Count * sizeof(float));
	oPoly.m_pTerm = (unsigned char*)pMalloc(&oMem_Mgr, iMax_Term_Count * oPoly.m_iElem_Count * sizeof(unsigned char));
	*poPoly = oPoly;
	return;
}
unsigned char* pGet_Term(Polynormial* poPoly, int iTerm)
{//ȡ�ڼ�������ʽ
	return &poPoly->m_pTerm[iTerm * poPoly->m_iElem_Count * sizeof(unsigned char)];
}

void Add_Poly_Term(Polynormial* poPoly, float fCoeff, int first, ...)
{
	if (!poPoly)
		return;

	Polynormial oPoly = *poPoly;
	unsigned char* pCur_Term;
	unsigned char Input[64];
	int i, j;

	//�����������еı�������
	va_list arg_list;
	va_start(arg_list, first);
	Input[0] = first;
	for (int i = 1; i < oPoly.m_iElem_Count; i++)
		Input[i] = va_arg(arg_list, int);

	//���ȱ���һ�β���ͬ������кϲ�
	for (i = 0; i < oPoly.m_iTerm_Count; i++)
	{
		pCur_Term = pGet_Term(&oPoly, i);
		for (j = 0; j < oPoly.m_iElem_Count; j++)
			if (pCur_Term[j] != Input[j])
				break;
		if (j == oPoly.m_iElem_Count)
		{
			oPoly.m_pCoeff[i] += fCoeff;
			return;
		}
	}

	if (oPoly.m_iTerm_Count == oPoly.m_iMax_Term_Count)
	{
		printf("Exceeds max term count:%d\n", oPoly.m_iMax_Term_Count);
		return;
	}
	pCur_Term = pGet_Term(&oPoly, oPoly.m_iTerm_Count);
	for (i = 0; i < oPoly.m_iElem_Count; i++)
		pCur_Term[i] = Input[i];

	//���ܻ����Ÿ���
	oPoly.m_pCoeff[oPoly.m_iTerm_Count] = fCoeff;
	oPoly.m_iTerm_Count++;
	*poPoly = oPoly;
	return;
}
void Disp(Polynormial oPoly, int iElem_No)
{//��ʾ����ʽ. ���iElem_No=-1�� ����ʾȫ��
	if (iElem_No < 0)
	{
		printf("����%dԪ����ʽ\n", oPoly.m_iElem_Count);
		for (int i = 0; i < oPoly.m_iTerm_Count; i++)
		{
			Disp(oPoly, i);
			if (i < oPoly.m_iTerm_Count - 1)
				printf(" + ");
			printf("\n");
		}
	}
	else
	{
		unsigned char* pCur_Term = pGet_Term(&oPoly, iElem_No);
		if (oPoly.m_pCoeff[iElem_No] != 1.f)
		{
			if (round(oPoly.m_pCoeff[iElem_No]) == oPoly.m_pCoeff[iElem_No])
				printf("(%d)", (int)oPoly.m_pCoeff[iElem_No]);
			else
				printf("(%f)", oPoly.m_pCoeff[iElem_No]);
		}

		for (int i = 0; i < oPoly.m_iElem_Count; i++)
		{
			if (pCur_Term[i])
				if (pCur_Term[i] > 1)
					printf("x%d^%d ", i, pCur_Term[i]);
				else
					printf("x%d ", i);
		}
	}
}

void Get_Derivation(Polynormial* poSource, int iElem, Polynormial* poDest)
{//����ʽ�󵼡������һԪ�����䵼��������Ƕ�Ԫ������ƫ��
//iElem�� ��ʾ���Ǹ�����������
	Polynormial oPoly = *poSource;
	int i;
	unsigned char* pCur_Term;
	Init_Polynormial(&oPoly, poSource->m_iElem_Count, poSource->m_iMax_Term_Count);
	Copy_Polynormial(poSource, &oPoly);

	for (i = 0; i < oPoly.m_iTerm_Count; )
	{
		pCur_Term = pGet_Term(&oPoly, i);
		if (pCur_Term[iElem] != 0)
		{
			oPoly.m_pCoeff[i] *= pCur_Term[iElem];
			pCur_Term[iElem]--;
			i++;
		}
		else
		{//������ iElem�������ĵ���ʽ����������ȥ��
			memcpy(pCur_Term, pCur_Term + oPoly.m_iElem_Count, (oPoly.m_iTerm_Count - i - 1) * oPoly.m_iElem_Count * sizeof(unsigned char));
			memcpy(&oPoly.m_pCoeff[i], &oPoly.m_pCoeff[i + 1], (oPoly.m_iTerm_Count - i - 1) * sizeof(float));
			oPoly.m_iTerm_Count--;
		}
	}
	if (poDest)
		Copy_Polynormial(&oPoly, poDest);
	Free_Polynormial(&oPoly);
}
float fGet_Polynormial_Value(Polynormial oPoly, float x[])
{//����x,����ú�����ֵ
	int i, j;
	unsigned char* pCur_Term;
	float fTotal = 0, fValue;
	for (i = 0; i < oPoly.m_iTerm_Count; i++)
	{
		pCur_Term = pGet_Term(&oPoly, i);
		fValue = oPoly.m_pCoeff[i];
		for (j = 0; j < oPoly.m_iElem_Count; j++)
		{
			if (pCur_Term[j])
				fValue *= (float)pow(x[j], pCur_Term[j]);
		}
		fTotal += fValue;
	}
	return fTotal;
}

template<typename _T> void Exp_Ref(_T A[], int n, _T B[], _T eps)
{//�������ָ��, AΪn�׾���, ���nҲ��n�׾���, ���������㷨���Ӷ���ֱ�Ӽ��㵽����Ϊֹ
	//���������ûŪ����
	//exp(A) = ��(n=0,n->��) (A^n)/n!
	int iSize = n * n * sizeof(_T);
	int i;
	_T* pTemp_1 = (_T*)pMalloc(&oMem_Mgr, iSize),
		* pTemp_2 = (_T*)pMalloc(&oMem_Mgr, iSize),
		* pSum = (_T*)pMalloc(&oMem_Mgr, iSize);
	
	//memset(pSum, 0, iSize);
	Gen_I_Matrix(pSum, n, n);
	Gen_I_Matrix(pTemp_1, n, n);

	unsigned long long iN_Factorial=0;
	for (i = 1;; i++)
	{
		/*if (i == 28)
			printf("Here");*/
		Matrix_Multiply(pTemp_1, n, n, A, n, pTemp_1);		
		//_T F = (_T)1.f / (_T)iFactorial(i);
		iN_Factorial = iFactorial(iN_Factorial, i);
		_T F = (_T)1.f / (_T)iN_Factorial;
		if (isinf(F))
		{
			F = 0;
			break;
		}
			
		Matrix_Multiply(pTemp_1, n, n, F, pTemp_2);
		//Disp(pTemp_1, 3, 3, "pTemp_1");
		if (fGet_Mod(pTemp_2, n * n) < eps)
			break;
		Matrix_Add(pTemp_2, pSum, n, pSum);
		//printf("i:%d %lf\n",i, pSum[0]);
	}

	memcpy(B, pSum, iSize);
	Free(&oMem_Mgr, pTemp_1);
	Free(&oMem_Mgr, pTemp_2);
	Free(&oMem_Mgr, pSum);
}
template<typename _T> void Gen_Ksi_by_Rotation_Vector_t(_T Rotation_Vector[4], _T t[3], _T Ksi[6], int n)
{//��һ��Rotation_Vector ��һ��t����һ��Ksi, nΪRotation_Vector��ά����3ά��4ά�Կ�
	_T Rotation_Vector_1[4], Ksi_1[6];
	union {
		_T J[3 * 3];
		_T J_Inv[3 * 3];
	};
	int i, iResult;
	if (n == 3)
	{
		Normalize(Rotation_Vector, 3, Rotation_Vector_1);
		Rotation_Vector_1[3] = fGet_Mod(Rotation_Vector, 3);
		memcpy(&Ksi_1[3], Rotation_Vector, 3 * sizeof(_T));
	}
	else
	{
		memcpy(Rotation_Vector_1, Rotation_Vector, 4 * sizeof(_T));
		for (i = 0; i < 3; i++)
			Ksi_1[3 + i] = Rotation_Vector[i] * Rotation_Vector[3];
	}

	Get_Jl_4(Rotation_Vector_1, J);
	Get_Inv_Matrix_Row_Op(J, J_Inv, 3, &iResult);
	if (!iResult)
	{
		printf("Rank <3 in Gen_Ksi_by_Rotation_Vector_t\n");
		return;
	}

	//��=J(-1)t
	Matrix_Multiply(J_Inv, 3, 3, t, 1, Ksi_1);
	memcpy(Ksi, Ksi_1, 6 * sizeof(_T));
	return;
}
template<typename _T> void Transpose_Multiply(_T A[], int m, int n, _T B[], int bAAt)
{//iFlag=0 ʱ B = A'A iFlag=1 ʱ B= AA'
	_T* At;
	_T* B1;
	if (bAAt == 0)
	{
		At = (_T*)pMalloc(&oMem_Mgr, n * m * sizeof(_T));
		Matrix_Transpose(A, m, n, At);
		//Matrix_Multiply(At, n, m, A, n, B);	//
		Transpose_Multiply(At, n, m, B,1);	//AtA
		Free(&oMem_Mgr, At);
	}else
	{//�˴������Ż�һ�£���A�������ƽ�
		int x, y;
		//if (abs((int)(A - B)) < m*m)
		if (A == B)
			B1 = (_T*)pMalloc(m * m * sizeof(_T));
		else
			B1 = B;

		//�˴�Ӧ�����öԳ��Լ��ټ�����
		for (y = 0; y < m; y++)
		{
			for (x = y; x < m; x++)
			{	//(y,x)ΪĿ���ַ
				int iPos_r = y * n,
					iPos_c = x * n;
				_T fValue = 0;
				for(int x1=0;x1<n;x1++)
					fValue += A[iPos_r++] * A[iPos_c++];
				B1[x * m + y] = B1[y * m + x] = fValue;
			}
		}
		//int iPos_r,iPos_Dest;
		//for (iPos_Dest = iPos_r = y = 0; y < m; y++, iPos_r += n)
		//{//���ƽ�
		//	int y1,iPos_c;
		//	for (iPos_c = y1 = 0; y1 < m; y1++, iPos_c += n, iPos_Dest++)
		//	{
		//		int iPos_r1 = iPos_r,
		//			iPos_c1 = iPos_c;
		//		_T fValue = 0;
		//		for (x = 0; x < n; x++,iPos_r1++,iPos_c1++)
		//			fValue += A[iPos_r1] * A[iPos_c1];
		//		B1[iPos_Dest] = fValue;
		//	}
		//}

		//if (abs((int)(A - B)) < m * m)
		if (A == B)
			Free(B1);
		//Matrix_Multiply(A, m, n, At, m, B);	//AAt
	}	
}

template<typename _T> void Gen_Translation_Matrix_2D(_T P_1[2], _T P_2[2], _T T[3 * 3])
{//��P_1��ƽ�Ƶ�P_2�㣬�����ƽ�ƾ���T
	T[0] = T[4] = T[8] = 1;
	T[1] = T[3] = T[6] = T[7] = 0;
	T[2] = P_2[0] - P_1[0];
	T[5] = P_2[1] - P_2[1];
	return;
}
template<typename _T> void Gen_Rotation_Matrix_2D(_T theta, _T R[3 * 3])
{//��ԭ����ʱ����ת�ȽǶ�
	//x' = x*cos�� - y*sin�� => x'	= cos�� -sin��	*	x
	//y' = x*sin�� + y*cos��	y'	  sin��  cos��		y
	R[0] = R[4] = cos(theta);
	R[3] = sin(theta), R[1] = -R[3];
	R[2] = R[5] = R[6] = R[7] = 0;
	R[8] = 1;
	return;
}

template<typename _T> void Gen_Roation_Matrix_2D(_T Rotation_Center[2], _T theta, _T R[3 * 3])
{//��ĳ����ת�ȽǶȣ�������ξ���
	//��һ������Rotation_Center�Ƶ�ԭ��
	_T T_1[3 * 3];	// = { 1, 0, -Rotation_Center[0],
	//0, 1, -Rotation_Center[1],
	//0, 0, 1 };
	_T O[2] = { 0.f,0.f };
	Gen_Translation_Matrix_2D(Rotation_Center, O, T_1);

	//�ڶ�����������ת�任
	_T R_2[3 * 3];
	Gen_Rotation_Matrix_2D(theta, R_2);

	//���������ƻ�ԭ��λ��
	_T T_3[3 * 3];
	Gen_Translation_Matrix_2D(O, Rotation_Center, T_3);

	//�������� T_3 * R_2 * T_1
	Matrix_Multiply(T_3, 3, 3, R_2, 3, R);
	Matrix_Multiply(R, 3, 3, T_1, 3, R);
	Disp(R, 3, 3, "R");
}
template<typename _T>void Get_R_t(_T T[4 * 4], _T R[3 * 3], _T t[3])
{//��4x4 ��ξ����г�ȡR��t
	if (R)
	{
		R[0] = T[0], R[1] = T[1], R[2] = T[2];
		R[3] = T[4], R[4] = T[5], R[5] = T[6];
		R[6] = T[8], R[7] = T[9], R[8] = T[10];
	}
	if (t)
		t[0] = T[3], t[1] = T[7], t[2] = T[11];
	return;
}
template<typename _T>void Gen_Rotation_Matrix_2D(_T R[2 * 2], float fTheta)
{//����һ����ת����,��Բ����ʱ����תfTheta���ȣ�����������P, ��P'=RxP;
	//ԭ����������	x' = x*cos(��) - y*sin(��) =	cos(��) -sin(��) *	x
	//					y' = x*sin(��) + y*cos(��)		sin(��) cos(��)		y
	R[3] = R[0] = cos(fTheta);
	R[2] = sin(fTheta);
	R[1] = -R[2];
	return;
}
template<typename _T>void Gen_Homo_Matrix_2D(_T R[2 * 2], _T t[2], _T T[3 * 3])
{
	T[0] = R[0], T[1] = R[1], T[2] = t[0];
	T[3] = R[2], T[4] = R[3], T[5] = t[1];
	T[6] = T[7] = 0, T[8] = 1;
}

template<typename _T> _T fGet_Value(Sparse_Matrix<_T>* poMatrix, int x, int y)
{
	typename Sparse_Matrix<_T>::Item* poItem;
	Get_Item(*poMatrix, x, y, poItem);
	return poItem ? poItem->m_fValue : 0;
}

template<typename _T>void Set_Value(Sparse_Matrix<_T>* poMatrix, int x, int y, _T fValue)
{//��ϡ�����ֵĳ��Ԫ
	/*Sparse_Matrix<_T>::Item* poItem;
	Get_Item(oMatrix, x, y, poItem);
	if (poItem)
	{
		poItem->m_fValue = fValue;
		return;
	}
	if (oMatrix.m_iCur_Item >= oMatrix.m_iItem_Count)
	{
		printf("Exceed max Item _Count:%d in Set_Value\n", oMatrix.m_iItem_Count);
		return;
	}
	poItem =&oMatrix.m_pBuffer[oMatrix.m_iCur_Item++];
	poItem->x = x;
	poItem->y = y;
	poItem->m_fValue = fValue;*/

	Sparse_Matrix<_T> oMatrix = *poMatrix;
	//���뵽��ӦΪֹ
	typename Sparse_Matrix<_T>::Item* poCur, * poPrevious, oItem;
	int iCur_Link_Item;
	oItem.x = x;
	oItem.y = y;
	oItem.m_fValue = fValue;
	//���ҵ���Ӧ����
	if (!(iCur_Link_Item = oMatrix.m_pRow[oItem.y]))
	{	//Ͱ��û��ָ��ֱ�Ӳ���
		oMatrix.m_pRow[oItem.y] = oMatrix.m_iCur_Item;
		oItem.m_iRow_Next = 0;
	}
	else {
		poCur = &oMatrix.m_pBuffer[iCur_Link_Item];
		poPrevious = NULL;
		while (poCur->x < oItem.x && poCur->m_iRow_Next)
		{
			poPrevious = poCur;
			poCur = &oMatrix.m_pBuffer[poCur->m_iRow_Next];
		}
		if (poCur->x < oItem.x)
		{
			oItem.m_iRow_Next = poCur->m_iRow_Next;
			poCur->m_iRow_Next = oMatrix.m_iCur_Item;
		}
		else if (poCur->x == oItem.x)
		{//��λ���Ѿ���ֵ����д����
			poCur->m_fValue = fValue;
			return;
		}else
		{
			if (poPrevious)
			{//����previous��cur֮��
				oItem.m_iRow_Next = poPrevious->m_iRow_Next;
				poPrevious->m_iRow_Next = oMatrix.m_iCur_Item;
			}
			else
			{
				oItem.m_iRow_Next = iCur_Link_Item;
				oMatrix.m_pRow[oItem.y] = oMatrix.m_iCur_Item;
			}
		}
	}

	//��˳����뵽��������
	if (!(iCur_Link_Item = oMatrix.m_pCol[oItem.x]))
	{	//Ͱ��û��ָ��ֱ�Ӳ���
		oMatrix.m_pCol[oItem.x] = oMatrix.m_iCur_Item;
		oItem.m_iCol_Next = 0;
	}
	else {
		poCur = &oMatrix.m_pBuffer[iCur_Link_Item];
		poPrevious = NULL;
		while (poCur->y < oItem.y && poCur->m_iCol_Next)
		{
			poPrevious = poCur;
			poCur = &oMatrix.m_pBuffer[poCur->m_iCol_Next];
		}
		if (poCur->y < oItem.y)
		{
			oItem.m_iCol_Next = poCur->m_iCol_Next;
			poCur->m_iCol_Next = oMatrix.m_iCur_Item;
		}
		else
		{
			if (poPrevious)
			{//����previous��cur֮��
				oItem.m_iCol_Next = poPrevious->m_iCol_Next;
				poPrevious->m_iCol_Next = oMatrix.m_iCur_Item;
			}
			else
			{
				oItem.m_iCol_Next = iCur_Link_Item;
				oMatrix.m_pCol[oItem.x] = oMatrix.m_iCur_Item;
			}
		}
	}

	if (oItem.y >= (unsigned int)oMatrix.m_iRow_Count)
		oMatrix.m_iRow_Count = oItem.y + 1;
	if (oItem.x >= (unsigned int)oMatrix.m_iCol_Count)
		oMatrix.m_iCol_Count = oItem.x + 1;

	if (oMatrix.m_iCur_Item >= oMatrix.m_iMax_Item_Count)
	{
		printf("Exceed max item count in Set_Value:%d %d\n", x, y);
		return;
	}
	oMatrix.m_pBuffer[oMatrix.m_iCur_Item++] = oItem;
	*poMatrix = oMatrix;
	return;
}

template<typename _T>void Matrix_Multiply(Sparse_Matrix<_T> A, Sparse_Matrix<_T> B, Sparse_Matrix<_T>* poC,int bCompact)
{
	int y, x, i0_Count = 0, iNew_Item_Count = 0;
	_T fValue;
	typename Sparse_Matrix<_T>::Item* poRow_Cur, * poCol_Cur, * poNew_Row_Pre = NULL, oNew_Item;//* poNew_Col_Pre = NULL,
	Sparse_Matrix<_T> oC;
	//pNew_Col_Preʵ������һ��Item, ���Դ���ϴθ��е����һ��Item,�Ա�ӿ��ٶ�
	typename Sparse_Matrix<_T>::Item** pNew_Col_Pre = (typename Sparse_Matrix<_T>::Item**)malloc(B.m_iCol_Count * sizeof(typename Sparse_Matrix<_T>::Item*));

	if (A.m_iCol_Count != B.m_iRow_Count || !poC)
	{
		printf("Mismatch matrix size in Matrix_Multiply\n");
		return;
	}
		
	if (!poC->m_iMax_Item_Count)
		Init_Sparse_Matrix(&oC, A.m_iRow_Count * B.m_iCol_Count, B.m_iCol_Count, A.m_iRow_Count);
		//Init_Sparse_Matrix(&oC, A.m_iRow_Count * B.m_iCol_Count, Max(A.m_iRow_Count, B.m_iCol_Count));
	else
	{
		oC = *poC;
		memset(oC.m_pRow, 0, oC.m_iRow_Count * sizeof(int));
		memset(oC.m_pCol, 0, oC.m_iCol_Count * sizeof(int));
	}
	for (y = 0; y < A.m_iRow_Count; y++)
	{
		for (x = 0; x < B.m_iCol_Count; x++)
		{
			fValue = 0;
			if (A.m_pRow[y] && B.m_pCol[x])
			{
				poRow_Cur = &A.m_pBuffer[A.m_pRow[y]];
				poCol_Cur = &B.m_pBuffer[B.m_pCol[x]];
				//int iCount = 0;
				do
				{
					if (poRow_Cur->x == poCol_Cur->y)
					{
						fValue += poRow_Cur->m_fValue * poCol_Cur->m_fValue;
						poRow_Cur = poRow_Cur->m_iRow_Next ? &A.m_pBuffer[poRow_Cur->m_iRow_Next] : NULL;
						poCol_Cur = poCol_Cur->m_iCol_Next ? &B.m_pBuffer[poCol_Cur->m_iCol_Next] : NULL;
					}else if (poRow_Cur->x < poCol_Cur->y)
						poRow_Cur = poRow_Cur->m_iRow_Next ? &A.m_pBuffer[poRow_Cur->m_iRow_Next] : NULL;
					else
						poCol_Cur = poCol_Cur->m_iCol_Next ? &B.m_pBuffer[poCol_Cur->m_iCol_Next] : NULL;
					//iCount++;
				} while (poRow_Cur && poCol_Cur);
				/*if (iCount > 3)
					printf("%d ",iCount);*/
			}

			if (fValue != 0)
			{//��Ŀ�����c�м�һ��Item
				oNew_Item.x = x;
				oNew_Item.y = y;
				oNew_Item.m_fValue = fValue;
				oNew_Item.m_iRow_Next = oNew_Item.m_iCol_Next = 0;
				oC.m_pBuffer[++iNew_Item_Count] = oNew_Item;

				//�ȼ���������
				if (!oC.m_pRow[y])
					oC.m_pRow[y] = iNew_Item_Count;
				else
					poNew_Row_Pre->m_iRow_Next = iNew_Item_Count;

				//�ټ���������
				if (!oC.m_pCol[x])
					oC.m_pCol[x] = iNew_Item_Count;
				else
					pNew_Col_Pre[x]->m_iCol_Next = iNew_Item_Count;
				pNew_Col_Pre[x] = poNew_Row_Pre = &oC.m_pBuffer[iNew_Item_Count];
			}else
				i0_Count++;
		}
	}
	oC.m_iRow_Count = A.m_iRow_Count;
	oC.m_iCol_Count = B.m_iCol_Count;

	oC.m_iCur_Item = iNew_Item_Count;
	if (bCompact)
		Compact_Sparse_Matrix(&oC);
	*poC = oC;
	free(pNew_Col_Pre);
	return;
}

template<typename _T>void Matrix_Multiply(Sparse_Matrix<_T> A, _T a)
{//C = aA
	int i;
	for (i = 1; i <= A.m_iCur_Item; i++)
		A.m_pBuffer[i].m_fValue *= a;
	return;
}
template<typename _T>void Matrix_Transpose_1(Sparse_Matrix<_T> A, Sparse_Matrix<_T>* poAt)
{
	Sparse_Matrix<_T> At = *poAt;
	typename Sparse_Matrix<_T>::Item* poCur;
	unsigned int iTemp;
	//Init_Sparse_Matrix(&At, A.m_iItem_Count, Max(A.m_iRow_Count, A.m_iCol_Count));
	if (!At.m_pBuffer || At.m_iRow_Count != A.m_iCol_Count || At.m_iCol_Count != A.m_iRow_Count)
	{
		printf("Invalid parameter in Matrix_Transpose_1\n");
		return;
	}
	At.m_iRow_Count = A.m_iCol_Count;
	At.m_iCol_Count = A.m_iRow_Count;

	memcpy(At.m_pBuffer + 1, A.m_pBuffer + 1, A.m_iMax_Item_Count * sizeof(struct Sparse_Matrix<_T>::Item));
	memcpy(At.m_pRow, A.m_pCol, At.m_iRow_Count * sizeof(unsigned int));
	memcpy(At.m_pCol, A.m_pRow, At.m_iCol_Count * sizeof(unsigned int));

	for (int i = 1; i <= At.m_iMax_Item_Count; i++)
	{
		poCur = &At.m_pBuffer[i];
		//������ָ�룬��ָ��
		iTemp = poCur->m_iRow_Next;
		poCur->m_iRow_Next = poCur->m_iCol_Next;
		poCur->m_iCol_Next = iTemp;

		//�����кţ��к�
		iTemp = poCur->x;
		poCur->x = poCur->y;
		poCur->y = iTemp;
	}
	At.m_iCur_Item = A.m_iCur_Item;
	*poAt = At;
	return;
}

//template<typename _T>void Init_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix, int iItem_Count, int iMax_Order)
//{
//	poMatrix->m_iMax_Item_Count = iItem_Count;
//	int iSize = poMatrix->m_iMax_Item_Count * sizeof(typename Sparse_Matrix<_T>::Item);
//	unsigned char* pBuffer = (unsigned char*)pMalloc(&oMem_Mgr, iSize + iMax_Order * 2 * sizeof(unsigned int));
//	poMatrix->m_pBuffer = (typename Sparse_Matrix<_T>::Item*)pBuffer;
//	if (!pBuffer)
//		return;
//	poMatrix->m_pBuffer--;	//�Ա��[1]��ʼ
//	poMatrix->m_pRow = (unsigned int*)(pBuffer + iSize);
//	memset(poMatrix->m_pRow, 0, iMax_Order * 2 * sizeof(unsigned int));
//	poMatrix->m_pCol = poMatrix->m_pRow + iMax_Order;
//	poMatrix->m_iRow_Count = poMatrix->m_iCol_Count = 0;
//}

template<typename _T>void Init_Sparse_Matrix(Sparse_Matrix<_T>* poMatrix, int iItem_Count, int w, int h)
{
	poMatrix->m_iMax_Item_Count = iItem_Count;
	unsigned int iSize = poMatrix->m_iMax_Item_Count * sizeof(struct Sparse_Matrix<_T>::Item);
	unsigned char* pBuffer = (unsigned char*)pMalloc(&oMem_Mgr, iSize + (w + h) * sizeof(unsigned int));
	if (!pBuffer)
	{
		printf("Fail to allocate mem in Init_Sparse_Matrix\n");
		return;
	}
	poMatrix->m_pRow = (unsigned int*)pBuffer;
	poMatrix->m_pCol = poMatrix->m_pRow + h;
	memset(poMatrix->m_pRow, 0, (w + h) * sizeof(unsigned int));
	pBuffer += (w + h) * sizeof(unsigned int);
	poMatrix->m_pBuffer = (typename Sparse_Matrix<_T>::Item*)pBuffer;
	poMatrix->m_pBuffer--;	//�Ա��[1]��ʼ

	poMatrix->m_iRow_Count = h;
	poMatrix->m_iCol_Count = w;
	poMatrix->m_iCur_Item = 1;
}
template<typename _T>void Re_Arrange_Sparse_Matrix(Sparse_Matrix<_T>* poA)
{//��������һ��Sparse_Matrix,ɾȥ0Ԫ��
	Sparse_Matrix<_T> oA = *poA, oB;
	int iMax_Item_Count = oA.m_iMax_Item_Count;
	//Disp_Mem(&oMem_Mgr, 0);
	Compact_Sparse_Matrix(&oA);
	//Disp_Mem(&oMem_Mgr, 0);
	Init_Sparse_Matrix(&oB, iMax_Item_Count, oA.m_iCol_Count, oA.m_iRow_Count);
	if (oB.m_pRow)
	{
		Copy_Sparse_Matrix(oA, &oB);
		Free_Sparse_Matrix(&oA);
	}
	*poA = oB;
	return;
}
template<typename _T>void Matrix_Minus(Sparse_Matrix<_T> oA, Sparse_Matrix<_T> oB, Sparse_Matrix<_T>* poC, int bCompact)
{//�˴���������ٶ���̫����������걸�� C= A+B, 1�� C�ٶ��Ѿ��㹻�ռ䡣������ʧ���˳�
	//2��A�����B����ͬ�ĳ���
	Sparse_Matrix<_T> oC;
	typename Sparse_Matrix<_T>::Item oNew_Item;

	if (oA.m_iRow_Count != oB.m_iRow_Count || oA.m_iCol_Count != oB.m_iCol_Count)
	{
		printf("Mis matched size in Matrix_Add\n");
		return;
	}
	Init_Sparse_Matrix(&oC, oA.m_iCur_Item + oB.m_iCur_Item, oA.m_iCol_Count, oA.m_iRow_Count);
	oC.m_iCur_Item = 1;
	int y;
	typename Sparse_Matrix<_T>::Item oB_Cur, oA_Cur;
	for (y = 0; y < oC.m_iRow_Count; y++)
	{
		if (oC.m_iCur_Item >= oC.m_iMax_Item_Count && (oA.m_pRow[y] || oB.m_pRow[y]))
		{
			printf("Exceed Item count in Matrix_Add\n");
			return;
		}
		if (oA.m_pRow[y] && oB.m_pRow[y])
		{
			oA_Cur = oA.m_pBuffer[oA.m_pRow[y]];
			oB_Cur = oB.m_pBuffer[oB.m_pRow[y]];
			oC.m_pRow[y] = oC.m_iCur_Item;
			while (1)
			{
				if (oA_Cur.x == oB_Cur.x)
				{//���Լ�
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue - oB_Cur.m_fValue, oC.m_iCur_Item + 1 };
					if (oA_Cur.m_iRow_Next)
						oA_Cur = oA.m_pBuffer[oA_Cur.m_iRow_Next];
					else
						oA_Cur.x = 0xFFFFFFFF;
					if (oB_Cur.m_iRow_Next)
						oB_Cur = oB.m_pBuffer[oB_Cur.m_iRow_Next];
					else
						oB_Cur.x = 0xFFFFFFFF;
				}
				else if (oA_Cur.x < oB_Cur.x)
				{//д��AԪ��
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue, oC.m_iCur_Item + 1 };
					if (oA_Cur.m_iRow_Next)
						oA_Cur = oA.m_pBuffer[oA_Cur.m_iRow_Next];
					else
						oA_Cur.x = 0xFFFFFFFF;
				}
				else if (oB_Cur.x < oA_Cur.x)
				{
					oNew_Item = { oB_Cur.x,oB_Cur.y,-oB_Cur.m_fValue, oC.m_iCur_Item + 1 };
					if (oB_Cur.m_iRow_Next)
						oB_Cur = oB.m_pBuffer[oB_Cur.m_iRow_Next];
					else
						oB_Cur.x = 0xFFFFFFFF;
				}
				oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
				if (oA_Cur.x == 0xFFFFFFFF && oB_Cur.x == 0xFFFFFFFF)
				{
					oC.m_pBuffer[oC.m_iCur_Item - 1].m_iRow_Next = 0;
					break;
				}
			}
		}else if (!oA.m_pRow[y] && oB.m_pRow[y])
		{//���������B���м��뵽C
			oB_Cur = oB.m_pBuffer[oB.m_pRow[y]];
			oC.m_pRow[y] = oC.m_iCur_Item;
			while (1)
			{
				if (oB_Cur.m_iRow_Next)
				{//������һ��Ԫ��
					oNew_Item = { oB_Cur.x,oB_Cur.y,-oB_Cur.m_fValue, oC.m_iCur_Item + 1 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					oB_Cur = oB.m_pBuffer[oB_Cur.m_iRow_Next];
				}else
				{
					oNew_Item = { oB_Cur.x,oB_Cur.y,-oB_Cur.m_fValue, 0 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					break;
				}
			}
		}else if (!oB.m_pRow[y] && oA.m_pRow[y])
		{//A���м���C
			oA_Cur = oA.m_pBuffer[oA.m_pRow[y]];
			oC.m_pRow[y] = oC.m_iCur_Item;
			while (1) {
				if (oA_Cur.m_iRow_Next)
				{//������һ��Ԫ��
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue, oC.m_iCur_Item + 1 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					oA_Cur = oA.m_pBuffer[oA_Cur.m_iRow_Next];
				}
				else
				{
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue, 0 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					break;
				}
			}
		}
	}
	Build_Link_Col(oC);
	if (poC->m_iMax_Item_Count)
		Free_Sparse_Matrix(poC);
	if (bCompact)
		Compact_Sparse_Matrix(&oC);
	*poC = oC;
	return;
}
template<typename _T>void Matrix_Add(Sparse_Matrix<_T> oA, Sparse_Matrix<_T> oB, Sparse_Matrix<_T>* poC,int bCompact)
{//�˴���������ٶ���̫����������걸�� C= A+B, 1�� C�ٶ��Ѿ��㹻�ռ䡣������ʧ���˳�
	//2��A�����B����ͬ�ĳ���
	Sparse_Matrix<_T> oC;
	typename Sparse_Matrix<_T>::Item oNew_Item;

	if (oA.m_iRow_Count != oB.m_iRow_Count || oA.m_iCol_Count != oB.m_iCol_Count)
	{
		printf("Mis matched size in Matrix_Add\n");
		return;
	}
	Init_Sparse_Matrix(&oC, oA.m_iCur_Item + oB.m_iCur_Item, oA.m_iCol_Count, oA.m_iRow_Count);
	oC.m_iCur_Item = 1;
	int y;
	typename Sparse_Matrix<_T>::Item oB_Cur, oA_Cur;
	for (y = 0; y < oC.m_iRow_Count; y++)
	{
		if (oC.m_iCur_Item >= oC.m_iMax_Item_Count && (oA.m_pRow[y] || oB.m_pRow[y]))
		{
			printf("Exceed Item count in Matrix_Add\n");
			return;
		}
		if (oA.m_pRow[y] && oB.m_pRow[y])
		{
			oA_Cur = oA.m_pBuffer[oA.m_pRow[y]];
			oB_Cur = oB.m_pBuffer[oB.m_pRow[y]];
			oC.m_pRow[y] = oC.m_iCur_Item;
			while (1)
			{
				if (oA_Cur.x == oB_Cur.x)
				{//���Լ�
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue + oB_Cur.m_fValue, oC.m_iCur_Item + 1 };
					if (oA_Cur.m_iRow_Next)
						oA_Cur = oA.m_pBuffer[oA_Cur.m_iRow_Next];
					else
						oA_Cur.x = 0xFFFFFFFF;
					if (oB_Cur.m_iRow_Next)
						oB_Cur = oB.m_pBuffer[oB_Cur.m_iRow_Next];
					else
						oB_Cur.x = 0xFFFFFFFF;
				}
				else if (oA_Cur.x < oB_Cur.x)
				{//д��AԪ��
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue, oC.m_iCur_Item + 1 };
					if (oA_Cur.m_iRow_Next)
						oA_Cur = oA.m_pBuffer[oA_Cur.m_iRow_Next];
					else
						oA_Cur.x = 0xFFFFFFFF;
				}
				else if (oB_Cur.x < oA_Cur.x)
				{
					oNew_Item = { oB_Cur.x,oB_Cur.y,oB_Cur.m_fValue, oC.m_iCur_Item + 1 };
					if (oB_Cur.m_iRow_Next)
						oB_Cur = oB.m_pBuffer[oB_Cur.m_iRow_Next];
					else
						oB_Cur.x = 0xFFFFFFFF;
				}
				oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
				if (oA_Cur.x == 0xFFFFFFFF && oB_Cur.x == 0xFFFFFFFF)
				{
					oC.m_pBuffer[oC.m_iCur_Item - 1].m_iRow_Next = 0;
					break;
				}
			}
		}
		else if (!oA.m_pRow[y] && oB.m_pRow[y])
		{//���������B���м��뵽C
			oB_Cur = oB.m_pBuffer[oB.m_pRow[y]];
			oC.m_pRow[y] = oC.m_iCur_Item;
			while (1)
			{
				if (oB_Cur.m_iRow_Next)
				{//������һ��Ԫ��
					oNew_Item = { oB_Cur.x,oB_Cur.y,oB_Cur.m_fValue, oC.m_iCur_Item + 1 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					oB_Cur = oB.m_pBuffer[oB_Cur.m_iRow_Next];
				}
				else
				{
					oNew_Item = { oB_Cur.x,oB_Cur.y,oB_Cur.m_fValue, 0 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					break;
				}
			}
		}
		else if (!oB.m_pRow[y] && oA.m_pRow[y])
		{//A���м���C
			oA_Cur = oA.m_pBuffer[oA.m_pRow[y]];
			oC.m_pRow[y] = oC.m_iCur_Item;
			while (1) {
				if (oA_Cur.m_iRow_Next)
				{//������һ��Ԫ��
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue, oC.m_iCur_Item + 1 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					oA_Cur = oA.m_pBuffer[oA_Cur.m_iRow_Next];
				}else
				{
					oNew_Item = { oA_Cur.x,oA_Cur.y,oA_Cur.m_fValue, 0 };
					oC.m_pBuffer[oC.m_iCur_Item++] = oNew_Item;
					break;
				}
			}
		}
	}
	Build_Link_Col(oC);
	if (poC->m_iMax_Item_Count)
		Free_Sparse_Matrix(poC);
	if (bCompact)
		Compact_Sparse_Matrix(&oC);
	*poC = oC;
	return;
}
template<typename _T>void Copy_Sparse_Matrix(Sparse_Matrix<_T> oA, Sparse_Matrix<_T>* poB)
{
	Sparse_Matrix<_T> oB = *poB;
	typename Sparse_Matrix<_T>::Item oOrg;
	int y;
	for (y = 0; y < oA.m_iRow_Count; y++)
	{
		if (!oA.m_pRow[y])
			continue;
		oOrg = oA.m_pBuffer[oA.m_pRow[y]];
		while (1)
		{
			if (oOrg.m_fValue != 0)
			{
				if (!oB.m_pRow[y])
					oB.m_pRow[y] = oB.m_iCur_Item;
				oB.m_pBuffer[oB.m_iCur_Item++] = { oOrg.x, oOrg.y,oOrg.m_fValue,oB.m_iCur_Item + 1 };
			}
			if (oOrg.m_iRow_Next)
				oOrg = oA.m_pBuffer[oOrg.m_iRow_Next];
			else
				break;
		}
		if (oB.m_pRow[y])
			oB.m_pBuffer[oB.m_iCur_Item - 1].m_iRow_Next = 0;
	}
	//Build_Link_Col(oB);
	Build_Link_Col_1(oB);

	*poB = oB;
}
template<typename _T>void Add_Col_Link(Sparse_Matrix<_T> oA, typename Sparse_Matrix<_T>::Item** pCol_All, typename Sparse_Matrix<_T>::Item* poItem)
{//����һ��Col�����һ��Item��������������
	typename Sparse_Matrix<_T>::Item *poCur = pCol_All[poItem->x],*poPre=NULL;
	//if (poItem->x == 6)
		//printf("y:%d\n", poCur->y);
	while (poCur->y <=poItem->y)
	{
		poPre = poCur;
		if (poCur->m_iCol_Next)
			poCur = &oA.m_pBuffer[poCur->m_iCol_Next];
		else
			break;
	}

	int iItem_Index;
	if (poPre)
	{
		iItem_Index =(int)(poItem - oA.m_pBuffer);
		poItem->m_iCol_Next = poPre->m_iCol_Next;
		poPre->m_iCol_Next = iItem_Index;
	}else
	{
		iItem_Index =(int)(pCol_All[poItem->x] - oA.m_pBuffer);
		poItem->m_iCol_Next = iItem_Index;
		//�˴���Ҫ����oA.m_pCol
		oA.m_pCol[poItem->x] =(int)(poItem - oA.m_pBuffer);
	}
	pCol_All[poItem->x] = poItem;
	return;
}
template<typename _T>void Test_Linear(_T A[], int iOrder, _T B[], _T X[],_T *pfError_Sum)
{//����Ax=b
	_T* B1 = (_T*)pMalloc(&oMem_Mgr, iOrder * sizeof(_T));
	_T fError_Sum;
	Matrix_Multiply(A, iOrder, iOrder, X, 1, B1);
	if (B)
		fError_Sum = fGet_Distance(B1, B, iOrder);
	else
		fError_Sum = fGet_Mod(B1, iOrder);
	//printf("Error sum:%f\n", fError_Sum);
	if (B1)
		Free(B1);
	if (pfError_Sum)
		*pfError_Sum += fError_Sum;
	return;
}
template<typename _T>void Test_Linear(Sparse_Matrix<_T> oA, _T B[],_T X[])
{
	//�˴��������㿴��
	Sparse_Matrix<_T> oB, oX;
	_T* B1 = (_T*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(_T));

	Init_Sparse_Matrix(&oB, oA.m_iRow_Count, 1, oA.m_iRow_Count);
	Init_Sparse_Matrix(&oX, oA.m_iRow_Count, 1, oA.m_iRow_Count);
	//Dense_2_Sparse(B, oA2.m_iRow_Count, 1, &oB);
	Dense_2_Sparse(X, oA.m_iRow_Count, 1, &oX);
	Matrix_Multiply(oA, oX, &oB);
	Sparse_2_Dense(oB, B1);
	printf("Error Sum:%f\n", fGet_Distance(B1, B, oB.m_iRow_Count));
	Free(&oMem_Mgr, B1);
	Free_Sparse_Matrix(&oB);
	Free_Sparse_Matrix(&oX);

	/*_T* A = (_T*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * oA.m_iCol_Count * sizeof(_T));
	double *A1= (double*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * oA.m_iCol_Count * sizeof(double));
	double* X1 = (double*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(double));
	double *B1= (double*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(double));
	int iResult;
	Sparse_2_Dense(oA, A);
	for (int i = 0; i < oA.m_iRow_Count * oA.m_iRow_Count; i++)
		A1[i] = A[i];
	for (int i = 0; i < oA.m_iRow_Count; i++)
		B1[i] = B[i];
	Solve_Linear_Gause(A1, oA.m_iRow_Count, B1, X1,&iResult);*/

	/*_T* A = (_T*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * oA.m_iCol_Count * sizeof(_T));
	_T* X1 = (_T*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(_T));
	int iResult;
	Sparse_2_Dense(oA, A);
	Solve_Linear_Gause(A, oA.m_iRow_Count, B, X1,&iResult);
	Free(&oMem_Mgr, A);
	Free(&oMem_Mgr, X1);*/
	return;
}

#pragma pack(4)
template<typename _T>struct Solve_Linear_Gause_Row_Item { //����������Ԫ�������Է���
	_T m_fValue;
	unsigned int x;
	union {
		unsigned int m_iNext;			//��ʵNode����һ��λ��
		unsigned int m_iFree_Pos;		//��Nodeָ�����λ��
	};	
};
template<typename _T>struct Solve_Linear_Gause_Index {       //����row, col index���㸴��Щ
	_T m_fHead_Value;
	unsigned long long m_iHead_Index : 24;			//һ�л�һ�еĶ�ͷ��Item Index
	unsigned long long m_iRow_Col_No : 20;          //ԭ�����к�
	unsigned long long m_iHead_Coord:20;				//������Index, �����ǿ�ʼԪ�ص�x �ͽ�Ԫ�ص�x
};
#pragma pack()

template<typename _T>void Solve_Linear_Gause_3(Sparse_Matrix<_T> oA, _T B[], _T X[], int* pbSuccess)
{//��С����̫��������ƭһƪר����ûɶ���ü�ֵ
	int x, y, bSuccess = 1;
	unsigned int iCur_Item, iFree_Item_End, iMax_Item_Count;
	Solve_Linear_Gause_Index<_T>* pRow;
	union {
		Solve_Linear_Gause_Index<_T> oMax_Row_Index;
		Solve_Linear_Gause_Index<_T> oRow_Index;
	};

	Solve_Linear_Gause_Row_Item<_T>* pBuffer,
		* pMax_Row, * poMax_Row_Item,
		* poItem_To_Elim;
	union {
		Solve_Linear_Gause_Row_Item<_T>* poRow_Item;
		Solve_Linear_Gause_Row_Item<_T>* poPre;
		unsigned int iSize;
	};
	pRow = (Solve_Linear_Gause_Index<_T>*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(Solve_Linear_Gause_Index<_T>));
	iMax_Item_Count = oA.m_iCur_Item * 20;
	iSize = iMax_Item_Count * sizeof(Solve_Linear_Gause_Row_Item<_T>);
	pBuffer = (Solve_Linear_Gause_Row_Item<_T>*)pMalloc(&oMem_Mgr, iSize);
	memset(pBuffer, 0, iSize);

	pBuffer--;
	typename Sparse_Matrix<_T>::Item* poItem;

	//��Aȫ������
	iCur_Item = 1;		//�˴�ָ����õ�λ��
	for (y = 0; y < oA.m_iRow_Count; y++)
	{
		oRow_Index.m_iHead_Index = iCur_Item;
		oRow_Index.m_iRow_No = y;
		if (!oA.m_pRow[y])
			oRow_Index.m_iHead_Coord = oA.m_iCol_Count;		//����ֻ�ܼ���B
		else
		{
			poItem = &oA.m_pBuffer[oA.m_pRow[y]];
			oRow_Index.m_iHead_Coord = poItem->x;
			oRow_Index.m_fHead_Value = poItem->m_fValue;
			while (1)
			{
				pBuffer[iCur_Item++] = { poItem->m_fValue, poItem->x,iCur_Item + 1 };
				if (poItem->m_iRow_Next)
					poItem = &oA.m_pBuffer[poItem->m_iRow_Next];
				else
					break;
			}
		}
		//�ټӶ�һ��BԪ��
		if (B[y])
			pBuffer[iCur_Item++] = { B[y], (unsigned int)oA.m_iCol_Count,0 };
		else
			pBuffer[iCur_Item - 1].m_iNext = 0;

		pRow[y] = oRow_Index;
	}
	iFree_Item_End = iCur_Item;

	_T fMax, fValue;
	const _T eps = (_T)1e-10;
	int iMax_Row;

	for (x = 0; x < oA.m_iCol_Count; x++)
	{
		//���ҵ�x�е�����Ԫ
		fMax = 0;
		for (y = x; y < oA.m_iRow_Count; y++)
		{
			oMax_Row_Index = pRow[y];
			if (oMax_Row_Index.m_iHead_Coord == x && Abs(oMax_Row_Index.m_fHead_Value) > fMax)
			{//ֻ����һ�еĿ�ʼλ��Ϊx���ж�
				fMax = Abs(oMax_Row_Index.m_fHead_Value);
				iMax_Row = y;
			}
		}
		if (fMax < eps)
		{
			printf("������ in Solve_Linear_Gause_2\n");
			bSuccess = 0;
			goto END;
		}
		fMax = pRow[iMax_Row].m_fHead_Value;

		//������Ԫ����ֵ
		oMax_Row_Index = pRow[iMax_Row];
		poMax_Row_Item = &pBuffer[oMax_Row_Index.m_iHead_Index];
		poMax_Row_Item->m_fValue = 1.f;
		if (poMax_Row_Item->m_iNext)
		{
			poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
			pMax_Row = poMax_Row_Item;
			while (1)
			{
				poMax_Row_Item->m_fValue /= oMax_Row_Index.m_fHead_Value;
				if (poMax_Row_Item->m_iNext)
					poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
				else
					break;
			}
			pRow[iMax_Row] = { pMax_Row->m_fValue, (unsigned int)(pMax_Row - pBuffer),oMax_Row_Index.m_iRow_No,(unsigned int)(pMax_Row->x) };
		}else
		{
			pRow[iMax_Row] = { poMax_Row_Item->m_fValue,0,oMax_Row_Index.m_iRow_No,0 };
			continue;   //B[x]=0 ������Ԫ��
		}
		y = x;

		//������Ԫ�е���Ԫ�ع黹
		if (iFree_Item_End >= iMax_Item_Count)
		{
			printf("Insufficient space of oA1 in Solve_Linear_Gause_1\n");
			return;
		}
		pBuffer[iFree_Item_End++].m_iFree_Pos = oMax_Row_Index.m_iHead_Index;

		//���л���ȥ
		swap(pRow[iMax_Row], pRow[y]);

		for (y = x + 1; y < oA.m_iRow_Count; y++)
		{//������Ԫѭ��

			Solve_Linear_Gause_Index<_T> oRow_To_Elim_Index = pRow[y];
			if (oRow_To_Elim_Index.m_iHead_Coord > (unsigned int)x)
				continue;
			poMax_Row_Item = pMax_Row;
			poPre = poItem_To_Elim = &pBuffer[oRow_To_Elim_Index.m_iHead_Index];
			fValue = poItem_To_Elim->m_fValue;
			poItem_To_Elim = poItem_To_Elim->m_iNext ? &pBuffer[poItem_To_Elim->m_iNext] : NULL;
			
			while (1)
			{
				if (iCur_Item >= iMax_Item_Count)
				{
					printf("Insufficient space of oA1 in Solve_Linear_Gause_1\n");
					bSuccess = 0;
					goto END;
				}
				//�ָ������
				if (!poItem_To_Elim)
				{//�ӽ�ȥ���˴������ò��ϣ���Ϊ��B�еĴ���
					pBuffer[iCur_Item] = { -fValue * poMax_Row_Item->m_fValue,poMax_Row_Item->x };
					poPre->m_iNext = iCur_Item++;
					if (!poMax_Row_Item->m_iNext)
						break;
					poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
					poPre = &pBuffer[poPre->m_iNext];
				}else if (poItem_To_Elim->x < poMax_Row_Item->x)
				{//������� ����Ԫ��Ԫ�����ҽ�һ��                 
					poPre = poItem_To_Elim;
					poItem_To_Elim = poItem_To_Elim->m_iNext ? &pBuffer[poItem_To_Elim->m_iNext] : NULL;
				}else if (poMax_Row_Item->x < poItem_To_Elim->x)
				{//��1��Ԫ�أ�����Ԫ��Ԫ��������һ��
					Solve_Linear_Gause_Row_Item<_T>* pNew_Item = &pBuffer[iCur_Item];
					if (pNew_Item->m_iFree_Pos)
					{
						int iFree_Pos = pNew_Item->m_iFree_Pos;
						pNew_Item = &pBuffer[iFree_Pos];
						*pNew_Item = { -fValue * poMax_Row_Item->m_fValue,poMax_Row_Item->x, poPre->m_iNext };
						poPre->m_iNext = iFree_Pos;
						if (!poMax_Row_Item->m_iNext)
							break;
						poPre = pNew_Item;
						poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];

						iFree_Item_End--;
						pBuffer[iCur_Item].m_iFree_Pos = pBuffer[iFree_Item_End].m_iFree_Pos;
						pBuffer[iFree_Item_End].m_iFree_Pos = 0;
					}else
					{
						/*if (iCur_Item == 806329)
							printf("Here");*/
						*pNew_Item = { -fValue * poMax_Row_Item->m_fValue,poMax_Row_Item->x, poPre->m_iNext };
						poPre->m_iNext = iCur_Item++;
						if (!poMax_Row_Item->m_iNext)
							break;
						poPre = &pBuffer[iCur_Item - 1];
						poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
					}					
				}else
				{//��ȣ�����ֵ��Ȼ����������ͬʱ�����ƶ�һ��
					poItem_To_Elim->m_fValue -= fValue * poMax_Row_Item->m_fValue;
					if (!poMax_Row_Item->m_iNext)
						break;
					poPre = poItem_To_Elim;
					poItem_To_Elim = poItem_To_Elim->m_iNext ? &pBuffer[poItem_To_Elim->m_iNext] : NULL;
					poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
				}
			}
			poItem_To_Elim = &pBuffer[oRow_To_Elim_Index.m_iHead_Index];
			//��һ�н�����ʱ�򣬰����е����ݹ黹����pBuffer�����
			if (iCur_Item > iFree_Item_End)
				iFree_Item_End = iCur_Item;
			pBuffer[iFree_Item_End++].m_iFree_Pos = oRow_To_Elim_Index.m_iHead_Index;
			
			if (poItem_To_Elim->m_iNext)
			{
				oRow_To_Elim_Index.m_iHead_Index = poItem_To_Elim->m_iNext;
				poItem_To_Elim = &pBuffer[oRow_To_Elim_Index.m_iHead_Index];
				oRow_To_Elim_Index.m_iHead_Coord = poItem_To_Elim->x;
				oRow_To_Elim_Index.m_fHead_Value = poItem_To_Elim->m_fValue;
				pRow[y] = oRow_To_Elim_Index;
			}else
			{
				printf("������ in Solve_Linear_Gause_2\n");
				bSuccess = 0;
				goto END;
			}
		}
	}
	//����������ʽ�ش�
	for (y = oA.m_iRow_Count - 1; y >= 0; y--)
	{
		oRow_Index = pRow[y];
		if (!oRow_Index.m_iHead_Index)
		{
			X[y] = 0;
			continue;
		}
		poRow_Item = &pBuffer[oRow_Index.m_iHead_Index];
		fValue = 0;

		while (1)
		{
			if (poRow_Item->x < (unsigned int)oA.m_iCol_Count)
			{
				fValue += poRow_Item->m_fValue * X[poRow_Item->x];
				if (poRow_Item->m_iNext)
					poRow_Item = &pBuffer[poRow_Item->m_iNext];
				else
					break;
			}
			else
				break;
		}
		if (poRow_Item->x == oA.m_iCol_Count)
			X[y] = poRow_Item->m_fValue - fValue;
		else
			X[y] = -fValue;
	}

END:
	if (pBuffer)
		Free(&oMem_Mgr, pBuffer + 1);
	if (pRow)
		Free(&oMem_Mgr, pRow);
	if (pbSuccess)
		*pbSuccess = bSuccess;
	return;
}
template<typename _T>void Solve_Linear_Gause_2(Sparse_Matrix<_T> oA, _T B[], _T X[], int* pbSuccess)
{//�����С���Ժ������������_1�汾
	int x, y, bSuccess = 1;
	unsigned int iCur_Item, iMax_Item_Count;

	Solve_Linear_Gause_Index<_T>* pRow;
	union {
		Solve_Linear_Gause_Index<_T> oMax_Row_Index;
		Solve_Linear_Gause_Index<_T> oRow_Index;
	};

	Solve_Linear_Gause_Row_Item<_T>* pBuffer,
		* pMax_Row, * poMax_Row_Item,
		* poItem_To_Elim;
	union {
		Solve_Linear_Gause_Row_Item<_T>* poRow_Item;
		Solve_Linear_Gause_Row_Item<_T>* poPre;
		unsigned int iSize;
	};
	pRow = (Solve_Linear_Gause_Index<_T>*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(Solve_Linear_Gause_Index<_T>));
	iMax_Item_Count = oA.m_iCur_Item * 20;
	iSize = iMax_Item_Count * sizeof(Solve_Linear_Gause_Row_Item<_T>);
	pBuffer = (Solve_Linear_Gause_Row_Item<_T>*)pMalloc(&oMem_Mgr, iSize);
	
	pBuffer--;
	typename Sparse_Matrix<_T>::Item* poItem;

	//��Aȫ������,��������
	iCur_Item = 1;
	for (y = 0; y < oA.m_iRow_Count; y++)
	{
		oRow_Index.m_iHead_Index = iCur_Item;
		oRow_Index.m_iRow_Col_No = y;
		if (!oA.m_pRow[y])
			oRow_Index.m_iHead_Coord = oA.m_iCol_Count;		//����ֻ�ܼ���B
		else
		{
			poItem = &oA.m_pBuffer[oA.m_pRow[y]];
			oRow_Index.m_iHead_Coord = poItem->x;
			oRow_Index.m_fHead_Value = poItem->m_fValue;
			while (1)
			{
				pBuffer[iCur_Item++] = { poItem->m_fValue, poItem->x,iCur_Item + 1 };
				if (poItem->m_iRow_Next)
					poItem = &oA.m_pBuffer[poItem->m_iRow_Next];
				else
					break;
			}
		}
		//�ټӶ�һ��BԪ��
		if (B[y])
			pBuffer[iCur_Item++] = { B[y], (unsigned int)oA.m_iCol_Count,0 };
		else
			pBuffer[iCur_Item - 1].m_iNext = 0;

		pRow[y] = oRow_Index;
	}

	_T fMax, fValue;
	const _T eps = (_T)1e-10;
	int iMax_Row;

	for (x = 0; x < oA.m_iCol_Count; x++)
	{
		//���ҵ�x�е�����Ԫ
		fMax = 0;
		for (y = x; y < oA.m_iRow_Count; y++)
		{
			oMax_Row_Index = pRow[y];
			if (oMax_Row_Index.m_iHead_Coord == x && Abs(oMax_Row_Index.m_fHead_Value) > fMax)
			{//ֻ����һ�еĿ�ʼλ��Ϊx���ж�
				fMax = Abs(oMax_Row_Index.m_fHead_Value);
				iMax_Row = y;
			}
		}
		if (fMax < eps)
		{
			printf("������ in Solve_Linear_Gause_2\n");
			bSuccess = 0;
			goto END;
		}
		fMax = pRow[iMax_Row].m_fHead_Value;

		//������Ԫ����ֵ
		oMax_Row_Index = pRow[iMax_Row];
		poMax_Row_Item = &pBuffer[oMax_Row_Index.m_iHead_Index];
		poMax_Row_Item->m_fValue = 1.f;
		if (poMax_Row_Item->m_iNext)
		{
			poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
			pMax_Row = poMax_Row_Item;
			while (1)
			{
				poMax_Row_Item->m_fValue /= oMax_Row_Index.m_fHead_Value;
				if (poMax_Row_Item->m_iNext)
					poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
				else
					break;
			}
			pRow[iMax_Row] = { pMax_Row->m_fValue, (unsigned int)(pMax_Row - pBuffer),oMax_Row_Index.m_iRow_Col_No,(unsigned int)(pMax_Row->x) };
		}else
		{
			pRow[iMax_Row] = { poMax_Row_Item->m_fValue,0,oMax_Row_Index.m_iRow_Col_No,0 };
			continue;   //B[x]=0 ������Ԫ��
		}
		y = x;

		//���л���ȥ
		swap(pRow[iMax_Row], pRow[y]);

		for (y = x + 1; y < oA.m_iRow_Count; y++)
		{//������Ԫѭ��

			Solve_Linear_Gause_Index<_T> oRow_To_Elim_Index = pRow[y];
			if (oRow_To_Elim_Index.m_iHead_Coord > (unsigned int)x)
				continue;
			poMax_Row_Item = pMax_Row;
			poPre = poItem_To_Elim = &pBuffer[oRow_To_Elim_Index.m_iHead_Index];
			fValue = poItem_To_Elim->m_fValue;
			poItem_To_Elim = poItem_To_Elim->m_iNext ? &pBuffer[poItem_To_Elim->m_iNext] : NULL;
			while (1)
			{
				if (iCur_Item >= iMax_Item_Count)
				{
					printf("Insufficient space of oA1 in Solve_Linear_Gause_1\n");
					bSuccess = 0;
					goto END;
				}
				//�ָ������
				if (!poItem_To_Elim)
				{//�ӽ�ȥ
					pBuffer[iCur_Item] = { -fValue * poMax_Row_Item->m_fValue,poMax_Row_Item->x };
					poPre->m_iNext = iCur_Item++;
					if (!poMax_Row_Item->m_iNext)
						break;
					poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
					poPre = &pBuffer[poPre->m_iNext];
				}else if (poItem_To_Elim->x < poMax_Row_Item->x)
				{//������� ����Ԫ��Ԫ�����ҽ�һ��                 
					poPre = poItem_To_Elim;
					poItem_To_Elim = poItem_To_Elim->m_iNext ? &pBuffer[poItem_To_Elim->m_iNext] : NULL;
				}else if (poMax_Row_Item->x < poItem_To_Elim->x)
				{//��1��Ԫ�أ�����Ԫ��Ԫ��������һ��
					pBuffer[iCur_Item] = { -fValue * poMax_Row_Item->m_fValue,poMax_Row_Item->x, poPre->m_iNext };
					poPre->m_iNext = iCur_Item++;
					if (!poMax_Row_Item->m_iNext)
						break;
					poPre = &pBuffer[iCur_Item - 1];
					poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
				}else
				{//��ȣ�����ֵ��Ȼ����������ͬʱ�����ƶ�һ��
					poItem_To_Elim->m_fValue -= fValue * poMax_Row_Item->m_fValue;
					if (!poMax_Row_Item->m_iNext)
						break;
					poPre = poItem_To_Elim;
					poItem_To_Elim = poItem_To_Elim->m_iNext ? &pBuffer[poItem_To_Elim->m_iNext] : NULL;
					poMax_Row_Item = &pBuffer[poMax_Row_Item->m_iNext];
				}
			}

			poItem_To_Elim = &pBuffer[oRow_To_Elim_Index.m_iHead_Index];
			if (poItem_To_Elim->m_iNext)
			{
				oRow_To_Elim_Index.m_iHead_Index = poItem_To_Elim->m_iNext;
				poItem_To_Elim = &pBuffer[oRow_To_Elim_Index.m_iHead_Index];
				oRow_To_Elim_Index.m_iHead_Coord = poItem_To_Elim->x;
				oRow_To_Elim_Index.m_fHead_Value = poItem_To_Elim->m_fValue;
				pRow[y] = oRow_To_Elim_Index;
			}else
			{
				printf("������ in Solve_Linear_Gause_2\n");
				bSuccess = 0;
				goto END;
			}
		}
	}

	//����������ʽ�ش�
	for (y = oA.m_iRow_Count - 1; y >= 0; y--)
	{
		oRow_Index = pRow[y];
		if (!oRow_Index.m_iHead_Index)
		{
			X[y] = 0;
			continue;
		}
		poRow_Item = &pBuffer[oRow_Index.m_iHead_Index];
		fValue = 0;

		while (1)
		{
			if (poRow_Item->x <(unsigned int)oA.m_iCol_Count)
			{
				fValue += poRow_Item->m_fValue * X[poRow_Item->x];
				if (poRow_Item->m_iNext)
					poRow_Item = &pBuffer[poRow_Item->m_iNext];
				else
					break;
			}else
				break;
		}
		if (poRow_Item->x == oA.m_iCol_Count)
			X[y] = poRow_Item->m_fValue - fValue;
		else
			X[y] = -fValue;
	}
	
END:
	if (pRow)
		Free(&oMem_Mgr, pRow);
	if (pBuffer)
		Free(&oMem_Mgr, pBuffer + 1);
	if (pbSuccess)
		*pbSuccess = bSuccess;
	//exit(0);
	return;
}

template<typename _T>void Solve_Linear_Gause_1(Sparse_Matrix<_T>oA2, _T B[], _T X[], int* pbSuccess)
{//����Ҫ����
	Sparse_Matrix<_T> oA;
	Init_Sparse_Matrix(&oA, oA2.m_iCur_Item * 20, oA2.m_iCol_Count + 1, oA2.m_iRow_Count);
	Copy_Sparse_Matrix(oA2, &oA);

	int x, y, * pQ = NULL, * pQ_Reverse = NULL,		//����������pQ_Reverse[3] ��ʾ��ʵ�еĵ�3����Q���е�������
		bSuccess = 1;
	const _T eps = (_T)1e-10;
	typename Sparse_Matrix<_T>::Item* poItem, ** pCol_All, oItem, * pMax_Row, * poMax_Row_Item, * poItem_To_Elim, * poPre;

	pQ = (int*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(int));
	pQ_Reverse = (int*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(int));
	pCol_All = (typename Sparse_Matrix<_T>::Item **)pMalloc(&oMem_Mgr,oA.m_iCol_Count * 8);
	_T* B1 = (_T*)pMalloc(&oMem_Mgr, oA.m_iCol_Count * sizeof(_T));
	for (y = 0; y < oA.m_iRow_Count; y++)
	{
		if (B[y] != 0)
		{
			oA.m_pBuffer[oA.m_iCur_Item] = { (unsigned int)oA.m_iCol_Count - 1,(unsigned int)y,B[y],0 };
			if (oA.m_pRow[y])
			{
				poItem = &oA.m_pBuffer[oA.m_pRow[y]];
				while (poItem->m_iRow_Next)
					poItem = &oA.m_pBuffer[poItem->m_iRow_Next];
				poItem->m_iRow_Next = oA.m_iCur_Item++;
			}
			else
				oA.m_pRow[y] = oA.m_iCur_Item++;//����û��Ԫ��
		}
		pQ_Reverse[y] = pQ[y] = y;
	}
	Build_Link_Col(oA);

	int iMax_Item, iCur_Item, iMax_Q_Row, iTemp;
	_T fMax, fValue;

	for (x = 0; x < oA.m_iRow_Count; x++)
	{//���ⲻ��Reverse Lookup
		if (!oA.m_pCol[x])
		{
			printf("������ in Solve_Linear_Gause_1\n");
			bSuccess = 0;
			goto END;
		}
		iCur_Item = iMax_Item = oA.m_pCol[x];
		oItem = oA.m_pBuffer[iMax_Item];

		//����һ��
		y = x;
		fMax = 0;
		while (1)
		{
			if (pQ_Reverse[oItem.y] >= y)   //ֻҪ��ǰʵ������Q[y]֮��Ϳ����ж�
			{
				if (Abs(oItem.m_fValue) > Abs(fMax))
				{
					iMax_Q_Row = pQ_Reverse[oItem.y];
					fMax = oItem.m_fValue;
					iMax_Item = iCur_Item;
				}
			}
			if (oItem.m_iCol_Next)
			{
				iCur_Item = oItem.m_iCol_Next;
				oItem = oA.m_pBuffer[oItem.m_iCol_Next];
			}else
				break;
		}

		if (abs(fMax) <= eps)
		{
			printf("������ in Solve_Linear_Gause_1\n");
			bSuccess = 0;
			goto END;
		}
		//���Խ���Q[y]��
		iTemp = pQ[y];
		pQ[y] = pQ[iMax_Q_Row];
		pQ[iMax_Q_Row] = iTemp;
		pQ_Reverse[pQ[y]] = y;
		pQ_Reverse[iTemp] = iMax_Q_Row;

		//������Ԫ�����и���
		poItem = &oA.m_pBuffer[iMax_Item];
		fValue = poItem->m_fValue;
		poItem->m_fValue = 1.f;

		//�˴�Ϊ�˾����ܼ���0Ԫ������������ֱ�ӽ�m_pRowָ���������ǰ���Ѿ���ȥ��Ԫ��
		oA.m_pRow[poItem->y] = iMax_Item;
		//�˴���ȫ����ȥ�����������ܸ���

		//��Ԫ�е���ֵ		
		//poPre = poItem;	//����Ϊ����0Ԫ��׼��

		if (poItem->m_iRow_Next)
			poItem =pMax_Row = &oA.m_pBuffer[poItem->m_iRow_Next];
		else
			continue;	//�˴������飬�о���ʹ���ﲻ�������������һB��Ϊ0��©����Ԫ
				
		while (1)
		{
			poItem->m_fValue /= fMax;

			//////��ʱ���ԣ�������0Ԫ
			//if (abs(poItem->m_fValue) < 0.00001f)
			//{
			//	printf("*");
			//	poPre->m_iRow_Next = poItem->m_iRow_Next;	//����
			//}	
			//else
			//	poPre = poItem;

			if (poItem->m_iRow_Next)
				poItem = &oA.m_pBuffer[poItem->m_iRow_Next];
			else
				break;
		}

		//����pCol_All
		for (int i = x + 1; i < oA.m_iCol_Count; i++)
			pCol_All[i] = &oA.m_pBuffer[oA.m_pCol[i]];

		//Disp(oA,"A");
		//Ȼ������һ�ж������н�����Ԫ����ͷ����
		poItem = &oA.m_pBuffer[oA.m_pCol[x]];
		while (1)
		{
			if (pQ_Reverse[poItem->y] > y)
			{//����Ҫ��Ԫ, ��Max_Row����
				fValue = poItem->m_fValue;
				poItem->m_fValue = 0.f;

				//�˴��鷳���������н�����
				poPre = poItem;
				if (poItem->m_iRow_Next)
				{
					poItem_To_Elim = &oA.m_pBuffer[poItem->m_iRow_Next];

					//Ϊ�˽�һ������0Ԫ���˴�ֱ�ӽ�m_pRow��ֵָ��poItem_To_Elim
					if (pQ_Reverse[poItem_To_Elim->y] > y && poItem_To_Elim->x == x + 1)
						oA.m_pRow[poItem->y] = poItem->m_iRow_Next;
					//�˴���ȫ����ȥ�����������ܸ���
				}else
					poItem_To_Elim = NULL;

				poMax_Row_Item = pMax_Row;
				while (1)
				{//���ѭ����һ����Ԫ��
					if (oA.m_iCur_Item >= oA.m_iMax_Item_Count)
					{
						printf("Insufficient space of oA1 in Solve_Linear_Gause_1\n");
						bSuccess = 0;
						goto END;
					}
					//�ָ������
					if (!poItem_To_Elim)
					{//�ӽ�ȥ
						oA.m_pBuffer[oA.m_iCur_Item] = { poMax_Row_Item->x,poPre->y, -fValue * poMax_Row_Item->m_fValue };
						Add_Col_Link(oA, pCol_All, &oA.m_pBuffer[oA.m_iCur_Item]);
						poPre->m_iRow_Next = oA.m_iCur_Item++;
						if (!poMax_Row_Item->m_iRow_Next)
							break;
						poMax_Row_Item = &oA.m_pBuffer[poMax_Row_Item->m_iRow_Next];
						poPre = &oA.m_pBuffer[poPre->m_iRow_Next];
					}else if (poItem_To_Elim->x < poMax_Row_Item->x)
					{//������� ����Ԫ��Ԫ�����ҽ�һ��                 
						poPre = poItem_To_Elim;
						poItem_To_Elim = poItem_To_Elim->m_iRow_Next ? &oA.m_pBuffer[poItem_To_Elim->m_iRow_Next] : NULL;
					}else if (poMax_Row_Item->x < poItem_To_Elim->x)
					{//��1��Ԫ�أ�����Ԫ��Ԫ��������һ��
						oA.m_pBuffer[oA.m_iCur_Item] = { poMax_Row_Item->x, poPre->y, -fValue * poMax_Row_Item->m_fValue,poPre->m_iRow_Next };
						Add_Col_Link(oA, pCol_All, &oA.m_pBuffer[oA.m_iCur_Item]);
						poPre->m_iRow_Next = oA.m_iCur_Item++;
						if (!poMax_Row_Item->m_iRow_Next)
							break;
						poPre = &oA.m_pBuffer[oA.m_iCur_Item - 1];
						poMax_Row_Item = &oA.m_pBuffer[poMax_Row_Item->m_iRow_Next];
					}else
					{//��ȣ�����ֵ��Ȼ����������ͬʱ�����ƶ�һ��
						poItem_To_Elim->m_fValue -= fValue * poMax_Row_Item->m_fValue;
						if (!poMax_Row_Item->m_iRow_Next)
							break;
						poPre = poItem_To_Elim;
						poItem_To_Elim = poItem_To_Elim->m_iRow_Next ? &oA.m_pBuffer[poItem_To_Elim->m_iRow_Next] : NULL;
						poMax_Row_Item = &oA.m_pBuffer[poMax_Row_Item->m_iRow_Next];
					}
				}
			}

			//Disp(oA, "A");
			if (poItem->m_iCol_Next)
			{
				//Ϊ������ϴ��룬��poItem->m_iCol_NextҲҪ����
				poPre = poItem;
				poItem = &oA.m_pBuffer[poItem->m_iCol_Next];
				poPre->m_iCol_Next = (int)NULL;
			}else
				break;
		}
	}
	Re_Arrange_Sparse_Matrix(&oA);
	if (!oA.m_pBuffer)
	{
		bSuccess = 0;
		goto END;
	}

	//��B����X
	memset(B1, 0, (oA.m_iCol_Count - 1) * sizeof(_T));
	if (oA.m_pCol[oA.m_iCol_Count - 1])
	{
		poItem = &oA.m_pBuffer[oA.m_pCol[oA.m_iCol_Count - 1]];
		while (1)
		{
			B1[poItem->y] = poItem->m_fValue;
			if (poItem->m_iCol_Next)
				poItem = &oA.m_pBuffer[poItem->m_iCol_Next];
			else
				break;
		}
	}

	//Disp(oA);		
	for (x = oA.m_iCol_Count - 2; x > 0; x--)
	{//���и�,�����ϲ����ж���ô�࣬ǰ���Ѿ���֤�����ȣ����Ժܶ��жϿ��Զ���
		poItem = &oA.m_pBuffer[oA.m_pCol[x]];
		//����һֱ�㵽poItem->y!=xΪֹ
		y = x;
		fValue = B1[pQ[y]];
		while (1)
		{
			if (pQ_Reverse[poItem->y] < x)
			{
				B1[poItem->y] += -fValue * poItem->m_fValue;
				poItem->m_fValue = 0;
			}	
			if (poItem->m_iCol_Next)
				poItem = &oA.m_pBuffer[poItem->m_iCol_Next];
			else
				break;
		}
	}
	for (y = 0; y < oA.m_iRow_Count; y++)
		X[y] = B1[pQ[y]];
	Test_Linear(oA2, B, X);
END:
	if (pbSuccess)
		*pbSuccess = bSuccess;
	Free_Sparse_Matrix(&oA);
	if (pQ)
		Free(&oMem_Mgr, pQ);
	if (pQ_Reverse)
		Free(&oMem_Mgr, pQ_Reverse);
	if (B1)
		Free(&oMem_Mgr, B1);
	if (pCol_All)
		Free(&oMem_Mgr, pCol_All);
		
	return;
}

template<typename _T>void Solve_Linear_Gause(Sparse_Matrix<_T> oA, _T B[], _T X[], int* pbSuccess)
{//ֻ��ϡ�����A��ϡ�������������������
	//�˺���������������ֻ���ڲο�
	typedef struct Reverse_Item {
		unsigned int m_iQ_Index : 31;	//�����Ѿ���Q������Ǹ�Index;
		unsigned int m_bDone : 1;		//�����Ƿ��Ѿ���������Ԫ�Ĵ���
	}Reverse_Item;
	int y, y1, x, x1, i;
	int iMax, iTemp, * pQ, bSuccess = 1, iRank = 0;
	Reverse_Item* pQ_Reverse;
	const _T eps = (_T)1e-10;
	typename Sparse_Matrix<_T>::Item* pMax_Row, * poMax_Row_Item = NULL;
	//Ҫ��һ��Q��Reverse look up table�������У���һ��pQ_Done����
	_T fMax, fValue;
	//�����־�������Ա�־ĳһ���Ƿ��Ѿ���������Ԫ�Ĵ���ϡ������ܷ���������Ԫ��������
	//unsigned char* pQ_Done = (unsigned char*)pMalloc(&oMem_Mgr, ((oA.m_iRow_Count + 7) >> 3));
	pQ = (int*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(int));
	pQ_Reverse = (Reverse_Item*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(Reverse_Item));

	Sparse_Matrix<_T> oA1;
	typename Sparse_Matrix<_T>::Item oItem, * poItem;
	//Init_Sparse_Matrix(&oA1, oA.m_iMax_Item_Count + oA.m_iRow_Count, oA.m_iCol_Count + 1, oA.m_iRow_Count);
	Init_Sparse_Matrix(&oA1, oA.m_iCur_Item*300, oA.m_iCol_Count + 1, oA.m_iRow_Count);
	if (/*!pQ_Done ||*/ !pQ || !oA1.m_pBuffer)
	{
		printf("Fail to allocate in Solve_Linear_Gause\n");
		goto END;
	}

	//��oA,B ����oA1��
	for (y = 0; y < oA1.m_iRow_Count; y++)
	{
		if (oA.m_pRow[y])
		{//������
			oItem = oA.m_pBuffer[oA.m_pRow[y]];
			oA1.m_pRow[y] = oA1.m_iCur_Item;
			while (1)
			{
				if (oItem.m_iRow_Next)
				{//����
					oA1.m_pBuffer[oA1.m_iCur_Item++] = { oItem.x,oItem.y,oItem.m_fValue,oA1.m_iCur_Item + 1 };
					oItem = oA.m_pBuffer[oItem.m_iRow_Next];
				}else
				{
					oA1.m_pBuffer[oA1.m_iCur_Item++] = { oItem.x,oItem.y,oItem.m_fValue,0 };
					break;
				}
			}
		}

		//�ٰ�B[y]�ӽ�ȥ
		if (B[y] != 0)
		{
			if (oA.m_pRow[y]) //��������У��ӵ���β
				oA1.m_pBuffer[oA1.m_iCur_Item - 1].m_iRow_Next = oA1.m_iCur_Item;
			else //û�о��Լ���һ��
				oA1.m_pRow[y] = oA1.m_iCur_Item;
			oA1.m_pBuffer[oA1.m_iCur_Item++] = { (unsigned int)oA.m_iRow_Count,(unsigned int)y,B[y] };
		}
		pQ[y] = y;	//ÿ����Ԫ���ڵ���
		pQ_Reverse[y] = { (unsigned int)y,0 };
	}
	Build_Link_Col(oA1);
		
	//Disp(oA1, "A1");
	for (x = 0; x < oA1.m_iCol_Count - 1; x++)
	{
		//if (x % 100 == 0)
			//printf("x:%d\n", x);
		if (!oA1.m_pCol[x])
		{
			printf("������\n");
			bSuccess = 0;
			goto END;
		}
		oItem = oA1.m_pBuffer[oA1.m_pCol[x]];
		int iBit;
		fMax = 0;

		while (1)
		{
			//�������Ƿ��Ѿ���������Ԫ�жϣ����겻��
			//iBit= bGet_Bit(pQ_Done, oItem.y);
			iBit = pQ_Reverse[oItem.y].m_bDone;
			if (!iBit)
			{
				if (abs(oItem.m_fValue) > abs(fMax))
				{
					fMax = oItem.m_fValue;
					iMax = oItem.y;
				}
			}
			if (oItem.m_iCol_Next)
				oItem = oA1.m_pBuffer[oItem.m_iCol_Next];
			else
				break;
		}
		if (abs(fMax) < eps)
		{
			printf("������\n");
			bSuccess = 0;
			goto END;
		}
		//printf("x:%d iMax:%d %f\n",x, iMax, fGet_Value(&oA1, 5, 5));
		//��ʱ��iMax��ʾ��ʵ�У�����pQ�����У�����Ҫת��ΪpQ������
		iMax = pQ_Reverse[iMax].m_iQ_Index;
		y = x;	//ע�⣬yֻ��pQ��λ�ã��Ѿ�������Ԫ������
		//�����ԪSWAP��Q�ĵ�ǰλ����		
		iTemp = pQ[y];
		pQ[y] = pQ[iMax];
		pQ[iMax] = iTemp;
		//Disp(pQ, 1, 3, "Q");
		//��Reverse Lookup tableҲ����
		pQ_Reverse[pQ[y]].m_iQ_Index = y;
		pQ_Reverse[pQ[iMax]].m_iQ_Index = iMax;
		iRank++;

		//��iMax���ڵ��н���ϵ�����㣬��ϵ��/=A[y][y]
		poItem = &oA1.m_pBuffer[oA1.m_pRow[pQ[y]]];
		x1 = y;	//��ʾ�㵽��һ��, ���²���Ҫ���и��ֱ߽��жϣ��������ʴ��󣬽ṹ�д�����ǰ��
		while (poItem->x != x1)
			poItem = &oA1.m_pBuffer[poItem->m_iRow_Next];
		poItem->m_fValue = 1.f;
		pMax_Row = poItem = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
		while (poItem)
		{
			poItem->m_fValue /= fMax;
			poItem = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
		}
		//�Ժ��������д���
		for (i = y + 1; i < oA1.m_iRow_Count; i++)
		{//i��ʾ��i��
			y1 = pQ[i];
			//if (pQ[y] == 2 && y1 == 5)
				//printf("%f\n", fGet_Value(&oA1, 12, 5));
			if (!oA1.m_pRow[y1])
				continue;
			//��y1�е�������Ԫ���Ժ�����޸�
			poItem = &oA1.m_pBuffer[oA1.m_pRow[y1]];
			//�ƶ�����Ԫ��
			while (poItem && (int)poItem->x < x)
				poItem = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
			if (!poItem || poItem->x != x)
				continue;
			//�Ѿ��ƶ�����Ԫ�л�֮��
			fValue = poItem->m_fValue;
			if (fValue == 0)
				continue;
			poItem->m_fValue = 0.f;

			//�˴����Խ�m_pRow[y1]ָ����һ������λ�ã��Ա����
			oA1.m_pRow[y1] = poItem->m_iRow_Next;

			//���Խ������������ʾ����10-20%
			//typename Sparse_Matrix<_T>::Item* poPre_Item = poItem;

			poItem = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
			poMax_Row_Item = pMax_Row;
			while (poMax_Row_Item)
			{
				//������Ҫ����
				if (!poItem)
				{//��poItem�����еĶ�Ӧ��λ���ϵ�������Ԫ��
					Set_Value(&oA1, poMax_Row_Item->x, y1, -fValue * poMax_Row_Item->m_fValue);
					poMax_Row_Item = poMax_Row_Item->m_iRow_Next ? &oA1.m_pBuffer[poMax_Row_Item->m_iRow_Next] : NULL;
					//poPre_Item = &oA1.m_pBuffer[oA1.m_iCur_Item - 1];
				}else if (poItem->x == poMax_Row_Item->x)
				{//���£�poItem �� poMax_Row_Item�������ƶ�
					poItem->m_fValue -= fValue * poMax_Row_Item->m_fValue;
					//poPre_Item = poItem;
					poItem = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
					poMax_Row_Item = poMax_Row_Item->m_iRow_Next ? &oA1.m_pBuffer[poMax_Row_Item->m_iRow_Next] : NULL;
				}else if (poItem->x < poMax_Row_Item->x)
				{//ͬ����poItem�����еĶ�Ӧ��λ���ϵ�������Ԫ�أ�poItem�����ƶ�
					//poPre_Item = poItem;
					poItem = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
				}else if (poItem->x > poMax_Row_Item->x)
				{//ͬ����poItem�����еĶ�Ӧ��λ���ϵ�������Ԫ�أ�poMax_Row_Item�����ƶ�
					//�˴�Ҫ�ص��Ż�������������
					Set_Value(&oA1, poMax_Row_Item->x, y1, -fValue * poMax_Row_Item->m_fValue);
					poMax_Row_Item = poMax_Row_Item->m_iRow_Next ? &oA1.m_pBuffer[poMax_Row_Item->m_iRow_Next] : NULL;
					//poPre_Item = &oA1.m_pBuffer[oA1.m_iCur_Item - 1];
				}
				if (oA1.m_iCur_Item >= oA1.m_iMax_Item_Count)
				{
					printf("Insufficient space in Solve_Linear_Gause\n");
					bSuccess = 0;
					goto END;
				}
			}
		}
		pQ_Reverse[pQ[y]].m_bDone = 1;
	}

	//�˴�������һ������ɾ��oA1�����еķ�0ֵ
	Re_Arrange_Sparse_Matrix(&oA1);
	typename Sparse_Matrix<_T>::Item* poB_Item;
	//Disp(oA1, "A1");

	//�ش�����Q[iOrder - 1]��ʼ�ش���������һ�����ϻش�
	for (y = oA1.m_iRow_Count - 1; y >= 0; y--)
	{
		poItem = poB_Item = &oA1.m_pBuffer[oA1.m_pRow[pQ[y]]];
		while (poB_Item->m_iRow_Next)
			poB_Item = &oA1.m_pBuffer[poB_Item->m_iRow_Next];
		fValue = poB_Item->x == oA1.m_iCol_Count - 1 ? poB_Item->m_fValue : 0;
		X[poItem->x] = fValue;
		if (fValue == 0)
			continue;
		//����������Ԫ����λ����Ԫ������
		poItem = &oA1.m_pBuffer[oA1.m_pCol[y]];
		while (1)
		{
			if (poItem->y != pQ[y])
			{
				//���Ź�ȥѰ��B��
				if (poItem->m_iRow_Next)
				{
					poB_Item = &oA1.m_pBuffer[poItem->m_iRow_Next];
					while (poB_Item->m_iRow_Next)
						poB_Item = &oA1.m_pBuffer[poB_Item->m_iRow_Next];
					if (poB_Item->x == oA1.m_iCol_Count - 1)
						poB_Item->m_fValue += -poItem->m_fValue * fValue;
					else//��һ����Item
						Set_Value(&oA1, poItem->x, poItem->y, -poItem->m_fValue * fValue);
				}else
					Set_Value(&oA1, poItem->x, poItem->y, -poItem->m_fValue * fValue);
				poItem->m_fValue = 0;
			}
			if (poItem->m_iCol_Next)
				poItem = &oA1.m_pBuffer[poItem->m_iCol_Next];
			else
				break;
		}
		//Disp(oA1,"A1");
	}
END:
	if (pQ)
		Free(&oMem_Mgr, pQ);
	if (pQ_Reverse)
		Free(&oMem_Mgr, pQ_Reverse);
	Free_Sparse_Matrix(&oA1);
	//if (pQ_Done)
		//Free(&oMem_Mgr, pQ_Done);
	if (pbSuccess)
		*pbSuccess = bSuccess;
	return;
}
template<typename _T>void Dense_2_Sparse(_T A[], int m, int n, Sparse_Matrix<_T>* poA,int bCompact)
{//�����ܿ콫һ�����ת��Ϊϡ������ʽ���ٶ�oA��֮ǰ�Ѿ���ʼ���ÿռ�
	int y, x, bLine_Add, iPos = 0;
	Sparse_Matrix<_T> oA = *poA;
	if (!oA.m_iMax_Item_Count)
	{
		printf("Sparse Matrix is not initialized yet\n");
		return;
	}
	Reset_Sparse_Matrix(&oA);	//�������
	for (y = 0; y < m; y++)
	{
		for (bLine_Add = x = 0; x < n; x++, iPos++)
		{
			if (A[iPos] != 0)
			{
				if (!bLine_Add)
				{//��δ������
					oA.m_pRow[y] = oA.m_iCur_Item;
					bLine_Add = 1;
				}
				oA.m_pBuffer[oA.m_iCur_Item++] = { (unsigned int)x,(unsigned int)y,A[iPos],oA.m_iCur_Item + 1 };
			}
		}
		if (bLine_Add)
			oA.m_pBuffer[oA.m_iCur_Item - 1].m_iRow_Next = 0;
	}
	Build_Link_Col(oA);
	if (bCompact)
		Compact_Sparse_Matrix(&oA);
	*poA = oA;
	//Disp(oA);
}
template<typename _T>void Sparse_2_Dense(Sparse_Matrix<_T> oA, _T B[])
{//��ϡ�����ת��Ϊ���ܣ�����������
	int i;
	typename Sparse_Matrix<_T>::Item oItem;
	if (oA.m_iRow_Count == 1)
	{//������
		if (oA.m_pRow[0])
		{
			oItem = oA.m_pBuffer[oA.m_pRow[0]];
			for (i = 0; i < oA.m_iCol_Count; i++)
			{
				if (oItem.x == i)
				{
					B[i] = oItem.m_fValue;
					if (oItem.m_iRow_Next)
						oItem = oA.m_pBuffer[oItem.m_iRow_Next];
					else
					{
						i++;
						break;
					}
				}else
					B[i] = 0;
			}
			memset(&B[i], 0, (oA.m_iCol_Count - i) * sizeof(_T));
		}else
			memset(B, 0, oA.m_iCol_Count * sizeof(_T));

	}else if (oA.m_iCol_Count == 1)
	{//������
		if (oA.m_pCol[0])
		{
			oItem = oA.m_pBuffer[oA.m_pCol[0]];
			for (i = 0; i < oA.m_iRow_Count; i++)
			{
				if (oItem.y == i)
				{
					B[i] = oItem.m_fValue;
					if (oItem.m_iCol_Next)
						oItem = oA.m_pBuffer[oItem.m_iCol_Next];
					else
					{
						i++;
						break;
					}
				}else
					B[i] = 0;
			}
			memset(&B[i], 0, (oA.m_iRow_Count - i) * sizeof(_T));
		}else
			memset(B, 0, oA.m_iRow_Count * sizeof(_T));
	}else
	{
		memset(B, 0, oA.m_iRow_Count * oA.m_iCol_Count * sizeof(_T));
		for (i = 0; i < oA.m_iRow_Count; i++)
		{
			if (oA.m_pRow[i])
			{//����������
				oItem = oA.m_pBuffer[oA.m_pRow[i]];
				while (1)
				{
					if (oItem.m_fValue != 0)
						B[oItem.y * oA.m_iCol_Count + oItem.x] = oItem.m_fValue;
					if (oItem.m_iRow_Next)
						oItem = oA.m_pBuffer[oItem.m_iRow_Next];
					else
						break;
				}
			}
		}
	}
		
	return;
}
template<typename _T>void Add_I_Matrix(_T A[], int iOrder,_T ramda)
{//ר��������-��䷽����������A�������һ��I��
	int i;
	for (i = 0; i < iOrder; i++)
		A[i * iOrder + i] += ramda;
	return;
}
template<typename _T>void Add_I_Matrix(Sparse_Matrix<_T>* poA, int* pbSuccess,_T ramda)
{//���������������Ҫ����Sigma_H����ʹ�ÿ��棬����һ��I����
	int y;
	Sparse_Matrix<_T> oA = *poA;
	typename Sparse_Matrix<_T>::Item* poItem;
	for (y = 0; y < oA.m_iRow_Count; y++)
	{
		if (!oA.m_pRow[y])
		{//����ȫ0���ü�һ��
			if (oA.m_iCur_Item >= oA.m_iMax_Item_Count)
			{
				printf("Exceed max Item count in Add_I_Matrix\n");
				return;
			}
			oA.m_pRow[y] = oA.m_iCur_Item;
			oA.m_pBuffer[oA.m_iCur_Item++] = { (unsigned int)y,(unsigned int)y,ramda };
		}else
		{//һֱ�Ƶ������λ��
			poItem = &oA.m_pBuffer[oA.m_pRow[y]];
			if (poItem->x > (unsigned int)y)
			{//�ӵ�����ͷ
				if (oA.m_iCur_Item >= oA.m_iMax_Item_Count)
				{
					printf("Exceed max Item count in Add_I_Matrix\n");
					return;
				}
				oA.m_pBuffer[oA.m_iCur_Item] = { (unsigned int)y,(unsigned int)y,ramda,(int)oA.m_pRow[y] };
				oA.m_pRow[y] = oA.m_iCur_Item++;
			}else
			{
				while (1)	// (poItem && poItem->x < y)
				{
					if (poItem->x == y)
					{
						poItem->m_fValue += ramda;
						break;
					}
					else if (!poItem->m_iRow_Next || (poItem->m_iRow_Next && (int)oA.m_pBuffer[poItem->m_iRow_Next].x > y))
					{//����
						if (oA.m_iCur_Item >= oA.m_iMax_Item_Count)
						{
							printf("Exceed max Item count in Add_I_Matrix\n");
							return;
						}
						oA.m_pBuffer[oA.m_iCur_Item] = { (unsigned int)y,(unsigned int)y,ramda,poItem->m_iRow_Next };
						poItem->m_iRow_Next = oA.m_iCur_Item++;
						break;
					}
					poItem = &oA.m_pBuffer[poItem->m_iRow_Next];
				}
			}
		}
	}
	Build_Link_Col(oA);
	*poA = oA;
}
template<typename _T>void Reset_Sparse_Matrix(Sparse_Matrix<_T>* poA)
{
	Sparse_Matrix<_T> oA = *poA;
	memset(oA.m_pRow, 0, oA.m_iRow_Count * sizeof(int));
	memset(oA.m_pCol, 0, oA.m_iCol_Count * sizeof(int));
	oA.m_iCur_Item = 1;
	*poA = oA;
}
template<typename _T>void Disp_Fillness(Sparse_Matrix<_T> oMatrix)
{
	int y, x;
	typename Sparse_Matrix<_T>::Item* poCur;
	for (y = 0; y < oMatrix.m_iRow_Count; y++)
	{
		if (!oMatrix.m_pRow[y])
			for (x = 0; x < oMatrix.m_iCol_Count; x++)
				printf(".");
		else
		{
			poCur = &oMatrix.m_pBuffer[oMatrix.m_pRow[y]];
			x = 0;
			do {
				for (; x < (int)poCur->x; x++)
					printf(".");
				if (poCur->m_fValue)
					printf("*");
				else
					printf(".");
				x++;
				if (poCur->m_iRow_Next)
					poCur = &oMatrix.m_pBuffer[poCur->m_iRow_Next];
				else
					break;
			} while (1);
			for (; x < oMatrix.m_iCol_Count; x++)
				printf(".");
		}
		printf("\n");
	}
}
template<typename _T>void Disp(Sparse_Matrix<_T> oMatrix, const char Caption[])
{
	int y, x;
	typename Sparse_Matrix<_T>::Item* poCur;
	if (Caption)
		printf("%s\n", Caption);
	for (y = 0; y < oMatrix.m_iRow_Count; y++)
	{
		if (!oMatrix.m_pRow[y])
			for (x = 0; x < oMatrix.m_iCol_Count; x++)
				printf("%.8f, ", 0.f);
		else
		{
			poCur = &oMatrix.m_pBuffer[oMatrix.m_pRow[y]];
			x = 0;
			do {
				for (; x <(int)poCur->x; x++)
					printf("%.8f, ", 0.f);
				printf("%.8f, ", poCur->m_fValue);
				x++;
				if (poCur->m_iRow_Next)
					poCur = &oMatrix.m_pBuffer[poCur->m_iRow_Next];
				else
					break;
			} while (1);
			for (; x < oMatrix.m_iCol_Count; x++)
				printf("%.8f, ", 0.f);
		}
		printf("\n");
	}
}
template<typename _T>void Get_Inv_Matrix_Row_Op(Sparse_Matrix<_T> oA, Sparse_Matrix<_T>* poA_Inv, int* pbSuccess)
{//ϡ��������棬��֪��û���ã��ȸ����˵,������Reverse_Q
	int* pQ, * pQ_Reverse, y, x, bSuccess = 1, iCol_Count_Minuse_1, iCur_Item = 0;
	Sparse_Matrix<_T> oA1, oA_Inv = *poA_Inv;
	const _T eps = (_T)1e-6;
	typename Sparse_Matrix<_T>::Item oItem, * poItem, * poPre, * pMax_Row, * poMax_Row_Item, * poItem_To_Elim=NULL;
	if (oA.m_iCol_Count != oA.m_iRow_Count)
	{
		printf("������ in Get_Inv_Matrix_Row_Op\n");
		return;
	}
	pQ = (int*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(int));
	pQ_Reverse = (int*)pMalloc(&oMem_Mgr, oA.m_iRow_Count * sizeof(int));

	Init_Sparse_Matrix(&oA1, oA.m_iCur_Item*60, oA.m_iCol_Count * 2, oA.m_iRow_Count);
	iCol_Count_Minuse_1 = oA1.m_iCol_Count - 1;
	//������������ұ�Ϊ��λ����
	for (y = 0; y < oA.m_iRow_Count; y++)
	{
		if (!oA.m_pRow[y])
		{
			printf("������ in Get_Inv_Matrix_Row_Op\n");
			bSuccess = 0;
			goto END;
		}			
		oItem = oA.m_pBuffer[oA.m_pRow[y]];
		while (1)
		{
			if (oItem.m_fValue != 0)
			{
				if (!oA1.m_pRow[y])
					oA1.m_pRow[y] = oA1.m_iCur_Item;
				oA1.m_pBuffer[oA1.m_iCur_Item++] = { oItem.x,oItem.y,oItem.m_fValue,oA1.m_iCur_Item + 1 };
			}
			if (oItem.m_iRow_Next)
				oItem = oA.m_pBuffer[oItem.m_iRow_Next];
			else
				break;
		}
		if (!oA1.m_pRow[y])
		{
			printf("������ in Get_Inv_Matrix_Row_Op\n");
			bSuccess = 0;
			goto END;
		}
		oA1.m_pBuffer[oA1.m_iCur_Item++] = { (unsigned int)(oA.m_iCol_Count + y),(unsigned int)y,1,0 };
		pQ_Reverse[y] = pQ[y] = y;
	}
	Build_Link_Col(oA1);

	int iMax_Item, iMax_Q_Row, iTemp;
	_T fMax, fValue;
	//Disp(oA1,"A1");
	for (x = 0; x < oA1.m_iRow_Count; x++)
	{//���ⲻ��Reverse Lookup
		if (!oA1.m_pCol[x])
		{
			printf("������ in Get_Inv_Matrix_Row_Op\n");
			bSuccess = 0;
			goto END;
		}
		iCur_Item = iMax_Item = oA1.m_pCol[x];
		oItem = oA1.m_pBuffer[iMax_Item];

		//����һ��
		y = x;
		fMax = 0;
		while (1)
		{
			if (pQ_Reverse[oItem.y] >= y)   //ֻҪ��ǰʵ������Q[y]֮��Ϳ����ж�
			{
				if (Abs(oItem.m_fValue) > Abs(fMax))
				{
					iMax_Q_Row = pQ_Reverse[oItem.y];
					fMax = oItem.m_fValue;
					iMax_Item = iCur_Item;
				}
			}
			if (oItem.m_iCol_Next)
			{
				iCur_Item = oItem.m_iCol_Next;
				oItem = oA1.m_pBuffer[oItem.m_iCol_Next];
			}else
				break;
		}
		if (abs(fMax) <= eps)
		{
			printf("������ in Get_Inv_Matrix_Row_Op\n");
			bSuccess = 0;
			goto END;
		}

		//���Խ���Q[y]��
		iTemp = pQ[y];
		pQ[y] = pQ[iMax_Q_Row];
		pQ[iMax_Q_Row] = iTemp;
		pQ_Reverse[pQ[y]] = y;
		pQ_Reverse[iTemp] = iMax_Q_Row;

		//������Ԫ�����и���
		poItem = &oA1.m_pBuffer[iMax_Item];
		fValue = poItem->m_fValue;
		poItem->m_fValue = 1.f;

		//�˴�Ϊ�˾����ܼ���0Ԫ������������ֱ�ӽ�m_pRowָ���������ǰ���Ѿ���ȥ��Ԫ��
		oA1.m_pRow[poItem->y] = iMax_Item;
		//�˴���ȫ����ȥ�����������ܸ���

		//��Ԫ�е���ֵ
		if (poItem->m_iRow_Next)
			poItem = pMax_Row = &oA1.m_pBuffer[poItem->m_iRow_Next];
		else
			continue;

		//poItem = pMax_Row = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
		while (poItem)
		{
			poItem->m_fValue /= fMax;
			poItem = poItem->m_iRow_Next ? &oA1.m_pBuffer[poItem->m_iRow_Next] : NULL;
		}
		//Disp(oA1, "A1");
		//Ȼ������һ�ж������н�����Ԫ����ͷ����
		poItem = &oA1.m_pBuffer[oA1.m_pCol[x]];
		while (1)
		{
			if (poItem->y != pQ[y] && poItem->m_fValue != 0)
			{//����Ҫ��Ԫ, ��Max_Row����                
				fValue = poItem->m_fValue;
				poItem->m_fValue = 0.f;

				//�˴��鷳���������н�����
				poPre = poItem;
				if (poItem->m_iRow_Next)
				{
					poItem_To_Elim = &oA1.m_pBuffer[poItem->m_iRow_Next];
					
					//Ϊ�˽�һ������0Ԫ���˴�ֱ�ӽ�m_pRow��ֵָ��poItem_To_Elim
					if ( pQ_Reverse[poItem_To_Elim->y] >y && poItem_To_Elim->x == x + 1)
						oA1.m_pRow[poItem->y] = poItem->m_iRow_Next;
					//�˴���ȫ����ȥ�����������ܸ���

				}else
					poItem_To_Elim = NULL;

				poMax_Row_Item = pMax_Row;
				while (1)
				{//���ѭ����һ����Ԫ��
					if (oA1.m_iCur_Item >= oA1.m_iMax_Item_Count)
					{
						printf("Insufficient space of oA1 in Get_Inv_Matrix_Row_Op\n");
						bSuccess = 0;
						goto END;
					}
					//�ָ������
					if (!poItem_To_Elim)
					{//�ӽ�ȥ
						oA1.m_pBuffer[oA1.m_iCur_Item] = { poMax_Row_Item->x,poPre->y, -fValue * poMax_Row_Item->m_fValue };
						poPre->m_iRow_Next = oA1.m_iCur_Item++;
						if (!poMax_Row_Item->m_iRow_Next)
							break;
						poMax_Row_Item = &oA1.m_pBuffer[poMax_Row_Item->m_iRow_Next];
						poPre = &oA1.m_pBuffer[poPre->m_iRow_Next];
					}else if (poItem_To_Elim->x < poMax_Row_Item->x)
					{//�ò��� ����Ԫ��Ԫ�����ҽ�һ��
						//oA1.m_pBuffer[oA1.m_iCur_Item]= { poMax_Row_Item->x, poMax_Row_Item->y, -fValue * poMax_Row_Item->m_fValue,poPre->m_iRow_Next };
						//poPre->m_iRow_Next = oA1.m_iCur_Item;                        
						poPre = poItem_To_Elim;
						poItem_To_Elim = poItem_To_Elim->m_iRow_Next ? &oA1.m_pBuffer[poItem_To_Elim->m_iRow_Next] : NULL;
					}else if (poMax_Row_Item->x < poItem_To_Elim->x)
					{//��1��Ԫ�أ�����Ԫ��Ԫ��������һ��
						oA1.m_pBuffer[oA1.m_iCur_Item] = { poMax_Row_Item->x, poPre->y, -fValue * poMax_Row_Item->m_fValue,poPre->m_iRow_Next };
						poPre->m_iRow_Next = oA1.m_iCur_Item++;
						if (!poMax_Row_Item->m_iRow_Next)
							break;
						poPre = &oA1.m_pBuffer[oA1.m_iCur_Item-1];
						poMax_Row_Item = &oA1.m_pBuffer[poMax_Row_Item->m_iRow_Next];
					}else
					{//��ȣ�����ֵ��Ȼ����������ͬʱ�����ƶ�һ��
						poItem_To_Elim->m_fValue -= fValue * poMax_Row_Item->m_fValue;
						if (!poMax_Row_Item->m_iRow_Next)
							break;
						poPre = poItem_To_Elim;
						poItem_To_Elim = poItem_To_Elim->m_iRow_Next ? &oA1.m_pBuffer[poItem_To_Elim->m_iRow_Next] : NULL;
						poMax_Row_Item = &oA1.m_pBuffer[poMax_Row_Item->m_iRow_Next];
					}
				}
			}
			//Disp(oA1, "A1");
			if (poItem->m_iCol_Next)
				poItem = &oA1.m_pBuffer[poItem->m_iCol_Next];
			else
				break;
		}
		Build_Link_Col(oA1);    //�˴������ٸĽ�������Ԫ�����Ҹ㣬ǰ�治��
	}

	//Disp(oA1, "A1");
	//��oA1�ұ߲��ֳ���poA_Inv�ϼ���
	oA_Inv = *poA_Inv;
	if(!oA_Inv.m_pBuffer)
		Init_Sparse_Matrix(&oA_Inv, oA1.m_iMax_Item_Count, oA.m_iCol_Count, oA.m_iRow_Count);
	for (y = 0; y < oA1.m_iRow_Count; y++)
	{
		oItem = oA1.m_pBuffer[oA1.m_pRow[pQ[y]]];
		//���ƶ�����������λ����
		while (1)
		{
			if (oItem.x >= (unsigned int)oA_Inv.m_iCol_Count || !oItem.m_iRow_Next)
				break;
			oItem = oA1.m_pBuffer[oItem.m_iRow_Next];
		}
		oA_Inv.m_pRow[y] = oA_Inv.m_iCur_Item;
		while (1)
		{
			if (oA_Inv.m_iCur_Item >= oA_Inv.m_iMax_Item_Count)
			{
				printf("Insufficient space of oA_Inv in Get_Inv_Matrix_Row_Op\n");
				bSuccess = 0;
				goto END;
			}
			if (oItem.m_fValue != 0)
				oA_Inv.m_pBuffer[oA_Inv.m_iCur_Item++] = { oItem.x - oA_Inv.m_iCol_Count,oItem.y,oItem.m_fValue,oA_Inv.m_iCur_Item + 1 };

			if (oItem.m_iRow_Next)
				oItem = oA1.m_pBuffer[oItem.m_iRow_Next];
			else
				break;
		}
		oA_Inv.m_pBuffer[oA_Inv.m_iCur_Item - 1].m_iRow_Next = 0;
	}
	Build_Link_Col(oA_Inv);
	//Disp(oA_Inv, "A_Inv");
	*poA_Inv = oA_Inv;
END:
	if (pbSuccess)
		*pbSuccess = bSuccess;
	Free_Sparse_Matrix(&oA1);
	if (pQ)
		Free(&oMem_Mgr, pQ);
	if (pQ_Reverse)
		Free(&oMem_Mgr, pQ_Reverse);
}
template<typename _T>void Rotation_Vector_4_2_3(_T V[], _T V1[])
{//��4ά��ת����ת��ά3ά��ת����

	for (int i = 0; i < 3; i++)
		V1[i] = V[i] * V[3];
}
template<typename _T>void Rotation_Vector_3_2_4(_T V[], _T V1[])
{//��3ά��ת����ת��Ϊ4ά��ת����
	_T V2[4];
	Normalize(V, 3, V2);
	V2[3] = fGet_Mod(V,3);
	memcpy(V1, V2, 4 * sizeof(_T));
}

template<typename _T>void Crop_Matrix(_T Source[], int m, int n, int x, int y, int w, int h, _T Dest[],int iDest_Stride)
{//����һ��m��n�еľ���Source, ��x,y���س�һ��h��w�е��ӿ����������Dest
//��iDest_Stride=0, ��Ŀ��Ĭ��δw*h
	if (iDest_Stride == 0)
		iDest_Stride = w;
	else if(w > iDest_Stride)
		w = iDest_Stride;
	_T* pSource_Cur = &Source[y * n + x],
		* pSource_End = pSource_Cur + h * n,
		* pDest_Cur = Dest;
	for (; pSource_Cur < pSource_End; pSource_Cur += n,pDest_Cur+=iDest_Stride)
		for (int i = 0; i < w; i++)
			pDest_Cur[i] = pSource_Cur[i];
	return;
}
template<typename _T>void Copy_Matrix_Partial(_T Source[], int iSource_Stride, int x, int y,_T Dest[], int m, int n)
{//�պ��෴���������һ���ֿ�����Ŀ�������
	int x1, y1;
	for (y1 = 0; y1 < m; y1++)
		for (x1 = 0; x1 < n; x1++)
			Dest[y1 * n + x1] = Source[(y1 + y) * iSource_Stride + x1 + x];
}
template<typename _T>void Matrix_Add_Partial(_T Source[], int iSource_Stride, int x1, int y1, int w, int h,
	_T Dest[], int iDest_Stride, int x2, int y2)
{//��Դ�ٳ�һ�鿽��Dest��Ӧ(x2,y2)��ʼ��λ����
	int i, j;
	for (i = 0; i < h; i++)
		for (j = 0; j < w; j++)
			Dest[(y2 + i) * iDest_Stride + x2 + j] += Source[(y1 + i) * iSource_Stride + x1 + j];
	return;
}
template<typename _T>void Copy_Matrix_Partial(_T Source[], int m, int n, _T Dest[], int iDest_Stride, int x, int y)
{//��Source���� Dest��Ӧλ����
	int x1, y1;
	//�˴����øĽ���Ҫ��Place_Image���������Ǹ������������̫�ֲ�
	for (y1 = 0; y1 < m; y1++)
		for (x1 = 0; x1 < n; x1++)
			Dest[(y1 + y) * iDest_Stride + x1 + x] = Source[y1 * n + x1];
	return;	
}

template<typename _T>void Disp_Fillness(_T A[], int m, int n)
{//��ʾһ������ķ������
	for (int y = 0; y < m; y++)
	{
		for (int x = 0; x < n; x++)
		{
			if (A[y * n + x] != 0)
				printf("*");
			else
				printf(".");
		}
		printf("\n");
	}
}

float fGet_e()
{//���ü�����ʽ���e�Ľ��ƽ�
	float e_pre=2, e=2;
	int n,n_perm=2;
	for (n = 2;;)
	{
		e += 1.f / n_perm;
		if (e == e_pre)
			return e;
		n++;
		n_perm *= n;
		e_pre = e;
	}
}
template<typename _T>void Solve_Linear_LLt(_T A[],int iOrder, _T* B, _T* X, int* pbSuccess)
{//��������������A�����Է��� Ax=b, ��LLt�ֽ�ֱ�ӷ���
//�����󵽣����������Ȼ�о������⣬ԭ���Ƿ�����Ԫ������ʹ��˹����Ԫ�������ھ������⣬����
//�������Ҫ�ŵ���ʵ�����в���
	int iResult;
	union {
		_T* L;
		_T* Lt;
	};
	L= (_T*)pMalloc(&oMem_Mgr, iOrder * iOrder * sizeof(_T));
	_T* y = (_T*)pMalloc(&oMem_Mgr, iOrder * sizeof(_T));

	Cholosky_Decompose(A, iOrder, L,&iResult);
	if(!iResult)
		goto END;
	//Disp(L, iOrder, iOrder, "L");
	//��ʱA = LLt ���� Ax=b => L*Lt*x = b, �� y= Lt*x, ���� Ly=b

	Solve_Linear_Gause(L, iOrder, B, y, &iResult);
	//��Ҫ�� Lt*x =y
	Matrix_Transpose(L, iOrder, iOrder, Lt);
	Solve_Linear_Gause(Lt, iOrder, y, X, &iResult);

END:
	if (L)
		Free(&oMem_Mgr, L);
	if (y)
		Free(&oMem_Mgr, y);
	if (pbSuccess)
		*pbSuccess = iResult;
	return;
}
template<typename _T>void Softmax(_T A[], int n, _T B[])
{//һάsoftmax
	_T fTotal = 0,fRecip;
	int i;
	for (i = 0; i < n; i++)
		fTotal+=(B[i] = exp(A[i]));
	fRecip = 1.f / fTotal;
	for (i = 0; i < n; i++)
		B[i] *= fRecip;
	return;
}

template<typename _T>void Matrix_2_R(_T Hr[3 * 3], _T R[3 * 3])
{//
	SVD_Info oSVD;
	int iResult;
	//Disp(H, 3, 3, "H");
	SVD_Alloc(3, 3, &oSVD, Hr);
	svd_3(Hr, oSVD, &iResult);
	//Test_SVD(Hr, oSVD, &iResult, 1e-7f);
	/*_T Ut[3*3], V[3*3];
	Matrix_Transpose((_T*)oSVD.U, 3, 3, Ut);
	Matrix_Transpose((_T*)oSVD.Vt, 3, 3, V);*/
	//R = VU'
	Matrix_Multiply((_T*)oSVD.U, 3, 3, (_T*)oSVD.Vt, 3, R);
	//Disp(R, 3, 3, "R");
	//printf("%d\n", bIs_Orthogonal(R, 3, 3));
	_T fDet = fGet_Determinant(R, 3);
	if (fDet < 0.f)
	{//����ϵ��Ҫ��V�ĵ�һ�г���-1
		printf("��������ϵ��Ҫ��");
		//ע�⣬Ӧ����V��������ȡ��������V'
		((_T*)oSVD.Vt)[6] *= -1, ((_T*)oSVD.Vt)[7] *= -1, ((_T*)oSVD.Vt)[8] *= -1;
		Matrix_Multiply((_T*)oSVD.U, 3, 3, (_T*)oSVD.Vt, 3, R);
		//Disp(R, 3, 3, "R");
	}
	//Disp((_T*)oSVD.S, 1, 3, "S");
	Free_SVD(&oSVD);
	return;
}
template<typename _T>int bTest_Eigen(_T A[], int n, _T fEigen_Value, _T Eigen_Vector[], _T eps)
{//����һ������ֵ�����������Ƿ�ΪA��
	_T* pTemp_1 = (_T*)pMalloc(n * sizeof(_T)),
		* pTemp_2 = (_T*)pMalloc(n * sizeof(_T));
	int bRet = 1;
	Matrix_Multiply(A, n, n, Eigen_Vector, 1, pTemp_1);
	Matrix_Multiply(Eigen_Vector, 1, n, fEigen_Value, pTemp_2);
	for (int i = 0; i < n; i++)
		pTemp_1[i] = Abs(pTemp_1[i]), pTemp_2[i] = Abs(pTemp_2[i]);

	if (sqrt(fGet_Distance(pTemp_1, pTemp_2, n)) > eps)
	{
		Disp(pTemp_1, 1, n, "Ax");
		Disp(pTemp_2, 1, n, "��x");
		printf("Error:%f\n", sqrt(fGet_Distance(pTemp_1, pTemp_2, n)));
		bRet = 0;
	}
	if (pTemp_1)Free(pTemp_1);
	if (pTemp_2)Free(pTemp_2);
	return bRet;
}

template<typename _T>static _T uniformRand(_T lowerBndr, _T upperBndr) {
	return lowerBndr +
		((_T)std::rand() / (RAND_MAX + 1.0)) * (upperBndr - lowerBndr);
}
template<typename _T> _T fGet_Random_No(_T mean, _T sigma) 
{//�����g20����������㷨���о��ܴ�
	_T y, r2;
	do {
		_T x =(_T)( - 1.0 + 2.0 * uniformRand(0.0, 1.0));
		y = (_T)(-1.0 + 2.0 * uniformRand(0.0, 1.0));
		r2 = x * x + y * y;
	} while (r2 > 1.0 || r2 == 0.0);
	return(_T)(mean + sigma * y * sqrt(-2.0 * log(r2) / r2));
}
template<typename _T>void Test_Inv_Matrix(_T A[], _T A_Inv[], int iOrder,_T *pfError_Sum)
{//����AA(-1)������I�Ա�
	int i;
	_T fError,* AA_Inv = (_T*)pMalloc(iOrder * iOrder * sizeof(_T));
	Matrix_Multiply(A, iOrder, iOrder, A_Inv, iOrder, AA_Inv);
	//Disp(AA_Inv, iOrder, iOrder, "AA_Inv");
	for (i = 0; i < iOrder; i++)
		AA_Inv[i * iOrder + i] -= 1;  
	fError = fGet_Mod(AA_Inv, iOrder * iOrder);
	if(pfError_Sum)
		*pfError_Sum +=fError;
	//printf("Error Sum:%lf\n",fError);
	Free(AA_Inv);
}
template<typename _T>void Matrix_Transpose(_T* A, int ma, int na, _T* At)
{//����ת�ã���ȫ������Դ���Ե���Ŀ��
	int y, x;
	_T* At_1;
	if (abs((int)(A - At)) < ma * na)
		At_1 = (_T*)pMalloc(&oMem_Mgr, ma * na * sizeof(_T));
	else
		At_1 = At;
	for (y = 0; y < ma; y++)
		for (x = 0; x < na; x++)
			At_1[x * ma + y] = A[y * na + x];
	if (abs((int)(A - At)) < ma * na)
	{
		memcpy(At, At_1, ma * na * sizeof(_T));
		Free(&oMem_Mgr, At_1);
	}
	return;
}

/**********************һЩ������λ�û���*******************/
int bGet_Matrix_xy(int iPos, int n, int* px, int* py)
{//����һ������Ľ�n,��һ��˳��λ�ã������ھ���������
	//��󲻵�Ϊ  (n + n-y-1)*y/2 <=N
	//������h, ��h-1�еĸ���Ϊ n-(h-1)
	//����Ϊ [n + n-(h-1)]*h/2 = (2n+1 -h)*h/2 <=N
	// (2n+1)*h - h^2 <=2N
	//h^2 - (2n+1)h >=-2N
	//h^2 - (2n+1)h + 2N >=0

	//���߿�������
	int a = 1, b = -(1 + 2 * n), c = 2 * iPos;
	float fDelta =(float)(b * b - 4 * c);
	float sqrt_b_sqr_4ac = sqrt(fDelta);
	if (fDelta < 0)
	{
		printf("��������");
		return 0;
	}

	float h[2] = { (-b + sqrt_b_sqr_4ac) / 2.f,	//������
		(-b - sqrt_b_sqr_4ac) / 2.f };

	//��������������y1[0]Ϊ�����y1[1]ΪС��
	//��������
	int h1 = (int)floor(h[1]);
	if (h1 > n || h1 < 0)
	{
		h1 =(int) std::floor(h[1]);
		if (h1 > n || h1 < 0)
			return 0;
	}

	//yΪʲô���ڸߣ���Ϊ����һ�п�ʼ
	int x, x_Start = h1, y = h1;
	int iLine_Start = (n + (n - y + 1)) * y / 2;
	if (iLine_Start == iPos)
		x = y;
	else
		x = iPos - iLine_Start + x_Start;
	*px = x, * py = y;
	return 1;
}
int iGet_Upper_Pos(int x, int y, int n)
{//����һ������Ľ���n������������(x,y)��Ӧ��λ��
	int iPos = (n + (n - y + 1)) * y / 2;
	return iPos + x - y;
}
/**********************һЩ������λ�û���*******************/
