#include <iostream>
#include "stdio.h"
#include "Common.h"
#include "Matrix.h"
const int LABEL_NOT_MARK = 0xFFFFFFF;
typedef struct GC {	//Graph Cut
	typedef struct Vertex {	//未必需要这个机构
		int m_iTo_Start,	//指向流出顶点的边集开始位置
			m_iFrom_Start;	//指向流入顶点的边集开始位置
		unsigned char m_iTo_Count,	//流出边（也是顶点）的数量，
			m_iFrom_Count;	//流入边的数量
	}Vertex;

	typedef struct Edge {
		int m_iVertex_From,	m_iVertex_To;	//看看能否省略From，只记录To
		float c;		//容量
		float f;		//流量，限制： f<=c
	}Edge;

	typedef struct Edge_1 {
		int m_iVertex_ID;
		float c;
		float f;
	}Edge_1;

	typedef struct Label {
		int m_iFrom;		//哪个顶点出发给它标的号，如果多个入点，按照图广度搜索入队先者进行处理
							//负数表示反向边
							// 0xFFFFFFF 表示未标号
		float a;			//可改进量，理论上初始化时给最大值
	}Label;

	Edge* m_pEdge;
	Edge_1* m_pTo, * m_pFrom;
	Label* m_pLabel;		//和Vertex一样多

	Vertex* m_pVertex;

	int m_iMax_Edge_Count, m_iMax_Vertex_Count;
	int m_iEdge_Count, m_iVertex_Count;
}GC;

typedef struct GC_Queue {
	int* m_pBuffer;	//还是装索引，此处为顶点索引
	int m_iHead;	//队头
	int m_iEnd;		//队尾
	int m_iCount;	//当前队列的元素个数
	int m_iBuffer_Size;	//当前队列可以装多少个元素，也是循环队列的模，一般来说，可以只装1/2的节点
}CN_Queue;

void Init_GC(GC *poGC,int iMax_Edge_Count,int iMax_Vertex_Count)
{
	GC oGC = {};
	int iSize;
	iSize = iMax_Edge_Count * sizeof(GC::Edge) + iMax_Vertex_Count * sizeof(GC::Vertex);
	oGC.m_iMax_Edge_Count = iMax_Edge_Count;
	oGC.m_iMax_Vertex_Count = iMax_Vertex_Count;

	oGC.m_pEdge = (GC::Edge*)pMalloc(iSize);
	memset(oGC.m_pEdge, 0, iSize);
	oGC.m_pVertex = (GC::Vertex*)(oGC.m_pEdge + iMax_Edge_Count);
	*poGC = oGC;
	return;
}

void Free_GC(GC* poGC)
{
	if (poGC->m_pEdge)
		Free(poGC->m_pEdge);
	if (poGC->m_pFrom)
		Free(poGC->m_pFrom);
	if (poGC->m_pLabel)
		Free(poGC->m_pLabel);
	*poGC = {};
}
void Add_Edge(GC* poGC, int iFrom, int iTo, float c)
{
	GC oGC = *poGC;
	int iMax = Max(iFrom, iTo)+1;
	oGC.m_iVertex_Count = Max(iMax, oGC.m_iVertex_Count);
	oGC.m_pEdge[oGC.m_iEdge_Count++] = { iFrom, iTo, c };

	oGC.m_pVertex[iFrom].m_iTo_Count++;
	oGC.m_pVertex[iTo].m_iFrom_Count++;

	*poGC = oGC;
	return;
}
void Set_GC_Index(GC* poGC)
{
	GC oGC = *poGC;
	oGC.m_pFrom = (GC::Edge_1*)pMalloc(oGC.m_iEdge_Count * 2 * sizeof(GC::Edge_1));
	memset(oGC.m_pFrom, 0, oGC.m_iEdge_Count * 2 * sizeof(GC::Edge_1));
	oGC.m_pTo = oGC.m_pFrom + oGC.m_iEdge_Count;
	for (int i = 1; i < oGC.m_iVertex_Count; i++)
	{//这步很难搞到GPU里
		oGC.m_pVertex[i].m_iTo_Start = oGC.m_pVertex[i-1].m_iTo_Start + oGC.m_pVertex[i - 1].m_iTo_Count;
		oGC.m_pVertex[i].m_iFrom_Start = oGC.m_pVertex[i - 1].m_iFrom_Start + oGC.m_pVertex[i - 1].m_iFrom_Count;

		//重置一下，后面接一下这个值
		oGC.m_pVertex[i - 1].m_iTo_Count = oGC.m_pVertex[i - 1].m_iFrom_Count = 0;
	}
	oGC.m_pVertex[5].m_iFrom_Count = oGC.m_pVertex[5].m_iTo_Count = 0;

	//剩下比较好办，可以GPU搞
	for (int i = 0; i < oGC.m_iEdge_Count; i++)
	{
		GC::Edge oEdge = oGC.m_pEdge[i];
		GC::Vertex oFrom = oGC.m_pVertex[oEdge.m_iVertex_From];
		GC::Vertex oTo = oGC.m_pVertex[oEdge.m_iVertex_To];
		//if(oGC.m_pVertex[5].m_iFrom_Count>0)
			//printf("Here");

		//int m_iVertex_ID;		float c;		float f;
		oGC.m_pTo[oFrom.m_iTo_Start + (oFrom.m_iTo_Count++)] = { oEdge.m_iVertex_To, oEdge.c, oEdge.f };
		oGC.m_pFrom[oTo.m_iFrom_Start + (oTo.m_iFrom_Count++)] = { oEdge.m_iVertex_From, oEdge.c,oEdge.f };
		//if (oGC.m_pTo[4].m_iVertex_ID == 1)
			//printf("Here");

		oGC.m_pVertex[oEdge.m_iVertex_From] = oFrom;
		oGC.m_pVertex[oEdge.m_iVertex_To] = oTo;
	}
	
	//此处可以做个内存整理

	*poGC = oGC;
	return;
}

void Disp_GC(GC oGC)
{
	printf("To:\n");
	for (int i = 0; i < oGC.m_iVertex_Count; i++)
	{
		GC::Vertex oVertex = oGC.m_pVertex[i];
		printf("Vertex:%d\n", i);
		for (int j = 0; j < oVertex.m_iTo_Count; j++)
		{
			GC::Edge_1 oEdge = oGC.m_pTo[oVertex.m_iTo_Start + j];
			printf("\t To:%d c:%f f:%f\n", oEdge.m_iVertex_ID, oEdge.c, oEdge.f);
		}
	}
	printf("From\n");
	for (int i = 0; i < oGC.m_iVertex_Count; i++)
	{
		GC::Vertex oVertex = oGC.m_pVertex[i];
		printf("Vertex:%d\n", i);
		for (int j = 0; j < oVertex.m_iFrom_Count; j++)
		{
			GC::Edge_1 oEdge = oGC.m_pFrom[oVertex.m_iFrom_Start + j];
			printf("\t From:%d c:%f f:%f\n", oEdge.m_iVertex_ID, oEdge.c, oEdge.f);
		}
	}
	return;
}
void Init_Label(GC* poGC)
{
	GC oGC = *poGC;
	oGC.m_pLabel = (GC::Label*)pMalloc(oGC.m_iVertex_Count * sizeof(GC::Label));
	for (int i = 0; i < oGC.m_iVertex_Count; i++)
		oGC.m_pLabel[i] = { LABEL_NOT_MARK, MAX_FLOAT};
	*poGC = oGC;
	return;
}
#define Init_Queue(oQueue, iBuffer_Size) \
{ \
	oQueue.m_iBuffer_Size = iBuffer_Size; \
	oQueue.m_iCount = oQueue.m_iEnd = oQueue.m_iHead = 0; \
	oQueue.m_pBuffer = (int*)pMalloc(sizeof(int) * iBuffer_Size); \
}
#define Free_Queue(oQueue) \
{ \
	Free(oQueue.m_pBuffer); \
}
#define In_Queue(oQueue, iIndex) \
{ \
	if ((oQueue).m_iEnd >= (oQueue).m_iBuffer_Size) \
	{ \
		printf("Error in In_Queue"); \
	}else \
	{ \
		(oQueue).m_pBuffer[(oQueue).m_iEnd++] = iIndex; \
		(oQueue).m_iCount++; \
	} \
}
#define Out_Queue(oQueue,iIndex) \
{ \
	iIndex=oQueue.m_pBuffer[oQueue.m_iHead++]; \
	oQueue.m_iCount--; \
}

void Add_To_Queue(GC_Queue*poQueue, GC oGC, int iVertex_ID)
{
	//对iVertex_ID后面的邻接顶点进行标号
	GC::Vertex oVertex = oGC.m_pVertex[iVertex_ID];
	GC::Edge_1 oEdge;
	
	for (int i = 0; i < oVertex.m_iTo_Count; i++)
	{
		oEdge = oGC.m_pTo[oVertex.m_iTo_Start+i];
		//对邻接顶点进行标号
		GC::Label *poLabel = &oGC.m_pLabel[oEdge.m_iVertex_ID];
		if (poLabel->m_iFrom != LABEL_NOT_MARK)
			continue;	//已经标号不必再搞

		poLabel->m_iFrom = iVertex_ID;
		poLabel->a = oEdge.c - oEdge.f;

		//将邻接顶点入队
		In_Queue(*poQueue, oEdge.m_iVertex_ID);
	}

	for (int i = 0; i < oVertex.m_iFrom_Count; i++)
	{
		oEdge = oGC.m_pTo[oVertex.m_iTo_Start + i];
		//对邻接顶点进行标号
		GC::Label* poLabel = &oGC.m_pLabel[oEdge.m_iVertex_ID];
		if (poLabel->m_iFrom != LABEL_NOT_MARK)
			continue;	//已经标号不必再搞

		poLabel->m_iFrom = -iVertex_ID;
		poLabel->a = oEdge.c - oEdge.f;

		//将邻接顶点入队
		In_Queue(*poQueue, oEdge.m_iVertex_ID);
	}

	//
	return;
}

void Disp_Label(GC oGC)
{//只显示标号
	printf("Label Status\n");
	for (int i = 0; i < oGC.m_iVertex_Count; i++)
	{
		GC::Label oLabel = oGC.m_pLabel[i];
		if (oLabel.m_iFrom == LABEL_NOT_MARK)
			continue;

		printf("Vertex:%d (%d,%f)\n ", i, oLabel.m_iFrom, oLabel.a);
	}

	return;
}
void Mark_Label(GC* poGC)
{//进行一次的标号
	GC_Queue oQueue;
	GC oGC = *poGC;
	int iVertex_ID;

	Init_Queue(oQueue, oGC.m_iVertex_Count);

	oGC.m_pLabel[0].m_iFrom = 0;

	In_Queue(oQueue, 0);
	while(oQueue.m_iCount)
	{
		Out_Queue(oQueue, iVertex_ID);
		Add_To_Queue(&oQueue,oGC, iVertex_ID);
		Disp_Label(oGC);
	}

	Free_Queue(oQueue);
	*poGC = oGC;
	return;
}
static void Test_1()
{
	GC oGC;
	Init_GC(&oGC, 100, 50);

	//e(0,1) = (8,0)
	Add_Edge(&oGC, 0, 1, 8);
	//e(0,2) = (4,0)
	Add_Edge(&oGC, 0, 2, 4);
	//e(1,3) = (2,0)
	Add_Edge(&oGC, 1, 3, 2);
	//e(1,4) = (2,0)
	Add_Edge(&oGC, 1, 4, 2);
	//e(2,1) = (4,0)
	Add_Edge(&oGC, 2, 1, 4);
	//e(2,3) = (1,0)
	Add_Edge(&oGC, 2, 3, 1);
	//e(2,4) = (4,0)
	Add_Edge(&oGC, 2, 4, 4);
	//e(3,4) = (6,0)
	Add_Edge(&oGC, 3, 4, 6);
	//e(3,5) = (9,0)
	Add_Edge(&oGC, 3, 5, 9);
	//e(4,5) = (7,0)
	Add_Edge(&oGC, 4, 5, 7);

	Set_GC_Index(&oGC);
	Init_Label(&oGC);

	Mark_Label(&oGC);

	Disp_GC(oGC);
	Free_GC(&oGC);
	return;
}

int main_1()
{
	Init_Env();
	Test_1();
	Free_Env();
	return 0;
}