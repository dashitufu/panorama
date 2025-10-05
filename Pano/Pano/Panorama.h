#pragma once

typedef struct Image_Pair {
    unsigned short m_iImage_A, m_iImage_B;
    short roi[2][2];
    int* m_pSub_dx_A, * m_pSub_dy_A,
        * m_pSub_dx_B, * m_pSub_dy_B;
    Image m_oImage_A_GPU, m_oImage_B_GPU;
    Image m_oImage_A, m_oImage_B;
}Image_Pair;

template <typename _T>struct K_Rinv {
    _T M[3 * 3];
};

typedef struct Overlap_Pair {
    unsigned short m_Pair[2];   //�����е�һ��Block_ID;
    short roi[2][2];            //��ʼλ�ã���ߣ��Ժ�����������
}Overlap_Pair;

template<typename _T> struct Sphere_Projector {
    _T scale;       //�˴�������� fx* seam_work_aspect_
    _T K[3 * 3];
    _T rinv[3 * 3];
    _T r_kinv[3 * 3];
    union {
        _T k_rinv[3 * 3];
        K_Rinv<_T> oK_Rinv;
    };    
    _T t[3];
};
template<typename _T>struct Camera {
    _T K[3 * 3], R[3 * 3], t[3], H[3 * 3];
};
typedef struct Compensator {
    int nr_feed = 1;
    float fSimilarity_Threshold = 1;
    double* m_pGain;    //Ҫע��Gain�Ĺ��죬��һ��һ��Image��Block����
                        //image 0: [h0][w0] Image 1: [h1][w1] ...

    double** m_pGain_Map;   //һͼһ����
}Compensator;
template<typename _T>struct Match_Item {
    short m_iImage_A, m_iImage_B;
    unsigned short m_iMatch_Count;	//�㹻����
    _T(*m_pPoint_1)[2],			//�ĳ�������ʽ�����ձ���
        (*m_pPoint_2)[2];
};
typedef struct Image_Size_In_Block {
    unsigned char m_Block_Per_Image[2]; //����ٸ����߶��ٸ�Block
    unsigned short m_iBlock_Start;      //һ��Warp_Image����Block���ĸ���ſ�ʼ
}Image_Size_In_Block;
template<typename _T>struct Stitch {
    //һЩ����
    _T seam_work_aspect = 0.57735026918962573;
    _T WORK_MEGAPIX = 0.3;
    _T registr_resol_ = WORK_MEGAPIX;
    _T seam_est_resol_ = 0.1;
    _T work_scale = 0.79056941504209488;
    _T seam_scale;
    _T warped_image_scale;
    _T compose_scale = 1.0;

    //ָʾwork_scale�Ƿ��Ѿ�����
    unsigned char m_bWork_Scale_Set = 0;
    unsigned char m_bSeam_Scale_set = 0;

    //����Size
    unsigned short m_Source_Size[2];
    unsigned short m_Seam_Size[2];

    Camera<_T>* m_pCamera;            //�ж�����ͼ���ж��ٸ����
    Match_Item<_T>* m_pImage_Match; //����ͼƥ��
    Image* m_pImage_Source;
    Image* m_pSeam_Est;             //�о��Ƿ�Ӵ�
    Image* m_pMask;
    Image* m_pImage_Warp;         //Ӧ��һ��image��һ��warp
    Image* m_pMasks_Warped;

    //��Image��ͷ������GPU��
    Image* m_pImage_Source_Header_GPU;  
    Image* m_pSeam_Est_Header_GPU; 
    Image* m_pMask_Header_GPU;;
    Image* m_pImage_Warp_Header_GPU;
    Image* m_pMask_Warp_Header_GPU;
    //unsigned char(*m_pBlock_Per_Image_GPU)[2];
    Image_Size_In_Block* m_pBlock_Per_Image_GPU;

    Image* m_pBlock_Image_Header_GPU;         //ÿһ��Image�ֳ����ٸ�Title
    unsigned char* m_pBlock_Image_Data_GPU;       //�����������

    int (*m_pCorner)[2][2];
    int (*m_pCorner_GPU)[2][2];

    int (*m_pSize)[2];              //װƴ��ʱ��ͼ���С��
    int m_iImage_Count;

    //#define Max_Image_Count  4

    Image_Size_In_Block *m_pBlock_Per_Image;     //ÿ��Images_warped����Ϊ���ٸ�32*32��
    
    int (*m_pBlock_Corner)[2];
    int (*m_pBlock_Corner_GPU)[2];

    Image* m_pBlock_Image;
    Image* m_pBlock_Masks;          //ÿ��Image�ֳ����ٸ�mask
    int m_iBlock_Count;             //�������е�Image_Block����

    Compensator m_oComp;

    _T *m_pKer_GPU;         //f = { 0.25, 0.5, 0.25,0,0 };
    int* m_pdx_dy_GPU;      //dx,dy���ٵĿռ�ָ��
    int** m_dx_dy;          //�����ͼ��dx,dyָ��
};

typedef struct Graph_Vertex {
    int m_iID;
    //float m_fSource_W, m_fSink_W;
    float weight;

    int parent;
    int first;
    int ts;
    int dist;

    unsigned char t;
    Graph_Vertex* next;
}Graph_Vertex;

typedef struct Graph_Edge {
    //int v0, v1; //�������������ֵ��
    //float w, revw;    
    float weight;
    int dst, next;
}Graph_Edge;

typedef struct GCGraph {    //ͼ��ָGraph Cut
    Graph_Vertex* m_Vertex;
    Graph_Edge* m_Edge;
    int m_iMax_Vertex_Count; //��������
    int m_iMax_Edge_Count;  //������
    int m_iCur_Edge;    //��ǰ�ӵ�������
    float flow;
}GCGraph;

typedef struct Blender {

    int dst_roi_final[2][2];
    int dst_roi[2][2];      //����һ��roi���� [0][0], [0][1] ���Ͻ����꣬ [1][0],[1][1]Ϊ���
    int m_iNum_Band;       //ȷ����Block�Ĵ�С��

    Image m_oMask;          //ͬ��ͼ��С
    //short (*m_pDst)[3];     //һ�Ŵ�ͼ��16λ��3ͨ��

    struct {
        Image::Part_1 dst_pyr_laplace[6],
            dst_band_weights[6];
    };

    struct {
        Image::Part_1* dst_pyr_laplace_Header_GPU,
            * dst_band_weights_Header_GPU;
    };

    Image m_oResult;

    unsigned char* m_pBuffer;
}Blender;

template<typename _T>void Init_Stitch(Stitch<_T>* poStitch, int iWidth, int iHeight);
template<typename _T>void Free_Stitch(Stitch<_T>* poStitch);

template<typename _T>void Resize_Seam_Image(Stitch<_T>* poStitch, Light_Ptr oPtr = {});

template<typename _T>void Warp_2(Image oImage, Image oMask, _T K[3 * 3], _T R[3 * 3], _T fScale,
    Image* poImage_Warped, Image* poMask_Warped, Interpolation_Flag iImage_Inter_Type, Interpolation_Flag iMask_Inter_Type,
    int Corner[2][2],Border_Type iImage_Border_Type = BORDER_REFLECT, Border_Type iMask_Border_Type = BORDER_CONSTANT, Point_Cloud<float>* poPC = NULL);

template<typename _T>void Warp_3(Image oImage, _T K[3 * 3], _T R[3 * 3], _T fScale,
    Image* poImage_Warped, Interpolation_Flag iImage_Inter_Type, Interpolation_Flag iMask_Inter_Type,
    int Corner[2][2], Border_Type iImage_Border_Type, Border_Type iMask_Border_Type, Point_Cloud<float>* poPC = NULL);

template<typename _T>void Feed(Stitch<_T>* poStitch);
template<typename _T>void Find(Stitch<_T>* poStitch);

template<typename _T> void Build_Map(int w, int h, _T K[], _T R[], _T fScale,
    /*_T** ppx_Map, _T** ppy_Map,*/
    int Dest_roi[2][2], Sphere_Projector<_T>* poProjector, Point_Cloud<float>* poPC = NULL);

template<typename _T>void Re_Map_3_GPU(Image oImage, Image oImage_Warped,
    int roi_x, int roi_y, Sphere_Projector<_T>* poProjector,
    Interpolation_Flag iImage_Inter_Type, Interpolation_Flag iMask_Inter_Type,
    Border_Type iImage_Border_Type, Border_Type iMask_Border_Type);

template<typename _T>void Block_Compensate(Stitch<_T>* poStitch,
    Image Image_Warp[], Image Image_Warp_Header_GPU[]);

void Dilate_GPU(Image Warp[], Image Warp_Header_GPU[], int iImage_Count);
template<typename _T>void Resize_Bitwise_And(Stitch<_T>* poStitch,
    Image Blend_Warp[], Image Blend_Warp_Header_GPU[]);

