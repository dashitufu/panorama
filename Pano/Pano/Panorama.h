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
    unsigned short m_Pair[2];   //交集中的一对Block_ID;
    short roi[2][2];            //开始位置，宽高，以后再慢慢整理
}Overlap_Pair;

template<typename _T> struct Sphere_Projector {
    _T scale;       //此处让相机的 fx* seam_work_aspect_
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
    double* m_pGain;    //要注意Gain的构造，是一个一个Image的Block构成
                        //image 0: [h0][w0] Image 1: [h1][w1] ...

    double** m_pGain_Map;   //一图一增益
}Compensator;
template<typename _T>struct Match_Item {
    short m_iImage_A, m_iImage_B;
    unsigned short m_iMatch_Count;	//足够大了
    _T(*m_pPoint_1)[2],			//改成这种形式更具普遍性
        (*m_pPoint_2)[2];
};
typedef struct Image_Size_In_Block {
    unsigned char m_Block_Per_Image[2]; //横多少个，高多少个Block
    unsigned short m_iBlock_Start;      //一个Warp_Image，其Block从哪个标号开始
}Image_Size_In_Block;
template<typename _T>struct Stitch {
    //一些常量
    _T seam_work_aspect = 0.57735026918962573;
    _T WORK_MEGAPIX = 0.3;
    _T registr_resol_ = WORK_MEGAPIX;
    _T seam_est_resol_ = 0.1;
    _T work_scale = 0.79056941504209488;
    _T seam_scale;
    _T warped_image_scale;
    _T compose_scale = 1.0;

    //指示work_scale是否已经设置
    unsigned char m_bWork_Scale_Set = 0;
    unsigned char m_bSeam_Scale_set = 0;

    //各种Size
    unsigned short m_Source_Size[2];
    unsigned short m_Seam_Size[2];

    Camera<_T>* m_pCamera;            //有多少张图就有多少个相机
    Match_Item<_T>* m_pImage_Match; //所有图匹配
    Image* m_pImage_Source;
    Image* m_pSeam_Est;             //感觉是缝接处
    Image* m_pMask;
    Image* m_pImage_Warp;         //应该一个image就一组warp
    Image* m_pMasks_Warped;

    //连Image的头都放在GPU内
    Image* m_pImage_Source_Header_GPU;  
    Image* m_pSeam_Est_Header_GPU; 
    Image* m_pMask_Header_GPU;;
    Image* m_pImage_Warp_Header_GPU;
    Image* m_pMask_Warp_Header_GPU;
    //unsigned char(*m_pBlock_Per_Image_GPU)[2];
    Image_Size_In_Block* m_pBlock_Per_Image_GPU;

    Image* m_pBlock_Image_Header_GPU;         //每一个Image分出多少个Title
    unsigned char* m_pBlock_Image_Data_GPU;       //放上面得内容

    int (*m_pCorner)[2][2];
    int (*m_pCorner_GPU)[2][2];

    int (*m_pSize)[2];              //装拼接时的图像大小？
    int m_iImage_Count;

    //#define Max_Image_Count  4

    Image_Size_In_Block *m_pBlock_Per_Image;     //每个Images_warped划分为多少个32*32块
    
    int (*m_pBlock_Corner)[2];
    int (*m_pBlock_Corner_GPU)[2];

    Image* m_pBlock_Image;
    Image* m_pBlock_Masks;          //每个Image分出多少个mask
    int m_iBlock_Count;             //以上所有的Image_Block数量

    Compensator m_oComp;

    _T *m_pKer_GPU;         //f = { 0.25, 0.5, 0.25,0,0 };
    int* m_pdx_dy_GPU;      //dx,dy开辟的空间指针
    int** m_dx_dy;          //具体各图的dx,dy指针
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
    //int v0, v1; //两个顶点的索引值？
    //float w, revw;    
    float weight;
    int dst, next;
}Graph_Edge;

typedef struct GCGraph {    //图像分割，Graph Cut
    Graph_Vertex* m_Vertex;
    Graph_Edge* m_Edge;
    int m_iMax_Vertex_Count; //顶点数量
    int m_iMax_Edge_Count;  //边数量
    int m_iCur_Edge;    //当前加到哪条边
    float flow;
}GCGraph;

typedef struct Blender {

    int dst_roi_final[2][2];
    int dst_roi[2][2];      //还是一个roi区域， [0][0], [0][1] 左上角坐标， [1][0],[1][1]为宽高
    int m_iNum_Band;       //确定了Block的大小？

    Image m_oMask;          //同大图大小
    //short (*m_pDst)[3];     //一张大图？16位，3通道

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

