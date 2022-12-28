#ifndef _SWITCH_HPP_
#define _SWITCH_HPP_
#include<iostream>
//是否使用多线程,这个开关只负责赵志豪自己写的类和函数
#define USE_OMP
#define SIN_THRE 0.2 //正弦值小于这个数字的不要
#define COL_THRE 0.95 //两个颜色的夹角正弦小于这个距离被视为不同的颜色
#define CORR_THRE 0.99f //若干颜色的互相关性的阈值
//当使用z平面策略来做的时候，对z直方图的区分度
//目前引入了高度的图片分布图，目前使用的是uchar来记录高度分布
//因此这里使用的数值应该尽量大于1/255，如果设置小于这个数值的精度值，那还需要修改load_img的代码
#define STEP_L1 (1.f/200) //一阶精度，在z上找坐标的时候每次步长是z范围的1%
//20211108是否使用欧氏距离，如果使用欧氏距离，则关于NCC的距离全部失效
//不过NCC的阈值仍然是有效的
//#define USE_EUCLID_DIS
//注意，当使用聚类算法的时候，下面的NCC阈值之类的仍然有效
//只不过使用聚类算法的时候，这里的阈值表示的是聚类结果的最大值与全部图片的比值
#define NCC_THRE 0.7 //两个图的颜色互相关性的阈值，也就是论文里面的NCC
//表示总共需要迭代的次数，这个次数会按照比例分配给每人高度阶层
#define FIX_MASK_TIMES 1000//按照一个固定的大小循环若干次，而不是依次递减的窗口大小
//将每个点具有可参考次数改为每个点具有的可被优化的次数
//每个点可被优化的次数是有限的，并且不会更新
//这个选项和可被参考次数这个选项是冲突的，不能同时打开
//#define USE_MAX_OPT_TIME
#ifndef USE_MAX_OPT_TIME
//是否限制每个点的被参考次数
#define USE_REFINE_TIME
#endif
//每个先验点可用次数的初始值
#define REF_TIME_INIT 1
//什么时候有必要开始判断结束循环的条件
#define CONSIDER_BREAK_TIME 400
//什么时候初始化所有的可用次数
inline bool NEED_INIT_REF(int id)
{
    if(id>=CONSIDER_BREAK_TIME)
    {
        return id%10==0;
    }
    //默认不更新可用次数
    return false;
}
//当global_help生效的时候，下面的这个宏是用来判断global_help具体在哪一次迭代的时候生效的
//下面是一个已经弃用的选项
#define SLIDER_BG_TIME 999 //这个次数之后，开始稠密滑动，即两次滑动窗口是挨着的,fix_mask_times有效的时候它才有效
//与迭代次数相关联的阈值判定,每次迭代的时候阈值都不一样
//这里输入的是迭代次数
inline double THRE_ID(int id)
{
    //判断如果是最后一个周期
    if(id>=CONSIDER_BREAK_TIME) return NCC_THRE-(id-CONSIDER_BREAK_TIME)*0.02;
    //每20个节点是一个小周期
    return NCC_THRE;
}
//可参考的阈值与真实阈值的差值
#define REF_THRE_DIS 0.4

//每次迭代的时候允许用到的最小的z阶层，用于保证高处的场景被优先重建
//不用的时候把它写成0就可以
inline int ALLOW_MIN_Z(int id)
{
    //在第1个周期，按照序号降低
    return 199.f-id*0.5;
}

//程序运行过程中输出的log文件信息的名称
//20211030，这是为了研究纯色区域传播较慢的问题特别设立的接口
inline std::string logFileName(int domX,int domY)
{
    //返回一个字符串
    return "/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/reconstruct/"
            +std::to_string(domX)+"-"+std::to_string(domY)+".txt";
}

//20211126 使用外部存储的方式存储直线信息
//#define SAVE_EXTERNAL_PROJECT_LINE
//每个DOM像素对应的二进制投影直线的文件路径
inline std::string projectBinFileName(int unitId)
{
    return "/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/projectLines/line_"
            +std::to_string(unitId)+".bin";
}

//20211029为了研究纯色区域传播较慢的问题，添加的一个用于debug的条件判断
inline bool debugStop(int xValue,int yValue)
{
    return false;
    if(239>=xValue && 239-xValue<=3 && 375>=yValue && 375-yValue<=3)
        return true;
    return false;
}

#define PRIOR_NCC 0.5//用于给先验点写入，让它一开始就成为一个可以使用的点
#define NCC_MIN_IMG 2 //计算某个位置的颜色时，至少需要有几个能看到该点的图片
#define NCC_RM_AVG //计算NCC分数的时候是否减去平均颜色
//20210930 为了让NCC的分数更合理而设置的一个选项，以前是减去每个面片各自的平均值
//这个选项会把算法改成减去平均面片的平均值
//这个选项与上面的选项不冲突，如果同时使用，会减去平均面片的平均颜色的通道平均颜色
//#define NCC_CENTER_AVG
#define NCC_NORM_1 //是否把颜色向量归一化成方差为1
//使用NCC分数的时候是否根据点到光心的连线与Z轴的夹角正弦来作为权值
#define USE_NCC_WEIGHT
//使用k均值算法的时候，这个选项无效
#define NCC_WRITE_AVG_COLOR //使用NCC算法写入颜色的时候，使用所有颜色的平均值来写入
#ifdef NCC_WRITE_AVG_COLOR
//只有写入平均颜色的时候下面这个选项才生效
//20211103
//选取中相机里面最接近投影点最接近主点的那个面片作为所谓的平均颜色
//这个选项经过验证具有特别良好的效果，基本上以后永远都会开着这个选项
#define USE_PRI_AS_AVG
#endif
//20211122 是否每个位置的颜色面片都需要使用权值参与运算
#ifdef USE_NCC_WEIGHT
//#define USE_EACH_NCC_WEIGHT
#endif
//NCC分数是否取绝对值
//#define NCC_ABS
//#define USE_NCC_MEDIUM //是否使用NCC数据的中位数作为比较依据，否则使用平均值来比较
#define MIN_CORR 4 //计算颜色的互相关性的时候最少需要有几个图片
#define MIN_MASK_PT 8 //在dom上遍历的时候，每次mask里面选取到的有效点不能少于这个数
#define MASK_STEP 1 //从dom图里面选取像素的时候使用的步长
//新的向量等于原始的单位向量加上0.4*新的随机单位向量
#define NORM_STEP 0.3 //法向量的更新幅度，原始的向量是单位向量
//20210926 法向量里面Z值的最小分量，为了保证法向量变换的过程中不要变得太离谱
#define MIN_Z_WEIGHT 0.5
#define NORM_OPT_TIMES 7 //法向量的最大更新次数,超过这个次数之后如果还没找到就返回了
#define Z_OPT_TIMES 10 //z值的随机更新次数,这个参数已经弃用了
#define RAND_Z 0.3 //对z做随机更新的时候使用的长度比例
#define RADI_L2 7 //在一阶精度的基础上，再寻找几次
#define SEARCH_UP_RATE 0.1 //比如说稀疏点云的高度范围是0~1,那么做稠密的时候搜索范围是-0.2~1.2
//遍历DOM的时候使用的矩形面片的大小
inline int MASK_SIZE(int iterId)
{
    //从3～5循环变化
    return 3;
}
#define SAVE_EACH_TIME //是否每运行完一次就保留一次运行的dom结果？
//图片路径的存储方式，传入的必须是字符串形式的id，不然可能会出问题
//对于从上到下重建的过程，这里传入的是z的阶层，而不是迭代的次数
#define SAVE_IMG_PATH(id) "/media/cvlab/data/workSpace/mainProject/topViewConstruct/workspace/reconstruct/"+id+".bmp"
#define NEED_SAVE(id) (id%20==0)
//#define GLOBAL_HELP //是否使用德劳内三角化提供的z值做辅助判断，当某些点距离先验点实在太远的时候使用
#define BUFFER_SIZE 0.4 //仅通过法向量来做优化的预计不成功的比例
#define USE_ONLY_MID_CAM //仅使用中间位置的图片
//传入一个图片的标号，判断这是否为一个中间位置的图片
inline bool MID_ID_JUDGE(int id)
{
    return id/100000==2;
}
//20211124 为了适应载入大场景DOM需求，把程序流程更改成动态载入图片
#define DYNAMIC_LOAD_IMAGE
//#define USE_ONLY_IMG //dom优先遍历的时候，仅使用投影最长的图片信息,全都遍历实在是太慢了 已经弃用了
#define DEBUG_PRINT //print调试的开关
//#define NEED_F //是否需要计算F矩阵
//使用的时候会在home目录检测是否存在"empty0"这个文件，如果存在才会保存图片
//#define SAVE_PROJ_PT //寻找颜色一致性的时候，是否保留投影点在图片上的位置，并把投影后的图片保存下来
//#define SET_COLOR_NORM 150.f //把所有颜色的模长都弄成这个数值，保持亮度一致性,不用可以注释掉
//聚类算法的各种参数
#define CLASS_NUM 3 //总共分类的个数
#define KMEAN_ITER_TIME 50 //聚类算法的迭代次数
//#define USE_KMEANS //在有歧义颜色的时候是否使用聚类算法，使用聚类算法就不用NCC了
#define USE_Z_PLANE //每个面片的初始位置使用平面，而不是根据每个高度使用斜面
//#define FIX_SAMPLE //固定采样的意思是，虽然掩膜很大，但只拿角的和各边中点以及整个正方开的中心点
#define OPT_REFINE //当待优化的点附近出现了已经优化过的点的时候，还要不要再优化它一下试试
//20210928 是否优化先验点,如果优化先验点的话，那稀疏点也有可能被改变
#define PRIOR_REFINE
//2021118 是否优化优化次数已经用完了的DOM像素
//#define REFINE_REFINE_ZERO
//#define USE_BRIGHT_CONSTRAIN //是否使用面片之间的亮度约束
//注意，use_img_z_constrain和上面的use_bright_constrain如果同时使用可能会出现问题
//#define USE_IMG_Z_CONSTRAIN //如果之前有一个更高的z点投影过图片上的某个像素，是否应该拒绝这个像素再受其它低位点的投影
#define SAVE_Z_MAP //使用Z值约束的情况下是否保存Z值约束形成的图片
//20210926 是否使用法向量传播的开关
//#define USE_NORM_PROG
//是否对倾斜遮挡的情况做判别 20211130
#ifdef USE_PRI_AS_AVG//只有采用选取主点坐标的方案，下面这个选项才有意义
//#define USE_OCCLUDE_DETECT
#endif
//差别遮挡时使用的步长
#define OCCLUDE_DETECT_SEG_NUM 5
//根据传入的移动步长与DOM像素的比值，判断是否需要在这个位置判断步长
inline bool needJudgeOcclude(double stepRate)
{
    return stepRate>1;
}
//判断某个周期是否需要判断遮挡问题
inline bool isIterateDetectOcclude(int iterTimes)
{
    return iterTimes>=100;
}
//是否使用复杂的遮挡判断方案
#define USE_COMPLEX_OCCLUDE_DETECTION
//检测遮挡的过程中，如果遇到的像素比当前的中心像素高，就直接判定为遮挡
//#define OCCLUDE_DETECT_ANY_HIGHER
//判断第2个高度是否比第1个高度更高
//这里面传入的是减去了最低值的高度
inline bool isDOMHigher(double height1,double height2)
{
    return height2-height1>(-25);
}
//是否需要在最后单纯做一次遮挡判断，用来处理颜色遮挡问题
//#define END_WITH_OCCLUDE_DETECTION
//2021-9-18
#define CONSTRAIN_Z_LINE //把上面的USE_IMG_Z_CONSTRAIN打开这个才有效，如果某个点被约束成功，那么下面所有的点都记录一下
//如果遇到了Z值约束，则仅仅是放弃把它作为中心像素，但仍然使用
#define BAN_CENTER_FOR_CONSTRAIN
//存储Z阶层时对Z阶层数据做的变换
//为了让相邻的Z值有一个接近的量化层次
inline int Z_LEVEL_TRANS(int zLevel)
{
    //如果小于0就什么都不做
    return zLevel;
}
//访问Z值图的时候做的变换，有时候可能需要把自身加上一个数字
//与上面的BAN_CENTER_FOR_CONSTRAIN配套使用
inline int Z_LEVEL_VISIT_TRANS(int zLevel)
{
    //如果小于0就什么都不做
    return zLevel+6;
}
#define BRIGHT_MAX_DIS 50 //当需要参考颜色的亮度的时候，亮度和主亮度的差别不能超过多少
//#define OVERLAP_RATE(imgNum) 4.f //图片的重叠度，用于判断是否需要把图片的分辨率降低的一个参数
//#define NCC_EVO //是否对NCC分数做与图片个数相关的根号处理


//下面是各种报错的错误码
#define ERROR_PLANE_PT_LESS_3 -1 //当计算一个平面方程的时候，传入的点少于三个的时候，返回这个值
#define ERROR_PLANE_COLLINEATION -2 //当计算一个平面方程的时候，如果传入的三个点是共线的，返回这个错误

//构造二维范围数据Range2D的时候，如果传入的是一个空的向量的话，会返回这个错误
#define ERROR_EMPTY_PT_LIST_MAKING_RANGE -3

//做分辨率转换的时候，像素长度为0是致命的
#define ERROR_ZERO_PIXEL_LEN -4

//画三角形的时候，如果检测到获取到的点不是三个，报这个错误
#define ERROR_TRI_PT_NOT_3 -5

//初始化一个图片的时候，使用了包含0项的size
#define ERROR_ZERO_IMG_SIZE -6

//法向量的模长是0
#define ERROR_ZERO_NORMAL -7

//传入两个数据，计算平面的另一个数据的时候，当有多解的时候返回这个错误
#define ERROR_MANY_VALUE -8

//计算平面数据的时候，如果传入的数据维度不正常，返回这个错误
#define ERROR_PLANE_DIM_OUT -9

//将一个vector转换为eigen的时候，如果指定的通道数与向量原来的维度不一样的话，报这个错
#define ERROR_EIGEN_VEC_DIFF_DIM -10

//当一个矩形面片里面的点个数太少的时候，就是出现这个错误
#define ERROR_MASK_PT_TOO_LESS -11

//DOM信息图的结果访问超过了索引位置
#define ERROR_DOMINFO_OUTRANGE -12

//图片访问范围的越界
#define ERROR_IMG_OUT_RANGE -13

//VectorInterface里面的两个操作单元的size不一致的时候会出这个错
#define ERROR_VEC_SIZE_DIFF -14

//计算NCC分数的时候，如果参与计算的量没有事先减去平均值，会throw这个东西
#define ERROR_NOT_AVG_ZERO -15

//对z做优化的时候需要把点传回面片，然后根据传回的点列表恢复面片里面每个点的坐标
//如果传回的面片size()异常，会返回这个错误
#define ERROR_POINT_DIFF_SIZE -16

//map里面不存在对应的索引但强行访问的时候出这个错误
#define ERROR_NO_MAP_IDX -17

//表示zRange的指针为空的情况
#define ERROR_EMPTY_ZRANGE -18

//强行访问了空的列表
#define ERROR_VISIT_EMPTY_LIST -19

//颜色序列的长度不是3的倍数，那说明这有问题
#define ERROR_NOT_3_TIMES -20

//访问vector的时候的访问越界
#define ERROR_VEC_OUT_RANGE -21

//获取平面法向量的时候取到了零向量
#define ERROR_GET_ZERO_NORM -22

//opencv读取图片失败的时候返回的错误
#define ERROR_IMREAD_FAIL -23

//如果z的阶层高度没有初始化过，会有这个错误
#define ERROR_NO_INIT_ZRANGE -24

//z优化器的个数和投影点的个数不同出现的问题
#define ERROR_PT_ZOPT_DIFF_SIZE -25

//对两个颜色面片的数据相加的时候，如果相加的两个颜色面片里面的颜色个数不一样，会出这个错
#define ERROR_ADD_PATCH_DIFF_SIZE -26

//如果颜色面片的投影点信息是空的但强行使用它，会出这个错
#define ERROR_PROJ_EMPTY -27

//使用主点坐标的时候，如果不存在中相机，则throw这个东西
#define ERROR_NO_MID_FOR_PRI -28

//投影直线打开失败的情况下报的错
#define ERROR_FAIL_OPEN_PROJECT_BIN -29

#define ERROR_INIT_PROJECT_LINES -30

#endif
