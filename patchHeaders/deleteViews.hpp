#pragma once

//根据某种规则删除哈希表里面的内容
template<class Views,
         class JudgeDeleteFunctor>
class DeleteViews
{
public:

    //用于判断一个view是否需要被删除
    JudgeDeleteFunctor needDelete;

    //遍历views然后按照指定的规则删除
    void operator()(Views& dstViews)
    {
        //遍历每个view
        for(auto iter=dstViews.begin();iter!=dstViews.end();)
        {
            //判断是否需要删除
            if(needDelete(*iter))
            {
                iter=dstViews.erase(iter);
            }
            else
            {
                ++iter;
            }
        }
    }
};
