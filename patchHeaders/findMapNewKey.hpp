#pragma once

//对于一个哈希表，随便寻找一个它还没加入过的键值
//键值默认都是unsigned,最好是整型
template<class Map,class Key=unsigned>
class FindMapNewKey
{
public:
    //对输入的map类型寻找新的键值
    Key operator()(const Map& targetMap)
    {
        //从size的位置开始找
        Key tempKey=targetMap.size();
        while(targetMap.count(tempKey))
        {
            ++tempKey;
        }
        //找到不存在的键值了就返回
        return tempKey;
    }
};
