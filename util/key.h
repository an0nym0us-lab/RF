#ifndef KEY
#define KEY
#include <vector>
#include <string>
#include <stdint.h>

struct Slice {
    std::string str;
    uint32_t ipaddr;
    double cost;
    //std::vector<uint8_t> hash_map_;
    bool operator< (const Slice &s) const{
        return cost > s.cost;
    }
    bool operator< (const Slice *s) const{
        return cost > s->cost;
    }
};

#endif