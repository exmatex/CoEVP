#ifndef included_ApproxNearestNeighbors_h
#define included_ApproxNearestNeighbors_h

#include <vector>
#define uint128_t unsigned __int128
#define uint128_t_undefined 0xffffffffffffffffffffffffffffffff
#define id_undefined -1

class ApproxNearestNeighbors {
public:
    virtual int insert(std::vector<double>& point, uint128_t& key) = 0;
    virtual void insert(std::vector<uint128_t> const& keys) = 0;
    virtual void remove(int id) = 0;
    virtual void knn(std::vector<double> const& x, int k, std::vector<int> &ids, std::vector<uint128_t> &keys, std::vector<double> &dists) = 0;
    virtual uint128_t getKey(int id) = 0;
};
#endif
