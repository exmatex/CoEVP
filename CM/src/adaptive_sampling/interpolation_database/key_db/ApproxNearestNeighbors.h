#ifndef included_ApproxNearestNeighbors_h
#define included_ApproxNearestNeighbors_h

#include <vector>
#define uint128_t unsigned __int128

class ApproxNearestNeighbors {
public:
    virtual void insert(std::vector<uint128_t> const& keys) = 0;
    virtual void knn(std::vector<double> const& x, int k, std::vector<uint128_t> &keys, std::vector<double> &dists) = 0;
};
#endif
