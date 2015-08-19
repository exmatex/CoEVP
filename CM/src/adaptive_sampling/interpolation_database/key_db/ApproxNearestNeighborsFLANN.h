#ifndef included_ApproxNearestNeighborsFLANN_h
#define included_ApproxNearestNeighborsFLANN_h
#include <functional>
#include "ApproxNearestNeighbors.h"

#ifdef FLANN
#include <flann/flann.hpp>
/*
 Generic KNN interface
*/

// TODO: remove dynamic memory allocations in insert() and knn_helper() for optimal performance
class ApproxNearestNeighborsFLANN : public ApproxNearestNeighbors {
public:
    int dim;
    std::function<double *(uint128_t)> database_pull_key;
    int n_trees, n_checks_default;
    
    flann::Index<flann::L2<double>> flann_index;
    bool is_empty = true;
    std::vector<uint128_t> keymap;
    
    ApproxNearestNeighborsFLANN(int dim, std::function<double *(uint128_t)> database_pull_key,
                                int n_trees, int n_checks_default):
        dim(dim), database_pull_key(database_pull_key), n_trees(n_trees), n_checks_default(n_checks_default),
        flann_index(flann::KDTreeIndexParams(n_trees))
    {
    }
    
    void insert(std::vector<uint128_t> const& keys);
    void knn_helper(std::vector<double> const& x, int k, int n_checks, std::vector<uint128_t> &keys, std::vector<double> &dists);
    void knn(std::vector<double> const& x, int k, std::vector<uint128_t> &keys, std::vector<double> &dists);
};
#endif
#endif
