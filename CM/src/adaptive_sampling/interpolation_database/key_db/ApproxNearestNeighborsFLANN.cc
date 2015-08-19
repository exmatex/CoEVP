#include "ApproxNearestNeighborsFLANN.h"

#ifdef FLANN
void ApproxNearestNeighborsFLANN::insert(std::vector<uint128_t> const& keys) {
        keymap.insert(keymap.end(), keys.begin(), keys.end());
        
        std::vector<double> raw_pts(keys.size() * dim);
        for (int i = 0; i < keys.size(); i++) {
            double *d = database_pull_key(keys[i]);
            for (int j = 0; j < dim; j++) {
                raw_pts[i*dim + j] = d[j];
            }
        }
        flann::Matrix<double> pts(raw_pts.data(), keys.size(), dim);
        if (is_empty) {
            flann_index.buildIndex(pts);
            is_empty = false;
        }
        else {
            flann_index.addPoints(pts);
        }
    }
    
    
void ApproxNearestNeighborsFLANN::knn_helper(std::vector<double> const& x, int k, int n_checks, std::vector<uint128_t> &keys, std::vector<double> &dists) {
        if (is_empty) {
            keys.resize(0);
            dists.resize(0);
        }
        else {
            flann::Matrix<double> query((double *)x.data(), 1, dim);
            
            std::vector<int> raw_indices(k);
            flann::Matrix<int> indices(raw_indices.data(), 1, k);
            
            std::vector<double> raw_dists(k);
            flann::Matrix<double> mdists(raw_dists.data(), 1, k);
            
            flann_index.knnSearch(query, indices, mdists, k, flann::SearchParams(n_checks));
            
            keys.resize(k);
            dists.resize(k);
            for (int i = 0; i < k; i++) {
                keys[i]  = keymap[(size_t)indices[0][i]];
                dists[i] = mdists[0][i];
            }
        }
    }
    
void ApproxNearestNeighborsFLANN::knn(std::vector<double> const& x, int k, std::vector<uint128_t> &keys, std::vector<double> &dists) {
        knn_helper(x, k, n_checks_default, keys, dists);
    }
#endif
