#include "ApproxNearestNeighborsFLANN.h"

#ifdef FLANN

void
ApproxNearestNeighborsFLANN::insert(std::vector<uint128_t> const& keys)
{
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

    
int
ApproxNearestNeighborsFLANN::insert(std::vector<double>& point,
                                    uint128_t& key)
{
   assert(point.size() == dim);

   int id = keymap.size();

   keymap.push_back(key);

   double * data = new double[dim];
   for (int i=0; i<dim; ++i) {
      data[i] = point[i];
   }

   points.push_back(data);

   flann::Matrix<double> pts(data, 1, dim);
   if (is_empty) {
      flann_index.buildIndex(pts);
      is_empty = false;
   }
   else {
      flann_index.addPoints(pts,10000);
   }

   return id;
}


void
ApproxNearestNeighborsFLANN::knn_helper(std::vector<double> const& x,
                                        int k,
                                        int n_checks,
                                        std::vector<int> &ids,
                                        std::vector<uint128_t> &keys,
                                        std::vector<double> &dists)
{
   if (is_empty) {
      ids.resize(0);
      keys.resize(0);
      dists.resize(0);
   }
   else {
      flann::Matrix<double> query((double *)x.data(), 1, dim);
            
      std::vector< std::vector<int> > indices;
      std::vector< std::vector<double> > mdists;

      flann_index.knnSearch(query, indices, mdists, k, flann::SearchParams(n_checks));
            
      int num_neighbors_found = indices[0].size();

      ids.resize(num_neighbors_found);
      keys.resize(num_neighbors_found);
      dists.resize(num_neighbors_found);
      for (int i = 0; i < num_neighbors_found; i++) {
         ids[i]  = indices[0][i];
         keys[i] = keymap[(size_t)ids[i]];
         if (keys[i] == uint128_t_undefined) {
            std::cout << "FLANN returned deleted model " << ids[i] << std::endl;
            exit(1);
         }
         dists[i] = mdists[0][i];
      }
   }
}

    
void
ApproxNearestNeighborsFLANN::knn(std::vector<double> const& x,
                                 int k,
                                 std::vector<int> &ids,
                                 std::vector<uint128_t> &keys,
                                 std::vector<double> &dists)
{
   knn_helper(x, k, n_checks_default, ids, keys, dists);
}


uint128_t
ApproxNearestNeighborsFLANN::getKey(int id)
{
   assert(id >= 0 && id < keymap.size());

   return keymap[id];
}


void
ApproxNearestNeighborsFLANN::remove(int id)
{
   assert(id >= 0 && id < keymap.size() && id < points.size());

   flann_index.removePoint(id);

   keymap[id] = uint128_t_undefined;

   delete points[id];
}


#endif
