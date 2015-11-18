#ifndef included_ApproxNearestNeighborsMTree_h
#define included_ApproxNearestNeighborsMTree_h

#include "ApproxNearestNeighbors.h"
//#include "mtreedb/MTree.h"
#include "MTree.h"

/*
  MTree KNN interface
*/

class ApproxNearestNeighborsMTree
   : public ApproxNearestNeighbors
   {
   public:

      ApproxNearestNeighborsMTree(int dim,
                                  const string& tree_name,
                                  const string& directory_name,
                                  ostream* error_log_stream = (ostream*)NULL,
                                  bool do_error_checking = false                               
                                  );
    
      ~ApproxNearestNeighborsMTree() {delete m_tree;}

      virtual int insert(std::vector<double>& point, uint128_t& key);

      virtual void insert(std::vector<uint128_t> const& keys);

      virtual void remove(int id);

      virtual void knn(std::vector<double> const& x, int k, std::vector<int> &ids, std::vector<uint128_t> &keys, std::vector<double> &dists);

      virtual uint128_t getKey(int id);

   private:

      int m_dim;
      MTree* m_tree;

};

#endif
