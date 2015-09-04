#include "ApproxNearestNeighborsMTree.h"
#include "MTreeKeyObjectFactory.h"

#include "base/ResponsePoint.h"


ApproxNearestNeighborsMTree::ApproxNearestNeighborsMTree( int dim,
                                                          const string& tree_name,
                                                          const string& directory_name,
                                                          ostream* error_log_stream,
                                                          bool do_error_checking )
   : m_dim(dim)
{
   m_tree = new MTree(tree_name, error_log_stream, do_error_checking);

   m_tree->initializeCreate(directory_name + "/" 
                            "kriging_model_database",
                            "krigcpl",
                            *(new krigcpl::MTreeKeyObjectFactory<uint128_t>));

   m_tree->setMaxNodeEntries(12);
}


int
ApproxNearestNeighborsMTree::insert(std::vector<double>& point,
                                    uint128_t&           key)
{
   krigcpl::ResponsePoint rpt(point.size(), point.data());

   krigcpl::MTreeKeyObject<uint128_t> mtree_object(key);

   m_tree->insertObject(mtree_object, rpt, 0.0);

   return mtree_object.getObjectId();
}


void
ApproxNearestNeighborsMTree::insert(std::vector<uint128_t> const& keys)
{
   std::cout << "ApproxNearestNeighborsMTree::insert(std::vector<uint128_t> const& keys) not implemented" << std::endl;
   exit(1);
}
    
    
void
ApproxNearestNeighborsMTree::remove(int id)
{
   m_tree->deleteObject(id);
}
    
    
void
ApproxNearestNeighborsMTree::knn(std::vector<double> const& x,
                                 int k_neighbors,
                                 std::vector<int> &ids,
                                 std::vector<uint128_t> &keys,
                                 std::vector<double> &dists)
{
   krigcpl::ResponsePoint point(x.size());

   for (int i=0; i<x.size(); ++i) {
      point[i] = x[i];
   }

   std::vector<MTreeSearchResult> searchResults;

   m_tree->searchKNN(searchResults,
                     point,
                     k_neighbors);

   int num_neighbors_found = searchResults.size();

   if (num_neighbors_found > 0) {

      ids.resize(num_neighbors_found);
      keys.resize(num_neighbors_found);
      dists.resize(num_neighbors_found);

      for (int i=0; i<num_neighbors_found; ++i) {

         const krigcpl::MTreeKeyObject<uint128_t> & mtree_object = 
            dynamic_cast<const krigcpl::MTreeKeyObject<uint128_t> &>(searchResults[i].getDataObject()); 

         ids[i]   = mtree_object.getObjectId();
         keys[i]  = mtree_object.getKey();
         dists[i] = searchResults[i].getDistanceToQueryPoint();
      }
   }
   else {
      ids.resize(0);
      keys.resize(0);
      dists.resize(0);
   }
}


uint128_t
ApproxNearestNeighborsMTree::getKey(int id)
{
   uint128_t key;

   const DBObjectPtr mtree_object_ptr = m_tree->getObject(id);

   if (mtree_object_ptr == NULL) {
      key = uint128_t_undefined;
   }
   else {

      const krigcpl::MTreeKeyObject<uint128_t> & mtree_object = 
         dynamic_cast<const krigcpl::MTreeKeyObject<uint128_t> &>(*mtree_object_ptr); 

      key = mtree_object.getKey();
   }

   return key;
}
