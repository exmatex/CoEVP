//
// File:        MTree.T
// Package:     MTree database
// 
// 
// 
// Description: Main Mtree index structure class.
//

template<typename DataPredicate>
inline bool
MTree::writeObjects(const DataPredicate & predicate) const
{

   return d_data_store.writeObjects(predicate);

}

