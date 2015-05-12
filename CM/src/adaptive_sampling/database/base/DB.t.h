//
// File:        DB.T
// 
// 
// 
// Description: Main DB index structure class.
//

template<typename DataPredicate>
inline bool
DB::writeObjects(const DataPredicate & predicate) const
{

   return d_data_store.writeObjects(predicate);

}

