//
// File:        MTreeDataStore.h
// Package:     MTree database
// 
// 
// 
// Description: Manager class for MTree node allocation and storage of data objects.
//

template<typename DataPredicate>
inline bool
MTreeDataStore::writeObjects(const DataPredicate & predicate)
{

   bool write_successful = false;

   //
   // iterate through all leaf nodes
   //

   unsigned int object_id = 0;
   int object_count = 0;
   while ( (object_id < d_object_info.size()) &&
           (object_count < d_num_objects) ) {
	
      if ( isValidObjectId(object_id) ) {
	  
         if ( !(d_object_info[object_id]->getInFile()) ) {
	
	    //
	    // get object pointer
	    //

	    DBObjectPtr objectPointer = getObjectPtr(object_id);

	    //
	    // write object if predicate is true
	    //

	    if (predicate(objectPointer) == true)
               write_successful |= writeDataObject(object_id);
	    
    
         }

         ++object_count;

      }

      ++object_id;
	
   }

   return write_successful;
   
}
    

