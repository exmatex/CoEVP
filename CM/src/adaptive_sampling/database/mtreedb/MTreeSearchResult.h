//
// File:        MTreeSearchResult.h
// Package:     MTree database
// 
// 
// 
// Description: Container for single data object resulting from MTree search
//

#ifndef included_mtreedb_MTreeSearchResult
#define included_mtreedb_MTreeSearchResult

#ifndef included_mtreedb_MTreePoint
#include "MTreePoint.h"
#endif
#ifndef included_mtreedb_MTreeEntry
#include "MTreeEntry.h"
#endif
#ifndef included_mtreedb_MTreeObjectFactory
#include "MTreeObjectFactory.h"
#endif
#ifndef included_mtreedb_MTreeDataStore
#include "MTreeDataStore.h"
#endif

namespace mtreedb {

/*!
 * @brief MTreeSearchResult is a container for information about
 * a single data object retrieved from an MTree during a search. 
 * 
 * @see mtreedb::MTreePoint 
 */

class MTreeSearchResult
{
public:
   friend class MTree;

   /*!
    * Default ctor for MTree search result sets result
    * information to undefined state.
    */
   MTreeSearchResult();

   /*!
    * Copy ctor for a MTree search result.
    */
   MTreeSearchResult(const MTreeSearchResult& result);
 
   /*!
    * Dtor for MTree search result objects.
    */
   virtual ~MTreeSearchResult();

   /*!
    * MTree search result copy assignment operator.
    */
   MTreeSearchResult& operator=(const MTreeSearchResult& rhs); 

   /*!
    * Return const reference to data object associated with search result.
    */
   const MTreeObject& getDataObject() const;

   /*!
    * Return distance of result to query point.
    */
   double getDistanceToQueryPoint() const;

   /*!
    * Return const reference to point representing data object.
    */
   const MTreePoint& getDataObjectPoint() const;

   /*!
    * Return radius of data object.
    */
   double getDataObjectRadius() const;

   /*!
    * Return const reference to query point associated with search result.
    */
   const MTreePoint& getQueryPoint() const;

   /*!
    * Less than operator to compare search results based on 
    * distance to query point.  This is used by the STL list
    * function sort() for sorting results of range search.
    */
   int operator< (const MTreeSearchResult& rhs) const;
   
private:
   /*
    * Private ctor for MTree search result class (called 
    * by MTree) sets query point and data object information 
    * from entry. Distance from data object to query point 
    * is not set.
    *
    * There is no error checking to make sure entry holds
    * good data object information.  However, when assertion 
    * checking is on, an assertion is thrown if a null pointer 
    * is passed.
    */
   MTreeSearchResult(MTreePointPtr query_point,
                     MTreeEntryPtr entry);

   /*!
    * Private function to set query point (called by MTree).
    */
   void setQueryPoint(MTreePointPtr query_point);

   /*
    * Private function to set distance of result to query point 
    * (called by MTree). There is no error checking.
    */
   void setDistanceToQueryPoint(double distance);

   /*
    * Private function to determine whether result object
    * contains a valid data object description (called by MTree).
    */
   bool isValidResult() const;
 
   /*
    * Private function to set object protection.
    */
   void setObjectProtection(bool make_safe);
 
   /*
    * Private function to finalize search result (called by MTree)
    * replaces data object and point with deep copies, if
    * the boolean argument is true.
    * This ensures integrity of data in MTree.  
    */
   void finalizeSearchResult(MTreeDataStore& data_store,
                             bool make_safe);

   MTreePointPtr   d_query_point;
   double          d_distance_to_query_point;

   MTreeObjectPtr  d_data_object;
   int             d_data_object_id;
   MTreePointPtr   d_data_object_point;
   double          d_data_object_radius;

   bool            d_is_valid_result;
};

}
#ifndef DEBUG_NO_INLINE
#include "MTreeSearchResult.I"
#endif
#endif
