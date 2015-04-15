//
// File:        MTreeSearchResult.h
// Package:     MTree database
// 
// 
// 
// Description: Container for single data object resulting from MTree search
//

#ifndef included_MTreeSearchResult
#define included_MTreeSearchResult

#ifndef included_DBSearchResult
#include <base/DBSearchResult.h>
#endif
#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif
#ifndef included_MTreeEntry
#include "MTreeEntry.h"
#endif
#ifndef included_DBObjectFactory
#include <base/DBObjectFactory.h>
#endif
#ifndef included_MTreeDataStore
#include "MTreeDataStore.h"
#endif

/*!
 * @brief MTreeSearchResult is a container for information about
 * a single data object retrieved from an MTree during a search. 
 * 
 * @see MetricSpacePoint 
 */

class MTreeSearchResult
  : public DBSearchResult
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
   //   MTreeSearchResult(const MTreeSearchResult& result);
 
   /*!
    * Dtor for MTree search result objects.
    */
   virtual ~MTreeSearchResult();

   /*!
    * Less than operator to compare search results based on 
    * distance to query point.  This is used by the STL list
    * function sort() for sorting results of range search.
    */
   virtual int operator< (const MTreeSearchResult& rhs) const;
   
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
   MTreeSearchResult(MetricSpacePointPtr query_point,
                     MTreeEntryPtr entry);

   /*!
    * Private function to set query point (called by MTree).
    */
   void setQueryPoint(MetricSpacePointPtr query_point);

   /*
    * Private function to set distance of result to query point 
    * (called by MTree). There is no error checking.
    */
   void setDistanceToQueryPoint(double distance);

   /*
    * Private function to determine whether result object
    * contains a valid data object description (called by MTree).
    */
   //   virtual bool isValidResult() const;
 
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


};

#ifndef DEBUG_NO_INLINE
#include "MTreeSearchResult.I"
#endif
#endif
