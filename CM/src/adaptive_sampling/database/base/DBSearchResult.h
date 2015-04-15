//
// File:        DBSearchResult.h
// 
// 
// 
// Description: Container for single data object resulting from DB search
//

#ifndef included_DBSearchResult
#define included_DBSearchResult

#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif
#ifndef included_DBObjectFactory
#include <base/DBObjectFactory.h>
#endif

/*!
 * @brief DBSearchResult is a container for information about
 * a single data object retrieved from a DB during a search. 
 * 
 * @see MetricSpacePoint 
 */

class DBSearchResult
{
public:

   /*!
    * Default ctor for DB search result sets result
    * information to undefined state.
    */
   DBSearchResult();

   /*!
    * Copy ctor for a DB search result.
    */
   DBSearchResult(const DBSearchResult& result);
 
   /*!
    * Dtor for DB search result objects.
    */
   virtual ~DBSearchResult();

   /*!
    * Less than operator to compare search results based on 
    * distance to query point.  This is used by the STL list
    * function sort() for sorting results of range search.
    */
   virtual int operator< (const DBSearchResult& rhs) const {return 0;}

   /*!
    * DB search result copy assignment operator.
    */
   virtual DBSearchResult& operator=(const DBSearchResult& rhs); 

   /*!
    * Return const reference to data object associated with search result.
    */
   virtual const DBObject& getDataObject() const;

   /*!
    * Return distance of result to query point.
    */
   virtual double getDistanceToQueryPoint() const;

   /*!
    * Return const reference to point representing data object.
    */
   virtual const MetricSpacePoint& getDataObjectPoint() const;

   /*!
    * Return radius of data object.
    */
   virtual double getDataObjectRadius() const;

   /*!
    * Return const reference to query point associated with search result.
    */
   virtual const MetricSpacePoint& getQueryPoint() const;

   bool isValidResult() const {return d_is_valid_result;}

   //   virtual void finalizeSearchResult(MTreeDataStore& data_store,
   //                                     bool make_safe) {;}

protected:

   double          d_distance_to_query_point;
   DBObjectPtr     d_data_object;
   double          d_data_object_radius;
   bool            d_is_valid_result;
   MetricSpacePointPtr   d_query_point;
   int             d_data_object_id;
   MetricSpacePointPtr   d_data_object_point;

private:
};

#ifndef DEBUG_NO_INLINE
#include "DBSearchResult.I"
#endif
#endif
