//
// File:        MTreeQuery.h
// Package:     MTree database
// 
// 
// 
// Description: Utility class used in MTree searches
//

#ifndef included_MTreeQuery
#define included_MTreeQuery

#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif
#ifndef included_MTreeEntry
#include "MTreeEntry.h"
#endif

class MTree;

/*!
 * @brief MTreeQuery is a simple utility class used in MTree
 * range searches and nearest-neighbor searches.  
 * 
 * This class is used in MTree search operations and should not
 * be used for other stuff.
 */
 
class MTreeQuery
{
public:
   /*!
    * Ctor for MTreeQuery sets query point and radius 
    * (no error checking), sets grade to zero, and caches
    * pointer to MTree (for statistics gathering).
    */
   MTreeQuery(MetricSpacePointPtr query_point,
              double query_radius,
              MTree* tree);

   /*!
    * Copy ctor for MTreeQuery.
    */
   MTreeQuery(const MTreeQuery& query);

   /*!
    * Dtor for MTreeQuery objects. 
    */
   virtual ~MTreeQuery();

   /*!
    * Return query point. 
    */
   MetricSpacePointPtr getPoint() const;

   /*!
    * Return query grade. 
    */
   double getGrade() const;

   /*!
    * Set query grade (no error checking).  
    */
   void setGrade(double grade);

   /*!
    * Return query radius. 
    */
   double getRadius() const;

   /*!
    * Set query radius to given value (no error checking). 
    */
   void setRadius(double radius);

   /*!
    * Return true if entry is consistent (i.e., potentially satisfies)
    * with query and false otherwise.
    */
   bool isConsistent(MTreeEntryPtr entry);

private:
   // The following are not implemented
   MTreeQuery();
   void operator=(const MTreeQuery&);
   
   MetricSpacePointPtr   d_query_point;
   double          d_radius;
   double          d_grade;
   MTree*          d_tree;
};

#ifndef DEBUG_NO_INLINE
#include "MTreeQuery.I"
#endif
#endif
