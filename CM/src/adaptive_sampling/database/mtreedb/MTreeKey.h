//
// File:        MTreeKey.h
// Package:     MTree database
// 
// 
// 
// Description: Description of metric space region associated with each MTree entry
//

#ifndef included_MTreeKey
#define included_MTreeKey

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif

/*!
 * @brief MTreeKey maintains the representation of the metric space region 
 * associated with each MTreeEntry in an MTree structure.  The metric space 
 * region is defined by a center point (an MetricSpacePoint object) and a radius.
 * 
 * @see MetricSpacePoint
 */

class MTreeKey
{
public:
   /*!
    * Static function to get initial undefined distance-to-parent
    * value common to all MTree key objects.
    */
   static double getUndefinedDistanceToParent();

   /*!
    * Default ctor for a MTree key object.
    */
   MTreeKey();

   /*!
    * Copy ctor for a MTree key object.
    */
   MTreeKey(const MTreeKey& key);

   /*!
    * MTree key ctor that specifies all data members.
    */
   MTreeKey(const MetricSpacePointPtr point,
            double radius, 
            double dist2parent = getUndefinedDistanceToParent());

   /*!
    * Dtor for MTree key objects.
    */
   virtual ~MTreeKey();

   /*!
    * MTree key copy assignment operator.
    */
   MTreeKey& operator=(const MTreeKey& rhs);

   /*!
    * Return const pointer to MetricSpacePoint associated with this 
    * key object.
    */
   const MetricSpacePointPtr getPoint() const;

   /*!
    * Return radius of region associated with this key object.
    */
   double getRadius() const;

   /*!
    * Set key radius to given value (no bounds checking!).
    */
   void setRadius(double radius);

   /*!
    * Return distance from center of region associated with this key object
    * to center of region associated with parent region.
    */
   double getDistanceToParent() const;

   /*!
    * Set distance from center of region associated with this key object
    * to center of region associated with parent region to given value
    * (no bounds checking!).
    */
   void setDistanceToParent(double distance);

   /*!
    * Compute and return distance between center of metric space region 
    * defined by this MTree key object and center of metric space region 
    * defined by argument MTree key object. 
    */
   double computeDistanceTo(const MTreeKey& other) const;

   /*!
    * Print MTree key data to the specified output stream.
    */
   void printClassData(ostream& stream) const;

private:
   /*
    * Undefined distance-to-parent shared by all key instances;
    * cannot be changed.
    */
   static double  s_undefined_distance_to_parent;

   /* 
    * Center point of metric space region defined by key;
    * null by default
    */
   MetricSpacePointPtr  d_point;

   /*
    * Radius of metric space region defined by key;
    * set to infinity by default.
    */
   double         d_radius;

   /*
    * Distance from center point of key region to center point 
    * of parent key region; set to negative max point distance 
    * by default.
    */
   double         d_distance_to_parent;

};

#ifndef DEBUG_NO_INLINE
#include "MTreeKey.I"
#endif
#endif
