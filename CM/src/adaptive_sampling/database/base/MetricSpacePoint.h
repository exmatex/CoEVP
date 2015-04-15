#ifndef included_MetricSpacePoint
#define included_MetricSpacePoint

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_toolbox_Database
#include "toolbox/database/Database.h"
#endif


class MetricSpacePoint;
typedef std::shared_ptr<MetricSpacePoint> MetricSpacePointPtr;

/*!
 * @brief MetricSpacePoint is an abstract base class declaring the interface for 
 * points (i.e., feature vectors) of a metric space.  
 */

class MetricSpacePoint 
{
public:
   /*!
    * Static member function that returns the maximum possible distance
    * between points in the metric space.  This value is either the value set 
    * by calling the setMaxDistance() function, or the default value 
    * value chosen as a reasonable conservative maximum.  The value
    * should be some reasonably large positive value.  It is used to 
    * optimize certain MTree traversal operations by reducing 
    * the number of distance computations performed.
    */
   static double getMaxDistance();

   /*!
    * Static member function that sets the maximum possible distance
    * between points in the metric space to the given value.  The given 
    * value is returned by subsequent calls to getMaxDistance().  
    *  
    * Note that there are no checks to ensure that this routine is 
    * not called excessively.  If it is used to reset the distance
    * after the getMaxDistance() function has been called, unanticipated
    * behavior may result.
    *
    * When assertion checking is active, an assertion is thrown if the 
    * value is <= 0.0.
    */
   static void setMaxDistance(double max_dist);

   /*!
    * Default ctor for a metric space point object.
    */
   MetricSpacePoint();

   /*!
    * Virtual dtor for metric space point objects.
    */
   virtual ~MetricSpacePoint();

   /*!
    * Pure virtual method to create and return smart pointer to a 
    * (deep) copy of this point object.
    */
   virtual MetricSpacePointPtr makeCopy() const = 0;

   /*!
    * Pure virtual method to compute and return distance between this point 
    * and another point given as a const reference.
    */ 
   virtual double computeDistanceTo(const MetricSpacePoint& other) const = 0;

   /*!
    * Pure virtual method to put point data to given database object.
    */
   virtual void putToDatabase(toolbox::Database& db) const = 0;

   /*!
    * Pure virtual method to get point data from given database object.
    */
   virtual void getFromDatabase(toolbox::Database& db) = 0;

   /*!
    * Method to compute and return distance between this point 
    * and another point given as an MetricSpacePointPtr.  Method is for
    * convenience and calls computeDistanceTo(const MetricSpacePoint& other) 
    * version.
    */ 
   double computeDistanceTo(MetricSpacePointPtr other) const;

   /*!
    * Virtual method to print concrete point object data to the specified 
    * output stream.   This is optional; a default no-op is supplied here.
    */
   virtual ostream& print(ostream& stream) const;

private:
   // The following are not implemented
   MetricSpacePoint(const MetricSpacePoint&);
   MetricSpacePoint& operator=(const MetricSpacePoint&);

   static double s_max_distance;

};


#ifndef DEBUG_NO_INLINE
#include "MetricSpacePoint.I"
#endif
#endif
