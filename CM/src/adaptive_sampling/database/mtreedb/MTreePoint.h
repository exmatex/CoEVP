#ifndef included_mtreedb_MTreePoint
#define included_mtreedb_MTreePoint

//#ifndef included_boost_shared_ptr
//#define included_boost_shared_ptr
//#include <boost/shared_ptr.hpp>
//#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_toolbox_Database
#include "toolbox/database/Database.h"
#endif

namespace mtreedb {

class MTreePoint;
//typedef boost::shared_ptr<MTreePoint> MTreePointPtr;
typedef std::shared_ptr<MTreePoint> MTreePointPtr;

/*!
 * @brief MTreePoint is an abstract base class declaring the interface for 
 * points (i.e., feature vectors) of entries indexed by an MTree structure.  
 * A representative point characterizes each data entry or routing entry 
 * in the tree.
 */

class MTreePoint 
{
public:
   /*!
    * Static member function that returns the maximum possible distance
    * between points in the MTree.  This value is either the value set 
    * by calling the setMaxDistance() function, or the default value 
    * value chosen as a reasonable conservative maximum.  The value
    * should be some reasonably large positive value.  It is used to 
    * optimize certain tree traversal operations by reducing 
    * the number of distance computations performed.
    */
   static double getMaxDistance();

   /*!
    * Static member function that sets the maximum possible distance
    * between points in the MTree to the given value.  The given 
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
    * Default ctor for an MTree point object.
    */
   MTreePoint();

   /*!
    * Virtual dtor for MTree point objects.
    */
   virtual ~MTreePoint();

   /*!
    * Pure virtual method to create and return smart pointer to a 
    * (deep) copy of this point object.
    */
   virtual MTreePointPtr makeCopy() const = 0;

   /*!
    * Pure virtual method to compute and return distance between this point 
    * and another point given as a const reference.
    */ 
   virtual double computeDistanceTo(const MTreePoint& other) const = 0;

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
    * and another point given as an MTreePointPtr.  Method is for
    * convenience and calls computeDistanceTo(const MTreePoint& other) 
    * version.
    */ 
   double computeDistanceTo(MTreePointPtr other) const;

   /*!
    * Vitual method to print concrete point object data to the specified 
    * output stream.   This is optional; a default no-op is supplied here.
    */
   virtual ostream& print(ostream& stream) const;

private:
   // The following are not implemented
   MTreePoint(const MTreePoint&);
   MTreePoint& operator=(const MTreePoint&);

   static double s_max_distance;

};


}
#ifndef DEBUG_NO_INLINE
#include "MTreePoint.I"
#endif
#endif
