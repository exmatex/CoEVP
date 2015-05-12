#ifndef included_krigcpl_ResponsePoint
#define included_krigcpl_ResponsePoint

#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif

#ifndef included_krigalg_Point
#include <base/Point.h>
#endif

#ifndef included_toolbox_Database
#include <toolbox/base/Database.h>
#endif

namespace krigcpl {

/*!
 * @brief A ResponsePoint object represents a point in the domain 
 *        of the response function; i.e., where the reponse is evaluated.  
 * 
 *        ResponsePoint provides an concrete implementation of the 
 *        MetricSpacePoint base class from which it is inherited.  
 *        Most ResponsePoint functionality is inherited from the 
 *        krigalg::Point class.
 * 
 * @see MetricSpacePoint
 * @see krigalg::Point
 */

class ResponsePoint :
   public MetricSpacePoint,
   public krigalg::Point 
{
public:
   /*!
    * Default ctor for response point.
    */
   ResponsePoint();
   
   /*!
    * Ctor for a response point object indicating dimension of response point.
    */
   ResponsePoint(int dimension);

   /*!
    * Ctor to create a response point object and initialize to given value.
    */
   ResponsePoint(int dimension, double value);

   /*!
    * Ctor to create a reposnse point object and initialize
    */
   ResponsePoint(int dimension, const double * coordinates);

   /*!
    * Copy ctor for a response point object.
    */
   ResponsePoint(const ResponsePoint& response_pt);

   /*!
    * Virtual destructor for response point objects.
    */
   virtual ~ResponsePoint();

   /*!
    * Ceate and return smart pointer to a (deep) copy of this response 
    * point object.
    */
   MetricSpacePointPtr makeCopy() const;

   /*!
    * Compute and return distance between this reponse point object and 
    * a given response point object.
    */
   double distance(const ResponsePoint& other_pt) const;

   /*!
    * Put point data to given database object.
    */
   void putToDatabase(toolbox::Database& db) const;

   /*!
    * Get point data from given database object.
    */
   void getFromDatabase(toolbox::Database& db);

   /*!
    * Print response point object data to the specified output stream.
    */
   friend ostream& operator<< (ostream& stream,
                               const ResponsePoint& point);
  
private:
   // This function is not implemented, but is declared here to avoid
   // having the compiler provide an implementation automatically.
   void operator=(const ResponsePoint&);

   /*!
    * Implementation of pure virtual distance function declared in 
    * MetricSpacePoint base class.
    */
   double computeDistanceTo(const MetricSpacePoint& other_pt) const;

   /*!
    * Print response point object data to the specified output stream.
    */
   ostream& print(ostream& stream) const;

};

}
#ifndef DEBUG_NO_INLINE
#include "ResponsePoint.I"
#endif
#endif
