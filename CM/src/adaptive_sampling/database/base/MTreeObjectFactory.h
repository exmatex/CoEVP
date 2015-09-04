//
// File:        MTreeObjectFactory.h
// 
// 
// 
// Description: Abstract base class for database object factory.
//

#ifndef included_MTreeObjectFactory
#define included_MTreeObjectFactory

#ifndef included_MTreeObject
#include "MTreeObject.h"
#endif

#ifndef included_toolbox_Database
#include "toolbox/Database.h"
#endif


/*!
 * @brief MTreeObjectFactory is an abstract base class used to allocate 
 *        database objects and a concrete instance must be provided
 *        for database objects.
 *
 * The separation of data object creation implemented by the factory 
 * class from the data object constructor allows concrete objects to 
 * be constructed without knowing the concrete type specifically.
 */

class MTreeObjectFactory
{
public:

   /*!
    * MTreeObjectFactory default ctor.
    */
   MTreeObjectFactory();

   /*!
    * Virtual dtor for MTreeObjectFactory.
    */
   virtual ~MTreeObjectFactory();

   /*!
    * Pure virtual method to create and return smart pointer to a
    * new data object and set its data members from the contents of
    * the given database. 
    */
   virtual MTreeObjectPtr allocateObject(toolbox::Database& db) const = 0;

   /*!
    * Virtual method to create an exact copy of the given data object
    * and return smart pointer to the copy.  
    *
    * The default implementation calls the MtreeObject method makeCopy().
    */
   virtual MTreeObjectPtr cloneObject(const MTreeObject& object) const;

private:
   // The following are not implemented
   MTreeObjectFactory(const MTreeObjectFactory&);
   void operator=(const MTreeObjectFactory&);

};

#endif
