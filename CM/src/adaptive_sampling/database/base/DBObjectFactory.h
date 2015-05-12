//
// File:        DBObjectFactory.h
// 
// 
// 
// Description: Abstract base class for database object factory.
//

#ifndef included_DBObjectFactory
#define included_DBObjectFactory

#ifndef included_DBObject
#include "DBObject.h"
#endif

#ifndef included_toolbox_Database
#include "toolbox/Database.h"
#endif


/*!
 * @brief DBObjectFactory is an abstract base class used to allocate 
 *        database objects and a concrete instance must be provided
 *        for database objects.
 *
 * The separation of data object creation implemented by the factory 
 * class from the data object constructor allows concrete objects to 
 * be constructed without knowing the concrete type specifically.
 */

class DBObjectFactory
{
public:

   /*!
    * DBObjectFactory default ctor.
    */
   DBObjectFactory();

   /*!
    * Virtual dtor for DBObjectFactory.
    */
   virtual ~DBObjectFactory();

   /*!
    * Pure virtual method to create and return smart pointer to a
    * new data object and set its data members from the contents of
    * the given database. 
    */
   virtual DBObjectPtr allocateObject(toolbox::Database& db) const = 0;

   /*!
    * Virtual method to create an exact copy of the given data object
    * and return smart pointer to the copy.  
    *
    * The default implementation calls the DBObject method makeCopy().
    */
   virtual DBObjectPtr cloneObject(const DBObject& object) const;

private:
   // The following are not implemented
   DBObjectFactory(const DBObjectFactory&);
   void operator=(const DBObjectFactory&);

};

#endif
