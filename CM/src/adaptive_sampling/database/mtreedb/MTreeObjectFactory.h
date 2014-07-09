//
// File:        MTreeObjectFactory.h
// Package:     MTree database
// 
// 
// 
// Description: Abstract base class for MTree data object factory.
//

#ifndef included_mtreedb_MTreeObjectFactory
#define included_mtreedb_MTreeObjectFactory

#ifndef included_mtreedb_MTreeObject
#include "MTreeObject.h"
#endif

#ifndef included_toolbox_Database
#include "toolbox/Database.h"
#endif

namespace mtreedb {

/*!
 * @brief MTreeObjectFactory is an abstract base class used to allocate 
 *        data objects indexed by an MTree structure and a concrete
 *        instance must be provided for data objects indexed by the MTree.
 *
 * The separation of data object creation implemented by the factory 
 * class from the data object constructor allows concrete objects to 
 * be constructed without knowing the concrete type specifically.
 */

class MTreeObjectFactory
{
public:
   friend class MTreeDataStore;

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
    * The default implementation calls the MTreeObject method makeCopy().
    */
   virtual MTreeObjectPtr cloneObject(const MTreeObject& object) const;

private:
   // The following are not implemented
   MTreeObjectFactory(const MTreeObjectFactory&);
   void operator=(const MTreeObjectFactory&);

};

}
#endif
