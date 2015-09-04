//
// File:        MTreeObject.h
// 
// 
// 
// Description: Abstract base class for database objects.
//

#ifndef included_MTreeObject
#define included_MTreeObject

#ifndef included_toolbox_Database
#include "toolbox/database/Database.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif


class MTreeObject;
typedef std::shared_ptr<MTreeObject> MTreeObjectPtr;

/*!
 * @brief MTreeObject is an abstract base class from which all database
 *        objects must be derived.
 *
 * By default, an object has an undefined object id.  For safety, and 
 * to preserve the uniqueness of object ids, an object id should only 
 * be set by the data store class.  Typically, this is done when an object 
 * is added to the data store.  The setObjectId() method is declared private.
 */

class MTreeObject
{
public:
   friend class MTreeDataStore;
   friend class MTreeSearchResult;

   /*!
    * Static function to get common undefined integer identifier for object.
    */
   static int getUndefinedId();

   /*!
    * Database object ctor.
    */
   MTreeObject();

   /*!
    * Virtual dtor for database objects.
    */
   virtual ~MTreeObject();

   /*!
    * Pure virtual method to create and return smart pointer to a
    * (deep) copy of this object.
    */
   virtual MTreeObjectPtr makeCopy() const = 0;

   /*!
    * Pure virtual method to write data members to given database.
    */
   virtual void putToDatabase(toolbox::Database& db) const = 0; 

   /*!
    * Vitual method to print concrete object data to the specified
    * output stream.  This is optional; a default no-op is supplied here.
    */
   virtual ostream& print(ostream& stream) const;

   /*!
    * Return integer identifier of object.
    */
   int getObjectId() const;

private:
   // The following are not implemented
   MTreeObject(const MTreeObject&);
   void operator=(const MTreeObject&);

   /*
    * Set object id to given integer value.  Generally, this method must only
    * be called once.  To avoid some (but not all) problems, this function
    * will throw an assertion (when assertion checking is active) if it is 
    * called and the object already has a valid object id.  If this occurs, 
    * there is a likely a problem somewhere.
    */
   void setObjectId(int id);

   void writeToDatabase(toolbox::Database& db) const;

   static int s_undefined_object_id;

   int d_object_id;

};

#ifndef DEBUG_NO_INLINE
#include "MTreeObject.I"
#endif
#endif
