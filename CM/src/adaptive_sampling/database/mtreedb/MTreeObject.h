//
// File:        MTreeObject.h
// Package:     MTree database
// 
// 
// 
// Description: Abstract base class for data objects indexed by an MTree.
//

#ifndef included_mtreedb_MTreeObject
#define included_mtreedb_MTreeObject

#ifndef included_toolbox_Database
#include "toolbox/database/Database.h"
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

//#ifndef included_boost_shared_ptr
//#define included_boost_shared_ptr
//#include <boost/shared_ptr.hpp>
//#endif

namespace mtreedb {

class MTreeObject;
//typedef boost::shared_ptr<MTreeObject> MTreeObjectPtr;
typedef std::shared_ptr<MTreeObject> MTreeObjectPtr;

/*!
 * @brief MTreeObject is an abstract base class from which all data 
 *        objects indexed by an MTree structure must be derived.
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
    * MTree object ctor.
    */
   MTreeObject();

   /*!
    * Virtual dtor for MTree objects.
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

}
#ifndef DEBUG_NO_INLINE
#include "MTreeObject.I"
#endif
#endif
