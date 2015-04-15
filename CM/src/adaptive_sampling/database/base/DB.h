//
// File:        DB.h
// 
// 
// 
// Description: Main DB structure class.
//

#ifndef included_DB
#define included_DB

#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_vector
#define included_vector
#include <vector>
using namespace std;
#endif

#ifndef included_list
#define included_list
#include <list>
using namespace std;
#endif

#ifndef included_DBObjectFactory
#include <base/DBObjectFactory.h>
#endif
#ifndef included_DBSearchResult
#include <base/DBSearchResult.h>
#endif
#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif

/*!
 * @brief DB is an abstract base class defining the capabilities of a
 *        database used in CoEVP's adaptive sampling algorithm.
 *        It was developed by abstracting the original MTree index
 *        structure (a dynamic, paged, metric tree) which is now
 *        a derived class.  This API supports insertion, deletion, and querying of 
 *        data objects that are described relative to each other in terms 
 *        of points in a general metric space.
 *
 * Typical usage of a DB object involves several operations:
 *
 * -# Create the DB by calling the ctor; e.g., DB my_db("MyDB").
 *
 * -# Initialize the DB in one of two ways. 
 *               To set up a new DB from scratch, call the 
 *               initializeCreate() method.  To set up a new DB
 *               and initialize it with data read from existing files, 
 *               call the initializeOpen() method.
 *
 * -# Insert data. Data objects are inserted using the insertObject() 
 *               method. Note that in addition to providing the data 
 *               object to the method, a point representing the object 
 *               and an object radius are required.  The object radius
 *               may be zero if the user desires to index points in the
 *               metric space rather than spherical regions.  However,
 *               each object must be specified by a unique point.
 * 
 * -# Query data. The data in the DB can be searched in two ways. 
 *               The searchKNN() method returns the k-nearest neighbors 
 *               of a given point, and the searchRange() method returns 
 *               all data objects within a given distance of a given point.
 *
 * -# Finalize the DB structure by calling the finalize() method.  
 *               This will write all DB index structure state and
 *               data object information to files for retrieval later.
 * 
 * -# Destroy the DB by calling the dtor explicitly or letting the
 *               DB go out of scope.
 * 
 * @see DBObject
 * @see MetricSpacePoint
 */

class DB
{
public:

   /*!
    * Ctor for DB object sets database name, turns on error checking
    * and sets error logging stream if specified.
    */
   DB(const string& db_name,
      ostream* error_log_stream = (ostream*)NULL,
      bool do_error_checking = false);

   /*!
    * Dtor for DB objects.
    */
   virtual ~DB();

   /*!
    * Return DB name string passed to the ctor.
    */
   string getName() const;

   //@{
   //! @name Methods for initializing/finalizing this DB object.

   /*!
    * Initialize DB to empty state.
    *
    * Either this method, or the initializeOpen() method must be called
    * before any operations involving data objects may be performed.
    * 
    * A DB object can only be initialized once.  Calling this method 
    * more than once or after calling he initializeOpen() method will 
    * result in an unrecoverable error.
    *
    * @param directory_name Const reference to string indicating the
    *                     directory into which all data files
    *                     associated with this DB will be written.
    * @param file_prefix  Const reference to string indicating the
    *                     prefix for all data files associated with
    *                     this DB.
    * @param obj_factory  Const reference to factory that creates 
    *                     concrete data objects to be indexed by this DB.
    */
   virtual void initializeCreate(const string& directory_name,
                                 const string& file_prefix,
                                 const DBObjectFactory& obj_factory) = 0;

   /*!
    * Initialize DB to state contained in existing data files.
    *
    * Either this method, or the initializeCreate() method must be called
    * before any operations involving data objects may be performed.
    *
    * A DB object can only be initialized once.  Calling this method
    * more than once or after calling he initializeCreate() method will
    * result in an unrecoverable error.
    *
    * @param directory_name Const reference to string indicating the
    *                     directory into which all data files associated 
    *                     with this DB will be read and written.
    * @param file_prefix  Const reference to string indicating the
    *                     prefix for all data files associated with
    *                     this DB. 
    * @param obj_factory  Const reference to factory that creates 
    *                     concrete data objects to be indexed by this DB.
    */
   virtual void initializeOpen(const string& directory_name,
                               const string& file_prefix,
                               const DBObjectFactory& obj_factory) = 0;

   /*!
    * Finalize DB structure.
    *
    * This method will write all DB structure state and
    * data object information to files for retrieval later (e.g.,
    * initialization of another DB object). 
    */
   virtual void finalize() = 0; 

   //@}
  
   //@{
   //! @name Methods for inserting, deleting and retrieving objects in this DB object.

   /*!
    * Insert object into DB and set the identifier of the object.
    * 
    * Note that in addition to providing the data object to the method, 
    * a point representing the object and an object radius are required.  
    * The object radius may be zero if the user desires to index points 
    * in the metric space rather than spherical regions.  However, each 
    * object must be specified by a unique point.
    *
    * To avoid external tampering with database contents, this method 
    * produces internal deep copies of the given point and data objects.
    *
    * @param object  Reference to data object.  The reference is non-const
    *                since the method sets the object identifier to match
    *                the identifier of the database copy.
    * @param point   Const reference to center point of object.
    * @param radius  Double radius of object about center point. When 
    *                assertion checking is on, assertion will result if 
    *                value is less than 0.
    */
   virtual void insertObject(DBObject& object,
                             const MetricSpacePoint& point,
                             double radius) = 0;
   
   /*!
    * Get copy of object indexed by DB given object identifier.
    *
    * @param object_id  Integer identifier of object to delete from DB.
    *                If this is not a valid id for an object indexed by
    *                the DB, the method will return a null pointer.
    */ 
   virtual DBObjectPtr getObject(int object_id) const = 0;

   /*!
    * Delete object from DB. 
    *
    *
    * @param object_id  Integer identifier of object to delete from DB.
    *                If this is not a valid id for an object indexed by
    *                the DB, the method will do nothing. 
    */
   virtual void deleteObject(int object_id) = 0;

   //@}
  
   /*!
    * Search DB for "k" nearest neighbors of given query point.
    * 
    * To avoid external tampering with database contents, this method 
    * returns deep copies of points and data objects within the set of
    * search results.
    *
    * @param results Reference to vector of DBSearchResult objects.
    *                The size of the vector will be the number of search
    *                results, up to a maximum of the given "K".  This
    *                vector will be cleared on entry to method so any
    *                pre-existing data it contains will be lost.  Also,
    *                the vector of results will be sorted with respect 
    *                to increaing distance from the query point.
    * @param query_point  Const reference to search query point.
    * @param k_neighbors  Integer maximum neighbors to return.  If 
    *                this value is <= 0, the result vector will be empty.
    * @param make_safe  Optional boolean indicating whether the search
    *                will return copies of objects in data base (true),
    *                or actual objects in the database (false).  The
    *                default value is false.  Passing a value of true
    *                prevents any possibility of data objects being
    *                modified outside of the database.  In the dafault
    *                case, if the user casts away the const-ness of 
    *                any data object in a search result, future database
    *                operations may produce unexpected behavior.
    */
   virtual void searchKNN(vector<DBSearchResult>& results,
                          const MetricSpacePoint& query_point,
                          int k_neighbors,
                          bool make_safe = false) = 0;

   /*!
    * Search the DB for all data objects within given distance of given 
    * query point.
    *
    * Note that if the given radius is zero, this method will return 
    * all data objects whose representative metric space regions contain
    * the given query point.  If the radius is greater than zero, this
    * method will return all data objects whose regions intersect the query
    * region.
    *
    * To avoid external tampering with database contents, this method
    * returns deep copies of points and data objects within the set of
    * search results.
    *
    * @param results Reference to list of DBSearchResult objects.
    *                The size of the list will be the number of search
    *                results.  This list will be cleared on entry to 
    *                method so any pre-existing data it contains will be 
    *                lost.  Also, the list of results will be sorted with
    *                respect to increasing distance from query point. 
    * @param query_point  Const reference to point at center of search.
    * @param radius  Double radius value indicating search region about 
    *                query point.  If this value is < 0, the result list 
    *                will be empty.
    * @param make_safe  Optional boolean indicating whether the search
    *                will return copies of objects in data base (true),
    *                or actual objects in the database (false).  The
    *                default value is false.  Passing a value of true
    *                prevents any possibility of data objects being
    *                modified outside of the database.  In the dafault
    *                case, if the user casts away the const-ness of 
    *                any data object in a search result, future database
    *                operations may produce unexpected behavior.
    */
   virtual void searchRange(list<DBSearchResult>& results,
                            const MetricSpacePoint& query_point,
                            double radius,
                            bool make_safe = false) = 0;

   /*!
    * Write objects to disk. Return boolean tru if write is successful
    * and false otherwise.
    */
   //   template<typename DataPredicate>
   //   bool writeObjects(const DataPredicate & predicate) const;

   virtual void outputStats(std::ostream & outputStream) = 0;

protected:

   /*
    * String name of DB used mainly in printing and error reporting
    */
   string d_db_name;

   /*
    * Error checking and reporing members; set in constructor
    * but can be changed via member functions.
    * Note: error checking is turned off by default.
    */
   bool     d_do_error_checking;
   ostream* d_error_log_stream;

private:
   // The following are not implemented
   DB(const DB&);
   void operator=(const DB&);
};

#ifndef DEBUG_NO_INLINE
#include "DB.I"
#endif

//#include "DB.t.h"

#endif
