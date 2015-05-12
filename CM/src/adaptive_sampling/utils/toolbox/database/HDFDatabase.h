//
// File:        HDFDatabase.h
// Package:     toolbox
// 
// 
// 
// Description: A database structure that stores HDF5 format data.
//
 
#ifndef included_toolbox_HDFDatabase
#define included_toolbox_HDFDatabase
 
/*
************************************************************************
*  THIS CLASS WILL BE UNDEFINED IF THE LIBRARY IS BUILT WITHOUT HDF5
************************************************************************
*/
#ifdef HAVE_PKG_hdf5

#ifndef included_toolbox_Database
#include "toolbox/database/Database.h"
#endif

#ifndef included_list
#define included_list
#include <list>
using namespace std;
#endif

#ifndef included_hdf5
#define included_hdf5
#undef RCSID
#include <hdf5.h>
#endif

namespace toolbox {
 
class HDFDatabase;
typedef std::shared_ptr<HDFDatabase> HDFDatabasePtr;
 
/*!
 * @brief HDFDatabase implements the interface of the Database
 * class to store data in the HDF5 (Hierarchical Data Format) data 
 * format.
 *
 * It is assumed that all processors will access the database in the same
 * manner.  Error reporting is done using the error 
 * reporting macros.
 *
 * @see toolbox::Database
 */

class HDFDatabase : public Database
{
public:
   /*!
    * The HDF database constructor creates an empty database with the
    * specified name.  By default the database will not be associated
    * with a file until it is mounted, using the mount() member function.
    * 
    * When assertion checking is active, the name string must be non-empty.
    */
   HDFDatabase(const string& name);

   /*!
    * The database destructor closes the HDF5 group or data file.
    */
   virtual ~HDFDatabase();

  /*!
    * Check whether key exists in database.
    *  
    * @return boolean true if the specified key exists in the database 
    * and false otherwise.
    *
    * @param key Key name to lookup.
    */
   virtual bool keyExists(const string& key);
 
   /*!
    * Retrieve vector of all keys in the database.
    */
   virtual void getAllKeys(vector<string>& keys);
 
   /*!
    * Return the size of the array associated with the key.  If the key is
    * associated with a scalar value, then one is returned.  If the key
    * does not exist, then zero is returned.
    *
    * @param key Key name in database.
    */
   virtual int getArraySize(const string& key);

   /*!
    * Check whether the specified key represents a database entry. 
    *  
    * @return boolean true if the key is associated with a database
    * entry; otherwise false.
    *
    * @param key Key name in database.
    */
   virtual bool isDatabase(const string& key);
 
   /*!
    * Create a new database with the specified key name and return
    * a smart pointer to it.  If the key already exists in the database,
    * then the old key record is deleted and the new one is silently
    * created in its place.
    *
    * @param key Key name in database.
    */
   virtual DatabasePtr putDatabase(const string& key);
 
   /*!
    * Get a pointer to the database with the specified key name.  If the
    * specified key does not exist in the database or it is not a database,
    * then an error message is printed and the program exits.
    *
    * @param key Key name in database.
    */
   virtual DatabasePtr getDatabase(const string& key);

   /*!
    * Return whether the specified key represents a boolean entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isBool(const string& key);
 
   /*!
    * Create a boolean scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putBool(const string& key, bool data);
 
   /*!
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first boolean data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void putBoolArray(
      const string& key, const bool* const data, const int nelements);

   /*!
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to vector with bool data to put in database.
    */
   virtual void putBoolArray(
      const string& key, const vector<bool>& data);
 
   /*!
    * Get a boolean scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * boolean scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual bool getBool(const string& key);

   /*!
    * Get a boolean array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place bool data.
    */
   virtual void getBoolArray(const string& key, vector<bool>& data);
 
   /*!
    * Get a boolean array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first boolean data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void getBoolArray(
      const string& key, bool* data, const int nelements);

   /*!
    * Return whether the specified key represents a character entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isChar(const string& key);
 
   /*!
    * Create a character scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putChar(const string& key, char data);
 
   /*!
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first char data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void putCharArray(
      const string& key, const char* const data, const int nelements);
 
   /*!
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to char vector with data to put in database.
    */
   virtual void putCharArray(
      const string& key, const vector<char>& data);
 
   /*!
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * character scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual char getChar(const string& key);

   /*!
    * Get a character array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first char data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void getCharArray(
      const string& key, char* data, const int nelements);
 
   /*!
    * Get a character array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place char data.
    */
   virtual void getCharArray(const string& key, vector<char>& data);

   /*!
    * Return whether the specified key represents a double entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isDouble(const string& key);
 
   /*!
    * Create a double scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putDouble(const string& key, double data);
 
   /*!
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first double data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void putDoubleArray(
      const string& key, const double* const data, const int nelements);
 
   /*!
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to double vector with data to put in database.
    */
   virtual void putDoubleArray(
      const string& key, const vector<double>& data);
 
   /*!
    * Get a double scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * double scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual double getDouble(const string& key);

   /*!
    * Get a double array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place double data.
    */
   virtual void getDoubleArray(const string& key, vector<double>& data);
 
   /*!
    * Get a double array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first double data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void getDoubleArray(
      const string& key, double* data, const int nelements);

   /*!
    * Return whether the specified key represents a float entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isFloat(const string& key);
 
   /*!
    * Create a float scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putFloat(const string& key, float data);
 
   /*!
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first float data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void putFloatArray(
      const string& key, const float* const data, const int nelements);
 
   /*!
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to vector with float data to put in database.
    */
   virtual void putFloatArray(
      const string& key, const vector<float>& data);
 
   /*!
    * Get a float scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * float scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual float getFloat(const string& key);

   /*!
    * Get a float array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place float data.
    */
   virtual void getFloatArray(const string& key, vector<float>& data);
 
   /*!
    * Get a float array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first float data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void getFloatArray(
      const string& key, float* data, const int nelements);

   /*!
    * Return whether the specified key represents a integer entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isInteger(const string& key);
 
   /*!
    * Create a integer scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putInteger(const string& key, int data);
 
   /*!
    * Create a integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first integer data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void putIntegerArray(
      const string& key, const int* const data, const int nelements);
 
   /*!
    * Create a integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to integer vector with data to put in database.
    */
   virtual void putIntegerArray(
      const string& key, const vector<int>& data);
 
   /*!
    * Get a integer scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * integer scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual int getInteger(const string& key);

   /*!
    * Get a integer array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place integer data.
    */
   virtual void getIntegerArray(const string& key, vector<int>& data);
 
   /*!
    * Get a integer array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first integer data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void getIntegerArray(
      const string& key, int* data, const int nelements);

   /*!
    * Return whether the specified key represents a string entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isString(const string& key);
 
   /*!
    * Create a string scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putString(const string& key, const string& data);
 
   /*!
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first string data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void putStringArray(
      const string& key, const string* const data, const int nelements);
 
   /*!
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to string vector with data to put in database.
    */
   virtual void putStringArray(
      const string& key, const vector<string>& data);
 
   /*!
    * Get a string scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * string scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual string getString(const string& key);

   /*!
    * Get a string array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place string data.
    */
   virtual void getStringArray(const string& key, vector<string>& data);
 
   /*!
    * Get a string array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.  The specified number of elements must match
    * exactly the number of elements in the array in the database.
    *
    * @param key       Key name in database.
    * @param data      Pointer to first string data value in array.
    * @param nelements Number of elements to write from array.
    */
   virtual void getStringArray(
      const string& key, string* data, const int nelements);

   /*!
    * Print contents of current database to the specified output stream.  
    * If no output stream is specified, then data is written to stream pout.
    * Note that none of the subdatabases contained in the current database 
    * will have their contents printed.  To view the contents of any other
    * database, you must call this print routine for that database. 
    */
   virtual void printClassData(ostream& os);

   /*!
    * Open a database file associated with a file name.  
    *
    * @return  integer indicating status of file open operation.  If 
    * the file does not exist when read only is specified, an error 
    * status is returned (<0).  If the file exists when opened for write, 
    * an error status is returned (<0).  Otherwise, a successful status 
    * (=1) is returned.
    * 
    * @param file_name Const reference to string name of the file to open.  
    * @param flags Const reference to string indicating the file access
    *              requested. Allowed values are: "R" (open existing file
    *              read only access), "W" (open existing file for 
    *              read/write access), "WN" (create new file with 
    *              read/write access -- if file already exists, its 
    *              contents will be cleared for new data to be written).
    *
    * When assertion checking is active, strings must be non-empty.
    */
   virtual int mount(const string& file_name,
                     const string& flags);

   /*!
    * Close the database file.
    */
   virtual void unmount();

private:
   HDFDatabase(const HDFDatabase&);   // not implemented
   void operator=(const HDFDatabase&);     // not implemented

   /*
    * Static function passed HDF5 iterator routine to look up database keys.
    */
   static herr_t iterateKeys(hid_t loc_id,
                             const char *name,
                             void *opdata);

   /*
    * Static member used to construct list of keys when searching for
    * database keys via the HDF5 iterator routine.
    */
   static void addKeyToList(const char *name,
                            int type,
                            void* database);

   /*
    * Private constructor used internally to create sub-databases.
    */
   HDFDatabase(const string& name, hid_t group_ID);

   /*
    * Private utility routine for inserting array data in the database
    */
   void insertArray(hid_t parent_id,
                    const char *name,
                    hsize_t offset,
                    int ndims,
                    const hsize_t dim[/*ndims*/],
                    const int *perm,
                    hid_t member_id) const;

   /*
    * Private utility routines for searching keys in database;
    */
   void performKeySearch();
   void cleanupKeySearch();

   /*!
    * @brief Write attribute for a given dataset.
    *
    * Currently only one attribute is kept for each dataset: its type.
    * The type attribute is used to determine what kind of data the
    * dataset represent.
    *
    * @param type_key Type identifier for the dataset
    * @param dataset_id The HDF dataset id
    */
    void writeAttribute(int type_key,
		        hid_t dataset_id);

   /*!
    * @brief Read attribute for a given dataset.
    *
    * Currently only one attribute is kept for each dataset: its type.
    * The type attribute is returned.
    *
    * @param dataset_id The HDF dataset id
    * @return type attribute
    */
   int readAttribute( hid_t dataset_id );

   struct hdf_complex {
      double re;
      double im;
   };

   /*
    * The following structure is used to store (key,type) pairs when
    * searching for keys in the database.
    */
   struct KeyData {
      string d_key;   // group or dataset name
      int    d_type;  // type of entry
   };

   static string s_top_level_search_group;
   static string s_group_to_search;
   static int s_still_searching;
   static int s_found_group;

   /*
    * HDF5 file and group id, boolean flag indicating whether database 
    * is associated with a mounted file, and name of this database object.
    */
   /*!
     @brief Whether database is mounted to a file
   */
   bool  d_is_file;

   /*!
     @brief ID of file attached to database

     Is either -1 (not mounted to file) or set to the return value from
     opening a file.
     Set to -1 on unmounting the file.
   */
   hid_t d_file_id;

   /*!
     @brief ID of group attached to database

     A database object is always attached to a group.
     The group id is set in the constructor when constructing from a group.
     If the object mounts a file, the group id is the file id.
   */
   hid_t d_group_id;

   /*
    * Name of this database object (passed into constructor, and list
    * of (key,type) pairs assembled when searching for keys.
    */ 
   string d_database_name;

   list<KeyData> d_keydata;

};


}

#endif

#endif
