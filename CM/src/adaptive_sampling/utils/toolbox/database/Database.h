//
// File:        Database.h
// Package:     toolbox
// 
// 
// 
// Description: An abstract base class for managing database objects in files
//
 
#ifndef included_toolbox_Database
#define included_toolbox_Database
 
#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
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
 
#include <memory>
 
namespace toolbox {
 
class Database;
typedef std::shared_ptr<Database> DatabasePtr;
 
/*!
 * @brief Database is an abstract base class defining a simple interface
 * for putting data into a file and getting data out of files.
 *
 * Each database stores (key,value) pairs in a hierarchy.
 * Each value may be another database, or a boolean, character, double,
 * float, integer, or string.
 *
 * Data is entered into the database through methods of the general form
 * putTYPE(key, TYPE) or putTYPEArray(key, TYPE array), where TYPE is the
 * type of value created.  If the specified key already exists in the
 * database, then the existing key and value are silently deleted.
 *
 * Data is extracted from the database through methods of the general form
 * TYPE = getTYPE(key), where TYPE is the type of value to be returned
 * from the database.
 */
 
class Database
{
public:
   /*!
    * The constructor for the database base class does nothing interesting.
    */
   Database();

   /*!
    * The virtual destructor for the database base class does nothing
    * interesting.
    */
   virtual ~Database();

   /*!
    * Check whether key exists in database.
    *
    * @return boolean true if the specified key exists in the database
    * and false otherwise.
    *
    * @param key Key name to lookup.
    */
   virtual bool keyExists(const string& key) = 0;

   /*!
    * Retrieve vector of all keys in the database.
    */
   virtual void getAllKeys(vector<string>& keys) = 0;
 
   /*!
    * Return the size of the array associated with the key.  If the key is
    * associated with a scalar value, then one is returned.  If the key
    * does not exist, then zero is returned.
    *
    * @param key Key name in database.
    */
   virtual int getArraySize(const string& key) = 0;
 
   /*!
    * Check whether the specified key represents a database entry.
    *
    * @return boolean true if the key is associated with a database
    * entry; otherwise false.
    *
    * @param key Key name in database.
    */
   virtual bool isDatabase(const string& key) = 0;

   /*!
    * Create a new database with the specified key name and return
    * a smart pointer to it.  If the key already exists in the database,
    * then the old key record is deleted and the new one is silently
    * created in its place.
    *
    * @param key Key name in database.
    */
   virtual DatabasePtr putDatabase(const string& key) = 0;
 
   /*!
    * Get a pointer to the database with the specified key name.  If the
    * specified key does not exist in the database or it is not a database,
    * then an error message is printed and the program exits.
    *
    * @param key Key name in database.
    */
   virtual DatabasePtr getDatabase(const string& key) = 0;
 
   /*!
    * Return whether the specified key represents a boolean entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isBool(const string& key) = 0;

   /*!
    * Create a boolean scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putBool(const string& key, bool data) = 0;
 
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
      const string& key, const bool* const data, const int nelements) = 0;

   /*!
    * Create a boolean array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to vector with bool data to put in database.
    */
   virtual void putBoolArray(
      const string& key, const vector<bool>& data) = 0;

   /*!
    * Get a boolean scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * boolean scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual bool getBool(const string& key) = 0;

   /*!
    * Get a boolean array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a boolean array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place bool data.
    */
   virtual void getBoolArray(const string& key, vector<bool>& data) = 0;

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
      const string& key, bool* data, const int nelements) = 0;
 
   /*!
    * Return whether the specified key represents a character entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isChar(const string& key) = 0;

   /*!
    * Create a character scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putChar(const string& key, char data) = 0;
 
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
      const string& key, const char* const data, const int nelements) = 0;

   /*!
    * Create a character array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to char vector with data to put in database.
    */
   virtual void putCharArray(
      const string& key, const vector<char>& data) = 0;

   /*!
    * Get a character entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * character scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual char getChar(const string& key) = 0;
 
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
      const string& key, char* data, const int nelements) = 0;

   /*!
    * Get a character array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a character array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place char data.
    */
   virtual void getCharArray(const string& key, vector<char>& data) = 0;

   /*!
    * Return whether the specified key represents a double entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isDouble(const string& key) = 0;
 
   /*!
    * Create a double scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putDouble(const string& key, double data) = 0;

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
      const string& key, const double* const data, const int nelements) = 0;
 
   /*!
    * Create a double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to double vector with data to put in database.
    */
   virtual void putDoubleArray(
      const string& key, const vector<double>& data) = 0; 

   /*!
    * Get a double scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * double scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual double getDouble(const string& key) = 0;
 
   /*!
    * Get a double array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place double data.
    */
   virtual void getDoubleArray(const string& key, vector<double>& data) = 0;

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
      const string& key, double* data, const int nelements) = 0;

   /*!
    * Return whether the specified key represents a float entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isFloat(const string& key) = 0;
 
   /*!
    * Create a float scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putFloat(const string& key, float data) = 0;

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
      const string& key, const float* const data, const int nelements) = 0;
 
   /*!
    * Create a float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to vector with float data to put in database.
    */
   virtual void putFloatArray(
      const string& key, const vector<float>& data) = 0;
 
   /*!
    * Get a float scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * float scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual float getFloat(const string& key) = 0;

   /*!
    * Get a float array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place float data.
    */
   virtual void getFloatArray(const string& key, vector<float>& data) = 0;
 
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
      const string& key, float* data, const int nelements) = 0;

   /*!
    * Return whether the specified key represents a integer entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isInteger(const string& key) = 0;
 
   /*!
    * Create a integer scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putInteger(const string& key, int data) = 0;
 
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
      const string& key, const int* const data, const int nelements) = 0;
 
   /*!
    * Create a integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to integer vector with data to put in database.
    */
   virtual void putIntegerArray(
      const string& key, const vector<int>& data) = 0;
 
   /*!
    * Get a integer scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * integer scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual int getInteger(const string& key) = 0;

   /*!
    * Get a integer array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place integer data.
    */
   virtual void getIntegerArray(const string& key, vector<int>& data) = 0;
 
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
      const string& key, int* data, const int nelements) = 0;

   /*!
    * Return whether the specified key represents a string entry.  If
    * the key does not exist, then false is returned.
    *
    * @param key Key name in database.
    */
   virtual bool isString(const string& key) = 0;
 
   /*!
    * Create a string scalar entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key Key name in database.
    * @param data Value to put into database.
    */
   virtual void putString(const string& key, const string& data) = 0;
 
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
      const string& key, const string* const data, const int nelements) = 0;
 
   /*!
    * Create a string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key  Key name in database.
    * @param data Const reference to string vector with data to put in database.
    */
   virtual void putStringArray(
      const string& key, const vector<string>& data) = 0;

   /*!
    * Get a string scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not a
    * string scalar, then an error message is printed and the program
    * exits.
    *
    * @param key Key name in database.
    */
   virtual string getString(const string& key) = 0;

   /*!
    * Get a string array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.
    *
    * @param key Key name in database.
    * @param data Reference to vector in which to place string data.
    */
   virtual void getStringArray(const string& key, vector<string>& data) = 0;
 
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
      const string& key, string* data, const int nelements) = 0;

   /*!
    * Get a bool scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * bool scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Reference to bool scalar to read from database.
    */
   void getScalar(const string& key, bool& scalar);

   /*!
    * Put a bool entry into the database with the specified key name.  
    *
    * @param key    Key name in database.
    * @param scalar Bool scalar to put into database.
    */
   void putScalar(const string& key, bool scalar);

   /*!
    * Get a bool array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a bool array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Reference to bool vector to read from database.
    */
   void getArray(const string& key, vector<bool>& array);

   /*!
    * Create an bool array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Const reference to bool vector to put into database.
    */
   void putArray(const string& key, const vector<bool>& array);

   /*!
    * Get a scalar char entry from the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * char scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Char scalar value read from database.
    */
   void getScalar(const string& key, char& scalar);

   /*!
    * Put a char entry into the database with the specified key name.
    *
    * @param key    Key name in database.
    * @param scalar Scalar char value to put into database.
    */
   void putScalar(const string& key, char scalar);

   /*!
    * Get a char array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a char array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Reference to char vector to read from database.
    */
   void getArray(const string& key, vector<char>& array);

   /*!
    * Create an char array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Const reference to char vector to put into database.
    */
   void putArray(const string& key, const vector<char>& array);

   /*!
    * Get a float scalar entry in the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * float scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Reference to float scalar to read from database.
    */
   void getScalar(const string& key, float& scalar);

   /*!
    * Put a float entry into the database with the specified key name.
    *
    * @param key    Key name in database.
    * @param scalar Value to put into database.
    */
   void putScalar(const string& key, float scalar);

   /*!
    * Get a float array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a float array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Reference to float vector to read from database.
    */
   void getArray(const string& key, vector<float>& array);

   /*!
    * Create an float array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Const reference to float vector to put into database.
    */
   void putArray(const string& key, const vector<float>& array);

   /*!
    * Get a double scalar entry from the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * double scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar Reference to double scalar to read from database.
    */
   void getScalar(const string& key, double& scalar);

   /*!
    * Put a double array entry from the database with the specified key name.
    *
    * @param key    Key name in database.
    * @param scalar Double scalar to put into database.
    */
   void putScalar(const string& key, double scalar);

   /*!
    * Get a double array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a double array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Reference to double vector to read from database.
    */
   void getArray(const string& key, vector<double>& array);

   /*!
    * Create an double array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Array to put into database.
    */
   void putArray(const string& key, const vector<double>& array);

   /*!
    * Get an integer scalar entry from the database with the specified key 
    * name.
    *
    * @param key    Key name in database.
    * @param scalar Integer scalar read from database.
    */
   void getScalar(const string& key, int& scalar);

   /*!
    * Put an integer scalar entry into the database with the given key name. 
    *
    * @param key    Key name in database.
    * @param scalar Scalar integer value to put into database.
    */
   void putScalar(const string& key, int scalar);

   /*!
    * Get an integer array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a integer array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Reference to integer vector to read from database.
    */
   void getArray(const string& key, vector<int>& array);

   /*!
    * Create an integer array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Const reference to integer vector to put into database.
    */
   void putArray(const string& key, const vector<int>& array);

   /*!
    * Get a scalar string entry from the database with the specified key name.
    * If the specified key does not exist in the database or is not an
    * string scalar, then an error message is printed and the program
    * exits.
    *
    * @param key    Key name in database.
    * @param scalar String scalar value read from database.
    */
   void getScalar(const string& key, string& scalar);
 
   /*!
    * Put a string entry into the database with the specified key name.
    *
    * @param key    Key name in database.
    * @param scalar Scalar string value to put into database.
    */
   void putScalar(const string& key, const string& scalar);

   /*!
    * Get a string array entry from the database with the specified key
    * name.  If the specified key does not exist in the database or
    * is not a string array, then an error message is printed and
    * the program exits.
    *
    * @param key    Key name in database.
    * @param array  Reference to string vector to read from database.
    */
   void getArray(const string& key, vector<string>& array);
 
   /*!
    * Create an string array entry in the database with the specified
    * key name.  If the key already exists in the database, then the old
    * key record is deleted and the new one is silently created in its place.
    *
    * @param key    Key name in database.
    * @param array  Const reference to string vector to put into database.
    */
   void putArray(const string& key, const vector<string>& array);

   /*!
    * Print the current database to the specified output stream.  If
    * no output stream is specified, then data is written to stream pout.
    * 
    * @param os Output stream.
    */
   virtual void printClassData(ostream& os) = 0;

};

}

#ifndef DEBUG_NO_INLINE
#include "Database.I"
#endif
#endif
