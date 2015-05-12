//
// File:        DBKeyObject.h
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: DB key object
//

#ifndef included_krigcpl_DBKeyObject_h
#define included_krigcpl_DBKeyObject_h

#include <base/DBObject.h>

namespace krigcpl {

    /*!
     * @brief Abstraction of a key object to be used with DB.
     */

    template <typename T>
    class DBKeyObject : public DBObject
    {
      
    public:

      DBKeyObject() {;}

      /*!
       * @brief Constructor for the DBKeyObject.
       * 
       * @param keyObject A handle to a key object to be stored on 
       *                    in a DB.
       */
      DBKeyObject(const T & keyObject);

      /*!
       * @brief Destructor for the DBKeyObject.
       */
      virtual ~DBKeyObject();

      /*!
       * @brief Concrete virtual method to create and return smart
       * pointer to a (deep) copy of this object.
       */
      DBObjectPtr makeCopy() const;

      /*!
       * @brief Concrete virtual method to write data members to given
       * database.
       */
      void putToDatabase(toolbox::Database & db) const;

      /*!
       * @brief Vitual method to print concrete object data to the specified
       * output stream.  This is optional; a default no-op is supplied
       * here.
       */
      ostream & print(ostream & stream) const;
      
      /*!
       * @brief Get a copy of the key object.
       */
      T getKey() const;

    private:
      //
      // not implemented
      //
      DBKeyObject(const DBKeyObject &);
      const DBKeyObject & operator=(const DBKeyObject &);

      //
      // data
      //
      T _keyObject;

    };

}

#include "DBKeyObject.I"

#endif // included_krigcpl_DBKeyObject_h
