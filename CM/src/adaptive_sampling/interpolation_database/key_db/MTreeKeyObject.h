//
// File:        MTreeKeyObject.h
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: DB key object
//

#ifndef included_krigcpl_MTreeKeyObject_h
#define included_krigcpl_MTreeKeyObject_h

#include <mtreedb/MTreeObject.h>

namespace krigcpl {

    /*!
     * @brief Abstraction of a key object to be used with DB.
     */

    template <typename T>
    class MTreeKeyObject : public MTreeObject
    {
      
    public:

      MTreeKeyObject() {;}

      /*!
       * @brief Constructor for the MTreeKeyObject.
       * 
       * @param keyObject A handle to a key object to be stored on 
       *                    in a DB.
       */
      MTreeKeyObject(const T & keyObject);

      /*!
       * @brief Destructor for the MTreeKeyObject.
       */
      virtual ~MTreeKeyObject();

      /*!
       * @brief Concrete virtual method to create and return smart
       * pointer to a (deep) copy of this object.
       */
      MTreeObjectPtr makeCopy() const;

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
      MTreeKeyObject(const MTreeKeyObject &);
      const MTreeKeyObject & operator=(const MTreeKeyObject &);

      //
      // data
      //
      T _keyObject;

    };

}

#include "MTreeKeyObject.I"

#endif // included_krigcpl_MTreeKeyObject_h
