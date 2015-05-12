//
// File:        DBModelObject.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: DB model object
//

#ifndef included_krigcpl_DBModelObject_h
#define included_krigcpl_DBModelObject_h

#include <iosfwd>

namespace krigcpl {

    /*!
     * @brief Abstraction of a model object to be used with DB.
     */

    template <typename T>
    class DBModelObject : public DBObject
    {
      
    public:
      /*!
       * @brief Constructor for the DBModelObject.
       * 
       * @param modelObject A handle to a model object to be stored on 
       *                    in a DB.
       */
      DBModelObject(const T & modelObject);

      /*!
       * @brief Destructor for the DBModelObject.
       */
      virtual ~DBModelObject();

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
       * @brief Get a copy of the model object.
       */
      T getModel() const;

    private:
      //
      // not implemented
      //
      DBModelObject(const DBModelObject &);
      const DBModelObject & operator=(const DBModelObject &);

      //
      // data
      //
      T _modelObject;

    };

}

#include "DBModelObject.I"

#endif // included_krigcpl_DBModelObject_h
