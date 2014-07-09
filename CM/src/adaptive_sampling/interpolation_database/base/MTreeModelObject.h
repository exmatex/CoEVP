//
// File:        MTreeModelObject.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: MTree model object
//

#ifndef included_krigcpl_MTreeModelObject_h
#define included_krigcpl_MTreeModelObject_h

#include <iosfwd>

namespace krigcpl {

    /*!
     * @brief Abstraction of a model object to be used with MTree DB.
     */

    template <typename T>
    class MTreeModelObject : public mtreedb::MTreeObject
    {
      
    public:
      /*!
       * @brief Constructor for the MTreeModelObject.
       * 
       * @param modelObject A handle to a model object to be stored on 
       *                    in an MTree DB.
       */
      MTreeModelObject(const T & modelObject);

      /*!
       * @brief Destructor for the MTreeModelObject.
       */
      virtual ~MTreeModelObject();

      /*!
       * @brief Concrete virtual method to create and return smart
       * pointer to a (deep) copy of this object.
       */
      mtreedb::MTreeObjectPtr makeCopy() const;

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
      MTreeModelObject(const MTreeModelObject &);
      const MTreeModelObject & operator=(const MTreeModelObject &);

      //
      // data
      //
      T _modelObject;

    };

}

#include "MTreeModelObject.I"

#endif // included_krigcpl_MTreeModelObject_h
