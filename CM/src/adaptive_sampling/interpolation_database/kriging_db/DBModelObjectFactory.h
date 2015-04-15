//
// File:        DBModelObjectFactory.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: DB model object factory
//

#ifndef included_krigcpl_DBModelObjectFactory_h
#define included_krigcpl_DBModelObjectFactory_h

#include "kriging/MultivariateDerivativeKrigingModel.h"

#include <base/DBObject.h>
#include <base/DBObjectFactory.h>

namespace krigcpl {

    /*!
     * @brief Abstraction of a model object factory fo the use with generic DB
     */

    template <typename T>
    class DBModelObjectFactory : public DBObjectFactory {
      
    public:
      
      /*!
       * @brief Constructor.
       */

      DBModelObjectFactory();
      
      /*!
       * @brief Destructor
       */

      virtual ~DBModelObjectFactory();

      /*!
       * @brief Allocate object and fill its contents from the database.
       *
       * @param db Handle to a database.
       *
       * @return Pointer to DBObject.
       */
      
      DBObjectPtr allocateObject(toolbox::Database& db) const;

    private:
      // The following are not implemented
      DBModelObjectFactory(const DBModelObjectFactory&);
      void operator=(const DBModelObjectFactory&);

    };

    //
    // template specializations
    //
    template<> DBObjectPtr 
      DBModelObjectFactory<krigalg::MultivariateDerivativeKrigingModel>::allocateObject(toolbox::Database& db) const;

}

#include "DBModelObjectFactory.I"

#endif // included_krigcpl_DBModelObjectFactory_h
