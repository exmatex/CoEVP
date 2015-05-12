//
// File:        DBKrigingModelObjectFactory.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: DB model object factory
//

#ifndef included_krigcpl_DBKrigingModelObjectFactory_h
#define included_krigcpl_DBKrigingModelObjectFactory_h

#include "kriging/MultivariateDerivativeKrigingModel.h"

#include <base/DBObject.h>
#include <base/DBObjectFactory.h>

namespace krigcpl {

    /*!
     * @brief Abstraction of a model object factory fo the use with generic DB
     */

    template <typename T>
    class DBKrigingModelObjectFactory : public DBObjectFactory {
      
    public:
      
      /*!
       * @brief Constructor.
       */

      DBKrigingModelObjectFactory();
      
      /*!
       * @brief Destructor
       */

      virtual ~DBKrigingModelObjectFactory();

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
      DBKrigingModelObjectFactory(const DBKrigingModelObjectFactory&);
      void operator=(const DBKrigingModelObjectFactory&);

    };

    //
    // template specializations
    //
    template<> DBObjectPtr 
      DBKrigingModelObjectFactory<krigalg::MultivariateDerivativeKrigingModel>::allocateObject(toolbox::Database& db) const;

}

#include "DBKrigingModelObjectFactory.I"

#endif // included_krigcpl_DBKrigingModelObjectFactory_h
