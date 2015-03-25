//
// File:        MTreeModelObjectFactory.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: MTree model object factory
//

#ifndef included_krigcpl_MTreeModelObjectFactory_h
#define included_krigcpl_MTreeModelObjectFactory_h

#include "kriging/MultivariateDerivativeKrigingModel.h"

#include <mtreedb/MTreeObject.h>
#include <mtreedb/MTreeObjectFactory.h>

namespace krigcpl {

    /*!
     * @brief Abstraction of a model object factory fo the use with
     * MTree DB.
     */

    template <typename T>
    class MTreeModelObjectFactory : public mtreedb::MTreeObjectFactory {
      
    public:
      
      /*!
       * @brief Constructor.
       */

      MTreeModelObjectFactory();
      
      /*!
       * @brief Destructor
       */

      virtual ~MTreeModelObjectFactory();

      /*!
       * @brief Allocate object and fill its contents from the database.
       *
       * @param db Handle to a database.
       *
       * @return Pointer to MTreeObject.
       */
      
      mtreedb::MTreeObjectPtr allocateObject(toolbox::Database& db) const;

    private:
      // The following are not implemented
      MTreeModelObjectFactory(const MTreeModelObjectFactory&);
      void operator=(const MTreeModelObjectFactory&);

    };

    //
    // template specializations
    //
    template<> mtreedb::MTreeObjectPtr 
      MTreeModelObjectFactory<krigalg::MultivariateDerivativeKrigingModel>::allocateObject(toolbox::Database& db) const;

}

#include "MTreeModelObjectFactory.I"

#endif // included_krigcpl_MTreeModelObjectFactory_h
