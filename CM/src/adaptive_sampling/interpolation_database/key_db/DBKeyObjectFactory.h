//
// File:        DBKeyObjectFactory.h
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: DB model object factory
//

#ifndef included_krigcpl_DBKeyObjectFactory_h
#define included_krigcpl_DBKeyObjectFactory_h

#include <base/DBObject.h>
#include <base/DBObjectFactory.h>

namespace krigcpl {

    /*!
     * @brief Abstraction of a model object factory fo the use with a key DB
     */

    template <typename T>
    class DBKeyObjectFactory : public DBObjectFactory {
      
    public:
      
      /*!
       * @brief Constructor.
       */

      DBKeyObjectFactory();
      
      /*!
       * @brief Destructor
       */

      virtual ~DBKeyObjectFactory();

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
      DBKeyObjectFactory(const DBKeyObjectFactory&);
      void operator=(const DBKeyObjectFactory&);

    };

    //
    // template specializations
    //
    template<> DBObjectPtr 
       DBKeyObjectFactory<std::string>::allocateObject(toolbox::Database& db) const;

}

#include "DBKeyObjectFactory.I"

#endif // included_krigcpl_DBKeyObjectFactory_h
