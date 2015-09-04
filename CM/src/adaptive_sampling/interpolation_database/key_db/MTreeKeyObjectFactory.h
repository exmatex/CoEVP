//
// File:        MTreeKeyObjectFactory.h
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: DB model object factory
//

#ifndef included_krigcpl_MTreeKeyObjectFactory_h
#define included_krigcpl_MTreeKeyObjectFactory_h

#include <base/DBObject.h>
#include <base/DBObjectFactory.h>

namespace krigcpl {

    /*!
     * @brief Abstraction of a model object factory fo the use with a key DB
     */

   template <typename T>
      class MTreeKeyObjectFactory : public DBObjectFactory {
      
    public:
      
      /*!
       * @brief Constructor.
       */

      MTreeKeyObjectFactory();
      
      /*!
       * @brief Destructor
       */

      virtual ~MTreeKeyObjectFactory();

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
      MTreeKeyObjectFactory(const MTreeKeyObjectFactory&);
      void operator=(const MTreeKeyObjectFactory&);

    };

    //
    // template specializations
    //
   template<> DBObjectPtr 
      MTreeKeyObjectFactory<std::string>::allocateObject(toolbox::Database& db) const;

}

#include "MTreeKeyObjectFactory.I"

#endif // included_krigcpl_MTreeKeyObjectFactory_h
