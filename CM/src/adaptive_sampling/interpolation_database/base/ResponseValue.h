#ifndef included_krigalg_ResponseValue
#define included_krigalg_ResponseValue

#ifndef included_krigalg_Value
#include <base/Value.h>
#endif

namespace krigcpl {

/*!
 * @brief A ResponseValue object represents the value of a response calculated at a
 * ResponsePoint. Most ResponsePoint functionality is inherited from the krigalg::Value class.
 *
 * @see krigalg::Value
 */

class ResponseValue : krigalg::Value
{
public:
   /*!
    * Ctor for a response value object indicating dimension of response value..
    */
   explicit ResponseValue(int dimension);

   /*!
    * Virtual destructor for response value objects.
    */
   virtual ~ResponseValue();

   /*!
    * Print response point object data to the specified output stream.
    */
   friend ostream& operator<< (ostream& stream,
                               const ResponseValue& value);

private:
   // These functions are not implemented, but are declared here to avoid
   // having the compiler provide default implementations.
   ResponseValue();
   ResponseValue(const ResponseValue&);
   void operator=(const ResponseValue&);

};

}
#ifndef DEBUG_NO_INLINE
#include "ResponseValue.I"
#endif
#endif
