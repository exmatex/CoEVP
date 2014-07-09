//
// File:        Utilities.h
// Package:     toolbox
// 
// 
// 
// Description: A collection of trivial utility functions such as min and max
//

#ifndef included_toolbox_Utilities
#define included_toolbox_Utilities

#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_toolbox_IOStream
#include "toolbox/stream/IOStream.h"
#endif

#ifndef included_sys_types
#include <sys/types.h>
#define included_sys_types
#endif

#ifndef included_sys_stat
#include <sys/stat.h>
#define included_sys_stat
#endif

namespace toolbox {

#ifdef _MSC_VER
#include <sys/types.h>
#include <sys/stat.h>
#include <direct.h>
typedef int mode_t;
#define  S_ISDIR(m)      (((m)&S_IFMT) == S_IFDIR)
#define S_IRUSR 2
#define S_IWUSR 2
#define S_IXUSR 2
#endif

/*!
 * @brief Utilities is a utility that defines simple utility functions 
 * for issuing error messages and warnings.  This class provides a single 
 * location of definition for these functions instead of sprinkling them 
 * throughout the code.
 */

struct Utilities
{
   /*!
    * Creates the directory specified in path.  Permissions are set 
    * by default to rwx by user.  The intermediate directories in the 
    * path are created if they do not already exist.  When 
    * only_node_zero_creates is true, only node zero creates the 
    * directories.  Otherwise, all nodes create the directories.
    */
   static void recursiveMkdir(const string& path, 
			      mode_t mode = (S_IRUSR|S_IWUSR|S_IXUSR),
			      bool only_node_zero_creates = true);

   /*!
    * Rename a file from old file name to new file name.
    */
   static void renameFile(const string& old_filename, 
                          const string& new_filename);

   /*!
    * Convert an non-negative (!) integer to a string, padded on the left
    * with zeros as needed so that the string contains the number of 
    * characters indicated by the width argument.  
    * 
    * This function performs no error checking on the relationship between
    * the value of num and the number of digits indicated by width.  However,
    * when assertion checking is on, an assertion will result if num
    * is < 0, width is <= 0, or width is too small to hold all the digits for num. 
    */
   static string intToString(int num, int width);

   /*!
    * Aborts the run after printing an error message with file and
    * linenumber information.
    */
   static void abort(const string &message, 
		     const string &filename,
		     const int line);

   /*!
    * Logs warning message with file & location.
    */
   static void printWarning(const string &message, 
		            const string &filename,
		            const int line);

};


/*!
 * A statement that does nothing, for insure++ make it something 
 * more complex than a simple C null statement to avoid a warning.
 */
#ifdef __INSURE__
#define NULL_STATEMENT if(0) int nullstatement=0
#else
#define NULL_STATEMENT
#endif

/*!
 * A null use of a variable, use to avoid GNU compiler 
 * warnings about unused variables.
 */
#define NULL_USE(variable) do { \
       if(0) {char *temp = (char *)&variable; temp++;} \
    } while (0)

   /*!
    * Throw an error exception from within any C++ source code.  The 
    * macro argument may be any standard ostream expression.  The file and
    * line number of the abort are also printed.
    */
#ifndef LACKS_SSTREAM
#define TBOX_ERROR(X) do {					\
      ostringstream tboxos;					\
      tboxos << X << ends;					\
      toolbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);\
} while (0)
#else
#define TBOX_ERROR(X) do {					\
      ostrstream tboxos;					\
      tboxos << X << ends;					\
      toolbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif

   /*!
    * Print a warning without exit.  Print file and line number of the warning.
    */
#ifndef LACKS_SSTREAM
#define TBOX_WARNING(X) do {					\
      ostringstream tboxos;					\
      tboxos << X << ends;					\
      toolbox::Utilities::printWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#else
#define TBOX_WARNING(X) do {					\
      ostrstream tboxos;					\
      tboxos << X << ends;					\
      toolbox::Utilities::printWarning(tboxos.str(), __FILE__, __LINE__);\
} while (0)
#endif

   /*!
    * Throw an error exception from within any C++ source code if the
    * given expression is not true.  This is a parallel-friendly version
    * of assert.
    * The file and line number of the abort are also printed.
    */
#ifndef LACKS_SSTREAM
#define TBOX_ASSERT(EXP) do {                                   \
      if ( !(EXP) ) {                                           \
         ostringstream tboxos;                                  \
         tboxos << "Failed assertion: " << ends;\
         toolbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);\
      }                                                         \
} while (0)
#else
#define TBOX_ASSERT(EXP) do {                                    \
      if ( !(EXP) ) {                                            \
         ostrstream tboxos;                                      \
         tboxos << "Failed assertion: " << ends;                 \
         toolbox::Utilities::abort(tboxos.str(), __FILE__, __LINE__);\
      }                                                          \
} while (0)
#endif


}

#endif
