//
// File:	PIO.h
// Package:	toolbox
// 
// 
// 
// Description:	Parallel I/O classes pout, perr, and plog and control class
//

#ifndef included_toolbox_PIO
#define included_toolbox_PIO

#ifndef included_fstream
#include <fstream>
#include <iostream>
#define included_fstream
using namespace std;
#endif

#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

namespace toolbox {

/*!
 * @brief PIO manages parallel stream I/O and logging.  Static member
 * function initialize() must be called before any of the parallel streams
 * pout, perr, or plog may be used.  Routine finalize() should also be called
 * before termination of the program.  
 *
 * By default, logging is disabled.  To enable logging, call one of the
 * routines logOnlyNodeZero() or logAllNodes().  Logging may be suspended
 * and resumed.
 */

struct PIO
{
   /*!
    * Initialize the parallel I/O streams.  This routine must be called
    * before using pout, perr, or plog.  This routine
    * must be called after the MPI routines have been initialized; e.g.,
    * after calling MPI::init().
    */
   static void initialize();

   /*!
    * Shut down the parallel I/O streams and close log files.  This routine
    * must be called before program termination.
    */
   static void finalize();

   /*!
    * Log messages for node zero only to the specified filename.  All output
    * to pout, perr, and plog on node zero will go to the log file.
    */
   static void logOnlyNodeZero(const string &filename);

   /*!
    * Log messages from all nodes.  The diagnostic data for processor XXXXX
    * will be sent to a file with the name filename.XXXXX, where filename is
    * the function argument.
    */
   static void logAllNodes(const string &filename);

   /*!
    * Temporarily suspend file logging.  Log file output will be discarded,
    * although the output file will not be closed.  Logging can be resumed
    * by calling member function resumeLogging().
    */
   static void suspendLogging();

   /*!
    * Resume logging after logging was suspended via member function
    * suspendLogging().
    */
   static void resumeLogging();

private:
   static void shutdownFilestream();	// shutdown the log filestream

   static int       s_rank;		// processor rank in MPI group
   static ofstream* s_filestream;	// NULL or log filestream
};

/*!
 * Parallel output stream pout writes to the standard output from node zero
 * only.  Output from other nodes is ignored.  If logging is enabled, then
 * output is mirrored to the log stream, as well.
 */
extern ostream pout;

/*!
 * Parallel output stream perr writes to the standard error from all nodes.
 * Output is prepended with the processor number.  If logging is enabled,
 * then output is mirrored to the log stream, as well.
 */
extern ostream perr;

/*!
 * Parallel output stream plog writes output to the log file.  When logging
 * from multiple processors, the processor number is appended to the filename.
 */
extern ostream plog;


}

#endif
