// DO-NOT-DELETE revisionify.begin() 
/*

                            Copyright (c) 2014.
               Lawrence Livermore National Security, LLC.
         Produced at the Lawrence Livermore National Laboratory
                             LLNL-CODE-656392.
                           All rights reserved.

This file is part of CoEVP, Version 1.0. Please also read this link -- http://www.opensource.org/licenses/index.php

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/
// DO-NOT-DELETE revisionify.end() 
//
// File:	PIO.cc
// Package:	toolbox
// 
// 
// 
// Description:	Parallel I/O classes pout, perr, and plog and control class
//

#include "toolbox/base/PIO.h"

#include <stdio.h>
#include "toolbox/parallel/MPI.h"
#include "toolbox/base/ParallelBuffer.h"

#ifndef NULL
#define NULL 0
#endif

namespace toolbox {

int       PIO::s_rank       = -1;
ofstream* PIO::s_filestream = NULL;

/*
*************************************************************************
*									*
* Define the parallel buffers and the associated ostream objects.	*
*									*
*************************************************************************
*/

static ParallelBuffer pout_buffer;
static ParallelBuffer perr_buffer;
static ParallelBuffer plog_buffer;

ostream pout(&pout_buffer);
ostream perr(&perr_buffer);
ostream plog(&plog_buffer);

/*
*************************************************************************
*									*
* Initialie the parallel I/O streams.  This routine must be called	*
* before pout, perr, and plog are used for output but after MPI has	*
* been initialized.  By default, logging is disabled.			*
*									*
*************************************************************************
*/

void PIO::initialize()
{
   s_rank       = MPI::getRank();
   s_filestream = NULL;
   
   /*
    * Initialize the standard parallel output stream
    */

   pout_buffer.setActive(s_rank == 0);
   pout_buffer.setPrefixString(string());
   pout_buffer.setOutputStream1(&cout);
   pout_buffer.setOutputStream2(NULL);

   /*
    * Initialize the error parallel output stream
    */

   char buffer[16];
   sprintf(buffer, "P=%05d:", s_rank);

   perr_buffer.setActive(true);
   perr_buffer.setPrefixString(buffer);
   perr_buffer.setOutputStream1(&cerr);
   perr_buffer.setOutputStream2(NULL);

   /*
    * Initialize the parallel log file (disabled by default)
    */

   plog_buffer.setActive(false);
   plog_buffer.setPrefixString(string());
   plog_buffer.setOutputStream1(NULL);
   plog_buffer.setOutputStream2(NULL);
}

/*
*************************************************************************
*									*
* Close the output streams.  Flush both cout and cerr.  If logging,	*
* then flush and close the log stream.					*
*									*
*************************************************************************
*/

void PIO::finalize()
{
   cout.flush();
   cerr.flush();
   shutdownFilestream();
}

/*
*************************************************************************
*									*
* If the log file stream is open, then shut down the filestream.  Close	*
* and flush the channel and disconnect the output stream buffers.	*
*									*
*************************************************************************
*/

void PIO::shutdownFilestream()
{
   if (s_filestream) {
      s_filestream->flush();
      s_filestream->close();

      delete s_filestream;
      s_filestream = NULL;

      pout_buffer.setOutputStream2(NULL);
      perr_buffer.setOutputStream2(NULL);
      plog_buffer.setOutputStream1(NULL);
      plog_buffer.setActive(false);
   }
}

/*
*************************************************************************
*									*
* Log messages for node zero only.  If a log stream was open, close	*
* it.  If this is node zero, then open a new log stream and set the	*
* appropriate buffer streams to point to the log file.			*
*									*
*************************************************************************
*/

void PIO::logOnlyNodeZero(const string &filename)
{
   /*
    * If the filestream was open, then close it and reset streams
    */

   shutdownFilestream();

   /*
    * If this is node zero, then open the log stream and redirect output
    */

   if (s_rank == 0) {
      s_filestream = new ofstream(filename.c_str());
      if (!(*s_filestream)) {
         delete s_filestream;
         s_filestream = NULL;
         perr << "PIO: Could not open log file ``" << filename.c_str() << "''\n";
      } else {
         pout_buffer.setOutputStream2(s_filestream);
         perr_buffer.setOutputStream2(s_filestream);
         plog_buffer.setOutputStream1(s_filestream);
         plog_buffer.setActive(true);
      }
   }
}

/*
*************************************************************************
*									*
* Log messages for all nodes.  If a log stream was open, the close it.	*
* Open a log stream on every processor.  The filename for the log file	*
* will be appended with the processor number.				*
*									*
*************************************************************************
*/

void PIO::logAllNodes(const string &filename)
{
   /*
    * If the filestream was open, then close it and reset streams
    */

   shutdownFilestream();

   /*
    * Open the log stream and redirect output
    */

   char *buffer = new char[filename.length() + 16];
   sprintf(buffer, "%s.%05d", filename.c_str(), s_rank);
   s_filestream = new ofstream(buffer);

   if (!(*s_filestream)) {
      delete s_filestream;
      s_filestream = NULL;
      perr << "PIO: Could not open log file ``" << buffer << "''\n";
   } else {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }

   delete [] buffer;
}

/*
*************************************************************************
*									*
* Suspend logging of data to the file stream.  This does not close the	*
* filestream (assuming it is open) but just disables logging.		*
*									*
*************************************************************************
*/

void PIO::suspendLogging()
{
   pout_buffer.setOutputStream2(NULL);
   perr_buffer.setOutputStream2(NULL);
   plog_buffer.setOutputStream1(NULL);
   plog_buffer.setActive(false);
}

/*
*************************************************************************
*									*
* Resume logging of the file stream (assuming it was open).  If the	*
* file stream is NULL, then do nothing.					*
*									*
*************************************************************************
*/

void PIO::resumeLogging()
{
   if (s_filestream) {
      pout_buffer.setOutputStream2(s_filestream);
      perr_buffer.setOutputStream2(s_filestream);
      plog_buffer.setOutputStream1(s_filestream);
      plog_buffer.setActive(true);
   }
}


}



