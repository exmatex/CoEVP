//
// File:	IOStream.h
// Package:	toolbox
// 
// 
// 
// Description:	Wrapper header file for standard IO stream classes
//

#ifndef included_Stream
#define included_Stream

#ifndef included_stdio
#define included_stdio
#include <stdio.h>
#endif

#ifndef included_iostream
#define included_iostream
#include <iostream>
#endif

#ifndef included_iomanip
#define included_iomanip
#include <iomanip>
#endif

#ifndef LACKS_SSTREAM
#ifndef included_sstream
#define included_sstream
#include <sstream>
#endif
#endif

#if 0   // strssream is deprecated, and we don't seem to need it anyway
#ifndef LACKS_STRSTREAM
#ifndef included_strstream
#define included_strstream
#include <strstream>
#endif
#endif
#endif
using namespace std;


#endif
