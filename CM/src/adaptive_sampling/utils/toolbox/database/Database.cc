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
// File:        Database.cc
// Package:     toolbox
// 
// 
// 
// Description: An abstract base class for managing database objects in files
//

#include "toolbox/database/Database.h"

#ifdef DEBUG_NO_INLINE
#include "toolbox/database/Database.I"
#endif

namespace toolbox {

Database::~Database()
{
}

/*  
 * Scalar/array functions for boolean type.
 */ 

void Database::getScalar(const string& key, bool& scalar)
{
   scalar = getBool(key);
}

void Database::putScalar(const string& key, bool scalar)
{
   putBool(key, scalar);
}

void Database::getArray(const string& key, vector<bool>& array)
{
   getBoolArray(key, array);
}

void Database::putArray(const string& key, const vector<bool>& array)
{
   putBoolArray(key, array);
}

/*  
 * Scalar/array functions for char type.
 */ 

void Database::getScalar(const string& key, char& scalar)
{
   scalar = getChar(key);
}

void Database::putScalar(const string& key, char scalar)
{
   putChar(key, scalar);
}

void Database::getArray(const string& key, vector<char>& array)
{
   getCharArray(key, array);
}

void Database::putArray(const string& key, const vector<char>& array)
{
   putCharArray(key, array);
}

/*  
 * Scalar/array functions for float type.
 */ 

void Database::getScalar(const string& key, float& scalar)
{
   scalar = getFloat(key);
}

void Database::putScalar(const string& key, float scalar)
{
   putFloat(key, scalar);
}

void Database::getArray(const string& key, vector<float>& array)
{
   getFloatArray(key, array);
}

void Database::putArray(const string& key, const vector<float>& array)
{
   putFloatArray(key, array);
}

/*  
 * Scalar/array functions for double type.
 */ 

void Database::getScalar(const string& key, double& scalar)
{
   scalar = getDouble(key);
}

void Database::putScalar(const string& key, double scalar)
{
   putDouble(key, scalar);
}

void Database::getArray(const string& key, vector<double>& array)
{
   getDoubleArray(key, array);
}

void Database::putArray(const string& key, const vector<double>& array)
{
   putDoubleArray(key, array);
}

/*  
 * Scalar/array functions for integer type.
 */ 

void Database::getScalar(const string& key, int& scalar)
{
   scalar = getInteger(key);
}

void Database::putScalar(const string& key, int scalar)
{
   putInteger(key, scalar);
}

void Database::getArray(const string& key, vector<int>& array)
{
   getIntegerArray(key, array);
}

void Database::putArray(const string& key, const vector<int>& array)
{
   putIntegerArray(key, array);
}

/* 
 * Scalar/array functions for string type.
 */
 
void Database::getScalar(const string& key, string& scalar)
{
   scalar = getString(key);
}
 
void Database::putScalar(const string& key, const string& scalar) 
{
   putString(key, scalar);
}
 
void Database::getArray(const string& key, vector<string>& array)
{
   getStringArray(key, array); 
}
 
void Database::putArray(const string& key, const vector<string>& array)
{
   putStringArray(key, array);
}


}




