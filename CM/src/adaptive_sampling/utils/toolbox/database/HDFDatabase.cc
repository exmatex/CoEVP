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
// File:        HDFDatabase.cc
// Package:     toolbox
// 
// 
// 
// Description: A database structure that stores HDF5 format data.
//

#include "toolbox/database/HDFDatabase.h"

#ifdef HAVE_PKG_hdf5

#include "toolbox/stream/IOStream.h"
#include "toolbox/base/Utilities.h"
#include "toolbox/base/MathUtilities.h"

#include <string.h>

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif
#endif

/*
*************************************************************************
*                                                                       *
* Integer keys for identifying types in HDF5 database.  Negative        *
* entries are used to distinguish arrays from scalars when printing     *
* key information.                                                      *
*                                                                       *
*************************************************************************
*/
#define KEY_DATABASE        (0)
#define KEY_BOOL_ARRAY      (1)
#define KEY_CHAR_ARRAY      (3)
#define KEY_DOUBLE_ARRAY    (5)
#define KEY_FLOAT_ARRAY     (6)
#define KEY_INT_ARRAY       (7)
#define KEY_STRING_ARRAY    (8)

#define KEY_BOOL_SCALAR     (-1)
#define KEY_CHAR_SCALAR     (-3)
#define KEY_DOUBLE_SCALAR   (-5)
#define KEY_FLOAT_SCALAR    (-6)
#define KEY_INT_SCALAR      (-7)
#define KEY_STRING_SCALAR   (-8)


/*
  Macros starting with H5T_MPTCOUPLER_ are for controlling the data
  type that is actually written to the file.  As long as
  these are not "native" types, the file should be portable.
*/

// Type used for writing simple (non-compound) data.
#define H5T_MPTCOUPLER_INT      H5T_STD_I32BE
#define H5T_MPTCOUPLER_FLOAT    H5T_IEEE_F32BE
#define H5T_MPTCOUPLER_DOUBLE   H5T_IEEE_F64BE
#define H5T_MPTCOUPLER_BOOL     H5T_STD_I8BE

// Type used for writing the data attribute key.
#define H5T_MPTCOUPLER_ATTR H5T_STD_I8BE


/*
*************************************************************************
*                                                                       *
* Macros to suppress the HDF5 messages sent to standard i/o; handle     *
* errors explicity within this code.                                    *
*                                                                       *
*************************************************************************
*/

#define BEGIN_SUPPRESS_HDF5_WARNINGS                   \
{                                                     \
   herr_t (*H5E_saved_efunc) (void*);                 \
   void *H5E_saved_edata;                             \
   H5Eget_auto(&H5E_saved_efunc, &H5E_saved_edata);   \
   H5Eset_auto(NULL, NULL);                              

#define END_SUPPRESS_HDF5_WARNINGS                     \
   H5Eset_auto(H5E_saved_efunc, H5E_saved_edata);     \
}



/*
*************************************************************************
* We may wish to assert HDF5 return values regardless of debug modes.   *
*************************************************************************
*/
#define ASSERT_HDF5_RETURN_VALUES
#ifdef ASSERT_HDF5_RETURN_VALUES
#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif
#endif



namespace toolbox {

string HDFDatabase::s_top_level_search_group = string();
string HDFDatabase::s_group_to_search = string();
int HDFDatabase::s_still_searching = 0;
int HDFDatabase::s_found_group = 0;

/*
*************************************************************************
*                                                                       *
* Static member function to iterate through the hdf5 data file and      *
* assemble a list of desired (key, type) pairs.                         *
*                                                                       *
*************************************************************************
*/

herr_t HDFDatabase::iterateKeys(
   hid_t loc_id,
   const char* name,
   void* opdata)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(name != (char*)NULL);
#endif

   if (s_still_searching) {

     H5G_stat_t statbuf;
     int type_key;
     herr_t errf;

     errf = H5Gget_objinfo(loc_id, name, 0, &statbuf);
#ifdef ASSERT_HDF5_RETURN_VALUES
     assert( errf >= 0 );
#endif

     switch (statbuf.type) {
     case H5G_GROUP: {
       if (s_top_level_search_group == "/") {
         addKeyToList(name, KEY_DATABASE, opdata);
       } else if ( !strcmp(name, s_group_to_search.c_str()) ) {
         hid_t grp;
         grp = H5Gopen(loc_id, name);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( grp >= 0 );
#endif
         s_found_group = true;
         s_still_searching =
           H5Giterate(grp, ".", NULL,
                      HDFDatabase::iterateKeys, opdata);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( s_still_searching >= 0 );
#endif
         s_found_group = false;
       } else {
         hid_t grp;
         grp = H5Gopen(loc_id, name);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( grp >= 0 );
#endif
         if (s_found_group) {
           addKeyToList(name, KEY_DATABASE, opdata);
         } else {
           errf = H5Giterate(grp, ".", NULL,
                             HDFDatabase::iterateKeys, opdata);
#ifdef ASSERT_HDF5_RETURN_VALUES
           assert( errf >= 0 );
#endif
         }
       }
       break;
     }

     case H5G_DATASET: {
       if (s_still_searching && s_found_group) {
         hid_t this_set;
         BEGIN_SUPPRESS_HDF5_WARNINGS
         this_set = H5Dopen(loc_id, name);
         END_SUPPRESS_HDF5_WARNINGS
         if (this_set > 0) {
            hid_t attr = H5Aopen_name(this_set, "Type");
#ifdef ASSERT_HDF5_RETURN_VALUES
            assert( attr >= 0 );
#endif
            errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
#ifdef ASSERT_HDF5_RETURN_VALUES
            assert( errf >= 0 );
#endif
            hid_t this_space = H5Dget_space(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
            assert( this_space >= 0 );
#endif
            hsize_t nsel = H5Sget_select_npoints(this_space);
            int array_size = int(nsel); 
            addKeyToList(name,
                         (array_size == 1 ? -type_key : type_key),
                         opdata);
            errf = H5Sclose(this_space);
#ifdef ASSERT_HDF5_RETURN_VALUES
            assert( errf >= 0 );
#endif
            errf = H5Aclose(attr);
#ifdef ASSERT_HDF5_RETURN_VALUES
            assert( errf >= 0 );
#endif
            errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
            assert( errf >= 0 );
#endif
         }
       }
       break;
     }

     default: {
       TBOX_ERROR("HDFDatabase key search error....\n"
                  << "   Unable to identify key = " << name 
                  << " as a known group or dataset" << endl);
     }
     }

   }
   return 0;
}

/*
*************************************************************************
*                                                                       *
* Static member function to add key to list for database associated     *
* with void* argument.                                                  *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::addKeyToList(
   const char* name,
   int type,
   void* database) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(name != (char*)NULL);
   assert(database != NULL);
#endif

   KeyData key_item;
   key_item.d_key  = name; 
   key_item.d_type = type; 

   ((HDFDatabase*)database)->d_keydata.push_back(key_item);
}

/*
*************************************************************************
*                                                                       *
* Public HDF database constructor creates an empty database with the    *
* specified name.  It sets the group_ID to a default value of -1.       *
* This data is used by member functions to track parent databases.      *
*                                                                       *
*************************************************************************
*/

HDFDatabase::HDFDatabase(const string& name) :
   d_is_file(false),
   d_file_id(-1),
   d_group_id(-1),
   d_database_name(name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
#endif
   d_keydata.clear(); 
}

/*
*************************************************************************
*                                                                       *
* Private HDF database constructor creates an empty database with the   *
* specified name.  The group_ID is used privately within                *
* the member functions to track parent databases.                       *
*                                                                       *
*************************************************************************
*/

HDFDatabase::HDFDatabase(
   const string& name, 
   hid_t group_ID) :
   d_is_file(false),
   d_file_id(-1),
   d_group_id(group_ID),
   d_database_name(name)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!name.empty());
#endif
   d_keydata.clear(); 
}

/*
*************************************************************************
*                                                                       *
* The database destructor closes the opened file or group.              *
*                                                                       *
*************************************************************************
*/

HDFDatabase::~HDFDatabase()
{
   herr_t errf;
   if (d_is_file) {
      unmount();
   } 
   if ( d_group_id != -1 ) {
      errf = H5Gclose(d_group_id);
#ifdef DEBUG_CHECK_ASSERTIONS
      assert( errf >= 0 );
#endif
   }
}

/*
*************************************************************************
*                                                                       *
* Return true if the key exists within the database; false otherwise.   *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::keyExists(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   bool key_exists = false;
   herr_t errf;
   
   hid_t this_set;
   BEGIN_SUPPRESS_HDF5_WARNINGS
   this_set = H5Dopen(d_group_id, key.c_str());
   END_SUPPRESS_HDF5_WARNINGS
   if (this_set > 0) {
      key_exists = true;
      errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
   }
   if (!key_exists) {
      hid_t this_group;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_group = H5Gopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_group > 0) {
         key_exists = true;
         errf = H5Gclose(this_group);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return key_exists;
}

/*
*************************************************************************
*                                                                       *
* Return all keys in the database.                                      *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::getAllKeys(vector<string>& keys)
{
   performKeySearch();
 
   keys.clear();
   keys.reserve(d_keydata.size());
 
   for (list<KeyData>::iterator ik = d_keydata.begin();
        ik != d_keydata.end(); ++ik) {
      keys.push_back((*ik).d_key);
   }
 
   cleanupKeySearch();
}

/*
*************************************************************************
*                                                                       *
* Return the size of the array associated with the key.  If the key     *
* does not exist, then zero is returned.                                *
* Array size is set based on the number of elements (points) within     *
* the dataspace defined by the named dataset (or key).                  *
*                                                                       *
*************************************************************************
*/

int HDFDatabase::getArraySize(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   herr_t errf;
   int array_size  = 0;

   hid_t this_set;
   BEGIN_SUPPRESS_HDF5_WARNINGS
   this_set = H5Dopen(d_group_id, key.c_str());
   END_SUPPRESS_HDF5_WARNINGS
   if (this_set > 0) {
      hid_t this_space = H5Dget_space(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( this_space >= 0 );
#endif
      hsize_t nsel = H5Sget_select_npoints(this_space);
      array_size = int(nsel);
      errf = H5Sclose(this_space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
   }

   return array_size;
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a database entry.  If the key does not exist, then false   *
* is returned.  The key represents a database (or hdf group) if the     *
* H5Gopen function on the key is successful.                            *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isDatabase(const string& key)
{
   bool is_database = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_group;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_group = H5Gopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_group > 0) {
         is_database = true;
         errf = H5Gclose(this_group);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return is_database;
}

/*
*************************************************************************
*                                                                       *
* Create a new database with the specified key name.                    *
*                                                                       *
*************************************************************************
*/

DatabasePtr HDFDatabase::putDatabase(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif

   //
   // check if a group exists and unlink it if it does
   //

   if (isDatabase(key.c_str()) == true) {

     herr_t unlink_return = H5Gunlink(d_group_id, key.c_str());
     
     if (unlink_return <  0)
       TBOX_ERROR("HDFDatabase::putDatabase() error in database "
		  << d_database_name << endl
		  << "could not unlink " << key << endl);
     
   }
   
   //
   //
   //

   string parent_name = d_database_name;
   hid_t  parent_id   = d_group_id;

   hid_t this_group = H5Gcreate(parent_id, key.c_str(), 0);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( this_group >= 0 );

#endif
   DatabasePtr new_database( new HDFDatabase(key, this_group) );

   return(new_database);
}

/*
************************************************************************
*                                                                      *
* Get the database with the specified key name.                        *
*                                                                      *
************************************************************************
*/

DatabasePtr HDFDatabase::getDatabase(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if (!isDatabase(key)) {
      TBOX_ERROR("HDFDatabase::getDatabase() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a database." << endl);
   }

   hid_t this_group = H5Gopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( this_group >= 0 );

#endif
   DatabasePtr database( new HDFDatabase(key, this_group) );

   return(database);
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a boolean entry.  If the key does not exist, then false    *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isBool(const string& key)
{
   bool is_boolean  = false;
   herr_t errf;
   
   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_BOOL_ARRAY) {
            is_boolean = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return is_boolean;
}

/*
*************************************************************************
*                                                                       *
* Create a boolean scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putBool(
   const string& key, 
   bool data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   putBoolArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a boolean array entry in the database with the specified       *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putBoolArray(
   const string& key, 
   const vector<bool>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif

   TBOX_ERROR("HDFDatabase::putBoolArray() error... "
              << "\n  Function not implemented properly and"
              << " needs to be fixed!" << endl);

#if 0
   if ( data.size() > 0 ) {
      putBoolArray(key, &data[0], data.size());
   } else {
      TBOX_ERROR("HDFDatabase::putBoolArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
#endif
}

/*
*************************************************************************
*                                                                       *
* Create a boolean array entry in the database with the specified       *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putBoolArray(
   const string& key, 
   const bool* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != (bool*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( space >= 0 );
#endif

      /*
        We cannot be sure exactly what bool is because it is
        represented differently on different platforms, and
        it may have been redefined, i.e., by the Boolean
        type.  We are unsure what the bool is so we convert it
        to the native int type (H5T_NATIVE_INT) before giving
        it to HDF.  When we write a bool, we write it the
        shortest integer type we can find, the H5T_MPTCOUPLER_BOOL
        type.
      */
      vector<int> data1( nelements );
      for ( int i=0; i<nelements; ++i ) data1[i] = data[i];

      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_MPTCOUPLER_BOOL,
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                      H5P_DEFAULT, &data1[0]);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_BOOL_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif

   } else {
      TBOX_ERROR("HDFDatabase::putBoolArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get boolean scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a boolean type.                                                 *
*                                                                      *
************************************************************************
*/

bool HDFDatabase::getBool(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   bool ret_val;
   getBoolArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Two routines to get boolean arrays from the database with the        *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a boolean type.                   *
*                                                                      *
************************************************************************
*/

void HDFDatabase::getBoolArray(const string& key, 
                               vector<bool>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if (!isBool(key)) {
      TBOX_ERROR("HDFDatabase::getBoolArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a bool array." << endl);
   }

   hid_t dset, dspace;
   hsize_t nsel;
   herr_t errf;

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   data.clear();
   data.resize(nsel);

   if (nsel > 0) {
      /*
        We cannot be sure exactly what bool is because it is
        represented differently on different platforms, and
        it may have been redefined, i.e., by the Boolean
        type.  So we read bools into native integer memory
        then convert.
      */
      vector<int> data1( nsel );
      int* locPtr = &data1[0];
      errf = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      // Convert what was just read in.
      for ( size_t i=0; i<nsel; ++i ) {
         data[i] = data1[i];
      }
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );

#endif
}

void HDFDatabase::getBoolArray(
   const string& key,
   bool* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   vector<bool> tmp;
   getBoolArray(key, tmp);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getBoolArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a char entry.  If the key does not exist, then false       *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isChar(const string& key)
{
   bool is_char  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_CHAR_ARRAY) {
            is_char = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return is_char;
}

/*
*************************************************************************
*                                                                       *
* Create a char scalar entry in the database with the specified         *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putChar(
   const string& key, 
   char data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   putCharArray(key, &data, 1);

}

/*
*************************************************************************
*                                                                       *
* Create a char array entry in the database with the specified          *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putCharArray(
   const string& key, 
   const vector<char>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if ( data.size() > 0 ) {
      putCharArray(key, &data[0], data.size());
   } else { 
      TBOX_ERROR("HDFDatabase::putCharArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a char array entry in the database with the specified          *
* key name. The charentry is defined by the hdf type H5T_C_S1.          *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putCharArray(
   const string& key,
   const char* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != (char*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hid_t atype, space, dataset;

      char* local_buf = new char[nelements];

      atype = H5Tcopy(H5T_C_S1);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( atype >= 0 );
#endif
      errf = H5Tset_size(atype, nelements);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Tset_strpad(atype, H5T_STR_NULLTERM);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      for (int i = 0; i < nelements; i++) {
         local_buf[i] = data[i];
      }

      space = H5Screate(H5S_SCALAR);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( space >= 0 );
#endif

      dataset = 
         H5Dcreate(d_group_id, key.c_str(), atype, space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( dataset >= 0 );
#endif

      errf = H5Dwrite(dataset, atype, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, local_buf);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_CHAR_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Tclose(atype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      delete[] local_buf;

   } else {
      TBOX_ERROR("HDFDatabase::putCharArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get char scalar entry from the database with the specified key       *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a char type.                                                    *
*                                                                      *
************************************************************************
*/

char HDFDatabase::getChar(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   char ret_val;
   getCharArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Two routines to get char arrays from the database with the           *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a char type.                      *
*                                                                      *
************************************************************************
*/

void HDFDatabase::getCharArray(const string& key,
                               vector<char>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if ( !isChar(key) ) {
      TBOX_ERROR("HDFDatabase::getCharArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a char array." << endl);
   } 

   hid_t   dset, dspace, dtype;
   size_t  nsel = 0;
   herr_t errf;

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dspace >= 0 );
#endif
   dtype  = H5Dget_type(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dtype >= 0 );
#endif
   nsel   = H5Tget_size(dtype);

   data.clear();
   data.resize(nsel);

   if (nsel > 0) {
      char* locPtr = &data[0];
      errf = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Tclose(dtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
}

void HDFDatabase::getCharArray(
   const string& key,
   char* data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   vector<char> tmp;
   getCharArray(key, tmp);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getCharArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a double entry.  If the key does not exist, then false     *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isDouble(const string& key)
{
   bool is_double = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_DOUBLE_ARRAY) {
            is_double = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return is_double;
}

/*
*************************************************************************
*                                                                       *
* Create a double scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDouble(
   const string& key, 
   double data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   putDoubleArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a double array entry in the database with the specified        *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDoubleArray(
   const string& key,
   const vector<double>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if ( data.size() > 0 ) {
      putDoubleArray(key, &data[0], data.size());
   } else {
      TBOX_ERROR("HDFDatabase::putDoubleArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a double array entry in the database with the specified        *
* key name.  The array type is based on the hdf type H5T_NATIVE_HDOUBLE.*
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putDoubleArray(
   const string& key,
   const double* const data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != (double*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( space >= 0 );

#endif
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_MPTCOUPLER_DOUBLE, 
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );

#endif
      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_DOUBLE_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );

#endif
   } else {
      TBOX_ERROR("HDFDatabase::putDoubleArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   } 
}

/*
************************************************************************
*                                                                      *
* Get double scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a double type.                                                  *
*                                                                      *
************************************************************************
*/

double HDFDatabase::getDouble(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   double ret_val;
   getDoubleArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Two routines to get double arrays from the database with the         *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a double type.                    *
*                                                                      *
************************************************************************
*/

void HDFDatabase::getDoubleArray(const string& key,
                                 vector<double>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   herr_t errf;
   if (!isDouble(key)) {
     TBOX_ERROR("HDFDatabase::getDoubleArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a double array." << endl);
   }

   hid_t   dset, dspace;
   hsize_t nsel;

   dset = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   data.clear();
   data.resize(nsel);

   if (nsel > 0) {
      double* locPtr = &data[0];
      errf = H5Dread(dset, H5T_NATIVE_DOUBLE, 
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
}

void HDFDatabase::getDoubleArray(
   const string& key, 
   double* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   vector<double> tmp;
   getDoubleArray(key, tmp);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getDoubleArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a float entry.  If the key does not exist, then false      *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isFloat(const string& key)
{
   bool is_float  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_FLOAT_ARRAY) {
            is_float = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return is_float;
}

/*
*************************************************************************
*                                                                       *
* Create a float scalar entry in the database with the specified        *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putFloat(
   const string& key, 
   float data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   putFloatArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a float array entry in the database with the specified         *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putFloatArray(
   const string& key,
   const vector<float>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if ( data.size() > 0 ) {
      putFloatArray(key, &data[0], data.size());
   } else {
      TBOX_ERROR("HDFDatabase::putFloatArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
  
}

/*
*************************************************************************
*                                                                       *
* Create a float array entry in the database with the specified         *
* key name.  The array type is based on the hdf type H5T_NATIVE_HFLOAT. *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putFloatArray(
   const string& key,
   const float* const data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != (float*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( space >= 0 );

#endif
      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_MPTCOUPLER_FLOAT, 
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( dataset >= 0 );
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_FLOAT_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif

   } else {
      TBOX_ERROR("HDFDatabase::putFloatArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get float scalar entry from the database with the specified key      *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a float type.                                                   *
*                                                                      *
************************************************************************
*/

float HDFDatabase::getFloat(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   float ret_val;
   getFloatArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Two routines to get float arrays from the database with the          *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a float type.                     *
*                                                                      *
************************************************************************
*/

void HDFDatabase::getFloatArray(const string& key,
                                vector<float>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   herr_t errf;
   if (!isFloat(key)) {
      TBOX_ERROR("HDFDatabase::getFloatArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a float array." << endl);
   }
 
   hid_t   dset, dspace;
   hsize_t nsel;

   dset = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   data.clear();
   data.resize(nsel);
 
   if (nsel > 0) {
      float* locPtr = &data[0];
      errf = H5Dread(dset, H5T_NATIVE_FLOAT, 
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
}

void HDFDatabase::getFloatArray(
   const string& key,
   float* data,
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   vector<float> tmp;
   getFloatArray(key, tmp);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getFloatArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a integer entry.  If the key does not exist, then false    *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isInteger(const string& key)
{
   bool is_int  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_INT_ARRAY) {
            is_int = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return is_int;
}

/*
*************************************************************************
*                                                                       *
* Create a integer scalar entry in the database with the specified      *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putInteger(
   const string& key,
   int data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   putIntegerArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create an integer array entry in the database with the specified      *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putIntegerArray(
   const string& key, 
   const vector<int>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if ( data.size() > 0 ) {
      putIntegerArray(key, &data[0], data.size());
   } else {
      TBOX_ERROR("HDFDatabase::putIntegerArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create an integer array entry in the database with the specified      *
* key name.  The array type is based on the hdf type H5T_NATIVE_HINT.   *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putIntegerArray(
   const string& key, 
   const int* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != (int*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert(space >= 0);
#endif

      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), H5T_MPTCOUPLER_INT, 
                                space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert(dataset >= 0);
#endif
      errf = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, data);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert(errf >= 0);
#endif

      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_INT_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif

   } else {
      TBOX_ERROR("HDFDatabase::putIntegerArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get integer scalar entry from the database with the specified key    *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a integer type.                                                 *
*                                                                      *
************************************************************************
*/

int HDFDatabase::getInteger(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   int ret_val;
   getIntegerArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Two routines to get integer arrays from the database with the        *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a integer type.                   *
*                                                                      *
************************************************************************
*/

void HDFDatabase::getIntegerArray(const string& key,
                                  vector<int>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   herr_t errf;
   if ( !isInteger(key) ) {
      TBOX_ERROR("HDFDatabase::getIntegerArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not an integer array." << endl);
   }

   hid_t   dset, dspace;
   hsize_t nsel;

   dset = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dspace >= 0 );
#endif
   nsel   = H5Sget_select_npoints(dspace);

   data.clear();
   data.resize(nsel); 

   if (nsel > 0) {
      int* locPtr = &data[0];
      errf = H5Dread(dset, H5T_NATIVE_INT, 
                     H5S_ALL, H5S_ALL, H5P_DEFAULT, locPtr);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );

#endif
}

void HDFDatabase::getIntegerArray(
   const string& key, 
   int* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   vector<int> tmp;
   getIntegerArray(key, tmp);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getIntegerArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }
}

/*
*************************************************************************
*                                                                       *
* Return true or false depending on whether the specified key           *
* represents a string entry.  If the key does not exist, then false     *
* is returned.                                                          *
*                                                                       *
*************************************************************************
*/

bool HDFDatabase::isString(const string& key)
{
   bool is_string  = false;
   herr_t errf;

   if (!key.empty()) {
      hid_t this_set;
      BEGIN_SUPPRESS_HDF5_WARNINGS
      this_set = H5Dopen(d_group_id, key.c_str());
      END_SUPPRESS_HDF5_WARNINGS
      if (this_set > 0) {
         int type_key = readAttribute(this_set);
         if (type_key == KEY_STRING_ARRAY) {
            is_string = true;
         }
         errf = H5Dclose(this_set);
#ifdef ASSERT_HDF5_RETURN_VALUES
         assert( errf >= 0 );
#endif
      }
   }

   return is_string;
}

/*
*************************************************************************
*                                                                       *
* Create a string scalar entry in the database with the specified       *
* key name.  A scalar entry is an array of one.                         *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putString(
   const string& key, 
   const string& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   putStringArray(key, &data, 1);
}

/*
*************************************************************************
*                                                                       *
* Create a string array entry in the database with the specified        *
* key name.                                                             *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putStringArray(
   const string& key, 
   const vector<string>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   if ( data.size() > 0 ) {
      putStringArray(key, &data[0], data.size());
   } else {
      TBOX_ERROR("HDFDatabase::putStringArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
*************************************************************************
*                                                                       *
* Create a double array entry in the database with the specified        *
* key name.  The array type is based on the hdf type H5T_C_S1.          *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::putStringArray(
   const string& key, 
   const string* const data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
   assert(data != (string*)NULL);
#endif
   herr_t errf;
   if (nelements > 0) {

      int maxlen = 0;
      int current, data_size;
      int i;
      for (i = 0; i < nelements; i++) {
         current = data[i].size(); 
         if ( current > maxlen ) maxlen = current;
      }

      char* local_buf = new char[nelements*(maxlen+1)];
      for (i = 0; i < nelements; i++) {
         strcpy(&local_buf[i*(maxlen+1)], data[i].c_str());
         data_size = data[i].size();
         if (data_size < maxlen) {
            memset(&local_buf[i*(maxlen+1)] + data_size + 1, 0, 
                   maxlen - data_size);
         }
      }

      hid_t atype = H5Tcopy(H5T_C_S1);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( atype >= 0 );
#endif
      errf = H5Tset_size(atype, maxlen+1);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Tset_strpad(atype, H5T_STR_NULLTERM);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif

      hsize_t dim[] = {nelements};
      hid_t space = H5Screate_simple(1, dim, NULL);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( space >= 0 );
#endif

      hid_t dataset = H5Dcreate(d_group_id, key.c_str(), 
                                atype, space, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( dataset >= 0 );
#endif

      errf = H5Dwrite(dataset, atype, H5S_ALL, H5S_ALL, 
                      H5P_DEFAULT, local_buf);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      // Write attribute so we know what kind of data this is.
      writeAttribute( KEY_STRING_ARRAY, dataset );

      errf = H5Sclose(space);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Tclose(atype);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      errf = H5Dclose(dataset);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      delete[] local_buf;

   } else {
      TBOX_ERROR("HDFDatabase::putStringArray() error in database "
         << d_database_name
         << "\n    Attempt to put zero-length array with key = "
         << key << endl);
   }
}

/*
************************************************************************
*                                                                      *
* Get string scalar entry from the database with the specified key     *
* name. An error message is printed and the program exits if the       *
* specified key does not exist in the database or is not associated    *
* with a string type.                                                  *
*                                                                      *
************************************************************************
*/

string HDFDatabase::getString(const string& key)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   string ret_val;
   getStringArray(key, &ret_val, 1);

   return(ret_val);
}

/*
************************************************************************
*                                                                      *
* Two routines to get string arrays from the database with the         *
* specified key name. In any case, an error message is printed and     *
* the program exits if the specified key does not exist in the         *
* database or is not associated with a string type.                    *
*                                                                      *
************************************************************************
*/

void HDFDatabase::getStringArray(const string& key,
                                 vector<string>& data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   herr_t errf;
   if (!isString(key)) {
      TBOX_ERROR("HDFDatabase::getStringArray() error in database "
         << d_database_name
         << "\n    Key = " << key << " is not a string array." << endl);
   }

   hsize_t nsel;
   size_t  dsize;
   hid_t   dset, dspace, dtype;
   char*   local_buf;

   dset   = H5Dopen(d_group_id, key.c_str());
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dset >= 0 );
#endif
   dspace = H5Dget_space(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dspace >= 0 );
#endif
   dtype  = H5Dget_type(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( dtype >= 0 );
#endif
   dsize  = H5Tget_size(dtype);
   nsel   = H5Sget_select_npoints(dspace);

   local_buf = new char[nsel*dsize];

   errf = H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, local_buf);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif

   data.clear();
   data.resize(nsel);

   for (int i = 0; i < (int)nsel; i++) {
      data[i] = &local_buf[i*dsize];
   }

   errf = H5Sclose(dspace);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Tclose(dtype);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Dclose(dset);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif

   delete[] local_buf;
}

void HDFDatabase::getStringArray(
   const string& key, 
   string* data, 
   const int nelements)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!key.empty());
#endif
   vector<string> tmp;
   getStringArray(key, tmp);
   const int tsize = tmp.size();

   if (nelements != tsize) {
      TBOX_ERROR("HDFDatabase::getStringArray() error in database "
         << d_database_name
         << "\n    Incorrect array size = " << nelements
         << " given for key = " << key
         << "\n    Actual array size = " << tsize << endl);
   }

   for (int i = 0; i < tsize; i++) {
      data[i] = tmp[i];
   }

}


void HDFDatabase::writeAttribute( int type_key,
                                       hid_t dataset_id
                                       )
{
   herr_t errf;
   hid_t attr_id = H5Screate(H5S_SCALAR);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( attr_id >= 0 );
#endif
   hid_t attr = H5Acreate(dataset_id, "Type", H5T_MPTCOUPLER_ATTR, 
                          attr_id, H5P_DEFAULT);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( attr >= 0 );
#endif
   errf = H5Awrite(attr, H5T_NATIVE_INT, &type_key);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Aclose(attr);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Sclose(attr_id);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
}


int HDFDatabase::readAttribute( hid_t dataset_id )
{
   herr_t errf;
   hid_t attr = H5Aopen_name(dataset_id, "Type");
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( attr >= 0 );
#endif
   int type_key;
   errf = H5Aread(attr, H5T_NATIVE_INT, &type_key);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Aclose(attr);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   return type_key;
}


/*
*************************************************************************
*                                                                       *
* Print contents of current database to the specified output stream.    *
* Note that contents of subdatabases will not be printed.  This must    *
* be done by iterating through all the subdatabases individually.       * 
*                                                                       *
*************************************************************************
*/

void HDFDatabase::printClassData(ostream& os)
{

   performKeySearch();

   if (d_keydata.empty()) {
      os << "Database named `"<< d_database_name 
         << "' has zero keys..." << endl;
   } else {
      os << "Printing contents of database named `" 
         << d_database_name << "'..." << endl;
   }

   for (list<KeyData>::iterator ik = d_keydata.begin();
        ik != d_keydata.end(); ++ik) {
      int t = (*ik).d_type;
      switch (MathUtilities<int>::Abs(t)) {
         case KEY_DATABASE: {
            os << "   Data entry `"<< (*ik).d_key << "' is"
               << " a database" << endl;   
            break;
         }
         case KEY_BOOL_ARRAY: {
            os << "   Data entry `"<< (*ik).d_key << "' is" << " a boolean ";
            os << ( (t < 0) ? "scalar" : "array") << endl;   
            break;
         }
         case KEY_CHAR_ARRAY: {
            os << "   Data entry `"<< (*ik).d_key << "' is" << " a char ";
            os << ( (t < 0) ? "scalar" : "array") << endl;   
            break;
         }
         case KEY_DOUBLE_ARRAY: {
            os << "   Data entry `"<< (*ik).d_key << "' is" << " a double ";
            os << ( (t < 0) ? "scalar" : "array") << endl;   
            break;
         }
         case KEY_FLOAT_ARRAY: {
            os << "   Data entry `"<< (*ik).d_key << "' is" << " a float ";
            os << ( (t < 0) ? "scalar" : "array") << endl;   
            break;
         }
         case KEY_INT_ARRAY: {
            os << "   Data entry `"<< (*ik).d_key << "' is" << " an integer ";
            os << ( (t < 0) ? "scalar" : "array") << endl;   
            break;
         }
         case KEY_STRING_ARRAY: {
            os << "   Data entry `"<< (*ik).d_key << "' is" << " a string ";
            os << ( (t < 0) ? "scalar" : "array") << endl;   
            break;
         }
         default: {
            TBOX_ERROR("HDFDatabase::printClassData error....\n"
               << "   Unable to identify key = " << (*ik).d_key
               << " as a known group or dataset" << endl);
         }
      }
   }

   cleanupKeySearch();

}

/*
*************************************************************************
*                                                                       *
* Open a HDF5 data-base file.  Flags can be "R" for read only or "W"    *
* for write.  If the file does not exist when read only is specified,   *
* an error status is returned (<0).  If the file exists when opened     *
* for write, an error status is returned (<0).                          * 
*                                                                       *
*************************************************************************
*/ 

int HDFDatabase::mount(
   const string& file_name, 
   const string& flags)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!file_name.empty());
   assert(!flags.empty());
#endif

   int status = 1;

   hid_t file_id = 0;

   if (flags == "R") {
      file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if ( file_id < 0 ) {
         TBOX_ERROR("HDFDatabase::mount error..."
                    << "\n   Unable to open HDF5 file " 
                    << file_name << endl);
      }
   } else if (flags == "W") {
      file_id = H5Fopen(file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      if ( file_id < 0 ) {
         TBOX_ERROR("HDFDatabase::mount error..."
                    << "\n   Unable to open HDF5 file " 
                    << file_name << endl);
      }
   } else if (flags == "WN") { 
      file_id = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, 
                          H5P_DEFAULT, H5P_DEFAULT);
      if ( file_id < 0 ) {
         TBOX_ERROR("HDFDatabase::mount error..."
                    << "\n   Unable to open or create HDF5 file " 
                    << file_name << endl);
      }
   } else {
      TBOX_ERROR("HDFDatabase::mount error...\n"
                 << "   database name is " << d_database_name
                 << "\n    file name = " << file_name
                 << "\n    unrecognized flag = " << flags << endl);  
   }

   if (file_id < 0) {
      status = (int)file_id;
   } else {
      d_is_file  = true;
      d_group_id = file_id;
      d_file_id  = file_id;
   }

   return status;

}

/*
*************************************************************************
*                                                                       *
* Close the open HDF data file specified by d_file_id.                  *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::unmount()
{
   herr_t errf;
   if (d_is_file) {
      errf = H5Fclose(d_file_id);
#ifdef ASSERT_HDF5_RETURN_VALUES
      assert( errf >= 0 );
#endif
      if ( d_group_id == d_file_id ) d_group_id = -1;
      d_file_id = -1;
      d_is_file = false;
   }
}

/*
*************************************************************************
*                                                                       *
* Private helper function for writing arrays in HDF5.  This function    *
* was deprecated in HDF5 1.4.  We replicate it here since it makes      *
* arrays easier to use in this database class.                          *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::insertArray(
   hid_t parent_id, 
   const char *name, 
   hsize_t offset, 
   int ndims, 
   const hsize_t dim[/*ndims*/], 
   const int *perm, 
   hid_t member_id) const
{
   herr_t errf;
#if (H5_VERS_MAJOR > 1) || ((H5_VERS_MAJOR == 1) && (H5_VERS_MINOR > 2))

   hid_t array = H5Tarray_create(member_id, ndims, dim, perm);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( array >= 0 );
#endif
   errf = H5Tinsert(parent_id, name, offset, array);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
   errf = H5Tclose(array);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
#else
   size_t newdim[H5S_MAX_RANK];
   for(int i = 0; i < ndims; i++) {
     newdim[i] = dim[i];
   }
    
   errf = H5Tinsert_array(parent_id, name, offset, ndims, newdim, perm, member_id);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
#endif
}

/*
*************************************************************************
*                                                                       *
* Private helper function for searching database keys.                  *
*                                                                       *
*************************************************************************
*/

void HDFDatabase::performKeySearch()
{
   herr_t errf;
   if (d_is_file) {
      s_group_to_search = "/";
      s_top_level_search_group = "/";
      s_found_group = 1;
   } else {
      s_group_to_search = d_database_name;
      s_top_level_search_group = string();
      s_found_group = 0;
   }

   s_still_searching = 1;

   errf = H5Giterate(d_group_id, "/", NULL,
                     HDFDatabase::iterateKeys, (void*)this);
#ifdef ASSERT_HDF5_RETURN_VALUES
   assert( errf >= 0 );
#endif
}

void HDFDatabase::cleanupKeySearch()
{
   s_top_level_search_group = string();
   s_group_to_search = string();
   s_still_searching = 0;
   s_found_group = 0;

   d_keydata.clear();
}

}

#endif



