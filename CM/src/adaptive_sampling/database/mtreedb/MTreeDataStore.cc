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
// File:        MTreeDataStore.cc
// Package:     MTree database
// 
// 
// 
// Description: Manager class for MTree node allocation and storage of data objects.
//

#ifndef included_MtreeDataStore_C
#define included_MTreeDataStore_C

#ifndef included_MTreeDataStore
#include "MTreeDataStore.h"
#endif

#ifdef HAVE_MPI
#ifndef included_toolbox_MPI
#include "toolbox/parallel/MPI.h"
#endif 
#endif // HAVE_MPI

#ifndef included_toolbox_Utilities
#include "toolbox/base/Utilities.h"
#endif 

//#ifdef DEBUG_CHECK_ASSERTIONS
//#ifndef included_cassert
//#define included_cassert
//#include <cassert>
//#endif
//#endif
#include <assert.h>

#ifndef included_MTree
#include "MTree.h"
#endif 

#ifndef included_MTreeEntry
#include "MTreeEntry.h"
#endif 

#ifdef DEBUG_NO_INLINE
#include "MTreeDataStore.I"
#endif

#if 0
#include <strstream>  // deprecated
#else
#include <sstream>     
#endif

      //
      // 
      //

      namespace {

	//
	// generate string representation of MPI rank
	//

	inline 
	std::string 
	getStringMPIRankRepresentation(int width = 6)
	{

	  //	  std::strstream representationStream;
	  std::stringstream representationStream;
	  representationStream << std::setw(width) << std::setfill('0')
#ifdef HAVE_MPI
			       << toolbox::MPI::getRank();
#else
	                       << 0;
#endif // HAVE_MPI
	  std::string representation;
	  representationStream >> representation;

	  return representation;

	}

      }

/*
*************************************************************************
*                                                                       *
* Default ctor and dtor.                                                *
*                                                                       *
*************************************************************************
*/

MTreeDataStore::MTreeDataStore()
: d_mtree((MTree*)NULL),
  d_object_factory((const DBObjectFactory*)NULL),
  d_is_initialized(false),
  d_is_open(false),
  d_num_leaf_nodes(0),
  d_num_objects(0),
  d_object_file_capacity(MTREE_DATA_STORE_OBJECT_FILE_CAPACITY),
  d_num_objects_in_files(0)
{
}

MTreeDataStore::~MTreeDataStore()
{

   d_mtree = (MTree*)NULL;
   d_object_factory = (const DBObjectFactory*)NULL;

   const int num_files = d_object_file_info.size();
   for (int file_index = 0; file_index < num_files; ++file_index) {
      if ( d_object_file_info[file_index] ) {
         delete d_object_file_info[file_index];
      }
   }
   d_object_file_info.clear();

   const int num_objects = d_object_info.size();
   for (int obj_index = 0; obj_index < num_objects; ++obj_index) {
      if ( d_object_info[obj_index] ) {
         delete d_object_info[obj_index];
      }
   }
   d_object_info.clear();

   const int num_leaf_nodes = d_leaf_node_info.size();
   for (int leaf_index = 0; leaf_index < num_leaf_nodes; ++leaf_index) {
      if ( d_leaf_node_info[leaf_index] ) {
         delete d_leaf_node_info[leaf_index];
      }
   }
   d_leaf_node_info.clear();

}

/*
*************************************************************************
*                                                                       *
* Initialize new data store.                                            * 
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::create(MTree* mtree,
                            const DBObjectFactory* obj_factory,
                            const string& directory_name,
                            const string& file_prefix)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(mtree != (MTree*)NULL);
   assert(obj_factory != (DBObjectFactory*)NULL);
   assert(!directory_name.empty());
   assert(!file_prefix.empty());
#endif

   if (!d_is_initialized) {

      d_mtree = mtree;
      d_object_factory = obj_factory;

      d_directory_name = directory_name + "_" + 
	getStringMPIRankRepresentation();

      d_file_prefix = file_prefix;

      d_mtree_index_file_name = getFullPathMTreeIndexFileName(); 

      d_data_store_file_name = getFullPathDataStoreFileName();

      d_object_file_prefix = getObjectFilePrefix();

      d_num_leaf_nodes = 0;
      d_leaf_node_info.reserve(LEAF_NODE_VECTOR_ALLOCATION_CHUNK);
      d_recycled_leaf_node_indices.clear();

      d_object_info.reserve(OBJECT_VECTOR_ALLOCATION_CHUNK);

#ifdef HAVE_PKG_hdf5
      toolbox::Utilities::recursiveMkdir(d_directory_name, 
					 S_IRWXU|S_IRWXG,
					 false);
#endif

      d_is_initialized = true;
      d_is_open        = true;

   } else {
      TBOX_WARNING("MTreeDataStore::create() warning"
                   << "\nCannot call create() method when "
                   << " store is already initialized." 
                   << "\n create() method call did nothing!" << endl);
   }

}

/*
*************************************************************************
*                                                                       *
* Open existing data store.                                             * 
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::open(MTree* mtree,
                          const DBObjectFactory* obj_factory,
                          const string& directory_name,
                          const string& file_prefix)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(mtree != (MTree*)NULL);
   assert(obj_factory != (DBObjectFactory*)NULL);
   assert(!directory_name.empty());
   assert(!file_prefix.empty());
#endif


   if (!d_is_initialized) {

      d_is_open = true;

      d_mtree = mtree;
      d_object_factory = obj_factory;

#ifdef HAVE_PKG_hdf5
      d_directory_name = directory_name + "_" + 
	getStringMPIRankRepresentation();

      d_file_prefix = file_prefix;

      d_mtree_index_file_name = getFullPathMTreeIndexFileName();

      d_data_store_file_name = getFullPathDataStoreFileName(); 

      d_object_file_prefix = getObjectFilePrefix();

      toolbox::HDFDatabasePtr 
         read_store_DB( new toolbox::HDFDatabase(d_data_store_file_name) );
      bool mount_successful =
         (read_store_DB->mount(d_data_store_file_name, "R") >= 0);

      bool read_successful = mount_successful;

      if ( read_successful) {
         read_successful = readObjectFileInfo(read_store_DB);
      }

      if ( read_successful) {
         read_successful = readObjectInfo(read_store_DB);
      }

      if ( read_successful) {
         bool rebuild_from_open = true;
         read_successful = readAllDataObjects(rebuild_from_open);
      }

      if ( mount_successful ) {
         read_store_DB->unmount();
      }
      read_store_DB.reset();

      if ( read_successful ) {
 
         toolbox::HDFDatabasePtr
            read_tree_DB( new toolbox::HDFDatabase(d_mtree_index_file_name) );
         mount_successful =
            (read_tree_DB->mount(d_mtree_index_file_name, "R") >= 0);
 
         read_successful = mount_successful;
 
         if (read_successful) {
            read_successful = d_mtree->getFromDatabase(read_tree_DB);
         }
 
         if ( mount_successful ) {
            read_tree_DB->unmount();
         }
         read_tree_DB.reset();
 
      }

      if ( !read_successful ) {
         TBOX_ERROR("MTreeDataStore::open() error..."
                    << "\nSome operation failed when reading data" 
                    << " from files. \nState of data store may not "
                    << " be initialized properly. " << endl);
      } else {
         d_is_initialized = true;
      }

#endif

   } else {
      TBOX_WARNING("MTreeDataStore::open() warning"
                   << "\nCannot call open() method when "
                   << " store is already initialized."
                   << "\n open() method call did nothing!" << endl);
   }

}

/*
*************************************************************************
*                                                                       *
* Close data store and write all unwritten data to files.               * 
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::close()
{
   if (d_is_open) {

#ifdef HAVE_PKG_hdf5
      bool write_successful = writeAllDataObjects();

      toolbox::HDFDatabasePtr 
         write_store_DB( new toolbox::HDFDatabase(d_data_store_file_name) );
      bool mount_successful =
         (write_store_DB->mount(d_data_store_file_name, "WN") >= 0);

      write_successful &= mount_successful;
      
      if (write_successful) {
         write_successful = writeObjectFileInfo(write_store_DB);
      }
  
      if (write_successful) {
         write_successful = writeObjectInfo(write_store_DB);
      }

      if ( mount_successful ) {
         write_store_DB->unmount();
      }
      write_store_DB.reset();

      if ( write_successful ) {

         toolbox::HDFDatabasePtr 
            write_tree_DB( new toolbox::HDFDatabase(d_mtree_index_file_name) );
         mount_successful =
            (write_tree_DB->mount(d_mtree_index_file_name, "WN") >= 0);

         write_successful = mount_successful;
      
         if (write_successful) {
            write_successful = d_mtree->putToDatabase(write_tree_DB);
         }

         if ( mount_successful ) {
            write_tree_DB->unmount();
         }
         write_tree_DB.reset();

      }

      if ( !write_successful ) {
         TBOX_ERROR("MTreeDataStore::close() error..."
                    << "\nSome operation failed when writing data" 
                    << " to files. \nState of data store may not "
                    << " be accurately preserved in files. " << endl);
      }

#endif

   }

   d_is_initialized = false;
   d_is_open        = false;
}

/*
*************************************************************************
*                                                                       *
* Add leaf node to store and set its leaf node identifier.              *
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::addLeafNode(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif

   if ( d_is_open ) { 

   if ( node->getLeafNodeId() == MTreeNode::getUndefinedId() ) {

      if (d_recycled_leaf_node_indices.empty()) {

         if (d_leaf_node_info.size() == d_leaf_node_info.capacity() ) {
            d_leaf_node_info.reserve( d_leaf_node_info.size() + 
                                      LEAF_NODE_VECTOR_ALLOCATION_CHUNK );
         }

         node->setLeafNodeId( d_leaf_node_info.size() );

         d_leaf_node_info.push_back( new LeafNodeInfo(node) );

      } else {

         node->setLeafNodeId( d_recycled_leaf_node_indices.front() );
         d_recycled_leaf_node_indices.pop_front();

         d_leaf_node_info[ node->getLeafNodeId() ] = new LeafNodeInfo(node); 

      }

      d_num_leaf_nodes++;

   } else {

      const int leaf_node_id = node->getLeafNodeId();
      const int leaf_node_info_size = d_leaf_node_info.size(); 

      if ( (leaf_node_id >= 0) && 
           (leaf_node_id < leaf_node_info_size) ) {

         if ( !d_leaf_node_info[leaf_node_id] ) {
            TBOX_ERROR("MTreeDataStore::addLeafNode() error..."
                       << "\nGiven node with id = " << node->getNodeId()
                       << " and leaf node id = " << leaf_node_id 
                       << " has leaf node id set but is not in data store." 
                       << "\nLeaf node id has probably been set outside of "
                       << " data store which will have unexpected behavior."
                       << endl);
         }
         if ( d_leaf_node_info[leaf_node_id]->getNode().get() != node.get() ) {
            TBOX_ERROR("MTreeDataStore::addLeafNode() error..."
                       << "\nGiven node with id = " << node->getNodeId()
                       << " and leaf node id = " << leaf_node_id
                       << " does not match leaf node in data store with" 
                       << " same id."
                       << "\nLeaf Node id has probably been set outside of "
                       << " data store which will have unexpected behavior."
                       << endl);
         }

      } else {

         TBOX_ERROR("MTreeDataStore::addLeafNode() error..."
                    << "\nGiven node with id = " << node->getNodeId()
                    << " has an unknown leaf node id = " << leaf_node_id
                    << "\nLeaf Node id has probably been set outside of "
                    << " data store which will have unexpected behavior."
                    << endl);

      }

   }

   } // if data store is open

}

/*
*************************************************************************
*                                                                       *
* Remove leaf node from data store.                                     *
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::removeLeafNode(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif

   if ( d_is_open ) { 

   const int leaf_node_id = node->getLeafNodeId();
   const int leaf_node_info_size = d_leaf_node_info.size(); 

   if ( (leaf_node_id >= 0) && 
        (leaf_node_id < leaf_node_info_size) &&
        d_leaf_node_info[leaf_node_id] ) {

       delete d_leaf_node_info[leaf_node_id];
       d_leaf_node_info[leaf_node_id] = (LeafNodeInfo*)NULL;
       d_recycled_leaf_node_indices.push_front(leaf_node_id);
       d_num_leaf_nodes--;

   }

   } // if data store is open

}

/*
*************************************************************************
*                                                                       *
* Determine whether given leaf node id is associated with a             *
* valid leaf node in data store.                                        * 
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::isValidLeafNodeId(int leaf_node_id) const
{
   bool ret_val = false;

   if ( d_is_open ) { 

   const int leaf_node_info_size = d_leaf_node_info.size(); 

   if ( (leaf_node_id >= 0) && 
        (leaf_node_id < leaf_node_info_size) &&
        d_leaf_node_info[leaf_node_id] ) {

      ret_val = true;

   }

   } // if data store is open

   return (ret_val); 
}

/*
*************************************************************************
*                                                                       *
* Add copy of data object to store and set identifier of orginal and    *
* copy if no potential problems.                                        *
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::addObject(DBObject& object)
{

   if ( d_is_open ) { 

   if (d_recycled_object_indices.empty()) {

      if (d_object_info.size() == d_object_info.capacity() ) {
         d_object_info.reserve( d_object_info.size() +
                                OBJECT_VECTOR_ALLOCATION_CHUNK );
      }

      DBObjectPtr store_obj( d_object_factory->cloneObject(object) );

      store_obj->setObjectId( d_object_info.size() );
      object.setObjectId( store_obj->getObjectId() );

      d_object_info.push_back( new ObjectInfo(store_obj) );

   } else {

      DBObjectPtr store_obj( d_object_factory->cloneObject(object) );

      store_obj->setObjectId( d_recycled_object_indices.front() );
      d_recycled_object_indices.pop_front();

      object.setObjectId( store_obj->getObjectId() );

      d_object_info[ store_obj->getObjectId() ] = new ObjectInfo(store_obj);

   }

   d_object_info[ object.getObjectId() ]->setInMemory(true);

   d_num_objects++;

   } // if data store is open

}

/*
*************************************************************************
*                                                                       *
* Remove object with given id from store.                               *
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::removeObject(int object_id)
{

   if ( d_is_open ) { 

   const int object_info_size = d_object_info.size();
   if ( (object_id >= 0) && 
        (object_id < object_info_size) &&
        d_object_info[object_id] &&
        !d_object_info[object_id]->getInFile() ) {

      delete d_object_info[object_id];
      d_object_info[object_id] = (ObjectInfo*)NULL;
      d_recycled_object_indices.push_front(object_id);
      d_num_objects--;

   }

   } // if data store is open

}

/*
*************************************************************************
*                                                                       *
* Determine whether object with given id lives in data store.           *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::isValidObjectId(int object_id) const
{
   bool ret_val = false;

   if ( d_is_open ) { 

   const int object_info_size = d_object_info.size();
   if ( (object_id >= 0) && 
        (object_id < object_info_size) &&
        d_object_info[object_id] ) {

      ret_val = true;

   }

   } // if data store is open

   return (ret_val); 
}

/*
*************************************************************************
*                                                                       *
* Return deep copy of data object in store with given identifier.       *
*                                                                       *
*************************************************************************
*/

DBObjectPtr MTreeDataStore::getObjectCopy(int object_id)
{
   DBObjectPtr ret_object;

   if ( d_is_open ) { 

   const int object_info_size = d_object_info.size(); 
   if ( (object_id >= 0) &&
        (object_id < object_info_size) &&
        d_object_info[object_id] ) {
     
     if ( d_object_info[object_id]->getInMemory() == false )
       readDataObject(object_id);
     
     // object has to be in memory now
     
     assert(d_object_info[object_id]->getInMemory() == true);

     if( d_object_info[object_id]->getInMemory() ) {
       
       ret_object = d_object_info[object_id]->getObject()->makeCopy();
       ret_object->setObjectId(object_id);
       
     }
   }
   
   } // if data store is open
   
   return(ret_object);
}

/*
*************************************************************************
*                                                                       *
* Return pointer to data object in store with given identifier.         *
*                                                                       *
*************************************************************************
*/

DBObjectPtr MTreeDataStore::getObjectPtr(int object_id)
{
   DBObjectPtr ret_object;

   if ( d_is_open ) { 

   const int object_info_size = d_object_info.size(); 
   if ( (object_id >= 0) &&
        (object_id < object_info_size) &&
        d_object_info[object_id] ){

     if ( d_object_info[object_id]->getInMemory() == false)
       readDataObject(object_id);

     // object has to be in memory now
     
     assert(d_object_info[object_id]->getInMemory() == true);
     
     if ( d_object_info[object_id]->getInMemory() ) {

       ret_object = d_object_info[object_id]->getObject();

     }
   }

   } // if data store is open

   return(ret_object);
}

/*
*************************************************************************
*                                                                       *
* Get node owning object with given identifier.                         *
*                                                                       *
*************************************************************************
*/

MTreeNodePtr MTreeDataStore::getNodeOwningObject(int object_id) const
{
   MTreeNodePtr ret_node;

   if ( d_is_open ) { 

   const int object_info_size = d_object_info.size();
   if ( (object_id >= 0) &&
        (object_id < object_info_size) &&
        d_object_info[object_id] ) {

      ret_node = 
         d_leaf_node_info[ 
            d_object_info[object_id]->getOwnerLeafNodeId() ]->getNode();

   }

   } // if data store is open

   return(ret_node);
}

/*
*************************************************************************
*                                                                       *
* Set information mapping data object to node.                          *
*                                                                       *
*************************************************************************
*/


bool MTreeDataStore::setNodeOwningObject(MTreeNodePtr node,
                                         MTreeEntryPtr entry)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
   assert(entry.get());
#endif
   bool operation_successful = false;

   if ( d_is_open ) { 

   const int object_info_size = d_object_info.size(); 
   const int leaf_node_info_size = d_leaf_node_info.size(); 

   const int obj_id = entry->getDataObjectId();
   const int leaf_node_id = node->getLeafNodeId();

   if ( node->isLeaf() &&
        entry->isDataEntry() &&
        (obj_id < object_info_size) &&
        d_object_info[obj_id] &&
        (leaf_node_id < leaf_node_info_size) && 
        d_leaf_node_info[leaf_node_id] ) {

      d_object_info[ obj_id ]->
         setOwnerLeafNodeInfo( leaf_node_id, 
                               entry->getPositionInNode() );

      d_leaf_node_info[ leaf_node_id ]->
         setOwnObjectId( obj_id, entry->getPositionInNode() );

      operation_successful = true;

   }

   } // if data store is open

   return(operation_successful);
}

/*
*************************************************************************
*                                                                       *
* Clean information mapping data objects to node.                       *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::clearOwnObjectIds(MTreeNodePtr node) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif
   bool operation_successful = false;

   if ( d_is_open ) { 

   const int leaf_node_info_size = d_leaf_node_info.size(); 
   const int leaf_node_id = node->getLeafNodeId();

   if ( node->isLeaf() &&
        (leaf_node_id < leaf_node_info_size) && 
        d_leaf_node_info[leaf_node_id] ) {

      vector<int>& object_ids = 
         d_leaf_node_info[ leaf_node_id ]->getObjectIds();

      for (unsigned int i = 0; i < object_ids.size(); ++i) {
         if ( ( object_ids[i] >= 0 ) &&
              d_object_info[ object_ids[i] ] ) {
            d_object_info[ object_ids[i] ]->
               setOwnerLeafNodeInfo( MTreeNode::getUndefinedId(), -1 );
         }
      }

      d_leaf_node_info[ leaf_node_id ]->clearOwnObjectIds();

      operation_successful = true;

   }

   } // if data store is open

   return(operation_successful);
}

/*
*************************************************************************
*                                                                       *
* Write all data object to HDF files.                                   *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::writeAllDataObjects()
{
   bool write_successful = false;

#ifdef HAVE_PKG_hdf5
   if ( d_is_open ) {

      write_successful = true;

      int last_file_index = -1;
      toolbox::HDFDatabasePtr write_DB;

      unsigned int object_id = 0;
      int object_count = 0;
      while ( (object_id < d_object_info.size()) &&
              (object_count < d_num_objects) ) {

         if ( isValidObjectId(object_id) ) {

            if ( !(d_object_info[object_id]->getInFile()) ) {

               string file_name;
               int file_index;
               bool new_file =
                  getFileNameAndIndexForObjectWrite(file_name,
                                                    file_index);
               string write_DB_name(d_object_file_prefix + file_name);

               if ( file_index != last_file_index ) {
 
                  if (write_DB.get()) {
                     write_DB->unmount();
                     write_DB.reset();
                  }
 
                  write_DB =
                     toolbox::HDFDatabasePtr(
                        new toolbox::HDFDatabase(write_DB_name) );
 
                  string flag = ( new_file ? "WN" : "W" );
                  write_successful &=
                     (write_DB->mount(write_DB_name, flag) >= 0);
 
               }  // if need to mount new file

               write_successful &=
                  writeDataObjectHelper(object_id,
                                        write_DB,
                                        file_index);
     
            }  // if object not already in file

            object_count++;
         } // if valid object id

         object_id++;
      }  // iterate over object ids

      if (write_DB.get()) {
         write_DB->unmount();
         write_DB.reset();
      }

      if (write_successful) {
         unsigned int leaf_node_id = 0;
         int node_count = 0;
         while ( (leaf_node_id < d_leaf_node_info.size()) &&
                 (node_count < d_num_leaf_nodes) ) {
            if ( isValidLeafNodeId(leaf_node_id) ) {
               d_leaf_node_info[leaf_node_id]->setInFile(true);
               d_leaf_node_info[leaf_node_id]->setInMemory(false);
               node_count++;
            }
            leaf_node_id++;
         }
      }

   }
#else
   TBOX_ERROR("MTreeDataStore error..."
              << "\n   Cannot write data to file..."
              << " code not compiled with HDF" << endl);
#endif

   return(write_successful);
}

/*
*************************************************************************
*                                                                       *
* Read all data objects from HDF files.                                 *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::readAllDataObjects(
   bool rebuild_from_open)
{
   bool read_successful = false;

#ifdef HAVE_PKG_hdf5
   if ( d_is_open ) {   

   read_successful = true;

   const int num_files = d_object_file_info.size();
   for (int file_index = 0; file_index < num_files; ++file_index) {

      const string read_DB_name( d_object_file_prefix + 
         d_object_file_info[ file_index ]->getFileName() );
 
      toolbox::HDFDatabasePtr read_DB( new toolbox::HDFDatabase(read_DB_name) );

      bool mount_successful = (read_DB->mount(read_DB_name, "R") >= 0); 
      read_successful &= mount_successful;

      const vector<int>& object_ids =
         d_object_file_info[ file_index ]->getObjectIds();

      const int num_objects = object_ids.size();

      if (rebuild_from_open) {
         for (int iobject = 0; iobject < num_objects; ++iobject) {
            const int object_id = object_ids[iobject];
            if ( d_object_info[ object_id ]->getInMemory() ) {
               read_successful &= readDataObjectHelper(object_id,
                                                       read_DB,
                                                       true);
            } 
         }
      } else {
         for (int iobject = 0; iobject < num_objects; ++iobject) {
            const int object_id = object_ids[iobject];
            if ( !(d_object_info[ object_id ]->getInMemory()) ) {
               read_successful &= readDataObjectHelper(object_id,
                                                       read_DB);
            }
         }
      }

      if ( mount_successful ) {
         read_DB->unmount();
      }
      read_DB.reset();
 
   }

   if ( !rebuild_from_open && read_successful ) {
      unsigned int leaf_node_id = 0;
      int node_count = 0;
      while ( (leaf_node_id < d_leaf_node_info.size()) &&
              (node_count < d_num_leaf_nodes) ) {
         if ( isValidLeafNodeId(leaf_node_id) ) {
            d_leaf_node_info[leaf_node_id]->setInMemory(true);
            node_count++;
         }
         leaf_node_id++;
      }
   }

   } // if data store is open
#else
   TBOX_ERROR("MTreeDataStore error..."
              << "\n   Cannot read data from files..."
              << " code not compiled with HDF" << endl);
#endif

   return(read_successful);
}

/*
*************************************************************************
*                                                                       *
* Write all data objects in leaf node to file                           *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::writeNodeDataObjects(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif

   bool write_successful = false;

#ifdef HAVE_PKG_hdf5
   if ( d_is_open ) {   

   const int leaf_node_id = node->getLeafNodeId(); 
   
   if ( isValidLeafNodeId(leaf_node_id) ) {

      if ( !(d_leaf_node_info[leaf_node_id]->getInFile()) ) {

         write_successful = true;

         const vector<int>& object_ids =
             d_leaf_node_info[leaf_node_id]->getObjectIds();

         int last_file_index = -1;
         toolbox::HDFDatabasePtr write_DB;

         for (unsigned int i = 0; i < object_ids.size(); ++i) {
            int object_id = object_ids[i];

            if ( isValidObjectId(object_id) &&
                 !(d_object_info[object_id]->getInFile()) ) {
         
               string file_name;
               int file_index;
               bool new_file =
                  getFileNameAndIndexForObjectWrite(file_name,
                                                    file_index); 
               string write_DB_name(d_object_file_prefix + file_name);

               if ( file_index != last_file_index ) { 

                  if (write_DB.get()) {
                     write_DB->unmount();
                     write_DB.reset();
                  }

                  write_DB = 
                     toolbox::HDFDatabasePtr(
                        new toolbox::HDFDatabase(write_DB_name) );

                  string flag = ( new_file ? "WN" : "W" );
                  write_successful &= 
                     (write_DB->mount(write_DB_name, flag) >= 0);

               }  // if need to mount new file

               write_successful &=
                  writeDataObjectHelper(object_id,
                                              write_DB,
                                              file_index);

            } // if valid object id and object not already written to file

         } // iterate over objects owned by leaf node

         if (write_DB.get()) {
            write_DB->unmount();
            write_DB.reset(); 
         }

         if (write_successful) {
            d_leaf_node_info[leaf_node_id]->setInFile(true); 
            d_leaf_node_info[leaf_node_id]->setInMemory(false); 
         }

      } else {  // leaf node objects already written to file
         write_successful = true;
      }

   } // if valid leaf node id and all leaf node objects 
     // not already written to file

   } // if data store is open
#else
   TBOX_ERROR("MTreeDataStore error..."
              << "\n   Cannot write data to file..."
              << " code not compiled with HDF" << endl);
#endif

   return(write_successful);
}

/*
*************************************************************************
*                                                                       *
* Read all data objects in leaf node from file                          *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::readNodeDataObjects(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif

   bool read_successful = false;

#ifdef HAVE_PKG_hdf5
   if ( d_is_open ) {   

   const int leaf_node_id = node->getLeafNodeId(); 
   
   if ( isValidLeafNodeId(leaf_node_id) ) {

      if ( !(d_leaf_node_info[leaf_node_id]->getInMemory()) ) {

         read_successful = true;

         const vector<int>& object_ids =
             d_leaf_node_info[leaf_node_id]->getObjectIds();

         int last_file_index = -1;
         toolbox::HDFDatabasePtr read_DB;

         for (unsigned int i = 0; i < object_ids.size(); ++i) {
            int object_id = object_ids[i];

            if ( isValidObjectId(object_id) &&
                 !(d_object_info[object_id]->getInMemory()) ) {

               const int file_index = d_object_info[object_id]->getFileIndex();

               if ( file_index != last_file_index ) { 

                  if (read_DB.get()) {
                     read_DB->unmount();
                     read_DB.reset();
                  }

                  const string read_DB_name( d_object_file_prefix + 
                     d_object_file_info[file_index]->getFileName() );
         
                  read_DB = 
                     toolbox::HDFDatabasePtr(
                        new toolbox::HDFDatabase(read_DB_name) );

                  read_successful &= (read_DB->mount(read_DB_name, "R") >= 0);

               }  // if need to mount new file

               read_successful &=
                  readDataObjectHelper(object_id,
                                       read_DB);

            } // if valid object id and object not already in memory

         } // iterate over objects owned by leaf node

         if (read_DB.get()) {
            read_DB->unmount();
            read_DB.reset(); 
         }

         if (read_successful) {
            d_leaf_node_info[leaf_node_id]->setInMemory(true); 
         }

      } else {  // leaf node objects already in memory
         read_successful = true;
      }

   } // if valid leaf node id and all leaf node objects 
     // do not already exist in memory

   } // if data store is open
#else
   TBOX_ERROR("MTreeDataStore error..."
              << "\n   Cannot read data from file..."
              << " code not compiled with HDF" << endl);
#endif

   return(read_successful);
}

/*
*************************************************************************
*                                                                       *
* Write single data object to file                                      *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::writeDataObject(int object_id)
{
   bool write_successful = false;

   if ( d_is_open ) {   
#ifdef HAVE_PKG_hdf5
//      std::cout << "[@" << toolbox::MPI::getRank() << "] " 
// 	       << "writing out object: " << object_id << std::endl;
   toolbox::HDFDatabasePtr dummy_write_DB;
   int                     dummy_file_index = -1;
   write_successful = writeDataObjectHelper(object_id,
                                            dummy_write_DB,
                                            dummy_file_index);
#else
   TBOX_ERROR("MTreeDataStore::writeDataObject() error..."
              << "\n   Cannot write data to file..."
              << " code not compiled with HDF" << endl);
#endif
   } // if data store is open

   return( write_successful );
}

/*
*************************************************************************
*                                                                       *
* Read single data object from HDF file.                                *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::readDataObject(int object_id)
{
   bool read_successful = false;

   if ( d_is_open ) {   
#ifdef HAVE_PKG_hdf5
     //     std::cout  << "[@" << toolbox::MPI::getRank() << "] " 
     //		<< "reading in object: " << object_id << std::endl;
   toolbox::HDFDatabasePtr dummy_read_DB;
   read_successful = readDataObjectHelper(object_id,
                                          dummy_read_DB);

#else
   TBOX_ERROR("MTreeDataStore::readDataObject() error..."
              << "\n   Cannot read data from file..."
              << " code not compiled with HDF" << endl);
#endif
   } // if data store is open

   return( read_successful );
}

/*
*************************************************************************
*                                                                       *
* Private method to write single data object to HDF file.               *
*                                                                       *
*************************************************************************
*/

#ifdef HAVE_PKG_hdf5

bool MTreeDataStore::writeDataObjectHelper(
   int object_id,
   toolbox::HDFDatabasePtr write_DB,
   int file_index)
{
   bool write_successful = false;

   if ( isValidObjectId(object_id) ) {

      if ( !(d_object_info[object_id]->getInFile()) ) {

         write_successful = true;

         bool local_DB = false;

         if ( !write_DB.get() || (file_index < 0) ) {

            string file_name; 
            bool new_file =
               getFileNameAndIndexForObjectWrite(file_name,
                                                 file_index);
            string write_DB_name(d_object_file_prefix + file_name);
             
            write_DB = 
               toolbox::HDFDatabasePtr( new toolbox::HDFDatabase(write_DB_name) );

            string flag = ( new_file ? "WN" : "W" );
            write_successful &= (write_DB->mount(write_DB_name, flag) >= 0);

            local_DB = true;

         }

         if ( write_successful ) {

            string object_db_name = getObjectDatabaseName(object_id);
            toolbox::DatabasePtr obj_db = 
               write_DB->putDatabase(object_db_name);

            d_object_info[object_id]->getObject()->
                writeToDatabase( *(obj_db.get()) );
 
            obj_db.reset();

            if (local_DB) { 
               write_DB->unmount();
               write_DB.reset(); 
            }

            d_object_info[ object_id ]->setInFile(true);
            d_object_info[ object_id ]->setFileIndex( file_index );
            d_object_info[ object_id ]->resetObjectPtr();
            d_object_info[ object_id ]->setInMemory(false);

            d_object_file_info[file_index]->addObjectId(object_id);

            d_num_objects_in_files++;

         }

      } else {
         write_successful = true;
      }

   }

   return(write_successful);
}

/*
*************************************************************************
*                                                                       *
* Private method to read single data object from HDF file.              *
*                                                                       *
*************************************************************************
*/
 
bool MTreeDataStore::readDataObjectHelper(
   int object_id,
   toolbox::HDFDatabasePtr read_DB,
   bool force_read)
{
   bool read_successful = false; 

   if ( isValidObjectId(object_id) ) {
 
      if ( force_read || (!d_object_info[object_id]->getInMemory()) ) {

         read_successful = true; 
 
         bool mount_successful = false;
 
         if ( !read_DB.get() ) {

            const int file_index = d_object_info[object_id]->getFileIndex();
 
            const string read_DB_name(d_object_file_prefix + 
               d_object_file_info[file_index]->getFileName() );
 
            read_DB = 
               toolbox::HDFDatabasePtr( new toolbox::HDFDatabase(read_DB_name) );
 
            mount_successful = (read_DB->mount(read_DB_name, "R") >= 0);
            read_successful = mount_successful; 
 
         }

         if ( read_successful ) {

            string object_db_name = getObjectDatabaseName(object_id);
            toolbox::DatabasePtr obj_db = read_DB->getDatabase(object_db_name);

            DBObjectPtr data_object = 
               d_object_factory->allocateObject( *(obj_db.get()) );

            data_object->setObjectId(object_id);
            d_object_info[ object_id ]->setObjectPtr(data_object);
            d_object_info[ object_id ]->setInMemory(true);
            d_object_info[ object_id ]->setInFile(false);

	    --d_num_objects_in_files;

         }

         if (mount_successful) {
            read_DB->unmount();
            read_DB.reset();
         } 

      }  // if object is in file

   } // if valid object id

   return(read_successful);
}

#endif

/*
*************************************************************************
*                                                                       *
* Private method to generate file name for next object write.           *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::getFileNameAndIndexForObjectWrite(
   string& file_name,
   int& file_index)
{
   bool new_file = false;


   file_index = d_num_objects_in_files / d_object_file_capacity;
   file_name = getObjectFileName(file_index);
   
   const int num_files = d_object_file_info.size();
   if ( file_index >= num_files ) {
     if ( d_object_file_info.size() == d_object_file_info.capacity() ) {
       d_object_file_info.reserve( d_object_file_info.size() +
				   OBJECT_FILE_VECTOR_ALLOCATION_CHUNK );
     }
     d_object_file_info.push_back( new ObjectFileInfo(file_name) ); 
     new_file = true;
   }

   return( new_file );
}

/*
*************************************************************************
*                                                                       *
* Private method to generate object file name (given file index),       *
* object database name (given object id), mtree file name,              *
* data store file name, and object file prefix.                         *
*                                                                       *
*************************************************************************
*/

string MTreeDataStore::getObjectFileName(int file_index) const
{
   return( d_file_prefix + "__data_objects." + 
           toolbox::Utilities::intToString(file_index, 8) );
}

string MTreeDataStore::getObjectDatabaseName(int object_id) const
{
   return( "data_object__" +
           toolbox::Utilities::intToString(object_id, 9) ); 
}

string MTreeDataStore::getFullPathMTreeIndexFileName() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_directory_name.empty());
   assert(!d_file_prefix.empty());
#endif
   return( d_directory_name + 
           "/" + d_file_prefix + 
           "__mtree_index_structure" );
}

string MTreeDataStore::getFullPathDataStoreFileName() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_directory_name.empty());
   assert(!d_file_prefix.empty());
#endif
   return( d_directory_name + 
           "/" + d_file_prefix + 
           "__data_store_summary" );
}

string MTreeDataStore::getObjectFilePrefix() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!d_directory_name.empty());
#endif
   return( d_directory_name + "/" );
}

/*
*************************************************************************
*                                                                       *
* Private method to write object file information to database.          *
*                                                                       *
*************************************************************************
*/

#ifdef HAVE_PKG_hdf5

bool MTreeDataStore::writeObjectFileInfo(
   toolbox::HDFDatabasePtr write_DB) const
{
   bool write_successful = false;
   write_successful = true;

   const int object_file_info_size = d_object_file_info.size();

   write_DB->putInteger("object_file_info_size", object_file_info_size);
   write_DB->putInteger("d_num_objects_in_files", d_num_objects_in_files);
   write_DB->putInteger("d_object_file_capacity", d_object_file_capacity);

   for (int findex = 0; findex < object_file_info_size; ++findex) {
      string file_db_name = d_object_file_info[ findex ]->getFileName();
      toolbox::DatabasePtr file_obj_db =
                           write_DB->putDatabase(file_db_name);
      file_obj_db->putIntegerArray( "object_ids", 
                                    d_object_file_info[ findex ]->
                                       getObjectIds() );
   }

   return( write_successful );
}

/*
*************************************************************************
*                                                                       *
* Private method to write object information to database.               *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::writeObjectInfo(
   toolbox::HDFDatabasePtr write_DB) const
{
   bool write_successful = false;
   write_successful = true;

   const unsigned int object_info_size = d_object_info.size();

   write_DB->putInteger("object_info_size", object_info_size);
   write_DB->putInteger("d_num_objects", d_num_objects);

   unsigned int object_id = 0;
   int object_count = 0;
   while ( (object_id < object_info_size) &&
           (object_count < d_num_objects) ) {
 
      if ( isValidObjectId(object_id) ) {

         string object_db_name = getObjectDatabaseName(object_id);
         toolbox::DatabasePtr obj_db =
            write_DB->putDatabase(object_db_name);

         obj_db->putInteger( "d_object_file_index", 
                             d_object_info[ object_id ]->getFileIndex() );

         obj_db->putBool( "d_in_memory", 
                           d_object_info[ object_id ]->getInMemory() );

         object_count++;
      } // if valid object id
 
      object_id++;
   }  // iterate over object ids
 
   return( write_successful );
}

#endif

/*
*************************************************************************
*                                                                       *
* Private method to read object file information from database.         *
*                                                                       *
*************************************************************************
*/

#ifdef HAVE_PKG_hdf5

bool MTreeDataStore::readObjectFileInfo(
   toolbox::HDFDatabasePtr read_DB)
{
   bool read_successful = false;
   read_successful = true;
 
   const int object_file_info_size = 
      read_DB->getInteger("object_file_info_size");
   d_num_objects_in_files = read_DB->getInteger("d_num_objects_in_files");
   d_object_file_capacity = read_DB->getInteger("d_object_file_capacity");

   d_object_file_info.clear();
   d_object_file_info.reserve(object_file_info_size);
 
   for (int findex = 0; findex < object_file_info_size; ++findex) {
      string file_db_name = getObjectFileName(findex);

      d_object_file_info.push_back( new ObjectFileInfo(file_db_name) );
      toolbox::DatabasePtr file_obj_db =
                           read_DB->getDatabase(file_db_name);

      vector<int> object_ids;
      file_obj_db->getIntegerArray( "object_ids", object_ids );

      const int num_ids = object_ids.size();
      for (int idi = 0; idi < num_ids; ++idi) {
         d_object_file_info[findex]->addObjectId( object_ids[idi] );
      }
   }

   return( read_successful );
}

/*
*************************************************************************
*                                                                       *
* Private method to read object information from database.              *
*                                                                       *
*************************************************************************
*/

bool MTreeDataStore::readObjectInfo(
   toolbox::HDFDatabasePtr read_DB)
{
   bool read_successful = false;
   read_successful = true;

   const int object_info_size = read_DB->getInteger("object_info_size");
 
   d_num_objects = read_DB->getInteger("d_num_objects");

   d_recycled_object_indices.clear();

   d_object_info.clear();
   d_object_info.resize(object_info_size);
   for (int obj_id = 0; obj_id < object_info_size; ++obj_id) {
      d_object_info[obj_id] = (ObjectInfo*)NULL;
   }
 
   int object_id = 0;
   int object_count = 0;
   while ( (object_id < object_info_size) &&
           (object_count < d_num_objects) ) {
 
      string object_db_name = getObjectDatabaseName(object_id);
      
      if ( read_DB->isDatabase(object_db_name) ) {
      
         toolbox::DatabasePtr obj_db =
            read_DB->getDatabase(object_db_name);

         DBObjectPtr dummy_obj_ptr;

         ObjectInfo* obj_info = new ObjectInfo(dummy_obj_ptr);

         obj_info->setFileIndex( obj_db->getInteger("d_object_file_index") );
         obj_info->setInMemory( obj_db->getBool("d_in_memory") );

         d_object_info[object_id] = obj_info;

         object_count++;

      } else {

         d_recycled_object_indices.push_front(object_id);   

      }
 
      object_id++;

   }  // iterate over object ids
 
   return( read_successful );
}

#endif

/*
*************************************************************************
*                                                                       *
* Methods for private leaf node info utility class.                     *
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::LeafNodeInfo::clearOwnObjectIds()
{
   const int num_objects = d_object_ids.size();
   for (int in = 0; in < num_objects; ++in) {
      d_object_ids[in] = DBObject::getUndefinedId();
   }
}

void MTreeDataStore::LeafNodeInfo::setOwnObjectId(int obj_id, int pos)
{
   const int object_capacity = d_object_ids.capacity();
   int i = d_object_ids.size();
   while ( (i <= pos) && (i < object_capacity) ) {
      d_object_ids.push_back(DBObject::getUndefinedId());
      i++;
   }
   d_object_ids[pos] = obj_id;
}

void MTreeDataStore::LeafNodeInfo::printClassData(ostream& stream) const
{
   if (d_node_ptr.get()) {
      stream << "node id = " << d_node_ptr->getNodeId() << endl;
      stream << "leaf node id = " << d_node_ptr->getLeafNodeId() << endl;
   } else {
      stream << "d_node_ptr is NULL" << endl;
   }
   stream << "d_in_file = " << d_in_file << endl;
   stream << "d_in_memory = " << d_in_memory << endl;
   for (unsigned int in = 0; in < d_object_ids.size(); ++in) {
      stream << "   d_object_ids[" << in << "] = "
             << d_object_ids[in] << endl;
   }
}

/*
*************************************************************************
*                                                                       *
* Methods for private object info utility class.                        *
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::ObjectInfo::printClassData(ostream& stream) const
{
   if (d_object_ptr.get()) {
      stream << "object id = " << d_object_ptr->getObjectId() << endl;
   } else {
      stream << "d_object_ptr is NULL" << endl;
   }
   stream << "d_owner_leaf_node_id = " << d_owner_leaf_node_id << endl;
   stream << "d_owner_leaf_node_entry_position = " 
          << d_owner_leaf_node_entry_position << endl;
   stream << "d_in_file = " << d_in_file << endl;
   stream << "d_in_memory = " << d_in_memory << endl;
   stream << "d_object_file_index = " << d_object_file_index << endl;
}

/*
*************************************************************************
*                                                                       *
* Methods for private object file info utility class.                   *
*                                                                       *
*************************************************************************
*/
 
void MTreeDataStore::ObjectFileInfo::printClassData(ostream& stream) const
{
   stream << "d_file_name = " << d_file_name << endl;
   for (unsigned int in = 0; in < d_object_ids.size(); ++in) {
      stream << "   d_object_ids[" << in << "] = "
             << d_object_ids[in] << endl;
   }
}

/*
*************************************************************************
*                                                                       *
* Write data store contents to given output stream.                     *
*                                                                       *
*************************************************************************
*/

void MTreeDataStore::printClassData(ostream& stream) const
{
   stream << "MTreeDataStore::printClassData()\n";
   stream << "--------------------------------------\n";
   stream << "this ptr = " << (MTreeDataStore*)this << endl;
   stream << "d_mtree = " << (MTree*)d_mtree << endl;
   stream << "d_object_factory = " 
          << (DBObjectFactory*)d_object_factory << endl;
   stream << "d_is_initialized = " << d_is_initialized << endl;
   stream << "d_is_open = " << d_is_open << endl;
   stream << "d_directory_name = " << d_directory_name << endl;
   stream << "d_mtree_index_file_name = " 
          << d_mtree_index_file_name << endl;
   stream << "d_data_store_file_name = " 
          << d_data_store_file_name << endl;
   stream << "d_object_file_prefix = " << d_object_file_prefix << endl;
   stream << "d_object_file_capacity = " << d_object_file_capacity << endl;
   stream << "d_num_objects_in_files = " << d_num_objects_in_files << endl;

   stream << "\n\n" << endl;
   stream << "d_num_leaf_nodes = " << d_num_leaf_nodes << endl;
   stream << "d_leaf_node_info size = " 
          << d_leaf_node_info.size() << endl;
   stream << "d_leaf_node_info capacity = " 
          << d_leaf_node_info.capacity() << endl;
   stream << "d_recycled_leaf_node_indices size = "
          << d_recycled_leaf_node_indices.size() << endl;
#if 0
   stream << "d_recycled_leaf_node_indices = \n    ";
   list<int>::iterator lnit(d_recycled_leaf_node_indices.begin());
   while ( lnit != d_recycled_leaf_node_indices.end() ) {
      stream << *lnit << " , ";
      ++lnit;
   }
#endif
   stream << "\n" << endl;

   stream << "Leaf Node Info ..." << endl;
   unsigned int lni = 0;
   int node_count = 0; 
   while ( (lni < d_leaf_node_info.size()) &&
           (node_count < d_num_leaf_nodes) ) {
      if (d_leaf_node_info[lni]) {
         stream << "   d_leaf_node_info[" << lni << "] = " << endl;
         d_leaf_node_info[lni]->printClassData(stream);
         stream << endl;
         node_count++;
      } else {
         stream << "   d_leaf_node_info[" << lni << "] = NULL\n" 
                << endl;
      }
      lni++; 
   }

   stream << "\n\n" << endl;
   stream << "d_num_objects = " << d_num_objects << endl;
   stream << "d_object_info size = " 
          << d_object_info.size() << endl;
   stream << "d_object_info capacity = " 
          << d_object_info.capacity() << endl;
   stream << "d_recycled_object_indices size = "
          << d_recycled_object_indices.size() << endl;
#if 0
   stream << "d_recycled_object_indices = \n    ";
   list<int>::iterator oit(d_recycled_object_indices.begin());
   while ( oit != d_recycled_object_indices.end() ) {
      stream << *oit << " , ";
      ++oit;
   }
#endif
   stream << "\n" << endl;

   stream << "Object Info ..." << endl;
   unsigned int oi = 0;
   int object_count = 0; 
   while ( (oi < d_object_info.size()) &&
           (object_count < d_num_objects) ) {
      if (d_object_info[oi]) {
         stream << "   d_object_info[" << oi << "] = " << endl;
         d_object_info[oi]->printClassData(stream);
         object_count++;
      } else {
         stream << "   d_object_info[" << oi << "] = NULL" << endl;
      }
      stream << endl;
      oi++; 
   }

   stream << "\n\n" << endl;
   stream << "d_object_file_info size = " 
          << d_object_file_info.size() << endl; 
   stream << "d_object_file_info capacity = " 
          << d_object_file_info.capacity() << endl; 
   stream << "\n" << endl;

   stream << "Object File Info ..." << endl;
   for (unsigned int fi = 0; fi < d_object_file_info.size(); ++fi) {
      stream << "   d_object_file_info[" << fi << "] = " << endl;
      d_object_file_info[fi]->printClassData(stream);
      stream << endl;
   }

}


#endif




