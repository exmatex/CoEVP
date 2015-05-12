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
// File:        MTree.cc
// Package:     MTree database
// Description: Main Mtree index structure class.
//

#ifndef included_MTree_C
#define included_MTree_C

#include "MTree.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif
#endif

#ifndef included_list
#define included_list
#include <list>
using namespace std;
#endif

#include "toolbox/base/Utilities.h"
#include "toolbox/base/MathUtilities.h"

#include "MTreeSearchNode.h"
#include "MTreeSearchQueue.h"

#ifdef DEBUG_NO_INLINE
#include "MTree.I"
#endif

/*
*************************************************************************
*                                                                       *
* MTree ctor and dtor.                                                  *
*                                                                       *
*************************************************************************
*/

MTree::MTree(const string& tree_name,
             ostream* error_log_stream,
             bool do_error_checking)
   : DB(tree_name, error_log_stream, do_error_checking),
     d_max_node_entries(DEFAULT_MAX_NODE_ENTRIES),
     d_root_node_promotion_method(MTreeNode::MIN_OVERLAP_PROMOTION),
     d_node_promotion_method(MTreeNode::MAX_SPREAD_DISTANCE_PROMOTION),
     d_node_partition_method(MTreeNode::HYPERPLANE_PARTITION),
     d_min_node_utilization(0.5)
{
// For calculating tree statistics, set up array of 
// node count per level sized to some large number 
// of levels we hope to never exceed 
   const int big_size = 50;
   d_number_nodes_in_level.resize(big_size); 
   for (int i = 0; i < big_size; ++i) {
      d_number_nodes_in_level[i] = 0;
   }
   clearLevelStatistics();
   clearOperationCountStatistics();
}

MTree::~MTree()
{
   finalize();

   clearLevelStatistics();
   clearOperationCountStatistics();

   d_root_node.reset();
}

/*
*************************************************************************
*                                                                       *
* MTree initialization routines: initializeCreate() creates a new       *
* empty data store and initializeOpen() initializes the data store to   *
* the contents of an exisiting data file.                               *
*                                                                       *
*************************************************************************
*/

void MTree::initializeCreate(const string& directory_name,
                             const string& file_prefix,
                             const DBObjectFactory& obj_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!directory_name.empty());
   assert(!file_prefix.empty());
#endif

   if (!d_data_store.isInitialized()) {
      d_data_store.create(this,
                          &obj_factory,
                          directory_name,
                          file_prefix);
   } else {
      TBOX_WARNING("MTree::initializeCreate() warning"
                   << " for tree named = " << d_db_name
                   << "\nCannot call initializeCreate() method on tree"
                   << " that is already initialized!" 
                   << "\n So calling initialization method did nothing!"
                   << endl);
   }

   clearLevelStatistics();
   clearOperationCountStatistics();

}

void MTree::initializeOpen(const string& directory_name,
                           const string& file_prefix,
                           const DBObjectFactory& obj_factory)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!directory_name.empty());
   assert(!file_prefix.empty());
#endif

   if (!d_data_store.isInitialized()) {
      d_data_store.open(this,
                        &obj_factory,
                        directory_name,
                        file_prefix);
   } else {
      TBOX_WARNING("MTree::initializeOpen() warning"
                   << " for tree named = " << d_db_name
                   << "\nCannot call initializeOpen() method on tree"
                   << " that is already initialized!"
                   << "\n So calling initialization method did nothing!"
                   << endl);
   }

}

/*
*************************************************************************
*                                                                       *
* MTree finalize routine.                                               *
*                                                                       *
*************************************************************************
*/

void MTree::finalize() 
{
   d_data_store.close();
}

/*
*************************************************************************
*                                                                       *
* Set maximum entries for each node in the MTree structure.             *
*                                                                       *
*************************************************************************
*/

void MTree::setMaxNodeEntries(int max_entries)
{
  if ( d_root_node || (max_entries < 2) ) {
     TBOX_WARNING("MTree::setMaxNodeEntries() warning"
                  << " for tree named = " << d_db_name
                  << "\nEither the given max number of entries = " 
                  << max_entries << " is less than 2,"
                  << "\nor some data object has been inserted already."
                  << "\nChanging the max number of entries under either"
                  << "\n of these situations is not allowed."
                  << "\nSo max node entries has not be changed from"
                  << "\npreviously set value = " << d_max_node_entries
                  << endl);
  } else {
     d_max_node_entries = max_entries; 
  }
}


/*
*************************************************************************
*                                                                       *
* Main routine for inserting an object into the MTree index structure.  *
* Note: deep copy of object and point are made for insertion into tree. *
*                                                                       *
*************************************************************************
*/
 
void MTree::insertObject(DBObject& object,
                         const MetricSpacePoint& point,
                         double radius)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(radius >= 0.0);
#endif

   if (!d_data_store.isInitialized()) {
      TBOX_ERROR("MTree object error: " << d_db_name
                 << "\n   an initialize method must be called"
                 << " before calling MTree::insertObject()" << endl);
   } else {

     d_num_distance_comps_in_last_insert = 0;
     d_current_distance_type = DISTANCE_INSERT;
     d_num_inserts++;

     if (!d_root_node) {
        createRootNode();
     }

     d_data_store.addObject(object);  

     MTreeKey key(point.makeCopy(), radius);

     MTreeEntryPtr entry( new MTreeEntry(key) );
     entry->setDataObjectId( object.getObjectId() );

     MTreeNodePtr node( pickNode(d_root_node, entry) );

     MTreeNode::insertEntry(entry, node);

     clearLevelStatistics();
 
     if (node->isOverfull()) {

        split(node);

        if (node->isOverfull()) {  // just to be anal....
           TBOX_ERROR("MTree::insertObject() error"
                      << " for tree named = " << d_db_name
                      << "\nSplit of overfull node id " 
                      << node->getNodeId()
                      << "\nfailed for object id " 
                      << entry->getDataObjectId()
                      << endl);
        }

     } else {

        mapDataObjectsToNode(node);

        node->resetRadiiUpToRoot();

     }
 
     if (d_do_error_checking) {
        *d_error_log_stream 
        << "Checking correctness from insert node up to root"
        << "\nafter inserting object # " << object.getObjectId()
        << " ..." << endl;
        if ( node->checkConsistencyUpToRoot(*d_error_log_stream) ) {
           *d_error_log_stream << "      Everything looks good!" << endl; 
        } else {
           *d_error_log_stream << "      Problems found!" << endl; 
        }
     }

     d_current_distance_type = DISTANCE_UNDEFINED;

   }

}


/*
*************************************************************************
*                                                                       *
* Routine for retrieving an object from the MTree by object id.         *
*                                                                       *
*************************************************************************
*/

DBObjectPtr MTree::getObject(int object_id) const
{
   DBObjectPtr ret_object( d_data_store.getObjectCopy(object_id) );
   return(ret_object);
}

/*
*************************************************************************
*                                                                       *
* Main routine for deleting an object from the MTree index structure.   *
*                                                                       *
*************************************************************************
*/
 
void MTree::deleteObject(int object_id)
{

   if ( d_data_store.isValidObjectId(object_id) ) {

      clearLevelStatistics();

      d_num_distance_comps_in_last_delete = 0;
      d_current_distance_type = DISTANCE_DELETE;
      d_num_deletes++;

      MTreeNodePtr node( d_data_store.getNodeOwningObject(object_id) );

      if ( node.get() && node->isLeaf() ) {

         MTreeEntryPtr entry;

         const int num_entries = node->getNumberEntries();
         bool found = false;
         int ie = 0;
         while ( !found && (ie < num_entries) ) {
            entry = node->getEntry(ie);
            if ( object_id == entry->getDataObjectId() ) {
               found = true;
            }  
            ie++;
         }

         if ( found ) {

            MTreeNode::deleteEntry(entry, node);
            d_data_store.removeObject(object_id);
            mapDataObjectsToNode(node);

            /* 
             * If the node is now empty but is not the root in the tree, 
             * we delete the node from the tree.
             * 
             * If the parent node contains only one other child node 
             * apart from the empty node, we attempt to move some entries 
             * from the other child to the empty node.  This helps to 
             * maintain some fanout in the tree.  If the other 
             * child node doesn't have enough entries to share (i.e.,
             * it has one or zero entries) then we cannot do this and
             * so the parent will have only one child node.
             *
             * If the parent node is now empty, we recurse and apply
             * the same algorithm to the parent.  We stop at the root
             * or the first non-empty parent node.
             */

            while ( !node->isRoot() &&
                    (node->getNumberEntries() == 0) ) {

               MTreeNodePtr  parent_node( node->getParentNode() );
               MTreeEntryPtr parent_entry( node->getParentEntry() );

               if ( parent_node->getNumberEntries() == 1 ) {

                  if ( node->isLeaf() ) {
                     d_data_store.removeLeafNode(node);
                     removeNodeFromLevelCount( node->getLevelInTree() );
                  }
                  MTreeNode::deleteEntry(parent_entry, parent_node);

               } else {

                  if ( parent_node->getNumberEntries() == 2 ) {  

                     MTreeNodePtr neighbor_node;
                     if ( parent_entry->getPositionInNode() == 0 ) {
                        neighbor_node = 
                           parent_node->getEntry(1)->getSubtreeNode();
                     } else {
                        neighbor_node = 
                           parent_node->getEntry(0)->getSubtreeNode();
                     }

                     if ( node->isLeaf() ) {
                        d_data_store.removeLeafNode(node);
                        removeNodeFromLevelCount( node->getLevelInTree() );
                     }
                     MTreeNode::deleteEntry(parent_entry, parent_node);

                     if ( neighbor_node->getNumberEntries() >= 2) {

                        parent_node->resetRadius();

                        bool use_root_promotion = false;
                        bool special_delete_split = true;
                        MTreeNodePtr new_node(
                           neighbor_node->split( use_root_promotion,
                                                 special_delete_split ) );

                        if (neighbor_node->isLeaf()) {
                           mapDataObjectsToNode(neighbor_node);
                           mapDataObjectsToNode(new_node);
                        } else {
                           new_node->
                              setLevelInTree( neighbor_node->getLevelInTree() );
                        }

                     } // if neighbor node can be split

                  } else { // parent node has more than 2 entries;
                           // simply remove empty child node

                     if ( node->isLeaf() ) { 
                        d_data_store.removeLeafNode(node);
                        removeNodeFromLevelCount( node->getLevelInTree() );
                     }
                     MTreeNode::deleteEntry(parent_entry, parent_node);

                  }

               }  // else parent node has at least two entries

               parent_node->resetRadius(); 

               if ( parent_node->getNumberEntries() == 0 ) {
                  parent_node->
                     setLevelInTree( parent_node->getLevelInTree() - 1 );
               }

               node = parent_node;

            } // while node to prune is empty but is not root node 

            if ( !node->isRoot() ) {
               node->resetRadiiUpToRoot();
            } 

         } // object's entry found in tree
                 
      } else { // object not assigned to node in store, 
               // so just attempt to remove object from data store.
         d_data_store.removeObject(object_id);
      }

      d_current_distance_type = DISTANCE_UNDEFINED;
      
   }  // if valid object id in data store

}

/*
*************************************************************************
*                                                                       *
* Search tree and return up to "K" nearest neighbors to query point.    * 
*                                                                       *
*************************************************************************
*/

void MTree::searchKNN(vector<DBSearchResult>& results,
                      const MetricSpacePoint& query_point,
                      int k_neighbors,
                      bool make_safe)
{
   if (k_neighbors > 0) {

      d_num_distance_comps_in_last_knn_query = 0;
      d_current_distance_type = DISTANCE_KNN_QUERY;
      d_num_knn_queries++;

      /*
       * Set up query to compare neighbors, results vector,
       * and queue of nodes to search.
       */

      const double max_dist = MetricSpacePoint::getMaxDistance();
      const int    klast    = k_neighbors - 1;

      MTreeQuery query( query_point.makeCopy(), max_dist, this );

      results.clear();
      results.reserve(k_neighbors);
      for (int irt = 0; irt < k_neighbors; ++irt) {
         MTreeSearchResult result;
         result.setQueryPoint( query.getPoint() );
         result.setDistanceToQueryPoint( max_dist );
         results.push_back(result); 
      }

      MTreeSearchQueue queue;

      /*
       * Start at root node and search for nearest neighbors.
       */

      MTreeNodePtr node( d_root_node );
      while ( node.get() ) {

         double old_dist = query.getGrade();

         const int num_entries = node->getNumberEntries();
         for (int ie = 0; ie < num_entries; ++ie) {
            MTreeEntryPtr entry(node->getEntry(ie));

            query.setGrade(old_dist);

            if ( query.isConsistent(entry) ) {

               if ( entry->isDataEntry() ) {

                  if ( query.getGrade() < 
                       results[klast].getDistanceToQueryPoint() ) {

                     int id = 0;
                     while ( query.getGrade() > 
                             results[id].getDistanceToQueryPoint() ) {
                        id++;
                     }

                     for (int ir = klast; ir > id; ir--) {
                        results[ir] = results[ir-1];
                     }

                     MTreeSearchResult result( query.getPoint(), entry );
                     result.setDistanceToQueryPoint( query.getGrade() );
                     results[id] = result;

                  }

                  query.setRadius( results[klast].getDistanceToQueryPoint() );

               } else {  // routing entry; insert child node into search queue  
                  double dmin = 
                     toolbox::MathUtilities<double>::Max( 0.0,
                        query.getGrade() - entry->getRadius() );
              
                  MTreeSearchNode search_node( entry->getSubtreeNode(),
                                               dmin,
                                               query.getGrade() );

                  queue.insert(search_node);

               } // else

            }  // if entry is consistent with query

         }  // iterate over entries in node

         /*
          * If queue is not empty, pull off the next node 
          * to search from the front of the queue.  If the
          * queue is empty, we are done and we terminate the 
          * search. 
          * 
          * The node we pull off the queue will be the subtree 
          * whose root is next closest to query point.  If this node 
          * is farther from the query point than the last item in  
          * the vector of results, we are done and we terminate the 
          * search.  Otherwise, we continue to search for data objects 
          * closer than the ones we have collected so far.
          */

         if ( queue.empty() ) { 
            node.reset();
         } else {
            MTreeSearchNode next_node = queue.getFirst();
            queue.removeFirst();
 
            if ( next_node.getBound() >= 
                 results[klast].getDistanceToQueryPoint() ) {
               queue.clear();
               node.reset();
            } else {
               node = next_node.getSearchNode();
               query.setGrade(next_node.getDistance());
            }
         } // else queue not empty

      }  // while still searching nodes

      /*
       * Search complete, finalize results to return.
       */
      for (int iresult = klast; iresult >= 0; iresult--) {
         if ( results[iresult].isValidResult() ) {
            ((MTreeSearchResult&)results[iresult]).finalizeSearchResult(d_data_store,
                                                  make_safe);
         } else {
            results.pop_back();
         }
      }

      d_current_distance_type = DISTANCE_UNDEFINED;

   } // if number of neighbors to find is > 0
}

/*
*************************************************************************
*                                                                       *
* Search tree and return all data entries whose regions have non-empty  *
* intersection with region given by query point and radius.             *
*                                                                       *
*************************************************************************
*/

void MTree::searchRange(list<DBSearchResult>& results,
                        const MetricSpacePoint& query_point,
                        double radius,
                        bool make_safe)
{
   if (radius >= 0.0) {

      d_num_distance_comps_in_last_range_query = 0;
      d_current_distance_type = DISTANCE_RANGE_QUERY;
      d_num_range_queries++;

      list<MTreeSearchResult> tmp_results;

      //
      // check if there exists a root node
      //

      if (d_root_node == NULL)
	return;

      //
      //
      //

      MTreeQuery query( query_point.makeCopy(), radius, this );

      searchRangeRecursive(tmp_results, query, d_root_node);

      /*
       * Search complete, finalize results to return.
       */
      results.clear();
      while ( !tmp_results.empty() ) {
          MTreeSearchResult result(tmp_results.front());
          tmp_results.pop_front();

          if ( result.isValidResult() ) {
             result.finalizeSearchResult(d_data_store,
                                         make_safe);
             results.push_back(result);
          }

      }

      results.sort();

      d_current_distance_type = DISTANCE_UNDEFINED;

   }
}

/*
*************************************************************************
*                                                                       *
* Private recursive member function used in rage search.                *
*                                                                       *
*************************************************************************
*/

void MTree::searchRangeRecursive(list<MTreeSearchResult>& results,
                                 const MTreeQuery& query,
                                 MTreeNodePtr node) const
{

   if ( node->isLeaf() ) {

      const int num_entries = node->getNumberEntries();
      for (int ie = 0; ie < num_entries; ++ie) {
         MTreeEntryPtr entry(node->getEntry(ie));
         MTreeQuery q_tmp(query); 
         if ( q_tmp.isConsistent(entry) ) {
            MTreeSearchResult result( q_tmp.getPoint(), entry );
            result.setDistanceToQueryPoint(q_tmp.getGrade());
            results.push_back(result);
         }
      }

   } else {

      const int num_entries = node->getNumberEntries();
      for (int ie = 0; ie < num_entries; ++ie) {
         MTreeEntryPtr entry(node->getEntry(ie));
         MTreeQuery q_tmp(query);
         if ( q_tmp.isConsistent(entry) ) {
            list<MTreeSearchResult> e_results;
            searchRangeRecursive(e_results,
                                 query,
                                 entry->getSubtreeNode()); 
            if (!e_results.empty()) {
               results.insert(results.end(), e_results.begin(), e_results.end()); 
               e_results.clear();
            }
         }
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Private member function to choose node in which to insert an entry.   *
* The method accepts a starting node and traverses the tree from that   *
* node to a suitable leaf node, which it returns.                       *
*                                                                       *
*************************************************************************
*/

MTreeNodePtr MTree::pickNode(MTreeNodePtr start_node,
                             MTreeEntryPtr entry) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(start_node.get());
   assert(entry.get());
#endif

   MTreeNodePtr insert_node(start_node);

   bool done = false;
   while (!done) {

      done = ( insert_node->isLeaf() );
 
      if ( !done ) {

         insert_node = insert_node->searchBestSubtreeForInsert(entry);  

         if ( !insert_node.get() ) {
            TBOX_ERROR("MTree::pickNode() error"
                       << " for tree named = " << d_db_name
                       << "\nFailed to find node to insert entry for object "
                       << entry->getDataObjectId() 
                       << endl);
         }

      }
 
   }

   return(insert_node);

}

/*
*************************************************************************
*                                                                       *
* Private function to split a given node if it is over-full. The split  *
* procedure adds a new entry and associated subtree node to the parent  *
* of the givne node.  Then, we travers the tree to the parent node and  *
* split it if it is also over-full.  This process continues up the tree *
* to the root node until no parent node requires splitting.  If we      *
* reach the root and it is over-full, the root node is split, which     *
* increments the height of the tree (i.e., the level of root node).     *
*                                                                       *
*************************************************************************
*/

void MTree::split(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif

   MTreeNodePtr node2split(node);

   while ( node2split->isOverfull() && 
           !node2split->isRoot() ) {

      MTreeNodePtr new_node(node2split->split());

      if (node2split->isLeaf()) {
         mapDataObjectsToNode(node2split);
         mapDataObjectsToNode(new_node);
      } else {
         new_node->setLevelInTree( node2split->getLevelInTree() ); 
      } 

      if (d_do_error_checking) {
         *d_error_log_stream 
            << "Checking correctness after splitting node # "
            << node2split->getNodeId() << endl; 
         *d_error_log_stream 
            << "   Checking split node..." << endl;
         if ( node2split->checkConsistency(*d_error_log_stream) ) {
            *d_error_log_stream << "      Everything looks good!" << endl; 
         } else {
            *d_error_log_stream << "      Problems found!" << endl; 
         }
         *d_error_log_stream 
            << "\n   Checking new node # " 
            <<  new_node->getNodeId() << endl;
         if ( new_node->checkConsistency(*d_error_log_stream) ) {
            *d_error_log_stream << "      Everything looks good!" << endl; 
         } else {
            *d_error_log_stream << "      Problems found!" << endl; 
         }
      }

      node2split = node2split->getParentNode();

   }

   if ( node2split->isRoot() ) {
      if ( node2split->isOverfull() ) { 
         splitRootNode();
      }
   } else {
      node2split->resetRadiiUpToRoot();
   }


}

/*
*************************************************************************
*                                                                       *
* Split root node if over-full.  Note that we first create a new node   *
* and transfer all entries in the tree to this new node.  Then, we      *
* randomly promote one of the entries to the root node.  This allows    *
* us to then call the standard node split routine on the new node.      *
*                                                                       *
*************************************************************************
*/

void MTree::splitRootNode()
{
   if ( d_root_node->isOverfull() ) {

      const bool root_was_leaf = d_root_node->isLeaf();

      MTreeNodePtr new_node( new MTreeNode(this, d_max_node_entries) ); 

      d_root_node->transferEntries(new_node);

      if (root_was_leaf) {
         d_data_store.removeLeafNode(d_root_node);
      }

      MTreeEntryPtr tmp_root_entry;
      new_node->promoteOne(tmp_root_entry, 
                           MTreeNode::CONFIRMED_RANDOM_PROMOTION);

      MTreeEntryPtr root_entry( new MTreeEntry(tmp_root_entry->getKey()) ); 
      root_entry->setSubtreeNode(new_node);

      MTreeNode::insertEntry(root_entry, d_root_node);
      new_node->setParentEntry(root_entry);

      bool use_root_promotion = true;
      MTreeNodePtr new_node2( new_node->split(use_root_promotion) );

      if (root_was_leaf) {
         mapDataObjectsToNode(new_node);
         mapDataObjectsToNode(new_node2);
      } 

      d_root_node->setLevelInTree( d_root_node->getLevelInTree()+1 ); 

      new_node->setLevelInTree( d_root_node->getLevelInTree()-1 );
      new_node2->setLevelInTree( d_root_node->getLevelInTree()-1 );

      if (d_do_error_checking) {
         *d_error_log_stream
            << "Checking correctness after splitting root node " << endl;
         *d_error_log_stream
            << "   Checking new root child node # " 
            << new_node->getNodeId() << "..." << endl;
         if ( new_node->checkConsistency(*d_error_log_stream) ) {
            *d_error_log_stream << "      Everything looks good!" << endl;
         } else {
            *d_error_log_stream << "      Problems found!" << endl;
         }
         *d_error_log_stream
            << "   Checking new root child node # " 
            << new_node2->getNodeId() << "..." << endl;
         if ( new_node2->checkConsistency(*d_error_log_stream) ) {
            *d_error_log_stream << "      Everything looks good!" << endl;
         } else {
            *d_error_log_stream << "      Problems found!" << endl;
         }
         *d_error_log_stream
            << "   Checking root node..." << endl;
         if ( d_root_node->checkConsistency(*d_error_log_stream) ) {
            *d_error_log_stream << "      Everything looks good!" << endl;
         } else {
            *d_error_log_stream << "      Problems found!" << endl;
         }
      }

   }
}

/*
*************************************************************************
*                                                                       *
* Allocate and initialize root node of tree.                            *
*                                                                       *
*************************************************************************
*/

void MTree::createRootNode()
{
   MTreeNodePtr root( new MTreeNode(this, d_max_node_entries) );
   d_root_node = root;
   d_root_node->setRootNode(true);
   d_root_node->setLevelInTree(0);
}

/*
*************************************************************************
*                                                                       *
* Set node ownership in data store for objects in node, if a leaf.      *
*                                                                       *
*************************************************************************
*/

void MTree::mapDataObjectsToNode(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif

   if ( node->isLeaf() ) {
      d_data_store.addLeafNode(node);
      bool operation_successful = d_data_store.clearOwnObjectIds(node);
      const int num_entries = node->getNumberEntries();
      for (int ie = 0; ie < num_entries; ++ie) {
         operation_successful &= 
            d_data_store.setNodeOwningObject(node, 
                                             node->getEntry(ie)); 
      }
      if (!operation_successful) {
         TBOX_ERROR("MTree::mapDataObjectsToNode() error"
                    << " for tree named = " << d_db_name
                    << "\nSome operation in mapping objects to node # " 
                    << node->getNodeId() << " failed!" 
                    << endl); 
      }
   } else {
      TBOX_ERROR("MTree::mapDataObjectsToNode() error"
                 << " for tree named = " << d_db_name
                 << "\nMethod called for node # " << node->getNodeId()
                 << "\nwhich is not marked as a leaf node in tree."
                 << endl);
   }

}

/*
*************************************************************************
*                                                                       *
* Check consistency of tree structure.                                  *
*                                                                       *
*************************************************************************
*/
 
bool MTree::checkConsistency(ostream& stream) const
{
   bool tree_is_consistent = true;

   stream << "\n Checking consistency of entire tree...." << endl;

   list<MTreeNodePtr> nodes2check; 
   nodes2check.push_front(d_root_node);

   while (!nodes2check.empty() ) {

      MTreeNodePtr check_node = nodes2check.front();
      nodes2check.pop_front();

      tree_is_consistent &= check_node->checkConsistency(stream);

      if ( !check_node->isLeaf() ) {
         for (int ie = check_node->getNumberEntries()-1; ie >= 0; ie--) {
            if ( check_node->getEntry(ie)->isRoutingEntry() ) {
               nodes2check.push_front(
                  check_node->getEntry(ie)->getSubtreeNode());
            } 
         }
      } 

   }

   if ( tree_is_consistent ) {
      stream << "   Everything looks good!" << endl;
   } else {
      stream << "   Problems found!" << endl;
   }
 
   return( tree_is_consistent );
}

/*
*************************************************************************
*                                                                       *
* Print node summary to given output stream.                            *
*                                                                       *
*************************************************************************
*/

void MTree::printNodeSummary(ostream& stream) const
{
   stream << "\nMTree::printNodeSummary()\n";
   stream << "d_db_name = " << d_db_name << endl;
   stream << "d_max_node_entries = " << d_max_node_entries << endl;
   const int number_levels_in_tree = getNumberLevels();
   stream << "number levels in tree = " << number_levels_in_tree << endl;
   for (int ln = 0; ln < number_levels_in_tree; ++ln) {
      stream << "   " << d_number_nodes_in_level[ln]
             << " nodes in level " << ln << endl;
   }
   stream << "--------------------------------------\n";

   if (number_levels_in_tree > 0) {

      list<MTreeNodePtr> nodes2print;
      nodes2print.push_front(d_root_node);

      while (!nodes2print.empty() ) {

         MTreeNodePtr print_node = nodes2print.front();
         nodes2print.pop_front();

         print_node->printSummary(stream);
         stream << endl;

         if ( !print_node->isLeaf() ) {
            for (int ie = print_node->getNumberEntries()-1; ie >= 0; ie--) {
               if ( print_node->getEntry(ie)->isRoutingEntry() ) {
                  nodes2print.push_front(
                     print_node->getEntry(ie)->getSubtreeNode());
               }
            }
         }

      } // end while

   }  // if > 0 levels in tree

}

/*
*************************************************************************
*                                                                       *
* Print entire tree structure to given output stream.                   *
*                                                                       *
*************************************************************************
*/

void MTree::printClassData(ostream& stream) const
{
   stream << "\nMTree::printClassData()\n";
   stream << "--------------------------------------\n";
   stream << "this ptr = " << (MTree*)this << endl;
   stream << "d_db_name = " << d_db_name << endl;
   stream << "d_max_node_entries = " << d_max_node_entries << endl;
   stream << "d_do_error_checking = " << d_do_error_checking << endl;
   stream << "d_root_node_promotion_method = " 
          << d_root_node_promotion_method << endl;
   stream << "d_node_promotion_method = " << d_node_promotion_method << endl;
   stream << "d_node_partition_method = " << d_node_partition_method << endl;
   stream << "d_min_node_utilization = " << d_min_node_utilization << endl;
   printAllTreeOperationStatistics(stream);

   const int number_levels_in_tree = getNumberLevels();
   stream << "number levels in tree = " << number_levels_in_tree << endl;
   for (int ln = 0; ln < number_levels_in_tree; ++ln) {
      stream << "   " << d_number_nodes_in_level[ln]
             << " nodes in level " << ln << endl;
   }

   if (number_levels_in_tree > 0) {

      stream << "\nPrinting nodes in tree..." << endl;

      list<MTreeNodePtr> nodes;
      nodes.push_front(d_root_node);

      while (!nodes.empty() ) {

         MTreeNodePtr print_node(nodes.front());
         nodes.pop_front();

         print_node->printClassData(stream);

         if ( !print_node->isLeaf() ) {
            for (int ie = print_node->getNumberEntries()-1; ie >= 0; ie--) {
               if ( print_node->getEntry(ie)->isRoutingEntry() ) {
                  nodes.push_front(
                     print_node->getEntry(ie)->getSubtreeNode());
               }
            }
         }

      }   // end while

   }  // if > 0 levels in tree
   
}

/*
*************************************************************************
*                                                                       *
* Print data store contents to given output stream.                     *
*                                                                       *
*************************************************************************
*/

void MTree::printDataStoreClassData(ostream& stream) const
{
   d_data_store.printClassData(stream);
}

/*
*************************************************************************
*                                                                       *
* Functions for computing and obtaining tree statistics.                * 
*                                                                       *
*************************************************************************
*/

void MTree::calculateLevelStatistics()
{

   clearLevelStatistics();

   /*
    * Allocate statistics information for current tree structure. 
    */

   const int num_levels = getNumberLevels();
   d_level_statistics.reserve(num_levels);
   
   for (int ln = 0; ln < num_levels; ++ln) {
      MTreeLevelStatistic* lstat = 
         new MTreeLevelStatistic(ln, d_number_nodes_in_level[ln]);
      d_level_statistics.push_back(lstat);
   }

   /*
    * Setup statistics information for current tree structure by level. 
    */

   if (num_levels > 0) {

      list<MTreeNodePtr> nodes;
      nodes.push_front(d_root_node);

      while (!nodes.empty() ) {

         MTreeNodePtr stat_node(nodes.front());
         nodes.pop_front();

         d_level_statistics[stat_node->getLevelInTree()]->
            addNodeStat(stat_node);

         if ( !stat_node->isLeaf() ) {
            for (int ie = stat_node->getNumberEntries()-1; ie >= 0; ie--) {
               if ( stat_node->getEntry(ie)->isRoutingEntry() ) {
                  nodes.push_front(
                     stat_node->getEntry(ie)->getSubtreeNode());
               }
            }
         }

      }   // end while

      /*
       * Accumulate total data object in each subtree.
       */

      for (int ln = 1; ln < num_levels; ++ln) {
         MTreeLevelStatistic* l_stats = d_level_statistics[ln];  
         const int n_nodes = l_stats->getNumberNodesOnLevel();

         MTreeLevelStatistic* sub_l_stats = d_level_statistics[ln-1];
     
         for (int ni = 0; ni < n_nodes; ++ni) {
            vector<int> sub_node_real_ids;
            l_stats->getSubtreeNodeRealIds(ni, sub_node_real_ids);
           
            int n_sub_objects = 0; 

            const int n_node_subtrees = sub_node_real_ids.size();
            for (int si = 0; si < n_node_subtrees; ++si) {
               n_sub_objects += 
                  sub_l_stats->
                     getTotalNumberDataObjectsInSubtreeForRealNodeId(sub_node_real_ids[si]);
            }

            l_stats->setTotalNumberDataObjectsInSubtree(ni, n_sub_objects);

         }

      }

   }  // if > 0 levels in tree

}

int MTree::getNumberLevels() const
{
   int num_levels_in_tree = 0;
   if (d_root_node) {
      num_levels_in_tree = d_root_node->getLevelInTree()+1;
   }
   return(num_levels_in_tree);
}

const MTreeLevelStatistic*
MTree::getLevelStatistics(int level_number) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level_number < getNumberLevels());
   assert(static_cast<unsigned int>(level_number) < d_level_statistics.size());
#endif
   return( d_level_statistics[level_number] );
}

void MTree::clearLevelStatistics()
{
   const int level_stat_size = d_level_statistics.size();
   if (level_stat_size > 0) {
      for (int ln = 0; ln < level_stat_size; ++ln) {
         if ( d_level_statistics[ln] ) {
            delete d_level_statistics[ln];
         }
      }
      d_level_statistics.clear();
   }
}

/*
*************************************************************************
*                                                                       *
* Print tree and tree level node statistics.                            *
*                                                                       *
*************************************************************************
*/

void MTree::printAllTreeStatistics(ostream& stream) const
{
   stream << "\nMTree::printAllTreeStatistics()\n";
   stream << "d_db_name = " << d_db_name << endl;
   printAllTreeOperationStatistics(stream);

   const int n_tree_levels = getNumberLevels();
   stream << "\n" << n_tree_levels << " levels in tree" << endl;

   const int level_stat_size = d_level_statistics.size();
   if (level_stat_size == 0) {
      stream << "\n   Level statistics not calculated."
             << " Must call calculateStatistics() first!" << endl;
   }

   if (level_stat_size < n_tree_levels) {
      stream << "\n   Level statistics not consistent with tree."
             << " Must call calculateStatistics() first!" << endl;
   }

   for (int ln = 0; ln < level_stat_size; ++ln) {
      if ( d_level_statistics[ln] ) {
         d_level_statistics[ln]->printClassData(stream);
      }
   }
}

void MTree::printTreeLevelStatistics(int level_number, 
                                     ostream& stream) const
{
   stream << "\nMTree::printTreeLevelStatistics()\n";

   const int level_stat_size = d_level_statistics.size();
   if ( (level_number >= level_stat_size) ||
        !d_level_statistics[level_number] ) {
      stream << "\n   Given level number = " << level_number
             << " is not a valid level number for calculated statistics!" << endl;
   } else {
      d_level_statistics[level_number]->printClassData(stream);
   }
}

/*
*************************************************************************
*                                                                       *
* Print tree operation count statistics.                                *
*                                                                       *
*************************************************************************
*/

void MTree::printAllTreeOperationStatistics(ostream& stream) const
{
   stream << "\nMTree::printAllTreeOperationStatistics()\n";
   stream << "------------------------------------------\n";
   stream << "d_num_range_queries = " << d_num_range_queries << endl;
   stream << "d_num_knn_queries = " << d_num_knn_queries << endl;
   stream << "d_num_inserts = " << d_num_inserts << endl;
   stream << "d_num_deletes = " << d_num_deletes << endl;
   stream << "d_total_distance_comps_in_range_queries = " 
          << d_total_distance_comps_in_range_queries << endl;
   stream << "d_total_distance_comps_in_knn_queries = " 
          << d_total_distance_comps_in_knn_queries << endl;
   stream << "d_total_distance_comps_in_inserts = " 
          << d_total_distance_comps_in_inserts << endl;
   stream << "d_total_distance_comps_in_deletes = " 
          << d_total_distance_comps_in_deletes << endl;
   stream << "d_num_distance_comps_in_last_range_query = " 
          << d_num_distance_comps_in_last_range_query << endl;
   stream << "d_num_distance_comps_in_last_knn_query = " 
          << d_num_distance_comps_in_last_knn_query << endl;
   stream << "d_num_distance_comps_in_last_insert = " 
          << d_num_distance_comps_in_last_insert << endl;
   stream << "d_num_distance_comps_in_last_delete = " 
          << d_num_distance_comps_in_last_delete << endl;
}

/*
*************************************************************************
*                                                                       *
* Private methods to build and clear tree operation count statistics.   * 
*                                                                       *
*************************************************************************
*/

void MTree::incrementDistanceComputeCount()
{
   switch(d_current_distance_type) {

      case DISTANCE_RANGE_QUERY: {
         d_total_distance_comps_in_range_queries++;
         d_num_distance_comps_in_last_range_query++;
         break;
      }

      case DISTANCE_KNN_QUERY: {
         d_total_distance_comps_in_knn_queries++;
         d_num_distance_comps_in_last_knn_query++;
         break;
      }

      case DISTANCE_INSERT: {
         d_total_distance_comps_in_inserts++;
         d_num_distance_comps_in_last_insert++;
         break;
      }

      case DISTANCE_DELETE: {
         d_total_distance_comps_in_deletes++;
         d_num_distance_comps_in_last_delete++;
         break;
      }

      default: {
      }
   }
}

void MTree::clearOperationCountStatistics()
{
  d_current_distance_type = DISTANCE_UNDEFINED;
  d_num_range_queries = 0; 
  d_num_knn_queries = 0; 
  d_num_inserts = 0; 
  d_num_deletes = 0; 
  d_total_distance_comps_in_range_queries = 0; 
  d_total_distance_comps_in_knn_queries = 0; 
  d_total_distance_comps_in_inserts = 0; 
  d_total_distance_comps_in_deletes = 0; 
  d_num_distance_comps_in_last_range_query = 0; 
  d_num_distance_comps_in_last_knn_query = 0; 
  d_num_distance_comps_in_last_insert = 0; 
  d_num_distance_comps_in_last_delete = 0;

}

void MTree::outputStats(std::ostream & outputStream)
{
   //
   // calculate statistics 
   //

   calculateLevelStatistics();

   //
   // print all tree statistics
   //

   printAllTreeStatistics(outputStream);

   //
   // get number of levels
   //

   const int numberLevels = getNumberLevels();

   //
   // output levels
   //

   for (int i = 0; i < numberLevels; ++i) {

      //
      // output level number
      //

      outputStream << "Level " << i << std::endl;

      //
      // get level stats
      //

      const MTreeLevelStatistic * levelStats = getLevelStatistics(i);

      //
      // get number of nodes at this level
      //

      const int numberNodes = levelStats->getNumberNodesOnLevel();

      //
      // iterate over nodes on this level
      //

      outputStream << "Node   Number entries   Number data leaf nodes"
                   << std::endl;

      for (int iNode = 0; iNode < numberNodes; ++iNode) {

         //
         // get node stat
         //

         const MTreeNodeStat & nodeStat = levelStats->getNodeStat(iNode);

         //
         // output 
         //
	    
         outputStream << iNode << " "
                      << nodeStat.getNumberEntries() << " "
                      << nodeStat.getTotalNumberDataObjectsInSubtree() << " "
                      << std::endl;
	  
      }	  

   }
	
   return;
}


#endif





