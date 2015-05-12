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
// File:        MTreeNode.cc
// Package:     MTree database
// 
// 
// 
// Description: Representation of node in an MTree.
//

#ifndef included_MtreeNode_C
#define included_MTreeNode_C

#ifndef included_MTreeNode
#include "MTreeNode.h"
#endif

#ifndef included_MTree
#include "MTree.h"
#endif

#ifndef included_toolbox_Utilities
#include "toolbox/Utilities.h"
#endif

#ifndef included_toolbox_MathUtilities
#include "toolbox/MathUtilities.h"
#endif 

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif
#endif

#ifdef DEBUG_NO_INLINE
#include "MTreeNode.I"
#endif


/*
*************************************************************************
*                                                                       *
* Node contructor.                                                      *
* FIXME: This used to an inline function.                               *
*                                                                       *
*************************************************************************
*/

MTreeNode::MTreeNode(MTree* tree,
                     int max_entries)
:  d_my_tree(tree),
   d_is_root_node(false),
   d_node_id( s_node_instance_counter++ ),
   d_leaf_node_id( MTreeNode::getUndefinedId() ),
   d_level_in_tree(0),
   d_max_entries(max_entries)
{
   d_entries.reserve(max_entries+1);
   d_my_tree->addNodeToLevelCount(d_level_in_tree);
}


/*
*************************************************************************
*                                                                       *
* Accessory function to get and set node id and get shared pointer      *
* to this object.                                                       *
* FIXME: This used to an inline function.                               *
*                                                                       *
*************************************************************************
*/
 
MTreeNodePtr MTreeNode::getMyself() const
{
   if (isRoot()) {
      return( d_my_tree->getRootNode() );
   } else {
#ifdef DEBUG_CHECK_ASSERTIONS
      assert( d_parent_entry );
#endif
      return( d_parent_entry->getSubtreeNode() );
   }
}

/*
*************************************************************************
*                                                                       *
* Accessory routines to check type of node and whether it is defined,   * 
* and also to set node level in tree.                                   *
* FIXME: This used to an inline function.                               *
*                                                                       *
*************************************************************************
*/

void MTreeNode::setLevelInTree(int level) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level >= 0);
#endif
   d_my_tree->removeNodeFromLevelCount(d_level_in_tree);
   d_level_in_tree = level;
   d_my_tree->addNodeToLevelCount(d_level_in_tree);
   if ( !isLeaf() ) {
      d_leaf_node_id = MTreeNode::getUndefinedId();
   }
}

/*
*************************************************************************
*                                                                       *
* File scope function for determining min entry of a vector<double>.    *
*                                                                       *
*************************************************************************
*/

static int findMin(const vector<double>& vec)
{
   int imin = -1;
  
   double min_val = toolbox::MathUtilities<double>::getMax(); 
   const int veclen = vec.size();
   for (int i = 0; i < veclen; ++i) {
      if (vec[i] < min_val) {
         imin = i;
         min_val = vec[i];
      }
   }

   return(imin); 
}

/*
*************************************************************************
*                                                                       *
* Initialization of and access to static data members.                  *
*                                                                       *
*************************************************************************
*/

int MTreeNode::s_undefined_node_id = -1;
int MTreeNode::s_node_instance_counter = 0;

int MTreeNode::getUndefinedId()
{
   return( s_undefined_node_id );
}

/*
*************************************************************************
*                                                                       *
* Dtor for MTreeNode.                                                   *
*                                                                       *
*************************************************************************
*/

MTreeNode::~MTreeNode()
{
   d_my_tree = (MTree*)NULL; 
   d_parent_entry.reset();
   const int nentries = d_entries.size();
   for (int ie = 0; ie < nentries; ++ie) {
      d_entries[ie].reset();
   }
   d_entries.clear();
}

/*
*************************************************************************
*                                                                       *
* Reset radii for all routing entries along path from this node to      *
* the root node.                                                        *
*                                                                       *
*************************************************************************
*/

void MTreeNode::resetRadiiUpToRoot()
{
   MTreeNodePtr node2reset( getMyself() ); 
   while ( !node2reset->isRoot() ) {
      node2reset->resetRadius();
      node2reset = node2reset->getParentNode();
   }
}

/*
*************************************************************************
*                                                                       *
* Reset radius for routing entry associated with this node.             *
*                                                                       *
* Method asssumes each node entry has correct radius and                *
* distance-to-parent value set.                                         *
*                                                                       *
*************************************************************************
*/

void MTreeNode::resetRadius()
{
   if ( !isRoot() ) {
      double max_radius = 0.0;
      const int num_entries = getNumberEntries();
      if ( num_entries > 0 ) {
         for (int ie = 0; ie < num_entries; ++ie) {
            MTreeEntryPtr entry(d_entries[ie]);
            max_radius = 
               toolbox::MathUtilities<double>::Max( max_radius,
                  entry->getDistanceToParent() + entry->getRadius() );
         }
      }
      d_parent_entry->setRadius(max_radius);
   }
}

/*
*************************************************************************
*                                                                       *
* Search this node for entry with best subtree for object insertion.    *
* This routine attampts to minimize the number of distance computations *
* be applying optimizations based on the triangle inequality.           *
*                                                                       *
*************************************************************************
*/

MTreeNodePtr 
MTreeNode::searchBestSubtreeForInsert(MTreeEntryPtr new_entry) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(new_entry.get());
#endif

   MTreeEntryPtr best_entry;

   double best_distance = toolbox::MathUtilities<double>::getMax();
   double best_penalty = toolbox::MathUtilities<double>::getMax();

   double new_dist2par = new_entry->getDistanceToParent();
   double new_radius = new_entry->getRadius();

   const int num_entries = getNumberEntries();

   vector<double> distances(getNumberEntries());

   for (int ie = 0; ie < num_entries; ++ie) {

      distances[ie] = toolbox::MathUtilities<double>::getMax();
      
      MTreeEntryPtr test_entry(d_entries[ie]);
      double test_dist2par = test_entry->getDistanceToParent();
      double test_radius = test_entry->getRadius();

      bool done_with_test = false;
      double tdistance = 0.0; 
      double tpenalty = toolbox::MathUtilities<double>::getMax();

      /*
       * If distance to parent for new entry has already been set, we 
       * attempt to use it to avoid distance computations by applying
       * the triangle inequality. 
       */

      if (new_dist2par > 0.0) {

         double rad_sum = test_radius + new_radius;
         double dist_penalty = 
            toolbox::MathUtilities<double>::Abs(new_dist2par - test_dist2par) - rad_sum;

         if (dist_penalty < 0.0) {
            /*
             * If distance penalty is negative, the hypothesis of Lemma 2 
             * doesn't apply.  Thus, we cannot conclude that the regions 
             * defined by the test entry and the new entry do not intersect.
             */
            dist_penalty += rad_sum - MetricSpacePoint::getMaxDistance(); 
         }

         if ( (best_penalty < toolbox::MathUtilities<double>::getMax()) &&
              (dist_penalty >= best_penalty) ) {
            /*
             * If the best penalty has been set (i.e., we have found a possible 
             * entry) and if the current distance penalty is greater than the 
             * best penalty, we are done with this test entry and we do not 
             * need to compute distance between the test entry and the new entry.
             */
            done_with_test = true;
         }

      }

      if (!done_with_test) {
        
         distances[ie] = test_entry->computeDistanceTo(new_entry);
         d_my_tree->incrementDistanceComputeCount(); 

         tdistance = distances[ie] + new_radius;
         tpenalty = tdistance - test_radius;

         if (tpenalty < 0.0) {
            /*
             * If tpenalty < 0.0, then the region defined by the new entry
             * fits entirely within the region defined by the test entry.
             * We set tpenalty so that it will be assigned to the nearest
             * such entry.
             */     
            tpenalty = tdistance - MetricSpacePoint::getMaxDistance();
         }

      }

      if ( !best_entry || (tpenalty < best_penalty) ) {
         best_entry = test_entry;
         best_distance = tdistance;
         best_penalty = tpenalty;
      }
   }

   /*
    * New entry will be inserted into region defined by best entry. 
    * Set distance of new entry to distance between it and best entry.
    * Then,
    */
   new_dist2par = distances[best_entry->getPositionInNode()];
   if ( !(new_dist2par < toolbox::MathUtilities<double>::getMax()) ) {
      new_dist2par = best_entry->computeDistanceTo(new_entry);
      d_my_tree->incrementDistanceComputeCount(); 
   }

   new_entry->setDistanceToParent(new_dist2par);

   return( best_entry->getSubtreeNode() );

}

/*
*************************************************************************
*                                                                       *
* Find proper position for entry in node entry array and insert.        *
*                                                                       *
* Preconditions:                                                        *
*    This method assumes -- entries in the node are properly ordered    *
*                           with respect to increasing distance from    *
*                           the parent entry of this node.              *
*                        -- the entry has a correct distance-to-parent  *
*                           value set with respect to the parent entry  *
*                           of the node.                                *
*                        -- the entry vector for the node is not        *
*                           full to capacity.                           *   
*                                                                       *
* Postconditions:                                                       *
*    - the entry is inserted into the entry vector for the node         *
*      at the proper position to preserve ordering with respect to      *
*      increasing distance from the parent entry of the node.           *
*    - the node of the entry is set to the node.                        *
*    - the position-in-node value for entry is set to inserted position.*
*                                                                       *
*************************************************************************
*/

void MTreeNode::insertEntry(MTreeEntryPtr entry,
                            MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(entry.get());
   assert(node.get());
   const int entries_cap = node.get()->d_entries.capacity();
   assert(node->getNumberEntries() < entries_cap);
#endif

   if ( node->isRoot() ) {

      entry->setNode( node );
      entry->setPositionInNode( node->getNumberEntries() );
      node.get()->d_entries.push_back(entry);

   } else {

      int position = 0;
      bool found = false;

      while ( !found && (position < node.get()->getNumberEntries()) ) {
         if ( node.get()->d_entries[position]->getDistanceToParent() >
              entry->getDistanceToParent() ) {
            found = true;
         } else {
            position++;
         }
      }

      MTreeNode::insertEntryAtPosition(entry, position, node);

   }

}

/*
*************************************************************************
*                                                                       *
* Insert entry at given position in node if position is correct.        *
* Otherwise find proper position for insertion.                         *
*                                                                       *
* Preconditions:                                                        *
*    This method assumes -- entries in the node are properly ordered    *
*                           with respect to increasing distance from    *
*                           the parent entry of the node.               *
*                        -- the entry has a correct distance-to-parent  *
*                           value set with respect to the parent entry  *
*                           of the node.                                *
*                        -- the entry vector for the node is not        *
*                           full to capacity.                           *   
*                                                                       *
* Postconditions:                                                       *
*    - the entry is inserted into the entry vector for the node         *
*      at the proper position to preserve ordering with respect to      *
*      increasing distance from the parent entry of the node.           *
*    - the node of the entry is set to the node.                        *
*    - the position-in-node value for entry is set to inserted position.*
*                                                                       *
*************************************************************************
*/

void MTreeNode::insertEntryAtPosition(MTreeEntryPtr entry,
                                      int position,
                                      MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(entry.get());
   assert(node.get());
   const int entries_cap = node.get()->d_entries.capacity();
   assert(node->getNumberEntries() < entries_cap);
#endif

   const int num_entries = node->getNumberEntries();

   bool position_correct =
        ( (position >= 0) && (position <= num_entries) );
   if ( position_correct && (num_entries > 0) ) {
      if (position > 0) {
         position_correct =
            ( node.get()->d_entries[position-1]->getDistanceToParent() <=
              entry->getDistanceToParent() );
      }
      if ( position_correct && (position < num_entries) ) {
         position_correct =
            ( entry->getDistanceToParent() <=
              node.get()->d_entries[position]->getDistanceToParent() );
      }
   } else {
      position_correct = (position == 0);
   }

   if (position_correct) {
      const int hi_move = node->getNumberEntries() - 1;
      const int lo_move = position;

      MTreeEntryPtr dummy;
      node.get()->d_entries.push_back(dummy);

      for (int ie = hi_move; ie >= lo_move; ie--) {
         node.get()->d_entries[ie+1] = node.get()->d_entries[ie];
         node.get()->d_entries[ie+1]->setPositionInNode(ie+1);
      }

      entry->setNode( node );
      entry->setPositionInNode(position);
      node.get()->d_entries[position] = entry;

   } else {
      MTreeNode::insertEntry(entry, node);
   }

}

/*
*************************************************************************
*                                                                       *
* Delete given entry from the node.                                     *
*                                                                       *
*************************************************************************
*/

void MTreeNode::deleteEntry(MTreeEntryPtr entry,
                            MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(entry.get());
   assert(node.get());
   assert(entry->getPositionInNode() < node->getNumberEntries());
#endif
   const int position = entry->getPositionInNode();

   if ( (entry->getNode().get() == node.get()) &&
        (entry.get() == node.get()->d_entries[position].get()) ) {

      node.get()->d_entries[position].reset();

      const int lo_move = position + 1;
      const int hi_move = node->getNumberEntries() - 1;

      for (int ie = lo_move; ie <= hi_move; ie++) {
         node.get()->d_entries[ie-1] = node.get()->d_entries[ie];
         node.get()->d_entries[ie-1]->setPositionInNode(ie-1);
      }

      node.get()->d_entries.pop_back();
   }

}

/*
*************************************************************************
*                                                                       *
* Simple transfer of entries from this node to argument node.           *
*                                                                       *
* This rotuine assumes argument node contains no entries.               *
*                                                                       *
* This routine assumes nothing about the order of the entries in this   *
* node nor whether their distance-to-parent values have been properly   *
* set.  On exit, this node has no entries, and the other node owns all  *
* entries.  The order of entries is preserved in the other node.        *
*                                                                       *
*************************************************************************
*/

void MTreeNode::transferEntries(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
   assert(node->getNumberEntries() == 0);
#endif
   d_entries.swap(node->d_entries);

   const int num_entries = node->getNumberEntries();
   for (int ie = 0; ie < num_entries; ++ie) {
      node->getEntry(ie)->setNode(node);
   }

}

/*
*************************************************************************
*                                                                       *
* Split this node and return the new node resulting from the split.     *
*                                                                       *
* If this node is the root of the tree, the method does nothing         *
* and returns a null pointer.  Otherwise, if this node is not           *
* over-full and the "special_delete_split" flag is false, the           *
* method does nothing and returns a null pointer.                       *
*                                                                       *
* Preconditions:                                                        *
*    This method assumes -- this node has at least two entries.         *
*                        -- all entries have correct radius and         *
*                           distance-to-parent information set.         *
*                                                                       *
* Postconditions:                                                       *
*    - this node remains a child of its parent node.                    *
*    - a new node is created and is a child of the parent of this node. *
*    - parent entries for this node and new node are promoted to the    *
*      parent node of this node.                                        *
*    - the entries originally in this node are partitioned into two     *
*      disjoint sets and assigned to this node and the new node.        *
*    - all distance-to-parent are set properly for all entries.         *
*    - the distance-to-parent values are also reset for the entries     *
*      promoted to the parent node.                                     *
*    - the radii of this node and the new node (i.e., promoted entry    *
*      radii) are reset to properly represent the regions defined       *
*      by the entries in each node.                                     *
*    - the radii of the parent node is reset, if parent is not root.    *
*                                                                       *
*************************************************************************
*/

MTreeNodePtr MTreeNode::split(bool use_root_promotion,
                              bool special_delete_split) 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(getNumberEntries() >= 2);
#endif

   if ( !isRoot() && (isOverfull() || special_delete_split) ) {

      MTreeNodePtr parent_node( getParentNode() ); 

      MTreeNodePtr node2( new MTreeNode( d_my_tree, d_max_entries) );

      /*
       * Choose two entries in this node to promote to routing 
       * entries in parent node.
       */ 
      
      MTreeEntryPtr tmp_entry1;
      MTreeEntryPtr tmp_entry2;
      bool new_parent_entry = promoteRoutingEntries(tmp_entry1, 
                                                    tmp_entry2, 
                                                    use_root_promotion,
                                                    special_delete_split);

      /*
       * Make copies of promoted entries for insertion into parent node.
       */ 
      
      MTreeEntryPtr entry1( new MTreeEntry(tmp_entry1->getKey()) );
      MTreeEntryPtr entry2( new MTreeEntry(tmp_entry2->getKey()) );

      /*
       * Partition entries between this node and node2.  Note that
       * entry1 will be assigned to this node and entry2 will be 
       * assigned to node2.
       */ 
      
      partitionEntries(tmp_entry1,
                       new_parent_entry, 
                       tmp_entry2,
                       node2); 

      /*
       * Insert promoted entries into parent node.
       */ 
      
      if (new_parent_entry) {
         entry1->setSubtreeNode( getMyself() );
         MTreeNode::deleteEntry(d_parent_entry, parent_node);
         setParentEntry(entry1);
         if (!parent_node->isRoot()) {
            entry1->setDistanceToParent( 
               entry1->computeDistanceTo( parent_node->getParentEntry() ) );
            d_my_tree->incrementDistanceComputeCount(); 
         } 
      }
      MTreeNode::insertEntry(entry1, parent_node);

      entry2->setSubtreeNode(node2);
      node2->setParentEntry(entry2);
      if (!parent_node->isRoot()) {
         entry2->setDistanceToParent(
               entry2->computeDistanceTo( parent_node->getParentEntry() ) );
         d_my_tree->incrementDistanceComputeCount(); 
      }
      MTreeNode::insertEntry(entry2, parent_node);

      /*
       * Reset radius for split node and node2 (and parent, if needed)
       * to preserve semantics of covering radii.
       */ 
      
      resetRadius();
      node2->resetRadius();

      if (!parent_node->isRoot()) {
         parent_node->resetRadius();
      }

      return( node2 );

   } else {
      MTreeNodePtr dummy;
      return( dummy );
   }

}

/*
*************************************************************************
*                                                                       *
* Identify two entries in the node for promotion based on the           *
* algorithm choice specified by the enum type MTreeNodePromotionMethod  *
* (see header file).  The boolean return value is true if neither of    *
* the promoted entries is the current node parent entry; otherwise it   *
* is false.                                                             *
*                                                                       *
* Preconditions:                                                        *
*    This method assumes -- this node has at least two entries.         *
*                        -- each entry has correct radius value set.    *
*                        -- each entry has correct distance-to-parent   *
*                           value set, except when root promotion       *
*                           method is used.                             *
*                                                                       *
* Postconditions:                                                       *
*    - entry1, entry2 will be set to the promoted entries.              *
*    - entry1 will be the promoted entry corresponding to this node.    *
*    - boolean return value will be true if entry1 IS NOT the current   *
*      parent entry of this node and false otherwise.                   *
*    - distance-to-parent and radius of each entry is unchanged.        *
*                                                                       *
*************************************************************************
*/

bool MTreeNode::promoteRoutingEntries(MTreeEntryPtr& entry1, 
                                      MTreeEntryPtr& entry2,
                                      bool use_root_promotion,
                                      bool special_delete_split) const
{

   bool new_parent_entry = false;

   if ( !isRoot() ) {

      if ( getNumberEntries() < 2 ) {

         TBOX_ERROR("MTreeNode::promoteRoutingEntries() error..."
                    << " Node # " << getNodeId()
                    << " has less than 2 entries!" << endl);

      } else {

         if ( isOverfull() || special_delete_split ) {
   
            if ( use_root_promotion ) {

               new_parent_entry = 
                  promoteTwo(entry1, 
                             entry2, 
                             d_my_tree->getRootNodePromotionMethod());

            } else {

               MTreeNodePromotionMethod promotion_method = 
                  d_my_tree->getNodePromotionMethod();

               if ( (promotion_method == CONFIRMED_RANDOM_PROMOTION) || 
                    (promotion_method == CONFIRMED_MAX_LB_DISTANCE_PROMOTION) ) {

                  new_parent_entry = false;
                  entry1 = d_parent_entry;
         
                  promoteOne(entry2, 
                             promotion_method);

               } else {
               
                  new_parent_entry = promoteTwo(entry1, 
                                                entry2, 
                                                promotion_method);

               }

            } // else do not use root promotion

         }  // if node is over-full

      } // else node has >= 2 entries

   } // if not root node

   return(new_parent_entry);

}

/*
*************************************************************************
*                                                                       *
* Identify two entries in node for promotion using the given method.    * 
*                                                                       *
* Preconditions:                                                        *
*    This method assumes -- this node has at least two entries.         *
*                        -- all entries have correct radius and         *
*                           distance-to-parent information set.         *
*                                                                       *
* Postconditions:                                                       *
*    - entry1, entry2 will be set to the promoted entries.              *
*    - entry1 will be the promoted entry corresponding to this node.    *
*    - boolean return value will be true if entry1 IS NOT the current   *
*      parent entry of this node and false otherwise.                   *
*    - distance-to-parent and radius of each entry is unchanged.        *
*                                                                       *
*************************************************************************
*/

bool MTreeNode::promoteTwo(MTreeEntryPtr& entry1,
                           MTreeEntryPtr& entry2,
                           MTreeNodePromotionMethod method) const 
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(getNumberEntries() >= 2);
#endif

   bool new_parent_entry = false;

   const int num_entries = getNumberEntries(); 

   int indx_entry1 = -1;
   int indx_entry2 = -1;

   switch(method) {

      case RANDOM_PROMOTION: {

         indx_entry1 = toolbox::MathUtilities<int>::Rand(0, num_entries);
         indx_entry2 = indx_entry1;
         while (indx_entry1 == indx_entry2) {
            indx_entry2 = 
               toolbox::MathUtilities<int>::Rand(0, num_entries);
         }

         break;
      }

      case MAX_SPREAD_DISTANCE_PROMOTION: {

         indx_entry1 = num_entries-1;
         indx_entry2 = num_entries-2;

         break;
      }

      case MIN_RADII_PROMOTION: {

         int ie;

         indx_entry1 = 0;
         double radius1 = d_entries[indx_entry1]->getRadius();

         indx_entry2 = 1;
         double radius2 = d_entries[indx_entry2]->getRadius();

         for (ie = 1; ie < num_entries; ++ie) {
            double eradius = d_entries[ie]->getRadius();
            if ( (eradius < radius1) ) {
               indx_entry1 = ie;
               radius1 = eradius;
            }
         }

         if ( indx_entry1 == 1 ) {
            indx_entry2 = 0;
            radius2 = d_entries[indx_entry2]->getRadius(); 
         }

         for (ie = 0; ie < num_entries; ++ie) {
            if ( ie != indx_entry1 ) {
               double eradius = d_entries[ie]->getRadius();
               if ( (eradius < radius1) ) {
                  indx_entry2 = ie;
                  radius2 = eradius;
               }
            }
         }

         break;
      }

      case MIN_OVERLAP_PROMOTION: {

         double min_overlap = -MetricSpacePoint::getMaxDistance();
 
         for (int ie1 = 0; ie1 < num_entries; ++ie1) {
            MTreeEntryPtr entry1(d_entries[ie1]);
            for (int ie2 = num_entries-1; ie2 > ie1; ie2--) {
               MTreeEntryPtr entry2(d_entries[ie2]);
               double overlap_val = 
                  entry2->computeDistanceTo(entry1) - 
                     ( entry1->getRadius() + entry2->getRadius() );
               d_my_tree->incrementDistanceComputeCount(); 
               if ( overlap_val > min_overlap ) {
                  min_overlap = overlap_val;
                  indx_entry1 = ie1;
                  indx_entry2 = ie2; 
               }
            }
         }
        
         break;
      } 

      default: {
         TBOX_ERROR("MTreeNode::promoteTwo() error..."
                    << " Promotion method " << method 
                    << " not recognized" << endl);
      }

   }

   if (indx_entry1 == indx_entry2) {
      TBOX_ERROR("MTreeNode::promoteTwo() error..."
                 << "\n  Promoted entries are the same, entry # " 
                 << indx_entry1 << endl);
   } else {
      int indx_tmp1 = toolbox::MathUtilities<int>::Min(indx_entry1,
                                                       indx_entry2);
      int indx_tmp2 = toolbox::MathUtilities<int>::Max(indx_entry1,
                                                       indx_entry2);
      entry1 = d_entries[indx_tmp1];
      entry2 = d_entries[indx_tmp2];
   }

   if ( entry2.get() == d_parent_entry.get() ) {
       MTreeEntryPtr etmp(entry1);
       entry1 = entry2;
       entry2 = etmp;
   }

   if ( entry1.get() != d_parent_entry.get() ) {
      new_parent_entry = true;
   }

   return(new_parent_entry);

}

/*
*************************************************************************
*                                                                       *
* Identify one entry in the node for promotion using the given method.  *
* The promoted entry will be different than the current node parent.    *
*                                                                       *
* Preconditions:                                                        *
*    This method assumes -- this node has at least two entries.         *
*                        -- all entries have correct radius and         *
*                           distance-to-parent information set.         *
*                                                                       *
* Postconditions:                                                       *
*    - entry2 will be set to the promoted entry.                        *
*    - entry2 will be different that the current parent entry           *
*      of this node.                                                    *
*    - distance-to-parent and radius of each entry is unchanged.        *
*                                                                       *
*************************************************************************
*/

void MTreeNode::promoteOne(MTreeEntryPtr& entry2,
                           MTreeNodePromotionMethod method) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(getNumberEntries() >= 2);
#endif

   const int num_entries = getNumberEntries(); 

   switch(method) {

      case CONFIRMED_RANDOM_PROMOTION: {

         int indx_entry2 = 
            toolbox::MathUtilities<int>::Rand(0, num_entries);
         entry2 = d_entries[indx_entry2];

         while ( entry2.get() == d_parent_entry.get() ) {
            indx_entry2 =
               toolbox::MathUtilities<int>::Rand(0, num_entries);
         }
         entry2 = d_entries[indx_entry2];

         break;
      }

      case CONFIRMED_MAX_LB_DISTANCE_PROMOTION: {

         entry2 = d_entries[num_entries-1]; 

         break;
      }

      default: {
         TBOX_ERROR("MTreeNode::promoteOne() error..."
                    << " Promotion method " << method 
                    << " not recognized" << endl);
      }

   }

}

/*
*************************************************************************
*                                                                       *
* Partition entries between this node and node2.                        *
*                                                                       *
* Preconditions:                                                        *
*    This method assumes -- entry1 and entry2 are proper entries to     *
*                           serve as parent entries of this node and    *
*                           node2, respectively.                        *
*                        -- this node owns all entries and node2        *
*                           has zero entries.                           *
*                        -- distance-to-parent values are correctly     *
*                           set for all entries.                        *
*                        -- if the boolean argument is false, this      *
*                           node does not have a new parent entry       *
*                           and so only distances to node2 parent       *
*                           entry must be computed.                     *
*                                                                       *
* Postconditions:                                                       *
*    - entry1 will be assigned to this node and end entry2 will         *
*      be assigned to node2.                                            *
*    - this node and node2 will own disjoint sets of entries whose      *
*      union is the set of entries in this node upon entry to method.   *
*    - distance-to-parent values will be correctly set for all          *
*      entries with respect to entry1 and entry2.                       *
*    - entries in this node and node2 are ordered correctly with        *
*      respect to increasing distance from their parent entries.        *
*                                                                       *
*************************************************************************
*/

void MTreeNode::partitionEntries(MTreeEntryPtr entry1,
                                 bool new_parent_entry,
                                 MTreeEntryPtr entry2,
                                 MTreeNodePtr& node2)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(entry1.get());
   assert(entry2.get());
   assert(node2.get());
   assert(getNumberEntries() >= 2);
   assert(node2->getNumberEntries() == 0);
#endif

   const int pos_entry1 = entry1->getPositionInNode();
   const int pos_entry2 = entry2->getPositionInNode();

#ifdef DEBUG_CHECK_ASSERTIONS
   assert(0 <= pos_entry1 && pos_entry1 < getNumberEntries());
   assert(0 <= pos_entry2 && pos_entry2 < getNumberEntries());
   assert(pos_entry1 != pos_entry2);
#endif

   if ( !isRoot() ) {

      /*
       * Both partition methods require distances from each
       * entry to node parent entries.  We compute only those
       * distances that must be computed and store all distances.
       */

      int ie;
      const int num_entries = getNumberEntries();

      vector<int> assignment(num_entries, -1);
 
      vector<double> distance1(num_entries,
                               toolbox::MathUtilities<double>::getMax());
      vector<double> distance2(num_entries,
                               toolbox::MathUtilities<double>::getMax());

      if (new_parent_entry) {
         for (ie = 0; ie < num_entries; ++ie) {
            distance1[ie] = d_entries[ie]->computeDistanceTo(entry1);
            d_my_tree->incrementDistanceComputeCount();
         }
      } else {
         for (ie = 0; ie < num_entries; ++ie) {
            distance1[ie] = d_entries[ie]->getDistanceToParent();
         }
      }

      for (ie = 0; ie < num_entries; ++ie) {
         distance2[ie] = d_entries[ie]->computeDistanceTo(entry2);
         d_my_tree->incrementDistanceComputeCount();
      }

      /*
       * Assign each of the given promoted entries to the proper node.
       */

      vector<double> tdist1(distance1);
      vector<double> tdist2(distance2);
 
      assignment[pos_entry1] = 1; 
      tdist1[pos_entry1] = toolbox::MathUtilities<double>::getMax();  

      assignment[pos_entry2] = 2; 
      tdist2[pos_entry2] = toolbox::MathUtilities<double>::getMax();  

      int num_entries1 = 1;
      int num_entries2 = 1;

      MTreeNodePartitionMethod partition_method =
         d_my_tree->getNodePartitionMethod();

      switch(partition_method) {

         case HYPERPLANE_PARTITION: {

            /*
             * 1) Set minimum number of node entries.
             * 2) Assign each entry to the node with the closest 
             *    parent entry. 
             * 3) If either node has fewer than minimum number
             *    of entries, move entries to the under-full node
             *    that are closest to the parent entry of the 
             *    under-full node.
             */

            double min_node_utilization =
               d_my_tree->getMinNodeUtilization();

            const int min_split_entries = 
               toolbox::MathUtilities<int>::Max( 1, 
                  static_cast<int>( getNumberEntries()*min_node_utilization) );

            /*
             * Assign rest of entries based on closest entry.  
             */

            for (ie = 0; ie < num_entries; ++ie) {
               if ( (ie != pos_entry1) && (ie != pos_entry2) ) {
                  if ( distance1[ie] < distance2[ie] ) {
                     assignment[ie] = 1;
                     num_entries1++;
                  } else {
                     assignment[ie] = 2;
                     num_entries2++;
                  }
               }
            }

            /*
             * Balance entries as needed.
             */

            if (num_entries1 < min_split_entries) {
               while (num_entries1 < min_split_entries) {
                  int indx_entry = findMin(tdist1);
                  if (assignment[indx_entry] != 1) {
                     assignment[indx_entry] = 1;
                     num_entries1++;
                     num_entries2--;
                  }
                  tdist1[indx_entry] = toolbox::MathUtilities<double>::getMax();
               }
            }

            /*
             * Assign rest of entries based on closest entry.  
             */

            if (num_entries2 < min_split_entries) {
               while (num_entries2 < min_split_entries) {
                  int indx_entry = findMin(tdist2);
                  if (assignment[indx_entry] != 2) {
                     assignment[indx_entry] = 2;
                     num_entries2++;
                     num_entries1--;
                  }
                  tdist2[indx_entry] = toolbox::MathUtilities<double>::getMax();
               }
            }

            break;
         }  // case HYPERPLANE_PARTITION

         case BALANCED_PARTITION: {

            /*
             * Alternately, assign unassigned entries to each node based on
             * the closest parent entry.
             */

            ie = 0;
            while (ie < num_entries) {

               if ( (ie != pos_entry1) && (ie != pos_entry2) ) {

                  int indx_entry = -1;
                  if ( ie % 2 ) {  // assign entry to node2
                     indx_entry = findMin(tdist2);
                     assignment[indx_entry] = 2;
                  } else { // assign entry to this node
                     indx_entry = findMin(tdist1);
                     assignment[indx_entry] = 1;
                  }
                  tdist1[indx_entry] = 
                     toolbox::MathUtilities<double>::getMax();
                  tdist2[indx_entry] = 
                     toolbox::MathUtilities<double>::getMax();

               }

               ie++;

            }

            break;
         }  // case BALANCED_PARTITION

         default: {
            TBOX_ERROR("MTreeNode::partitionEntries() error..."
                       << " Partition method " << partition_method
                       << " not recognized" << endl);
         }

      } // switch(partition_method)

      vector<MTreeEntryPtr> tmp_entries(d_entries);
      for (ie = 0; ie < num_entries; ++ie) {
         d_entries.pop_back();
      }

      for (ie = 0; ie < num_entries; ++ie) {
         if (assignment[ie] == 1) {
            tmp_entries[ie]->setDistanceToParent(distance1[ie]);
            MTreeNode::insertEntry( tmp_entries[ie], getMyself() );
         } else {
            tmp_entries[ie]->setDistanceToParent(distance2[ie]);
            MTreeNode::insertEntry( tmp_entries[ie], node2 );
         }
      }


   }  // if not root node

}

/*
*************************************************************************
*                                                                       *
* Check consistency of this node and all others along its path to root. *
*                                                                       *
*************************************************************************
*/

bool MTreeNode::checkConsistencyUpToRoot(ostream& stream) const
{
   bool nodes_are_consistent = true;

   MTreeNodePtr node2check( getMyself() );
   nodes_are_consistent &= node2check->checkConsistency(stream);
   while ( !node2check->isRoot() ) {
      node2check = node2check->getParentNode();
      nodes_are_consistent &= node2check->checkConsistency(stream);
   }

   return( nodes_are_consistent );
}

/*
*************************************************************************
*                                                                       *
* Check consistency of this node.                                       *
*                                                                       *
*************************************************************************
*/
 
bool MTreeNode::checkConsistency(ostream& stream) const 
{
   bool node_is_consistent = true;

   if ( d_node_id <= MTreeNode::getUndefinedId() ) {
      stream << "MTREE NODE ERROR: Bad node id!" << endl;
      node_is_consistent = false;
   }

   if ( isLeaf() && d_level_in_tree != 0 ) {
      stream << "MTREE NODE ERROR: Leaf node has non-zero level in tree!" << endl;
      node_is_consistent = false;
   }

   if ( isRoot() && d_parent_entry.get() ) {
      stream << "MTREE NODE ERROR: Root node has non-null parent entry!" << endl;
      node_is_consistent = false;
   }

   if ( getNumberEntries() == 0 ) {
      stream << "MTREE NODE ERROR: Node has zero entries!" << endl;
      node_is_consistent = false;
   }

   if ( getNumberEntries() > d_max_entries ) {
      stream << "MTREE NODE ERROR: Node has more than max entries!" << endl;
      node_is_consistent = false;
   }

   int ie;

   bool entries_are_consistent = true;
   for (ie = 0; ie < getNumberEntries(); ++ie) {
      if ( d_entries[ie]->getPositionInNode() != ie ) {
         stream << "MTREE NODE ERROR: Entry " << ie 
                << " has position set to " 
                << d_entries[ie]->getPositionInNode() << endl;
         entries_are_consistent = false;
      }     
      entries_are_consistent &= 
         d_entries[ie]->checkConsistency(getMyself(),
                                         stream); 
   }
   if ( !entries_are_consistent ) {
      stream << "MTREE NODE ERROR: Node has bad entries!" << endl;
      node_is_consistent = false;
   }

   for (ie = 0; ie < getNumberEntries()-1; ++ie) {
      if ( d_entries[ie]->getDistanceToParent() >
           d_entries[ie+1]->getDistanceToParent() ) {
         stream << "MTREE NODE ERROR: Entry " << ie << " out of order!" << endl;
         stream << " entry[" << ie << "] distance = "
                << d_entries[ie]->getDistanceToParent()
                << " : entry[" << ie+1 << "] distance = "
                << d_entries[ie+1]->getDistanceToParent() << endl;
         node_is_consistent = false;
      }
   }

   if ( !node_is_consistent ) {
      printClassData(stream);
   }
 
   return( node_is_consistent );
}

/*
*************************************************************************
*                                                                       *
* Print node summary data to given output stream.                       *
*                                                                       *
*************************************************************************
*/

void MTreeNode::printSummary(ostream& stream) const
{
   stream << "   Node # " << d_node_id << " : " << endl;
   if ( isRoot() ) {
      stream << "      Root node" << endl;
   } 
   if ( isLeaf() ) {
      stream << "      Leaf node" << endl;
      stream << "      Leaf Node # " << d_leaf_node_id << " : " << endl;
   }
   stream << "      d_level_in_tree = " << d_level_in_tree << endl;
   stream << "      number of entries = " << getNumberEntries() << endl;
   if ( isLeaf() ) {
      for (int ie = 0; ie < getNumberEntries(); ++ie) {
          stream << "entry # " << ie << ": " << endl;
          stream << "     point = ";
                 getEntry(ie)->getPoint()->print(stream);    
          stream << "     radius = " << getEntry(ie)->getRadius() << endl;    
          stream << "     owns object " 
                 << getEntry(ie)->getDataObjectId() << endl;
      }
   } else {
      for (int ie = 0; ie < getNumberEntries(); ++ie) {
          stream << "entry # " << ie << ": " << endl;
          stream << "     point = ";
                 getEntry(ie)->getPoint()->print(stream);    
          stream << "     radius = " << getEntry(ie)->getRadius() << endl;    
          stream << "     owns subtree node " 
                 << getEntry(ie)->getSubtreeNode()->getNodeId() << endl;
      }
   }
   
}

/*
*************************************************************************
*                                                                       *
* Print all node data to given output stream.                           *
*                                                                       *
*************************************************************************
*/

void MTreeNode::printClassData(ostream& stream) const
{
   stream << "MTreeNode::printClassData()\n";
   stream << "--------------------------------------\n";
   stream << "this ptr = " << (MTreeNode*)this << endl;
   stream << "d_my_tree = " << (MTree*)d_my_tree << endl;
   stream << "isDefined() = " << isDefined() << endl;
   stream << "d_is_root_node = " << d_is_root_node << endl;
   stream << "isLeaf() = " << isLeaf() << endl;
   stream << "d_node_id = " << d_node_id << endl;
   stream << "d_level_in_tree = " << d_level_in_tree << endl;
   stream << "d_parent_entry = " << (MTreeEntry*)d_parent_entry.get() << endl;
   stream << "d_max_entries = " << d_max_entries << endl;
   stream << "number of entries = " << getNumberEntries() << endl;

   stream << "\nPrinting entries in node..." << endl;

   for (int ie = 0; ie < getNumberEntries(); ++ie) {
       stream << "entry # " << ie << endl;
       getEntry(ie)->printClassData(stream);
   }
   stream << "\n" << endl;
}


#endif




