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
// File:        MTreeEntry.cc
// Package:     MTree database
// 
// 
// 
// Description: Representation of data entry or routing entry in MTree
//

#ifndef included_MtreeEntry_C
#define included_MTreeEntry_C

#include "MTreeEntry.h"

#include "toolbox/base/MathUtilities.h"
#include "MTreeNode.h"

#ifdef DEBUG_NO_INLINE
#include "MTreeEntry.I"
#endif

#include <assert.h>

/*
*************************************************************************
*                                                                       *
* Ctor for MTreeEntry initializes key and set entry type and position   *
* in node to be undefined.                                              *
*                                                                       *
*************************************************************************
*/
 
MTreeEntry::MTreeEntry(const MTreeKey& key)
: d_entry_type(UNDEFINED_ENTRY),
  d_key(key),
  d_my_node_id(MTreeNode::getUndefinedId()),
  d_my_position_in_node(-1),
  d_subtree_node_id(MTreeNode::getUndefinedId()),
  d_data_object_id(DBObject::getUndefinedId())
{
}

/*
*************************************************************************
*                                                                       *
* Dtor for entry.                                                       *
*                                                                       *
*************************************************************************
*/

MTreeEntry::~MTreeEntry()
{
   d_my_node.reset();
   d_subtree_node.reset();
}

/*
*************************************************************************
*                                                                       *
* Accessory routines to set node owning entry, position in node,        *
* subtree node, and object of entry.                                    *
*                                                                       *
*************************************************************************
*/
 
void MTreeEntry::setNode(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node.get());
#endif
   d_my_node = node;
   d_my_node_id = node->getNodeId();
   if (d_my_node->isRoot()) {
      setDistanceToParent( MTreeKey::getUndefinedDistanceToParent() );
   }
}

void MTreeEntry::setPositionInNode(int pos)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_my_node.get());
   assert(0 <= pos && pos <= d_my_node->getMaxEntries());
#endif
   d_my_position_in_node = pos;
}

void MTreeEntry::setSubtreeNode(MTreeNodePtr subtree)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(subtree.get());
#endif
   d_entry_type = ROUTING_ENTRY;
   d_subtree_node = subtree;
   d_subtree_node_id = subtree->getNodeId();
   d_data_object_id = DBObject::getUndefinedId();
}

void MTreeEntry::setDataObjectId(int object_id)
{
   d_entry_type = DATA_ENTRY;
   d_data_object_id = object_id;
   d_subtree_node.reset();
   d_subtree_node_id = MTreeNode::getUndefinedId();
}

/*
*************************************************************************
*                                                                       *
* Check consistency of this entry.                                      *
*                                                                       *
*************************************************************************
*/
 
bool MTreeEntry::checkConsistency(MTreeNodePtr my_node,
                                  ostream& stream) const
{
   assert(my_node.get());

   bool entry_is_consistent = true;

   if ( d_entry_type == UNDEFINED_ENTRY ) {
      stream << "MTREE ENTRY ERROR: Entry type is undefined!" << endl;
      entry_is_consistent = false;
   }

   if ( !d_my_node.get() ) {
      stream << "MTREE ENTRY ERROR: "
             << "Entry has null pointer to owning node!" << endl;
      entry_is_consistent = false;
   }

   if ( d_my_node.get() != my_node.get() ) {
      stream << "MTREE ENTRY ERROR: "
             << "Entry owning node not same as node "
             << "passed to checkConsistency() method!" << endl;
      entry_is_consistent = false;
   }

   if ( !d_my_node.get() ) {
      stream << "MTREE ENTRY ERROR: "
             << "Entry has null pointer to owning node!" << endl;
      entry_is_consistent = false;
   }
   
   if ( d_my_node->isRoot() ) {
      if ( getDistanceToParent() != 
           MTreeKey::getUndefinedDistanceToParent() ) {
         stream << "MTREE ENTRY ERROR: root node entry distance-to-parent " 
                << " is not undefined!" << endl;
         entry_is_consistent = false;
      } 
   } else {
      double my_dist2par = getDistanceToParent();
      double actual_dist2par = 
         computeDistanceTo( d_my_node->getParentEntry() );
      if ( !toolbox::MathUtilities<double>::equalEps(
            my_dist2par, actual_dist2par) ) {
         stream << "MTREE ENTRY ERROR: distance-to-parent incorrect!" 
                << "\n  my_dist2par = " << my_dist2par << endl
                << "  actual_dist2par = " << actual_dist2par << endl;
         entry_is_consistent = false;
      }
   }

   if ( d_my_node->isLeaf() ) {
      
      if ( d_entry_type != DATA_ENTRY ) {
         stream << "MTREE ENTRY ERROR: "
                << "Leaf node entry type not set to data entry!" << endl;
         entry_is_consistent = false;
      } 

      if ( d_subtree_node.get() ) {
         stream << "MTREE ENTRY ERROR: "
                << "Subtree node set for leaf node entry!" << endl;
         entry_is_consistent = false;
      }

      if ( d_data_object_id == DBObject::getUndefinedId() ) {
         stream << "MTREE ENTRY ERROR: "
                << "Data object id not set for leaf node entry!" << endl;
         entry_is_consistent = false;
      }

   } else {  // internal node entry

      if ( d_entry_type != ROUTING_ENTRY ) {
         stream << "MTREE ENTRY ERROR: "
                << "Non-leaf node entry type not set to routing entry!" 
                << endl;
         entry_is_consistent = false;
      } 

      if ( !d_subtree_node.get() ) {
         stream << "MTREE ENTRY ERROR: "
                << "Subtree node not set for non-leaf node entry!"
                << endl;
         entry_is_consistent = false;
      }

      if ( d_data_object_id != DBObject::getUndefinedId() ) {
         stream << "MTREE ENTRY ERROR: "
                << "Data object id defined for non-leaf node entry!"
                << endl;
         entry_is_consistent = false;
      }

      double subtree_radius = 0.0;
      for (int ie = 0; ie < d_subtree_node->getNumberEntries(); ++ie) {
         subtree_radius = 
            toolbox::MathUtilities<double>::Max(
               subtree_radius, 
               d_subtree_node->getEntry(ie)->getDistanceToParent() +
               d_subtree_node->getEntry(ie)->getRadius() );
      }
      if ( getRadius() < subtree_radius ) {
         stream << "MTREE ENTRY ERROR: Covering radius too small!" << endl;
         stream << "  my_radius = " << getRadius() << endl;
         stream << "  subtree_radius = " << subtree_radius << endl;
         entry_is_consistent = false;
      }

   }

   if ( !entry_is_consistent ) {
      printClassData(stream);
   }

   return( entry_is_consistent );
}

/*
*************************************************************************
*                                                                       *
* Print entry to given output stream.                                   *
*                                                                       *
*************************************************************************
*/

void MTreeEntry::printClassData(ostream& stream) const
{
   stream << "MTreeEntry::printClassData()\n";
   stream << "--------------------------------------\n";
   stream << "this ptr = " << (MTreeEntry*)this << endl;
   stream << "d_entry_type = " << d_entry_type << endl;
   stream << "isDefined() = " << isDefined() << endl;
   stream << "isDataEntry() = " << isDataEntry() << endl;
   stream << "isRoutingEntry() = " << isRoutingEntry() << endl;
   stream << "d_my_node = " << (MTreeNode*)d_my_node.get() << endl;
   stream << "d_my_node_id = " << d_my_node_id << endl;
   stream << "d_my_position_in_node = " << d_my_position_in_node << endl;
   stream << "d_subtree_node = " << (MTreeNode*)d_subtree_node.get() << endl;
   stream << "d_subtree_node_id = " << d_subtree_node_id << endl;
   stream << "d_data_object_id = " << d_data_object_id << endl;

   stream << "\nPrinting entry key..." << endl;
   d_key.printClassData(stream);

   stream << "\n" << endl;
}


#endif




