/* DO-NOT-DELETE revisionify.begin() */
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
/* DO-NOT-DELETE revisionify.end() */
//
// File:        MTreeNode.I
// Package:     MTree database
// 
// 
// 
// Description: Representation of node in an MTree.
//

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif
#endif

#ifndef included_toolbox_Utilities
#include "toolbox/base/Utilities.h"
#endif

#ifndef included_MTree
#include "MTree.h"
#endif

#ifdef DEBUG_NO_INLINE
#define inline
#endif

/*
*************************************************************************
*                                                                       *
* Accessory functions to get and set node id and get shared pointer     *
* to this object.                                                       *
*                                                                       *
*************************************************************************
*/

inline
int MTreeNode::getNodeId() const
{
   return(d_node_id);
}

/*
*************************************************************************
*                                                                       *
* Accessory routines to check type of node and whether it is defined,   * 
* and also to set node level in tree.                                   *
*                                                                       *
*************************************************************************
*/

inline
bool MTreeNode::isDefined() const
{
   return( d_node_id > MTreeNode::getUndefinedId() );
}

inline
bool MTreeNode::isRoot() const
{
   return( d_is_root_node );
}

inline
bool MTreeNode::isLeaf() const
{
   return( d_level_in_tree == 0 );
}

inline
int MTreeNode::getMaxEntries() const
{
   return( d_max_entries );
}

inline
void MTreeNode::setRootNode(bool is_root)
{
   d_is_root_node = is_root;
   d_parent_entry.reset();
}

inline
int MTreeNode::getLevelInTree() const
{
   return( d_level_in_tree );
}

/*
*************************************************************************
*                                                                       *
* Accessory functions to check number of entries in node and to         *
* retrieve the entry at a given position.                               *
*                                                                       *
*************************************************************************
*/

inline
int MTreeNode::getNumberEntries() const
{
   return( d_entries.size() );
}

inline
bool MTreeNode::isOverfull() const
{
   const int num_entries = d_entries.size();
   return( num_entries > d_max_entries );
}

inline
MTreeEntryPtr MTreeNode::getEntry(int position) const
{
   MTreeEntryPtr ret_entry;
   if ( (position >= 0) && (position < getNumberEntries()) ) {
      ret_entry = d_entries[position];
   }
   return( ret_entry );
}

/*
*************************************************************************
*                                                                       *
* Accessory function to get parent node of this node if it exists.      *
*                                                                       *
*************************************************************************
*/

inline
MTreeNodePtr MTreeNode::getParentNode() const
{
   MTreeNodePtr ret_node;
   if ( d_parent_entry.get() ) {
      ret_node = d_parent_entry->getNode(); 
   }
   
   return( ret_node );
}

/*
*************************************************************************
*                                                                       *
* Accessory function to set and get parent entry of this node.          *
*                                                                       *
*************************************************************************
*/

inline
void MTreeNode::setParentEntry(MTreeEntryPtr parent)
{
   if ( isRoot() ) {
      TBOX_ERROR("MTreeNode::setParentEntry() error..."
                 << " Cannot set parent entry of root node!" << endl);
   }
   d_parent_entry = parent;
}

inline
MTreeEntryPtr MTreeNode::getParentEntry() const
{
   return( d_parent_entry );
}

/*
*************************************************************************
*                                                                       *
* Private accessory function to set and get leaf node id.               *
*                                                                       *
*************************************************************************
*/

inline
int MTreeNode::getLeafNodeId() const
{
   return(d_leaf_node_id);
}

inline
void MTreeNode::setLeafNodeId(int id)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(d_leaf_node_id == MTreeNode::getUndefinedId());
   assert(id >= 0);
#endif
   d_leaf_node_id = id;
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif




