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
// File:        MTreeNodeStat.cc
// Package:     MTree database
// 
// 
// 
// Description: Simple class for holding MTree statistic data for a node in tree
//

#ifndef included_MtreeNodeStat_C
#define included_MTreeNodeStat_C

#include "MTreeNodeStat.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif
#endif

/*
*************************************************************************
*                                                                       *
* Ctors for MTreeNodeStat object.                                       *
*                                                                       *
*************************************************************************
*/

MTreeNodeStat::MTreeNodeStat()
{
}

MTreeNodeStat::MTreeNodeStat(MTreeNodePtr node)
:
   d_node(node),
   d_total_objects_in_subtree(0)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node);
#endif

   if ( d_node->isLeaf() ) {
      const int nentries = d_node->getNumberEntries();
      d_node_object_ids.reserve(nentries);
      for (int ie = 0; ie < nentries; ++ie) {
         d_node_object_ids.push_back(d_node->getEntry(ie)->getDataObjectId());
      }
      d_total_objects_in_subtree = nentries;
   }
}

MTreeNodeStat::MTreeNodeStat(const MTreeNodeStat& nstat)
:
   d_node(nstat.d_node),
   d_total_objects_in_subtree(nstat.d_total_objects_in_subtree),
   d_node_object_ids(nstat.d_node_object_ids)
{
}

/*
*************************************************************************
*                                                                       *
* Dtor for MTreeNodeStat object.                                        *
*                                                                       *
*************************************************************************
*/

MTreeNodeStat::~MTreeNodeStat()
{
   d_node.reset();
   d_node_object_ids.clear();
}

/*
*************************************************************************
*                                                                       *
* Copy assignment operatoy for MTreeNodeStat object.                    *
*                                                                       *
*************************************************************************
*/

MTreeNodeStat& MTreeNodeStat::operator=(const MTreeNodeStat& rhs)
{
   d_node                     = rhs.d_node;
   d_total_objects_in_subtree = rhs.d_total_objects_in_subtree;
   d_node_object_ids          = rhs.d_node_object_ids;
   return(*this);
}

/*
*************************************************************************
*                                                                       *
* Public member functions to access node statistic data.                *
*                                                                       *
*************************************************************************
*/

int MTreeNodeStat::getNodeId() const 
{
   return( d_node->getNodeId() );
}

int MTreeNodeStat::getLevelNumber() const
{
   return( d_node->getLevelInTree() );
}

int MTreeNodeStat::getNumberEntries() const
{
   return( d_node->getNumberEntries() );
}

bool MTreeNodeStat::isRoot() const
{
   return( d_node->isRoot() );
}

bool MTreeNodeStat::isLeaf() const
{
   return( d_node->isLeaf() );
}

double MTreeNodeStat::getCoveringRadius() const
{
   double radius = -1.0;
   if ( !d_node->isRoot() ) {
      radius = d_node->getParentEntry()->getRadius();
   }
   return( radius );
}

MetricSpacePointPtr MTreeNodeStat::getCenterPoint() const
{
   MetricSpacePointPtr point;
   if ( !d_node->isRoot() ) {
      point = d_node->getParentEntry()->getPoint();
   }
   return( point );
}

int MTreeNodeStat::getTotalNumberDataObjectsInSubtree() const
{
   return( d_total_objects_in_subtree );
}

const vector<int>& MTreeNodeStat::getDataObjectIds() const
{
   return( d_node_object_ids ); 
}

void MTreeNodeStat::printClassData(ostream& stream) const
{
   stream << "        Node id = " << getNodeId() << endl;
   stream << "        Num entries = " << getNumberEntries() << endl;
   stream << "        Num data objects in subtree = "
          << getTotalNumberDataObjectsInSubtree() << endl;
   stream << "        Root node = " << isRoot() << endl;
   stream << "        Leaf node = " << isLeaf() << endl;

   if ( !isRoot() ) {
      stream << "        Covering radius = " 
             << getCoveringRadius() << endl;
      stream << "        Center point = "; 
                getCenterPoint()->print(stream);
   }

   if ( isLeaf() ) {
      stream << "        Owns objects ";
      unsigned int ie = 0;
      for ( ; ie < d_node_object_ids.size() - 1; ++ie) {
         stream << d_node_object_ids[ie] << " , ";
      }
      if (ie > 0) {
         stream << d_node_object_ids[ie] << endl;
      }
   }
}

void MTreeNodeStat::setTotalNumberDataObjectsInSubtree(int n)
{
   d_total_objects_in_subtree = n;
}

#endif




