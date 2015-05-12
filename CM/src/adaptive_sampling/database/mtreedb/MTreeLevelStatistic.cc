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
// File:        MTreeLevelStatistic.cc
// Package:     MTree database
// 
// 
// 
// Description: Simple class for holding MTree statistic for a level in the tree
//

#ifndef included_MtreeLevelStatistic_C
#define included_MTreeLevelStatistic_C

#include "MTreeLevelStatistic.h"

#ifdef DEBUG_CHECK_ASSERTIONS
#ifndef included_cassert
#define included_cassert
#include <cassert>
#endif
#endif

/*
*************************************************************************
*                                                                       *
* Ctor for MTree level statistic object.                                *
*                                                                       *
*************************************************************************
*/

MTreeLevelStatistic::MTreeLevelStatistic(
   int level_number,
   int max_num_nodes) 
:
   d_level_number(level_number)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(level_number >= 0);
   assert(max_num_nodes >= 0);
#endif
   d_num_nodes_added = 0;
   d_node_stats.reserve(max_num_nodes);
}


/*
*************************************************************************
*                                                                       *
* Dtor for MTree level statistic object.                                *
*                                                                       *
*************************************************************************
*/

MTreeLevelStatistic::~MTreeLevelStatistic()
{
   d_level_number = -1;
   d_node_stats.clear();
}

/*
*************************************************************************
*                                                                       *
* Public member functions to access level statistic data.               *
*                                                                       *
*************************************************************************
*/

int MTreeLevelStatistic::getLevelNumber() const
{
   return(d_level_number);
}

int MTreeLevelStatistic::getNumberNodesOnLevel() const
{
   return( static_cast<int>(d_node_stats.size()) );
}

const MTreeNodeStat& MTreeLevelStatistic::getNodeStat(int indx) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(0 <= indx && static_cast<unsigned int>(indx) < d_node_stats.size());
#endif
   return(d_node_stats[indx]);
}

void MTreeLevelStatistic::printClassData(ostream& stream) const
{
   stream << "\nNode statistic data for tree level " 
          << d_level_number << endl;
   stream << "   " << d_num_nodes_added << " nodes..." << endl;
   for (int n = 0; n < d_num_nodes_added; ++n) {
      stream << "\n     Node entry # " << n << endl;
      d_node_stats[n].printClassData(stream);
   }
}

/*
*************************************************************************
*                                                                       *
* Private function to add MTreeNodeStat to this level statistic object. *
*                                                                       *
*************************************************************************
*/

void MTreeLevelStatistic::addNodeStat(MTreeNodePtr node)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node);
   assert(node->getLevelInTree() == d_level_number);
   assert(d_node_stats.capacity() > static_cast<unsigned int>(d_num_nodes_added));
#endif
   MTreeNodeStat nstat(node);
   d_node_stats.push_back(nstat);
   d_num_nodes_added++;
}

/*
*************************************************************************
*                                                                       *
* Private functions for setting number of objects in subtree of         *
* MTreeNodeStat objects.                                                *
*                                                                       *
*************************************************************************
*/

void MTreeLevelStatistic::getSubtreeNodeRealIds(
   int node_indx,
   vector<int>& sub_indx)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node_indx < d_num_nodes_added);
   assert(d_level_number > 0);
#endif
   const MTreeNodeStat& nstat = d_node_stats[node_indx];
   const int n_entries = nstat.d_node->getNumberEntries();
    
   sub_indx.clear();
   sub_indx.reserve(n_entries);

   for (int ie = 0; ie < n_entries; ++ie) {
      sub_indx.push_back( nstat.d_node->
                             getEntry(ie)->
                                getSubtreeNode()->
                                   getNodeId() );
   }

}

int 
MTreeLevelStatistic::getTotalNumberDataObjectsInSubtreeForRealNodeId(
   int real_node_id)
{
   int n_objects = 0;

   bool found = false;
   int n = 0; 
   while ( !found && (n < d_num_nodes_added) ) { 
      if ( d_node_stats[n].getNodeId() == real_node_id ) {
         n_objects = 
            d_node_stats[n].getTotalNumberDataObjectsInSubtree(); 
         found = true;
      }
      n++;
   }
   return( n_objects );
}

void MTreeLevelStatistic::setTotalNumberDataObjectsInSubtree(
   int node_indx,
   int num_objects)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(node_indx < d_num_nodes_added);
#endif
    d_node_stats[node_indx].
       setTotalNumberDataObjectsInSubtree(num_objects);
 
}

#endif




