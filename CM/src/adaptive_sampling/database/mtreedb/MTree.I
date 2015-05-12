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
// File:        MTree.I
// Package:     MTree database
// Description: Main Mtree index structure class.
//

#ifdef DEBUG_NO_INLINE
#define inline
#endif

/*
*************************************************************************
*                                                                       *
* inline methods to set node parameters.                                * 
*                                                                       *
*************************************************************************
*/

inline 
void MTree::setNodePromotionMethod(
   MTreeNode::MTreeNodePromotionMethod method)
{
   d_node_promotion_method = method;
}

inline 
void MTree::setRootNodePromotionMethod(
   MTreeNode::MTreeNodePromotionMethod method)
{
   switch(method) {

      case MTreeNode::RANDOM_PROMOTION:
      case MTreeNode::MIN_RADII_PROMOTION:
      case MTreeNode::MIN_OVERLAP_PROMOTION: {
         d_root_node_promotion_method = method;
         break;
      }  // valid root promotion options

      default: {
      }  // use default method if not acceptable

   }
}

inline 
void MTree::setNodePartitionMethod(
   MTreeNode::MTreeNodePartitionMethod method,
   double min_utilization)
{
   d_node_partition_method = method;
   if ( (0.0 < min_utilization) && (min_utilization <= 0.5) ) {
      d_min_node_utilization = min_utilization;
   }
}

/*
*************************************************************************
*                                                                       *
* inline methods to get insert operation count information.             *
*                                                                       *
*************************************************************************
*/

inline
int MTree::getTotalInsertCount() const
{
   return(d_num_inserts);
}

inline
int MTree::getTotalInsertDistanceCount() const
{
   return(d_total_distance_comps_in_inserts);
}

inline
int MTree::getLastInsertDistanceCount() const
{
   return(d_num_distance_comps_in_last_insert);
}

/*
*************************************************************************
*                                                                       *
* inline methods to get delete operation count information.             *
*                                                                       *
*************************************************************************
*/

inline
int MTree::getTotalDeleteCount() const
{
   return(d_num_deletes);
}

inline
int MTree::getTotalDeleteDistanceCount() const
{
   return(d_total_distance_comps_in_deletes);
}

inline
int MTree::getLastDeleteDistanceCount() const
{
   return(d_num_distance_comps_in_last_delete);
}

/*
*************************************************************************
*                                                                       *
* inline methods to get KNN search operation count information.         *
*                                                                       *
*************************************************************************
*/

inline
int MTree::getTotalKNNSearchCount() const
{
   return(d_num_knn_queries);
}

inline
int MTree::getTotalKNNSearchDistanceCount() const
{
   return(d_total_distance_comps_in_knn_queries);
}

inline
int MTree::getLastKNNSearchDistanceCount() const
{
   return(d_num_distance_comps_in_last_knn_query);
}

/*
*************************************************************************
*                                                                       *
* inline methods to get range search operation count information.       *
*                                                                       *
*************************************************************************
*/

inline
int MTree::getTotalRangeSearchCount() const
{
   return(d_num_range_queries);
}

inline
int MTree::getTotalRangeSearchDistanceCount() const
{
   return(d_total_distance_comps_in_range_queries);
}

inline
int MTree::getLastRangeSearchDistanceCount() const
{
   return(d_num_distance_comps_in_last_range_query);
}


/*
*************************************************************************
*                                                                       *
* Private inline methods to get node parameters.                        * 
*                                                                       *
*************************************************************************
*/

inline
MTreeNode::MTreeNodePromotionMethod 
MTree::getRootNodePromotionMethod() const
{
   return( d_root_node_promotion_method );
}

inline
MTreeNode::MTreeNodePromotionMethod 
MTree::getNodePromotionMethod() const
{
   return( d_node_promotion_method );
}

inline
MTreeNode::MTreeNodePartitionMethod 
MTree::getNodePartitionMethod() const
{
   return( d_node_partition_method );
}

inline
double MTree::getMinNodeUtilization() const
{
   return( d_min_node_utilization );
}

/*
*************************************************************************
*                                                                       *
* Private inline members for incrementing and decrementing              *
* level node counts.                                                    *
*                                                                       *
*************************************************************************
*/

inline
void MTree::addNodeToLevelCount(int level)
{
   d_number_nodes_in_level[level]++;
}

inline
void MTree::removeNodeFromLevelCount(int level)
{
   d_number_nodes_in_level[level]--;
}


#ifdef DEBUG_NO_INLINE
#undef inline
#endif




