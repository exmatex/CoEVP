//
// File:        MTreeLevelStatistic.h
// Package:     MTree database
// 
// 
// 
// Description: Simple class for holding MTree statistic for a level in the tree
//

#ifndef included_MTreeLevelStatistic
#define included_MTreeLevelStatistic

#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif
#ifndef included_MTreeNode
#include "MTreeNode.h"
#endif
#ifndef included_MTreeNodeStat
#include "MTreeNodeStat.h"
#endif

#ifndef included_vector
#define included_vector
#include <vector>
using namespace std;
#endif


class MTreeNodeStat;

/*!
 * @brief MTreeLevelStatistic is a simple class for storing
 * statistics for a single MTree level.
 * 
 * @see MetricSpacePoint 
 */

class MTreeLevelStatistic
{
public:
   friend class MTree;

   /*!
    * Default ctor for MTreeLevelStatistic.
    */
   MTreeLevelStatistic(int level_number,
                       int max_num_nodes);

   /*!
    * Dtor for MTree search result objects.
    */
   virtual ~MTreeLevelStatistic();

   /*!
    * Return level number for this level statistics object.
    */
   int getLevelNumber() const; 

   /*!
    * Return number of node statistics associated with this 
    * level statistics object.
    */
   int getNumberNodesOnLevel() const;

   /*!
    * Return const reference to MTreeNodeStat object with given
    * index.  
    *
    * Note that valid indices are in range 0 to getNumberNodesOnLevel()-1.
    */
   const MTreeNodeStat& getNodeStat(int indx) const;

   /*!
    * Print level node statistic data to the specified output stream.
    */
   void printClassData(ostream& stream) const;

private:
   // The following are not implemented
   MTreeLevelStatistic();
   MTreeLevelStatistic(const MTreeLevelStatistic&);
   void operator=(const MTreeLevelStatistic&);

   /*
    * Private member functions used by MTree class for setting level 
    * statistics.
    */
   void addNodeStat(MTreeNodePtr node);
   void getSubtreeNodeRealIds(int node_indx,
                              vector<int>& sub_indx);
   int getTotalNumberDataObjectsInSubtreeForRealNodeId(int real_node_id);
   void setTotalNumberDataObjectsInSubtree(int node_indx,
                                           int num_objects);

   int d_level_number;
   int d_num_nodes_added;
   vector<MTreeNodeStat> d_node_stats;
};

#endif
