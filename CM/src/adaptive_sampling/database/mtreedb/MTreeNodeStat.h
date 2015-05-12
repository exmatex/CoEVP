//
// File:        MTreeNodeStat.h
// Package:     MTree database
// 
// 
// 
// Description: Simple class for holding MTree statistic data for a node in tree
//

#ifndef included_MTreeNodeStat
#define included_MTreeNodeStat

#ifndef included_MetricSpacePoint
#include <base/MetricSpacePoint.h>
#endif
#ifndef included_MTreeNode
#include "MTreeNode.h"
#endif

#ifndef included_vector
#define included_vector
#include <vector>
using namespace std;
#endif


/*!
 * @brief MTreeNodeStat is a simple class for accessing
 * statistical informatio for a single MTree node.
 * 
 * @see MTreeNode 
 * @see MetricSpacePoint 
 */

class MTreeNodeStat
{
public:
   friend class MTreeLevelStatistic;

   /*!
    * Default ctor for MTreeNodeStat.
    */
   MTreeNodeStat();

   /*!
    * Ctor for MTreeNodeStat.
    */
   MTreeNodeStat(MTreeNodePtr node);

   /*!
    * Copy ctor for a MTreeNodeStat.
    */
   MTreeNodeStat(const MTreeNodeStat& nstat);

   /*!
    * Dtor for MTreeNodeStat.
    */
   virtual ~MTreeNodeStat();

   /*!
    * MTreeNodeStat copy assignment operator.
    */
   MTreeNodeStat& operator=(const MTreeNodeStat& rhs);

   /*!
    * Return id of actual node in tree associated with this MTreeNodeStat.
    */
   int getNodeId() const;

   /*!
    * Return level number for node associated with this MTreeNodeStat.
    */
   int getLevelNumber() const; 

   /*!
    * Return number of entries for node associated with this MTreeNodeStat.
    */
   int getNumberEntries() const; 

   /*!
    * Return true if node associated with this MTreeNodeStat is root 
    * node of the tree; otherwise, return false.
    */
   bool isRoot() const;

   /*!
    * Return true if node associated with this MTreeNodeStat is a leaf 
    * node in the tree; otherwise, return false.
    * Note: isLeaf() == ( (getLevelNumber() == 0) ).
    */
   bool isLeaf() const;

   /*!
    * Return covering radius of node associated with this MTreeNodeStat.
    * Note: if node is root of tree, it has no covering radius and 
    * an invalid value ( < 0.0 ) is returned.
    */
   double getCoveringRadius() const;

   /*!
    * Return smart pointer to center point of node associated with 
    * this MTreeNodeStat.  Note: if node is root of tree, it has no 
    * center point and a null pointer is returned.
    */
   MetricSpacePointPtr getCenterPoint() const;

   /*!
    * Return total number of data objects in all leaf nodes of node
    * associated with this MTreeNodeStat.
    */
   int getTotalNumberDataObjectsInSubtree() const;

   /*!
    * Return vector of data object ids in node associated with
    * this MTreeNodeStat.  Note: if node is not a leaf tree, it owns no
    * data objects directly and an empty vector is returned.
    */
   const vector<int>& getDataObjectIds() const;

   /*!
    * Print level node statistic data to the specified output stream.
    */
   void printClassData(ostream& stream) const;

private:
   void setTotalNumberDataObjectsInSubtree(int n);

   MTreeNodePtr  d_node;
   int           d_total_objects_in_subtree;
   vector<int>   d_node_object_ids;
};

#endif
