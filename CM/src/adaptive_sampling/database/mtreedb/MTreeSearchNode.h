//
// File:        MTreeSearchNode.h
// Package:     MTree database
// 
// 
// 
// Description: Utility class used in MTree searches
//

#ifndef included_MTreeSearchNode
#define included_MTreeSearchNode

#ifndef included_MTreeNode
#include "MTreeNode.h"
#endif

/*!
 * @brief MTreeSearchNode keeps track of an MTreeNode to search 
 * during a nearest-neighbor search. 
 * 
 * This class is used in MTree search operations and should not
 * be used for other stuff.
 */
 
class MTreeSearchNode
{
public:
   friend class MTreeSearchQueue;

   /*!
    * Ctor for MTreeSearchNode sets data members to given 
    * values (no error checking).
    */
   MTreeSearchNode(MTreeNodePtr node,
                   double bound,
                   double distance);

   /*!
    * Copy ctor for MTreeSearchNode (no error checking).
    */
   MTreeSearchNode(const MTreeSearchNode& node);

   /*!
    * Copy assignment operator for MTreeSearchNode (no error checking).
    */
   MTreeSearchNode operator=(const MTreeSearchNode& rhs);
 
   /*!
    * Dtor for MTreeSearchNode. 
    */
   virtual ~MTreeSearchNode();

   /*!
    * Return pointer to node to search.
    */
   MTreeNodePtr getSearchNode() const;

   /*!
    * Return bound on search distance.
    */
   double getBound() const;

   /*!
    * Return search distance.
    */
   double getDistance() const;

private:
   // The following is not implemented
   MTreeSearchNode();
   
   MTreeNodePtr    d_search_node;
   double          d_bound;
   double          d_distance;
 
};

#ifndef DEBUG_NO_INLINE
#include "MTreeSearchNode.I"
#endif
#endif
