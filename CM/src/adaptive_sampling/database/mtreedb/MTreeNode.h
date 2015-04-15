//
// File:        MTreeNode.h
// Package:     MTree database
// 
// 
// 
// Description: Representation of node in an MTree.
//

#ifndef included_MTreeNode
#define included_MTreeNode

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_vector
#define included_vector
#include <vector>
using namespace std;
#endif

#ifndef included_MTreeEntry
#include "MTreeEntry.h"
#endif


class MTree;

class MTreeNode;
typedef std::shared_ptr<MTreeNode> MTreeNodePtr;

/*!
 * @brief MTreeNode represents a node in an MTree.  Each node holds a 
 * set of (data or routing) entries in the tree.
 *
 * Upon creation, a unique node id is assigned to the node and the node 
 * is marked as a leaf node in the tree.  A node is changed from a leaf 
 * node by calling the method setLevelInTree() with an integer value > 0. 
 *
 * A node is set to be the root node of the MTree by calling the 
 * setRootNode() method with a true boolean value argument.  
 *
 * @see MTreeEntry
 */

class MTreeNode
{
public:
   friend class MTreeDataStore;

   /*!
    * Static function to get integer identifier for undefined node
    * common to all node objects.
    */
   static int getUndefinedId();

   /*!
    * Constructor for MTreeNode class. 
    * Upon creation, a unique node id is assigned and node is marked 
    * as a leaf node in the tree.  
    *
    * @param tree  Bare pointer to MTree owning this node.
    * @param max_entries  Integer maximum number of entries that node 
    *                     may own.  When assertion checking is on, an 
    *                     assertion will result if value is < 2.
    */
   MTreeNode(MTree* tree, 
             int max_entries);

   /*!
    * Enumerated type for entry promotion method options.  These options 
    * are the simplest and least expensive computationally among those 
    * described in the literature (except for MIN_OVERLAP_PROMOTION, 
    * which is not found in the literature).  Other, more sophisticated and 
    * expensive options may provide better tree searching performance.
    * 
    * "Confirmed" promotion indicates that only one object will be promoted 
    * and the split node will keep its parent entry.  Otherwise, two objects 
    * will be promoted and the split node may receive a new parent entry.  
    * Thus, "confirmed" promotion and MAX_SPREAD_DISTANCE_PROMOTION are not 
    * allowed in the root node, since it has no parent entries before split. 
    *
    * The choices are:
    * 
    * -# RANDOM_PROMOTION  
    *      Choose two entries in a node randomly to promote.
    * -# MAX_SPREAD_DISTANCE_PROMOTION  
    *      [default method for all non-root nodes!]
    *      Choose the two entries in a node having largest distances to
    *      parent entry of the node. 
    * -# MIN_RADII_PROMOTION 
    *      Choose the two entries in a node having the two minimum radii
    *      among all entries in the node.
    * -# MIN_OVERLAP_PROMOTION
    *      [default method for root node!]
    *      Choose the two entries E_i, E_j that maximize the value 
    *      d_ij - (rad_i + rad_j), where d_ij is the distance between 
    *      E_i, E_j and rad_i, rad_j are the radii of E_i, E_j.
    * -# CONFIRMED_RANDOM_PROMOTION
    *      Keep the current parent of a node and choose the other entry to
    *      promote randomly.
    * -# CONFIRMED_MAX_LB_DISTANCE_PROMOTION
    *      Keep the current parent of a node and choose the other entry to
    *      promote having the largest distances to parent entry of the node.
    */
   enum MTreeNodePromotionMethod { RANDOM_PROMOTION = 0,
                                   MAX_SPREAD_DISTANCE_PROMOTION = 1,
                                   MIN_RADII_PROMOTION = 2,
                                   MIN_OVERLAP_PROMOTION = 3,
                                   CONFIRMED_RANDOM_PROMOTION = 4,
                                   CONFIRMED_MAX_LB_DISTANCE_PROMOTION = 5 };

   /*!
    * Enumerated type for partitioning node entries during a node split. 
    * Others are possible.  Only these options are implemented.  
    *
    * The choices are:
    *
    * -# HYPERPLANE_PARTITION [default method!]
    *      Assign each entry to the node with the nearest promoted entry.  
    *      Generally, this results in an unbalanced split.  However, 
    *      a minimum node utilization fraction can be applied to control 
    *      the amount of imbalance.
    * -# BALANCED_PARTITION
    *      Assign the entry closest to the first promoted entry to the first 
    *      node.  Then, assign the entry closest to the second promoted 
    *      entry to the second node.  Continue until all entries 
    *      are assigned.  The result is a balanced split. 
    */
   enum MTreeNodePartitionMethod { HYPERPLANE_PARTITION = 0,
                                   BALANCED_PARTITION = 1 };

   /*!
    * Dtor for MTree node objects.
    */
   virtual ~MTreeNode();

   /*!
    * Function to get proper shared pointer to this MTreeNode object.
    */
   MTreeNodePtr getThisPtr() const;

   /*!
    * Return integer identifier of node.  If node is undefined, the value
    * provided by MTreeNode::getUndefinedId() is returned.
    */
   int getNodeId() const;

   /*!
    * Return leaf node integer identifier.  If node is undefined, the value
    * provided by MTreeNode::getUndefinedId() is returned.  Note: the 
    * leaf node id is distinct from the node id and should only be used by the 
    * data store.
    */
   int getLeafNodeId() const;

   /*!
    * Return true if the node has a valid node id (i.e., node id is not 
    * equal to MTreeNode::getUndefinedId()).
    */
   bool isDefined() const;

   /*!
    * Return true if node is the root of the tree; otherwise, return false.
    */
   bool isRoot() const;

   /*!
    * Return true if node is a leaf node in the tree; otherwise, return false.
    * Note: isLeaf() == ( (getLevelInTree() == 0) ).
    */
   bool isLeaf() const;

   /*!
    * Return max number of entries allowed in node; set in constructor.
    */
   int getMaxEntries() const;

   /*!
    * Set node to be the root node or not the root node of the tree
    * depending on value of boolean argument.
    */
   void setRootNode(bool is_root);

   /*!
    * Set level in tree to given level number.  Note: a level number of zero
    * indicates a leaf node, a level number > zero indicates a routing node.
    *
    * When assrtion checking is on, an assertion is thrown when level < 0.
    */
   void setLevelInTree(int level);

   /*!
    * Return level of node in treer.  Note: a level number of zero indicates 
    * a leaf node, a level number > zero indicates a routing node.
    */
   int getLevelInTree() const;

   /*!
    * Return pointer to parent node of this node.  If this node is root node,
    * then a null pointer is returned.
    */
   MTreeNodePtr getParentNode() const;

   /*!
    * Set parent entry of this node to given node.  The parent entry should
    * only be set for non-root nodes.  If this method is called on the root
    * node, an error results.
    */
   void setParentEntry(MTreeEntryPtr parent);

   /*!
    * Return pointer to parent entry of this node.  If this node is root 
    * node, then a null pointer is returned.
    */
   MTreeEntryPtr getParentEntry() const;

   /*!
    * Return integer number of entries currently in node.
    */
   int getNumberEntries() const;

   /*!
    * Return true if node has more entries than max specified in
    * call to MTreeNode::createNode(); otherwise, return false.
    * A return value of true indicates that the node must be split.
    */
   bool isOverfull() const;

   /*!
    * Return pointer to node entry at given position.  If position is 
    * invalid for node, then a null pointer is returned.
    */ 
   MTreeEntryPtr getEntry(int position) const;

   /*!
    * Search entries of this node for one best suited for
    * insertion of entry.  
    * 
    * @return Pointer to best subtree node among entries in this node.  
    *         Note that this method should not be called on a leaf node 
    *         in the tree.  If it is, it will return a null pointer.
    *
    * @param  entry Pointer to new entry.  When assertion checking is on,
    *         passing a null pointer will throw an assertion.
    */
   MTreeNodePtr searchBestSubtreeForInsert(MTreeEntryPtr entry) const;

   /*!
    * Determine position in node for new entry and insert.  This routine
    * is static so that pointer to node for insertion can be passed
    * (circumventing need to obtain shared pointer to this pointer).
    *
    * Note that this routine does not re-compute radius information for 
    * the node or other nodes that may be affected by the insertion.
    *
    * @param  entry Pointer to new entry.  When assertion checking 
    *         is on, passing a null pointer will throw an assertion.
    * @param  node Pointer to node into which to insert entry.  When 
    *         assertion checking is on, passing a null pointer will 
    *         throw an assertion.
    */
   static void insertEntry(MTreeEntryPtr entry,
                           MTreeNodePtr node);
   
   /*!
    * Insert entry in node at specified position.  This routine
    * is static so that pointer to node for insertion can be passed
    * (circumventing need to obtain shared pointer to this pointer).
    *
    * Note that this routine does not re-compute radius information for 
    * the node or other nodes that may be affected by the insertion.
    *
    * @param  entry Pointer to new entry.  When assertion checking 
    *         is on, passing a null pointer will throw an assertion.
    * @param  position Integer position in node of new entry.  When 
    *         assertion checking is on, this routine throws an assertion 
    *         if inserting the entry at the given position will break 
    *         the node ordering in any way.
    * @param  node Pointer to node into which to insert entry.  When 
    *         assertion checking is on, passing a null pointer will 
    *         throw an assertion.
    */
   static void insertEntryAtPosition(MTreeEntryPtr entry, 
                                     int position,
                                     MTreeNodePtr node);

   /*!
    * Delete given entry from node.  This routine is static so that 
    * pointer to node for entry deletion can be passed
    * (circumventing need to obtain shared pointer to this pointer).
    *
    * If the entry is not assigned to the node, then the method does 
    * nothing.
    *
    * Note that this routine does not re-compute radius information for 
    * the node or other nodes that may be affected by the entry deletion.
    *
    * @param  entry Pointer to entry.  When assertion checking 
    *         is on, passing a null pointer will throw an assertion.
    * @param  node Pointer to node from which to delet entry.  When 
    *         assertion checking is on, passing a null pointer will 
    *         throw an assertion.
    */
   static void deleteEntry(MTreeEntryPtr entry, 
                           MTreeNodePtr node);

   /*!
    * Transfer entries from this node to another node.
    *
    * This rotuine assumes argument node contains no entries. 
    *
    * This routine assumes nothing about the order of the entries in this
    * node nor whether their distance-to-parent values have been properly
    * set.  On exit, this node has no entries, and the other node owns all
    * entries.  The order of entries is preserved in the other node.  No
    * state of entries has been changed.
    *
    * @param  node Pointer to node to which all entries are transferred.
    */
   void transferEntries(MTreeNodePtr node);

   /*!
    * Split this node and return the new node resulting from the split.
    *
    * If this node is the root of the tree, the method does nothing 
    * and returns a null pointer.  Otherwise, if this node is not
    * over-full and the "special_delete_split" flag is false, the
    * method does nothing and returns a null pointer.
    * 
    * When this node is split, its entries are partitioned between it and
    * the new node.  Parent entries for this node and the new node are 
    * promoted to the parent node of this node.  All distance-to-parent 
    * values are reset for all entries in the node when the method 
    * was called.  The distance-to-parent values are also reset for the 
    * entries promoted to the parent node.
    *
    * @return Pointer to new node resulting from the split.
    * 
    * @param use_root_promotion Optional boolean indicating whether 
    *        root node promotion method should be used for split.  
    *        Dafault is false.
    * @param special_delete_split Optional boolean forcing a split for
    *        a non-root node that is not over-full.  Dafault is false.
    *        This flag should be set to true with care and is generally
    *        used for the special case of moving entries to an empty
    *        node after deleting objects from the tree.
    */
   MTreeNodePtr split(bool use_root_promotion = false,
                      bool special_delete_split = false);

   /*!
    * Identify two entries in the node for promotion using the
    * given method.
    * 
    * @return Boolean false if entry1 is the same as the parent entry of
    *         the node; otherwise, true (new parent entry).
    *
    * @param  entry1  Pointer to promoted entry to represent this node.
    * @param  entry2  Pointer to promoted entry to represent another node
    *                 (e.g., during a node split).
    * @param  method  enum type MTreeNodePromotionMethod value indicating
    *                 promotion method to use.
    */
   bool promoteTwo(MTreeEntryPtr& entry1,
                   MTreeEntryPtr& entry2,
                   MTreeNodePromotionMethod method) const;

   /*!
    * Identify one entry in the node for promotion using the
    * given method.
    *
    * @param  entry2  Pointer to promoted entry.  It will be different
    *                 than parent of node.
    * @param  method  enum type MTreeNodePromotionMethod value indicating
    *                 promotion method to use.
    */
   void promoteOne(MTreeEntryPtr& entry2,
                   MTreeNodePromotionMethod method) const;

   /*!
    * Reset radii for all routing entries along path from this node 
    * up to the root node.
    */
   void resetRadiiUpToRoot();

   /*!
    * Reset radius for routing entry associated with this node.
    *
    * Method asssumes each node entry has correct radius and
    * distance-to-parent value set.
    */
   void resetRadius();

   /*!
    * Check consistency of nodes and their entries starting
    * at this node and continuing along path up to root node.
    *
    * @return false if no inconsistencies found; otherwise true.
    *
    * @param stream Output stream to which all inconsistencies
    *               will be printed.
    */
   bool checkConsistencyUpToRoot(ostream& stream) const;

   /*!
    * Check consistency of this node and its entries.
    *
    * @return false if no inconsistencies found; otherwise true.
    *
    * @param stream Output stream to which all inconsistencies
    *               will be printed.
    */
   bool checkConsistency(ostream& stream) const;

   /*!
    * Print node summary to the given output stream.
    */
   void printSummary(ostream& stream) const;

   /*!
    * Print all MTree node data to the specified output stream. 
    */
   void printClassData(ostream& stream) const;

private:
   // The following are not implemented
   MTreeNode();
   MTreeNode(const MTreeNode&);
   void operator=(const MTreeNode&);

   /*
    * Private method to identify two entries in the node for promotion.
    *
    * The boolean return value is false if entry1 is the same as the 
    * parent entry of the node and true otherwise (new parent entry).
    */ 
   bool promoteRoutingEntries(MTreeEntryPtr& entry1,
                              MTreeEntryPtr& entry2,
                              bool use_root_promotion,
                              bool special_delete_split) const;

   /*
    * Private method to partition entries between this node and node2.
    * 
    * Partition method used depends on the value of s_partition_method.
    */
   void partitionEntries(MTreeEntryPtr entry1,
                         bool new_parent_entry,
                         MTreeEntryPtr entry2,
                         MTreeNodePtr& node2);

   /*
    * Set leaf node id to given integer value.  Generally, this method  
    * should only be called once.  To avoid some (but not all) problems, 
    * this function will throw an assertion if it is called and the node 
    * already has a valid leaf node id.  If this occurs, there is a 
    * likely a problem somewhere.  Note: the leaf node id is distinct 
    * from the node id and is only used by the data store.
    */
   void setLeafNodeId(int id);

   /*
    * Get smart pointer to myself.
    */
   MTreeNodePtr getMyself() const; 

   /*
    * Undefined node id shared by all node instances;
    * cannot be changed.
    */
   static int s_undefined_node_id;

   /*
    * Static node instance counter shared by all node instances;
    * used for generating unique node id.
    */
   static int s_node_instance_counter;

   /*
    * Pointer to MTree used to access node split algorithm parameters
    * and update statistics.
    */
   MTree* d_my_tree;

   /*
    * true if node is root of tree; otherwise false.
    * false by default.
    */   
   bool   d_is_root_node;            

   /*
    * Unique node identifier set when node is constructed.
    */
   int    d_node_id;

   /*
    * Unique leaf node identifier set by data store when data 
    * objects are assigned to leaf node; undefined by default.
    */
   int    d_leaf_node_id;

   /*
    * Level of node in tree (zero for leaf nodes and > zero for 
    * non-leaf nodes); zero by default.
    * Note: node level decreases from root to leaf.
    */
   int    d_level_in_tree;           

   /*
    * pointer to entry holding this node as subtree;
    * null by default. 
    * note should be null only when root node.
    */
   MTreeEntryPtr   d_parent_entry;   

   /*
    * Maximum number of entries allowed in node (set in constructor)
    * and set of entries owned by this node.
    */
   int    d_max_entries;            
   vector<MTreeEntryPtr> d_entries;  

};

#ifndef DEBUG_NO_INLINE
#include "MTreeNode.I"
#endif
#endif
