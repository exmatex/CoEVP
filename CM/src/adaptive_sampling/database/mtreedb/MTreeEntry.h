//
// File:        MTreeEntry.h
// Package:     MTree database
// 
// 
// 
// Description: Representation of data entry or routing entry in MTree
//

#ifndef included_MTreeEntry
#define included_MTreeEntry

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif

#ifndef included_MTreeKey
#include "MTreeKey.h"
#endif
#ifndef included_DBObject
#include <base/DBObject.h>
#endif

class MTreeEntry;
typedef std::shared_ptr<MTreeEntry> MTreeEntryPtr;

// This forward declaration is kinda gross, but it eliminates a
// circular dependency.
class MTreeNode;
typedef std::shared_ptr<MTreeNode> MTreeNodePtr;

/*!
 * @brief MTreeEntry represents an MTree entry as a key-object pair stored 
 * in an MTreeNode.  The key (an MTreeKey object) defines a region in the 
 * MTree. The object is an actual data item in the tree if the node owning 
 * the entry is a leaf.  Otherwise, the node owning the entry is an internal
 * node in the tree and the object is a node holding the subtree associated 
 * with the entry.
 *
 * When the entry holds a data object and lives in a leaf node, the member
 * function isDataEntry() will return true.  Otherwise, it will return false.
 *
 * When the entry holds a subtree node, the member function isRoutingEntry() 
 * will return true.  Otherwise, it will return false.
 *
 * The member function isDefined() will return true if either 
 * isRoutingEntry() or isDataEntry() would return true.  Otherwise, it 
 * will return false and the entry is undefined.
 * 
 * @see MTreeKey
 */

class MTreeEntry
{
public:
   /*!
    * Ctor for a MTree entry object sets entry key and entry type undefined.
    */
   MTreeEntry(const MTreeKey& key);

   /*!
    * Dtor for MTree entry objects.
    */
   virtual ~MTreeEntry();

   /*!
    * Return true if the entry holds a data object and lives in a leaf 
    * node of a tree (isDataEntry() would return true) or if the entry 
    * holds a subtree node (isRoutingEntry() would return true); 
    * otherwise, return false.
    */
   bool isDefined() const;

   /*!
    * Return true if the entry holds a data object and lives in a leaf 
    * node of a tree; otherwise, return false.
    */
   bool isDataEntry() const;

   /*!
    * Return true if the entry holds a subtree node; otherwise, return false.
    */
   bool isRoutingEntry() const;

   /*!
    * Return const reference to key associated with this entry.
    */
   const MTreeKey& getKey() const;

   /*!
    * Set pointer to node holding this entry.
    *
    * When assertion checking is on, throws an assertion if pointer is null.
    */
   void setNode(MTreeNodePtr node); 

   /*!
    * Return pointer to node owning this entry.
    */
   MTreeNodePtr getNode() const; 

   /*!
    * Set position of entry in node that owns it.
    *
    * When assertion checking is on, throws an assertion if node pointer 
    * is null or if position is out of range. 
    */
   void setPositionInNode(int pos); 

   /*!
    * Return position of entry in node that owns it.
    */
   int getPositionInNode() const;

   /*!
    * Set subtree node associated with this entry, set entry type to
    * ROUTING_ENTRY and reset data object pointer and id to undefined.
    *
    * When assertion checking is on, throws an assertion if pointer is null.
    */
   void setSubtreeNode(MTreeNodePtr subtree);

   /*!
    * Return pointer to root node of subtree of this entry.
    *
    * If the entry is in a leaf node (i.e., the entry is a data entry), 
    * then a null pointer is returned.  When assertion checking is active, 
    * an assertion is thrown if the entry is not a routing entry.
    */
   MTreeNodePtr getSubtreeNode() const; 

   /*!
    * Set data object id associated with this entry, set entry type to
    * DATA_ENTRY and reset subtree node and id to undefined.
    *
    * Note: No error checking on validity of node id.
    */
   void setDataObjectId(int object_id);

   /*!
    * Return id of data object of this entry.
    *
    * If the entry is not in a leaf node (i.e., the entry is a routine 
    * entry), then an undefined object id is returned. 
    */
   int getDataObjectId() const; 

   /*!
    * Check consistency of entry.
    *
    * @return false if no inconsistencies found; otherwise true.
    *
    * @param my_node Pointer to node owning this entry.  If 
    *                null an assertion will be thrown.
    * @param stream Output stream to which all inconsistencies
    *               will be printed.
    */
   bool checkConsistency(MTreeNodePtr my_node,
                         ostream& stream) const;

   /*!
    * Print MTree entry data to the specified output stream. 
    */
   void printClassData(ostream& stream) const;

   //@{
   //! @name Methods that make operations on key more convenient.

   /*!
    * Compute and return distance between this entry and given entry.
    */
   double computeDistanceTo(MTreeEntryPtr other) const;

   /*!
    * Set distance to parent of this entry to given value.
    */
   void setDistanceToParent(double distance);

   /*!
    * Return distance between entry objects.
    */
   double getDistanceToParent() const;

   /*!
    * Set radius of this entry to given value.
    */
   void setRadius(double radius);

   /*!
    * Return radius of this entry.
    */
   double getRadius() const;

   /*!
    * Return point of this entry.
    */
   MetricSpacePointPtr getPoint() const;

   //@}

private:
   // The following are not implemented
   MTreeEntry();
   MTreeEntry(const MTreeEntry&);
   void operator=(const MTreeEntry&);

   /*
    * Enumerated type for possible entry types.
    */
   enum MTreeEntryType { UNDEFINED_ENTRY = -1,
                         ROUTING_ENTRY = 0,
                         DATA_ENTRY = 1 };

   /*
    * Enumerated type indicating type of entry;
    * entry is undefined by default, it is changed
    * when either data entry or subtree node is set
    */
   MTreeEntryType  d_entry_type;        

   /*
    * Metric space region for entry; initialized in constructor.
    */
   MTreeKey      d_key;                 

   /*
    * Pointer to node containing this entry; null by default.
    */
   MTreeNodePtr  d_my_node;
   int           d_my_node_id;

   /*
    * Position of this entry among entries in owning node; 
    * Invalid position by default; set when entry is added to node. 
    */
   int           d_my_position_in_node; 

   /*
    * Pointer to node holding subtree of entry; null by default
    *  - should be non-null only when entry is in an internal node
    *    (i.e., entry type is ROUTING_ENTRY)
    *  - should be null otherwise (i.e., when entry type is 
    *    UNDEFINED_ENTRY or entry is in a leaf node
    *    and entry type is DATA_ENTRY)
    */ 
   MTreeNodePtr  d_subtree_node;        
   int           d_subtree_node_id;

   /*
    * Identifier of data object associated with entry; 
    * it is undefined by default.
    *  - should be defined only when entry is in a leaf node
    *    (i.e., entry type is DATA_ENTRY)
    *  - should be undefined otherwise (i.e., when entry type is 
    *    UNDEFINED_ENTRY or entry is in an internal node
    *    and entry type is ROUTING_ENTRY)
    */
   int            d_data_object_id;

};

#ifndef DEBUG_NO_INLINE
#include "MTreeEntry.I"
#endif
#endif
