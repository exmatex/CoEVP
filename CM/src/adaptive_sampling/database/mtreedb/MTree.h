//
// File:        MTree.h
// Package:     MTree database
// 
// 
// 
// Description: Main Mtree index structure class.
//

#ifndef included_MTree
#define included_MTree

#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif

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

#ifndef included_DB
#include <base/DB.h>
#endif
#ifndef included_MTreeDataStore
#include "MTreeDataStore.h"
#endif
#ifndef included_MTreeEntry
#include "MTreeEntry.h"
#endif
#ifndef included_MTreeNode
#include "MTreeNode.h"
#endif
#ifndef included_DBObject
#include <base/DBObject.h>
#endif
#ifndef included_MTreeSearchResult
#include "MTreeSearchResult.h"
#endif
#ifndef included_MTreeQuery
#include "MTreeQuery.h"
#endif
#ifndef included_MTreeLevelStatistic
#include "MTreeLevelStatistic.h"
#endif

#ifndef NULL
#define NULL (0)
#endif


class MTreeLevelStatistic;
class MTreeSearchResult;

/*!
 * @brief MTree is the main class providing the capabilities of an 
 *        MTree index structure, which is a dynamic, paged, metric tree.  
 *        Specifically, it supports insertion, deletion, and querying of 
 *        data objects that are described relative to each other in terms 
 *        of points in a general metric space.
 *
 * Typical usage of an MTree object involves several operations:
 *
 * -# Create the tree calling the ctor; e.g., MTree my_tree("MyTree").
 *
 * -# Initialize the tree in one of two ways. 
 *               To set up a new tree from scratch, call the 
 *               initializeCreate() method.  To set up a new tree 
 *               and initialize it with data read from existing files, 
 *               call the initializeOpen() method.
 *
 * -# Insert data. Data objects are inserted using the insertObject() 
 *               method. Note that in addition to providing the data 
 *               object to the method, a point representing the object 
 *               and an object radius are required.  The object radius
 *               may be zero if the user desires to index points in the
 *               metric space rather than spherical regions.  However,
 *               each object must be specified by a unique point.
 * 
 * -# Query data. The data in the tree can be searched in two ways. 
 *               The searchKNN() method returns the k-nearest neighbors 
 *               of a given point, and the searchRange() method returns 
 *               all data objects within a given distance of a given point.
 *
 * -# Finalize the tree structure by calling the finalize() method.  
 *               This will write all MTree index structure state and
 *               data object information to files for retrieval later.
 * 
 * -# Destroy the tree by calling the dtor explicitly or letting the
 *               tree go out of scope.
 * 
 * @see DBObject
 * @see MetricSpacePoint
 * @see MTreeNode
 * @see MTreeDataStore
 */

class MTree
   : public DB
{
public:
   friend class MTreeNode;
   friend class MTreeQuery;

   /*!
    * Ctor for MTree object sets tree name and default node size and
    * turns on error checking and sets error logging stream if specified.
    */
   MTree(const string& tree_name,
         ostream* error_log_stream = (ostream*)NULL,
         bool do_error_checking = false);

   /*!
    * Dtor for MTree objects.
    */
   virtual ~MTree();

   //@{
   //! @name Methods for initializing/finalizing this MTree object.

   /*!
    * Initialize MTree to empty state (containing no data objects).
    *
    * Either this method, or the initializeOpen() method must be called
    * before any tree operations involving data objects may be performed.
    * 
    * An MTree object can only be initialized once.  Calling this method 
    * more than once or after calling he initializeOpen() method will 
    * result in an unrecoverable error.
    *
    * @param directory_name Const reference to string indicating the
    *                     directory into which all data files
    *                     associated with this tree will be written.
    * @param file_prefix  Const reference to string indicating the
    *                     prefix for all data files associated with
    *                     this tree.
    * @param obj_factory  Const reference to factory that creates 
    *                     concrete data objects to be indexed by tree.
    */
   virtual void initializeCreate(const string& directory_name,
                                 const string& file_prefix,
                                 const DBObjectFactory& obj_factory);

   /*!
    * Initialize MTree to state contained in existing data files.
    *
    * Either this method, or the initializeCreate() method must be called
    * before any tree operations involving data objects may be performed.
    *
    * An MTree object can only be initialized once.  Calling this method
    * more than once or after calling he initializeCreate() method will
    * result in an unrecoverable error.
    *
    * @param directory_name Const reference to string indicating the
    *                     directory into which all data files associated 
    *                     with this tree will be read and written.
    * @param file_prefix  Const reference to string indicating the
    *                     prefix for all data files associated with
    *                     this tree. 
    * @param obj_factory  Const reference to factory that creates 
    *                     concrete data objects to be indexed by tree.
    */
   virtual void initializeOpen(const string& directory_name,
                               const string& file_prefix,
                               const DBObjectFactory& obj_factory);

   /*!
    * Finalize MTree index structure.
    *
    * This method will write all MTree index structure state and
    * data object information to files for retrieval later (e.g.,
    * initialization of another tree object). 
    */
   virtual void finalize(); 

   //@}
  
   //@{
   //! @name Methods for setting parameters for this MTree object.

   /*!
    * Set maximum number of node entries for all nodes in MTree.  
    *
    * If this method is not called, the default value given by
    * DEFAULT_MAX_NODE_ENTRIES is used. 
    *
    * Note that we constrain the tree so that each node has the 
    * the same maximum number of entries.  To successfully the 
    * max number of entries, this method must be called before 
    * any data object is interted into the tree.  If this method 
    * is called after an object has been inserted, the max 
    * number of entries will not be changed.   
    * 
    * Also, the given max number of entries must be at least two; 
    * if not, the max number of entries will not be changed.
    * 
    * @param max_entries  Integer maximum number of node entries.
    */
   void setMaxNodeEntries(int max_entries);

   /*!
    * Set promotion method for all MTreeNodes in this tree, except 
    * the root node.  If this method is not called, the default
    * will be used (see description of MTreeNodePromotionMethod enum 
    * type in MTreeNode header file).  If this method is called after 
    * an object has been inserted, the node promotion method 
    * will not be changed.
    *
    * @param method enum type MTreeNode::MTreeNodePromotionMethod value.
    */
   void setNodePromotionMethod(
      MTreeNode::MTreeNodePromotionMethod method);

   /*!
    * Set promotion method for the root node.  If this method is not called, 
    * the default will be used (see description of MTreeNodePromotionMethod 
    * type in MTreeNode header file).  If this method is called after an 
    * object has been inserted, the node promotion method will not be 
    * changed. 
    *
    * @param method enum type MTreeNode::MTreeNodePromotionMethod value.  
    *               Method choice for root node is more restrictive than 
    *               for other nodes. If a disallowed value is given, the 
    *               default method is used (see description of 
    *               MTreeNodePromotionMethod enum type in MTreeNode 
    *               header file).
    */
   void setRootNodePromotionMethod(
       MTreeNode::MTreeNodePromotionMethod method);

   /*!
    * Set partition method for all MTreeNodes in this tree.  If this method 
    * is not called, the default will be used (see description of 
    * MTreeNodePartitionMethod type in MTreeNode header file).  If this 
    * method is called after an object has been inserted, the node partition 
    * method will not be changed.
    *
    * @param method enum type MTree::MTreeNodePartitionMethod value.
    * @param min_utilization Optional double value indicating the minimum 
    *        node utilization fraction.  This value is used only for 
    *        unbalanced partition strategies. When partitioning, entries 
    *        will be moved to a new node until the fraction of total 
    *        entries in the new node is greater than or equal to this value.
    *        The default value is 0.5 which will result in a balanced 
    *        split. Valid values are in the range (0.0, 0.5]. If a value 
    *        outside of this range is given, the default is used.
    */
   void setNodePartitionMethod(
      MTreeNode::MTreeNodePartitionMethod method,
      double min_utilization = 0.5);

   //@}
  
   //@{
   //! @name Methods for inserting, deleting and retrieving objects in this MTree object.

   /*!
    * Insert object into tree and set the identifier of the object.
    * 
    * Note that in addition to providing the data object to the method, 
    * a point representing the object and an object radius are required.  
    * The object radius may be zero if the user desires to index points 
    * in the metric space rather than spherical regions.  However, each 
    * object must be specified by a unique point.
    *
    * To avoid external tampering with database contents, this method 
    * produces internal deep copies of the given point and data objects.
    *
    * @param object  Reference to data object.  The reference is non-const
    *                since the method sets the object identifier to match
    *                the identifier of the database copy.
    * @param point   Const reference to center point of object.
    * @param radius  Double radius of object about center point. When 
    *                assertion checking is on, assertion will result if 
    *                value is less than 0.
    */
   virtual void insertObject(DBObject& object,
                             const MetricSpacePoint& point,
                             double radius);
   
   /*!
    * Get copy of object indexed by tree given object identifier.
    *
    * @param object_id  Integer identifier of object to delete from tree.
    *                If this is not a valid id for an object indexed by
    *                the tree, the method will return a null pointer.
    */ 
   virtual DBObjectPtr getObject(int object_id) const;

   /*!
    * Delete object from tree. 
    *
    * This method attempts to minimize the impact of the object deletion on 
    * the tree index structure.  However, the method will move around entries 
    * in the tree to avoid empty nodes and to limit the number of nodes that
    * have only one entry.   Depending on pattern of object deletions (e.g., 
    * frequency and in which regions of the tree), they can result in 
    * significant tree restructuring which can be expensive and may 
    * negatively impact the performance of other operations such as 
    * object insertions and searches.
    *
    * @param object_id  Integer identifier of object to delete from tree.
    *                If this is not a valid id for an object indexed by
    *                the tree, the method will do nothing. 
    */
   virtual void deleteObject(int object_id);

   //@}
  
   /*!
    * Search tree for "k" nearest neighbors of given query point.
    * 
    * To avoid external tampering with database contents, this method 
    * returns deep copies of points and data objects within the set of
    * search results.
    *
    * @param results Reference to vector of MTreeSearchResult objects.
    *                The size of the vector will be the number of search
    *                results, up to a maximum of the given "K".  This
    *                vector will be cleared on entry to method so any
    *                pre-existing data it contains will be lost.  Also,
    *                the vector of results will be sorted with respect 
    *                to increaing distance from the query point.
    * @param query_point  Const reference to search query point.
    * @param k_neighbors  Integer maximum neighbors to return.  If 
    *                this value is <= 0, the result vector will be empty.
    * @param make_safe  Optional boolean indicating whether the search
    *                will return copies of objects in data base (true),
    *                or actual objects in the database (false).  The
    *                default value is false.  Passing a value of true
    *                prevents any possibility of data objects being
    *                modified outside of the database.  In the dafault
    *                case, if the user casts away the const-ness of 
    *                any data object in a search result, future database
    *                operations may produce unexpected behavior.
    */
   virtual void searchKNN(vector<DBSearchResult>& results,
                          const MetricSpacePoint& query_point,
                          int k_neighbors,
                          bool make_safe = false);

   /*!
    * Search tree for all data objects within given distance of given 
    * query point.
    *
    * Note that if the given radius is zero, this method will return 
    * all data objects whose representative metric space regions contain
    * the given query point.  If the radius is greater than zero, this
    * method will return all data objects whose regions intersect the query
    * region.
    *
    * To avoid external tampering with database contents, this method
    * returns deep copies of points and data objects within the set of
    * search results.
    *
    * @param results Reference to list of MTreeSearchResult objects.
    *                The size of the list will be the number of search
    *                results.  This list will be cleared on entry to 
    *                method so any pre-existing data it contains will be 
    *                lost.  Also, the list of results will be sorted with
    *                respect to increasing distance from query point. 
    * @param query_point  Const reference to point at center of search.
    * @param radius  Double radius value indicating search region about 
    *                query point.  If this value is < 0, the result list 
    *                will be empty.
    * @param make_safe  Optional boolean indicating whether the search
    *                will return copies of objects in data base (true),
    *                or actual objects in the database (false).  The
    *                default value is false.  Passing a value of true
    *                prevents any possibility of data objects being
    *                modified outside of the database.  In the dafault
    *                case, if the user casts away the const-ness of 
    *                any data object in a search result, future database
    *                operations may produce unexpected behavior.
    */
   virtual void searchRange(list<DBSearchResult>& results,
                            const MetricSpacePoint& query_point,
                            double radius,
                            bool make_safe = false);

   //@}
  
   //@{
   //! @name Methods for checking MTree consistency and obtaining MTree level statistics.

   /*!
    * Check consistency of entire tree structure including all 
    * nodes and entries.
    * @return false if no inconsistencies found; otherwise true.
    *
    * @param stream Output stream to which all inconsistencies
    *               will be printed.
    */
   bool checkConsistency(ostream& stream) const;

   /*!
    * Return number of levels in this MTree object.
    */
   int getNumberLevels() const;

   /*!
    * Calculate statistics associated with levels in MTree structure.
    */
   void calculateLevelStatistics();

   /*!
    * Return statistics object associated with given MTree level.
    * Note that levels numbers increase from leaf nodes (level == zero)
    * up to root of tree (level == getNumberLevels()).  
    */
   const MTreeLevelStatistic* 
   getLevelStatistics(int level_number) const; 

   //@}

   //@{
   //! @name Methods for obtaining MTree insertion statistics.

   /*!
    * Get total number of object insertions performed by this tree;
    *     i.e., number of times insertObject() function was called.
    */
   int getTotalInsertCount() const;

   /*!
    * Get total number of computations of distance between objects
    *     during all object insertions performed by this tree.
    */
   int getTotalInsertDistanceCount() const;

   /*!
    * Get number of computations of distance between objects
    *     during last object insertion performed by this tree.
    */
   int getLastInsertDistanceCount() const;

   //@}
    
   //@{
   //! @name Methods for obtaining MTree deletion statistics.

   /*!
    * Get total number of object deletions performed by this tree;
    *     i.e., number of times deleteObject() function was called.
    */
   int getTotalDeleteCount() const;

   /*!
    * Get total number of computations of distance between objects
    *     during all object deletions performed by this tree.
    */
   int getTotalDeleteDistanceCount() const;

   /*!
    * Get number of computations of distance between objects
    *     during last object deletion performed by this tree.
    */
   int getLastDeleteDistanceCount() const;

   //@}
    
   //@{
   //! @name Methods for obtaining MTree nearest neighbor search statistics.

   /*!
    * Get total number of nearest neighbor searches performed by this tree;
    *     i.e., number of times searchKNN() function was called.
    */
   int getTotalKNNSearchCount() const;

   /*!
    * Get total number of computations of distance between objects
    *     during all nearest neighbor searches performed by this tree.
    */
   int getTotalKNNSearchDistanceCount() const;

   /*!
    * Get number of computations of distance between objects
    *     during last nearest neighbor search performed by this tree.
    */
   int getLastKNNSearchDistanceCount() const;

   //@}
    
   //@{
   //! @name Methods for obtaining MTree range search statistics.

   /*!
    * Get total number of range searches performed by this tree;
    *     i.e., number of times searchRange() function was called.
    */
   int getTotalRangeSearchCount() const;

   /*!
    * Get total number of computations of distance between objects
    *     during all range searches performed by this tree.
    */
   int getTotalRangeSearchDistanceCount() const;

   /*!
    * Get number of computations of distance between objects
    *     during last range search performed by this tree.
    */
   int getLastRangeSearchDistanceCount() const;

   //@}
    
   //@{
   //! @name Methods for printing information about this MTree.

   /*!
    * Print calculated set of all MTree statistics.  To ensure
    * that printed results correspond to current state of tree,
    * the calculateStatistics() member function should be called
    * before this print routine.
    */
   void printAllTreeStatistics(ostream& stream) const;

   /*!
    * Print calculated MTree statistics for given level number.  
    * To ensure that printed results correspond to current state of tree,
    * the calculateStatistics() member function should be called
    * before this print routine.
    */
   void printTreeLevelStatistics(int level_number,
                                 ostream& stream) const;

   /*!
    * Print all MTree operation count statistics.  This includes
    * insertions, deletions, searches, and distance computations.
    */
   void printAllTreeOperationStatistics(ostream& stream) const;

   /*!
    * Print summary of MTree node configuration to given output stream. 
    */
   void printNodeSummary(ostream& stream) const;

   /*!
    * Print entire MTree data structure to the specified output stream. 
    */
   void printClassData(ostream& stream) const;

   /*!
    * Print MTree data store contents to the specified output stream.
    */
   void printDataStoreClassData(ostream& stream) const;

   //@}

   //@{
   //! @name Methods for writing/reading MTree to/from database file.

   /*!
    * Write entire MTree index structure to given database. 
    * Return boolean true if write is successful; false otherwise.
    * 
    * In general, this method should only be called by the data
    * store object associated with this tree object.
    */
   bool putToDatabase(toolbox::DatabasePtr write_tree_DB) const 
   { return true; }

   /*!
    * Read entire MTree index structure from given database. 
    * Return boolean true if read is successful; false otherwise.
    * 
    * In general, this method should only be called by the data
    * store object associated with this tree object.
    */
   bool getFromDatabase(toolbox::DatabasePtr read_tree_DB) 
   { return true; }

   /*!
    * Write objects to disk. Return boolean tru if write is successful
    * and false otherwise.
    */
   template<typename DataPredicate>
   bool writeObjects(const DataPredicate & predicate) const;
  
   //@}

   //@{
   //! @name Methods for testing...will be removed in future.

   // Testing .....
   bool testWriteObject(int object_id)
   {
      return( d_data_store.writeDataObject(object_id) );
   }

   // Testing .....
   bool testWriteAllObjects()
   {
      return( d_data_store.writeAllDataObjects() );
   }

   // Testing .....
   bool testReadObject(int object_id)
   {
      return( d_data_store.readDataObject(object_id) );
   }

   // Testing .....
   bool testReadAllObjects()
   {
      return( d_data_store.readAllDataObjects() );
   }

   virtual void outputStats(std::ostream & outputStream);

   //@}

private:
   // The following are not implemented
   MTree(const MTree&);
   void operator=(const MTree&);

   /*
    * Private methods to allow node objects to access algorithm
    * parameters.
    */
   MTreeNode::MTreeNodePromotionMethod getRootNodePromotionMethod() const;
   MTreeNode::MTreeNodePromotionMethod getNodePromotionMethod() const;
   MTreeNode::MTreeNodePartitionMethod getNodePartitionMethod() const;
   double getMinNodeUtilization() const;

   /*
    * Private method to recurse to child nodes used in range search.
    */
   void searchRangeRecursive(list<MTreeSearchResult>& results,
                             const MTreeQuery& query,
                             MTreeNodePtr node) const;

   /*
    * Private method to select node into which entry will be inserted.
    * The method accepts a starting node and traverses the tree from that
    * node to a suitable leaf node, which it returns.
    */
   MTreeNodePtr pickNode(MTreeNodePtr start_node,
                         MTreeEntryPtr entry) const;

   /*
    * Private method to split given node, if over-full, and continue to
    * split nodes nodes up to the root node as needed.
    */
   void split(MTreeNodePtr node);

   /*
    * Private method to split root node, if over-full.
    */
   void splitRootNode();

   /*
    * Private method to create and initialize root node of tree.
    */
   void createRootNode();

   /*
    * Private method to set node ownership in data store for 
    * objects in a leaf node.
    */
   void mapDataObjectsToNode(MTreeNodePtr node);

   /*
    * Private method to get smart pointer to root node.
    */
   MTreeNodePtr getRootNode() {return d_root_node;}

   /*
    * Default value for max number of entries in a node.  Not sure
    * what this should be; I pulled this value out of my butt.
    */
   enum { DEFAULT_MAX_NODE_ENTRIES = 7 };

   /*
    * Data store that holds nodes and data objects for tree
    */
   mutable MTreeDataStore  d_data_store;

   /*
    * Root node of tree; null by default.  Set when first object inserted.
    */
   MTreeNodePtr  d_root_node;

   /*
    * Max number of entries allowed in each node of tree;
    * Set to enum value by default.  
    * Note: Can be changed from default via member function
    *       only before first object is inserted.
    */
   int d_max_node_entries;

   /*
    * Root node and non-root node promotion and partition methods 
    * shared by all nodes in this tree; set to defaults in constructor.
    * Changed via methods described above.
    */
   MTreeNode::MTreeNodePromotionMethod d_root_node_promotion_method;
   MTreeNode::MTreeNodePromotionMethod d_node_promotion_method;
   MTreeNode::MTreeNodePartitionMethod d_node_partition_method;
   double d_min_node_utilization;

   /*
    * Data members and functions for gathering tree statistics.
    */
   void addNodeToLevelCount(int level);
   void removeNodeFromLevelCount(int level);
   void clearLevelStatistics();
   void incrementDistanceComputeCount();
   void clearOperationCountStatistics();

   enum MTreeDistanceType { DISTANCE_RANGE_QUERY = 0,
                            DISTANCE_KNN_QUERY = 1,
                            DISTANCE_INSERT = 2, 
                            DISTANCE_DELETE = 3,
                            DISTANCE_UNDEFINED = 5 };

   MTreeDistanceType d_current_distance_type;

   int d_num_range_queries;
   int d_num_knn_queries;
   int d_num_inserts;
   int d_num_deletes;

   int d_total_distance_comps_in_range_queries;
   int d_total_distance_comps_in_knn_queries;
   int d_total_distance_comps_in_inserts;
   int d_total_distance_comps_in_deletes;

   int d_num_distance_comps_in_last_range_query;
   int d_num_distance_comps_in_last_knn_query;
   int d_num_distance_comps_in_last_insert;
   int d_num_distance_comps_in_last_delete;

   vector<int> d_number_nodes_in_level;
   vector<MTreeLevelStatistic*> d_level_statistics;

};

#ifndef DEBUG_NO_INLINE
#include "MTree.I"
#endif

#include "MTree.t.h"

#endif
