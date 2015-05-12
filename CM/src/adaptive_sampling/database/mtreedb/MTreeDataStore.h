//
// File:        MTreeDataStore.h
// Package:     MTree database
// 
// 
// 
// Description: Manager class for MTree node allocation and storage of data objects.
//

#ifndef included_MTreeDataStore
#define included_MTreeDataStore

#ifndef included_iostream
#define included_iostream
#include <iostream>
using namespace std;
#endif
#ifndef included_String
#include <string>
using namespace std;
#define included_String
#endif
#ifndef included_list
#define included_list
#include <list>
using namespace std;
#endif
#ifndef included_vector
#define included_vector
#include <vector>
using namespace std;
#endif


#ifdef HAVE_PKG_hdf5
#ifndef included_toolbox_HDFDatabase
#include "toolbox/database/HDFDatabase.h"
#endif
#endif

#ifndef included_DBObject
#include <base/DBObject.h>
#endif
#ifndef included_DBObjectFactory
#include <base/DBObjectFactory.h>
#endif

#ifndef NULL
#define NULL (0)
#endif


class MTree;
class MTreeNode;
typedef std::shared_ptr<MTreeNode> MTreeNodePtr;
class MTreeEntry;
typedef std::shared_ptr<MTreeEntry> MTreeEntryPtr;

/*!
 * @brief MTreeDataStore maintains data objects and nodes for an MTree 
 *        index structure.  Specifically, it stores a mapping between 
 *        data objects and leaf nodes in the tree.  It also manages writing
 *        and reading an mtree index structure and its object data to and 
 *        from HDF5 files.  Note that the data store is designed to be 
 *        used by the MTree class and is not typically used directly.
 * 
 * Typical usage of the data store by the MTree index structure involves
 * several operations.  Note that these operations are invoked by the
 * MTree object associated with the data store and are not called explicitly
 * by a user.
 *
 * -# Creation.  The tree creates a data store object that will managing 
 *               its data using the default MTreeDataStore constructor.
 *
 * -# Initialization.  The tree initializes its data store in one of two 
 *               ways. To set up a new empty data store from scratch, 
 *               call the create() method.  To set up a new data store 
 *               and initialize it with data  read from existing files, 
 *               call the open() method.
 * 
 * -# Manipulate data. Nodes and objects are added, retrieved, removed, 
 *               etc. from the data store using the appropriate methods. 
 *
 * -# Finish.    Call the close() method when done with the data store 
 *               to write all unwritten data as well as the associated
 *               MTree index structure to files.
 * 
 * @see MTree
 * @see MTreeNode
 * @see DBObject
 * @see DBObjectFactory
 */

class MTreeDataStore
{
public:
   /*!
    * Default ctor for MTreeDataStore object.
    */
   MTreeDataStore();

   /*!
    * Dtor for MTreeDataStore objects.
    */
   virtual ~MTreeDataStore();

   /*!
    * Return true if data store object has been initialized, 
    * and false otherwise.
    *
    * The data store is initialized by calling either of the
    * methods create() or open(). 
    */
   bool isInitialized() const;

   /*!
    * Initialize empty data store and set it up to write data
    * to specified directory and files.
    * 
    * If the data store is already initialized, then method issues
    * a warning and does nothing.
    *
    * @param mtree        Bare pointer to Mtree for which 
    *                     data store will manage data.
    * @param obj_factory  Const bare pointer to factory used to create 
    *                     copies of data objects.
    * @param directory_name Const reference to string indicating the
    *                     directory into which all data files 
    *                     associated with this data store will be
    *                     written.  
    * @param file_prefix  Const reference to string indicating the
    *                     prefix for all data files associated with
    *                     this data store.
    *
    * When assertion checking is on, passing an empty string or a null
    * pointer will throw an assertion.
    */ 
   void create(MTree* mtree,
               const DBObjectFactory* obj_factory,
               const string& directory_name,
               const string& file_prefix);

   /*!
    * Initialize data store to contents of data files associated
    * with given prefix in specified directory.
    * 
    * If the data store is already initialized, then method issues
    * a warning and does nothing.
    *
    * @param mtree        Bare pointer to Mtree for which 
    *                     data store will manage data.
    * @param obj_factory  Const bare pointer to factory used to create 
    *                     copies of data objects.
    * @param directory_name Const reference to string indicating the
    *                     directory into which all data files
    *                     associated with this data store will be
    *                     written.
    * @param file_prefix  Const reference to string indicating the
    *                     prefix for all data files associated with
    *                     this data store.
    *
    * When assertion checking is on, passing an empty string or a null
    * pointer will throw an assertion.
    *
    * If either the data files or the directory do not exist, the 
    * an unrecoverable error will result.
    */ 
   void open(MTree* mtree,
             const DBObjectFactory* obj_factory,
             const string& directory_name,
             const string& file_prefix); 

   /*!
    * Write all un-written data to files and close data store.
    * 
    * It is important to note that after calling this method, no
    * other methods will do anything until the data store open()
    * or create method is called.
    */
   void close();

   /*!
    * Add leaf node information to data store and set its 
    * leaf node identifier.
    *
    * If the node has its leaf node identifier set, and the node 
    * matches the node in the data store with that identifier, the
    * routine does nothing.  If the node has its leaf node identifier 
    * set, and the node does not match the leaf node in the data store 
    * with that identifier, then an unrecoverable error results.  Such
    * a situation indicates a potentially substantial problem.
    * 
    * This method does not assign to the node data objects in the data 
    * store that the leaf node may own. Each data object must be explicitly 
    * mapped to the node that owns it by calling the setNodeOwningObject() 
    * method.
    * 
    * @param node  Pointer to leaf node to be added.  When assertion
    *              checking is on, an assertion is thrown if the
    *              pointer is null.
    */ 
   void addLeafNode(MTreeNodePtr node);

   /*!
    * Remove leaf node from data store. 
    *
    * If the node is not a leaf node, or does not have its node 
    * identifier set, the routine does nothing.
    * 
    * @param node  Pointer to leaf node to be removed from data store.
    *              When assertion checking is on, an assertion is thrown 
    *              if the pointer is null.
    */ 
   void removeLeafNode(MTreeNodePtr node);

   /*!
    * Determine whether given leaf node id is associated with a  
    * valid leaf node in data store. 
    *
    * @return Boolean true if id is associated with a leaf node in 
    *         the data store; otherwise false.
    *
    * @param  leaf_node_id  Integer identifier of leaf node.
    */
   bool isValidLeafNodeId(int leaf_node_id) const;

   /*!
    * Create copy of given data object, set object identifier for
    * original object and copy and add the copy to the data store.
    *
    * If the given object has its identifier set, and the object 
    * matches an object in the data store with that identifier, the
    * routine does nothing.  If the object has its identifier set,
    * and the object does not match the object in the data store with
    * that identifier, then an unrecoverable error results.  Such
    * a situation indicates a potentially substantial problem.
    *
    * This method does not assign the data object to the node that
    * owns it if the object is owned by a node.  Each data object must 
    * be explicitly mapped to the node that owns it  by calling 
    * the setNodeOwningObject() method.
    *
    * @param object  Reference to object to be added.  The reference
    *                is non-const since method sets object identifier.
    */
   void addObject(DBObject& object);

   /*!
    * Remove object from data store.
    *
    * @param object_id  Integer id of object to remove from data store.
    *              If object id does not match the id of some object in
    *              the data store, then this routine does nothing.
    */
   void removeObject(int object_id); 

   /*!
    * Determine whether object with given id lives in data store. 
    *
    * @return Boolean true if id is associated with an object in data store;
    *         otherwise false.
    *
    * @param  object_id  Integer identifier of data object.
    */
   bool isValidObjectId(int object_id) const;

   /*!
    * Return deep copy of data object in store with given identifier.
    *
    * @return Pointer to data object.  This will be null if identifier 
    *         does not correspond to any data object in the data store.
    * 
    * @param  object_id  Integer identifier of data object.
    */
   DBObjectPtr getObjectCopy(int object_id);

   /*!
    * Return pointer to data object in store with given identifier.
    *
    * @return Pointer to data object.  This will be null if identifier 
    *         does not correspond to any data object in the data store.
    * 
    * @param  object_id  Integer identifier of data object.
    *
    * Note that it is the responsibility of the user calling this
    * routine to not modify the data object returned.  If it is
    * altered in any way, then unexpected behavior can result.
    */
   DBObjectPtr getObjectPtr(int object_id);

   /*!
    * Get node owning object with given identifier. 
    *
    * Successful execution of this method requires that the object with
    * the given id has been added to the data store and that it is mapped
    * to some node in the data store.  If neither of these conditions has
    * been met, the returned node pointer will be null.
    *
    * @return      Pointer to node owning object.  Will be null if object
    *              with given id is not in data store or object is in 
    *              store but not mapped to a node.
    * 
    * @param  object_id  Integer identifier of data object.
    */
   MTreeNodePtr getNodeOwningObject(int object_id) const;

   /*!
    * Set information mapping data object to node.
    *
    * Successful execution of this method requires that the node and 
    * data object have already been added to the data store using the 
    * addNode() and addObject() methods, respectively. If the given node 
    * does not own the given entry, if either the node or object have not
    * been added to the store, or if the entry does not contain a data 
    * object (e.g., the node is not a leaf node and the entry is a routing
    * entry), then the method does nothing.
    *
    * @return Boolean true if operation succeeded; otherwise false.
    *
    * @param node  Pointer to node to which data object will be assigned.
    * @param entry Pointer to entry holding data object.
    *
    * When assertion checking is on, an assertion is thrown if either 
    * pointer is null.
    */
   bool setNodeOwningObject(MTreeNodePtr node,
                            MTreeEntryPtr entry);

   /*!
    * Clear information mapping data objects to node.
    *
    * @return Boolean true if operation succeeded; otherwise false.
    *
    * This routing should be called for each affected node before 
    * remapping objects to nodes.  Note that this operation is
    * distinct from simply adding a new object to a node.
    *
    * @param node  Pointer to node to which data objects will be unassigned.
    *
    * When assertion checking is on, an assertion is thrown if the node
    * pointer is null.
    */
   bool clearOwnObjectIds(MTreeNodePtr node);

   /*!
    * Write each data object in data store to proper HDF5 file.
    * 
    * Note that this routine will only write those data objects 
    * that have not already been written to a file.
    *
    * @return Boolean true if write operation succeeded; otherwise false.
    */
   bool writeAllDataObjects();

   /*!
    * Write data objects associated with given leaf node to proper HDF5 file.
    * 
    * Note that this routine will only write those data objects 
    * that have not already been written to a file.
    *
    * @return Boolean true if write operation succeeded; otherwise false.
    *
    * @param  node  Pointer to leaf node.  If not a valid leaf node, 
    *         function does nothing.  When assertion checking is on, 
    *         an assertion will result if pointer is null.
    */
   bool writeNodeDataObjects(MTreeNodePtr node);

   /*!
    * Write data objects satisfying given predicate.
    *
    * @return Boolean true if write operation succeeded; otherwise false.
    *
    * @param Unary predicate.
    */
   template<typename DataPredicate>
   bool writeObjects(const DataPredicate & predicate);
   
   /*!
    * Write object with given integer identifier to HDF5 file if it does
    * not currently exist in a file.
    *
    * @return Boolean true if write operation succeeded; otherwise false.
    *
    * @param  object_id  Integer identifier of data object. If not a valid
    *         data object id, function does nothing.
    */
   bool writeDataObject(int object_id);

   /*!
    * Read each data object in data store from proper HDF5 file.
    *
    * Note that this routine will only read those data objects 
    * that do not currently exist in memory.
    *
    * @return Boolean true if read operation succeeded; otherwise false.
    *
    * @param rebuild_from_open Boolean true if reading objects during
    *        initialization of data store from exisiting data files 
    *        (i.e., called from open()).  False by default.  When true
    *        only those objects marked as residing in memory (e.g., 
    *        when data store object from which files were constructed
    *        wxisted) are read from files.  When true, all objects 
    *        are read from files.
    */
   bool readAllDataObjects(bool rebuild_from_open = false);

   /*!
    * Read data objects associated with given leaf node from proper HDF5 file.
    * 
    * Note that this routine will only read those data objects 
    * that do not currently exist in memory.
    *
    * @return Boolean true if read operation succeeded; otherwise false.
    *
    * @param  node  Pointer to leaf node.  If not a valid leaf node, 
    *         function does nothing.  When assertion checking is on, 
    *         an assertion will result if pointer is null.
    */
   bool readNodeDataObjects(MTreeNodePtr node);

   /*!
    * Read object with given integer identifier from HDF5 file if it does
    * not currently exist in memory.
    *
    * @return Boolean true if read operation succeeded; otherwise false.
    *
    * @param  object_id  Integer identifier of data object. If not a valid
    *         data object id, function does nothing.
    */
   bool readDataObject(int object_id);

   /*!
    * Print all data store data (except the actual data objects 
    * themselves) to the specified output stream.
    */
   void printClassData(ostream& stream) const;

private:
   // The following are not implemented
   MTreeDataStore(const MTreeDataStore&);
   void operator=(const MTreeDataStore&);

#ifdef HAVE_PKG_hdf5
   /*
    * Private method to write object with given integer identifier 
    * to given file index and database if it does not currently exist 
    * in a file.  Return boolean true if write operation succeeded; 
    * otherwise false.
    */
   bool writeDataObjectHelper(
      int object_id,
      toolbox::HDFDatabasePtr write_DB, 
      int file_index);

   /*
    * Private method to read object with given integer identifier 
    * from given database if it currently exists in a file, or if
    * boolean flag to force read is true.  Return boolean true if 
    * read operation succeeded; otherwise false.
    */
   bool readDataObjectHelper(
      int object_id,
      toolbox::HDFDatabasePtr read_DB,
      bool force_read = false);
#endif

   /*
    * Private method to generate file name and file index 
    * for next object write.  Returns boolean true if name
    * and index refer to a new (non-existent file); otherwise, false.
    */
   bool getFileNameAndIndexForObjectWrite(string& file_name,
                                          int& file_index);

   /*
    * Private method to generate file name string for given index.
    */
   string getObjectFileName(int file_index) const;

   /*
    * Private method to generate database name for given object id.
    */
   string getObjectDatabaseName(int object_id) const;

   /*
    * Private method to generate full path mtree index file name.
    */
   string getFullPathMTreeIndexFileName() const;

   /*
    * Private method to generate full path data store file name.
    */
   string getFullPathDataStoreFileName() const;

   /*
    * Private method to generate object file prefix.
    */
   string getObjectFilePrefix() const;

#ifdef HAVE_PKG_hdf5
   /*
    * Private methods to read/write object file information
    * to/from file. 
    */
   bool writeObjectFileInfo(toolbox::HDFDatabasePtr write_DB) const;
   bool readObjectFileInfo(toolbox::HDFDatabasePtr read_DB);

   /*
    * Private methods to read/write object information
    * to/from file. 
    */
   bool writeObjectInfo(toolbox::HDFDatabasePtr write_DB) const;
   bool readObjectInfo(toolbox::HDFDatabasePtr read_DB);
#endif

   /*
    * Private utility class for managing information about leaf
    * nodes owning data objects maintained by data store.
    */   
   class LeafNodeInfo 
   {
      public:
         LeafNodeInfo(MTreeNodePtr node);

         ~LeafNodeInfo();

         MTreeNodePtr getNode() const;

         vector<int>& getObjectIds();

         void clearOwnObjectIds(); 

         void setOwnObjectId(int obj_id, int pos); 

         void setInFile(bool value);
         bool getInFile() const;

         void setInMemory(bool value);
         bool getInMemory() const;

         void printClassData(ostream& stream) const;

      private:
         LeafNodeInfo(const LeafNodeInfo&);
         LeafNodeInfo& operator = (const LeafNodeInfo&);
  
         MTreeNodePtr d_node_ptr;
         bool d_in_file;
         bool d_in_memory;

         vector<int> d_object_ids;
   };

   /*
    * Private utility class for managing information about data
    * objects maintained by data store.
    */
   class ObjectInfo
   {
      public:
         ObjectInfo(DBObjectPtr object);

         ~ObjectInfo();

         DBObjectPtr getObject() const;

         void resetObjectPtr();

         void setObjectPtr(DBObjectPtr obj);

         void setOwnerLeafNodeInfo(int leaf_node_id,
                                   int entry_position);

         int getOwnerLeafNodeId() const;
         int getOwnerLeafNodeEntryPosition() const;

         void setInFile(bool value);
         bool getInFile() const;

         void setInMemory(bool value);
         bool getInMemory() const;

         void setFileIndex(int index);
         int getFileIndex() const;

         void printClassData(ostream& stream) const;

      private:
         ObjectInfo(const ObjectInfo&);
         ObjectInfo& operator = (const ObjectInfo&);

         DBObjectPtr d_object_ptr;
         int    d_owner_leaf_node_id;
         int    d_owner_leaf_node_entry_position;
         bool   d_in_file;
         bool   d_in_memory;
         int    d_object_file_index;
   };

   /*
    * Private utility class for managing object information 
    * in files.
    */   
   class ObjectFileInfo 
   {
      public:
         ObjectFileInfo(const string& file_name);

         ~ObjectFileInfo();

         void addObjectId(int object_id);

         const string& getFileName() const;

         const vector<int>& getObjectIds() const;

         bool fileIsEmpty() const;
         bool fileIsFull() const;

         void printClassData(ostream& stream) const;

      private:
         ObjectFileInfo(const ObjectFileInfo&);
         ObjectFileInfo& operator = (const ObjectFileInfo&);
  
         string d_file_name;
         vector<int> d_object_ids;
   };

   /*
    * enumerated types used for internal data management.
    */
   enum { MTREE_DATA_STORE_FILE_FORMAT_VERSION = 1,
          MTREE_DATA_STORE_OBJECT_FILE_CAPACITY = 25,
          LEAF_NODE_VECTOR_ALLOCATION_CHUNK = 100,
          OBJECT_VECTOR_ALLOCATION_CHUNK = 200 ,
          OBJECT_FILE_VECTOR_ALLOCATION_CHUNK = 50 };

   MTree* d_mtree;
   const DBObjectFactory* d_object_factory;

   bool      d_is_initialized;
   bool      d_is_open;

   /*
    * Top-level directory name and file prefix used to write/read
    * data to/from files.
    */
   string    d_directory_name;
   string    d_file_prefix;

   /*
    * File name for MTree index structure.
    */
   string    d_mtree_index_file_name;

   /*
    * File name for data store.
    */
   string    d_data_store_file_name;

   /*
    * Data members used to map leaf nodes of the tree to data objects. 
    */
   int                   d_num_leaf_nodes;
   vector<LeafNodeInfo*> d_leaf_node_info;
   list<int>             d_recycled_leaf_node_indices;

   /*
    * Data members used to manage data objects.
    */
   int                 d_num_objects;
   vector<ObjectInfo*> d_object_info;
   list<int>           d_recycled_object_indices;

   /*
    * Data members used for managing data object files.
    */
   string                   d_object_file_prefix;
   int                      d_object_file_capacity;
   int                      d_num_objects_in_files;
   vector<ObjectFileInfo*>  d_object_file_info;

};

#ifndef DEBUG_NO_INLINE
#include "MTreeDataStore.I"
#endif

#include "MTreeDataStore.t.h"

#endif
