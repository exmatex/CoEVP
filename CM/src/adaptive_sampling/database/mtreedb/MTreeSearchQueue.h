//
// File:        MTreeSearchQueue.h
// Package:     MTree database
// 
// 
// 
// Description: Utility class used in MTree searches
//

#ifndef included_MTreeSearchQueue
#define included_MTreeSearchQueue

#ifndef included_list
#define included_list
#include <list>
using namespace std;
#endif

#ifndef included_MTreeSearchNode
#include "MTreeSearchNode.h"
#endif

/*!
 * @brief MTreeSearchQueue implements a simple ordered list 
 * of MTreeSearchNode objects used in MTree nearest-neighbor searches.
 * 
 * This class is used in MTree search operations and should not
 * be used for other stuff.
 */
 
class MTreeSearchQueue
{
public:
   /*!
    * Default ctor creates an empty queue. 
    */
   MTreeSearchQueue();

   /*!
    * Dtor for MTreeSearchQueue. 
    */
   virtual ~MTreeSearchQueue();

   /*!
    * Return true if queue is empty; otherwise return false.
    */
   bool empty() const;

   /*!
    * Return first MTreeSearchNode object in queue.
    */
   MTreeSearchNode getFirst() const;
 
   /*!
    * Remove first MTreeSearchNode object in queue.
    */
   void removeFirst();
 
   /*!
    * Return last MTreeSearchNode object in queue.
    */
   MTreeSearchNode getLast() const;

   /*!
    * Remove last MTreeSearchNode object in queue.
    */
   void removeLast();
 
   /*!
    * Insert MTreeSearchNode object in queue in proper location.
    */
   void insert(MTreeSearchNode& search);
 
   /*!
    * Clear queue of all MTreeSearchNode objects.
    */
   void clear();
 
private:
   // The following are not implemented
   MTreeSearchQueue(const MTreeSearchQueue&);
   void operator=(const MTreeSearchQueue&);

   list<MTreeSearchNode> d_queue;
 
};

#ifndef DEBUG_NO_INLINE
#include "MTreeSearchQueue.I"
#endif
#endif
