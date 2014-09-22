//
// File:        InterpolationDataBase.h
// Package:     kriging coupler
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: Abstract class encapsulating interpolation database
//

#ifndef included_krigcpl_InterpolationDataBase_h
#define included_krigcpl_InterpolationDataBase_h

#include <iosfwd>
#include <string>
#include <vector>

namespace krigcpl {
    
    /*!
     * @brief InterpolationDataBase is an abstract class encapsulating
     * operations of obtaining the value at a given point location.
     *
     */
    
    class InterpolationDataBase {

      //
      // types
      //

    public:

      enum { LOST_HINT_FLAG = 0,
	     USED_HINT_FLAG,
	     MODEL_SIZE_LIMIT_FLAG,
	     MODEL_INSERT_LIMIT_FLAG,
	     MODEL_NEW_START_FLAG,
	     NEW_HINT_FLAG,
	     NUMBER_FLAGS
      };

      //
      // methods
      //
    public:
      /*!
       * Construction
       *
       * @param pointDimension The dimension of the point space.
       * @param valueDimension The dimension of the value space.
       */
      InterpolationDataBase(int pointDimension,
			    int valueDimension);
      
      /*!
       * Construction
       *
       * @param pointDimension The dimension of the point space.
       * @param valueDimension The dimension of the value space.
       * @param fileName File name to be used for seeding the database.
       */
      InterpolationDataBase(int                 pointDimension,
			    int                 valueDimension,
			    const std::string & fileName);
      
      /*!
       * Virtual destructor.
       */
      virtual ~InterpolationDataBase() = 0;
      
      /*!
       * Compute interpolated value at a point.
       *
       * @param value Pointer for storing the value. Size of at least
       *              _valueDimension assumed.
       * @param hint  Reference to integer. This variable may be used to 
       *              provide a hint to the database. May be updated upon 
       *              return. An example of use would be the case in which 
       *              the database contains a collection of models. hint could
       *              then contain an index to the appropriate model.
       * @param point Pointer for accesing the point. Needs to have the size
       *              of at least _pointDimension.
       * @param flags Handle to a container for storing flags related
       *              to the inner workings of the interpolation database.
       * @param error_estimate Error estimate
       *
       * @return true if the interpolation successful; false otherwise. 
       */
      virtual bool interpolate(double            * value,
			       int               & hint,
			       const double      * point,
			       std::vector<bool> & flags,
                               double            & error_estimate) = 0;
      /*!
       * Compute interpolated value at a point.
       *
       * @param value Pointer for storing the value. Size of at least
       *              _valueDimension assumed.
       * @param hintList Pointer to an array of hints.
       * @param numberHints Size of hintList array.
       * @param oVIndexForMin Integer index of the value component used to
       *                      check for the best hint-based model.
       * @param hintUsed Reference to integer id of a hint actually used.
       * @param hint  Reference to integer. This variable may be used to 
       *              provide a hint to the database. May be updated upon 
       *              return. An example of use would be the case in which 
       *              the database contains a collection of models. hint could
       *              then contain an index to the appropriate model.
       * @param point Pointer for accesing the point. Needs to have the size
       *              of at least _pointDimension.
       * @param flags Handle to a container for storing flags related
       *              to the inner workings of the interpolation database.
       * @param error_estimate Error estimate
       *
       * @return true if the interpolation successful; false otherwise. 
       */      
      virtual bool interpolate(double            * value,
			       const int         * hintList,
			       int                 numberHints, 
			       int                 oVIndexForMin, 
			       int                & hintUsed,
			       const double       * point,
			       std::vector<bool>  & flags,
                               double             & error_estimate) = 0;

      /*!
       * Compute interpolated value at a point.
       *
       * @param value Pointer for storing the value. Size of at least
       *              _valueDimension assumed.
       * @param gradient Pointer for storing gradient of the value wrt.
       *                 point evaluated at the point.
       * @param hint  Reference to integer. This variable may be used to 
       *              provide a hint to the database. May be updated upon 
       *              return. An example of use would be the case in which 
       *              the database contains a collection of models. hint could
       *              then contain an index to the appropriate model.
       * @param point Pointer for accesing the point. Needs to have the size
       *              of at least _pointDimension.
       * @param flags Handle to a container for storing flags related
       *              to the inner workings of the interpolation database.
       * @param error_estimate Error estimate
       * @return true if the interpolation successful; false otherwise. 
       */
      virtual bool interpolate(double            * value,
			       double            * gradient,
			       int               & hint,
			       const double      * point,
			       std::vector<bool> & flags,
                               double            & error_estimate) = 0;
      
      /*!
       * Compute interpolated value at a point.
       *
       * @param value Pointer for storing the value. Size of at least
       *              _valueDimension assumed.
       * @param gradient Pointer for storing gradient of the value wrt.
       *                 point evaluated at the point.
       * @param hintList Pointer to an array of hints.
       * @param numberHints Size of hintList array.
       * @param oVIndexForMin Integer index of the value component used to
       *                      check for the best hint-based model.
       * @param hintUsed Reference to integer id of a hint actually used.
       * @param hint  Reference to integer. This variable may be used to 
       *              provide a hint to the database. May be updated upon 
       *              return. An example of use would be the case in which 
       *              the database contains a collection of models. hint could
       *              then contain an index to the appropriate model.
       * @param point Pointer for accesing the point. Needs to have the size
       *              of at least _pointDimension.
       * @param flags Handle to a container for storing flags related
       *              to the inner workings of the interpolation database.
       * @param error_estimate Error estimate
       *
       * @return true if the interpolation successful; false otherwise. 
       */      
      virtual bool interpolate(double            * value,
			       double            * gradient,
			       const int         * hintList,
			       int                 numberHints, 
			       int                 oVIndexForMin, 
			       int                & hintUsed,
			       const double       * point,
			       std::vector<bool>  & flags,
                               double             & error_estimate) = 0;

      /*!
       * Insert the point-value pair into the database.
       *
       * @param hint   A hint for the database.
       * @param point  Pointer to point data. Needs to have the size of
       *               at least _pointDimension.
       * @param value  Pointer to value data. Needs to have the size of 
       *               at least _valueDimension
       * @param gradient Pointer to gradient of the value wrt. point.
       * @param flags Handle to a container for storing flags related
       *              to the inner workings of the interpolation database.
       *
       */
      virtual void insert(int               & hint,
			  const double      * point,
			  const double      * value,
			  const double      * gradient,
			  std::vector<bool> & flags) = 0;
      
      /*!
       * Insert the point-value pair into the database.
       *
       * @param hintUsed Reference for integer hint used in the insertion.
       * @param point  Pointer to point data. Needs to have the size of
       *               at least _pointDimension.
       * @param value  Pointer to value data. Needs to have the size of 
       *               at least _valueDimension
       * @param gradient Pointer to gradient of the value wrt. point.
       * @param hintList Pointer to an array of hints.
       * @param numberHints Size of hintList array.
       * @param forceInsert A flag to force creation of a new model
       *                    containing a single point-value pair.
       * @param flags Handle to a container for storing flags related
       *              to the inner workings of the interpolation database.
       */
      virtual void insert(int               & hintUsed,
			  const double      * point,
			  const double      * value,
			  const double      * gradient,
			  const int         * hintList,
			  int                 numberHints,
			  bool                forceInsert,
			  std::vector<bool> & flags) = 0;
      
      virtual double interpolateSpecificModel(double            * value,
                                              double            * gradient,
                                              int               & model,
                                              const double      * point,
                                              std::vector<bool> & flags) = 0;

      /*! 
       * Get the point dimension
       *
       * @return Dimensionality of the point space.
       */ 
      int getPointDimension() const;
      
      /*! 
       * Get the value dimension
       *
       * @return Dimensionality of the value space.
       */
      int getValueDimension() const;

      /*!
       * Get the number of performance statistic data collected
       *
       * @return Number of data collected.
       */
      virtual int getNumberStatistics() const = 0;

      /*!
       * Provide performance statistic data collected so far
       *
       * @param stats A handle to an array.
       * @param size  Size of the stats array. 
       */
      virtual void getStatistics(double * stats,
				 int      size) const = 0;
      
      /*! 
       * Provide string descriptions of statistics data.
       *
       * @return An STL-vector of strings.
       */
      virtual std::vector<std::string> getStatisticsNames() const = 0;

      /*!
       * Print DB stats
       * 
       * @param outputStream Stream to be used for output.
       */

      virtual void printDBStats(std::ostream & outputStream) = 0;

      /*!
       * Swap out some objects in order to free up memory
       */

      virtual void swapOutObjects() const = 0;

    protected:

    private:
      // Not implemented
      InterpolationDataBase(const InterpolationDataBase &);
      const InterpolationDataBase & operator=(const InterpolationDataBase&);


      //
      // data
      //
      
    public:

    protected:

    private:
      int _pointDimension;
      int _valueDimension;

    };

}

#endif // included_krigcpl_InterpolationDataBase_h
