//
// File:        KrigingInterpolationKeyDB.h
// 
// Revision:    $Revision$
// Modified:    $Date$
// Description: Kriging interpolation using a keyed database.
//

#ifndef included_krigcpl_KrigingInterpolationKeyDB_h
#define included_krigcpl_KrigingInterpolationKeyDB_h

#include <unordered_map>

#ifndef included_krigcpl_InterpolationDataBase_h
#include "base/InterpolationDataBase.h"
#endif 

#ifndef included_krigalg_InterpolationModelFactory_h
#include <base/InterpolationModelFactory.h>
#endif

#ifndef included_base_DB
#include <base/DB.h>
#endif // included_base_DB

namespace krigcpl {

    typedef std::unordered_map<std::string, krigalg::InterpolationModelPtr> InterpolationModelDataBase;

    /*!
     * @brief Concrete implementation of the InterpolationDataBase class using
     * kriging as the interpolation model
     */

    class KrigingInterpolationKeyDB :
      public InterpolationDataBase {

    public:

      /*!
       * Construction.
       * 
       * @param pointDimension The dimension of the point space.
       * @param valueDimension The dimension of the value space.
       * @param modelFactory   Handle to a factory for creating interpolation models
       * @param keyDB          Handle to a database for keys indexing the interpolation database
       * @param modelDB        Handle to an interpolation model database
       * @param maxKrigingModelSize Maximum number of point/value pairs in
       *                            a single kriging model.
       * @param maxNumberSearchModels Maximum number of kriging models to be
       *                              pulled out of the database during the
       *                              search for the best kriging model.
       * @param useHint Whether or not to make use of hint data.
       * @param meanErrorFactor The value of the coefficient multiplying
       *                        the mean square error.
       * @param tolerance Requested tolerance.
       * @param maxQueryPointModelDistance The maximum distance between the 
       *                                   query point and the model for which
       *                                   interpolation is still attempted. 
       * @param agingThreshold Time threshold for object aging.
       */
      KrigingInterpolationKeyDB(int    pointDimension,
                                int    valueDimension,
                                const krigalg::InterpolationModelFactoryPointer  & modelFactory,
                                DB&    keyDB,
                                InterpolationModelDataBase& modelDB,
                                int    maxKrigingModelSize,
                                int    maxNumberSearchModels,
                                bool   useHint,
                                double meanErrorFactor,
                                double tolerance,
                                double maxQueryPointModelDistance,
                                int    agingThreshold);
      /*!
       * Construction.
       * 
       * @param pointDimension The dimension of the point space.
       * @param valueDimension The dimension of the value space.
       * @param modelFactory   Handle to a factory for creating interpolation models
       * @param keyDB          Handle to a database for keys indexing the interpolation database
       * @param modelDB        Handle to an interpolation model database
       * @param maxKrigingModelSize Maximum number of point/value pairs in
       *                            a single kriging model.
       * @param maxNumberSearchModels Maximum number of kriging models to be
       *                              pulled out of the database during the
       *                              search for the best kriging model.
       * @param useHint Whether or not to make use of hint data.
       * @param meanErrorFactor The value of the coefficient multiplying
       *                        the mean square error.
       * @param tolerance Requested tolerance.
       * @param maxQueryPointModelDistance The maximum distance between the 
       *                                   query point and the model for which
       *                                   interpolation is still attempted.
       * @param agingThreshold Time threshold for object aging.
       * @param directoryName Name of the directory to use for storage
       *                           of disk DB data.
       * @param fileName File name to be used for seeding the database.
       */
      KrigingInterpolationKeyDB(int    pointDimension,
                                int    valueDimension,
                                const krigalg::InterpolationModelFactoryPointer  & modelFactory,
                                DB&    keyDB,
                                InterpolationModelDataBase& modelDB,
                                int    maxKrigingModelSize,
                                int    maxNumberSearchModels,
                                bool   useHint,
                                double meanErrorFactor,
                                double tolerance,
                                double maxQueryPointModelDistance,
                                int    agingThreshold,
                                const std::string & directoryName,
                                const std::string & fileName);
      
      /*!
       * Destruction.
       */
      ~KrigingInterpolationKeyDB();

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
                               double            & error_estimate);

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
                               double             & error_estimate);

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
       *
       * @return true if the interpolation successful; false otherwise. 
       */
      virtual bool interpolate(double            * value,
			       double            * gradient,
			       int               & hint,
			       const double      * point,
			       std::vector<bool> & flags,
                               double            & error_estimate);	

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
                               double             & error_estimate);

      virtual double interpolateSpecificModel(double            * value,
                                              double            * gradient,
                                              int               & model,
                                              const double      * point,
                                              std::vector<bool> & flags);

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
       */
      virtual void insert(int               & hint,
			  const double      * point,
			  const double      * value,
			  const double      * gradient,
			  std::vector<bool> & flags);

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
			  std::vector<bool> & flags);
      
      /*!
       * Get the number of performance statistic data collected
       *
       * @return Number of data collected.
       */
      virtual int getNumberStatistics() const;

      /*!
       * Provide performance statistic data collected so far
       *
       * @param stats A handle to an array.
       * @param size  Size of the stats array. 
       */
      virtual void getStatistics(double * stats,
				 int      size) const;

      /*! 
       * Provide string descriptions of statistics data.
       *
       * @return An STL-vector of strings.
       */
      virtual std::vector<std::string> getStatisticsNames() const;

      /*!
       * Print DB stats
       * 
       * @param outputStream Stream to be used for output.
       */

      virtual void printDBStats(std::ostream & outputStream);

#if 0
      /*!
       * Output object given predicate.
       *
       * @param predicate Unary predicate to test for objects.
       */
      template<typename DataPredicate>
      bool writeObjects(const DataPredicate & predicate) const
      {
	return _krigingModelDB.writeObjects(predicate);
      }
#endif

      /*!
       * Swap out some objects in order to free up memory
       */

      virtual void swapOutObjects() const;

      /*!
       * Perform a query for k-closest interpolants.
       *
       * @param numberInterpolants The number of closest interpolants 
       *        requested (k).
       *
       * @param point The point.
       *
       * @return An STL-vector containg kriging interpolants.
       */

      std::vector<krigalg::InterpolationModelPtr>
	getKrigingModels(int            numberModels,
			 const double * point);

      /*!
       * Perform a range query given a point.
       *
       * @param radius Radius of the range query.
       * @param point  Point to use.
       *
       * @return An STL-vector containg kriging interpolants.
       */

      std::vector<krigalg::InterpolationModelPtr>
	getKrigingModels(double         radius,
			 const double * point);

    protected:

    private:
      // Not implemented
      KrigingInterpolationKeyDB(const KrigingInterpolationKeyDB &);
      const KrigingInterpolationKeyDB & operator=(const KrigingInterpolationKeyDB&);

      //
      // data
      //

    public:

    protected:
      
      krigalg::InterpolationModelFactoryPointer _modelFactory;

      DB&                         _keyDB;
      InterpolationModelDataBase& _modelDB;

      int    _maxKrigingModelSize;
      int    _maxNumberSearchModels;
      bool   _useHint;
      double _meanErrorFactor;
      double _tolerance;
      double _maxQueryPointModelDistance;
      
      //
      // various statistics
      //

      int _numberKrigingModels;
      int _numberPointValuePairs;

      //
      // kriging model aging threshold in seconds. kriging models that
      // have not been access in the last threshold seconds will be
      // moved to a disk.
      //

      int _agingThreshold;

    private:

    };

}

#endif // included_krigcpl_KrigingInterpolationKeyDB_h
