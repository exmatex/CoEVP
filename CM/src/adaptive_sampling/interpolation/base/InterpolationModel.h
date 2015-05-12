//
// File:        InterpolationModel.h
// Package:     kriging algorithm
// 
// 
// 
// Description: Class implementing kriging interpolation.
//
// $Id: KrigingModel.h,v 1.2 2005/08/24 18:32:33 knap2 Exp $
//
// $Log: KrigingModel.h,v $
// Revision 1.2  2005/08/24 18:32:33  knap2
// Added data members to cache various arrays. Added getError() member
// function.
//
// Revision 1.1  2005/08/23 21:12:40  knap2
// Initial source.
//
//

#if !defined(included_krigalg_InterpolationModel)
#define included_krigalg_InterpolationModel

#include <memory>
 
#ifndef included_vector
#define included_vector
#include <vector>
#endif // included_vector

#ifndef included_iosfwd
#define included_iosfwd
#include <iosfwd>
using namespace std;
#endif

#ifndef included_krigalg_Value
#include "Value.h"
#endif

#ifndef included_krigalg_Point
#include "Point.h"
#endif

#include "TimeRecorder.h"

  //
  // forward declarations
  //
 
  namespace toolbox {
 
    class Database;
 
  }
 
  //
  //
  //

  namespace krigalg {

  //
  // forward declarations
  //

  class InterpolationModel;

  //
  // local types
  //
  
  typedef std::shared_ptr<InterpolationModel> InterpolationModelPtr;

  /*!
   * @brief Interface for a generic interpolation model. The model
   * is a multivariate interpolation model and covers
   * the univariate interpolation as a special case.
   */

  class InterpolationModel {

    //
    // methods
    //

  public:
    /*!
     * Destruction.
     */
    virtual ~InterpolationModel() = 0;

    //
    // features of the interpolation model.
    //

    /*!
     * Check whether the interpolation model provides error. 
     *
     * @return True if the interpolation model provides error 
     *         and false otherwise.
     */
    
    virtual bool hasError() const = 0;

    /*!
     * Check whether the interpolation model provides gradient 
     * estimate along with the interpolated value.
     *
     * @return True if the interpolation model provides gradient
     *         and false otherwise.
     */

    virtual bool hasGradient() const = 0;

    //
    // meta-methods
    //

    /*!
     * @brief Clone the object.
     *
     * @return Pointer to a copy of InterpolationModel.
     */

    virtual InterpolationModel * clone() const = 0;

    /*!
     * Add point/value pair to the interpolation model. In general,
     * after this operation the interpolation model may not be ready 
     * to carry out interpolation. InterpolationModel::build() will
     * have to be called explicitly after the last point has been added.
     *
     * @param point Point to be added.
     * @param values Values associated with the Point..
     *
     * @return boolean true if the point can be added into the model and
     * false otherwise.
     */
    
    virtual bool addPoint(const Point              & point,
			  const std::vector<Value> & values) = 0;
    
    /*!
     * Get number of points in the model.
     *
     * @return integer number of points currently in the model.
     */
       
    virtual int getNumberPoints() const = 0;
 
    /*!
     * Get handle to model points.
     *
     * @return const reference to an STL vector holding the point data.
     */
 
    virtual const std::vector<Point> & getPoints() const = 0;

    /*!
     * Get number of values in the model.
     *
     * @return integer number of values.
     */
 
    virtual int getNumberValues() const = 0;
 
    /*!
     * Get point dimension.
     *
     * @return integer point dimension.
     */
                                                                               
    virtual int getPointDimension() const = 0;
                                                                               
    /*!
     * Get value dimension.
     *
     * @return integer value dimension.
     */
                                                                               
    virtual int getValueDimension() const = 0;

    /*!
     * Check if a model is valid.
     *
     * @return bool value; true if a model is valid and false
     * otherwise.
     */
    
    virtual bool isValid() const = 0;

    /*!
     * Interpolate value at a point.
     *
     * @param valueId index of the value to interpolate.
     * @param point reference to a point to interpolate at.
     *
     * @return value at the point.
     */
    
    virtual Value interpolate(int           valueId,
			      const Point & point) const = 0;

    /*!
     * Estimate error at a point.
     *
     * @param valueId index of the value to interpolate.
     * @param point reference to a point to interpolate at.
     *
     * @return mean squared error at the point.
     */

    virtual Value getMeanSquaredError(int           valueId,
				      const Point & point) const = 0;
    

    /*!
     * Get time recorded in TimeRecorder.
     *
     * @return TimeRecorder with the last recorded time.
     */

    TimeRecorder getTime() const;
  
    /*!
     * InterpolationModel output.
     */

    virtual void putToDatabase(toolbox::Database & db) const = 0;
    virtual void getFromDatabase(toolbox::Database & db) = 0;
    
    virtual void pack(std::vector<double> & packedContainer) const = 0;
    virtual void unpack(const std::vector<double> & packedContainer) = 0;


    friend std::ostream & 
      operator<<(std::ostream             & outputStream,
		 const InterpolationModel & interpolationModel);

  protected:
    /*!
     * Construction.
     */
    InterpolationModel();

    /*!
     * Build the model using accumulated points.
     */
    virtual void build() = 0;    

    /*!
     * @brief print the contents of InterpolationModel to a stream.
     *
     * @param outputStream the stream to use for output.
     */
    virtual void print(std::ostream & outputStream) const = 0;
    

  private:
    //
    // copy construction/assignment
    //
    
    // InterpolationModel(const InterpolationModel &);
    // const InterpolationModel & operator=(const InterpolationModel &);

    //
    // data
    //

  public:

  protected:

    //
    // store last access time
    //
    
    mutable TimeRecorder                                   _timeRecorder;

  private:

  };

}

#endif // included_krigalg_InterpolationModel
