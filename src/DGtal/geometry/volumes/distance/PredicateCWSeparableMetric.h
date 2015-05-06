/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file PredicateCWSeparableMetric.h
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2015/04/17
 *
 * Header file for module PredicateCWSeparableMetric.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(PredicateCWSeparableMetric_RECURSES)
#error Recursive header files inclusion detected in PredicateCWSeparableMetric.h
#else // defined(PredicateCWSeparableMetric_RECURSES)
/** Prevents recursive inclusion of headers. */
#define PredicateCWSeparableMetric_RECURSES

#if !defined PredicateCWSeparableMetric_h
/** Prevents repeated inclusion of headers. */
#define PredicateCWSeparableMetric_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include <cmath>
#include "DGtal/math/BasicMathFunctions.h"
#include "DGtal/kernel/CInteger.h"
#include "DGtal/kernel/CSpace.h"
#include "DGtal/kernel/CInteger.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class PredicateCWSeparableMetric
  /**
   * Description of template class 'PredicateCWSeparableMetric' <p>
   * \brief Aim:
   */
  template <typename TSpace, typename TWeight,
    typename TPromoted=DGtal::int64_t>
    class PredicateCWSeparableMetric
    {
    // ----------------------- Standard services ------------------------------
  public:
  
  
  ///Copy the space type
  typedef TSpace Space;
  BOOST_CONCEPT_ASSERT(( concepts::CSpace<TSpace> ));
  
  ///Type for points
  typedef typename Space::Point Point;
  ///Type for points
  typedef typename Point::Coordinate Abscissa;
  ///Type for vectors
  typedef typename Space::Vector Vector;
  
  ///Type for internal (integer) distance values
  typedef TPromoted Promoted;
  
  // Type of the weights
  //typedef TPromoted Weight;
  typedef TWeight Weight;
  
  // Type for "real" distance values
  typedef long double Distance;

  BOOST_CONCEPT_ASSERT(( concepts::CInteger<Promoted> ));
  
  ///Type for Value (alias)
  typedef TPromoted Value;
  
  /**
   * Constructor.
   */
  PredicateCWSeparableMetric();

  
  /**
   * Destructor.
   */
  ~PredicateCWSeparableMetric();

  /**
   * Copy constructor.
   */
  PredicateCWSeparableMetric ( const PredicateCWSeparableMetric & /*other*/ ) {}

  /**
   * Assignment.
   * @return a reference on 'this'.
   */
  PredicateCWSeparableMetric & operator= ( const PredicateCWSeparableMetric & /*other*/ ) { return *this;}
  

  
    // ----------------------- Interface --------------------------------------
  public:

  // ----------------------- CPowerMetric --------------------------------------
  
  
  /**
   *  Return the restricted CW distance of a point @a aPoint and a weighted
   *  point (@a aQ,@a aWq): returns 0 if @a aPoint is inside the ball of center @a aQ and radius @a aWq[0], 1 if @a aPoint is outside the ball of center @a aQ and radius @a aWq[1], the CW distance otherwise.
   *
   * @param aPoint a point
   * @param aQ a second point
   * @param aWq weights of the second point (a vector of two numbers such that aWq[0] <= aWq[1])
   *
   * @return the CW distance between aPoint and (Q,WQ)
   */
  Distance restrictedCWDistance(const Point &aPoint,
		    const Point &aQ,
		    const Weight &aWq) const;


  
  /**
   *  Return the CW distance of a point @a aPoint and a weighted
   *  point (@a aQ,@a aWq)
   *
   * @param aPoint a point
   * @param aQ a second point
   * @param aWq weights of the second point (a vector of two numbers such that aWq[0] <= aWq[1])
   *
   * @return the CW distance between aPoint and (Q,WQ)
   */
  Distance CWDistance(const Point &aPoint,
		    const Point &aQ,
		    const Weight &aWq) const;
  
  
  /**
   * Given an origin and two points, this method decides which one
   * is closest to the origin. This method should be faster than
   * comparing distance values.
   *
   * @param origin the origin
   * @param first  the first point
   * @param wF the first point weights
   * @param second the second point
   * @param wS the second point weights
   *
   * @return a Closest enum: FIRST, SECOND or BOTH.
     */
  DGtal::Closest closestCW(const Point &origin,
			      const Point &first,
			      const Weight &wF,
			      const Point &second,
			      const Weight &wS) const;
  
  
  // ----------------------- CPowerSeparableMetric --------------------------------------
  
  /**
     * Given three weighted sites (u,v,w) and a straight segment
     * [startingPoint,endPoint] along dimension dim, we detect if the
     * cells of @a u and @a w @e strictly hide the cell of @a v on the
     * straight line according to the CW metric.
     *
     * This method is in @f$ O(log(n))@f$ if @a n is the size of the
     * straight segment. For @f$ l_2@f$ metric (p=2), the method is in
     * @f$ O(1)@f$.
     *
     * @pre u,v and w must be such that u[dim] < v[dim] < w[dim]
     *
     * @param u a site
     * @param wu weights
     * @param v a site
     * @param wv weights
     * @param w a site
     * @param ww a weights
     * @param startingPoint starting point of the segment
     * @param endPoint end point of the segment
     * @param dim direction of the straight line
     *
     * @return true if (u,w) hides v.
     */
    bool hiddenByCW(const Point &u,
		    const Weight &wu,
		    const Point &v,
		    const Weight &wv,
		    const Point &w,
		    const Weight &ww,
		    const Point &startingPoint,
		    const Point &endPoint,
		    const typename Point::UnsignedComponent dim) const;
  
  

  

  /**
   * Writes/Displays the object on an output stream.
   * @param out the output stream where the object is written.
   */
  void selfDisplay ( std::ostream & out ) const;
  
  /**
   * Checks the validity/consistency of the object.
   * @return 'true' if the object is valid, 'false' otherwise.
   */
  bool isValid() const;
  
  // ------------------------- Protected Datas ------------------------------
  private:
  
  /**
   * Compute the Lp distance without the computation of the power
   * 1/p. I.e. only @f$ \sum |p_i- q_i|^p@f$ is given.
   *
   * @param aP a first point
   * @param aQ a second point
   *
   * @return the power p of the l_p distance between aP and aQ.
   */
  Promoted DistanceRepresentation(const Point &aP, const Point &aQ) const;
  
  /**
   * Perform a binary search on the interval [lower,upper] to
   * detect the mid-point between u and v according to the weighted l_p
   * distance.
   *
   * @param udim coordinate of u along dimension dim
   * @param vdim coordinate of v along dimension dim
   * @param nu  partial distance of u (sum of |xj-x_i|^p) discarding
   * the term along the dimension dim
   * @param nv partial distance of v (sum of |xj-x_i|^p) discarding
   * the term along the dimension dim
   * @param lower interval lower bound
   * @param upper interval upper bound
   *
   * @return the Voronoi boundary point coordinates along dimension dim.
   */
  Abscissa binarySearchHidden(const Abscissa &udim,
			      const Abscissa &vdim,
			      const Promoted &nu,
			      const Weight &wu,
			      const Promoted &nv,
			      const Weight &wv,
			      const Abscissa &lower,
			      const Abscissa &upper) const;
  
  

  // ------------------------- Private Datas --------------------------------
  private:
  
  // ------------------------- Hidden services ------------------------------
  /* protected: */
  
  /* /\** */
  /*  * Constructor. */
  /*  * Forbidden by default (protected to avoid g++ warnings). */
  /*  *\/ */
  /* PredicateCWSeparableMetric(); */
  
  /* private: */
  
  /* /\** */
  /*  * Copy constructor. */
  /*  * @param other the object to clone. */
  /*  * Forbidden by default. */
  /*  *\/ */
  /* PredicateCWSeparableMetric ( const PredicateCWSeparableMetric & other ); */
  
  /* /\** */
  /*  * Assignment. */
  /*  * @param other the object to copy. */
  /*  * @return a reference on 'this'. */
  /*  * Forbidden by default. */
  /*  *\/ */
  /* PredicateCWSeparableMetric & operator= ( const PredicateCWSeparableMetric & other ); */
  
    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class PredicateCWSeparableMetric


  /**
   * Overloads 'operator<<' for displaying objects of class 'PredicateCWSeparableMetric'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'PredicateCWSeparableMetric' to write.
   * @return the output stream after the writing.
   */
  /* template <typename T> */
  /* std::ostream& */
  /* operator<< ( std::ostream & out, const PredicateCWSeparableMetric<T> & object ); */

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/volumes/distance/PredicateCWSeparableMetric.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined PredicateCWSeparableMetric_h

#undef PredicateCWSeparableMetric_RECURSES
#endif // else defined(PredicateCWSeparableMetric_RECURSES)
